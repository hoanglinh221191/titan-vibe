"""Read-only source-energy identity and exact segment-local XTC alignment.

Outputs are audit tables, not DHAM transition counts. No interpolation, nearest
energy substitution, cross-segment pairing or physical-rate claim is made.
"""
import argparse,bisect,collections,concurrent.futures,gzip,hashlib,itertools,json,math,re,subprocess
from pathlib import Path
from xtc_segment_inventory import headers

def sha(path):
 h=hashlib.sha256()
 with Path(path).open('rb') as f:
  for chunk in iter(lambda:f.read(1024*1024),b''):h.update(chunk)
 return h.hexdigest()

def read_xvg(path):
 rows=[];labels={}
 for line in Path(path).read_text().splitlines():
  m=re.search(r'@\s+s(\d+)\s+legend\s+"(.*?)"',line)
  if m:labels[int(m[1])]=m[2]
  if not line.strip() or line.startswith(('#','@')):continue
  row=tuple(map(float,line.split()))
  if len(row)!=4 or not all(map(math.isfinite,row)):raise ValueError('Invalid energy row')
  if rows and row[0]<=rows[-1][0]:raise ValueError('Non-increasing energy time')
  rows.append(row)
 if labels!={0:'Potential',1:'Temperature',2:'Volume'}:raise ValueError(f'Unexpected energy fields: {labels}')
 return rows

def segment_identity(source,cat_lookup):
 matches=[];offsets=collections.Counter();unmatched=ambiguous=0
 for row in source:
  hits=cat_lookup.get(row[1:],[])
  if not hits:unmatched+=1
  elif len(hits)!=1:ambiguous+=1
  else:
   index,global_time=hits[0];delta=round(global_time-row[0],6)
   offsets[delta]+=1;matches.append((row[0],index,global_time))
 if unmatched or ambiguous or len(offsets)!=1:
  raise ValueError(f'Energy identity failed: missing={unmatched} ambiguous={ambiguous} offsets={dict(offsets)}')
 return next(iter(offsets)),matches

def determinant(m):
 return m[0]*(m[4]*m[8]-m[5]*m[7])-m[1]*(m[3]*m[8]-m[5]*m[6])+m[2]*(m[3]*m[7]-m[4]*m[6])

def closed_no_exchange(times,start,end):
 return bisect.bisect_right(times,end)==bisect.bisect_left(times,start)

def log_parameters(path):
 with Path(path).open(errors='replace') as f:head=''.join(itertools.islice(f,400))
 keys=['integrator','dt','nstenergy','nstxout-compressed','tcoupl','pcoupl','pcoupltype','free-energy']
 params={key:(re.search(r'^\s*'+re.escape(key)+r'\s*=\s*(.+)$',head,re.M).group(1).strip() if re.search(r'^\s*'+re.escape(key)+r'\s*=\s*(.+)$',head,re.M) else None) for key in keys}
 m=re.search(r'^\s*ref-t:\s*(.+)$',head,re.M);params['ref_t_K']=list(map(float,m[1].split())) if m else []
 pressure=re.findall(r'ref-p\[\s*\d+\]=\{([^}]+)\}',head)
 params['ref_p_matrix_bar']=[list(map(float,line.split(','))) for line in pressure]
 params['replex_steps']=int(re.search(r'-replex\s+(\d+)',head)[1])
 params['version']=re.search(r'^GROMACS version:\s*(.+)',head,re.M)[1].strip()
 params['header_sha256']=hashlib.sha256(head.encode()).hexdigest()
 return params

def audit_slot(slot,dataset,work,gmx,clock,events_root):
 w=work/f'a{slot}';w.mkdir();folder=dataset/f'a{slot}'
 catpath=folder/f'energy_{slot}.xvg';cat=read_xvg(catpath);lookup=collections.defaultdict(list)
 for i,row in enumerate(cat):lookup[row[1:]].append((i,row[0]))
 segments=[];energy_by_segment={};matched_by_segment={};changed={};covered=set()
 for seg in range(1,21):
  source=folder/f'{seg}.edr';out=w/f'{seg}.xvg'
  digest=sha(source)
  r=subprocess.run([gmx,'energy','-f',str(source),'-o',str(out)],input='Potential\nTemperature\nVolume\n0\n',cwd=w,capture_output=True,text=True,timeout=120)
  (w/f'{seg}_extract.log').write_text(r.stdout+r.stderr);r.check_returncode()
  rows=read_xvg(out);discarded_start=None
  try:offset,matches=segment_identity(rows,lookup)
  except ValueError:
   # Concatenation may retain the previous segment's final sample instead of
   # this segment's distinct t=0 sample. Accept only this precisely located
   # omission, with the retained row already verified against the previous EDR.
   if seg==1 or rows[0][0]!=0:raise
   offset,matches=segment_identity(rows[1:],lookup)
   retained=[(i,row) for i,row in enumerate(cat) if abs(row[0]-offset)<1e-6]
   if len(retained)!=1 or retained[0][0] not in covered:raise ValueError('Unverified boundary energy omission')
   discarded_start={'source_local_row':rows[0],'retained_cat_row_index':retained[0][0],'retained_cat_row':retained[0][1],'reason':'source t=0 omitted; cat retains a row already verified from preceding segment'}
  if any(abs(b[0]-a[0]-20)>1e-6 for a,b in zip(rows,rows[1:])):raise ValueError('Unexpected local energy cadence')
  energy_by_segment[seg]={row[0]:row for row in rows};matched_by_segment[seg]={t:(idx,global_t) for t,idx,global_t in matches}
  covered.update(idx for _,idx,_ in matches)
  params=log_parameters(folder/f'{seg}.log')
  if float(params['dt'])!=.002 or int(params['nstenergy'])!=10000 or int(params['nstxout-compressed'])!=25000 or params['replex_steps']!=2500:raise ValueError('Unexpected timing protocol')
  if not params['ref_t_K'] or len(set(params['ref_t_K']))!=1:raise ValueError('Unknown/multiple temperature groups')
  ev=json.loads((events_root/f'{seg}_exchange_events.json').read_text())
  changed[seg]=[time for step,time,swaps in ev if any(slot-1 in pair for pair in swaps)]
  segments.append(dict(segment=seg,source_edr_sha256=digest,extracted_xvg_sha256=sha(out),source_rows=len(rows),discarded_source_start=discarded_start,energy_global_offset_ps=offset,clock_energy_minus_xtc_ps=offset-clock['blocks'][seg-1]['xtc_offset_ps'],parameters=params))
 if len(covered)!=len(cat):raise ValueError(f'Concatenated energy contains {len(cat)-len(covered)} unaccounted rows')
 # All frames retained with explicit null energy for non-coincident output times.
 counts=collections.Counter();previous=None;maxvol=0;boundary_rows=0
 framehash=hashlib.sha256();xtc=folder/'cat.xtc';xtcstat=xtc.stat()
 with gzip.open(w/'alignment.jsonl.gz','wt',encoding='utf-8') as target:
  blocks=clock['blocks'];segidx=0
  for index,frame in enumerate(headers(xtc)):
   while index>blocks[segidx]['last_frame_zero_based']:segidx+=1
   block=blocks[segidx];seg=segidx+1;local=frame['step']*.002
   if abs(local+block['xtc_offset_ps']-frame['time_ps'])>1e-5 or frame['natoms']!=1254:raise ValueError('Current XTC differs from verified clock map')
   framehash.update(json.dumps(frame,sort_keys=True,separators=(',',':')).encode())
   en=energy_by_segment[seg].get(local);match=matched_by_segment[seg].get(local)
   volume=determinant(frame['cell_nm']);delta=abs(volume-en[3]) if en else None
   if en:maxvol=max(maxvol,delta)
   # XVG volume is printed with 6 decimals; XTC cell entries are float32.
   volume_ok=en is not None and delta<=max(1e-4,2e-6*abs(en[3]))
   boundary=previous is None or previous['segment']!=seg
   candidate=(not boundary and abs(local-previous['local_time_ps']-50)<1e-6 and closed_no_exchange(changed[seg],previous['local_time_ps'],local))
   boundary_rows+=boundary;counts['frames']+=1;counts['exact_energy' if en else 'no_exact_energy']+=1
   counts['volume_pass' if volume_ok else 'volume_mismatch' if en else 'volume_unavailable']+=1
   counts['closed_50ps_exchange_free_candidates']+=candidate
   row=dict(slot_one_based=slot,segment=seg,frame_zero_based=index,xtc_byte_offset=frame['offset'],local_step=frame['step'],local_time_ps=local,xtc_global_time_ps=frame['time_ps'],energy_cat_row_zero_based=match[0] if match else None,energy_global_time_ps=match[1] if match else None,potential_kj_mol=en[1] if en else None,temperature_observed_K=en[2] if en else None,volume_nm3=en[3] if en else None,xtc_volume_nm3=volume,volume_absolute_error_nm3=delta,volume_match=volume_ok if en else None,alignment_status='EXACT_LOCAL_TIME' if en else 'NO_EXACT_ENERGY_SAMPLE',segment_start=boundary,closed_50ps_no_exchange_from_previous=bool(candidate),physical_transition_certified=False)
   target.write(json.dumps(row,separators=(',',':'))+'\n');previous=row
 assert counts['frames']==21240 and boundary_rows==20
 result=dict(status='PASS_IDENTITY_AND_CLOCK' if not counts['volume_mismatch'] else 'IDENTITY_PASS_VOLUME_MISMATCH',slot=slot,cat_energy_sha256=sha(catpath),cat_energy_rows=len(cat),all_cat_energy_rows_accounted=True,segments=segments,counts=dict(counts),max_volume_difference_nm3=maxvol,xtc_bytes=xtcstat.st_size,xtc_mtime_ns=xtcstat.st_mtime_ns,xtc_header_sha256=framehash.hexdigest(),alignment_sha256=sha(w/'alignment.jsonl.gz'),exchange_history_scope='a1/a48 previously verified witnesses; shared global history, not a new all-log audit',physical_dham_authorized=False)
 (w/'audit.json').write_text(json.dumps(result,indent=2)+'\n');print(f"a{slot}: {result['status']} {dict(counts)} max_dV={maxvol:.9g}",flush=True)
 return result

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--dataset',type=Path,required=True);p.add_argument('--work-root',type=Path,required=True);p.add_argument('--clock-map',type=Path,required=True);p.add_argument('--events-root',type=Path,required=True);p.add_argument('--gmx',default='gmx');p.add_argument('--slots',default=','.join(map(str,range(1,49))));p.add_argument('--workers',type=int,default=3);a=p.parse_args()
 w=a.work_root.resolve();w.mkdir(parents=True,exist_ok=False);clocks=json.loads(a.clock_map.read_text());slots=[int(x) for x in a.slots.split(',')]
 if len(set(slots))!=len(slots) or any(x<1 or x>48 for x in slots):raise ValueError('Invalid slots')
 records=[];failures=[]
 with concurrent.futures.ThreadPoolExecutor(max_workers=a.workers) as pool:
  fs={pool.submit(audit_slot,k,a.dataset.resolve(),w,a.gmx,clocks['replicas'][k-1],a.events_root.resolve()):k for k in slots}
  for f in concurrent.futures.as_completed(fs):
   try:records.append(f.result())
   except Exception as error:failures.append({'slot':fs[f],'error':repr(error)});print('FAILED',fs[f],repr(error),flush=True)
 summary={'status':'PASS_IDENTITY_AND_CLOCK' if not failures and all(r['status']=='PASS_IDENTITY_AND_CLOCK' for r in records) else 'NEEDS_REVIEW','slots':sorted(slots),'completed':len(records),'failures':failures,'clock_map_sha256':sha(a.clock_map),'script_sha256':sha(__file__),'records':[{'slot':r['slot'],'status':r['status'],'counts':r['counts'],'audit_sha256':sha(w/f"a{r['slot']}/audit.json")} for r in sorted(records,key=lambda r:r['slot'])]}
 (w/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
 print(summary['status'],len(records),'slots',len(failures),'failures',flush=True)
 if failures:raise SystemExit(1)

if __name__=='__main__':main()
