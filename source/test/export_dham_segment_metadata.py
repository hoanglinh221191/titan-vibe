"""Export verified alignment to the production DHAM counts-only sidecar format."""
import argparse,gzip,hashlib,json
from pathlib import Path
from dimer_segment_alignment import sha
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--eligibility',type=Path,required=True);p.add_argument('--audit-roots',nargs='+',type=Path,required=True);p.add_argument('--events-root',type=Path,required=True);p.add_argument('--work-root',type=Path,required=True);a=p.parse_args();w=a.work_root.resolve();w.mkdir(parents=True,exist_ok=False)
j=json.loads(a.eligibility.read_text());assert j['status']=='TECHNICAL_ALIGNMENT_PASS_PHYSICAL_DHAM_NOT_QUALIFIED'
paths={};audits={};inputs={str(a.eligibility):sha(a.eligibility)}
for root in a.audit_roots:
 for audit in root.glob('a*/audit.json'):
  r=json.loads(audit.read_text());k=r['slot'];assert k not in paths;paths[k]=audit.parent;audits[k]=r
  key=next(q for q in j['inputs_sha256'] if Path(q).resolve()==audit.resolve());assert j['inputs_sha256'][key]==sha(audit)
assert sorted(paths)==list(range(1,49))
events={}
for seg in range(1,21):
 path=a.events_root/f'{seg}_exchange_events.json';key=next(q for q in j['inputs_sha256'] if Path(q).resolve()==path.resolve());assert sha(path)==j['inputs_sha256'][key];events[seg]=json.loads(path.read_text());inputs[str(path)]=sha(path)
out=w/'dham_segments.txt';nframes=nexchange=0
with out.open('x',encoding='utf-8') as f:
 f.write('TITAN_DHAM_SEGMENTS_V1\n# replica indices zero based; segment IDs positive; ps, kJ/mol, nm3, K, bar\n')
 for k in sorted(paths):
  params=audits[k]['segments'][0]['parameters'];f.write(f"replica {k-1} {params['ref_t_K'][0]:.12g} 1\n")
  for seg in range(1,21):
   ev=events[seg];assert ev and ev[0][1]==5
   f.write(f'segment {k-1} {seg} 0 {ev[-1][1]:.12g}\n')
  align=paths[k]/'alignment.jsonl.gz';assert sha(align)==audits[k]['alignment_sha256'];inputs[str(align)]=sha(align)
  with gzip.open(align,'rt',encoding='utf-8') as g:
   for line in g:
    row=json.loads(line);exact=row['alignment_status']=='EXACT_LOCAL_TIME';energy=format(row['potential_kj_mol'],'.17g') if exact else 'NA'
    volume=row['volume_nm3'] if exact else row['xtc_volume_nm3']
    f.write(f"frame {k-1} {row['frame_zero_based']} {row['xtc_global_time_ps']:.12g} {row['segment']} {row['local_time_ps']:.12g} {energy} {volume:.17g} {'exact' if exact else 'absent'}\n");nframes+=1
  for seg,ev in events.items():
   for _,time,pairs in ev:
    if any(k-1 in pair for pair in pairs):f.write(f'exchange {k-1} {seg} {time:.12g}\n');nexchange+=1
(w/'temperatures.txt').write_text(''.join(str(audits[k]['segments'][0]['parameters']['ref_t_K'][0])+'\n' for k in sorted(paths)))
(w/'manifest.json').write_text(json.dumps({'status':'PASS_EXPORT_COUNTS_AUDIT_ONLY','frames':nframes,'exchange_slot_events':nexchange,'metadata_sha256':sha(out),'temperatures_sha256':sha(w/'temperatures.txt'),'inputs_sha256':inputs,'script_sha256':sha(__file__),'physical_rates_validated':False},indent=2)+'\n')
print('PASS export',nframes,'frames',nexchange,'slot exchange events')
