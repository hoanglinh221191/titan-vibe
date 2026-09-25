"""Fingerprint effective topology and protocol in each stored TPR; no MD."""
import argparse,concurrent.futures,hashlib,json,re,subprocess,threading
from pathlib import Path
from dimer_segment_alignment import sha

# Output schedules and run identity are not Hamiltonian differences. All other
# effective input-record fields remain in the protocol fingerprint.
EXEMPT={'nsteps','tinit','init-step','simulation-part','ld-seed','ref-t','nstxout','nstvout','nstfout','nstlog','nstcalcenergy','nstenergy','nstxout-compressed','compressed-x-precision'}
def fingerprint(text):
 start=text.index('\ninputrec:\n')+1;stop=text.index('\nheader:\n',start)
 head=text[start:stop];topstart=text.index('\ntopology:\n')+1;topend=text.index('\nbox (3x3):',topstart)
 topology=text[topstart:topend]
 canonical=[];excluded={}
 for line in head.splitlines():
  line=' '.join(line.split());m=re.match(r'([^ =:]+)\s*(?:=|:)\s*(.*)',line)
  if m and m[1] in EXEMPT:excluded[m[1]]=m[2];continue
  canonical.append(line)
 protocol='\n'.join(canonical)+'\n';topology='\n'.join(' '.join(line.split()) for line in topology.splitlines())+'\n'
 digest=lambda s:hashlib.sha256(s.encode()).hexdigest()
 return {'topology_sha256':digest(topology),'protocol_sha256':digest(protocol),'excluded_run_fields':excluded},topology,protocol

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--dataset',type=Path,required=True);p.add_argument('--work-root',type=Path,required=True);p.add_argument('--gmx',default='gmx');p.add_argument('--workers',type=int,default=3);a=p.parse_args();w=a.work_root.resolve();w.mkdir(parents=True,exist_ok=False);lock=threading.Lock();seen=set()
 def inspect(k,seg):
  path=a.dataset/f'a{k}/{seg}.tpr';digest=sha(path)
  r=subprocess.run([a.gmx,'dump','-s',str(path)],capture_output=True,text=True,timeout=120);r.check_returncode()
  result,top,protocol=fingerprint(r.stdout);result.update(slot=k,segment=seg,tpr_sha256=digest)
  with lock:
   for label,value in [('topology',top),('protocol',protocol)]:
    h=result[label+'_sha256']
    if h not in seen:(w/f'{label}_{h}.txt').write_text(value);seen.add(h)
  return result
 records=[];failures=[]
 with concurrent.futures.ThreadPoolExecutor(max_workers=a.workers) as pool:
  fs={pool.submit(inspect,k,s):(k,s) for k in range(1,49) for s in range(1,21)}
  for f in concurrent.futures.as_completed(fs):
   try:records.append(f.result())
   except Exception as e:failures.append({'slot_segment':fs[f],'error':repr(e)});print('FAIL',fs[f],repr(e),flush=True)
   if (len(records)+len(failures))%60==0:print('TPR',len(records),'verified',len(failures),'failures',flush=True)
 groups={kind:{} for kind in ['topology','protocol']}
 for r in records:
  for kind in groups:groups[kind].setdefault(r[kind+'_sha256'],[]).append([r['slot'],r['segment']])
 result={'status':'PASS_COMMON_TOPOLOGY_AND_PROTOCOL' if not failures and all(len(g)==1 for g in groups.values()) else 'DIFFERENCES_REQUIRE_ASSESSMENT','records':sorted(records,key=lambda r:(r['slot'],r['segment'])),'groups':groups,'failures':failures,'excluded_fields':sorted(EXEMPT),'scope':'Equality at GROMACS dump precision of complete topology and all input-record fields except declared run/output/temperature fields; coordinates, velocities, box and binary-container headers are not compared. This is not bitwise equality of all TPR floating-point values.','script_sha256':sha(__file__)}
 (w/'audit.json').write_text(json.dumps(result,indent=2)+'\n');print(result['status'],{k:len(v) for k,v in groups.items()},flush=True)
 if failures:raise SystemExit(1)

if __name__=='__main__':main()
