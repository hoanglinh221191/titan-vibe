"""Separate output-group differences from force-field and protocol differences."""
import argparse,collections,hashlib,json
from pathlib import Path
from dimer_segment_alignment import sha

def forcefield_without_output_group(text):
 lines=[]
 for line in text.splitlines():
  if line.startswith('grp[Compressed X]'):continue
  if line.startswith('allocated '):
   values=line.split()[1:];assert len(values)==10;values[7]='0';line='allocated '+' '.join(values)
  elif line.startswith('groupnr['):
   left,right=line.split('=',1);values=right.split();assert len(values)==10;values[7]='0';line=left+'= '+' '.join(values)
  lines.append(line)
 return '\n'.join(lines)+'\n'

def fields(text):
 out={}
 for line in text.splitlines():
  if '=' in line:k,v=line.split('=',1);out[k.strip()]=v.strip()
  elif line.startswith(('tau-t:','ref-t:')):k,v=line.split(':',1);out[k]=v.strip()
 return out

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--audit-root',type=Path,required=True);p.add_argument('--work-root',type=Path,required=True);a=p.parse_args();w=a.work_root.resolve();w.mkdir(parents=True,exist_ok=False);raw=json.loads((a.audit_root/'audit.json').read_text());assert len(raw['records'])==960 and not raw['failures']
 maps={}
 for h in raw['groups']['topology']:
  text=forcefield_without_output_group((a.audit_root/f'topology_{h}.txt').read_text());digest=hashlib.sha256(text.encode()).hexdigest();maps[h]=digest
  if not (w/f'forcefield_{digest}.txt').exists():(w/f'forcefield_{digest}.txt').write_text(text)
 assert len(set(maps.values()))==1,'Topology changes beyond compressed-output group'
 protocol={h:fields((a.audit_root/f'protocol_{h}.txt').read_text()) for h in raw['groups']['protocol']}
 keys=['cutoff-scheme','coulombtype','coulomb-modifier','rcoulomb-switch','rcoulomb','epsilon-r','epsilon-rf','vdwtype','vdw-modifier','rvdw-switch','rvdw','DispCorr','ewald-rtol','ewald-rtol-lj','lj-pme-comb-rule']
 familykeys=['pme-order','fourierspacing','tcoupl','nsttcouple','tau-t','pcoupl','tau-p']
 families={};nonbonded={};summary=[]
 for r in raw['records']:
  par=protocol[r['protocol_sha256']];family={k:par.get(k) for k in familykeys};label=json.dumps(family,sort_keys=True);nb=json.dumps({k:par.get(k) for k in keys},sort_keys=True)
  families.setdefault(label,[]).append([r['slot'],r['segment']]);nonbonded.setdefault(nb,[]).append([r['slot'],r['segment']])
  summary.append({'slot':r['slot'],'segment':r['segment'],'topology_sha256':r['topology_sha256'],'normalized_forcefield_sha256':maps[r['topology_sha256']],'protocol_sha256':r['protocol_sha256'],'family':family,'fourier_grid':[par.get(k) for k in ['fourier-nx','fourier-ny','fourier-nz']],'rlist':par.get('rlist'),'nstlist':par.get('nstlist')})
 result={'status':'COMMON_FORCEFIELD_WITH_DIFFERENT_NUMERICAL_AND_DYNAMICS_PROTOCOLS','source_audit_sha256':sha(a.audit_root/'audit.json'),'normalized_forcefield_sha256':next(iter(maps.values())),'normalization':'Only compressed-X output group declaration and its column in allocated/groupnr removed; all force-field and other group fields retained.','raw_topology_groups':len(maps),'nonbonded_parameter_groups':[{'parameters':json.loads(k),'members':v} for k,v in nonbonded.items()],'protocol_families':[{'parameters':json.loads(k),'members':v} for k,v in families.items()],'records':summary,'decision':'Do not infer identical numerical Hamiltonian or kinetics from common force-field topology. Stratify protocol regimes; cross-evaluate numerical energy/force differences and test conformational equilibration before pooling.','script_sha256':sha(__file__)}
 (w/'classification.json').write_text(json.dumps(result,indent=2)+'\n')
 print(json.dumps({'status':result['status'],'nonbonded_groups':len(nonbonded),'families':[{'parameters':json.loads(k),'segments':sorted(set(s for _,s in v)),'slots':len(set(k for k,_ in v)),'members':len(v)} for k,v in families.items()]},indent=2))

if __name__=='__main__':main()
