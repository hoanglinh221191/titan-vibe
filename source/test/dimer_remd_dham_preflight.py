"""Read-only qualification of stored dimer REMD metadata for DHAM validation.

Does not generate physical transition counts or fit a state model. A temperature
slot trajectory can jump at exchanges; demultiplexing alone changes temperature
and is insufficient to establish target-temperature physical kinetics.
"""
import argparse,hashlib,json,math,re
from pathlib import Path
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--dataset',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
a=p.parse_args();rows=[]
for folder in sorted((x for x in a.dataset.iterdir() if x.is_dir() and re.fullmatch(r'a\d+',x.name)),key=lambda x:int(x.name[1:])):
 mdps=list(folder.glob('md_*.mdp'));mdp=next((x for x in mdps if x.stem=='md_'+folder.name[1:]),mdps[0] if mdps else None)
 params={}
 if mdp:
  for line in mdp.read_text(errors='replace').splitlines():
   line=line.split(';')[0]
   if '=' in line:k,v=line.split('=',1);params[k.strip().lower().replace('_','-')]=v.strip()
 row={'replica':folder.name,'path':str(folder),'trajectory_present':(folder/'cat.xtc').exists(),'params':{k:params.get(k) for k in ['dt','nstenergy','nstxout-compressed','ref-t','ref-p','pcoupl','tcoupl']}}
 if mdp:row['mdp_sha256']=hashlib.sha256(mdp.read_bytes()).hexdigest()
 energy=folder/('energy_'+folder.name[1:]+'.xvg');row['energy_file']=str(energy);row['energy_present']=energy.exists()
 if energy.exists():
  labels={};count=0;first=last=None;steps=set();bad=0;duplicate=0;nonmonotonic=0
  for line in energy.read_text().splitlines():
   m=re.search(r'@\s+s(\d+)\s+legend\s+"(.*?)"',line)
   if m:labels[int(m[1])]=m[2]
   if not line.strip() or line.startswith(('#','@')):continue
   try:values=list(map(float,line.split()))
   except ValueError:bad+=1;continue
   if len(values)!=4 or not all(map(math.isfinite,values)):bad+=1;continue
   t=values[0]
   if last is not None:
    steps.add(round(t-last,6));duplicate+=int(t==last);nonmonotonic+=int(t<last)
   if first is None:first=t
   last=t;count+=1
  row['energy']={'labels':labels,'frames':count,'first_ps':first,'last_ps':last,'intervals_ps':sorted(steps),'invalid_rows':bad,'duplicate_adjacent_times':duplicate,'time_reversals':nonmonotonic,'sha256':hashlib.sha256(energy.read_bytes()).hexdigest()}
 rows.append(row)
log=a.dataset/'a1/10.log';text=log.read_text(errors='replace') if log.exists() else ''
match=re.search(r'-replex\s+(\d+)',text);exchange_steps=int(match[1]) if match else None
result={'scope':__doc__,'dataset':str(a.dataset),'replica_count':len(rows),'reference_log_sha256':hashlib.sha256(log.read_bytes()).hexdigest() if log.exists() else None,'exchange_steps':exchange_steps,'constant_pressure_remd':'Using Constant Pressure REMD' in text,'replicas':rows,
 'decision':'QUALIFIED_METADATA_ONLY_NO_PHYSICAL_DHAM_COUNTS',
 'requirements':['Verify concatenation segment time identity and exchange history against coordinate frames.',
 'Use U+pV for NPT reduced potentials with matched frame times; do not drop the volume term.',
 'Coordinate cadence is coarser than exchange cadence; require exchange-free constant-temperature lag intervals. Do not count slot jumps or infer physical rates from demultiplexed changing-temperature paths.',
 'State-energy DHAM remains a small-bias approximation. Use a separate equilibrium reweighting validation for equilibrium populations.']}
with a.output.open('x',encoding='utf-8') as f:json.dump(result,f,indent=2);f.write('\n')
print('QUALIFIED metadata',len(rows),'replicas; physical DHAM transitions NOT QUALIFIED')
