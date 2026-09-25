"""Energy diagnostics stratified by TPR-defined regimes, not energy-selected cuts."""
import argparse,collections,gzip,json
from pathlib import Path
import numpy as np
from dimer_segment_alignment import sha
from dimer_equilibrium_diagnostics import overlap
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--classification',type=Path,required=True);p.add_argument('--samples',type=Path,required=True);p.add_argument('--work-root',type=Path,required=True);a=p.parse_args();w=a.work_root.resolve();w.mkdir(parents=True,exist_ok=False)
c=json.loads(a.classification.read_text());families=[];lookup={}
for i,f in enumerate(c['protocol_families']):
 segs=sorted(set(s for k,s in f['members']));families.append({'id':i,'segments':segs,'parameters':f['parameters']})
 for key in f['members']:assert tuple(key) not in lookup;lookup[tuple(key)]=i
data=collections.defaultdict(list)
with gzip.open(a.samples,'rt') as f:
 for line in f:
  r=json.loads(line);data[lookup[r['slot'],r['segment']],r['slot']].append(r['configurational_enthalpy_kj_mol'])
per_slot=[];adjacent=[]
for family in families:
 i=family['id']
 for k in range(1,49):
  v=np.array(data[i,k]);n=len(v);sd=v.std(ddof=1)
  per_slot.append({'regime':i,'slot':k,'samples':n,'mean_H_kj_mol':float(v.mean()),'sd_H_kj_mol':float(sd),'half_split_shift_sd':float((v[n//2:].mean()-v[:n//2].mean())/sd)})
  if k<48:
   neighbor=np.array(data[i,k+1]);adjacent.append({'regime':i,'slots':[k,k+1],'overlap_64':overlap(v,neighbor,64),'overlap_128':overlap(v,neighbor,128)})
summary=[{'regime':r['id'],'segments':r['segments'],'samples_per_slot':next(v['samples'] for v in per_slot if v['regime']==r['id']),'min_adjacent_overlap_64':min(v['overlap_64'] for v in adjacent if v['regime']==r['id']),'min_adjacent_overlap_128':min(v['overlap_128'] for v in adjacent if v['regime']==r['id']),'max_abs_half_split_shift_sd':max(abs(v['half_split_shift_sd']) for v in per_slot if v['regime']==r['id'])} for r in families]
result={'status':'STRATIFIED_ENERGY_DIAGNOSTICS_NOT_EQUILIBRATION_CERTIFICATE','regimes':families,'summary':summary,'per_slot':per_slot,'adjacent':adjacent,'classification_sha256':sha(a.classification),'samples_sha256':sha(a.samples),'script_sha256':sha(__file__),'decision':'Energy overlap is present within each TPR-defined regime; conformational overlap/equilibration, joint effective sample size and numerical equivalence across regimes remain unverified. No fit or burn-in selection.'}
(w/'diagnostics.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(summary,indent=2))
