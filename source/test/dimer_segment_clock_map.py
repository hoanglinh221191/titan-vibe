"""Verify segment clock resets in every replica without changing source trajectories."""
import argparse,json,hashlib
from pathlib import Path
from xtc_segment_inventory import headers
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--dataset',type=Path,required=True);p.add_argument('--energy-audit',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
energy=json.loads(a.energy_audit.read_text())['segments'];offsets={r['segment']:float(next(iter(r['energy_offsets_ps']))) for r in energy if len(r['energy_offsets_ps'])==1}
replicas=[]
for replica in range(1,49):
 path=a.dataset/f'a{replica}/cat.xtc';blocks=[];previous=None
 for index,r in enumerate(headers(path)):
  local=r['step']*.002;offset=r['time_ps']-local
  if previous is None or r['step']<=previous['step']:
   blocks.append(dict(segment=len(blocks)+1,first_frame_zero_based=index,last_frame_zero_based=index,first_local_ps=local,last_local_ps=local,xtc_offset_ps=offset,atom_count=r['natoms']))
  block=blocks[-1];assert abs(offset-block['xtc_offset_ps'])<1e-5 and r['natoms']==1254
  block['last_frame_zero_based']=index;block['last_local_ps']=local;previous=r
 assert len(blocks)==20
 for block in blocks:
  block['energy_offset_ps_a1_verified']=offsets[block['segment']]
  block['energy_minus_xtc_clock_ps']=offsets[block['segment']]-block['xtc_offset_ps']
 replicas.append({'replica':replica,'file':str(path),'blocks':blocks});print('PASS',replica,'20 segment blocks',flush=True)
assert all(r['blocks']==replicas[0]['blocks'] for r in replicas),'Replica segment clocks differ'
record={'status':'PASS_CLOCK_STRUCTURE','scope':__doc__,'dt_ps':.002,'replicas':replicas,'energy_value_identity_scope':'a1 EDR versus concatenated XVG only; other slots have identical coordinate clock structures but source-energy values were not compared individually.','physical_kinetics_certified':False,'energy_audit_sha256':hashlib.sha256(a.energy_audit.read_bytes()).hexdigest()}
with a.output.open('x',encoding='utf-8') as f:json.dump(record,f,indent=2);f.write('\n')
