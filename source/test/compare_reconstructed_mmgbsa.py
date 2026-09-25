"""Component comparison and SASA resolution sweep against fixed AmberTools references."""
import argparse,json,hashlib,math,subprocess
from pathlib import Path
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--reference',type=Path,required=True);p.add_argument('--driver',type=Path,required=True);p.add_argument('--work-root',type=Path,required=True);a=p.parse_args();w=a.work_root.resolve();w.mkdir(parents=True,exist_ok=False)
ref=json.loads(a.reference.read_text());records=[]
for model in ref['models']:
 refs=[r for r in ref['records'] if r['model']==model]
 fixture=a.reference.parent/refs[0]['fixture']
 for points in [384,1536,6144]:
  run=subprocess.run([str(a.driver.resolve()),str(fixture.resolve()),str(points)],capture_output=True,text=True,check=True)
  (w/f'model_{model}_points_{points}.log').write_text(run.stdout+run.stderr)
  rows=[list(map(float,line.split())) for line in run.stdout.splitlines() if len(line.split())==9]
  assert len(rows)==2
  for values,r in zip(rows,refs):
   salt,ee,ev,gb,np,total,ar,al,ac=values;assert abs(salt-r['salt'])<1e-10
   actual=dict(elec=ee,vdw=ev,gb=gb,nonpolar=np,total_mmgbsa=total,sasa_a2=ac-ar-al)
   terms=[]
   for name in ['elec','vdw','gb','sasa_a2','nonpolar','total_mmgbsa']:
    expected=r['delta'][name];error=abs(actual[name]-expected)
    sasa_tol=max(ref['sasa_delta_tolerance_abs_a2'],ref['sasa_delta_tolerance_rel']*abs(r['delta']['sasa_a2']))
    energy_tol=ref['mmgb_tolerance_abs_kj']+ref['mmgb_tolerance_rel']*abs(expected)
    tol=sasa_tol if name=='sasa_a2' else sasa_tol*ref['gamma_kcal_a2']*4.184 if name=='nonpolar' else energy_tol+sasa_tol*ref['gamma_kcal_a2']*4.184 if name=='total_mmgbsa' else energy_tol
    terms.append(dict(component=name,titan=actual[name],amber=expected,absolute_error=error,tolerance=tol,status='PASS' if math.isfinite(actual[name]) and error<=tol else 'FAIL'))
   record=dict(model=model,points=points,salt=salt,terms=terms,titan_endpoint_sasa_a2=[ar,al,ac],amber_endpoint_sasa_a2=r['endpoint_sasa_a2']);records.append(record)
   print(model,points,salt,'; '.join(t['component']+':'+t['status'] for t in terms),flush=True)
status='PASS' if all(t['status']=='PASS' for r in records if r['points']==6144 for t in r['terms']) else 'FAIL'
result=dict(status=status,scope=ref['scope'],high_resolution_gate=6144,records=records,reference_sha256=hashlib.sha256(a.reference.read_bytes()).hexdigest(),driver_sha256=hashlib.sha256(a.driver.read_bytes()).hexdigest())
(w/'comparison.json').write_text(json.dumps(result,indent=2)+'\n');raise SystemExit(status!='PASS')
