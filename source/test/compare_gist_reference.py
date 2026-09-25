"""Compare production GIST with independent arithmetic and CPPTRAJ."""
import argparse,csv,hashlib,json,math,subprocess
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--reference',type=Path,required=True);p.add_argument('--driver',type=Path,required=True);p.add_argument('--work-root',type=Path,required=True);a=p.parse_args();w=a.work_root.resolve();w.mkdir(parents=True,exist_ok=False)
ref=json.loads((a.reference/'reference.json').read_text());assert ref['status']=='REFERENCE_COMPLETE'
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
for name,digest in ref['inputs_sha256'].items():assert sha(a.reference/name)==digest
prefix=w/'titan';r=subprocess.run([str(a.driver.resolve()),str((a.reference/'titan_fixture.txt').resolve()),str(prefix)],cwd=w,capture_output=True,text=True,encoding='utf-8',errors='replace');(w/'driver.log').write_text(r.stdout+r.stderr,encoding='utf-8');r.check_returncode()
titan=list(csv.DictReader((w/'titan_numerical.csv').open()));assert len(titan)==4
lines=(a.reference/'amber-output.dat').read_text().splitlines();header=next(line.split() for line in lines if line.lstrip().startswith('voxel '));amber={}
for line in lines:
 words=line.split()
 if len(words)!=len(header) or not words[0].isdigit():continue
 row=dict(zip(header,map(float,words)))
 if row['population']>0:amber[int(row['voxel'])]=row
assert len(amber)==4
checks=[]
for row in titan:
 xyz=[float(row[k]) for k in ['x','y','z']];expected=min(ref['numerical_reference'],key=lambda r:sum((a-b)**2 for a,b in zip(xyz,r['center_A'])));ar=amber[int(row['voxel'])];assert int(row['population'])==ref['frames']==ar['population']
 # CPPTRAJ stores T delta S, whereas TITAN explicitly exports minus T delta S.
 pairs=[('Eww_arithmetic',float(row['Eww_kj']),expected['Eww_kj'],2e-8),('Sorient_arithmetic',float(row['minus_TSorient_kj']),expected['minus_TSorient_kj'],2e-7),('Eww_CPPTRAJ',float(row['Eww_kj']),ar['Eww-norm-unref(kcal/mol)']*4.184,.002),('Sorient_CPPTRAJ',float(row['minus_TSorient_kj']),-ar['dTSorient-norm(kcal/mol)']*4.184,.002)]
 for name,actual,reference,tol in pairs:
  error=abs(actual-reference);checks.append({'voxel':int(row['voxel']),'term':name,'titan':actual,'reference':reference,'error':error,'tolerance':tol,'status':'PASS' if math.isfinite(actual) and error<=tol else 'FAIL'})
status='PASS' if all(r['status']=='PASS' for r in checks) else 'FAIL';result={'status':status,'scope':'Synthetic four-water/ghost-solute fixture, 1024 orientations per site, all pairs inside cutoff; no physical equilibration or bulk hydration free energy claim','checks':checks,'reference_sha256':sha(a.reference/'reference.json'),'cpptraj_output_sha256':sha(a.reference/'amber-output.dat'),'driver_sha256':sha(a.driver),'script_sha256':sha(__file__)}
(w/'comparison.json').write_text(json.dumps(result,indent=2)+'\n');print(status,[(r['term'],r['error']) for r in checks]);raise SystemExit(status!='PASS')
