"""Compare the native MM+GB endpoint kernel to frozen AmberTools references."""
import argparse,hashlib,json,subprocess
from pathlib import Path
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--reference',type=Path,required=True);p.add_argument('--driver',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
a=p.parse_args();reference=json.loads(a.reference.read_text());rows=[]
for case in reference['references']:
 fixture=a.reference.parent/('igb'+str(case['igb']))/'fixture.txt'
 rgbmax=reference.get('rgbmax_angstrom',999)
 command=[str(a.driver.resolve()),str(fixture.resolve())]
 # A finite manifest rgbmax must reach TITAN too, otherwise the two sides
 # integrate the descreening sum over different ranges.
 if rgbmax is not None and float(rgbmax)<900: command.append(repr(float(rgbmax)))
 r=subprocess.run(command,capture_output=True,text=True,check=True)
 values=list(map(float,r.stdout.split()));assert len(values)==3
 terms=[]
 for name,value in zip(['elec','vdw','gb'],values):
  expected=case['delta_kj'][name];error=abs(value-expected)
  tolerance=reference['fixed_tolerance_kj']+reference['fixed_relative_tolerance']*abs(expected)
  terms.append(dict(component=name,titan_kj=value,amber_kj=expected,absolute_error_kj=error,tolerance_kj=tolerance,status='PASS' if error<=tolerance else 'FAIL'))
 rows.append(dict(igb=case['igb'],terms=terms,fixture_sha256=hashlib.sha256(fixture.read_bytes()).hexdigest()))
status='PASS' if all(t['status']=='PASS' for r in rows for t in r['terms']) else 'FAIL'
record=dict(status=status,scope=reference['scope'],rgbmax_angstrom=reference.get('rgbmax_angstrom'),rows=rows,reference_sha256=hashlib.sha256(a.reference.read_bytes()).hexdigest(),driver_sha256=hashlib.sha256(a.driver.read_bytes()).hexdigest())
with a.output.open('x',encoding='utf-8') as f:json.dump(record,f,indent=2);f.write('\n')
print(json.dumps(record,indent=2));raise SystemExit(status!='PASS')
