"""Compare four nonbonded energy components and all atom forces independently."""
import argparse,json,subprocess,hashlib
from pathlib import Path
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--reference',type=Path,required=True)
p.add_argument('--driver',type=Path,required=True)
p.add_argument('--output',type=Path,required=True)
p.add_argument('--match-single-precision-input',action='store_true')
a=p.parse_args()
ref=json.loads(a.reference.read_text());results=[]
for case in ref['cases']:
    top=a.reference.parent/case['name']/'processed.top'
    assert hashlib.sha256(top.read_bytes()).hexdigest()==case['topology_sha256']
    command=[str(a.driver.resolve()),str(top.resolve())]
    if a.match_single_precision_input:command.append('float-coordinates')
    run=subprocess.run(command,capture_output=True,text=True,check=True)
    actual=list(map(float,run.stdout.split()))
    values=dict(zip(case['energy_columns'],case['energy_row']))
    expected=[values[k] for k in ['Coulomb (SR)','LJ (SR)','Coulomb-14','LJ-14']]+case['force_row'][1:]
    assert len(actual)==len(expected)==19
    defects=[abs(x-y)/max(1,abs(y)) for x,y in zip(actual,expected)]
    passed=max(defects)<2e-5
    results.append({'case':case['name'],'status':'PASS' if passed else 'FAIL','worst_scaled_error':max(defects),
                    'titan':actual,'gromacs':expected})
    print(f'{results[-1]["status"]} {case["name"]}: worst error={max(defects):.4g}',flush=True)
record={'status':'PASS' if all(r['status']=='PASS' for r in results) else 'FAIL',
        'tolerance':2e-5,'comparison':'each SR/14 energy component and all forces; GROMACS single precision and text rounding',
        'single_precision_input_coordinates':a.match_single_precision_input,
        'driver_sha256':hashlib.sha256(a.driver.read_bytes()).hexdigest(),'cases':results}
a.output.write_text(json.dumps(record,indent=2)+'\n')
if record['status']!='PASS':raise SystemExit(1)
