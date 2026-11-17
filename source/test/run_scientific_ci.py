"""Run bounded scientific and performance gates with per-stage CI evidence."""
import argparse,json,subprocess,sys,os
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--compiler',required=True);p.add_argument('--binary',required=True);a=p.parse_args();w=Path('test/tmp/scientific_ci');w.mkdir(parents=True,exist_ok=False)
os.environ['PYTHONUTF8']='1'
os.environ['PYTHONIOENCODING']='utf-8'
cases=[('cli','tools/test_cli_validation.py',['--binary',a.binary]),('pbe_p0','test/pbe_p0_oracles.py',[]),('science_p1','test/science_p1_oracles.py',['--with-engine','--verify-regression']),('ene02','test/ene02_oracles.py',['--with-engine']),('dham','test/dham_segment_oracles.py',['--with-engine']),('paper','test/llps_paper_geometry_oracles.py',['--with-engine']),('llps_cli','test/llps_phase_cli_smoke.py',['--binary',a.binary]),('bie_route','test/pbe_bie_resolution_ses_smoke.py',['--binary',a.binary]),('gist_geometry','test/gist_geometry_oracles.py',[]),('gist_thermo','test/gist_thermo_oracles.py',[]),('alpb','test/alpb_zero_site_oracles.py',[]),('gist_performance','test/gist_performance_gate.py',[])]
cases.insert(next(i for i,c in enumerate(cases) if c[0]=='alpb'),('gist_cli','test/gist_cli_reference_smoke.py',['--binary',a.binary,'--driver','test/tmp/gist_thermo_ci/engine'+('.exe' if a.binary.endswith('.exe') else '')]))
records=[]
cases.insert(0,('abi_profile','test/science_build_profile_oracles.py',[]))
for name,script,options in cases:
 cmd=[sys.executable,script]+options
 if name!='cli':cmd+=['--work-root',str(Path('test/tmp')/(name+'_ci'))]
 if name not in ['cli','llps_cli','gist_cli','bie_route']:cmd+=['--compiler',a.compiler]
 r=subprocess.run(cmd,capture_output=True,text=True,encoding='utf-8',errors='replace');(w/(name+'.log')).write_text(r.stdout+r.stderr,encoding='utf-8');records.append({'case':name,'status':'PASS' if r.returncode==0 else 'FAIL','exit_code':r.returncode});(w/'summary.json').write_text(json.dumps({'status':'IN_PROGRESS' if r.returncode==0 else 'FAIL','cases':records},indent=2)+'\n');print(name,records[-1]['status'],flush=True)
 if r.returncode:print((r.stdout+r.stderr)[-6000:]);raise SystemExit(r.returncode)
(w/'summary.json').write_text(json.dumps({'status':'PASS','cases':records},indent=2)+'\n')
