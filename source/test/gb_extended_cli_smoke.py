"""Exercise two-frame GBn2/GBNSR6 through the complete TITAN executable.

Uses the fixture produced by gbsa_perf_entropy_smoke.ps1; --fixture-root points
to that run. Each run uses a fresh directory and preserves failures.
"""
import argparse, hashlib, json, shutil, subprocess
from pathlib import Path
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--binary',type=Path,required=True)
p.add_argument('--fixture-root',type=Path,required=True)
p.add_argument('--work-root',type=Path,required=True)
a=p.parse_args();w=a.work_root.resolve();w.mkdir(parents=True,exist_ok=False)
base=(a.fixture_root/'arg.in').read_text().splitlines();rows=[]
for model in [8,66]:
 for nt in [1,2]:
  folder=w/f'igb{model}_nt{nt}';folder.mkdir()
  opts=[l for l in base if not l.startswith('-nt=')]
  opts += [f'-nt={nt}',f'-igb={model}'] if model==8 else [f'-nt={nt}','-gbnsr6','-gbnsr6-spacing=0.25']
  inp=folder/'arg.in';inp.write_text('\n'.join(opts)+'\n')
  result=subprocess.run([str(a.binary.resolve()),'-a','arg.in'],cwd=folder,capture_output=True,text=True,encoding='utf-8',errors='replace')
  (folder/'console.log').write_text(result.stdout+result.stderr,encoding='utf-8');result.check_returncode()
  summary=(folder/'gbsa_summary.dat').read_text()
  assert f'gbsa_igb {model}' in summary
  snapshots=(folder/'traj_0_gbsa_snapshots.dat').read_text().splitlines()[1:]
  assert len(snapshots)==2
  if model==66:
   assert 'gbsa_model_status experimental_SES_quadrature' in summary
   assert 'alpb_coulomb_scaling none_GBNSR6' in summary
   assert 'gbnsr6_spacing_angstrom 0.25' in summary
   assert (folder/'gbsa_alpb_sizes.csv').exists()
  rows.append(dict(model=model,threads=nt,snapshots=snapshots,status='PASS'))
  print(model,nt,'PASS',flush=True)
 for one,two in zip(rows[-2]['snapshots'],rows[-1]['snapshots']):assert one==two,'thread-dependent energy'

# Config errors must be rejected, rather than silently ignoring incompatible options.
for label,extra in [('cutoff',['-gbsa-rgbmax=25']),('spacing',['-gbnsr6-spacing=0']),('cas',['-cas'])]:
 folder=w/('reject_'+label);folder.mkdir();inp=folder/'arg.in'
 inp.write_text('\n'.join(base+['-gbnsr6']+extra)+'\n')
 result=subprocess.run([str(a.binary.resolve()),'-a','arg.in'],cwd=folder,capture_output=True,text=True,encoding='utf-8',errors='replace')
 text=result.stdout+result.stderr;(folder/'console.log').write_text(text,encoding='utf-8')
 assert result.returncode!=0 and 'GBNSR6 requires' in text,(label,text)
report={'status':'PASS','records':rows,'rejections':['cutoff','spacing','cas'],
 'binary_sha256':hashlib.sha256(a.binary.read_bytes()).hexdigest()}
(w/'verification.json').write_text(json.dumps(report,indent=2)+'\n')
