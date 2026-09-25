"""GIST registration/output gates with targeted rotation and frozen-energy mutations."""
import argparse,hashlib,json,os,subprocess
from pathlib import Path
from science_build_profile import engine_arch_flags
root=Path(__file__).resolve().parents[1]
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--compiler',default='clang++' if os.name=='nt' else 'g++-13');p.add_argument('--objects-root',type=Path,default=root/'build/native/src');p.add_argument('--sanitize-gist',action='store_true');p.add_argument('--work-root',type=Path,required=True);a=p.parse_args();w=a.work_root.resolve();w.mkdir(parents=True,exist_ok=False)
def run(tag,cmd,cwd,expected_error=None):
 r=subprocess.run(cmd,cwd=cwd,capture_output=True,text=True,encoding='utf-8',errors='replace');text=r.stdout+r.stderr;(w/(tag+'.log')).write_text(text,encoding='utf-8')
 if expected_error:
  if r.returncode==0 or expected_error not in text:raise RuntimeError('Mutation did not fail at expected oracle: '+text[-1000:])
 else:r.check_returncode()
 print(tag,'EXPECTED_FAIL' if expected_error else 'PASS',flush=True)
objects=sorted(p.resolve() for p in a.objects_root.rglob('*.o') if p.name.lower() not in ['main.o','titan_gist.o']);assert objects
rsp=w/'objects.rsp';rsp.write_text('\n'.join('"'+p.as_posix()+'"' for p in objects))
flags=['-std=c++23','-O1' if a.sanitize_gist else '-O2','-fno-access-control','-Isrc','-Isrc/eigen','-Isrc/gmxtraj/include',*engine_arch_flags(a.objects_root),'-fopenmp','-pthread']
if os.name=='nt':flags+=['-fuse-ld=lld']
if a.sanitize_gist:flags+=['-fsanitize=address,undefined','-fno-omit-frame-pointer']
suffix='.exe' if os.name=='nt' else '';source='test/titan_gist_geometry_oracle.cpp';exe=w/('oracle'+suffix)
run('build',[a.compiler]+flags+[source,'@'+str(rsp),'-o',str(exe)],root)
case=w/'checked';case.mkdir();run('run',[str(exe)],case)
j=json.loads((case/'gist_checked_observables.json').read_text());assert j['processed_frames']==2 and not j['free_energy_available'] and j['E_ww_available'] and j['orientational_entropy_available']
mutations={}
if not a.sanitize_gist:
 text=(root/'src/titan_gist.cpp').read_text()
 for label,anchor,replacement,error in [
  ('rotation','return U * D * V.transpose();','return V * D * U.transpose();','rigid image/rotation occupancy invariance'),
  ('frozen_receptor','for(size_t j=0;j<n_rec_heavy_;++j) {','for(size_t j=n_rec_heavy_;j<n_rec_heavy_;++j) {','current receptor energy')]:
  assert text.count(anchor)==1
  mutant=w/(label+'.cpp');mutant.write_text(text.replace(anchor,replacement));bad=w/(label+'_oracle'+suffix)
  run(label+'_build',[a.compiler]+flags+['-DTITAN_GIST_IMPLEMENTATION="'+mutant.as_posix()+'"',source,'@'+str(rsp),'-o',str(bad)],root)
  case=w/label;case.mkdir();run(label+'_run',[str(bad)],case,error);mutations[label]='EXPECTED_FAIL'
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
sources=['src/titan_gist.cpp','src/titan_gist.h','test/titan_gist_geometry_oracle.cpp','test/gist_geometry_oracles.py']
(w/'verification.json').write_text(json.dumps({'status':'PASS','sanitizer_gist_and_test_only':a.sanitize_gist,'mutations':mutations,'source_sha256':{p:sha(root/p) for p in sources},'object_sha256':{str(p.relative_to(root)):sha(p) for p in objects},'oracle_sha256':sha(exe)},indent=2)+'\n')
