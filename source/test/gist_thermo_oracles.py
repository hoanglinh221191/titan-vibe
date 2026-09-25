"""SO(3) analytic/regression checks and production Eww/entropy gates."""
import argparse,hashlib,json,os,subprocess
from pathlib import Path
from science_build_profile import engine_arch_flags
root=Path(__file__).resolve().parents[1]
p=argparse.ArgumentParser();p.add_argument('--compiler',default='g++-13');p.add_argument('--objects-root',type=Path,default=root/'build/native/src');p.add_argument('--work-root',type=Path,required=True);p.add_argument('--sanitize',action='store_true');a=p.parse_args();w=a.work_root.resolve();w.mkdir(parents=True,exist_ok=False)
flags=['-std=c++23','-O1' if a.sanitize else '-O2','-fno-access-control','-Isrc','-Isrc/eigen','-Isrc/gmxtraj/include',*engine_arch_flags(a.objects_root),'-fopenmp','-pthread']
if os.name=='nt':flags+=['-fuse-ld=lld']
if a.sanitize:flags+=['-fsanitize=address,undefined','-fno-omit-frame-pointer']
suffix='.exe' if os.name=='nt' else '';env=dict(os.environ,ASAN_OPTIONS='detect_leaks=1:halt_on_error=1',UBSAN_OPTIONS='halt_on_error=1:print_stacktrace=1')
def run(name,cmd,cwd,expected=False):
 r=subprocess.run(cmd,cwd=cwd,capture_output=True,text=True,encoding='utf-8',errors='replace',env=env);text=r.stdout+r.stderr;(w/(name+'.log')).write_text(text,encoding='utf-8')
 if expected:assert r.returncode!=0 and 'Haar ball pi normalization' in text,text[-1000:]
 else:r.check_returncode()
 print(name,'EXPECTED_FAIL' if expected else 'PASS',flush=True)
source='test/titan_gist_thermo_oracle.cpp';unit=w/('unit'+suffix);run('unit_build',[a.compiler]+flags+[source,'-o',str(unit)],root);run('unit_run',[str(unit)],w)
objects=sorted(p.resolve() for p in a.objects_root.rglob('*.o') if p.name.lower()!='main.o');assert len(objects)>200
rsp=w/'objects.rsp';rsp.write_text('\n'.join('"'+p.as_posix()+'"' for p in objects));exe=w/('engine'+suffix);run('engine_build',[a.compiler]+flags+['-DTITAN_TEST_ENGINE',source,'@'+str(rsp),'-o',str(exe)],root);run('engine_run',[str(exe)],w)
mut=w/'mutant';mut.mkdir();text=(root/'src/titan_gist_entropy.h').read_text();anchor='return volume/std::acos(-1.);';assert text.count(anchor)==1;(mut/'titan_gist_entropy.h').write_text(text.replace(anchor,'return volume/(2*std::acos(-1.));'))
bad=w/('mutant_unit'+suffix);run('mutation_build',[a.compiler,'-I'+str(mut)]+flags+[source,'-o',str(bad)],root);run('mutation_run',[str(bad)],w,True)
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
(w/'verification.json').write_text(json.dumps({'status':'PASS','sanitized':a.sanitize,'mutation':'EXPECTED_FAIL','source_sha256':{p:sha(root/p) for p in ['src/titan_gist_entropy.h','src/titan_gist.cpp','src/titan_gist.h',source,'test/gist_thermo_oracles.py']},'engine_sha256':sha(exe),'objects_sha256':{str(p.relative_to(root)):sha(p) for p in objects}},indent=2)+'\n')
