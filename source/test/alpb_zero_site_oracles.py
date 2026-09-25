"""Zero-radius ALPB volume sites and production failure propagation."""
import argparse,hashlib,json,os,subprocess
from pathlib import Path
from science_build_profile import engine_arch_flags
root=Path(__file__).resolve().parents[1]
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--compiler',default='g++-13');p.add_argument('--objects-root',type=Path,default=root/'build/native/src');p.add_argument('--work-root',type=Path,required=True);a=p.parse_args();w=a.work_root.resolve();w.mkdir(parents=True,exist_ok=False)
objects=sorted(p.resolve() for p in a.objects_root.rglob('*.o') if p.name.lower()!='main.o');rsp=w/'objects.rsp';rsp.write_text('\n'.join('"'+p.as_posix()+'"' for p in objects));suffix='.exe' if os.name=='nt' else ''
flags=['-std=c++23','-O2','-Isrc','-Isrc/eigen','-Isrc/gmxtraj/include',*engine_arch_flags(a.objects_root),'-fopenmp','-pthread']
if os.name=='nt':flags+=['-fuse-ld=lld']
def run(tag,cmd,cwd,fail=False):
 r=subprocess.run(cmd,cwd=cwd,capture_output=True,text=True,encoding='utf-8',errors='replace');text=r.stdout+r.stderr;(w/(tag+'.log')).write_text(text)
 if fail:assert r.returncode!=0 and 'zero-radius spectator' in text
 else:r.check_returncode()
 print(tag,'EXPECTED_FAIL' if fail else 'PASS',flush=True)
exe=w/('oracle'+suffix);source='test/titan_alpb_zero_site_oracle.cpp';run('build',[a.compiler]+flags+['-DTITAN_TEST_ENGINE',source,'@'+str(rsp),'-o',str(exe)],root);run('run',[str(exe)],w)
mut=w/'mutant';mut.mkdir();text=(root/'src/titan_gb_model.h').read_text();anchor='!std::isfinite(radius) || radius < 0.0)';assert text.count(anchor)==1;(mut/'titan_gb_model.h').write_text(text.replace(anchor,'!std::isfinite(radius) || radius <= 0.0)'))
bad=w/('mutant_oracle'+suffix);run('mutation_build',[a.compiler,'-I'+str(mut)]+flags+[source,'-o',str(bad)],root);run('mutation_run',[str(bad)],w,True)
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
(w/'verification.json').write_text(json.dumps({'status':'PASS','mutation':'EXPECTED_FAIL','source_sha256':{p:sha(root/p) for p in ['src/titan_gb_model.h','src/gbsa_intrin.cpp',source,'test/alpb_zero_site_oracles.py']},'oracle_sha256':sha(exe),'object_sha256':{str(p.relative_to(root)):sha(p) for p in objects}},indent=2)+'\n')
