"""Paper profile lattice geometry, existing fractal regressions and pipeline."""
import argparse,hashlib,json,os,subprocess
from pathlib import Path
from science_build_profile import engine_arch_flags
root=Path(__file__).resolve().parents[1]
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--compiler',default='clang++' if os.name=='nt' else 'g++-13');p.add_argument('--objects-root',type=Path,default=root/'build/native/src');p.add_argument('--with-engine',action='store_true');p.add_argument('--sanitize',action='store_true');p.add_argument('--work-root',type=Path,required=True);a=p.parse_args();w=a.work_root.resolve();w.mkdir(parents=True,exist_ok=False)
def run(tag,cmd):
 r=subprocess.run(cmd,cwd=root if 'build' in tag else w,capture_output=True,text=True,encoding='utf-8',errors='replace');(w/(tag+'.log')).write_text(r.stdout+r.stderr,encoding='utf-8');r.check_returncode();print(tag,'PASS',flush=True)
flags=['-std=c++23','-O1' if a.sanitize else '-O2','-Isrc','-Isrc/eigen']
if a.sanitize:flags+=['-fsanitize=address,undefined','-fno-omit-frame-pointer']
suffix='.exe' if os.name=='nt' else '';exe=w/('geometry'+suffix)
run('geometry_build',[a.compiler]+flags+['test/titan_llps_paper_geometry_oracle.cpp','src/pbctopo/lattice_math.cpp','-o',str(exe)]);run('geometry_run',[str(exe)])
if not a.sanitize:
 for source in ['llps_paper_fractal_smoke','llps_paper_fractal_integration_smoke']:
  exe=w/(source+suffix);cmd=[a.compiler]+flags+['-DTITAN_LLPS_TESTING','-ffunction-sections','-fdata-sections',f'test/{source}.cpp','src/titan_llps_paper.cpp','src/pbctopo/lattice_math.cpp','-Wl,--gc-sections','-o',str(exe)]
  if os.name=='nt':cmd.insert(1,'-fuse-ld=lld')
  run(source+'_build',cmd);run(source+'_run',[str(exe)])
objects=[]
if a.with_engine:
 if a.sanitize:raise ValueError('Engine sanitizer needs separately instrumented objects')
 objects=sorted(x.resolve() for x in a.objects_root.rglob('*.o') if x.name.lower()!='main.o');rsp=w/'objects.rsp';rsp.write_text('\n'.join('"'+x.as_posix()+'"' for x in objects));exe=w/('pipeline'+suffix)
 cmd=[a.compiler]+flags+['-DTITAN_TEST_ENGINE','-fno-access-control','-Isrc/gmxtraj/include',*engine_arch_flags(a.objects_root),'-fopenmp','-pthread','test/titan_llps_paper_geometry_oracle.cpp','@'+str(rsp),'-o',str(exe)]
 if os.name=='nt':cmd.insert(1,'-fuse-ld=lld')
 run('pipeline_build',cmd);run('pipeline_run',[str(exe)])
sources=['src/titan_llps_profile_geometry.h','src/titan_llps_paper.cpp','src/titan_llps_cluster.cpp','src/titan_llps_cluster.h','src/impl_engine.cpp','test/titan_llps_paper_geometry_oracle.cpp','test/llps_paper_geometry_oracles.py']
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
(w/'verification.json').write_text(json.dumps({'status':'PASS','sanitizer':a.sanitize,'with_engine':a.with_engine,'source_sha256':{p:sha(root/p) for p in sources},'object_sha256':{str(p.relative_to(root)):sha(p) for p in objects}},indent=2)+'\n')
