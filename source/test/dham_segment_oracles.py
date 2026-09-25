"""Native DHAM metadata gates and production postprocessing integration."""
import argparse,hashlib,json,os,subprocess
from pathlib import Path
from science_build_profile import engine_arch_flags
root=Path(__file__).resolve().parents[1]
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--compiler',default='clang++' if os.name=='nt' else 'g++-13');p.add_argument('--objects-root',type=Path,default=root/'build/native/src');p.add_argument('--with-engine',action='store_true');p.add_argument('--sanitize',action='store_true');p.add_argument('--real-metadata',type=Path);p.add_argument('--work-root',type=Path,required=True);a=p.parse_args();w=a.work_root.resolve();w.mkdir(parents=True,exist_ok=False)
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def run(tag,cmd,cwd,expect_failure=False):
 r=subprocess.run(cmd,cwd=cwd,capture_output=True,text=True,encoding='utf-8',errors='replace');(w/(tag+'.log')).write_text(r.stdout+r.stderr,encoding='utf-8')
 if expect_failure:
  if r.returncode==0 or 'transition oracle' not in r.stdout+r.stderr:raise RuntimeError('Mutation did not fail at the intended transition check')
 else:r.check_returncode()
 print(tag, 'EXPECTED_FAIL' if expect_failure else 'PASS',flush=True)
suffix='.exe' if os.name=='nt' else ''
flags=['-std=c++23','-O1' if a.sanitize else '-O2','-Isrc']
if a.sanitize:flags+=['-fsanitize=address,undefined','-fno-omit-frame-pointer']
exe=w/('oracle'+suffix);run('oracle_build',[a.compiler]+flags+['test/titan_dham_segments_oracle.cpp','-o',str(exe)],root)
case=w/'unit';case.mkdir();run('unit',[str(exe)],case)
mut=w/'mutant';mut.mkdir();s=(root/'src/titan_dham_segments.h').read_text();anchor='Decision transition(int replica,const Frame& first,const Frame& last,double lag_ps) const {';assert s.count(anchor)==1
(mut/'titan_dham_segments.h').write_text(s.replace(anchor,anchor+'\n  return Decision::eligible; // regression: unfiltered legacy counts\n'))
bad=w/('mutant_oracle'+suffix);run('mutation_build',[a.compiler,'-I'+str(mut)]+flags+['test/titan_dham_segments_oracle.cpp','-o',str(bad)],root);run('mutation',[str(bad)],mut,True)
if a.real_metadata:run('real_metadata',[str(exe),str(a.real_metadata.resolve())],w)
objects=[]
if a.with_engine:
 if a.sanitize:raise ValueError('Engine objects need separate instrumentation; use header-only sanitizer gate')
 objects=sorted(x.resolve() for x in a.objects_root.rglob('*.o') if x.name.lower()!='main.o');assert objects
 rsp=w/'objects.rsp';rsp.write_text('\n'.join('"'+x.as_posix()+'"' for x in objects))
 linkflags=['-Isrc/eigen','-Isrc/gmxtraj/include',*engine_arch_flags(a.objects_root),'-fopenmp','-pthread']
 if os.name=='nt':linkflags+=['-fuse-ld=lld']
 pipeline=w/('pipeline'+suffix);run('pipeline_build',[a.compiler]+flags+linkflags+['test/titan_dham_segments_pipeline.cpp','@'+str(rsp),'-o',str(pipeline)],root)
 case=w/'pipeline_run';case.mkdir();run('pipeline',[str(pipeline)],case)
sources=['src/titan_dham_segments.h','src/engine_postproc.cpp','src/arguments.h','src/arguments.cpp','src/cli/parse_arguments.cpp','test/titan_dham_segments_oracle.cpp','test/titan_dham_segments_pipeline.cpp',str(Path(__file__).relative_to(root))]
record={'status':'PASS','mutation':'EXPECTED_FAIL','with_engine':a.with_engine,'sanitizer':a.sanitize,'source_sha256':{p:sha(root/p) for p in sources},'object_sha256':{str(p.relative_to(root)):sha(p) for p in objects},'oracle_sha256':sha(exe)}
if a.real_metadata:record['real_metadata_sha256']=sha(a.real_metadata)
(w/'verification.json').write_text(json.dumps(record,indent=2)+'\n')
