"""LLPS geometry and pipeline gate. Keep each execution in a fresh directory."""
import argparse,hashlib,json,os,subprocess
from pathlib import Path
from science_build_profile import engine_arch_flags
root=Path(__file__).resolve().parents[1]
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--compiler',default='clang++' if os.name=='nt' else 'g++-13')
p.add_argument('--work-root',type=Path,required=True)
p.add_argument('--objects-root',type=Path,default=root/'build/native/src')
p.add_argument('--with-engine',action='store_true');p.add_argument('--sanitize',action='store_true')
a=p.parse_args();work=a.work_root.resolve();work.mkdir(parents=True,exist_ok=False)
def run(tag,cmd):
 r=subprocess.run(cmd,cwd=root if 'build' in tag else work,capture_output=True,text=True,encoding='utf-8',errors='replace')
 (work/(tag+'.log')).write_text(r.stdout+r.stderr,encoding='utf-8')
 if 'run' in tag or r.returncode: print(r.stdout+r.stderr,end='',flush=True)
 r.check_returncode()
flags=['-std=c++23','-O1' if a.sanitize else '-O2','-mavx2','-mfma','-Isrc','-Isrc/eigen']
if a.sanitize:flags+=['-g','-fsanitize=address,undefined','-fno-omit-frame-pointer']
suffix='.exe' if os.name=='nt' else ''
exe=work/('geometry'+suffix)
run('geometry_build',[a.compiler]+flags+['test/titan_llps_geometry_oracle.cpp','src/pbctopo/lattice_math.cpp','-o',str(exe)])
run('geometry_run',[str(exe)])
record={'geometry':'PASS','sanitizers':a.sanitize}
if a.with_engine:
 if a.sanitize:raise ValueError('Use a separately instrumented engine build for pipeline sanitizers')
 objects=sorted(x.resolve() for x in a.objects_root.rglob('*.o') if x.name.lower()!='main.o')
 rsp=work/'objects.rsp';rsp.write_text('\n'.join('"'+x.as_posix()+'"' for x in objects),encoding='utf-8')
 exe=work/('pipeline'+suffix)
 cmd=[a.compiler]+[f for f in flags if f not in ('-mavx2','-mfma')]+engine_arch_flags(a.objects_root)+['-fno-access-control','-pthread','-fopenmp','-Isrc/gmxtraj/include','test/titan_llps_pipeline_oracle.cpp','@'+str(rsp),'-o',str(exe)]
 if os.name=='nt':cmd.insert(1,'-fuse-ld=lld')
 run('pipeline_build',cmd);run('pipeline_run',[str(exe)]);record['pipeline']='PASS'
 record['objects_sha256']={str(x.relative_to(root)):hashlib.sha256(x.read_bytes()).hexdigest() for x in objects}
sources=list((root/'src').glob('titan_llps*.cpp'))+[root/'src/titan_llps.h',root/'src/titan_llps_geometry.h',*sorted((root/'src/cli').glob('parse_*.cpp')),root/'src/arguments.h',*sorted((root/'src/engine').iterdir()),root/'src/pbctopo/lattice_math.cpp']
record['sources_sha256']={str(x.relative_to(root)):hashlib.sha256(x.read_bytes()).hexdigest() for x in sources}
(work/'verification.json').write_text(json.dumps(record,indent=2)+'\n',encoding='utf-8')
