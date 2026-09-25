"""Run DBSCAN/LJ independent geometry and executor regression gates."""
import argparse,hashlib,json,os,subprocess
from pathlib import Path
from science_build_profile import engine_arch_flags
root=Path(__file__).resolve().parents[1]
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--compiler',default='clang++' if os.name=='nt' else 'g++-13');p.add_argument('--work-root',type=Path,required=True);p.add_argument('--objects-root',type=Path,default=root/'build/native/src');p.add_argument('--verify-regression',action='store_true');a=p.parse_args();w=a.work_root.resolve();w.mkdir(parents=True,exist_ok=False)
def run(tag,cmd,expected=0):
 r=subprocess.run(cmd,cwd=root if 'build' in tag else w,capture_output=True,text=True,encoding='utf-8',errors='replace');(w/(tag+'.log')).write_text(r.stdout+r.stderr,encoding='utf-8')
 if 'run' in tag or r.returncode!=expected:print(r.stdout+r.stderr,end='',flush=True)
 if expected==0:r.check_returncode()
 elif r.returncode==0 or 'DBSCAN neighbor' not in r.stdout+r.stderr:raise RuntimeError('Worker callback mutation did not trigger intended regression')
objects=sorted(x.resolve() for x in a.objects_root.rglob('*.o') if x.name.lower()!='main.o')
rsp=w/'objects.rsp';rsp.write_text('\n'.join('"'+x.as_posix()+'"' for x in objects),encoding='utf-8')
flags=['-std=c++23','-O2',*engine_arch_flags(a.objects_root),'-pthread','-fopenmp','-fno-access-control','-Isrc','-Isrc/eigen','-Isrc/gmxtraj/include']
if os.name=='nt':flags+=['-fuse-ld=lld']
exe=w/('oracle.exe' if os.name=='nt' else 'oracle');source='test/titan_dbscan_lj_geometry_oracle.cpp'
run('build',[a.compiler]+flags+[source,'@'+str(rsp),'-o',str(exe)]);run('run',[str(exe)])
record={'status':'PASS','binary_sha256':hashlib.sha256(exe.read_bytes()).hexdigest()}
if a.verify_regression:
 s=(root/'src/titan_llps_cluster_dbscan.cpp').read_text(encoding='utf-8');anchor='[&](std::size_t /*worker*/, std::size_t begin,\n                            std::size_t end)'
 assert s.count(anchor)==1
 s=s.replace(anchor,'[&](std::size_t begin, std::size_t end,\n                            std::size_t /*worker*/)')
 mutant=w/'old_callback.cpp';mutant.write_text(s,encoding='utf-8')
 mrsp=w/'mutant.rsp';mrsp.write_text('\n'.join('"'+x.as_posix()+'"' for x in objects if x.name!='titan_llps_cluster_dbscan.o'),encoding='utf-8')
 bad=w/('regression.exe' if os.name=='nt' else 'regression');run('regression_build',[a.compiler]+flags+[source,str(mutant),'@'+str(mrsp),'-o',str(bad)]);run('regression_run',[str(bad)],1);record['old_worker_callback']='EXPECTED_FAIL'
files=['src/titan_exact_periodic_box.h','src/titan_lj.cpp','src/titan_llps_cluster.cpp','src/titan_llps_cluster_dbscan.cpp','src/titan_llps_cluster_dynamics.cpp','src/titan_llps_cluster_rna_compaction.cpp',source]
record['source_sha256']={f:hashlib.sha256((root/f).read_bytes()).hexdigest() for f in files};record['object_sha256']={str(f.relative_to(root)):hashlib.sha256(f.read_bytes()).hexdigest() for f in objects}
(w/'verification.json').write_text(json.dumps(record,indent=2)+'\n',encoding='utf-8')
