"""Hash-bound GIST scaling benchmark with exact-output gates.

Relative timing floors are deliberately loose gross-regression sentinels, not
portable speed guarantees. Keep raw medians, CPU and compiler with each run.
"""
import argparse,csv,hashlib,json,os,platform,subprocess
from pathlib import Path
from science_build_profile import engine_arch_flags
root=Path(__file__).resolve().parents[1]
p=argparse.ArgumentParser();p.add_argument('--compiler',default='g++-13');p.add_argument('--objects-root',type=Path,default=root/'build/native/src');p.add_argument('--work-root',type=Path,required=True);a=p.parse_args();w=a.work_root.resolve();w.mkdir(parents=True,exist_ok=False)
objects=sorted(p.resolve() for p in a.objects_root.rglob('*.o') if p.name.lower()!='main.o');assert len(objects)>200;rsp=w/'objects.rsp';rsp.write_text('\n'.join('"'+p.as_posix()+'"' for p in objects));exe=w/('benchmark.exe' if os.name=='nt' else 'benchmark')
flags=['-std=c++23','-O2','-fno-access-control','-Isrc','-Isrc/eigen','-Isrc/gmxtraj/include',*engine_arch_flags(a.objects_root),'-fopenmp','-pthread']
if os.name=='nt':flags+=['-fuse-ld=lld']
r=subprocess.run([a.compiler]+flags+['test/titan_gist_benchmark.cpp','@'+str(rsp),'-o',str(exe)],cwd=root,capture_output=True,text=True,encoding='utf-8',errors='replace');(w/'build.log').write_text(r.stdout+r.stderr,encoding='utf-8');r.check_returncode()
r=subprocess.run([str(exe)],cwd=w,capture_output=True,text=True,encoding='utf-8',errors='replace');(w/'run.log').write_text(r.stdout+r.stderr,encoding='utf-8');r.check_returncode();lines=[l for l in r.stdout.splitlines() if l.startswith(('kind,','orientation,','Eww,'))];rows=list(csv.DictReader(lines));assert len(rows)==6
for row in rows:
 for key in list(row)[1:]:row[key]=float(row[key])
 assert row['max_absolute_error']<=1e-8
large=[r for r in rows if (r['kind'],int(r['n'])) in [('orientation',16384),('Eww',2744)]]
# Fixed before timing: reject catastrophic reversals, tolerate noisy shared CI.
assert all(r['speedup']>=.5 for r in large),large
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
record={'status':'PASS','rows':rows,'timing_method':'one correctness warmup, median of three runs per algorithm, single calling thread','timing_gate':'largest workloads must be at least 0.5x exhaustive speed; gross regression sentinel only','host':platform.platform(),'processor':platform.processor(),'compiler_version':subprocess.check_output([a.compiler,'--version'],text=True).splitlines()[0],'flags':flags,'binary_sha256':sha(exe),'sources_sha256':{p:sha(root/p) for p in ['src/titan_gist.cpp','src/titan_gist.h','src/titan_gist_entropy.h','test/titan_gist_benchmark.cpp','test/gist_performance_gate.py']},'objects_sha256':{str(p.relative_to(root)):sha(p) for p in objects}}
if Path('/proc/cpuinfo').exists():record['cpu_model']=next(l.split(':',1)[1].strip() for l in Path('/proc/cpuinfo').read_text().splitlines() if l.startswith('model name'))
(w/'verification.json').write_text(json.dumps(record,indent=2)+'\n');print(json.dumps(rows,indent=2),flush=True)
