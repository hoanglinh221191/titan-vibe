"""Topology nonbonded regression gate, optionally with native dispatch/ASan."""
import argparse,hashlib,json,os,subprocess
from pathlib import Path
from science_build_profile import engine_arch_flags
root=Path(__file__).resolve().parents[1]
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--compiler',default='clang++' if os.name=='nt' else 'g++')
p.add_argument('--work-root',type=Path,default=root/'test/tmp/ene02_oracles')
p.add_argument('--objects-root',type=Path,default=root/'build/native/src')
p.add_argument('--with-engine',action='store_true')
p.add_argument('--sanitize',action='store_true')
a=p.parse_args();work=a.work_root.resolve();work.mkdir(parents=True,exist_ok=False)
suffix='.exe' if os.name=='nt' else ''
def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def run(tag,cmd,cwd):
    result=subprocess.run(cmd,cwd=cwd,capture_output=True,text=True,encoding='utf-8',errors='replace')
    (work/(tag+'.log')).write_text(result.stdout+result.stderr,encoding='utf-8')
    if result.returncode:print(result.stdout+result.stderr);result.check_returncode()
    if tag.endswith('run'):print(result.stdout+result.stderr,end='',flush=True)
flags=['-std=c++23','-O1' if a.sanitize else '-O2','-mavx2','-mfma','-Isrc','-Isrc/eigen']
if a.sanitize:flags+=['-g','-fsanitize=address,undefined','-fno-omit-frame-pointer']
exe=work/('oracle'+suffix)
sources=['test/titan_nonbonded_topology_oracle.cpp','src/topology/nonbonded_topology.cpp','src/topology/nonbonded_row_kernel.cpp','src/titan_topology_energy.cpp']
run('unit_build',[a.compiler]+flags+sources+['-o',str(exe)],root)
run('unit_run',[str(exe)],work)
record={'unit_status':'PASS','sanitizers':a.sanitize,'unit_binary_sha256':sha(exe)}
if a.with_engine:
    if a.sanitize:raise ValueError('Engine objects require their own sanitized build; use unit sanitizer separately')
    objects=sorted(x.resolve() for x in a.objects_root.rglob('*.o') if x.name.lower()!='main.o')
    required={'nonbonded_topology.o','titan_topology_energy.o'}
    if not required.issubset({x.name for x in objects}):raise RuntimeError('Build current TITAN first')
    rsp=work/'objects.rsp';rsp.write_text('\n'.join('"'+x.as_posix()+'"' for x in objects),encoding='utf-8')
    driver=work/('dispatch'+suffix)
    cmd=[a.compiler]+[f for f in flags if f not in ('-mavx2','-mfma')]+engine_arch_flags(a.objects_root)+['-Isrc/gmxtraj/include','-pthread','-fopenmp',
         'test/titan_topology_dispatch_oracle.cpp','@'+str(rsp),'-o',str(driver)]
    if os.name=='nt':cmd.insert(1,'-fuse-ld=lld')
    run('dispatch_build',cmd,root);run('dispatch_run',[str(driver),str(work/'chain.top')],work)
    record['dispatch_status']='PASS';record['dispatch_binary_sha256']=sha(driver)
    record['objects_sha256']={str(x.relative_to(root)):sha(x) for x in objects}
bound=sources+['src/topology/nonbonded_topology.h','src/titan_topology_energy.h',
 'src/energy.cpp','src/energy_intrin.cpp','src/energy_intrin_nofilter.cpp','src/force.cpp',
 'src/force_intrin.cpp','src/cache_sigma_epsilon.h',*sorted(p.relative_to(root).as_posix() for p in (root/'src/engine').iterdir()),'src/titan_energy_decomposition.cpp',
 'test/titan_topology_dispatch_oracle.cpp']
record['sources_sha256']={s:sha(root/s) for s in bound}
(work/'verification.json').write_text(json.dumps(record,indent=2)+'\n',encoding='utf-8')
