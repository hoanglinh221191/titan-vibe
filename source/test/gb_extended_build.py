"""Build the GB scientific driver against a freshly built TITAN engine.

Build TITAN first using make all PROFILE=native (or portable). This driver uses
the default Makefile configuration with the RF model disabled.
"""
import argparse, hashlib, json, os, subprocess
from pathlib import Path
from science_build_profile import engine_arch_flags
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--compiler',default='clang++' if os.name=='nt' else 'g++-13')
p.add_argument('--objects-root',type=Path,default=Path('build/native/src'))
p.add_argument('--work-root',type=Path,required=True)
a=p.parse_args();root=Path(__file__).resolve().parents[1]
w=a.work_root.resolve();w.mkdir(parents=True,exist_ok=False)
objects=sorted(f.resolve() for f in a.objects_root.rglob('*.o')
               if f.name.lower() not in ['main.o','titan_rf_model.o'])
assert any(f.name=='titan_rf_model_closed.o' for f in objects),'default RF-closed build required'
rsp=w/'objects.rsp';rsp.write_text('\n'.join('"'+f.as_posix()+'"' for f in objects))
exe=w/('driver.exe' if os.name=='nt' else 'driver')
cmd=[a.compiler,'-std=c++23','-O2','-Isrc','-Isrc/eigen','-Isrc/gmxtraj/include',
 *engine_arch_flags(a.objects_root),'-fopenmp','-pthread','test/titan_gb_extended_driver.cpp',
 '@'+str(rsp),'-o',str(exe)]
if os.name=='nt':cmd+=['-fuse-ld=lld']
with (w/'build.log').open('w') as log:subprocess.run(cmd,cwd=root,stdout=log,stderr=subprocess.STDOUT,check=True)
sha=lambda f:hashlib.sha256(f.read_bytes()).hexdigest()
(w/'build.json').write_text(json.dumps({'command':cmd,'driver_sha256':sha(exe),
 'objects_sha256':{str(f.relative_to(root)):sha(f) for f in objects},
 'source_sha256':{str(f.relative_to(root)):sha(f) for f in [root/'src/gbsa_intrin.cpp',
 root/'src/titan_gbn2.h',root/'src/titan_gbn2_tables.h',root/'src/titan_gbnsr6.h',
 root/'test/titan_gb_extended_driver.cpp']}},indent=2)+'\n')
print(exe)
