"""Check Eigen ABI profile matching with two separately compiled translation units.

The probe only returns compile-time integer constants; it executes no SIMD
instructions, so AVX-512 alignment mismatch can be tested on an AVX2 host.
"""
import argparse,hashlib,json,os,subprocess
from pathlib import Path
from science_build_profile import engine_arch_flags
root=Path(__file__).resolve().parents[1]
p=argparse.ArgumentParser();p.add_argument('--compiler',required=True);p.add_argument('--objects-root',type=Path,default=root/'build/native/src');p.add_argument('--work-root',type=Path,required=True);a=p.parse_args();w=a.work_root.resolve();w.mkdir(parents=True,exist_ok=False)
engine=w/'engine.cpp';engine.write_text('#include <Eigen/Core>\nextern "C" int engine_alignment() { return EIGEN_DEFAULT_ALIGN_BYTES; }\n')
driver=w/'driver.cpp';driver.write_text('#include <Eigen/Core>\n#include <cstdio>\nextern "C" int engine_alignment();\nint main(){const int a=EIGEN_DEFAULT_ALIGN_BYTES,b=engine_alignment();std::printf("Eigen driver=%d engine=%d\\n",a,b);return a!=b;}\n')
records=[]
for name,driver_flags,engine_flags,expected in [
 ('current_profile',engine_arch_flags(a.objects_root),engine_arch_flags(a.objects_root),0),
 ('portable',['-mavx2','-mfma'],['-mavx2','-mfma'],0),
 ('avx512_match',['-mavx512f','-mfma'],['-mavx512f','-mfma'],0),
 ('mixed_avx2_avx512',['-mavx2','-mfma'],['-mavx512f','-mfma'],1)]:
 obj=w/(name+'.o');exe=w/(name+('.exe' if os.name=='nt' else ''));common=['-std=c++23','-O0','-Isrc/eigen']
 cmds=[[a.compiler]+common+engine_flags+['-c',str(engine),'-o',str(obj)],[a.compiler]+common+driver_flags+[str(driver),str(obj),'-o',str(exe)]]
 for command in cmds:
  r=subprocess.run(command,cwd=root,capture_output=True,text=True,encoding='utf-8',errors='replace');r.check_returncode()
 r=subprocess.run([str(exe)],capture_output=True,text=True,encoding='utf-8');assert r.returncode==expected,(name,r.stdout,r.stderr)
 if expected:assert 'driver=32 engine=64' in r.stdout
 (w/(name+'.log')).write_text(r.stdout+r.stderr,encoding='utf-8');records.append({'case':name,'status':'EXPECTED_FAIL' if expected else 'PASS','driver_flags':driver_flags,'engine_flags':engine_flags,'output':r.stdout.strip()})
fake=w/'invalid_profile';(fake/'generated').mkdir(parents=True);(fake/'generated/titan_build_manifest.json').write_text('{"profile":"unknown"}')
try:engine_arch_flags(fake/'src')
except ValueError:pass
else:raise AssertionError('Unknown engine ABI profile was accepted')
sha=lambda path:hashlib.sha256(Path(path).read_bytes()).hexdigest()
(w/'verification.json').write_text(json.dumps({'status':'PASS','cases':records,'unknown_profile':'REJECTED','helper_sha256':sha(root/'test/science_build_profile.py'),'script_sha256':sha(__file__)},indent=2)+'\n',encoding='utf-8');print('PASS Eigen ABI profile guard',records,flush=True)
