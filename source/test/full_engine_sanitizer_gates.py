"""Representative cross-module gates linked to an entirely ASan/UBSan build.

This instruments all engine translation units, but does not imply execution or
validation of every engine path. Run after an isolated SANITIZERS build.
"""
import argparse,hashlib,json,os,subprocess
from pathlib import Path
from science_build_profile import engine_arch_flags
root=Path(__file__).resolve().parents[1]
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--objects-root',type=Path,required=True);p.add_argument('--binary',type=Path,required=True);p.add_argument('--compiler',default='g++-13');p.add_argument('--work-root',type=Path,required=True);a=p.parse_args();w=a.work_root.resolve();w.mkdir(parents=True,exist_ok=False)
objects=sorted(p.resolve() for p in a.objects_root.rglob('*.o') if p.name.lower()!='main.o');assert len(objects)>200
env=dict(os.environ,ASAN_OPTIONS='detect_leaks=1:halt_on_error=1',UBSAN_OPTIONS='halt_on_error=1:print_stacktrace=1')
def run(tag,cmd,cwd):
 r=subprocess.run(cmd,cwd=cwd,capture_output=True,text=True,encoding='utf-8',errors='replace',env=env);text=r.stdout+r.stderr;(w/(tag+'.log')).write_text(text,encoding='utf-8');r.check_returncode()
 if 'ERROR: AddressSanitizer' in text or 'runtime error:' in text:raise RuntimeError('Sanitizer diagnostic in '+tag)
 print(tag,'PASS',flush=True)
flags=['-std=c++23','-O1','-g1','-fsanitize=address,undefined','-fno-omit-frame-pointer','-fno-access-control',*engine_arch_flags(a.objects_root),'-fopenmp','-pthread','-Isrc','-Isrc/eigen','-Isrc/gmxtraj/include']
cases=[('dham_segment','test/titan_dham_segments_pipeline.cpp',[]),('paper_geometry','test/titan_llps_paper_geometry_oracle.cpp',['-DTITAN_TEST_ENGINE']),('gist_geometry','test/titan_gist_geometry_oracle.cpp',[]),('gist_thermo','test/titan_gist_thermo_oracle.cpp',['-DTITAN_TEST_ENGINE']),('alpb_zero_site','test/titan_alpb_zero_site_oracle.cpp',['-DTITAN_TEST_ENGINE']),('energy_p1','test/titan_energy_p1_oracle.cpp',[]),('dbscan_lj','test/titan_dbscan_lj_geometry_oracle.cpp',[]),('kinetics_p1','test/titan_science_p1_oracle.cpp',[])]
evidence=[]
run('abi_profile',['python3','test/science_build_profile_oracles.py','--compiler',a.compiler,'--objects-root',str(a.objects_root.resolve()),'--work-root',str(w/'abi_profile')],root)
for name,source,extra in cases:
 rsp=w/(name+'.rsp');rsp.write_text('\n'.join('"'+p.as_posix()+'"' for p in objects if name!='gist_geometry' or p.name!='titan_gist.o'))
 exe=w/name;run(name+'_build',[a.compiler]+flags+extra+[source,'@'+str(rsp),'-o',str(exe)],root)
 case=w/(name+'_run');case.mkdir();run(name+'_run',[str(exe)],case);evidence.append({'case':name,'status':'PASS'})
run('cli_validation',['python3','tools/test_cli_validation.py','--binary',str(a.binary.resolve())],root)
run('llps_cli',['python3','test/llps_phase_cli_smoke.py','--binary',str(a.binary.resolve()),'--work-root',str(w/'llps_cli')],root)
run('gist_cli',['python3','test/gist_cli_reference_smoke.py','--binary',str(a.binary.resolve()),'--driver',str(w/'gist_thermo'),'--work-root',str(w/'gist_cli')],root)
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
record={'status':'PASS','scope':__doc__,'cases':evidence,'cli_validation':'PASS','llps_cli_cases':7,'binary_sha256':sha(a.binary),'objects_sha256':{str(p.relative_to(root)):sha(p) for p in objects},'compiler_flags':flags,'ASAN_OPTIONS':env['ASAN_OPTIONS'],'UBSAN_OPTIONS':env['UBSAN_OPTIONS'],'script_sha256':sha(__file__)}
(w/'verification.json').write_text(json.dumps(record,indent=2)+'\n')
