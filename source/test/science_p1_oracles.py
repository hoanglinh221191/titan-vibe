"""Independent synthetic oracles for the September 2026 scientific corrections."""
import argparse
import hashlib
import json
import os
from pathlib import Path
from science_build_profile import engine_arch_flags
import subprocess

root = Path(__file__).resolve().parents[1]
p = argparse.ArgumentParser(description=__doc__)
p.add_argument('--compiler', default='clang++' if os.name == 'nt' else 'g++')
p.add_argument('--work-root', type=Path, default=root / 'test/tmp/science_p1_oracles')
p.add_argument('--verify-regression', action='store_true', help='Restore the old lag-buffer overwrite in a temporary source and require the oracle to fail')
p.add_argument('--case', default='all')
p.add_argument('--with-engine', action='store_true', help='Link integration oracles against a freshly built native engine')
p.add_argument('--objects-root', type=Path, default=root / 'build/native/src')
a = p.parse_args()
a.work_root.mkdir(parents=True, exist_ok=True)
exe = a.work_root.resolve() / ('science_p1.exe' if os.name == 'nt' else 'science_p1')
cmd = [a.compiler, '-std=c++23', '-O2', '-mavx2', '-mfma',
       '-Isrc', '-Isrc/eigen', '-pthread',
       'test/titan_science_p1_oracle.cpp', str(root / 'src/titan_kinetics.cpp'),
       str(root / 'src/titan_dham.cpp'), '-o', str(exe)]
subprocess.run(cmd, cwd=root, check=True)
result = subprocess.run([str(exe), a.case], cwd=a.work_root, capture_output=True, text=True)
(a.work_root / 'result.log').write_text(result.stdout + result.stderr, encoding='utf-8')
print(result.stdout + result.stderr, end='')
result.check_returncode()

mapping_exe=a.work_root.resolve()/('dham_mapping.exe' if os.name=='nt' else 'dham_mapping')
subprocess.run([a.compiler,'-std=c++23','-O2','-mavx2','-mfma','-pthread','-Isrc','-Isrc/eigen',
                'test/titan_dham_mapping_oracle.cpp','src/titan_dham.cpp','-o',str(mapping_exe)],cwd=root,check=True)
mapping=subprocess.run([str(mapping_exe)],cwd=a.work_root,capture_output=True,text=True)
(a.work_root/'dham_mapping.log').write_text(mapping.stdout+mapping.stderr,encoding='utf-8')
print(mapping.stdout+mapping.stderr,end='');mapping.check_returncode()
def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()
record={'core_binary_sha256':sha(exe),'core_status':'PASS',
        'dham_mapping_status':'PASS','dham_mapping_binary_sha256':sha(mapping_exe),
        'source_sha256':{str(f.relative_to(root)):sha(f) for f in
          [root/'src/titan_kinetics.cpp',root/'src/titan_dham.cpp',
           root/'src/titan_markov_estimator.h',root/'src/titan_qha_model.h',
           root/'src/box_type.h',root/'src/titan_neighbor_grid.h',
           root/'test/titan_science_p1_oracle.cpp']}}
if a.verify_regression:
    source=(root/'src/titan_kinetics.cpp').read_text(encoding='utf-8')
    anchor='    Eigen::Map<const Eigen::VectorXd> x_lagged(circ_buf_.row(circ_head_).data(),'
    if source.count(anchor)!=1: raise RuntimeError('Lag-buffer regression anchor changed')
    broken=a.work_root.resolve()/'old_lag_buffer.cpp'
    broken.write_text(source.replace(anchor,'    circ_buf_.row(circ_head_) = x.transpose();\n'+anchor),encoding='utf-8')
    broken_exe=a.work_root.resolve()/('old_lag_buffer.exe' if os.name=='nt' else 'old_lag_buffer')
    mutation_cmd=[str(broken) if part==str(root/'src/titan_kinetics.cpp') else
                  str(broken_exe) if part==str(exe) else part for part in cmd]
    subprocess.run(mutation_cmd,cwd=root,check=True)
    old=subprocess.run([str(broken_exe),'tica'],cwd=a.work_root,capture_output=True,text=True)
    (a.work_root/'old_lag_buffer_result.log').write_text(old.stdout+old.stderr,encoding='utf-8')
    if old.returncode==0 or 'P0 TICA lag buffer' not in old.stdout+old.stderr:
        raise RuntimeError('Old lag-buffer mutation did not fail the intended oracle')
    record['old_lag_buffer']='EXPECTED_FAIL'
    print('PASS regression sensitivity: old TICA lag-buffer overwrite rejected')
if a.with_engine:
    objects=[x.resolve() for x in a.objects_root.rglob('*.o') if x.name.lower()!='main.o']
    if not objects: raise RuntimeError('Build TITAN first')
    rsp=a.work_root.resolve()/'engine_objects.rsp'
    rsp.write_text('\n'.join('"'+x.as_posix()+'"' for x in objects),encoding='utf-8')
    exe=a.work_root.resolve()/('energy_p1.exe' if os.name=='nt' else 'energy_p1')
    cmd=[a.compiler,'-std=c++23','-O2',*engine_arch_flags(a.objects_root),'-fno-access-control',
         '-Isrc','-Isrc/eigen','-Isrc/gmxtraj/include','-pthread','-fopenmp',
         'test/titan_energy_p1_oracle.cpp','@'+str(rsp),'-o',str(exe)]
    if os.name=='nt':cmd.insert(1,'-fuse-ld=lld')
    subprocess.run(cmd,cwd=root,check=True)
    result=subprocess.run([str(exe)],cwd=a.work_root,capture_output=True,text=True)
    (a.work_root/'engine_result.log').write_text(result.stdout+result.stderr,encoding='utf-8')
    print(result.stdout+result.stderr,end='');result.check_returncode()
    record['engine_status']='PASS'
    record['engine_binary_sha256']=sha(exe)
    record['linked_objects_sha256']={str(x.relative_to(root)):sha(x) for x in objects}
(a.work_root/'verification.json').write_text(json.dumps(record,indent=2)+'\n',encoding='utf-8')
