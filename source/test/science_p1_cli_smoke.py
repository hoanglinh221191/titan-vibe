"""Native CLI controls for scientific output guards on a tiny synthetic fixture."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import subprocess

root = Path(__file__).resolve().parents[1]
p = argparse.ArgumentParser(description=__doc__)
p.add_argument('--binary', type=Path, default=root / ('titan.exe' if (root/'titan.exe').exists() else 'titan'))
p.add_argument('--work-root', type=Path, default=root/'test/tmp/science_p1_cli')
a = p.parse_args()
binary = a.binary.resolve()
work = a.work_root.resolve()
work.mkdir(parents=True, exist_ok=True)
atoms = [('N',-4,0,0),('CA',-2.55,0,0),('C',-1.95,1.35,0),
         ('O',-2.55,2.35,0),('CB',-2.05,-0.75,1.25)]
pdb=['CRYST1   30.000   30.000   30.000  90.00  90.00  90.00 P 1           1']
for frame in (1,2):
    pdb += [f'TITLE     Synthetic energy control t= {frame-1:.5f}',f'MODEL     {frame:4d}']
    for chain,shift in [('A',0),('B',6.0)]:
        for index,(name,x,y,z) in enumerate(atoms,1):
            serial=index+(5 if chain=='B' else 0)
            residue=1 if chain=='A' else 2
            pdb.append(f'ATOM  {serial:5d} {name:^4s} ALA {chain}{residue:4d}    {x+shift:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           {name[0]}')
        pdb.append('TER')
    pdb.append('ENDMDL')
pdb.append('END')
(work/'fixture.pdb').write_text('\n'.join(pdb)+'\n',encoding='ascii')
(work/'fixture.ndx').write_text('[ System ]\n1 2 3 4 5 6 7 8 9 10\n[ Receptor ]\n1 2 3 4 5\n[ Ligand ]\n6 7 8 9 10\n',encoding='ascii')
# Deliberately synthetic parameters: this is a control-flow fixture, not a
# calibrated force field. A complete monolithic topology avoids local FF assets.
topology=['[ defaults ]','1 2 yes 0.5 0.833333333333',
          '[ atomtypes ]','X 6 12.0 0.0 A 0.3 0.2',
          '[ moleculetype ]','; name\tnrexcl','Fixture 3','[ atoms ]']
for i in range(1,11):
    topology.append(f'{i} X {1 if i<=5 else 2} ALA {atoms[(i-1)%5][0]} {i} {0.1 if i<=5 else -0.1} 12.0')
topology += ['[ system ]','Synthetic control','[ molecules ]','Fixture 1']
(work/'fixture.top').write_text('\n'.join(topology)+'\n',encoding='ascii')
lib=work.as_posix()+'/'
record={'binary_sha256':hashlib.sha256(binary.read_bytes()).hexdigest(),'cases':{}}
for name,group,extra,expected in [
    ('distinct_control','1-2',[],None),
    ('self_energy','1-1',[],None),
    ('self_force','1-1',['-cal-force'],None),
    ('self_decomposition','1-1',['-edecomp=true'],None),
    ('unsupported_self','1-1',[],'unsupported_topology_self_energy'),
    ('self_gbsa','1-1',['-gbsa=true'],'invalid_gbsa_identical_selections')]:
    case=work/name
    case.mkdir()
    options=[f'-f={(work/"fixture.pdb").as_posix()}','-1pdb=true','-ff=AMBER99SB',
             f'-n={(work/"fixture.ndx").as_posix()}',f'-g={group}',
             '-mol-1-type=protein','-mol-2-type=lig',f'-path-a={lib}',f'-path-n={lib}',
             f'-common_topo={(work/"fixture.top").as_posix()}','-nt=1']+extra
    if name=='unsupported_self':
        incomplete='\n'.join(topology).replace('1 2 yes 0.5 0.833333333333','1 2 no 0.5 0.833333333333')
        incomplete=incomplete.replace('[ system ]','[ pairs ]\n1 4 1\n[ system ]')
        (case/'incomplete.top').write_text(incomplete+'\n',encoding='utf-8')
        options=[f'-common_topo={(case/"incomplete.top").as_posix()}' if opt.startswith('-common_topo=') else opt for opt in options]
    (case/'arg.in').write_text('\n'.join(options)+'\n',encoding='utf-8')
    result=subprocess.run([str(binary),'-a',str(case/'arg.in')],cwd=case,
                          capture_output=True,text=True,encoding='utf-8',errors='replace',timeout=60)
    output=result.stdout+result.stderr
    (case/'result.log').write_text(output,encoding='utf-8')
    if expected is None:
        if result.returncode!=0: raise RuntimeError(f'{name} failed: {output[-2000:]}')
    elif result.returncode==0 or expected not in output:
        raise RuntimeError(f'{name} did not reject for {expected}: {output[-2000:]}')
    record['cases'][name]={'exit_code':result.returncode,'status':'PASS','diagnostic':expected}
    print(f'PASS {name}',flush=True)
for override in (False,True):
    name='gromos_override' if override else 'gromos_mixing_control'
    case=work/name;case.mkdir()
    gromos_top=topology.copy()
    gromos_top[1]='1 1 no 1.0 1.0'
    gromos_top[3]='X 6 12.0 0.0 A 0.001 0.000001'
    if override:gromos_top[4:4]=['[ nonbond_params ]','X X 1 0.004 0.000012']
    (case/'fixture.top').write_text('\n'.join(gromos_top)+'\n',encoding='ascii')
    options=[f'-f={(work/"fixture.pdb").as_posix()}','-1pdb=true','-ff=GROMOS54A7',
             f'-n={(work/"fixture.ndx").as_posix()}','-g=1-2',
             '-mol-1-type=protein','-mol-2-type=lig',f'-path-a={lib}',f'-path-n={lib}',
             f'-common_topo={(case/"fixture.top").as_posix()}','-nt=1']
    (case/'arg.in').write_text('\n'.join(options)+'\n',encoding='utf-8')
    result=subprocess.run([str(binary),'-a',str(case/'arg.in')],cwd=case,
                          capture_output=True,text=True,encoding='utf-8',errors='replace',timeout=60)
    output=result.stdout+result.stderr;(case/'result.log').write_text(output,encoding='utf-8')
    if result.returncode!=0:raise RuntimeError(f'{name}: {output[-2000:]}')
    c6,c12=(0.004,0.000012) if override else (0.001,0.000001)
    expected=sum(c12/(0.1*math.dist((x,y,z),(u+6,v,w)))**12-
                 c6/(0.1*math.dist((x,y,z),(u+6,v,w)))**6
                 for _,x,y,z in atoms for _,u,v,w in atoms)
    energies=[]
    for line in (case/'cal.xvg').read_text(encoding='utf-8').splitlines():
        if not line.strip() or line.lstrip().startswith(('#','@')):continue
        energies.append(float(line.split()[2]))
    if len(energies)!=2 or any(abs(e-expected)>1e-6*max(1,abs(expected)) for e in energies):
        raise RuntimeError(f'{name}: energy {energies} != independent pair sum {expected}')
    record['cases'][name]={'status':'PASS','vdw_kj_mol':energies,'expected_vdw_kj_mol':expected}
    print(f'PASS {name} parser -> type IDs -> pair table -> energy',flush=True)
(work/'verification.json').write_text(json.dumps(record,indent=2)+'\n',encoding='utf-8')
