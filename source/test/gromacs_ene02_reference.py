"""Generate independent GROMACS single-frame nonbonded references (no MD)."""
import argparse
import hashlib
import json
import re
from pathlib import Path
import subprocess

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--work-root',type=Path,required=True)
p.add_argument('--gmx',default='gmx')
a=p.parse_args();work=a.work_root.resolve();work.mkdir(parents=True,exist_ok=False)
def run(directory,tag,arguments,stdin=None):
    result=subprocess.run([a.gmx]+arguments,input=stdin,cwd=directory,capture_output=True,text=True,timeout=60)
    (directory/(tag+'.log')).write_text(result.stdout+result.stderr,encoding='utf-8')
    if result.returncode:raise RuntimeError(f'{tag} failed in {directory}; see log')
    return result.stdout+result.stderr
positions=[(1.,1.,1.),(1.15,1.,1.),(1.30,1.,1.),(1.45,1.,1.),(1.10,1.15,1.02)]
charges=[.2,.1,-.1,-.4,.5]
mdp='''integrator = md
dt = 0.001
nsteps = 0
pbc = xyz
constraints = none
cutoff-scheme = Verlet
nstlist = 10
rlist = 4.0
verlet-buffer-tolerance = -1
rcoulomb = 4.0
rvdw = 4.0
coulombtype = Cut-off
coulomb-modifier = None
vdwtype = Cut-off
vdw-modifier = None
epsilon-r = 2
DispCorr = no
tcoupl = no
pcoupl = no
comm-mode = None
nstcalcenergy = 1
nstenergy = 1
nstfout = 1
nstxout = 1
nstvout = 0
nstlog = 1
'''
cases=[]
for rule in (1,2,3):
  for kind in ('generated','pairtypes','explicit','nbfix_generated'):
    case=work/f'rule{rule}_{kind}';case.mkdir()
    types=['X 6 12 0 A 0.004 0.000012','Y 6 12 0 A 0.009 0.000027'] if rule==1 else ['X 6 12 0 A 0.3 0.2','Y 6 12 0 A 0.4 0.8']
    top=['[ defaults ]',f'1 {rule} yes 0.5 0.8333333333333333','[ atomtypes ]']+types
    v,w=(0.02,0.00007) if rule==1 else (0.32,0.6)
    if kind=='pairtypes':top+=['[ pairtypes ]',f'X Y 1 {v} {w}']
    if kind=='nbfix_generated':top+=['[ nonbond_params ]',f'X Y 1 {v} {w}']
    top+=['[ moleculetype ]','Chain 3','[ atoms ]']
    for i in range(1,6):top.append(f'{i} {"X" if i in (1,2,5) else "Y"} 1 T C{i} {i} {charges[i-1]} 12')
    top+=['[ bonds ]','1 2 1 0.15 0','2 3 1 0.15 0','3 4 1 0.15 0','[ pairs ]',
          f'1 4 1 {v} {w}' if kind=='explicit' else '1 4 1',
          '[ system ]','Synthetic nonbonded oracle','[ molecules ]','Chain 1']
    (case/'topol.top').write_text('\n'.join(top)+'\n')
    gro=['Synthetic oracle','5']
    for i,(x,y,z) in enumerate(positions,1):gro.append(f'{1:5d}{"T":<5s}{"C"+str(i):>5s}{i:5d}{x:8.3f}{y:8.3f}{z:8.3f}')
    gro.append('10.0 10.0 10.0');(case/'coords.gro').write_text('\n'.join(gro)+'\n')
    (case/'rerun.mdp').write_text(mdp)
    run(case,'grompp',['grompp','-f','rerun.mdp','-c','coords.gro','-p','topol.top','-pp','processed.top','-o','rerun.tpr'])
    run(case,'rerun',['mdrun','-s','rerun.tpr','-rerun','coords.gro','-deffnm','rerun','-nt','1'])
    run(case,'energy',['energy','-f','rerun.edr','-o','energies.xvg'],
        'Coulomb-(SR)\nLJ-(SR)\nCoulomb-14\nLJ-14\n0\n')
    # Request per-atom forces, not the sum over a group.
    run(case,'force',['traj','-s','rerun.tpr','-f','rerun.trr','-of','forces.xvg','-xvg','none'],'0\n')
    energy_lines=(case/'energies.xvg').read_text().splitlines()
    energy=[float(x) for line in energy_lines if line.strip() and not line.startswith(('#','@')) for x in line.split()]
    columns=['Time']+[m[1] for line in energy_lines if (m:=re.match(r'@ s\d+ legend "([^"]+)"',line))]
    force=[float(x) for x in (case/'forces.xvg').read_text().splitlines()[0].split()]
    cases.append({'name':case.name,'energy_columns':columns,'energy_row':energy,'force_row':force,
                  'topology_sha256':hashlib.sha256((case/'processed.top').read_bytes()).hexdigest()})
    print('PASS GROMACS reference '+case.name,flush=True)
(work/'reference.json').write_text(json.dumps({'cases':cases,'positions_nm':positions,'dielectric':2,
    'reference_scope':'single-frame plain cutoff nonbonded terms; zero bonded force constants; no dynamics'},indent=2)+'\n')
