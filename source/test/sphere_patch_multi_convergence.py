"""Run a bounded multi-atom convergence campaign, preserving failed cases."""
import argparse
import json
import pathlib
import subprocess
import time
import shutil

p=argparse.ArgumentParser()
p.add_argument('--driver',type=pathlib.Path,required=True)
p.add_argument('--output',type=pathlib.Path,required=True)
p.add_argument('--ligand',type=pathlib.Path,required=True)
a=p.parse_args();a.output.mkdir(parents=True,exist_ok=False)
cases={
    'disjoint2':[(0,0,0,1.7,1),(0,0,5,1.7,-1)],
    'disjoint3':[(0,0,0,1.7,1),(0,0,5,1.7,-.5),(0,0,10,1.7,.3)],
    'overlap2':[(0,0,0,1.7,1),(2.5,0,0,1.7,-1)],
    'overlap3':[(0,0,0,1.7,1),(2.5,0,0,1.7,-1),(1.25,2,0,1.7,.5)],
}
for name,atoms in cases.items():
    (a.output/(name+'.atoms')).write_text(str(len(atoms))+'\n'+'\n'.join(' '.join(map(str,b)) for b in atoms)+'\n')
shutil.copyfile(a.ligand,a.output/'ligand.atoms')
(a.output/'thresholds.json').write_text(json.dumps(dict(
    residual=2e-8,timeout_seconds=300,tree_theta=.25,
    tree_direct_relative=1e-4,finest_step_relative=.01,
    independent_oracle_relative=.01,
    note='mesh convergence is not an independent accuracy certificate'),indent=2))
rows=[]
schedule=[]
for name in ('disjoint2','disjoint3','overlap2','overlap3'):
    for points in (96,384,1536):schedule.append((name,points,'direct' if points<1536 else 'treecode'))
    schedule.append((name,384,'treecode'))
for points in (24,96,384,1536):schedule.append(('ligand',points,'treecode'))
for index,(name,points,mode) in enumerate(schedule):
    command=[str(a.driver.resolve()),str((a.output/(name+'.atoms')).resolve()),str(points),mode,'1e-8','.25']
    stem=a.output/f'{index:02d}_{name}_{points}_{mode}'
    stem.with_suffix('.command.json').write_text(json.dumps(command,indent=2))
    print('START',name,points,mode,flush=True);start=time.perf_counter()
    row=dict(fixture=name,points=points,mode=mode)
    with stem.with_suffix('.log').open('wb') as output:
        try:
            proc=subprocess.run(command,stdout=output,stderr=subprocess.STDOUT,timeout=300)
            row['returncode']=proc.returncode
        except subprocess.TimeoutExpired:row['timeout']=True
    text=stem.with_suffix('.log').read_text(errors='replace')
    for line in text.splitlines():
        if line.startswith('RESULT '):row.update(json.loads(line[7:]))
    row['elapsed_seconds']=time.perf_counter()-start
    rows.append(row);(a.output/'rows.json').write_text(json.dumps(rows,indent=2))
    print('DONE',json.dumps(row),flush=True)
complete=all(r.get('returncode')==0 and r.get('converged') and
             r.get('residual',1)>0 and r['residual']<2e-8 for r in rows)
(a.output/'completion.json').write_text(json.dumps(dict(
    all_scheduled_solves_complete=complete,
    note='Completion is not an energy-accuracy certificate; compare the independent oracle and refinement steps.'),indent=2))
raise SystemExit(0 if complete else 1)
