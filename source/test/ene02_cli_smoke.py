"""Check native self-energy/decomposition against the frozen GROMACS references."""
import argparse,csv,json,re,subprocess
from pathlib import Path
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--reference',type=Path,required=True)
p.add_argument('--binary',type=Path,required=True)
p.add_argument('--work-root',type=Path,required=True)
a=p.parse_args();work=a.work_root.resolve();work.mkdir(parents=True,exist_ok=False)
ref=json.loads(a.reference.read_text());binary=a.binary.resolve();results=[]
for entry in ref['cases']:
    for mode in ('stream','cache'):
        case=work/(entry['name']+'_'+mode);case.mkdir()
        atoms=[]
        for frame in (1,2):
            atoms+=['CRYST1  100.000  100.000  100.000  90.00  90.00  90.00 P 1           1',
                    f'TITLE     Synthetic t= {frame-1:.5f}',f'MODEL     {frame:4d}']
            for i,(x,y,z) in enumerate(ref['positions_nm'],1):
                atoms.append(f'ATOM  {i:5d} {"C"+str(i):^4s} {"T":>3s} A{1:4d}    {10*x:8.3f}{10*y:8.3f}{10*z:8.3f}  1.00  0.00           C')
            atoms+=['TER','ENDMDL']
        atoms.append('END');(case/'coords.pdb').write_text('\n'.join(atoms)+'\n')
        (case/'selection.ndx').write_text('[ System ]\n1 2 3 4 5\n[ Selected ]\n1 2 3 4 5\n')
        top=(a.reference.parent/entry['name']/'processed.top').resolve()
        ff={'rule1':'GROMOS54A7','rule2':'AMBER99SB','rule3':'OPLSAA'}[entry['name'].split('_')[0]]
        args=[f'-f={(case/"coords.pdb").as_posix()}','-1pdb=true',f'-ff={ff}',
              f'-common_topo={top.as_posix()}',f'-n={(case/"selection.ndx").as_posix()}',
              '-g=1-1','-epsilon_in=2','-mol-1-type=protein','-mol-2-type=lig',
              f'-path-a={case.as_posix()}/',f'-path-n={case.as_posix()}/',
              '-edecomp=true','-cal-force','-nt=2',f'-cache={"no" if mode=="stream" else "yes"}']
        (case/'arg.in').write_text('\n'.join(args)+'\n')
        run=subprocess.run([str(binary),'-a',str(case/'arg.in')],cwd=case,capture_output=True,text=True,encoding='utf-8',errors='replace',timeout=60)
        (case/'result.log').write_text(run.stdout+run.stderr,encoding='utf-8')
        if run.returncode:raise RuntimeError(f'{case.name} runtime failed; see result.log')
        allocation=re.search(r'Pair Cache Entries\s*:\s*(\d+)',(case/'log.txt').read_text(encoding='utf-8'))
        if not allocation or int(allocation[1])!=(0 if mode=='stream' else 25):
            raise RuntimeError(f'{case.name} did not use the requested allocation mode')
        values=dict(zip(entry['energy_columns'],entry['energy_row']))
        expected_e=values['Coulomb (SR)']+values['Coulomb-14']
        expected_v=values['LJ (SR)']+values['LJ-14']
        energies=[list(map(float,line.split())) for line in (case/'cal.xvg').read_text().splitlines() if line.strip() and not line.lstrip().startswith(('#','@'))]
        if len(energies)!=2:raise RuntimeError('Expected two frames')
        # CLI output is rounded and GROMACS inputs/computation use single precision.
        for row in energies:
            if abs(row[1]-expected_e)>2e-5*max(1,abs(expected_e)) or abs(row[2]-expected_v)>2e-5*max(1,abs(expected_v)):
                raise RuntimeError(f'{case.name} energy differs: {row} vs {expected_e},{expected_v}')
        decomposition=[]
        for suffix,expected in [('elec',energies[0][1]),('vdw',energies[0][2])]:
            matches=list(case.rglob(f'edecomp_mm_residue_matrix_{suffix}.csv'))
            if len(matches)!=1:raise RuntimeError(f'{case.name} missing matrix {suffix}')
            rows=list(csv.DictReader(matches[0].read_text().splitlines()))
            total=sum(float(row['signed_mean']) for row in rows)
            if abs(total-expected)>1e-6*max(1,abs(expected)):raise RuntimeError(f'{case.name} decomposition closure: {total} vs {expected}')
            decomposition.append(total)
        results.append({'case':case.name,'status':'PASS','energy':energies,'decomposition':decomposition})
        print('PASS native '+case.name,flush=True)
(work/'verification.json').write_text(json.dumps({'status':'PASS','cases':results},indent=2)+'\n')
