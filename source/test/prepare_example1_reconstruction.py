"""Reconstruct the unavailable example/1 protein topology without changing protonation observed in PDB."""
import argparse,hashlib,itertools,json,subprocess
from pathlib import Path
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--pdb',type=Path,required=True);p.add_argument('--work-root',type=Path,required=True);p.add_argument('--gmx',default='gmx');a=p.parse_args();w=a.work_root.resolve();w.mkdir(parents=True,exist_ok=False)
lines=[]
for line in a.pdb.read_text().splitlines():
 if line.startswith(('ATOM','HETATM')) and line[17:20].strip()!='MOL':lines.append(line[:26]+' '+line[27:])
 elif line.startswith('TER'):lines.append(line)
 elif line.startswith('ENDMDL'):break
result=[];choices=[]
for key,group in itertools.groupby(lines,key=lambda line:line[17:27] if line.startswith(('ATOM','HETATM')) else line):
 group=list(group)
 if group[0].startswith(('ATOM','HETATM')) and group[0][17:20]=='HIS':
  names={x[12:16].strip() for x in group};d='HD1' in names;e='HE2' in names
  if not(d or e):raise ValueError('Histidine protonation cannot be recovered from observed hydrogen names')
  label='HIP' if d and e else 'HID' if d else 'HIE';choices.append({'residue':key,'template':label});group=[line[:17]+label+line[20:] for line in group]
 result+=group
(w/'protein.pdb').write_text('\n'.join(result)+'\nEND\n')
cmd=[a.gmx,'pdb2gmx','-f','protein.pdb','-o','protein_rebuilt.pdb','-p','protein.top','-ff','amber99sb-ildn','-water','none','-ignh']
r=subprocess.run(cmd,cwd=w,capture_output=True,text=True);(w/'pdb2gmx.log').write_text(r.stdout+r.stderr);r.check_returncode()
(w/'reconstruction.json').write_text(json.dumps({'scope':__doc__,'original_topology_available':False,'source_sha256':hashlib.sha256(a.pdb.read_bytes()).hexdigest(),'command':cmd,'histidine_choices':choices,'coordinates_policy':'Use input coordinates after atom-name/order and bond checks; generated hydrogen positions are not used by the comparison.'},indent=2)+'\n')
