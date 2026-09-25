"""Read all 48 protein XTCs and bind descriptors to prior segment alignment."""
import argparse,concurrent.futures,gzip,hashlib,json,subprocess,time
from pathlib import Path
import numpy as np
root=Path(__file__).resolve().parents[1]
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(4*1024*1024),b''):h.update(b)
 return h.hexdigest()
def main():
 p=argparse.ArgumentParser();p.add_argument('--dataset',type=Path,required=True);p.add_argument('--work-root',type=Path,required=True);p.add_argument('--workers',type=int,default=4);a=p.parse_args();w=a.work_root.resolve();w.mkdir(parents=True,exist_ok=False)
 # Selection is atom-name matched against the actual source GRO and topology.
 selections=[];atom_names=[];offset=0
 for c,name in enumerate(['topol_Protein_chain_A.itp','topol_Protein_chain_A2.itp']):
  section='';atoms=[]
  for line in (a.dataset/name).read_text().splitlines():
   line=line.split(';')[0].strip()
   if line.startswith('['):section=line.strip('[] ').strip();continue
   if section=='atoms' and line:
    f=line.split();atoms.append((int(f[0]),int(f[2]),f[4]));atom_names.append((f[3],f[4]))
  assert len(atoms)==627 and [i for i,r,n in atoms]==list(range(1,628))
  selections.extend((c,offset+i-1) for i,r,n in atoms if n=='CA');offset+=len(atoms)
 gro=(a.dataset/'a1/1.gro').read_text().splitlines()[2:1256]
 assert [(l[5:10].strip(),l[10:15].strip()) for l in gro]==atom_names
 sel=w/'ca_selection.txt';sel.write_text(''.join(f'{c} {i}\n' for c,i in selections));assert len(selections)==84
 exe=w/'features';cmd=['g++-13','-std=c++23','-O3','-Isrc','-Isrc/gmxtraj/include','test/dimer_conformation_features.cpp','src/pbctopo/lattice_math.cpp','src/gmxtraj/xdrfile.cpp','src/gmxtraj/xdrfile_xtc.cpp','-o',str(exe)]
 # The XDR sources in this tree use .cpp; compile the actual production reader.
 with (w/'build.log').open('w') as f:subprocess.run(cmd,cwd=root,stdout=f,stderr=subprocess.STDOUT,check=True)
 def slot(k):
  trajectory=a.dataset/f'a{k}/cat.xtc';auditroot=root/'test/tmp'/('dimer_alignment_pilot_r2_20260909' if k==2 else 'dimer_alignment_remaining_20260909')/f'a{k}'
  audit=json.loads((auditroot/'audit.json').read_text());alignment=auditroot/'alignment.jsonl.gz';assert sha(alignment)==audit['alignment_sha256'];assert trajectory.stat().st_size==audit['xtc_bytes']
  before=sha(trajectory);out=w/f'a{k}.csv';start=time.perf_counter();r=subprocess.run([str(exe),str(trajectory),str(sel),str(out)],capture_output=True,text=True);(w/f'a{k}.log').write_text(r.stdout+r.stderr);r.check_returncode()
  data=np.genfromtxt(out,delimiter=',',names=True);assert len(data)==21240
  rows=[json.loads(l) for l in gzip.open(alignment,'rt')];assert all(int(x['frame'])==row['frame_zero_based'] and int(x['step'])==row['local_step'] and abs(x['time_ps']-row['xtc_global_time_ps'])<1e-6 for x,row in zip(data,rows))
  # Align exact samples without interpolation; retain absent energies as NaN.
  np.savez_compressed(w/f'a{k}.npz',features=np.column_stack([data[n] for n in data.dtype.names]),segment=np.array([r['segment'] for r in rows]),local_ps=np.array([r['local_time_ps'] for r in rows]),energy=np.array([r['potential_kj_mol'] if r['potential_kj_mol'] is not None else np.nan for r in rows]),volume=np.array([r['xtc_volume_nm3'] for r in rows]))
  assert sha(trajectory)==before
  record={'slot':k,'frames':len(data),'valid_frames':int(data['valid'].sum()),'elapsed_s':time.perf_counter()-start,'trajectory_sha256':before,'alignment_sha256':sha(alignment),'csv_sha256':sha(out),'npz_sha256':sha(w/f'a{k}.npz')};print(k,record['valid_frames'],flush=True);return record
 with concurrent.futures.ThreadPoolExecutor(max_workers=a.workers) as pool:records=list(pool.map(slot,range(1,49)))
 result={'status':'PASS_EXTRACTION','scope':'CA geometry descriptors only; not state/Markov/equilibrium certification','columns':list(np.genfromtxt(w/'a1.csv',delimiter=',',names=True,max_rows=1).dtype.names),'records':records,'selection_sha256':sha(sel),'driver_sha256':sha(exe),'script_sha256':sha(__file__),'driver_source_sha256':sha(root/'test/dimer_conformation_features.cpp')}
 (w/'manifest.json').write_text(json.dumps(result,indent=2)+'\n')
if __name__=='__main__':main()
