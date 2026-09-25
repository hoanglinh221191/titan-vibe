"""Independent GROMACS decode and exhaustive-lattice geometry spot checks."""
import argparse,hashlib,itertools,json,re,subprocess
from pathlib import Path
import numpy as np
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def descriptors(x,box,selection):
 inv=np.linalg.inv(box);shifts=np.array(list(itertools.product(range(-2,3),repeat=3)))
 def mic(d):
  center=np.rint(d@inv);images=d-(center+shifts)@box;norm=np.sum(images*images,axis=1);order=np.argsort(norm);assert norm[order[1]]-norm[order[0]]>1e-10;return images[order[0]]
 com=[];rg=[];end=[];maximum=0.
 for ids in selection:
  chain=x[ids];lift=np.zeros_like(chain);lift[0]=chain[0]
  for i in range(1,len(chain)):
   d=mic(chain[i]-chain[i-1]);maximum=max(maximum,float(np.linalg.norm(d)));assert .2<np.linalg.norm(d)<.5;lift[i]=lift[i-1]+d
  center=lift.mean(axis=0);com.append(center);rg.append(np.sqrt(np.mean(np.sum((lift-center)**2,axis=1))));end.append(np.linalg.norm(lift[-1]-lift[0]))
 return np.array([np.linalg.norm(mic(com[1]-com[0])),min(rg),max(rg),min(end),max(end),maximum])
def main():
 p=argparse.ArgumentParser();p.add_argument('--dataset',type=Path,required=True);p.add_argument('--features',type=Path,required=True);p.add_argument('--gmx',type=Path,required=True);p.add_argument('--work-root',type=Path,required=True);a=p.parse_args();w=a.work_root.resolve();w.mkdir(parents=True,exist_ok=False)
 sel=np.loadtxt(a.features/'ca_selection.txt',dtype=int);selection=[sel[sel[:,0]==c,1] for c in range(2)];records=[]
 for k in [1,48]:
  data=np.genfromtxt(a.features/f'a{k}.csv',delimiter=',',names=True)
  for index in [0,10620,21239]:
   name=f'a{k}_f{index}';out=w/(name+'.gro');cmd=[str(a.gmx),'trjconv','-f',str(a.dataset/f'a{k}/cat.xtc'),'-s',str(a.dataset/f'a{k}/1.tpr'),'-dump',str(data['time_ps'][index]),'-o',str(out),'-ndec','6']
   r=subprocess.run(cmd,input='Protein\n',capture_output=True,text=True);(w/(name+'.log')).write_text(r.stdout+r.stderr);r.check_returncode()
   lines=out.read_text().splitlines();assert int(lines[1])==1254;assert abs(float(re.search(r't=\s*([\d.+-]+)',lines[0])[1])-data['time_ps'][index])<1e-4
   # GROMACS may retain the input precision despite -ndec. Infer fixed field
   # width from decimal positions, as required by the GRO format.
   dots=[i for i,c in enumerate(lines[2][20:]) if c=='.'];width=dots[1]-dots[0]
   x=np.array([[float(line[20+width*j:20+width*(j+1)]) for j in range(3)] for line in lines[2:1256]])
   cell=list(map(float,lines[-1].split()));assert len(cell) in [3,9]
   if len(cell)==3:box=np.diag(cell)
   else:box=np.array([[cell[0],cell[3],cell[4]],[cell[5],cell[1],cell[6]],[cell[7],cell[8],cell[2]]])
   result=descriptors(x,box,selection);expected=np.array([data[n][index] for n in data.dtype.names[4:]])
   error=np.abs(result-expected);assert error.max()<1e-4,(name,error)
   shifted=x.copy();rng=np.random.default_rng(831);shifted+=rng.integers(-3,4,size=(len(x),3))@box;assert np.max(np.abs(descriptors(shifted,box,selection[::-1])-result))<1e-9
   records.append({'slot':k,'frame':index,'max_error_nm':float(error.max()),'gro_sha256':sha(out),'raw_decode':'GROMACS trjconv','geometry_reference':'independent NumPy exhaustive lattice enumeration; per-CA image and chain-swap invariance'})
 (w/'verification.json').write_text(json.dumps({'status':'PASS','scope':'Six decoded frames at low/high temperature; not all-frame independent coordinate identity','records':records,'feature_manifest_sha256':sha(a.features/'manifest.json'),'script_sha256':sha(__file__)},indent=2)+'\n')
 print('PASS six independent geometry checks',max(r['max_error_nm'] for r in records))
if __name__=='__main__':main()
