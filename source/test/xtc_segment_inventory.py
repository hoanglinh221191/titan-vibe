"""Read XTC frame headers without decoding or modifying the trajectory.

Offsets, atom counts, simulation steps, timestamps and cell matrices are retained.
Coordinate identity requires a separate atom-matched decoded comparison.
"""
import argparse,json,struct,hashlib,collections
from pathlib import Path
def headers(path):
 size=path.stat().st_size
 with path.open('rb') as f:
  while f.tell()<size:
   offset=f.tell();raw=f.read(52)
   if len(raw)!=52:raise ValueError(f'{path}: truncated frame header at {offset}')
   magic,natoms,step,time,*cell=struct.unpack('>iiif9f',raw)
   if magic!=1995 or natoms<=0:raise ValueError(f'{path}: invalid XTC header at {offset}')
   count=struct.unpack('>i',f.read(4))[0]
   if count!=natoms:raise ValueError(f'{path}: atom count mismatch')
   if count<=9:f.seek(count*12,1)
   else:
    compressed=f.read(36)
    if len(compressed)!=36:raise ValueError('truncated compressed header')
    precision,*others=struct.unpack('>f8i',compressed);nbytes=others[-1]
    if nbytes<0:raise ValueError('negative XTC payload size')
    f.seek((nbytes+3)//4*4,1)
   if f.tell()>size:raise ValueError('truncated XTC payload')
   yield dict(offset=offset,end=f.tell(),natoms=natoms,step=step,time_ps=time,cell_nm=cell)
def summarize(path):
 rows=list(headers(path));dt=collections.Counter(round(b['time_ps']-a['time_ps'],6) for a,b in zip(rows,rows[1:]))
 return dict(path=str(path),bytes=path.stat().st_size,nframes=len(rows),atom_counts=sorted(set(x['natoms'] for x in rows)),first=rows[0] if rows else None,last=rows[-1] if rows else None,intervals_ps=dict(dt),nonincreasing=sum(b['time_ps']<=a['time_ps'] for a,b in zip(rows,rows[1:])),header_sha256=hashlib.sha256(json.dumps(rows,sort_keys=True).encode()).hexdigest(),frames=rows)
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--folder',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
 result=[]
 for path in sorted(a.folder.glob('*.xtc')):
  row=summarize(path);result.append(row)
  print(path.name,row['atom_counts'],row['nframes'],row['first']['time_ps'],row['last']['time_ps'],row['intervals_ps'],flush=True)
 with a.output.open('x',encoding='utf-8') as f:json.dump(result,f);f.write('\n')
