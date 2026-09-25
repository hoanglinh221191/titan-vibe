"""Match stored concatenations to their source frames by compressed-coordinate identity."""
import argparse,hashlib,json,collections
from pathlib import Path
from xtc_segment_inventory import headers
def fingerprints(path):
 with path.open('rb') as f:
  for row in headers(path):
   f.seek(row['offset']);raw=f.read(row['end']-row['offset'])
   # Exclude the rewritten time/step, retain atom count, box and coordinate block.
   row['fingerprint']=hashlib.sha256(raw[4:8]+raw[16:]).hexdigest();yield row
def match(source,target):
 rows=list(fingerprints(target));index=collections.defaultdict(list)
 for i,r in enumerate(rows):index[r['fingerprint']].append(i)
 result=[]
 for path in source:
  offsets=collections.Counter();matched=0;missing=0;ambiguous=0;examples=[];first=None;last=None
  for r in fingerprints(path):
   if first is None:first=r['time_ps']
   last=r['time_ps'];hits=index.get(r['fingerprint'],[])
   if len(hits)==1:matched+=1;offsets[round(rows[hits[0]]['time_ps']-r['time_ps'],6)]+=1
   elif not hits:missing+=1
   else:ambiguous+=1
   if len(examples)<3:examples.append({'source_ps':r['time_ps'],'target_frames':hits})
  entry=dict(source=str(path),target=str(target),matched=matched,missing=missing,ambiguous=ambiguous,time_offsets_ps=dict(offsets),first_ps=first,last_ps=last,examples=examples);result.append(entry)
  print(path.name,'->',target.name,'matched',matched,'missing',missing,'offsets',dict(offsets),flush=True)
 return result
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--folder',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
 result=match([a.folder/'cat1_5.xtc',a.folder/'pro_cat6_20.xtc'],a.folder/'cat.xtc')
 result+=match(sorted((x for x in a.folder.glob('*.xtc') if x.stem.isdigit()),key=lambda x:int(x.stem)),a.folder/'cat6_20.xtc')
 with a.output.open('x',encoding='utf-8') as f:json.dump(result,f,indent=2);f.write('\n')
