"""Audit reset segment clocks and replica exchanges; no kinetic fit or data modification."""
import argparse,subprocess,json,re,hashlib,collections
from pathlib import Path
def exchange_log(path):
 events=[];pending=None
 with path.open(errors='replace') as f:
  for line in f:
   m=re.search(r'Replica exchange at step (\d+) time ([\d.]+)',line)
   if m:pending=(int(m[1]),float(m[2]))
   elif line.startswith('Repl ex') and pending:
    tokens=line[7:].split();swaps=[]
    for i,t in enumerate(tokens):
     if t=='x':swaps.append((int(tokens[i-1]),int(tokens[i+1])))
    events.append((pending[0],pending[1],swaps));pending=None
 return events
def xvg(path):
 return [tuple(map(float,line.split())) for line in path.read_text().splitlines() if line and not line.startswith(('#','@'))]
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--dataset',type=Path,required=True);p.add_argument('--work-root',type=Path,required=True);p.add_argument('--xtc-identity',type=Path,required=True);a=p.parse_args();w=a.work_root.resolve();w.mkdir(parents=True,exist_ok=False)
 cat=xvg(a.dataset/'a1/energy_1.xvg');lookup=collections.defaultdict(list)
 for row in cat:lookup[(row[1],row[3])].append(row[0])
 identity=json.loads(a.xtc_identity.read_text());stem=lambda r:Path(r['source'].replace('\\','/')).stem
 parent_offset=next(float(next(iter(r['time_offsets_ps']))) for r in identity if stem(r)=='pro_cat6_20' and len(r['time_offsets_ps'])==1)
 xtc_offsets={int(stem(r)):float(next(iter(r['time_offsets_ps'])))+parent_offset for r in identity if stem(r).isdigit() and len(r['time_offsets_ps'])==1}
 summary=[]
 for segment in range(1,21):
  source=a.dataset/f'a1/{segment}.edr';out=w/f'{segment}.xvg'
  run=subprocess.run(['gmx','energy','-f',str(source),'-o',str(out)],input='Potential\nVolume\n0\n',capture_output=True,text=True,cwd=w)
  (w/f'{segment}_energy.log').write_text(run.stdout+run.stderr);run.check_returncode()
  energy=xvg(out);offsets=collections.Counter();missing=0
  for row in energy:
   hits=lookup[(row[1],row[2])]
   if len(hits)==1:offsets[round(hits[0]-row[0],6)]+=1
   else:missing+=1
  events=exchange_log(a.dataset/f'a1/{segment}.log');witness=exchange_log(a.dataset/f'a48/{segment}.log')
  digest=lambda rows:hashlib.sha256(json.dumps(rows,separators=(',',':')).encode()).hexdigest()
  visits=[0]*48;accepted=[0]*48
  for step,t,swaps in events:
   for i,j in swaps:accepted[i]+=1;accepted[j]+=1
  noexchange=0;intervals=0
  # Count only interior coordinate intervals (t,t+50]; exchange at the start
  # does not invalidate an interval if the stored coordinate is post-exchange.
  changed_by_interval=collections.defaultdict(set)
  for _,t,swaps in events:
   start=int((t-1e-6)//50)*50
   changed_by_interval[start].update(v for pair in swaps for v in pair)
  for start in range(50,int(energy[-1][0])-49,50):
   changed=changed_by_interval[start]
   noexchange+=48-len(changed);intervals+=48
  record=dict(segment=segment,energy_first_ps=energy[0][0],energy_last_ps=energy[-1][0],energy_offsets_ps=dict(offsets),unmatched_energy_rows=missing,
    exchange_events=len(events),exchange_first_ps=events[0][1] if events else None,exchange_last_ps=events[-1][1] if events else None,
    exchange_digest_a1=digest(events),exchange_digest_a48=digest(witness),witness_exact=events==witness,
    nonincreasing_exchange_times=sum(b[1]<=c[1] for c,b in zip(events,events[1:])),accepted_per_slot=accepted,
    candidate_50ps_intervals=intervals,exchange_free_50ps_intervals=noexchange,xtc_global_offset_ps=xtc_offsets.get(segment))
  if len(offsets)==1 and segment in xtc_offsets:record['energy_minus_xtc_clock_ps']=float(next(iter(offsets)))-xtc_offsets[segment]
  summary.append(record);(w/f'{segment}_exchange_events.json').write_text(json.dumps(events))
  print(segment,'energy offset',dict(offsets),'clock mismatch',record.get('energy_minus_xtc_clock_ps'),'exchange witness',record['witness_exact'],flush=True)
 (w/'audit.json').write_text(json.dumps({'scope':__doc__,'segments':summary,'decision':'NO_DIRECT_GLOBAL_TIMESTAMP_ALIGNMENT; segment-aware mapping required; exchange-free candidates do not certify Markov kinetics.'},indent=2)+'\n')
