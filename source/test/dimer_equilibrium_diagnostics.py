"""Energy-only NPT overlap/drift/correlation diagnostics, without reweighting fit.

Predeclared: 64/128 shared-range histogram bins; per-segment half-split drift;
energy ESS from an initial-positive, monotone autocorrelation sequence at 100 ps.
These diagnostics do not establish conformational equilibration or joint ESS.
"""
import argparse,collections,gzip,hashlib,json
from pathlib import Path
import numpy as np

def inefficiency(values):
 x=np.asarray(values,dtype=float);x=x-x.mean();n=len(x)
 if n<4 or np.dot(x,x)<=0:return float(n)
 size=1<<(2*n-1).bit_length();f=np.fft.rfft(x,n=size);c=np.fft.irfft(f*f.conj(),n=size)[:n];c/=np.arange(n,0,-1);rho=c/c[0]
 total=0.;previous=float('inf')
 for i in range(1,n-1,2):
  pair=float(rho[i]+rho[i+1])
  if pair<=0:break
  pair=min(pair,previous);total+=pair;previous=pair
 return min(float(n),max(1.,1.+2.*total))

def overlap(a,b,bins):
 edges=np.linspace(min(a.min(),b.min()),max(a.max(),b.max()),bins+1)
 if edges[0]==edges[-1]:return 1.
 x=np.histogram(a,bins=edges)[0]/len(a);y=np.histogram(b,bins=edges)[0]/len(b)
 return float(np.minimum(x,y).sum())

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--samples',type=Path,required=True);p.add_argument('--work-root',type=Path,required=True);a=p.parse_args();w=a.work_root.resolve();w.mkdir(parents=True,exist_ok=False)
 rng=np.random.default_rng(20260909);iid=rng.normal(size=100000);noise=rng.normal(size=100000);ar=np.empty_like(noise);ar[0]=noise[0]
 for i in range(1,len(ar)):ar[i]=.8*ar[i-1]+noise[i]
 gi=inefficiency(iid);ga=inefficiency(ar);assert 1<=gi<1.4 and abs(ga-9)<2
 assert overlap(np.array([0.,1.]),np.array([0.,1.]),64)==1 and overlap(np.array([0.,1.]),np.array([3.,4.]),64)==0
 data=collections.defaultdict(list);times=collections.defaultdict(list);temps={}
 with gzip.open(a.samples,'rt') as f:
  for line in f:
   r=json.loads(line);key=(r['slot'],r['segment']);data[key].append(r['configurational_enthalpy_kj_mol']);times[key].append(r['local_time_ps']);temps[r['slot']]=r['ref_temperature_K']
 assert len(data)==960
 persegment=[];perslot=[];phases={}
 for k in range(1,49):
  ess=0;maxdrift=0
  for s in range(1,21):
   values=np.asarray(data[k,s]);assert np.allclose(np.diff(times[k,s]),100,atol=1e-7)
   g=inefficiency(values);n=len(values);sd=float(values.std(ddof=1));drift=float((values[n//2:].mean()-values[:n//2].mean())/sd)
   persegment.append({'slot':k,'segment':s,'n':n,'mean_H_kj_mol':float(values.mean()),'sd_H_kj_mol':sd,'half_split_mean_shift_in_sd':drift,'energy_g':g,'energy_ess_diagnostic':n/g});ess+=n/g;maxdrift=max(maxdrift,abs(drift))
  allvalues=np.concatenate([data[k,s] for s in range(1,21)])
  for name,segs in [('segments_1_5',range(1,6)),('segments_6_20',range(6,21))]:phases[k,name]=np.concatenate([data[k,s] for s in segs])
  a1,a2=phases[k,'segments_1_5'],phases[k,'segments_6_20'];shift=float((a2.mean()-a1.mean())/allvalues.std(ddof=1))
  perslot.append({'slot':k,'temperature_K':temps[k],'n_exact':len(allvalues),'mean_H_kj_mol':float(allvalues.mean()),'sd_H_kj_mol':float(allvalues.std(ddof=1)),'sum_segment_energy_ess_diagnostic':ess,'max_abs_within_segment_half_shift_sd':maxdrift,'phase_6_20_minus_1_5_shift_sd':shift})
 adjacent=[]
 for phase in ['segments_1_5','segments_6_20']:
  for k in range(1,48):adjacent.append({'phase':phase,'slots':[k,k+1],'overlap_64_bins':overlap(phases[k,phase],phases[k+1,phase],64),'overlap_128_bins':overlap(phases[k,phase],phases[k+1,phase],128)})
 out={'status':'ENERGY_DIAGNOSTICS_ONLY_NOT_EQUILIBRATION_CERTIFICATE','oracles':{'iid_g':gi,'ar1_rho_0_8_g':ga,'ar1_expected_g':9,'histogram_identical_and_disjoint':'PASS'},'sample_sha256':hashlib.sha256(a.samples.read_bytes()).hexdigest(),'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'slots':perslot,'segments':persegment,'adjacent_overlap':adjacent,'limitations':['Energy overlap does not establish overlap of relevant conformations.','Segment ESS sums ignore dependence across segment boundaries and across exchanging slots; not a joint ESS or an uncertainty bound.','No burn-in time selected from these diagnostics; no fitted equilibrium weights/FES or rates.','Protocol changes and state-specific equilibration require separate assessment.']}
 (w/'diagnostics.json').write_text(json.dumps(out,indent=2)+'\n')
 print(json.dumps({'status':out['status'],'oracles':out['oracles'],'overlap64_min':min(r['overlap_64_bins'] for r in adjacent),'overlap128_min':min(r['overlap_128_bins'] for r in adjacent),'max_phase_shift_sd':max(abs(r['phase_6_20_minus_1_5_shift_sd']) for r in perslot),'max_half_shift_sd':max(r['max_abs_within_segment_half_shift_sd'] for r in perslot)},indent=2))

if __name__=='__main__':main()
