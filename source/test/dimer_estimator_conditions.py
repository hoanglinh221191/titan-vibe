"""Predeclared conformation/censoring/state-bias diagnostics; no physical fit.

States are geometric screens (12/36 bins), never certified MSM states. Drift,
overlap, ESS and finite-sample uncertainties are descriptive, not p-values from
independent replicas. Temperature slots participate in the same exchange run.
"""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
root=Path(__file__).resolve().parents[1]
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def js(a,b):
 a=np.asarray(a,dtype=float);b=np.asarray(b,dtype=float)
 if not a.sum() or not b.sum():return None
 a/=a.sum();b/=b.sum();m=(a+b)/2
 return float(.5*sum(np.sum(p[p>0]*np.log(p[p>0]/m[p>0])) for p in [a,b]))
def gtime(x):
 x=np.asarray(x,float);x=x-x.mean();n=len(x)
 if n<4 or not np.dot(x,x):return 1.
 m=1<<(2*n-1).bit_length();f=np.fft.rfft(x,m);ac=np.fft.irfft(f*f.conj(),m)[:n];ac/=np.arange(n,0,-1);ac/=ac[0]
 pairs=[]
 for i in range(1,min(n-1,n//2),2):
  pair=float(ac[i]+ac[i+1])
  if pair<=0:break
  pairs.append(min(pair,pairs[-1]) if pairs else pair)
 return max(1.,1+2*sum(pairs))
def state_labels(x,fine):
 distance=np.digitize(x[:,4],[1.,1.5,2.,2.5,3.] if fine else [1.,2.,3.])
 rg=np.digitize(x[:,6],[.8,1.,1.2,1.5,2.] if fine else [1.,1.5])
 return distance*(6 if fine else 3)+rg
def lme(x):
 top=float(x.max());return float(top+np.log(np.exp(x-top).mean()))
def main():
 p=argparse.ArgumentParser();p.add_argument('--features',type=Path,required=True);p.add_argument('--work-root',type=Path,required=True);a=p.parse_args();w=a.work_root.resolve();w.mkdir(parents=True,exist_ok=False)
 manifest=json.loads((a.features/'manifest.json').read_text());assert manifest['status']=='PASS_EXTRACTION'
 eligibility=json.loads((root/'test/tmp/dimer_dham_eligibility_20260909/eligibility.json').read_text());temps={r['slot']:r['ref_temperature_K'] for r in eligibility['slot_records']}
 events={};event_hashes={}
 for s in range(1,21):
  path=root/f'test/tmp/dimer_exchange_energy_r2_20260908/{s}_exchange_events.json';j=json.loads(path.read_text());event_hashes[str(path.relative_to(root))]=sha(path)
  # Retain the previously audited a1/a48 global exchange witnesses.
  for k in range(1,49):events[k,s]=np.array([t for _,t,pairs in j if any(k-1 in pair for pair in pairs)])
 rng=np.random.default_rng(271828);iid=rng.normal(size=50000);ar=np.zeros(100000);noise=rng.normal(size=len(ar))
 for i in range(1,len(ar)):ar[i]=.8*ar[i-1]+noise[i]
 oracle={'iid_g':gtime(iid),'ar08_g':gtime(ar),'identical_js':js([1,2],[1,2]),'disjoint_js':js([1,0],[0,1])};assert abs(oracle['iid_g']-1)<.15 and abs(oracle['ar08_g']-9)<2 and abs(oracle['disjoint_js']-np.log(2))<1e-12
 records=[];histograms={};bias=[];censor=[];total_frames=0
 for rec in manifest['records']:
  k=rec['slot'];path=a.features/f'a{k}.npz';assert sha(path)==rec['npz_sha256'];d=np.load(path);x=d['features'];segments=d['segment'];local=d['local_ps'];valid=x[:,3]==1;total_frames+=len(x)
  H=d['energy']+.0602214076*d['volume'];dbeta=1/(.00831446261815324*temps[k])-1/(.00831446261815324*300)
  for regime,(lo,hi) in enumerate([(1,9),(10,13),(14,20)]):
   mask=(segments>=lo)&(segments<=hi)&valid;y=x[mask];assert len(y)>20
   desc=[]
   for col in range(4,9):
    z=y[:,col];half=len(z)//2;sd=float(z.std());parts=np.array_split(z,4)
    # ESS estimates are per-segment only; do not sum into a joint replica ESS.
    ess=[{'segment':s,'n':int(np.sum((segments==s)&valid)),'g':gtime(x[(segments==s)&valid,col])} for s in range(lo,hi+1)]
    desc.append({'column':manifest['columns'][col],'mean':float(z.mean()),'sd':sd,'half_mean_shift_over_sd':float(abs(z[:half].mean()-z[half:].mean())/sd) if sd else 0.,'quarter_means':[float(p.mean()) for p in parts],'segment_correlation':ess})
   records.append({'slot':k,'temperature_K':temps[k],'regime':regime,'segments':[lo,hi],'n':len(y),'descriptors':desc})
   for fine in [False,True]:
    ns=36 if fine else 12;labels=state_labels(x,fine);z=labels[mask];half=len(z)//2;hist=np.bincount(z,minlength=ns);histograms[k,regime,ns]=hist
    records[-1].setdefault('state_screens',[]).append({'states':ns,'occupied':int(np.sum(hist>0)),'half_js_nats':js(np.bincount(z[:half],minlength=ns),np.bincount(z[half:],minlength=ns))})
    for state in np.flatnonzero(hist):
     es=H[mask&(labels==state)&np.isfinite(H)]
     if len(es)<20:continue
     centered=-dbeta*(es-es.mean());bias.append({'slot':k,'regime':regime,'state_count':ns,'state':int(state),'n_exact':len(es),'reduced_bias_sd':float(centered.std()),'logmeanexp_minus_mean_bias':lme(centered)})
    for lag in [50,100,200,500]:
     before=np.zeros(ns,dtype=int);retained=np.zeros(ns,dtype=int);exact=np.zeros(ns,dtype=int);counts=np.zeros((ns,ns),dtype=int)
     for s in range(lo,hi+1):
      ids=np.flatnonzero((segments==s)&valid);times=local[ids];ends=times+lag;idx=np.searchsorted(times,ends);exists=idx<len(ids);start=ids[exists];end=ids[idx[exists]];ok=np.abs(local[end]-local[start]-lag)<1e-7;start=start[ok];end=end[ok]
      ev=events[k,s];keep=np.searchsorted(ev,local[end],side='right')==np.searchsorted(ev,local[start],side='left');both=np.isfinite(H[start])&np.isfinite(H[end]);np.add.at(before,labels[start],1);np.add.at(retained,labels[start[keep]],1);np.add.at(exact,labels[start[keep&both]],1);np.add.at(counts,(labels[start[keep]],labels[end[keep]]),1)
     censor.append({'slot':k,'regime':regime,'state_count':ns,'lag_ps':lag,'before':before.tolist(),'retained_no_exchange':retained.tolist(),'both_exact_retained':exact.tolist(),'state_distribution_js_nats':js(before,retained),'offdiagonal_retained':int(counts.sum()-np.trace(counts)),'retained_counts':counts.tolist()})
  print('diagnosed',k,flush=True)
 overlap=[]
 for regime in range(3):
  for ns in [12,36]:
   for k in range(1,48):
    a1=histograms[k,regime,ns];b=histograms[k+1,regime,ns];overlap.append({'regime':regime,'states':ns,'slots':[k,k+1],'js_nats':js(a1,b),'histogram_overlap':float(np.minimum(a1/a1.sum(),b/b.sum()).sum())})
 result={'status':'CONDITIONS_ASSESSED_PHYSICAL_DHAM_NO_GO','frames':total_frames,'valid_frames':sum(r['valid_frames'] for r in manifest['records']),'state_definition':'CA COM MIC distance and maximum CA Rg, fixed 12/36 geometric bins; chain permutation invariant; not a Markov state model','oracles':oracle,'records':records,'adjacent_temperature_conformation_overlap':overlap,'state_mean_bias_checks':bias,'censoring':censor,'decisions':{'physical_dham_rates':'NO_GO: no correction/likelihood for state-dependent exchange-survival selection; no lag plateau or validated kinetic state model','legacy_mean_energy_dham':'NO_GO until within-state reduced bias variation and target dynamics assumptions are validated','conformation_equilibrium':'Descriptive within-regime/time-block screens only; one coupled REMD run does not supply independent replicate uncertainty','protocol_pooling':'Disabled: three protocol regimes evaluated separately','equilibrium_alternative':'Per-configuration equilibrium reweighting may be assessed within one numerical protocol after conformation convergence and correlated-sample uncertainty checks; no weights or FES fit here','new_md_performed':False},'features_manifest_sha256':sha(a.features/'manifest.json'),'classification_sha256':sha(root/'test/tmp/dimer_protocol_classification_20260909/classification.json'),'exchange_witnesses_sha256':event_hashes,'script_sha256':sha(__file__)}
 (w/'conditions.json').write_text(json.dumps(result,indent=2)+'\n')
if __name__=='__main__':main()
