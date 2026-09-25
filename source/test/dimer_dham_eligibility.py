"""Verify audit artifacts and quantify DHAM data eligibility without fitting.

Predeclared lags: 50, 100, 200, 500 ps. Closed-interval no-exchange selection
is reported as a conservative availability count, never an unbiased estimator.
"""
import argparse,bisect,collections,gzip,hashlib,json,math
from pathlib import Path
from dimer_segment_alignment import sha,closed_no_exchange

LAGS=(50,100,200,500)
PV_KJ_PER_BAR_NM3=0.0602214076

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--audit-roots',nargs='+',type=Path,required=True);p.add_argument('--events-root',type=Path,required=True);p.add_argument('--work-root',type=Path,required=True);a=p.parse_args();w=a.work_root.resolve();w.mkdir(parents=True,exist_ok=False)
 slots={};inputs={}
 for folder in a.audit_roots:
  sm=folder/'summary.json';summary=json.loads(sm.read_text());assert summary['status']=='PASS_IDENTITY_AND_CLOCK'
  inputs[str(sm)]=sha(sm)
  for rec in summary['records']:
   k=rec['slot'];assert k not in slots;path=folder/f'a{k}'
   assert sha(path/'audit.json')==rec['audit_sha256'];slots[k]=path
 assert sorted(slots)==list(range(1,49))
 # Only previously verified global histories are used. This is not a re-audit
 # of every replica log; their provenance is explicitly retained.
 changes={}
 for seg in range(1,21):
  path=a.events_root/f'{seg}_exchange_events.json';events=json.loads(path.read_text());inputs[str(path)]=sha(path)
  for k in slots:changes[k,seg]=[t for _,t,pairs in events if any(k-1 in pair for pair in pairs)]
 totals=collections.Counter();slot_records=[];lag_totals={lag:collections.Counter() for lag in LAGS};all_segments=[]
 with gzip.open(w/'exact_npt_samples.jsonl.gz','wt',encoding='utf-8') as output:
  for k,path in sorted(slots.items()):
   audit=json.loads((path/'audit.json').read_text());alignment=path/'alignment.jsonl.gz';assert sha(alignment)==audit['alignment_sha256'];inputs[str(path/'audit.json')]=sha(path/'audit.json');inputs[str(alignment)]=sha(alignment)
   source={};temperatures=[];offsets=[]
   for s in audit['segments']:
    seg=s['segment'];enfile=path/f'{seg}.xvg';assert sha(enfile)==s['extracted_xvg_sha256']
    # Independent numeric reader for exact alignment verification.
    erows=[tuple(map(float,line.split())) for line in enfile.read_text().splitlines() if line.strip() and not line.startswith(('#','@'))]
    source[seg]={r[0]:r[1:] for r in erows};params=s['parameters'];temps=params['ref_t_K'];assert len(set(temps))==1
    assert params['ref_p_matrix_bar']==[[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]]
    assert params['tcoupl']=='V-rescale' and params['pcoupl']=='C-rescale' and params['free-energy']=='no'
    temperatures.append(temps[0]);offsets.append(s['clock_energy_minus_xtc_ps'])
    all_segments.append({'slot':k,'segment':seg,'energy_offset_ps':s['energy_global_offset_ps'],'energy_minus_xtc_ps':s['clock_energy_minus_xtc_ps'],'omitted_source_t0':s['discarded_source_start'] is not None,'ref_temperature_K':temps[0]})
   assert len(set(temperatures))==1
   groups=collections.defaultdict(list);count=collections.Counter();max_volume=0
   with gzip.open(alignment,'rt',encoding='utf-8') as f:
    for idx,line in enumerate(f):
     row=json.loads(line);assert row['slot_one_based']==k and row['frame_zero_based']==idx
     seg=row['segment'];local=row['local_time_ps'];assert abs(local-row['local_step']*.002)<1e-7
     exact=source[seg].get(local);assert bool(exact)==(row['alignment_status']=='EXACT_LOCAL_TIME')
     if exact:
      assert (row['potential_kj_mol'],row['temperature_observed_K'],row['volume_nm3'])==exact
      assert row['energy_cat_row_zero_based'] is not None and row['volume_match']
      assert abs(row['energy_global_time_ps']-row['xtc_global_time_ps']-offsets[seg-1])<1e-6
      max_volume=max(max_volume,row['volume_absolute_error_nm3']);count['exact_npt_samples']+=1
      sample={'slot':k,'segment':seg,'frame_zero_based':idx,'local_step':row['local_step'],'local_time_ps':local,'ref_temperature_K':temperatures[0],'potential_kj_mol':exact[0],'volume_nm3':exact[2],'pressure_bar':1.0,'configurational_enthalpy_kj_mol':exact[0]+PV_KJ_PER_BAR_NM3*exact[2],'xtc_byte_offset':row['xtc_byte_offset'],'physical_transition_certified':False}
      output.write(json.dumps(sample,separators=(',',':'))+'\n')
     else:assert row['potential_kj_mol'] is None and row['energy_global_time_ps'] is None;count['missing_exact_energy']+=1
     groups[seg].append((local,bool(exact)));count['frames']+=1
   assert count['frames']==21240 and len(groups)==20
   lag_record={}
   for lag in LAGS:
    stat=collections.Counter()
    for seg,rows in groups.items():
     bytime=dict(rows);changed=changes[k,seg]
     for start,exact in rows:
      end=start+lag
      if end not in bytime:continue
      stat['within_segment_windows']+=1
      good=closed_no_exchange(changed,start,end)
      both=exact and bytime[end]
      stat['closed_no_exchange']+=good;stat['both_endpoints_exact_energy']+=both;stat['closed_no_exchange_and_both_exact']+=good and both
    lag_record[str(lag)]=dict(stat);lag_totals[lag].update(stat)
   assert lag_record['50']['closed_no_exchange']==audit['counts']['closed_50ps_exchange_free_candidates']
   totals.update(count);slot_records.append({'slot':k,'ref_temperature_K':temperatures[0],'counts':dict(count),'max_volume_difference_nm3':max_volume,'lags_ps':lag_record})
   print('VERIFIED',k,dict(count),flush=True)
 temps=[r['ref_temperature_K'] for r in slot_records];assert all(b>a for a,b in zip(temps,temps[1:]))
 result={'status':'TECHNICAL_ALIGNMENT_PASS_PHYSICAL_DHAM_NOT_QUALIFIED','counts':dict(totals),'energy_slots_verified':48,'energy_segments_verified':len(all_segments),'discarded_source_start_count':sum(s['omitted_source_t0'] for s in all_segments),'temperature_range_K':[temps[0],temps[-1]],'lags_ps':{str(lag):dict(stat) for lag,stat in lag_totals.items()},'slot_records':slot_records,'segment_records':all_segments,'endpoint_energy_gate_note':'Requiring exact energy at both window endpoints is a conservative availability diagnostic, not a universal DHAM requirement. State-mean-energy methods can use a subset, but must validate within-state sampling and bias approximation.','decision':{'exact_equilibrium_sample_table':'Available as candidate input; equilibration, common Hamiltonian and statistical overlap/ESS still require checks.','current_titan_dham_direct_input':'NOT_QUALIFIED: flat per-slot transition counting lacks exchange history and explicit segment IDs; nearest-energy alignment accepts noncoincident samples.','physical_rates':'NOT_QUALIFIED: exchange discontinuities, state-dependent censoring, lag/Markov/state-bias validity untested.','fit_performed':False},'npt':{'pressure_bar':1.0,'observable':'U + P*V','conversion_kj_per_bar_nm3':PV_KJ_PER_BAR_NM3,'note':'Use reference temperatures, not instantaneous kinetic temperatures. TITAN already has an NPT mode; this table does not invoke its estimator.'},'inputs_sha256':inputs,'exact_npt_samples_sha256':sha(w/'exact_npt_samples.jsonl.gz'),'script_sha256':sha(__file__)}
 (w/'eligibility.json').write_text(json.dumps(result,indent=2)+'\n')
 print(json.dumps({'status':result['status'],'counts':result['counts'],'lags':result['lags_ps']}),flush=True)

if __name__=='__main__':main()
