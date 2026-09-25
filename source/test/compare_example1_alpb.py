"""ALPB component gates with explicit EEL convention and previously qualified SASA."""
import argparse,hashlib,json,math,subprocess
from pathlib import Path
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--reference',type=Path,required=True);p.add_argument('--base-reference',type=Path,required=True);p.add_argument('--base-comparison',type=Path,required=True);p.add_argument('--driver',type=Path,required=True);p.add_argument('--work-root',type=Path,required=True);a=p.parse_args();w=a.work_root.resolve();w.mkdir(parents=True,exist_ok=False)
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
ref=json.loads(a.reference.read_text());base=json.loads(a.base_reference.read_text());previous=json.loads(a.base_comparison.read_text());assert previous['status']=='PASS' and previous['reference_sha256']==sha(a.base_reference)
records=[]
for model in [0,25,51]:
 br=[r for r in base['records'] if r['model']==model];fixture=a.base_reference.parent/br[0]['fixture']
 r=subprocess.run([str(a.driver.resolve()),str(fixture.resolve())],capture_output=True,text=True);(w/f'model_{model}.log').write_text(r.stdout+r.stderr);r.check_returncode()
 lines=[list(map(float,line.split())) for line in r.stdout.splitlines() if len(line.split())==7];assert len(lines)==2
 for values in lines:
  salt,ee,vdw,gpol,ar,al,ac=values;rr=next(r for r in ref['records'] if r['model']==model and r['salt']==salt);old=next(r for r in previous['records'] if r['model']==model and r['salt']==salt and r['points']==6144);bb=next(r for r in br if r['salt']==salt)
  # ALPB changes only the polar term. Reuse the already qualified 6144-point
  # SASA result with its reference/input hashes; the 32-point driver SASA is
  # deliberately not part of this total-energy comparison.
  np_titan=next(t['titan'] for t in old['terms'] if t['component']=='nonpolar');np_ref=bb['delta']['nonpolar']
  # TITAN now applies the same 1/(1+alpha*beta) factor to the direct Coulomb
  # term that sander applies for alpb=1, so the comparison is against the raw
  # sander EEL. The vacuum-convention column stays in the reference for audit.
  actual={'EEL':ee,'VDWAALS':vdw,'EGB':gpol};terms=[]
  for name,value in actual.items():
   expected=rr['delta_kj_mol'][name];tol=.01+1e-4*abs(expected);error=abs(value-expected)
   terms.append({'component':name,'titan':value,'amber':expected,'absolute_error':error,'tolerance':tol,'status':'PASS' if math.isfinite(value) and error<=tol else 'FAIL'})
  sizes=[]
  for name,value in zip(['receptor','ligand','complex'],[ar,al,ac]):
   expected=rr['arad_A'][name];error=abs(value-expected);sizes.append({'endpoint':name,'titan_A':value,'elsize_A':expected,'absolute_error_A':error,'tolerance_A':.000501,'status':'PASS' if error<=.000501 else 'FAIL'})
  total_titan=ee+vdw+gpol+np_titan;total_ref=sum(rr['delta_kj_mol'][n] for n in ['EEL_vacuum_convention','VDWAALS','EGB'])+np_ref
  sasa_tol=max(5,.01*abs(bb['delta']['sasa_a2']))*.0072*4.184;total_tol=.01+1e-4*abs(total_ref)+sasa_tol
  total={'titan':total_titan,'amber_vacuum_convention':total_ref,'amber_raw_alpb_total':sum(rr['delta_kj_mol'][n] for n in ['EEL','VDWAALS','EGB'])+np_ref,'absolute_error':abs(total_titan-total_ref),'tolerance':total_tol,'status':'PASS' if abs(total_titan-total_ref)<=total_tol else 'FAIL'}
  records.append({'model':model,'salt':salt,'terms':terms,'sizes':sizes,'total_with_previously_qualified_sasa':total});print(model,salt,[(t['component'],t['status'],t['absolute_error']) for t in terms],total['status'],flush=True)
status='PASS' if all(t['status']=='PASS' for r in records for t in r['terms']+r['sizes']+[r['total_with_previously_qualified_sasa']]) else 'FAIL'
(w/'comparison.json').write_text(json.dumps({'status':status,'scope':'Reconstructed topology; independent ALPB polar and A_det gates. Vacuum-MM convention explicitly aligned; raw Sander ALPB totals are retained but are not asserted equal. SASA reused from qualified 6144-point calculations.','records':records,'reference_sha256':sha(a.reference),'base_reference_sha256':sha(a.base_reference),'base_comparison_sha256':sha(a.base_comparison),'driver_sha256':sha(a.driver),'script_sha256':sha(__file__)},indent=2)+'\n');raise SystemExit(status!='PASS')
