"""Replay deposited fixed inputs; compare coordinates and scientific report fields."""
import argparse,csv,hashlib,json,subprocess,time
from concurrent.futures import ThreadPoolExecutor,as_completed
from pathlib import Path

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--all',action='store_true',help='Run all 480 held-out inputs; default: 20 stratified smoke cases')
    p.add_argument('--output',required=True,type=Path,help='New output directory; existing directories are refused')
    p.add_argument('--binary',type=Path,default=Path(__file__).resolve().parent/'titan.exe')
    p.add_argument('--jobs',type=int,default=4,help='Independent one-thread processes (default 4)')
    a=p.parse_args();root=Path(__file__).resolve().parent;out=a.output.resolve()
    out.mkdir(parents=True,exist_ok=False)
    cases=json.loads((root/'examples/heldout480/manifest.json').read_text())
    results=[]
    binary_hash=hashlib.sha256(a.binary.read_bytes()).hexdigest()
    def run_case(case):
        target=out/case['case_id'];target.mkdir(parents=True)
        inp=root/'examples/heldout480'/case['case_id']
        text=(inp/'arg.template').read_text().replace('@FF@',(root/'forcefields').as_posix()).replace('@CASE@',inp.as_posix())
        (target/'arg.in').write_text(text)
        t=time.perf_counter()
        with (target/'stdout.txt').open('wb') as so,(target/'stderr.txt').open('wb') as se:
            run=subprocess.run([str(a.binary.resolve()),'-a','arg.in'],cwd=target,stdout=so,stderr=se,timeout=180)
        result={'case_id':case['case_id'],'exit_code':run.returncode,'elapsed_s':time.perf_counter()-t}
        pdb=target/'oracle_v2_pbctopo.pdb';report=target/'pbctopo_frame_report.csv'
        result['coordinate_hash_matches']=pdb.is_file() and hashlib.sha256(pdb.read_bytes()).hexdigest()==case['expected_pdb_sha256']
        if report.is_file():
            current=next(csv.DictReader(report.open()))
            result['field_mismatches']={k:{'expected':v,'actual':current.get(k)} for k,v in case['expected_fields'].items() if current.get(k)!=v}
        else:result['field_mismatches']={'report':'missing'}
        result['pass']=run.returncode==0 and result['coordinate_hash_matches'] and not result['field_mismatches']
        return result
    selected=[c for c in cases if a.all or c['smoke']]
    with ThreadPoolExecutor(max_workers=max(1,a.jobs)) as pool:
        pending=[pool.submit(run_case,c) for c in selected]
        for task in as_completed(pending):
            result=task.result();results.append(result)
            (out/'REPLAY_RESULT.json').write_text(json.dumps({'purpose':'Implementation replay against archived coordinates and scientific report fields; not a new independent truth assessment','binary_sha256':binary_hash,'all_pass':all(x['pass'] for x in results),'completed':len(results),'expected':len(selected),'results':sorted(results,key=lambda r:r['case_id'])},indent=2))
            print(len(results),result['case_id'],'PASS' if result['pass'] else 'FAIL',flush=True)
    return 0 if all(x['pass'] for x in results) else 1

if __name__=='__main__':raise SystemExit(main())
