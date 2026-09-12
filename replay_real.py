"""Run the five systems directly from the 11 September deposited data archive."""
from pathlib import Path
import argparse,csv,json,subprocess,time,hashlib

FIELDS=('status','mode','component_count','evidence_state','fallback_reason','hard_feasible','hard_certified','evidence_graph_connected','soft_observed_contacts_checked','soft_observed_selected_compatible_contact_count','soft_observed_hypothesis_count','soft_observed_compatible_hypothesis_count','soft_observed_ambiguous_relation_count','soft_observed_no_support_relation_count','soft_observed_alternative_hypothesis_contact_count','feasible_assignment_set_enumerated','feasible_assignments_within_domain','wrapped_coordinate_hash')
def coordinates(path,first_only):
    rows=[]
    for line in path.read_text().splitlines():
        if line.startswith(('ATOM  ','HETATM')):rows.append(line[:54])
        elif line.startswith('ENDMDL') and rows and first_only:break
    return rows

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--data-root',required=True,type=Path,help='Extracted vibe_data_update_20260911 directory containing INPUT_PATH_MAP.json')
    p.add_argument('--output',required=True,type=Path,help='New directory; existing directories are refused')
    p.add_argument('--all-frames',action='store_true',help='All 367 comparison frames; default first frame of each system')
    p.add_argument('--only-extras',action='store_true',help='Only the two additional Figure 1 frames')
    p.add_argument('--binary',type=Path,default=Path(__file__).resolve().parent/'titan.exe')
    a=p.parse_args();root=a.data_root.resolve();out=a.output.resolve();out.mkdir(parents=True,exist_ok=False)
    bindings=json.loads((root/'INPUT_PATH_MAP.json').read_text());records={}
    data=root/'jcim_revision_controls_20260911'
    jobs=[]
    for sid in sorted(bindings):
        if not a.only_extras:jobs.append((sid,'A'))
        if a.all_frames or a.only_extras:
            jobs.extend((sid,d.name) for d in sorted((data/'A_vibe'/sid).glob('t*')) if d.is_dir())
    for sid,tag in jobs:
        b=bindings[sid]
        archived=data/'A_vibe'/sid/tag;target=out/sid/tag;target.mkdir(parents=True)
        first_only=tag=='A' and not a.all_frames
        values={'-s':root/b['pdb'],'-n':root/b['ndx'],'-common_topo':root/b['top'],'-traj':data/'A_inputs'/sid/f'full_{tag}.xtc'}
        for f in values.values():
            if not f.is_file():raise FileNotFoundError(f)
        lines=[]
        for line in (archived/'arg.in').read_text().splitlines():
            k,sep,v=line.partition('=')
            if k in values:v=values[k].as_posix()
            elif k=='-e' and first_only:v='0'
            lines.append(k+sep+v)
        (target/'arg.in').write_text('\n'.join(lines)+'\n')
        t=time.perf_counter()
        with (target/'stdout.txt').open('wb') as so,(target/'stderr.txt').open('wb') as se:
            run=subprocess.run([str(a.binary.resolve()),'-a','arg.in'],cwd=target,stdout=so,stderr=se,timeout=900)
        result={'exit_code':run.returncode,'elapsed_s':time.perf_counter()-t}
        name=next(s.split('=',1)[1] for s in lines if s.startswith('-opdb='))
        if run.returncode==0 and (target/name).exists():
            x=coordinates(target/name,False);y=coordinates(archived/name,first_only)
            result.update({'atom_records':len(x),'coordinates_identical_at_PDB_precision':bool(x) and x==y})
            current=list(csv.DictReader((target/'pbctopo_frame_report.csv').open()))
            expected=list(csv.DictReader((archived/'pbctopo_frame_report.csv').open()))
            if first_only:expected=expected[:1]
            result['frames']=len(current)
            result['scientific_fields_identical']=len(current)==len(expected) and all(x.get(k)==y.get(k) for x,y in zip(current,expected) for k in FIELDS)
        result['pass']=run.returncode==0 and result.get('coordinates_identical_at_PDB_precision',False) and result.get('scientific_fields_identical',False)
        records[sid+'/'+tag]=result
        (out/'REAL_REPLAY_RESULT.json').write_text(json.dumps({'binary_sha256':hashlib.sha256(a.binary.read_bytes()).hexdigest(),'all_frames_requested':a.all_frames,'all_pass':all(v['pass'] for v in records.values()),'systems_completed':len(records),'results':records},indent=2))
        print(sid,tag,result,flush=True)
    return 0 if all(x['pass'] for x in records.values()) else 1
if __name__=='__main__':raise SystemExit(main())
