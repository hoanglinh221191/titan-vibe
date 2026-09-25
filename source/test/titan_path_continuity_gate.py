"""Real TITAN CLI no-MD regression for explicit route state continuity."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess

ROOT=Path(__file__).resolve().parents[1]
ARG=ROOT/'test/titan_path_g3_example1/arg_nav_p0_nt1.in'

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--titan',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    titan=args.titan.resolve(strict=True)
    output=args.output.resolve()
    output.mkdir(parents=True,exist_ok=False)
    cases=[]
    def run(name,extra=(),expected=None):
        work=output/name
        work.mkdir()
        command=[str(titan),'-a',str(ARG),*map(str,extra)]
        r=subprocess.run(command,cwd=work,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,encoding='utf-8',errors='replace')
        (work/'console.txt').write_text(r.stdout,encoding='utf-8')
        if expected is None:
            if r.returncode!=0: raise RuntimeError(r.stdout[-2000:])
            raw=(work/'pull_direction_status.txt').read_text()
            status=dict(line.split('=',1) for line in raw.splitlines() if '=' in line)
            if status.get('STATUS')!='CONTINUE': raise RuntimeError('unexpected terminal '+str(status))
        else:
            if r.returncode==0 or expected not in r.stdout:
                raise RuntimeError('expected rejection '+expected+'\n'+r.stdout[-2000:])
            if (work/'pull_direction_status.txt').exists():
                raise RuntimeError('uncertain projection emitted a decision')
            status=None
        cases.append({'name':name,'exit_code':r.returncode,'expected_diagnostic':expected,
                      'console_sha256':sha(work/'console.txt')})
        return work,status
    first,state=run('initial')
    previous=first/'pull_direction_status.txt'
    _,bounded=run('bounded',['-pathway-previous-status',previous,'-pathway-continuity-window','1'])
    assert bounded['PROJECTION_CONTINUITY']=='BOUNDED_PREVIOUS_PROGRESS'
    assert abs(float(bounded['PATH_PROGRESS_A'])-float(state['PATH_PROGRESS_A']))<1e-8
    assert bounded['MODE']==state['MODE']
    raw=previous.read_text()
    variants={
        'wrong_identity':(raw.replace(state['ROUTE_ID'],'UNTRUSTED'), 'identity/state mismatch'),
        'future_time':(raw.replace('SNAPSHOT_TIME_PS='+state['SNAPSHOT_TIME_PS'],'SNAPSHOT_TIME_PS=1e99'),'progress/length/time mismatch'),
        'duplicate_mode':(raw+'MODE=Recovery\n','duplicate field'),
        'clipped_window':(raw.replace('PATH_PROGRESS_A='+state['PATH_PROGRESS_A'],'PATH_PROGRESS_A=4'),'PROJECTION_CONTINUITY_WINDOW_EXHAUSTED'),
    }
    for name,(content,error) in variants.items():
        p=output/(name+'_previous.txt')
        p.write_text(content)
        run(name,['-pathway-previous-status',p,'-pathway-continuity-window','1'],error)
    # Argument-file path normalization must resolve relative to its file, not cwd.
    lines=[]
    path_flags={'-topo','-topo2','-common_topo','-path-a','-path-n','-f','-n','-pathway-file'}
    for line in ARG.read_text().splitlines():
        key,sep,value=line.partition('=')
        if sep and key in path_flags: line=key+'='+str((ARG.parent/value).resolve())
        lines.append(line)
    relative=output/'relative.in'
    relative.write_text('\n'.join(lines)+'\n-pathway-previous-status=initial/pull_direction_status.txt\n-pathway-continuity-window=1\n')
    _,relative_state=run('relative_arg',['-a',relative])
    assert relative_state['PROJECTION_CONTINUITY']=='BOUNDED_PREVIOUS_PROGRESS'
    forbidden=[str(p) for p in output.rglob('*') if p.suffix in ('.tpr','.cpt','.xtc','.trr','.edr')]
    assert not forbidden
    report={'status':'PASS_NO_MD_CONTINUITY','cases':cases,'titan_sha256':sha(titan),'md_executed':False}
    (output/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print('PASS TITAN_PATH_CONTINUITY_GATE cases='+str(len(cases))+' no_md=true')

if __name__=='__main__': main()
