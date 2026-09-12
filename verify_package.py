from pathlib import Path
import hashlib
root=Path(__file__).resolve().parent
count=0
for line in (root/'SHA256SUMS.txt').read_text().splitlines():
    expected,name=line.split('  ',1)
    p=root/name
    if not p.is_file():raise SystemExit('Missing file: '+name)
    h=hashlib.sha256()
    with p.open('rb') as f:
        for data in iter(lambda:f.read(4194304),b''):h.update(data)
    if h.hexdigest()!=expected:raise SystemExit('Checksum mismatch: '+name)
    count+=1
print('PASS:',count,'files match the package manifest')
