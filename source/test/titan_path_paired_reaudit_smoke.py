"""No-MD recovery safety gates: path isolation, identity, failure reporting."""
import hashlib
import io
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'tools'))
import titan_path_paired_reaudit as r


class Safety(unittest.TestCase):
    def test_archive_names(self):
        for name in ('../x','/x','a/../../x','C:/x','a\\x','a//x','./x',''):
            with self.subTest(name=name),self.assertRaises(ValueError):r.safe_name(name)
        self.assertEqual(str(r.safe_name('generated/a.json')),'generated/a.json')

    def test_md_and_build_commands_forbidden(self):
        with patch.object(r.subprocess,'run',side_effect=AssertionError('must not execute')):
            for args in (['mdrun'],['grompp'],['smd'],['sbatch'],['bash']):
                with self.assertRaisesRegex(ValueError,'NO_MD_COMMAND'):r.run_native('/gmx',args,Path('/tmp'),'test')

    def test_bad_archive_hash_no_extraction(self):
        with tempfile.TemporaryDirectory() as temp:
            base=Path(temp);p=base/'bad.tar.gz';p.write_bytes(b'not archive')
            with self.assertRaisesRegex(ValueError,'FAILED_ARCHIVE_HASH'):r.extract_verified(p,base/'retained')
            self.assertFalse((base/'retained').exists())

    def test_runtime_exclusive(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp);(root/'generated').mkdir()
            with patch.object(r,'extract_verified',side_effect=AssertionError('must not extract')),self.assertRaises(FileExistsError):
                r.run(root,Path('/input'),Path('/gmx'))

    def test_failed_marker_is_not_pass(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp);(root/'downloads').mkdir();(root/'downloads/FAILED.json').write_text('{"status":"FAILED"}')
            with patch('sys.stdout',new_callable=io.StringIO) as out,self.assertRaisesRegex(ValueError,'REAUDIT_NOT_PASS'):r.audit(root)
            self.assertIn('FAILED',out.getvalue())

    def test_static_hash_and_extra_file(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp);(root/'a.py').write_bytes(b'original')
            (root/'STATIC.json').write_text(json.dumps({'a.py':hashlib.sha256(b'original').hexdigest()}));(root/'STATIC.json.sha256').write_text('outer lock')
            r.verify_package(root)
            (root/'a.py').write_bytes(b'changed')
            with self.assertRaisesRegex(ValueError,'STATIC_HASH'):r.verify_package(root)
            (root/'extra').write_text('extra')
            with self.assertRaisesRegex(ValueError,'STATIC_EXACT_SET'):r.verify_package(root)


if __name__=='__main__':unittest.main()
