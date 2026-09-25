"""Synthetic-only tests of immutable equilibration subset collection."""
from pathlib import Path
import contextlib
import io
import sys
import tarfile
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/"tools"))
import titan_path_nav6c_frame_inputs as collector


class CollectionSmoke(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix="nav6c-collector-SYNTHETIC-")
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)/"pilot"
        self.root.mkdir()
        self.output = Path(self.temp.name)/"collected"
        self.expected = {}
        for name in collector.INPUTS:
            path = self.root/name; path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(("SYNTHETIC NOT MD: "+name).encode())
            self.expected[name] = {"path": name, "bytes": str(path.stat().st_size), "sha256": collector.sha(path)}
        manifest = self.root/collector.MANIFEST; manifest.parent.mkdir()
        manifest.write_text("SYNTHETIC MANIFEST, NOT REAL EVIDENCE\n")
        self.manifest_hash = collector.sha(manifest)

    def collect_synthetic(self):
        with contextlib.redirect_stdout(io.StringIO()):
            return collector._collect(self.root, self.output, self.expected, self.manifest_hash)

    def test_exact_archive_and_unchanged_inputs(self):
        result = self.collect_synthetic()
        self.assertFalse(result["md_executed"])
        self.assertFalse(result["job_submitted"])
        with tarfile.open(self.output/collector.BUNDLE) as archive:
            self.assertEqual(set(archive.getnames()), set(collector.INPUTS)|{collector.MANIFEST})
        collector.verify_inputs(self.root, self.expected)

    def test_public_entry_rejects_synthetic_manifest(self):
        with self.assertRaisesRegex(ValueError, "manifest hash mismatch"):
            collector.collect(self.root, self.output)
        self.assertFalse(self.output.exists())

    def test_missing_input_rejected_before_output(self):
        (self.root/collector.INPUTS[0]).unlink()
        with self.assertRaisesRegex(ValueError, "Missing or changed"):
            self.collect_synthetic()
        self.assertFalse(self.output.exists())

    def test_corruption_rejected_before_output(self):
        (self.root/collector.INPUTS[0]).write_text("WRONG")
        with self.assertRaisesRegex(ValueError, "Missing or changed"):
            self.collect_synthetic()
        self.assertFalse(self.output.exists())

    def test_output_inside_original_forbidden(self):
        self.output = self.root/"new"
        with self.assertRaisesRegex(ValueError, "disjoint"):
            self.collect_synthetic()

    def test_no_overwrite_or_silent_retry(self):
        self.collect_synthetic()
        before = collector.sha(self.output/collector.BUNDLE)
        with self.assertRaisesRegex(ValueError, "Output exists"):
            self.collect_synthetic()
        self.assertEqual(before, collector.sha(self.output/collector.BUNDLE))

    def test_path_traversal_forbidden(self):
        with self.assertRaisesRegex(ValueError, "Unsafe"):
            collector.member(self.root, "../outside")

    def test_drive_path_forbidden(self):
        with self.assertRaisesRegex(ValueError, "Unsafe"):
            collector.member(self.root, "C:/outside")

    def test_changed_manifest_retains_partial_not_ready(self):
        self.manifest_hash = "0"*64
        with self.assertRaisesRegex(ValueError, "Archive payload hash"):
            self.collect_synthetic()
        self.assertTrue((self.output/(collector.BUNDLE+".partial")).is_file())
        self.assertFalse((self.output/collector.BUNDLE).exists())


if __name__ == "__main__":
    unittest.main()
