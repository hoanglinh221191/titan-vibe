#!/usr/bin/env python3
"""Focused fail-closed tests for the VIBE annotated candidate-tag gate."""

from __future__ import annotations

import hashlib
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "tools/vibe_evidence_lock"))

from freeze_candidate import (  # noqa: E402
    audit_release_contract,
    audit_required_paths_against_revision,
    decide_freeze_verdict,
)


MANUAL_PHRASES = [
    "Optional contacts are observational evidence only.",
    "H0 is weak and underdetermined.",
]
VALID_MANUAL = "\n".join(MANUAL_PHRASES) + "\n"
SPEC_PATH = "docs/spec.md"
SPEC_BYTES = b"immutable scientific contract\n"


def run(repo: Path, *arguments: str) -> None:
    subprocess.run(
        ["git", *arguments],
        cwd=repo,
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )


class CandidateTagGateTest(unittest.TestCase):
    def make_repo(self, manual: str = VALID_MANUAL) -> tuple[tempfile.TemporaryDirectory, Path]:
        temporary = tempfile.TemporaryDirectory(prefix="vibe_tag_gate_")
        repo = Path(temporary.name)
        run(repo, "init", "--initial-branch=main")
        run(repo, "config", "user.name", "VIBE Test")
        run(repo, "config", "user.email", "vibe-test@example.invalid")
        (repo / "required.txt").write_bytes(b"required\n")
        (repo / "manual.txt").write_bytes(manual.encode("utf-8"))
        (repo / "docs").mkdir()
        (repo / SPEC_PATH).write_bytes(SPEC_BYTES)
        run(repo, "add", "required.txt", "manual.txt", SPEC_PATH)
        run(repo, "commit", "-m", "candidate")
        return temporary, repo

    def audit(self, repo: Path, tag: str = "vibe-candidate", required=None, spec_hash=None):
        return audit_release_contract(
            repo,
            tag,
            required or ["required.txt"],
            [
                {
                    "path": SPEC_PATH,
                    "sha256": spec_hash or hashlib.sha256(SPEC_BYTES).hexdigest(),
                }
            ],
            "manual.txt",
            MANUAL_PHRASES,
        )

    def test_annotated_tag_at_clean_head_passes(self) -> None:
        temporary, repo = self.make_repo()
        self.addCleanup(temporary.cleanup)
        run(repo, "tag", "-a", "vibe-candidate", "-m", "candidate tag")
        result = self.audit(repo)
        self.assertTrue(result["pass"])
        self.assertTrue(result["annotated_tag_exists"])
        self.assertTrue(result["tag_resolves_to_head"])
        self.assertTrue(result["worktree_clean"])
        self.assertTrue(result["immutable_release_inputs_valid"])
        self.assertTrue(result["manual_semantics_valid"])
        self.assertTrue(result["tag_object_sha"])
        self.assertEqual(result["commit_sha"], result["head_sha"])
        self.assertTrue(result["tree_sha"])
        self.assertTrue(result["tagged_manifest_digest"])

    def test_missing_tag_fails(self) -> None:
        temporary, repo = self.make_repo()
        self.addCleanup(temporary.cleanup)
        result = self.audit(repo)
        self.assertFalse(result["pass"])
        self.assertIn("candidate_tag_missing", result["failure_reasons"])
        self.assertFalse(result["immutable_release_inputs_valid"])

    def test_lightweight_tag_fails(self) -> None:
        temporary, repo = self.make_repo()
        self.addCleanup(temporary.cleanup)
        run(repo, "tag", "vibe-candidate")
        result = self.audit(repo)
        self.assertFalse(result["pass"])
        self.assertIn("candidate_tag_is_not_annotated", result["failure_reasons"])

    def test_wrong_commit_tag_fails(self) -> None:
        temporary, repo = self.make_repo()
        self.addCleanup(temporary.cleanup)
        run(repo, "tag", "-a", "vibe-candidate", "-m", "old candidate")
        (repo / "later.txt").write_bytes(b"later\n")
        run(repo, "add", "later.txt")
        run(repo, "commit", "-m", "later")
        result = self.audit(repo)
        self.assertFalse(result["pass"])
        self.assertIn("candidate_tag_does_not_resolve_to_head", result["failure_reasons"])

    def test_dirty_checkout_fails(self) -> None:
        temporary, repo = self.make_repo()
        self.addCleanup(temporary.cleanup)
        run(repo, "tag", "-a", "vibe-candidate", "-m", "candidate tag")
        (repo / "unrelated.tmp").write_bytes(b"dirty\n")
        result = self.audit(repo)
        self.assertFalse(result["pass"])
        self.assertIn("originating_worktree_is_dirty", result["failure_reasons"])

    def test_clean_crlf_checkout_matches_tag_through_git_clean_filter(self) -> None:
        temporary, repo = self.make_repo()
        self.addCleanup(temporary.cleanup)
        run(repo, "tag", "-a", "vibe-candidate", "-m", "candidate tag")
        run(repo, "config", "core.autocrlf", "true")
        for relative in ("required.txt", "manual.txt", SPEC_PATH):
            (repo / relative).unlink()
        run(repo, "checkout", "--", ".")
        raw_manual = (repo / "manual.txt").read_bytes()
        tag_manual = subprocess.run(
            ["git", "show", "vibe-candidate:manual.txt"],
            cwd=repo,
            check=True,
            stdout=subprocess.PIPE,
        ).stdout
        self.assertIn(b"\r\n", raw_manual)
        self.assertNotEqual(raw_manual, tag_manual)
        result = self.audit(repo)
        self.assertTrue(result["pass"])
        self.assertTrue(result["required_worktree_git_clean_content_matches_tag"])
        head_audit = audit_required_paths_against_revision(
            repo, "HEAD", ["required.txt", "manual.txt", SPEC_PATH]
        )
        self.assertTrue(head_audit["reconstructible_from_revision"])
        self.assertEqual(head_audit["uncommitted_required_files"], [])

    def test_full_release_inventory_blocks_missing_manual_and_spec(self) -> None:
        temporary = tempfile.TemporaryDirectory(prefix="vibe_release_inventory_")
        self.addCleanup(temporary.cleanup)
        repo = Path(temporary.name)
        run(repo, "init", "--initial-branch=main")
        run(repo, "config", "user.name", "VIBE Test")
        run(repo, "config", "user.email", "vibe-test@example.invalid")
        (repo / "overlay.txt").write_bytes(b"committed overlay\n")
        run(repo, "add", "overlay.txt")
        run(repo, "commit", "-m", "overlay only")
        (repo / "manual.txt").write_text(VALID_MANUAL, encoding="utf-8")
        (repo / "docs").mkdir()
        (repo / SPEC_PATH).write_bytes(SPEC_BYTES)

        overlay_audit = audit_required_paths_against_revision(
            repo, "HEAD", ["overlay.txt"]
        )
        release_audit = audit_required_paths_against_revision(
            repo, "HEAD", ["overlay.txt", "manual.txt", SPEC_PATH]
        )
        self.assertTrue(overlay_audit["reconstructible_from_revision"])
        self.assertFalse(release_audit["reconstructible_from_revision"])
        self.assertEqual(
            release_audit["unversioned_required_files"],
            [SPEC_PATH, "manual.txt"],
        )
        self.assertEqual(
            decide_freeze_verdict(
                "preparation",
                True,
                False,
                len(release_audit["uncommitted_required_files"]),
            ),
            "BLOCKED_PENDING_VERSIONING",
        )

    def test_required_file_mismatch_fails(self) -> None:
        temporary, repo = self.make_repo()
        self.addCleanup(temporary.cleanup)
        run(repo, "tag", "-a", "vibe-candidate", "-m", "candidate tag")
        (repo / "required.txt").write_bytes(b"mutated\n")
        result = self.audit(repo)
        self.assertFalse(result["pass"])
        self.assertEqual(result["worktree_git_clean_mismatches"], ["required.txt"])

    def test_overlay_only_required_file_fails(self) -> None:
        temporary, repo = self.make_repo()
        self.addCleanup(temporary.cleanup)
        run(repo, "tag", "-a", "vibe-candidate", "-m", "candidate tag")
        (repo / "overlay-only.txt").write_bytes(b"overlay\n")
        result = self.audit(repo, required=["required.txt", "overlay-only.txt"])
        self.assertFalse(result["pass"])
        self.assertEqual(result["missing_from_tag"], ["overlay-only.txt"])

    def test_wrong_immutable_spec_hash_fails(self) -> None:
        temporary, repo = self.make_repo()
        self.addCleanup(temporary.cleanup)
        run(repo, "tag", "-a", "vibe-candidate", "-m", "candidate tag")
        result = self.audit(repo, spec_hash="0" * 64)
        self.assertFalse(result["pass"])
        self.assertFalse(result["immutable_release_inputs_valid"])
        self.assertNotEqual(
            decide_freeze_verdict("release", True, result["pass"], 0), "READY"
        )

    def test_immutable_spec_missing_from_tag_fails_with_both_reasons(self) -> None:
        temporary, repo = self.make_repo()
        self.addCleanup(temporary.cleanup)
        run(repo, "tag", "-a", "vibe-candidate", "-m", "candidate tag")
        missing_path = "docs/missing-spec.md"
        result = audit_release_contract(
            repo,
            "vibe-candidate",
            ["required.txt"],
            [{"path": missing_path, "sha256": "0" * 64}],
            "manual.txt",
            MANUAL_PHRASES,
        )
        self.assertFalse(result["pass"])
        self.assertIn(missing_path, result["missing_from_tag"])
        self.assertIn(
            "required_paths_missing_from_candidate_tag", result["failure_reasons"]
        )
        self.assertIn(
            "immutable_release_input_hash_mismatch", result["failure_reasons"]
        )

    def test_missing_manual_semantics_fails(self) -> None:
        temporary, repo = self.make_repo(manual="PBCTOPO options only\n")
        self.addCleanup(temporary.cleanup)
        run(repo, "tag", "-a", "vibe-candidate", "-m", "candidate tag")
        result = self.audit(repo)
        self.assertFalse(result["pass"])
        self.assertFalse(result["manual_semantics_valid"])
        self.assertEqual(result["manual_semantics_missing_phrases"], MANUAL_PHRASES)
        self.assertNotEqual(
            decide_freeze_verdict("release", True, result["pass"], 0), "READY"
        )

    def test_verdict_boundary_is_release_only_and_fail_closed(self) -> None:
        self.assertEqual(
            decide_freeze_verdict("preparation", True, True, 0),
            "BLOCKED_MISSING_TAG",
        )
        self.assertEqual(
            decide_freeze_verdict("preparation", True, True, 3),
            "BLOCKED_PENDING_VERSIONING",
        )
        self.assertEqual(
            decide_freeze_verdict("release", True, False, 0),
            "BLOCKED_TAG_INTEGRITY",
        )
        self.assertEqual(
            decide_freeze_verdict("release", False, True, 0),
            "BLOCKED_TECHNICAL",
        )
        self.assertEqual(
            decide_freeze_verdict("release", True, True, 0),
            "READY",
        )


if __name__ == "__main__":
    unittest.main(verbosity=2)
