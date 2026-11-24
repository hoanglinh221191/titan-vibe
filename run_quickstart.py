"""Run the quick-start examples and compare them with the recorded results.

  python run_quickstart.py --output quickstart_runs

Every example folder under examples/quickstart is copied into a new output directory and run there
(the package is not written to). For every run the frame-report fields in EXPECTED.json are compared
exactly, frame by frame, and the centroid of every molecule in the output PDB is compared within
--tolerance (default 0.01 A): coordinates are printed to 0.001 A, so builds may differ in the last
digit, while a molecule written at another lattice image moves by a cell length.
"""
import argparse
import csv
import json
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path


def default_binary():
    """titan-vibe.exe in this package on Windows; source/titan-vibe (the Linux build) elsewhere."""
    root = Path(__file__).resolve().parent
    names = ['titan-vibe.exe'] if os.name == 'nt' else ['titan-vibe']
    candidates = [root / n for n in names] + [root / 'source' / n for n in names]
    if os.name != 'nt':
        candidates.reverse()
    return next((c for c in candidates if c.is_file()), candidates[0])


def read_models(pdb):
    """Atom coordinates of every model of a multi-model PDB."""
    models, cur = [], []
    with open(pdb, encoding="ascii", errors="replace") as fh:
        for line in fh:
            if line.startswith(("ATOM", "HETATM")):
                cur.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
            elif line.startswith("ENDMDL") and cur:
                models.append(cur)
                cur = []
    if cur:
        models.append(cur)
    return models


def centroids(model, ranges):
    out = []
    for a, b in ranges:
        n = b - a
        out.append([sum(model[i][k] for i in range(a, b)) / n for k in range(3)])
    return out


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--output", required=True, type=Path, help="New directory; an existing one is refused")
    p.add_argument("--binary", type=Path, default=default_binary(),
                   help="Default: titan-vibe.exe (Windows) or source/titan-vibe (Linux) in this package")
    p.add_argument("--tolerance", type=float, default=0.01, help="Centroid tolerance in Angstrom")
    a = p.parse_args()
    root = Path(__file__).resolve().parent
    qs = root / "examples" / "quickstart"
    out = a.output.resolve()
    out.mkdir(parents=True, exist_ok=False)
    ff = (root / "forcefields").as_posix()
    failures = 0
    expected_files = sorted(qs.rglob("EXPECTED.json"))
    if not expected_files:
        print("RESULT: FAIL (no EXPECTED.json under examples/quickstart)")
        sys.exit(1)
    for expected_path in expected_files:
        example = expected_path.parent
        rel = example.relative_to(qs)
        expected = json.loads(expected_path.read_text(encoding="utf-8"))
        work = out / rel
        shutil.copytree(example, work)
        for arg, exp in expected["runs"].items():
            text = (work / arg).read_text()
            for depth in ("../../../../forcefields", "../../../forcefields"):
                text = text.replace(depth, ff)
            (work / arg).write_text(text)
            t0 = time.perf_counter()
            with open(work / (arg + ".stdout.txt"), "wb") as so, open(work / (arg + ".stderr.txt"), "wb") as se:
                run = subprocess.run([str(a.binary.resolve()), "-a", arg], cwd=work, stdout=so, stderr=se)
            secs = time.perf_counter() - t0
            problems = []
            if run.returncode != 0:
                problems.append(f"exit code {run.returncode}")
            report = work / "pbctopo_frame_report.csv"
            rows = list(csv.DictReader(open(report, newline=""))) if report.is_file() else []
            if len(rows) != len(exp["frames"]):
                problems.append(f"{len(rows)} report rows, expected {len(exp['frames'])}")
            for row, want in zip(rows, exp["frames"]):
                for k, v in want.items():
                    if row.get(k) != v:
                        problems.append(f"frame {want.get('frame')}: {k} = {row.get(k)!r}, expected {v!r}")
            pdb = work / exp["output_pdb"]
            models = read_models(pdb) if pdb.is_file() else []
            if len(models) != len(exp["centroids_A"]):
                problems.append(f"{len(models)} models in {exp['output_pdb']}, expected {len(exp['centroids_A'])}")
            worst = 0.0
            for model, want in zip(models, exp["centroids_A"]):
                got = centroids(model, exp["molecule_ranges"])
                worst = max([worst] + [abs(g - w) for gc, wc in zip(got, want) for g, w in zip(gc, wc)])
            if worst > a.tolerance:
                problems.append(f"molecule centroid differs by {worst:.4f} A")
            # keep every run's outputs when an example has several argument files
            keep = work / ("out_" + Path(arg).stem)
            keep.mkdir()
            for f in list(work.glob("pbctopo_*")) + [pdb, work / "log.txt", work / "cal.txt"]:
                if f.is_file():
                    shutil.move(str(f), str(keep / f.name))
            status = "PASS" if not problems else "FAIL"
            failures += bool(problems)
            print(f"{status}  {rel.as_posix():<32} {arg:<15} {len(rows):>3} frames  {secs:6.1f} s  "
                  f"max centroid diff {worst:.4f} A")
            for msg in problems[:10]:
                print("      " + msg)
    print("RESULT:", "PASS" if failures == 0 else f"FAIL ({failures} runs)")
    sys.exit(1 if failures else 0)


if __name__ == "__main__":
    main()
