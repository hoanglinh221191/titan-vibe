"""Geometry policy regression using sphere_patch_surface_audit.cpp.

Usage: python test/sphere_patch_policy_regression.py --after AUDIT_EXE
       --before BASELINE_EXE --fixtures ATOMS_DIRECTORY --output NEW_DIRECTORY
The ligand and its translated copy are optional; sphere/two are generated.
Build both executables with the architecture flags of their engine objects.
"""
import argparse
import csv
import io
import json
import math
import pathlib
import subprocess


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--after", required=True, type=pathlib.Path)
    p.add_argument("--before", required=True, type=pathlib.Path)
    p.add_argument("--fixtures", type=pathlib.Path)
    p.add_argument("--output", required=True, type=pathlib.Path)
    a = p.parse_args()
    a.output.mkdir(parents=True, exist_ok=False)
    fixtures = {}
    for name, content in {
        "sphere": "1\n0 0 0 1.7 1\n",
        "two": "2\n0 0 0 1.7 1\n2.5 0 0 1.7 -1\n",
        "sphere_shift": "1\n100 200 300 1.7 1\n",
    }.items():
        fixtures[name] = a.output / (name + ".atoms")
        fixtures[name].write_text(content)
    if a.fixtures:
        for name in ("lig", "lig_shift"):
            target = a.fixtures / (name + ".atoms")
            if target.exists():
                fixtures[name] = target
    checks, rows, cache = [], [], {}

    def run(version, name, n, mode, surface="vdw"):
        key = (version, name, n, mode, surface)
        if key in cache:
            return cache[key]
        exe = a.after if version == "after" else a.before
        r = subprocess.run([str(exe.resolve()), str(fixtures[name].resolve()),
                            str(n), mode, surface], text=True,
                           capture_output=True, timeout=60, check=True)
        (a.output / ("_".join(map(str, key)) + ".log")).write_text(r.stdout)
        row = list(csv.DictReader(io.StringIO(r.stdout)))[-1]
        cache[key] = row
        rows.append(dict(version=version, fixture=name, **row))
        return row

    def check(name, condition):
        checks.append(dict(name=name, passed=bool(condition)))

    def metrics(row):
        return {k: v for k, v in row.items()
                if k not in ("wall_ms", "mesh_mode", "points")}

    for name in ("sphere", "two", *(["lig"] if "lig" in fixtures else [])):
        for n in (96, 384, 1536):
            auto = run("after", name, n, "auto")
            ico = run("after", name, n, "icosphere")
            check(f"{name}/{n}: auto equals explicit icosphere", metrics(auto) == metrics(ico))
            for mode in ("latlong", "icosphere"):
                if name == "sphere":
                    check(f"{name}/{n}/{mode}: isolated explicit mode preserved",
                          metrics(run("before", name, n, mode)) ==
                          metrics(run("after", name, n, mode)))
                else:
                    check(f"{name}/{n}/{mode}: clipped vertices and centroids exposed",
                          float(run("after", name, n, mode)["burial_depth_max"]) < 1e-10)
        check(f"{name}: SAS auto preserves legacy geometry",
              metrics(run("before", name, 96, "auto", "sas")) ==
              metrics(run("after", name, 96, "auto", "sas")))
    for n in (96, 384, 1536, 65536):
        sphere = run("after", "sphere", n, "auto")
        check(f"sphere/{n}: exact spherical area", abs(float(sphere["area_rel_error"])) < 1e-11)
        check(f"sphere/{n}: exact radius and outward normals",
              abs(float(sphere["sphere_rad_min"]) - 1.7) < 1e-12 and
              abs(float(sphere["sphere_rad_max"]) - 1.7) < 1e-12 and
              float(sphere["sphere_normal_cos_min"]) > 1 - 1e-12)
    lat = run("after", "sphere", 3969, "latlong")
    for n in (6144, 16384, 65536):
        check(f"latlong/{n}: documented saturation",
              metrics(lat) == metrics(run("after", "sphere", n, "latlong")))
    check("icosphere: refinement beyond latlong cap",
          int(run("after", "sphere", 65536, "auto")["panels"]) == 81920 and
          int(run("after", "sphere", 6144, "auto")["panels"]) == 20480)
    for name in ("sphere", *(["lig"] if "lig_shift" in fixtures else [])):
        r, s = run("after", name, 384, "auto"), run("after", name + "_shift", 384, "auto")
        # Circle cuts depend on atom-coordinate differences. Translating the
        # decimal fixture perturbs those differences at floating-point scale;
        # retain bit equality for a sphere and a 64-ulp bound for a molecule.
        area_equal = r["area"] == s["area"] if name == "sphere" else (
            abs(float(r["area"])-float(s["area"])) <= 64*math.ulp(float(r["area"])))
        check(f"{name}: translation area and burial count invariant within roundoff",
              area_equal and r["burial_points"] == s["burial_points"] and r["panels"] == s["panels"])
    for name in ("two", *(["lig"] if "lig" in fixtures else [])):
        depths = [float(run("after", name, n, "auto")["burial_depth_max"])
                  for n in (96, 384, 1536)]
        check(f"{name}: burial absent throughout refinement sequence",
              max(depths) < 1e-10)
    check("two spheres: closed-form union area within 1.5% at finest tested mesh",
          abs(float(run("after", "two", 1536, "auto")["area_rel_error"])) < 0.015)
    (a.output / "results.json").write_text(json.dumps(dict(checks=checks, rows=rows), indent=2))
    failed = [x["name"] for x in checks if not x["passed"]]
    print(json.dumps(dict(passed=len(checks)-len(failed), total=len(checks), failed=failed)))
    return bool(failed)


if __name__ == "__main__":
    raise SystemExit(main())
