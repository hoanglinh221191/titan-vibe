"""End-to-end gate for the BIE route MM/PBSA uses.

The route: analytic SES accepted at grid resolution (patch-ray-resolution)
and remeshed onto the SES, curved PN node-patch, classic treecode with the
fast far field, classic ILU(0), energy far-eta 16 and sample order 2.

1. Kirkwood sphere: R = 1.7 A, unit charge at 0.5 R, eps 1/80, 0.15 M,
   against the screened analytic series.
   - h = 1.0 within 1%, h = 0.5 within 0.2%, the error falling with h.
   - The marching mesh (remesh off) at h = 0.5 within 0.5%.
   - The exact energy projection (far-eta 0, sample order 4) within 2e-4 of
     the MM/PBSA defaults.
2. 1SH1 (696 atoms, h = 1.0): the regression value within 1e-4, and the
   panel and iteration counts near theirs.
3. 1A63 (2,065 atoms, h = 0.5, surface only): the build that one
   uncertified ray row used to abort must succeed, with no vertex left to
   the morphology predicate and the volume on its h trend.
4. A solve that runs out of GMRES iterations must not pass for a result:
   exit code non-zero, status failed, NaN energies, no energy in
   pbe_summary.dat.
5. The p1 FMM route with the curved node-patch cell (P1_ROUTE, the TITAN
   BIE-PB paper protocol): the Kirkwood sphere at h = 0.5 within 2e-3 on the
   curved cell, and at theta 0.05, where no pair is left to the far field,
   the accelerated operator equal to the dense audit rows to 1e-12 - the
   packed near loop, its corrections and the near kernel are one operator.
6. Charge-adaptive remeshing on P1_ROUTE: a charge 0.34 A below the sphere
   surface at h = 1.0 within 2.5e-3 (the uniform mesh is 5.6% out).
7. The guarded local p4 far field on P1_ROUTE (1SH1, h = 1.0): sampled
   operator audit below 1e-4.
"""
import argparse
import csv
import hashlib
import json
import math
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

# MM/PBSA settings of docs/TITAN_MMPBSA_RESOLUTION_SES_2026-10-06.md.
ROUTE = [
    ("-pbe", "bie"),
    ("-pbe-surface", "vdw"),
    ("-pbe-bie-primary-surface", "analytic-ses"),
    ("-pbe-bie-analytic-ses-engine", "patch-ray-resolution"),
    ("-pbe-bie-ses-grid-spacing", "1.0"),
    ("-pbe-bie-discretization", "node-patch"),
    ("-pbe-bie-sum", "treecode"),
    ("-pbe-bie-formulation", "jcp2013-treecode"),
    ("-pbe-bie-fmm-m2l", "off"),
    ("-pbe-bie-tree-builder", "topdown-aabb"),
    ("-pbe-bie-tree-leaf-size", "8"),
    ("-pbe-bie-leaf-mac", "on"),
    ("-pbe-bie-opening-theta", "0.7"),
    ("-pbe-bie-taylor-order", "4"),
    ("-pbe-bie-precond", "classic-nearfield-ilu0"),
    ("-pbe-bie-precond-ilu-neighbor-cap", "64"),
    ("-pbe-bie-precond-memory-budget-mb", "4096"),
    ("-pbe-bie-curved-pair-cache-mb", "64"),
    ("-pbe-bie-node-patch-operator-audit", "off"),
    ("-pbe-bie-energy-far-eta", "16"),
    ("-pbe-bie-curved-sample-order", "2"),
    ("-pbe-bie-gmres-restart", "32"),
    ("-pbe-max-iter", "240"),
    ("-pbe-rms", "1e-4"),
    ("-pbe-dphi", "1e-4"),
    ("-pbe-int-diel", "1"),
    ("-pbe-ext-diel", "80"),
    ("-pbe-salt", "0.15"),
    ("-pbe-temp", "298.15"),
    ("-pbe-probe", "1.4"),
    ("-pbe-timing", "detail"),
]

# The p1 FMM route of the TITAN BIE-PB paper: degree-one Cartesian multipole
# far field, FMM-near ILU(0), curved node-patch with its near zone.
P1_ROUTE = [
    ("-pbe", "bie"),
    ("-pbe-surface", "vdw"),
    ("-pbe-bie-primary-surface", "analytic-ses"),
    ("-pbe-bie-analytic-ses-engine", "patch-ray-resolution"),
    ("-pbe-bie-ses-grid-spacing", "1.0"),
    ("-pbe-bie-discretization", "node-patch"),
    ("-pbe-bie-fast-operator", "p1"),
    ("-pbe-bie-tree-leaf-size", "64"),
    ("-pbe-bie-opening-theta", "0.45"),
    ("-pbe-bie-near-integration", "analytic_singular"),
    ("-pbe-bie-quad-buffer", "off"),
    ("-pbe-bie-precond", "fmm-nearfield-ilu0"),
    ("-pbe-bie-precond-ilu-neighbor-cap", "64"),
    ("-pbe-bie-precond-memory-budget-mb", "4096"),
    ("-pbe-bie-matvec-audit", "sampled"),
    ("-pbe-bie-matvec-audit-sample-count", "128"),
    ("-pbe-bie-gmres-restart", "128"),
    ("-pbe-max-iter", "240"),
    ("-pbe-rms", "1e-6"),
    ("-pbe-dphi", "1e-6"),
    ("-pbe-int-diel", "1"),
    ("-pbe-ext-diel", "80"),
    ("-pbe-salt", "0.15"),
    ("-pbe-temp", "298.15"),
    ("-pbe-probe", "1.4"),
    ("-pbe-timing", "detail"),
]

# The guarded degree-four accuracy control over P1_ROUTE.
P1_LOCAL_P4 = {
    "-pbe-bie-fmm-local": "gmres-p4",
    "-pbe-bie-fmm-local-refine": "operator-tail-static",
    "-pbe-bie-fmm-local-refine-tol": "0.002",
    "-pbe-bie-fmm-local-cache": "green-jet-compressed",
    "-pbe-bie-fmm-local-cache-budget-mb": "256",
    "-pbe-bie-fmm-local-fallback": "packed-avx2",
    "-pbe-bie-fmm-local-plan-builder": "metric-cache",
    "-pbe-bie-fmm-local-refine-probe-count": "4",
    "-pbe-bie-fmm-local-refine-probe-mode": "rademacher",
    "-pbe-bie-fmm-local-geometry-guard": "strict-geometry-mac",
    "-pbe-bie-fmm-local-geometry-theta": "0.5",
    "-pbe-bie-taylor-order": "4",
    "-pbe-bie-matvec-audit": "sampled-detail",
    "-pbe-bie-matvec-audit-sample-count": "128",
    "-pbe-rms": "1e-4",
    "-pbe-dphi": "1e-4",
}

# 1SH1 at h = 1.0 with ROUTE, recorded 2026-10-09 (clang, Windows), after the
# resolution projection stopped moving vertices onto another sheet than their
# own (3 refused here; -725.8708031 kcal/mol and 5773 panels before). A change
# beyond the tolerance is a change of the route's numerics, not rounding.
# Re-recorded 2026-10-14 when node-patch cells took the analytic SES normal at
# their vertex instead of the area-weighted star normal (-725.729078 before).
# On CHARMM 1SH1 over h = 1.25..0.5 the two normals converge alike and differ
# by at most 0.03% at any one spacing.
# Re-recorded 2026-10-18: the PN lift is proven unfolded and apart from every
# other patch except at a common vertex or edge, every near pair takes the
# admissible subdivision, and a bad patch is mended by its vertex normals
# (the star normal for patches and collocation alike) instead of being
# flattened, so the patches keep one tangent plane at each vertex.  The
# vertices of the few triangles still flat (mesh folds) are corners, whose
# cells collocate at a smooth point of a patch, and every self block is
# converged per kernel in the vertex frame (-725.4504757 before, 2.4e-4).
SH1_ENERGY_KCAL = -725.6301695
SH1_PANELS = 5770
SH1_ITERATIONS = 11

# 1A63 at h = 0.5, surface only, recorded 2026-10-10 (clang, Windows). The
# neighbouring spacings give 17749.0 (h = 0.6) and 17755.9 (h = 0.45).
A63_H05_VOLUME_A3 = 17752.9


def kirkwood_kcal(radius, depth, eps_in, eps_out, kappa, lmax=80):
    """Reaction energy (kcal/mol) of a unit charge at `depth` from the centre
    of a sphere of `radius`, ions excluded from the sphere, screened outside.
    Uses modified spherical Bessel functions k_l; their normalisation cancels
    in k_l'/k_l."""
    coulomb = 332.063713299
    x = kappa * radius
    k = [math.exp(-x) / x, math.exp(-x) * (1 + x) / x ** 2]
    for l in range(1, lmax + 1):
        k.append(k[l - 1] + (2 * l + 1) / x * k[l])
    energy = 0.0
    for l in range(lmax + 1):
        dk = -k[1] if l == 0 else -k[l - 1] - (l + 1) / x * k[l]
        ratio = kappa * dk / k[l]
        num = (l + 1) / radius ** (l + 2) + (eps_out / eps_in) * ratio / radius ** (l + 1)
        den = eps_in * l * radius ** (l - 1) - eps_out * ratio * radius ** l
        energy += 0.5 * coulomb * depth ** l * num / den * depth ** l
    return energy


def titan_kappa(salt_m, temp_k, eps_out):
    """Debye kappa (1/A) exactly as TITAN computes it (debye_kappa_inv_a)."""
    return 0.3291 * math.sqrt(salt_m) * math.sqrt(298.15 / temp_k) * math.sqrt(78.5 / eps_out)


def write_inputs(work):
    # A neutral four-atom driver; the PB source is the PQR.
    (work / "driver.pdb").write_text(
        "CRYST1   50.000   50.000   50.000  90.00  90.00  90.00 P 1           1\n"
        "ATOM      1  N   VAL A   1      -2.000   0.000   0.000  1.00  0.00           N\n"
        "ATOM      2  CA  VAL A   1      -0.600   0.100   0.000  1.00  0.00           C\n"
        "ATOM      3  C   VAL A   1       0.000   1.450   0.000  1.00  0.00           C\n"
        "ATOM      4  O   VAL A   1      -0.600   2.430   0.000  1.00  0.00           O\n"
        "END\n")
    (work / "driver.top").write_text(
        "[ defaults ]\n1 2 yes 0.5 0.8333333333\n[ atomtypes ]\nT 12.0 0.0 A 0.3 0.1\n"
        "[ moleculetype ]\nTEST 3\n[ atoms ]\n"
        "1 T 1 VAL N 1 0.0 12.0\n2 T 1 VAL CA 2 0.0 12.0\n"
        "3 T 1 VAL C 3 0.0 12.0\n4 T 1 VAL O 4 0.0 12.0\n"
        "[ system ]\nNeutral four-atom driver; PB source is the PQR\n[ molecules ]\nTEST 1\n")
    # Neutral sphere of 1.7 A with the unit charge inside it at 0.85 A.
    (work / "kirkwood.pqr").write_text(
        "ATOM      1 C   MOL A   1       0.000   0.000   0.000  0.0000 1.7000\n"
        "ATOM      2 C   MOL A   2       0.850   0.000   0.000  1.0000 0.1000\n"
        "END\n")
    # The same sphere with the charge 0.34 A below its surface.
    (work / "kirkwood_shallow.pqr").write_text(
        "ATOM      1 C   MOL A   1       0.000   0.000   0.000  0.0000 1.7000\n"
        "ATOM      2 C   MOL A   2       1.360   0.000   0.000  1.0000 0.1000\n"
        "END\n")
    (work / "lib").mkdir(exist_ok=True)


def read_csv_last(path):
    rows = list(csv.DictReader(path.open(encoding="utf-8")))
    if not rows:
        raise RuntimeError("empty " + str(path))
    return rows[-1]


def read_summary(path):
    out = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if line and not line.startswith("#") and " " in line:
            key, value = line.split(" ", 1)
            out[key] = value.strip()
    return out


def run_case(binary, work, name, pqr, threads, overrides=(), route=ROUTE):
    case = work / name
    case.mkdir(exist_ok=True)
    options = dict(route)
    options.update(overrides)
    lines = [
        "-f=" + (work / "driver.pdb").as_posix(),
        "-1pdb=true",
        "-common_topo=" + (work / "driver.top").as_posix(),
        "-ff=AMBER99SB",
        "-g=1&1",
        "-path-a=" + (work / "lib").as_posix() + "/",
        "-path-n=" + (work / "lib").as_posix() + "/",
        "-pbe-test-pqr=" + pqr.as_posix(),
        "-nt=" + str(threads),
    ] + [key + "=" + value for key, value in options.items()]
    (case / "arg.in").write_text("\n".join(lines) + "\n")
    result = subprocess.run([str(binary), "-a", "arg.in"], cwd=case, capture_output=True,
                            text=True, encoding="utf-8", errors="replace", timeout=600)
    output = result.stdout + result.stderr
    (case / "result.log").write_text(output, encoding="utf-8")
    energy = profile = summary = None
    if (case / "pbe_energy.csv").exists():
        energy = read_csv_last(case / "pbe_energy.csv")
    if (case / "pbe_tabi_profile.csv").exists():
        profile = read_csv_last(case / "pbe_tabi_profile.csv")
    if (case / "pbe_summary.dat").exists():
        summary = read_summary(case / "pbe_summary.dat")
    return result.returncode, output, energy, profile, summary


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def solved(name, returncode, output, energy, profile, summary):
    """A successful solve: exit 0, status ok, a finite energy that the summary
    counts, and the remesh actually ran unless it was switched off."""
    require(returncode == 0, f"{name}: exit {returncode}\n{output[-3000:]}")
    require(energy is not None and profile is not None and summary is not None,
            f"{name}: missing outputs")
    require(energy["status"] == "ok", f"{name}: status {energy['status']}")
    value = float(energy["E_rxn_kcalmol"])
    require(math.isfinite(value), f"{name}: non-finite energy")
    require(summary.get("energy_frames") == "1", f"{name}: summary energy_frames "
            + str(summary.get("energy_frames")))
    surface = re.search(r"TABI SES patch-ray-resolution: (.*)", output)
    require(surface is not None, f"{name}: no patch-ray-resolution surface line")
    return value, surface.group(1)


def remesh_iterations(surface_line):
    found = re.search(r"remesh_iterations=(\d+)", surface_line)
    return int(found.group(1)) if found else 0


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--work-root", type=Path,
                        default=ROOT / "test/tmp/pbe_bie_resolution_ses")
    parser.add_argument("--threads", type=int, default=2)
    args = parser.parse_args()
    binary = args.binary.resolve()
    work = args.work_root.resolve()
    work.mkdir(parents=True, exist_ok=True)
    write_inputs(work)
    kirkwood = work / "kirkwood.pqr"
    reference = kirkwood_kcal(1.7, 0.85, 1.0, 80.0, titan_kappa(0.15, 298.15, 80.0))
    records = []

    def record(name, **fields):
        records.append(dict(case=name, status="PASS", **fields))
        print("PASS BIE resolution SES", name,
              " ".join(f"{k}={v}" for k, v in fields.items()), flush=True)

    # 1. Kirkwood.
    errors, values = {}, {}
    for h, gate in (("1.0", 1.0e-2), ("0.5", 2.0e-3)):
        name = "kirkwood_h" + h
        r = run_case(binary, work, name, kirkwood, args.threads,
                     {"-pbe-bie-ses-grid-spacing": h})
        value, surface = solved(name, *r)
        require(remesh_iterations(surface) > 0, f"{name}: the remesh did not run")
        values[h] = value
        errors[h] = (value - reference) / abs(reference)
        require(abs(errors[h]) < gate,
                f"{name}: {value} vs Kirkwood {reference}, error {errors[h]:.3e} > {gate}")
        record(name, energy_kcal=value, reference_kcal=round(reference, 8),
               rel_error=f"{errors[h]:.3e}")
    require(abs(errors["0.5"]) < abs(errors["1.0"]),
            f"Kirkwood error did not fall with h: {errors}")

    name = "kirkwood_h0.5_marching"
    r = run_case(binary, work, name, kirkwood, args.threads,
                 {"-pbe-bie-ses-grid-spacing": "0.5", "-pbe-bie-ses-remesh-edge": "off"})
    value, surface = solved(name, *r)
    require(remesh_iterations(surface) == 0, f"{name}: remesh ran although off")
    marching_error = (value - reference) / abs(reference)
    require(abs(marching_error) < 5.0e-3, f"{name}: error {marching_error:.3e}")
    record(name, energy_kcal=value, rel_error=f"{marching_error:.3e}")

    name = "kirkwood_h0.5_exact_projection"
    r = run_case(binary, work, name, kirkwood, args.threads,
                 {"-pbe-bie-ses-grid-spacing": "0.5", "-pbe-bie-energy-far-eta": "0",
                  "-pbe-bie-curved-sample-order": "4"})
    exact, _ = solved(name, *r)
    gap = abs(exact - values["0.5"]) / abs(exact)
    require(gap < 2.0e-4, f"{name}: far-eta 16 / order 2 differ from exact by {gap:.3e}")
    record(name, energy_kcal=exact, gap_to_defaults=f"{gap:.3e}")

    # 2. 1SH1 regression.
    name = "1sh1_h1.0"
    pqr = ROOT / "test/pbe_inputs_public/1sh1_full_amber_ph7.pqr"
    r = run_case(binary, work, name, pqr, args.threads)
    value, surface = solved(name, *r)
    profile = r[3]
    drift = (value - SH1_ENERGY_KCAL) / abs(SH1_ENERGY_KCAL)
    require(abs(drift) < 1.0e-4,
            f"{name}: {value} vs recorded {SH1_ENERGY_KCAL} kcal/mol, drift {drift:.3e}")
    panels = int(profile["panels"])
    iterations = int(profile["iterations"])
    require(abs(panels - SH1_PANELS) <= 0.01 * SH1_PANELS,
            f"{name}: {panels} panels vs {SH1_PANELS}")
    require(iterations <= 2 * SH1_ITERATIONS,
            f"{name}: {iterations} GMRES iterations vs {SH1_ITERATIONS}")
    record(name, energy_kcal=value, drift=f"{drift:.3e}", panels=panels,
           iterations=iterations)

    # 3. 1A63 at h = 0.5: one ray row's torus quartic stays uncertified (with
    # clang on Windows). The surface must still be built, its vertices
    # classified by certified parity nearby rather than by the sampled
    # morphology predicate, and its volume on the h trend.
    name = "1a63_h0.5_surface"
    pqr = ROOT / "test/pbe_inputs_public/1a63_full_amber_ph7.pqr"
    returncode, output, _, _, _ = run_case(
        binary, work, name, pqr, args.threads,
        {"-pbe-bie-ses-grid-spacing": "0.5", "-pbe-bie-surface-only": "true"})
    require(returncode == 0, f"{name}: exit {returncode}\n{output[-3000:]}")
    surface = re.search(r"TABI SES patch-ray-resolution: (.*)", output)
    require(surface is not None, f"{name}: no patch-ray-resolution surface line")
    counters = {key: int(value) for key, value in re.findall(
        r"(uncertified_ray_rows|nearby_parity_vertices|fallback_sign_vertices)=(\d+)",
        surface.group(1))}
    require(counters.get("fallback_sign_vertices") == 0,
            f"{name}: {counters} - a vertex reached the morphology predicate")
    quality = read_csv_last(work / name / "pbe_tabi_surface_quality.csv")
    volume = float(quality["enclosed_volume"])
    require(abs(volume - A63_H05_VOLUME_A3) <= 1.0e-3 * A63_H05_VOLUME_A3,
            f"{name}: enclosed volume {volume} vs {A63_H05_VOLUME_A3} A^3")
    record(name, volume_A3=volume, **counters)

    # 4. A failed solve.
    name = "kirkwood_gmres_exhausted"
    returncode, output, energy, profile, summary = run_case(
        binary, work, name, kirkwood, args.threads,
        {"-pbe-bie-ses-grid-spacing": "0.5", "-pbe-max-iter": "2"})
    require(returncode != 0, f"{name}: a failed solve exited 0")
    require("did not succeed" in output, f"{name}: no failure message")
    require(energy is not None and energy["status"] == "failed",
            f"{name}: status {energy and energy['status']}")
    for column in ("E_rxn_kT", "E_rxn_kJmol", "E_rxn_kcalmol"):
        require(math.isnan(float(energy[column])),
                f"{name}: {column}={energy[column]} for a failed solve")
    require(summary is not None and summary.get("energy_frames") == "0",
            f"{name}: summary counted the failed frame's energy")
    require(math.isnan(float(summary.get("energy_kcalmol_mean", "0"))),
            f"{name}: summary mean {summary.get('energy_kcalmol_mean')}")
    record(name, exit_code=returncode, reason=profile["convergence_reason"])

    # 5. The p1 FMM route with the curved cell.
    name = "p1_curved_kirkwood_h0.5"
    r = run_case(binary, work, name, kirkwood, args.threads,
                 {"-pbe-bie-ses-grid-spacing": "0.5"}, P1_ROUTE)
    value, _ = solved(name, *r)
    profile = r[3]
    require(profile["node_patch_representation"] == "curved_pn",
            f"{name}: representation {profile['node_patch_representation']}")
    error = (value - reference) / abs(reference)
    require(abs(error) < 2.0e-3,
            f"{name}: {value} vs Kirkwood {reference}, error {error:.3e} > 2e-3")
    record(name, energy_kcal=value, rel_error=f"{error:.3e}",
           policy=profile["near_kernel_policy"])

    # With no far field the accelerated operator is the dense one, in the
    # pre-GMRES probe and in the manuscript's audit of the solution
    # (-pbe-bie-matvec-audit=targets: min(512, N) targets, both equations),
    # whose direct-row residual is then the solve's own residual.
    name = "p1_curved_operator_theta0.05"
    r = run_case(binary, work, name, kirkwood, args.threads,
                 {"-pbe-bie-ses-grid-spacing": "0.5",
                  "-pbe-bie-opening-theta": "0.05",
                  "-pbe-bie-matvec-audit": "targets",
                  "-pbe-bie-matvec-audit-sample-count": "512"}, P1_ROUTE)
    solved(name, *r)
    profile = r[3]
    probe = float(profile["fmm_cartesian_pre_audit_aggregate_rel"])
    audit = float(profile["tree_vs_direct_matvec_rel"])
    rows = int(profile["matvec_audit_rows"])
    nodes = int(r[2]["solute_points"])
    direct = float(profile["fmm_local_sampled_direct_true_residual"])
    true_residual = float(profile["true_residual_rel"])
    require(probe < 1.0e-12 and audit < 1.0e-12,
            f"{name}: accelerated and dense rows differ by {probe:.3e} (probe), "
            f"{audit:.3e} (solution) with no far field")
    require(rows == 2 * min(512, nodes),
            f"{name}: {rows} audit rows for {nodes} nodes")
    require(0.0 < direct < 1.0e-3 and direct < 100 * max(true_residual, 1e-12),
            f"{name}: sampled direct residual {direct:.3e}, true residual {true_residual:.3e}")
    record(name, operator_audit_rel=f"{audit:.3e}", probe_rel=f"{probe:.3e}",
           audit_rows=rows, direct_residual=f"{direct:.3e}")

    # 6. Charge-adaptive remeshing: a charge 0.34 A below the surface at
    # h = 1.0, which the uniform mesh misses by 5.6%.
    name = "p1_curved_charge_refine_shallow_h1.0"
    shallow = kirkwood_kcal(1.7, 1.36, 1.0, 80.0, titan_kappa(0.15, 298.15, 80.0))
    r = run_case(binary, work, name, work / "kirkwood_shallow.pqr", args.threads,
                 {"-pbe-bie-ses-grid-spacing": "1.0",
                  "-pbe-bie-ses-charge-refine": "0.35",
                  "-pbe-bie-ses-charge-refine-min-edge": "0.1"}, P1_ROUTE)
    value, surface = solved(name, *r)
    refined = re.search(r"refined_vertices=(\d+)", surface)
    require(refined is not None and int(refined.group(1)) > 0,
            f"{name}: no vertex was refined")
    error = (value - shallow) / abs(shallow)
    require(abs(error) < 2.5e-3,
            f"{name}: {value} vs Kirkwood {shallow}, error {error:.3e} > 2.5e-3")
    record(name, energy_kcal=value, rel_error=f"{error:.3e}",
           refined_vertices=int(refined.group(1)))

    # 7. The guarded local p4 far field on the same curved operator: its
    # packed fallback must carry the node-patch near corrections too, or the
    # sampled audit stays at the 1e-2 level instead of 1e-5.
    name = "1sh1_h1.0_p1_curved_local_p4"
    r = run_case(binary, work, name,
                 ROOT / "test/pbe_inputs_public/1sh1_full_amber_ph7.pqr",
                 args.threads, P1_LOCAL_P4, P1_ROUTE)
    value, _ = solved(name, *r)
    audit = float(r[3]["tree_vs_direct_matvec_rel"])
    require(audit < 1.0e-4, f"{name}: sampled operator audit {audit:.3e} >= 1e-4")
    record(name, energy_kcal=value, operator_audit_rel=f"{audit:.3e}")

    (work / "verification.json").write_text(json.dumps(dict(
        binary_sha256=hashlib.sha256(binary.read_bytes()).hexdigest(),
        kirkwood_reference_kcal=reference, cases=records), indent=2) + "\n")
    print("PASS pbe_bie_resolution_ses_smoke", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
