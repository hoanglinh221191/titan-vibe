#!/usr/bin/env python3
"""Build matched TITAN/APBS surface-potential assets for PyMOL and VTK tools."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path


def read_vtk(path: Path):
    lines = path.read_text(encoding="utf-8").splitlines()
    points = []
    triangles = []
    scalars = {}
    normals = []
    i = 0
    point_count = 0
    while i < len(lines):
        fields = lines[i].split()
        if fields[:1] == ["POINTS"]:
            point_count = int(fields[1])
            for row in lines[i + 1 : i + 1 + point_count]:
                values = [float(value) for value in row.split()]
                points.append(tuple(values[:3]))
            i += point_count + 1
            continue
        if fields[:1] == ["POLYGONS"]:
            polygon_count = int(fields[1])
            for row in lines[i + 1 : i + 1 + polygon_count]:
                values = [int(value) for value in row.split()]
                if values[0] != 3:
                    raise ValueError("Only triangular TITAN PolyData is supported")
                triangles.append(tuple(values[1:4]))
            i += polygon_count + 1
            continue
        if fields[:1] == ["POINT_DATA"]:
            if int(fields[1]) != point_count:
                raise ValueError("TITAN VTK POINT_DATA does not match POINTS")
            i += 1
            while i < len(lines):
                fields = lines[i].split()
                if fields[:1] == ["SCALARS"]:
                    name = fields[1]
                    if i + 1 >= len(lines) or not lines[i + 1].startswith("LOOKUP_TABLE"):
                        raise ValueError(f"Missing LOOKUP_TABLE for {name}")
                    values = [float(value) for value in lines[i + 2 : i + 2 + point_count]]
                    if len(values) != point_count:
                        raise ValueError(f"Incomplete scalar field {name}")
                    scalars[name] = values
                    i += point_count + 2
                    continue
                if fields[:1] == ["VECTORS"]:
                    normals = [
                        tuple(float(value) for value in row.split()[:3])
                        for row in lines[i + 1 : i + 1 + point_count]
                    ]
                    i += point_count + 1
                    continue
                i += 1
            break
        i += 1
    if not points or not triangles or "total_potential_kT_e" not in scalars:
        raise ValueError(
            "TITAN VTK lacks connected POINT_DATA total_potential_kT_e")
    if not normals:
        normals = [(0.0, 0.0, 1.0)] * len(points)
    return points, triangles, normals, scalars


def read_dx(path: Path):
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    counts = None
    origin = None
    deltas = []
    values = []
    reading = False
    expected = None
    for line in lines:
        fields = line.split()
        if len(fields) >= 8 and fields[:5] == ["object", "1", "class", "gridpositions", "counts"]:
            counts = tuple(int(value) for value in fields[5:8])
        elif fields[:1] == ["origin"] and len(fields) >= 4:
            origin = tuple(float(value) for value in fields[1:4])
        elif fields[:1] == ["delta"] and len(fields) >= 4:
            deltas.append(tuple(float(value) for value in fields[1:4]))
        elif "data follows" in line and "class array" in line:
            if "items" not in fields:
                raise ValueError("OpenDX array lacks item count")
            expected = int(fields[fields.index("items") + 1])
            reading = True
        elif reading:
            if line.startswith(("attribute", "object", "component")):
                reading = False
                continue
            for value in fields:
                values.append(float(value))
    if counts is None or origin is None or len(deltas) != 3:
        raise ValueError("OpenDX grid geometry is incomplete")
    if expected is None or len(values) < expected:
        raise ValueError("OpenDX potential array is incomplete")
    values = values[:expected]
    if math.prod(counts) != len(values):
        raise ValueError("OpenDX item count does not match grid dimensions")
    return counts, origin, deltas, values


def sample_dx(point, counts, origin, deltas, values):
    # APBS writes an axis-aligned grid and uses z as the fastest data index.
    spacing = []
    for axis, delta in enumerate(deltas):
        spacing.append(delta[axis])
        if abs(delta[(axis + 1) % 3]) > 1.0e-12 or abs(delta[(axis + 2) % 3]) > 1.0e-12:
            raise ValueError("Only axis-aligned APBS OpenDX grids are supported")
    u = [(point[axis] - origin[axis]) / spacing[axis] for axis in range(3)]
    base = [math.floor(value) for value in u]
    frac = [u[axis] - base[axis] for axis in range(3)]
    if any(base[axis] < 0 or base[axis] + 1 >= counts[axis] for axis in range(3)):
        return math.nan

    def at(ix, iy, iz):
        return values[(ix * counts[1] + iy) * counts[2] + iz]

    total = 0.0
    for dx in (0, 1):
        wx = frac[0] if dx else 1.0 - frac[0]
        for dy in (0, 1):
            wy = frac[1] if dy else 1.0 - frac[1]
            for dz in (0, 1):
                wz = frac[2] if dz else 1.0 - frac[2]
                total += wx * wy * wz * at(base[0] + dx, base[1] + dy, base[2] + dz)
    return total


def write_comparison_vtk(path, points, triangles, normals, titan, apbs_total,
                         apbs_reference, apbs_reaction):
    with path.open("w", encoding="ascii", newline="\n") as out:
        out.write("# vtk DataFile Version 3.0\n")
        out.write("6LZG TITAN and APBS potential on matched analytic-SES nodes; unit=kT/e\n")
        out.write("ASCII\nDATASET POLYDATA\n")
        out.write(f"POINTS {len(points)} double\n")
        for xyz in points:
            out.write("{:.17g} {:.17g} {:.17g}\n".format(*xyz))
        out.write(f"POLYGONS {len(triangles)} {4 * len(triangles)}\n")
        for triangle in triangles:
            out.write("3 {} {} {}\n".format(*triangle))
        out.write(f"POINT_DATA {len(points)}\n")
        for name, data in (
            ("titan_total_potential_kT_e", titan),
            ("apbs_total_potential_kT_e", apbs_total),
            ("apbs_reference_potential_kT_e", apbs_reference),
            ("apbs_reaction_potential_kT_e", apbs_reaction),
            ("titan_minus_apbs_total_kT_e", [a - b for a, b in zip(titan, apbs_total)]),
            ("titan_minus_apbs_reaction_kT_e", [a - b for a, b in zip(titan, apbs_reaction)]),
        ):
            out.write(f"SCALARS {name} double 1\nLOOKUP_TABLE default\n")
            for value in data:
                out.write(f"{value:.17g}\n")
        out.write("VECTORS normal double\n")
        for normal in normals:
            out.write("{:.17g} {:.17g} {:.17g}\n".format(*normal))


def write_pymol_loader(path: Path, vtk_name: str, pqr_name: str,
                       titan_limit: float, apbs_limit: float,
                       apbs_field: str):
    template = r'''from pathlib import Path
from pymol import cmd
from pymol.cgo import BEGIN, COLOR, END, NORMAL, TRIANGLES, VERTEX

ROOT = Path.cwd()
if not (ROOT / __VTK__).exists():
    ROOT = Path(__FALLBACK_ROOT__)
VTK = ROOT / __VTK__
PQR = ROOT / __PQR__
TITAN_LIMIT = __TITAN_LIMIT__
APBS_LIMIT = __APBS_LIMIT__

def read_vtk(path):
    lines = path.read_text().splitlines()
    points, triangles, normals, fields = [], [], [], {}
    n = 0
    i = 0
    while i < len(lines):
        f = lines[i].split()
        if f[:1] == ['POINTS']:
            n = int(f[1]); points = [tuple(map(float, r.split()[:3])) for r in lines[i+1:i+1+n]]; i += n + 1; continue
        if f[:1] == ['POLYGONS']:
            m = int(f[1]); triangles = [tuple(map(int, r.split()[1:4])) for r in lines[i+1:i+1+m]]; i += m + 1; continue
        if f[:1] == ['SCALARS']:
            name = f[1]; fields[name] = [float(v) for v in lines[i+2:i+2+n]]; i += n + 2; continue
        if f[:1] == ['VECTORS']:
            normals = [tuple(map(float, r.split()[:3])) for r in lines[i+1:i+1+n]]; i += n + 1; continue
        i += 1
    return points, triangles, normals, fields

def bwr(value, limit):
    x = max(-1.0, min(1.0, value / limit))
    if x < 0.0:
        return (1.0, 1.0 + x, 1.0 + x)
    return (1.0 - x, 1.0 - x, 1.0)

def mesh_cgo(points, triangles, normals, values, shift, limit):
    obj = [BEGIN, TRIANGLES]
    for tri in triangles:
        for idx in tri:
            x, y, z = points[idx]
            nx, ny, nz = normals[idx]
            r, g, b = bwr(values[idx], limit)
            obj += [COLOR, r, g, b, NORMAL, nx, ny, nz, VERTEX, x + shift, y, z]
    obj.append(END)
    return obj

points, triangles, normals, fields = read_vtk(VTK)
span = max(p[0] for p in points) - min(p[0] for p in points)
offset = 0.60 * span + 8.0
cmd.load_cgo(mesh_cgo(points, triangles, normals, fields['titan_total_potential_kT_e'], -offset, TITAN_LIMIT), 'TITAN_analytic_SES')
cmd.load_cgo(mesh_cgo(points, triangles, normals, fields[__APBS_FIELD__], offset, APBS_LIMIT), 'APBS_on_TITAN_SES')
cmd.load(str(PQR), 'structure_left')
cmd.create('structure_right', 'structure_left')
cmd.translate([-offset, 0.0, 0.0], 'structure_left')
cmd.translate([offset, 0.0, 0.0], 'structure_right')
cmd.hide('everything', 'structure_left or structure_right')
cmd.show('cartoon', 'structure_left or structure_right')
cmd.color('gray70', 'structure_left or structure_right')
cmd.set('cgo_transparency', 0.0)
cmd.set('orthoscopic', 1)
cmd.set('ray_opaque_background', 0)
cmd.bg_color('white')
cmd.zoom('all', buffer=4.0)
'''
    text = (
        template.replace("__VTK__", repr(vtk_name))
        .replace("__PQR__", repr(pqr_name))
        .replace("__TITAN_LIMIT__", repr(float(titan_limit)))
        .replace("__APBS_LIMIT__", repr(float(apbs_limit)))
        .replace("__FALLBACK_ROOT__", repr(str(path.parent.resolve())))
        .replace("__APBS_FIELD__", repr(apbs_field))
    )
    path.write_text(text, encoding="utf-8")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--titan-vtk", required=True, type=Path)
    parser.add_argument("--apbs-dx", required=True, type=Path)
    parser.add_argument("--apbs-reference-dx", required=True, type=Path)
    parser.add_argument("--pqr", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--color-limit-kt-e", type=float, default=5.0)
    parser.add_argument("--color-mode", choices=("independent", "common"), default="independent")
    parser.add_argument("--apbs-field", choices=("total", "reaction"), default="total")
    args = parser.parse_args()

    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    points, triangles, normals, scalars = read_vtk(args.titan_vtk.resolve())
    counts, origin, deltas, dx_values = read_dx(args.apbs_dx.resolve())
    ref_counts, ref_origin, ref_deltas, ref_values = read_dx(args.apbs_reference_dx.resolve())
    if counts != ref_counts or origin != ref_origin or deltas != ref_deltas:
        raise ValueError("APBS solvated and reference OpenDX grids differ")
    titan = scalars["total_potential_kT_e"]
    apbs_total = [sample_dx(point, counts, origin, deltas, dx_values) for point in points]
    apbs_reference = [sample_dx(point, counts, origin, deltas, ref_values) for point in points]
    apbs_reaction = [a - b for a, b in zip(apbs_total, apbs_reference)]
    outside = sum(not math.isfinite(value) for value in apbs_total + apbs_reference)
    if outside:
        raise ValueError(f"{outside} TITAN surface nodes fall outside the APBS grid")

    vtk_path = output / "6lzg_titan_vs_apbs_surface.vtk"
    csv_path = output / "6lzg_titan_vs_apbs_surface.csv"
    pqr_path = output / "6lzg_complex_charmm_ph7.pqr"
    write_comparison_vtk(vtk_path, points, triangles, normals, titan,
                         apbs_total, apbs_reference, apbs_reaction)
    pqr_path.write_bytes(args.pqr.resolve().read_bytes())
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(("node_id", "x_A", "y_A", "z_A", "titan_total_potential_kT_e",
                         "apbs_total_potential_kT_e", "apbs_reference_potential_kT_e",
                         "apbs_reaction_potential_kT_e", "titan_minus_apbs_total_kT_e",
                         "titan_minus_apbs_reaction_kT_e"))
        for index, (point, titan_value, total_value, reference_value, reaction_value) in enumerate(
                zip(points, titan, apbs_total, apbs_reference, apbs_reaction)):
            writer.writerow((index, *point, titan_value, total_value, reference_value,
                             reaction_value, titan_value - total_value,
                             titan_value - reaction_value))

    loader_path = output / "view_6lzg_titan_vs_apbs.py"
    selected_name = ("apbs_total_potential_kT_e" if args.apbs_field == "total"
                     else "apbs_reaction_potential_kT_e")
    selected_values = apbs_total if args.apbs_field == "total" else apbs_reaction
    def robust_limit(values):
        ordered = sorted(abs(value) for value in values if math.isfinite(value))
        return max(1.0e-12, ordered[round(0.95 * (len(ordered) - 1))])
    if args.color_mode == "common":
        titan_limit = args.color_limit_kt_e
        apbs_limit = args.color_limit_kt_e
    else:
        titan_limit = robust_limit(titan)
        apbs_limit = robust_limit(selected_values)
    write_pymol_loader(loader_path, vtk_path.name, pqr_path.name,
                       titan_limit, apbs_limit, selected_name)
    total_delta = [a - b for a, b in zip(titan, apbs_total)]
    reaction_delta = [a - b for a, b in zip(titan, apbs_reaction)]
    total_rms = math.sqrt(sum(value * value for value in total_delta) / len(total_delta))
    reaction_rms = math.sqrt(sum(value * value for value in reaction_delta) / len(reaction_delta))
    with (output / "6lzg_surface_comparison_summary.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(("nodes", "triangles", "apbs_field", "color_mode",
                         "titan_color_limit_kT_e", "apbs_color_limit_kT_e",
                         "titan_min_kT_e", "titan_max_kT_e", "apbs_total_min_kT_e",
                         "apbs_total_max_kT_e", "apbs_reaction_min_kT_e",
                         "apbs_reaction_max_kT_e", "rms_vs_total_kT_e",
                         "rms_vs_reaction_kT_e"))
        writer.writerow((len(points), len(triangles), args.apbs_field,
                         args.color_mode, titan_limit, apbs_limit,
                         min(titan), max(titan),
                         min(apbs_total), max(apbs_total), min(apbs_reaction),
                         max(apbs_reaction), total_rms, reaction_rms))
    print(f"comparison_ready nodes={len(points)} triangles={len(triangles)} "
          f"rms_total_kT_e={total_rms:.6g} rms_reaction_kT_e={reaction_rms:.6g}")


if __name__ == "__main__":
    main()
