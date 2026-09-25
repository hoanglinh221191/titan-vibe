#!/usr/bin/env python3
"""Generate Figures 1-4 for the ILU(0) robustness manuscript story."""

from __future__ import annotations

import argparse
import csv
import math
import subprocess
from pathlib import Path
from xml.sax.saxutils import escape


COLORS = {
    "ink": "#222222",
    "muted": "#666666",
    "grid": "#D9D9D9",
    "blue": "#0072B2",
    "light_blue": "#56B4E9",
    "orange": "#D55E00",
    "green": "#009E73",
    "magenta": "#CC79A7",
    "gray": "#777777",
    "light": "#F7F7F7",
}
LABELS = {
    "diagonal": "Diagonal",
    "leaf_block_jacobi": "Leaf block Jacobi",
    "fmm_near_graph_ras": "FMM-near RAS",
    "distance_knn_graph_ras": "Distance-kNN RAS",
    "distance_knn_ilu0": "Distance-kNN ILU(0)",
    "fmm_near_graph_ilu0": "FMM-near ILU(0)",
}
METHOD_COLORS = {
    "diagonal": COLORS["gray"],
    "leaf_block_jacobi": COLORS["orange"],
    "fmm_near_graph_ras": COLORS["blue"],
    "distance_knn_graph_ras": COLORS["light_blue"],
    "distance_knn_ilu0": COLORS["magenta"],
    "fmm_near_graph_ilu0": COLORS["green"],
}


def args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--matrix", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


def t(x: float, y: float, value: str, size: int = 24, anchor: str = "middle", weight: str = "normal", color: str = "#222222") -> str:
    return f'<text x="{x:.1f}" y="{y:.1f}" text-anchor="{anchor}" font-family="Arial" font-size="{size}" font-weight="{weight}" fill="{color}">{escape(value)}</text>'


def box(x: float, y: float, w: float, h: float, fill: str, stroke: str, title: str, lines: list[str]) -> list[str]:
    out = [f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="8" fill="{fill}" stroke="{stroke}" stroke-width="3"/>']
    out.append(t(x + w / 2, y + 34, title, 25, weight="bold"))
    for i, line in enumerate(lines):
        out.append(t(x + w / 2, y + 72 + i * 30, line, 21, color=COLORS["muted"]))
    return out


def arrow(x1: float, y1: float, x2: float, y2: float, color: str = "#444444") -> str:
    return f'<line x1="{x1}" y1="{y1}" x2="{x2}" y2="{y2}" stroke="{color}" stroke-width="4" marker-end="url(#arrow)"/>'


def base_svg(width: int, height: int, mm_height: int) -> list[str]:
    return [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="180mm" height="{mm_height}mm" viewBox="0 0 {width} {height}">',
        '<defs><marker id="arrow" markerWidth="10" markerHeight="10" refX="9" refY="3" orient="auto" markerUnits="strokeWidth"><path d="M0,0 L0,6 L9,3 z" fill="#444444"/></marker></defs>',
        '<rect width="100%" height="100%" fill="white"/>',
    ]


def export(svg: list[str], prefix: Path) -> None:
    prefix.parent.mkdir(parents=True, exist_ok=True)
    svg_path = prefix.with_suffix(".svg")
    svg_path.write_text("\n".join(svg + ["</svg>"]), encoding="ascii")
    subprocess.run(["inkscape", str(svg_path), "--export-type=pdf", f"--export-filename={prefix.with_suffix('.pdf')}"] , check=True)
    subprocess.run(["inkscape", str(svg_path), "--export-type=png", f"--export-filename={prefix.with_suffix('.png')}", "--export-dpi=600"], check=True)


def figure1(output: Path) -> None:
    svg = base_svg(2200, 760, 66)

    # Panel A: analytic patch-ray SES. The drawing is intentionally schematic:
    # its role is to show the three analytic patch classes and ray-to-mesh path.
    svg.append(t(35, 48, "A", 31, anchor="start", weight="bold"))
    svg.append(t(90, 48, "Analytic patch-ray SES", 28, anchor="start", weight="bold"))
    svg.append('<rect x="28" y="72" width="650" height="620" rx="8" fill="#FBFBFB" stroke="#B9B9B9" stroke-width="2"/>')

    atoms = [(195, 330, 86, "#9AB6CC"), (350, 265, 72, "#D7A45E"), (470, 370, 98, "#A4B8CE")]
    for x, y, radius, color in atoms:
        svg.append(f'<circle cx="{x}" cy="{y}" r="{radius}" fill="{color}" fill-opacity="0.82" stroke="#555555" stroke-width="2.5"/>')
    contour = "M92,360 C95,230 170,145 292,165 C345,118 438,150 470,220 C585,225 626,315 590,408 C560,520 430,555 345,505 C252,550 120,500 92,360 Z"
    svg.append(f'<path d="{contour}" fill="#DDF3EA" fill-opacity="0.35" stroke="{COLORS["green"]}" stroke-width="7"/>')
    svg.append(f'<path d="M280,168 C326,195 349,210 385,187" fill="none" stroke="{COLORS["orange"]}" stroke-width="7"/>')
    svg.append(f'<path d="M386,188 C411,213 437,225 470,220" fill="none" stroke="{COLORS["magenta"]}" stroke-width="7"/>')
    svg.append(t(145, 165, "contact", 21, weight="bold", color=COLORS["green"]))
    svg.append(t(330, 120, "toroidal", 21, weight="bold", color=COLORS["orange"]))
    svg.append(t(495, 155, "re-entrant", 21, weight="bold", color=COLORS["magenta"]))

    ray_y = (270, 355, 440)
    ray_hits = ((94, 582), (92, 600), (115, 565))
    for y, (x1, x2) in zip(ray_y, ray_hits):
        svg.append(f'<line x1="55" y1="{y}" x2="635" y2="{y}" stroke="#6D6D6D" stroke-width="2" stroke-dasharray="9 8"/>')
        svg.append(f'<circle cx="{x1}" cy="{y}" r="6" fill="white" stroke="{COLORS["blue"]}" stroke-width="3"/>')
        svg.append(f'<circle cx="{x2}" cy="{y}" r="6" fill="white" stroke="{COLORS["blue"]}" stroke-width="3"/>')
    svg.append(t(545, 475, "ray intersections", 20, weight="bold", color=COLORS["blue"]))
    mesh_nodes = [(110, 520), (190, 535), (270, 515), (350, 505), (430, 525), (520, 500), (585, 445)]
    for i, (x, y) in enumerate(mesh_nodes):
        svg.append(f'<circle cx="{x}" cy="{y}" r="5" fill="white" stroke="{COLORS["green"]}" stroke-width="2.5"/>')
        if i:
            px, py = mesh_nodes[i - 1]
            svg.append(f'<line x1="{px}" y1="{py}" x2="{x}" y2="{y}" stroke="{COLORS["green"]}" stroke-width="2.5"/>')
        if i >= 2:
            px, py = mesh_nodes[i - 2]
            svg.append(f'<line x1="{px}" y1="{py}" x2="{x}" y2="{y}" stroke="#91BFAE" stroke-width="1.6"/>')
    svg.append(t(350, 640, "closed oriented mesh", 22, weight="bold", color=COLORS["green"]))

    # Panel B: graph construction from the same coupled surface operator.
    svg.append(t(735, 48, "B", 31, anchor="start", weight="bold"))
    svg.append(t(790, 48, "Coupled near graph", 28, anchor="start", weight="bold"))
    svg.append('<rect x="728" y="72" width="650" height="620" rx="8" fill="#FBFBFB" stroke="#B9B9B9" stroke-width="2"/>')
    svg += box(775, 120, 250, 135, "#F4FAFF", COLORS["blue"], "FMM near graph", ["tree-derived edges"])
    svg += box(775, 345, 250, 135, "#FFF5FA", COLORS["magenta"], "Distance-kNN", ["geometric edges"])
    svg += box(1090, 220, 245, 220, "#F7F7F7", "#555555", "Sparse operator", ["2 x 2 K1-K4 blocks", "self edge", "neighbor cap 64"])
    svg += [arrow(1025, 188, 1090, 275), arrow(1025, 412, 1090, 385)]
    svg.append(t(1050, 560, "same coupled local operator", 22, weight="bold", color=COLORS["ink"]))

    # Panel C: move the former Figure S5 algebra comparison into the main figure.
    svg.append(t(1435, 48, "C", 31, anchor="start", weight="bold"))
    svg.append(t(1490, 48, "Preconditioner algebra", 28, anchor="start", weight="bold"))
    svg.append('<rect x="1428" y="72" width="744" height="620" rx="8" fill="#FBFBFB" stroke="#B9B9B9" stroke-width="2"/>')
    svg += box(1470, 115, 285, 155, "#FFF8F2", COLORS["orange"], "Local action", ["block Jacobi / RAS", "no edge elimination"])
    svg += box(1470, 385, 285, 175, "#F2FBF7", COLORS["green"], "Block ILU(0)", ["zero fill on graph", "triangular solves"])
    svg += box(1835, 115, 285, 155, "#F5F5F5", COLORS["gray"], "System dependent", ["24-720 iterations", "stagnation observed"])
    svg += box(1835, 385, 285, 175, "#F2FBF7", COLORS["green"], "Stable in tests", ["10-48 iterations", "both graph sources"])
    svg += [arrow(1755, 192, 1835, 192), arrow(1755, 472, 1835, 472)]
    svg.append(t(1800, 635, "factorization is the stabilizing step", 23, weight="bold", color=COLORS["green"]))
    export(svg, output / "figure1_ilu_workflow")


def figure2(output: Path) -> None:
    svg = base_svg(1800, 720, 72)
    svg.append(t(220, 55, "Independent graph constructions", 27, weight="bold"))
    svg += box(55, 100, 330, 150, "#F4FAFF", COLORS["blue"], "FMM near-source graph", ["tree-derived neighborhoods"])
    svg += box(55, 330, 330, 150, "#FFF5FA", COLORS["magenta"], "Distance-kNN graph", ["global geometric neighbors"])
    svg += box(480, 185, 350, 220, "#F7F7F7", "#555555", "Coupled sparse operator", ["2 x 2 K1-K4 edge blocks", "self edge + cap 64", "fixed sparsity pattern"])
    svg += [arrow(385, 175, 480, 250), arrow(385, 405, 480, 340)]

    svg.append(t(1090, 55, "Algebra applied to the graph", 27, weight="bold"))
    svg += box(900, 100, 375, 160, "#FFF8F2", COLORS["orange"], "Non-factorized local action", ["block Jacobi or RAS", "no elimination across edges"])
    svg += box(900, 350, 375, 185, "#F2FBF7", COLORS["green"], "Block ILU(0)", ["A approximately L U", "zero fill on graph", "forward / backward solves"])
    svg += [arrow(830, 295, 900, 180), arrow(830, 295, 900, 440)]
    svg += box(1390, 100, 350, 160, "#FFF8F2", COLORS["gray"], "System-dependent", ["22-720 iterations", "frequent stagnation"])
    svg += box(1390, 350, 350, 185, "#F2FBF7", COLORS["green"], "Convergence in tested cases", ["9-48 iterations", "both graph sources"])
    svg += [arrow(1275, 180, 1390, 180), arrow(1275, 440, 1390, 440)]
    svg.append(t(1090, 645, "Within the tested cases, factorized coupled algebra stabilizes both graph constructions.", 24, weight="bold", color=COLORS["green"]))
    export(svg, output / "figure2_ilu_mechanism")


def iteration_fill(iterations: int) -> str:
    if iterations >= 720:
        return "#D9D9D9"
    if iterations <= 50:
        return "#BFE7D6"
    if iterations <= 200:
        return "#FDE3B0"
    return "#F6C4A7"


def figure3(output: Path) -> None:
    svg = base_svg(1800, 760, 76)
    svg.append(t(55, 62, "A", 31, anchor="start", weight="bold"))
    svg.append(t(105, 62, "Voxel-derived watertight SES", 28, anchor="start", weight="bold"))
    svg.append(t(930, 62, "B", 31, anchor="start", weight="bold"))
    svg.append(t(980, 62, "Near-graph block ILU(0)", 28, anchor="start", weight="bold"))

    # Panel A: atoms, voxel morphology and a closed triangulated SES contour.
    for gx in range(80, 850, 38):
        svg.append(f'<line x1="{gx}" y1="100" x2="{gx}" y2="650" stroke="#E8E8E8" stroke-width="1"/>')
    for gy in range(118, 651, 38):
        svg.append(f'<line x1="70" y1="{gy}" x2="850" y2="{gy}" stroke="#E8E8E8" stroke-width="1"/>')
    atoms = [(250, 330, 90), (390, 250, 78), (520, 355, 105), (650, 270, 72)]
    atom_colors = ["#A7B6C8", "#D8A15D", "#7CAAC7", "#B5A0C9"]
    for (x, y, r), color in zip(atoms, atom_colors):
        svg.append(f'<circle cx="{x}" cy="{y}" r="{r}" fill="{color}" fill-opacity="0.75" stroke="#555" stroke-width="2"/>')
    contour = "M135,360 C135,235 210,160 330,168 C405,100 520,135 565,190 C680,155 770,225 770,330 C820,430 735,535 625,520 C550,600 430,570 370,525 C255,565 150,495 135,360 Z"
    svg.append(f'<path d="{contour}" fill="#DDF3EA" fill-opacity="0.55" stroke="{COLORS["green"]}" stroke-width="7"/>')
    nodes = [(145,350),(175,245),(255,180),(345,165),(430,140),(535,175),(625,170),(720,225),(770,320),(785,410),(720,495),(625,520),(530,555),(430,565),(335,525),(235,530),(160,460)]
    for i, (x, y) in enumerate(nodes):
        x2, y2 = nodes[(i + 1) % len(nodes)]
        svg.append(f'<line x1="{x}" y1="{y}" x2="{x2}" y2="{y2}" stroke="{COLORS["green"]}" stroke-width="3"/>')
        if i % 2 == 0:
            svg.append(f'<line x1="{x}" y1="{y}" x2="460" y2="350" stroke="#91BFAE" stroke-width="1.5"/>')
        svg.append(f'<circle cx="{x}" cy="{y}" r="5" fill="white" stroke="{COLORS["green"]}" stroke-width="2.5"/>')
    svg.append(f'<circle cx="745" cy="140" r="28" fill="white" stroke="{COLORS["orange"]}" stroke-width="4"/>')
    svg.append(arrow(745, 170, 710, 215, COLORS["orange"]))
    svg.append(t(745, 105, "probe", 20, color=COLORS["orange"]))
    svg.append(t(455, 690, "closed, oriented triangular boundary", 23, weight="bold", color=COLORS["green"]))

    # Panel B: cross-leaf couplings and block sparse factorization.
    svg.append(f'<line x1="1180" y1="115" x2="1180" y2="625" stroke="#777" stroke-width="3" stroke-dasharray="10 8"/>')
    svg.append(t(1180, 650, "tree-leaf boundary", 20, color=COLORS["muted"]))
    left_nodes = [(1010,180),(1060,245),(1085,330),(1060,420),(1015,500)]
    right_nodes = [(1340,175),(1295,245),(1270,330),(1300,420),(1350,505)]
    for side, color in ((left_nodes,COLORS["blue"]),(right_nodes,COLORS["light_blue"])):
        for i,(x,y) in enumerate(side):
            if i:
                px,py=side[i-1]
                svg.append(f'<line x1="{px}" y1="{py}" x2="{x}" y2="{y}" stroke="{color}" stroke-width="5"/>')
            svg.append(f'<circle cx="{x}" cy="{y}" r="10" fill="white" stroke="{color}" stroke-width="4"/>')
    for li, ri in ((1,1),(2,2),(3,3),(2,1),(2,3)):
        x1,y1=left_nodes[li];x2,y2=right_nodes[ri]
        svg.append(f'<line x1="{x1}" y1="{y1}" x2="{x2}" y2="{y2}" stroke="{COLORS["green"]}" stroke-width="5"/>')
    svg.append(t(1180, 92, "strong cross-leaf couplings", 21, weight="bold", color=COLORS["green"]))
    svg.append(arrow(1370, 340, 1450, 340, COLORS["green"]))

    mx, my, cell = 1470, 185, 48
    pattern = {(0,0),(0,1),(1,0),(1,1),(1,2),(2,1),(2,2),(2,3),(3,2),(3,3),(3,4),(4,3),(4,4)}
    for r in range(5):
        for c in range(5):
            fill = COLORS["green"] if (r,c) in pattern else "white"
            svg.append(f'<rect x="{mx+c*cell}" y="{my+r*cell}" width="{cell-5}" height="{cell-5}" fill="{fill}" fill-opacity="0.82" stroke="#BBBBBB"/>')
    svg.append(t(mx+2.5*cell, my-25, "2 x 2 block sparsity", 21, weight="bold"))
    svg.append(t(mx+2.5*cell, my+5*cell+42, "ILU(0): zero fill", 22, weight="bold", color=COLORS["green"]))
    svg.append(t(1350, 690, "factorization propagates local coupling across the graph", 23, weight="bold", color=COLORS["green"]))
    export(svg, output / "figure3_surface_ilu_illustration")


def figure_s3(rows: list[dict[str, str]], output: Path) -> None:
    svg = base_svg(1800, 920, 92)
    methods = list(LABELS)
    fixtures = [
        ("sphere-patch", "barnase", "Barnase"), ("sphere-patch", "1a63", "1A63"), ("sphere-patch", "1sh1", "1SH1"),
        ("grid-ses", "barnase", "Barnase"), ("grid-ses", "1a63", "1A63"), ("grid-ses", "1sh1", "1SH1"),
    ]
    x0, y0, label_w, cell_w, cell_h = 90, 210, 210, 240, 105
    legend = [("<=50 iterations", "#BFE7D6"), ("51-200", "#FDE3B0"), ("201-719", "#F6C4A7"), ("720 / stalled", "#D9D9D9")]
    for i, (label, color) in enumerate(legend):
        x = 350 + i * 330
        svg.append(f'<rect x="{x}" y="35" width="28" height="22" fill="{color}" stroke="#777"/>')
        svg.append(t(x + 40, 54, label, 21, anchor="start"))
    for ci, method in enumerate(methods):
        svg.append(t(x0 + label_w + ci * cell_w + cell_w / 2, y0 - 35, LABELS[method], 21, weight="bold"))
    for ri, (surface, fixture, fixture_label) in enumerate(fixtures):
        y = y0 + ri * cell_h
        if ri in (0, 3):
            group = "Sphere-patch surface" if surface == "sphere-patch" else "Voxel-derived watertight SES"
            group_y = y + 1.5 * cell_h
            svg.append(f'<text x="28" y="{group_y}" text-anchor="middle" font-family="Arial" font-size="22" font-weight="bold" fill="{COLORS["muted"]}" transform="rotate(-90 28,{group_y})">{escape(group)}</text>')
        svg.append(t(x0 + label_w - 20, y + 61, fixture_label, 24, anchor="end", weight="bold"))
        for ci, method in enumerate(methods):
            matches = [r for r in rows if r["surface_protocol"] == surface and r["fixture"] == fixture and r["preconditioner_label"] == method]
            if len(matches) != 1:
                raise RuntimeError(f"Missing matrix row: {surface} {fixture} {method}")
            iterations = int(float(matches[0]["iterations"]))
            x = x0 + label_w + ci * cell_w
            svg.append(f'<rect x="{x+5}" y="{y+6}" width="{cell_w-10}" height="{cell_h-12}" rx="6" fill="{iteration_fill(iterations)}" stroke="white" stroke-width="3"/>')
            svg.append(t(x + cell_w / 2, y + 52, str(iterations), 29, weight="bold"))
            svg.append(t(x + cell_w / 2, y + 79, "stall" if iterations >= 720 else "pass", 19, color=COLORS["muted"]))
    export(svg, output / "figureS3_iteration_matrix")


def figure4(rows: list[dict[str, str]], output: Path) -> None:
    svg = base_svg(2400, 820, 61)
    methods = list(LABELS)
    panels = [
        ("sphere-patch", "Sphere-patch stress surface"),
        ("grid-ses", "Watertight Grid-SES"),
        ("analytic-ses", "Analytic patch-ray SES"),
    ]
    colors = METHOD_COLORS
    x_min, x_max = 0.5, 500.0
    legend = [("Total wall time", "circle"), ("Preconditioner setup", "square")]
    for i, (label, shape) in enumerate(legend):
        x = 900 + i * 420
        if shape == "circle":
            svg.append(f'<circle cx="{x}" cy="38" r="8" fill="#333"/>')
        else:
            svg.append(f'<rect x="{x-8}" y="30" width="16" height="16" fill="white" stroke="#333" stroke-width="3"/>')
        svg.append(t(x + 22, 48, label, 28, anchor="start"))
    for pi, (surface, title) in enumerate(panels):
        px0 = 275 + pi * 710
        px1, py0, py1 = px0 + 570, 145, 700
        svg.append(t(px0, 108, chr(65 + pi), 33, anchor="start", weight="bold"))
        svg.append(t(px0 + 48, 108, title, 29, anchor="start", weight="bold"))
        for tick in (1, 10, 100, 500):
            x = px0 + math.log10(tick / x_min) / math.log10(x_max / x_min) * (px1 - px0)
            svg.append(f'<line x1="{x}" y1="{py0}" x2="{x}" y2="{py1}" stroke="#DDDDDD"/>')
            svg.append(t(x, py1 + 33, str(tick), 26))
        for mi, method in enumerate(methods):
            y = py0 + 42 + mi * 84
            matches = [r for r in rows if r["surface_protocol"] == surface and r["fixture"] == "1a63" and r["preconditioner_label"] == method]
            if len(matches) != 1:
                raise RuntimeError(f"Missing 1A63 timing row: {surface} {method}")
            row = matches[0]
            total, setup, iterations = float(row["wall_seconds"]), float(row["setup_seconds"]), int(float(row["iterations"]))
            x_total = px0 + math.log10(total / x_min) / math.log10(x_max / x_min) * (px1 - px0)
            svg.append(f'<line x1="{px0}" y1="{y}" x2="{x_total}" y2="{y}" stroke="{colors[method]}" stroke-width="4"/>')
            svg.append(f'<circle cx="{x_total}" cy="{y}" r="8" fill="{colors[method]}"/>')
            if setup > 0:
                x_setup = px0 + math.log10(setup / x_min) / math.log10(x_max / x_min) * (px1 - px0)
                svg.append(f'<rect x="{x_setup-7}" y="{y-7}" width="14" height="14" fill="white" stroke="{colors[method]}" stroke-width="3"/>')
            if pi == 0:
                svg.append(t(px0 - 18, y + 8, LABELS[method], 26, anchor="end"))
            suffix = "; stall" if iterations >= 720 else f"; {iterations} iter"
            if iterations >= 720:
                svg.append(t(x_total - 10, y - 15, f"{total:.1f} s; stall", 22, anchor="end", color=COLORS["muted"]))
            elif x_total > px1 - 150:
                svg.append(t(x_total - 10, y - 15, f"{total:.1f} s{suffix}", 22, anchor="end", color=COLORS["muted"]))
            else:
                svg.append(t(x_total + 10, y + 8, f"{total:.1f} s{suffix}", 22, anchor="start", color=COLORS["muted"]))
        svg.append(t((px0 + px1) / 2, 777, "Time (s, log scale)", 27))
    export(svg, output / "figure4_1a63_setup_total_time")


def main() -> None:
    parsed = args()
    with parsed.matrix.open(newline="", encoding="utf-8-sig") as handle:
        rows = list(csv.DictReader(handle))
    rows = [
        r
        for r in rows
        if r["fixture"] in {"barnase", "1a63", "1sh1"}
        and r["surface_protocol"] in {"sphere-patch", "grid-ses", "analytic-ses"}
    ]
    figure1(parsed.output_dir)
    figure2(parsed.output_dir)
    figure3(parsed.output_dir)
    figure_s3(rows, parsed.output_dir)
    figure4(rows, parsed.output_dir)


if __name__ == "__main__":
    main()
