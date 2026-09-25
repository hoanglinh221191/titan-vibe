#!/usr/bin/env python3
"""Create the manuscript scaling companion without third-party Python packages."""

from __future__ import annotations

import argparse
import csv
import math
import subprocess
from pathlib import Path
from xml.sax.saxutils import escape


PRECONDITIONERS = {
    "diagonal": ("Diagonal", "#555555"),
    "leaf_block_jacobi": ("Leaf block Jacobi", "#D55E00"),
    "fmm_near_graph_ras": ("FMM-near RAS", "#0072B2"),
    "fmm_near_graph_ilu0": ("FMM-near ILU(0)", "#009E73"),
}
GRAPH_COLORS = {
    "distance_knn_ilu0": ("Distance-kNN graph", "#CC79A7"),
    "fmm_near_graph_ilu0": ("FMM-near graph", "#009E73"),
}
SURFACE_LABELS = {
    "sphere-patch": "Sphere-patch surface",
    "grid-ses": "Voxel-derived watertight SES",
    "analytic-ses": "Analytic patch-ray SES",
}
FIXTURE_LABELS = {"barnase": "Barnase", "1a63": "1A63", "1sh1": "1SH1"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--matrix", type=Path, required=True)
    parser.add_argument("--output-prefix", type=Path, required=True)
    return parser.parse_args()


def f(row: dict[str, str], name: str) -> float:
    return float(row[name])


def marker(x: float, y: float, surface: str, color: str, size: float = 7.0) -> list[str]:
    if surface == "sphere-patch":
        return [f'<circle cx="{x:.2f}" cy="{y:.2f}" r="{size:.2f}" fill="white" stroke="{color}" stroke-width="3"/>']
    if surface == "grid-ses":
        points = f"{x:.2f},{y-size:.2f} {x-size:.2f},{y+size:.2f} {x+size:.2f},{y+size:.2f}"
    else:
        points = f"{x:.2f},{y-size:.2f} {x+size:.2f},{y:.2f} {x:.2f},{y+size:.2f} {x-size:.2f},{y:.2f}"
    return [f'<polygon points="{points}" fill="{color}" stroke="{color}" stroke-width="3"/>']


def line_key(x: float, y: float, color: str) -> str:
    return f'<line x1="{x}" y1="{y}" x2="{x+38}" y2="{y}" stroke="{color}" stroke-width="6" stroke-linecap="round"/>'


def main() -> None:
    args = parse_args()
    with args.matrix.open(newline="", encoding="utf-8-sig") as handle:
        rows = [
            row for row in csv.DictReader(handle)
            if row["fixture"] in {"barnase", "1a63", "1sh1"}
            and row["surface_protocol"] in SURFACE_LABELS
        ]
    if not rows:
        raise SystemExit("No graph-preconditioner rows found")

    width, height = 1800, 760
    left, right, top, bottom, gap = 125, 55, 175, 100, 120
    panel_width = (width - left - right - gap) / 2
    panel_height = height - top - bottom
    svg: list[str] = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="180mm" height="72mm" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
    ]

    for i, (label, color) in enumerate(PRECONDITIONERS.values()):
        x, y = 55 + i * 425, 35
        svg.append(line_key(x, y, color))
        svg.append(f'<text x="{x+50}" y="{y+8}" font-family="Arial" font-size="23">{escape(label)}</text>')
    surface_x = (55, 500, 985)
    for x, surface in zip(surface_x, SURFACE_LABELS):
        y = 88
        svg += marker(x + 9, y, surface, "#333333", 7)
        svg.append(f'<text x="{x+32}" y="{y+8}" font-family="Arial" font-size="23">{escape(SURFACE_LABELS[surface])}</text>')
    x, y = 1460, 88
    svg.append(line_key(x, y, GRAPH_COLORS["distance_knn_ilu0"][1]))
    svg.append(f'<text x="{x+50}" y="{y+8}" font-family="Arial" font-size="23">kNN ILU(0)</text>')

    panels = [
        ("GMRES iterations", 7.0, 900.0),
        ("Peak RSS per 10^3 panels (MB)", 5.5, 10.0),
    ]
    x_min, x_max = 5.0, 40.0

    for pi, (ylabel, ymin, ymax) in enumerate(panels):
        x0 = left + pi * (panel_width + gap)
        x1, y0, y1 = x0 + panel_width, top, top + panel_height
        svg.append(f'<line x1="{x0}" y1="{y1}" x2="{x1}" y2="{y1}" stroke="black" stroke-width="1.5"/>')
        svg.append(f'<line x1="{x0}" y1="{y0}" x2="{x0}" y2="{y1}" stroke="black" stroke-width="1.5"/>')
        for tick in (10, 20, 30, 40):
            x = x0 + (tick - x_min) / (x_max - x_min) * panel_width
            svg.append(f'<line x1="{x}" y1="{y1}" x2="{x}" y2="{y1+7}" stroke="black"/>')
            svg.append(f'<text x="{x}" y="{y1+31}" text-anchor="middle" font-family="Arial" font-size="20">{tick}</text>')

        if pi == 0:
            yticks = (10, 30, 100, 300, 720)
            def ymap(value: float) -> float:
                return y1 - (math.log10(value) - math.log10(ymin)) / (math.log10(ymax) - math.log10(ymin)) * panel_height
        else:
            yticks = (6, 7, 8, 9, 10)
            def ymap(value: float) -> float:
                return y1 - (value - ymin) / (ymax - ymin) * panel_height

        for tick in yticks:
            y = ymap(float(tick))
            svg.append(f'<line x1="{x0}" y1="{y}" x2="{x1}" y2="{y}" stroke="#dddddd"/>')
            svg.append(f'<text x="{x0-12}" y="{y+7}" text-anchor="end" font-family="Arial" font-size="20">{tick}</text>')

        svg.append(f'<text x="{x0-20}" y="{y0-18}" font-family="Arial" font-size="27" font-weight="bold">{chr(65+pi)}</text>')
        svg.append(f'<text x="{(x0+x1)/2}" y="{height-25}" text-anchor="middle" font-family="Arial" font-size="22">Panels (10^3)</text>')
        svg.append(f'<text x="{x0-66}" y="{(y0+y1)/2}" text-anchor="middle" font-family="Arial" font-size="22" transform="rotate(-90 {x0-66},{(y0+y1)/2})">{escape(ylabel)}</text>')

        if pi == 0:
            offsets = {name: shift for name, shift in zip(PRECONDITIONERS, (-0.45, -0.15, 0.15, 0.45))}
            for key, (_, color) in PRECONDITIONERS.items():
                for surface in SURFACE_LABELS:
                    for row in [r for r in rows if r["preconditioner_label"] == key and r["surface_protocol"] == surface]:
                        xv = f(row, "panels") / 1000.0 + offsets[key]
                        x = x0 + (xv - x_min) / (x_max - x_min) * panel_width
                        y = ymap(f(row, "iterations"))
                        svg += marker(x, y, surface, color)
                        if key == "fmm_near_graph_ilu0" and surface == "analytic-ses":
                            anchor = "end" if row["fixture"] == "1a63" else "start"
                            label_x = x - 10 if anchor == "end" else x + 10
                            svg.append(f'<text x="{label_x}" y="{y-12}" text-anchor="{anchor}" font-family="Arial" font-size="19" fill="#333333">{FIXTURE_LABELS[row["fixture"]]}</text>')
        else:
            metric = "peak_rss_mb"
            for gi, (key, (_, color)) in enumerate(GRAPH_COLORS.items()):
                for surface in SURFACE_LABELS:
                    for row in [r for r in rows if r["preconditioner_label"] == key and r["surface_protocol"] == surface]:
                        xv = f(row, "panels") / 1000.0 + (-0.18 if gi == 0 else 0.18)
                        value = 1000.0 * f(row, metric) / f(row, "panels")
                        x = x0 + (xv - x_min) / (x_max - x_min) * panel_width
                        svg += marker(x, ymap(value), surface, color)

    svg.append("</svg>")
    args.output_prefix.parent.mkdir(parents=True, exist_ok=True)
    svg_path = args.output_prefix.with_suffix(".svg")
    svg_path.write_text("\n".join(svg), encoding="ascii")
    subprocess.run(["inkscape", str(svg_path), "--export-type=pdf", f"--export-filename={args.output_prefix.with_suffix('.pdf')}"] , check=True)
    subprocess.run(["inkscape", str(svg_path), "--export-type=png", f"--export-filename={args.output_prefix.with_suffix('.png')}", "--export-dpi=600"], check=True)


if __name__ == "__main__":
    main()
