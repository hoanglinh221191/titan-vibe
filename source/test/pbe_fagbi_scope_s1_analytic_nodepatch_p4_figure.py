from __future__ import annotations

import argparse
import csv
import math
import shutil
import subprocess
from pathlib import Path
from xml.sax.saxutils import escape


KJ_PER_KCAL = 4.184
BLUE = "#3B6EA5"
ORANGE = "#D95F59"
INK = "#202020"


def number(value: str) -> float:
    try:
        parsed = float(value)
    except (TypeError, ValueError):
        return math.nan
    return parsed if math.isfinite(parsed) else math.nan


def nice_ticks(lo: float, hi: float, target_count: int = 6) -> tuple[float, float, list[float]]:
    if hi <= lo:
        return lo, hi, [lo]
    raw_step = (hi - lo) / max(1, target_count - 1)
    magnitude = 10.0 ** math.floor(math.log10(raw_step))
    normalized = raw_step / magnitude
    if normalized <= 1.0:
        nice = 1.0
    elif normalized <= 2.0:
        nice = 2.0
    elif normalized <= 2.5:
        nice = 2.5
    elif normalized <= 5.0:
        nice = 5.0
    else:
        nice = 10.0
    step = nice * magnitude
    axis_lo = math.floor(lo / step) * step
    axis_hi = math.ceil(hi / step) * step
    count = int(round((axis_hi - axis_lo) / step))
    return axis_lo, axis_hi, [axis_lo + index * step for index in range(count + 1)]


def text(x: float, y: float, value: str, size: int = 24, anchor: str = "start",
         weight: str = "normal", rotate: int | None = None) -> str:
    transform = f' transform="rotate({rotate} {x:.2f} {y:.2f})"' if rotate is not None else ""
    return (
        f'<text x="{x:.2f}" y="{y:.2f}" font-family="Arial" font-size="{size}" '
        f'font-weight="{weight}" text-anchor="{anchor}" fill="{INK}"{transform}>'
        f'{escape(value)}</text>'
    )


def rect(x: float, y: float, width: float, height: float, fill: str) -> str:
    return f'<rect x="{x:.2f}" y="{y:.2f}" width="{width:.2f}" height="{height:.2f}" fill="{fill}"/>'


def line(x1: float, y1: float, x2: float, y2: float, dash: str = "") -> str:
    extra = f' stroke-dasharray="{dash}"' if dash else ""
    return f'<line x1="{x1:.2f}" y1="{y1:.2f}" x2="{x2:.2f}" y2="{y2:.2f}" stroke="{INK}" stroke-width="2"{extra}/>'


def draw_grouped_panel(parts: list[str], panel: tuple[int, int, int, int], letter: str,
                       labels: list[str], first: list[float], second: list[float],
                       ylabel: str, log_scale: bool = False, legend: bool = False) -> None:
    x0, y0, width, height = panel
    left, right, top, bottom = 112, 24, 58, 88
    px, py = x0 + left, y0 + top
    pw, ph = width - left - right, height - top - bottom
    values = [value for value in first + second if math.isfinite(value) and value > 0] if log_scale else [value for value in first + second if math.isfinite(value)]
    if log_scale:
        lo = 10 ** math.floor(math.log10(min(values)))
        hi = 10 ** math.ceil(math.log10(max(values)))
        mapper = lambda value: py + ph * (math.log10(hi) - math.log10(value)) / (math.log10(hi) - math.log10(lo))
        tick_values = [10 ** exponent for exponent in range(int(math.log10(lo)), int(math.log10(hi)) + 1)]
        baseline = py + ph
    else:
        lo = min(min(values), 0.0)
        hi = max(max(values), 0.0)
        lo, hi, tick_values = nice_ticks(lo, hi)
        mapper = lambda value: py + ph * (hi - value) / (hi - lo)
        baseline = mapper(0.0)
    parts.extend([line(px, py, px, py + ph), line(px, baseline, px + pw, baseline)])
    for value in tick_values:
        yy = mapper(value)
        parts.append(line(px - 6, yy, px, yy))
        label = f"{value:g}" if not log_scale else f"{int(value):d}"
        parts.append(text(px - 12, yy + 8, label, 20, "end"))
    group_width = pw / max(1, len(labels))
    bar_width = min(52.0, group_width * 0.25)
    for index, label in enumerate(labels):
        center = px + group_width * (index + 0.5)
        for offset, value, color in ((-bar_width / 2, first[index], BLUE), (bar_width / 2, second[index], ORANGE)):
            if not math.isfinite(value) or (log_scale and value <= 0):
                continue
            yy = mapper(value)
            top_y, bar_h = (yy, baseline - yy) if yy <= baseline else (baseline, yy - baseline)
            parts.append(rect(center + offset - bar_width / 2, top_y, bar_width, max(1.0, bar_h), color))
        parts.append(text(center, py + ph + 34, label, 21, "middle"))
    parts.append(text(x0 + 20, y0 + 35, letter, 34, weight="bold"))
    parts.append(text(x0 + 28, py + ph / 2, ylabel, 23, "middle", rotate=-90))
    if legend:
        lx, ly = px + pw - 255, y0 + 22
        parts.extend([
            rect(lx, ly, 25, 16, BLUE), text(lx + 34, ly + 15, "TITAN p4", 20),
            rect(lx + 132, ly, 25, 16, ORANGE), text(lx + 166, ly + 15, "APBS d257", 20),
        ])


def draw_single_panel(parts: list[str], panel: tuple[int, int, int, int], letter: str,
                      labels: list[str], values: list[float], ylabel: str) -> None:
    x0, y0, width, height = panel
    left, right, top, bottom = 112, 24, 58, 88
    px, py = x0 + left, y0 + top
    pw, ph = width - left - right, height - top - bottom
    _, hi, tick_values = nice_ticks(0.0, max(1.0, max(values) * 1.08), target_count=6)
    mapper = lambda value: py + ph * (hi - value) / hi
    baseline = py + ph
    parts.extend([line(px, py, px, baseline), line(px, baseline, px + pw, baseline)])
    for value in tick_values:
        yy = mapper(value)
        parts.append(line(px - 6, yy, px, yy))
        parts.append(text(px - 12, yy + 8, f"{value:g}", 20, "end"))
    threshold = mapper(1.0)
    parts.append(line(px, threshold, px + pw, threshold, "8 6"))
    group_width = pw / max(1, len(labels))
    bar_width = min(80.0, group_width * 0.42)
    for index, label in enumerate(labels):
        center = px + group_width * (index + 0.5)
        yy = mapper(values[index])
        parts.append(rect(center - bar_width / 2, yy, bar_width, baseline - yy, "#E69F36"))
        parts.append(text(center, py + ph + 34, label, 21, "middle"))
        parts.append(text(center, yy - 10, f"{values[index]:.2f}", 18, "middle"))
    parts.append(text(x0 + 20, y0 + 35, letter, 34, weight="bold"))
    parts.append(text(x0 + 28, py + ph / 2, ylabel, 23, "middle", rotate=-90))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--matrix", required=True, type=Path)
    parser.add_argument("--png", required=True, type=Path)
    parser.add_argument("--svg", required=True, type=Path)
    parser.add_argument("--summary", required=True, type=Path)
    args = parser.parse_args()
    with args.matrix.open(newline="", encoding="utf-8-sig") as handle:
        rows = list(csv.DictReader(handle))

    records: list[dict[str, float | str]] = []
    for fixture in sorted({row["fixture"] for row in rows if row.get("fixture")}):
        titan = [row for row in rows if row.get("fixture") == fixture and row.get("backend") == "titan_guarded_p4" and row.get("status") == "ok"]
        apbs = [row for row in rows if row.get("fixture") == fixture and row.get("backend") == "apbs" and row.get("status") == "ok"]
        if not titan or len(apbs) < 2:
            continue
        t = titan[-1]
        surface = t.get("primary_surface", "").replace("-", "_")
        discretization = t.get("discretization", "").replace("-", "_")
        if surface != "titan_analytic_patch_ray_ses" or discretization != "node_patch":
            raise RuntimeError(f"{fixture}: expected analytic patch-ray SES/node-patch, got {surface}/{discretization}")
        apbs.sort(key=lambda row: int(row["dime"]))
        coarse, fine = apbs[0], apbs[-1]
        coarse_energy = number(coarse["apbs_energy_kj_mol"])
        fine_energy = number(fine["apbs_energy_kj_mol"])
        records.append({
            "fixture": fixture.upper(),
            "titan_energy": number(t["titan_energy_kj_mol"]) / KJ_PER_KCAL,
            "apbs_energy": fine_energy / KJ_PER_KCAL,
            "apbs_change": 100.0 * abs(fine_energy - coarse_energy) / max(abs(fine_energy), 1.0e-30),
            "titan_rss": number(t["peak_rss_mb"]), "apbs_rss": number(fine["peak_rss_mb"]),
            "titan_wall": number(t["wall_seconds"]), "apbs_wall": number(fine["wall_seconds"]),
            "audit_rms": number(t["sampled_audit_rms"]), "panels": number(t["panels"]),
            "iterations": number(t["iterations"]), "apbs_coarse_dime": coarse["dime"], "apbs_fine_dime": fine["dime"],
        })
    if not records:
        raise RuntimeError("No complete analytic-SES/node-patch p4 plus APBS pairs were found")

    labels = [str(record["fixture"]) for record in records]
    parts = ['<svg xmlns="http://www.w3.org/2000/svg" width="1800" height="1400" viewBox="0 0 1800 1400">', '<rect width="1800" height="1400" fill="white"/>']
    panels = [(35, 35, 850, 640), (915, 35, 850, 640), (35, 715, 850, 640), (915, 715, 850, 640)]
    draw_grouped_panel(parts, panels[0], "A", labels, [float(r["titan_energy"]) for r in records], [float(r["apbs_energy"]) for r in records], "ΔGsolv (kcal mol−1)", legend=True)
    draw_single_panel(parts, panels[1], "B", labels, [float(r["apbs_change"]) for r in records], "APBS refinement change (%)")
    draw_grouped_panel(parts, panels[2], "C", labels, [float(r["titan_rss"]) for r in records], [float(r["apbs_rss"]) for r in records], "Peak RSS (MB)", log_scale=True)
    draw_grouped_panel(parts, panels[3], "D", labels, [float(r["titan_wall"]) for r in records], [float(r["apbs_wall"]) for r in records], "Wall time (s)")
    parts.append('</svg>')
    args.svg.parent.mkdir(parents=True, exist_ok=True)
    args.svg.write_text("\n".join(parts), encoding="utf-8")

    columns = list(records[0].keys())
    with args.summary.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader(); writer.writerows(records)

    inkscape = shutil.which("inkscape") or shutil.which("inkscape.exe")
    if not inkscape:
        known = Path(r"C:\Program Files\Inkscape\bin\inkscape.exe")
        inkscape = str(known) if known.exists() else ""
    if not inkscape:
        raise RuntimeError("Inkscape is required to export the Figure S1 PNG")
    subprocess.run([inkscape, str(args.svg), "--export-type=png", f"--export-filename={args.png}", "--export-dpi=600"], check=True)


if __name__ == "__main__":
    main()
