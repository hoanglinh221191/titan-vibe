#!/usr/bin/env python3
"""Combine per-model PDB2PQR outputs into a whitespace-safe capsid PQR."""

from __future__ import annotations

import argparse
from pathlib import Path


def normalize_atom_line(line: str) -> str:
    if len(line) < 55:
        raise ValueError(f"Malformed PDB2PQR atom line: {line.rstrip()}")
    record = line[0:6].strip()
    serial = line[6:11].strip()
    atom_name = line[12:16].strip()
    residue = line[17:20].strip()
    chain = line[21:22].strip() or "_"
    residue_id = line[22:26].strip()
    x = line[30:38].strip()
    y = line[38:46].strip()
    z = line[46:54].strip()
    tail = line[54:].split()
    if len(tail) < 2:
        raise ValueError(f"Missing charge/radius: {line.rstrip()}")
    return (
        f"{record:<6}{serial} {atom_name} {residue} {chain} {residue_id} "
        f"{x} {y} {z} {tail[0]} {tail[1]}\n"
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--model-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--models", type=int, required=True)
    args = parser.parse_args()

    atom_count = 0
    partial = args.output.with_suffix(args.output.suffix + ".partial")
    with partial.open("w", encoding="ascii", newline="\n") as writer:
        for model in range(1, args.models + 1):
            source = args.model_root / f"model_{model:03d}.pqr"
            if not source.is_file():
                raise FileNotFoundError(f"Missing model PQR: {source}")
            with source.open("r", encoding="ascii", errors="strict") as reader:
                for line in reader:
                    if line.startswith(("ATOM  ", "HETATM")):
                        writer.write(normalize_atom_line(line))
                        atom_count += 1
        writer.write("END\n")
    partial.replace(args.output)
    print(f"combined_atoms={atom_count}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
