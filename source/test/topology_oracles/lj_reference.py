"""Intermolecular Lennard-Jones energy in plain Python, the reference the
kernels are compared against (kJ/mol, distances in nm, no periodic images:
TITAN's pair energies use direct distances unless -pbcn is given).

Per-atom parameters (v, w) are sigma/epsilon for combination rules 2 and 3
and C6/C12 for rule 1 (GROMOS). A pair listed in `pairs` (by atom type,
[ nonbond_params ]) takes those values instead of the combination rule."""
import math


def mix(rule, vi, wi, vj, wj):
    if rule == 2:
        return 0.5 * (vi + vj), math.sqrt(wi * wj)
    return math.sqrt(vi * vj), math.sqrt(wi * wj)


def pair_energy(rule, v, w, r):
    if rule == 1:
        r6 = r ** 6
        return w / (r6 * r6) - v / r6
    if v <= 0.0 or w == 0.0:  # charmm::lj_energy
        return 0.0
    sr6 = (v / r) ** 6
    return 4.0 * w * (sr6 * sr6 - sr6)


def distances(xyz_a, xyz_b):
    """Row-major list of |a_i - b_j| for every pair."""
    out = []
    for ax, ay, az in xyz_a:
        for bx, by, bz in xyz_b:
            out.append(math.sqrt((ax - bx) ** 2 + (ay - by) ** 2 + (az - bz) ** 2))
    return out


def interaction(atoms_a, atoms_b, dist, rule, pairs=None, cutoff=None):
    """Sum over a x b. atoms_*: [(type, v, w)]; dist from distances()."""
    pairs = pairs or {}
    total, k = 0.0, 0
    for ta, va, wa in atoms_a:
        for tb, vb, wb in atoms_b:
            r = dist[k]
            k += 1
            if cutoff is not None and r > cutoff:
                continue
            vw = pairs.get((ta, tb))
            v, w = vw if vw is not None else mix(rule, va, wa, vb, wb)
            total += pair_energy(rule, v, w, r)
    return total


def type_pair_energies(atoms_a, atoms_b, dist, rule, cutoff=None):
    """{(type_a, type_b): energy with the combination rule}, to choose the
    pairs an override will change most."""
    out = {}
    k = 0
    for ta, va, wa in atoms_a:
        for tb, vb, wb in atoms_b:
            r = dist[k]
            k += 1
            if cutoff is not None and r > cutoff:
                continue
            v, w = mix(rule, va, wa, vb, wb)
            key = (ta, tb) if ta <= tb else (tb, ta)
            out[key] = out.get(key, 0.0) + pair_energy(rule, v, w, r)
    return out
