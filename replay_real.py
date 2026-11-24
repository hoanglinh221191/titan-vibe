"""Run the five systems directly from the 11 September deposited data archive.

The reconstructed atom records are compared with the deposited outputs at PDB precision. The
scientific report fields are compared with validation/REAL367_FIELDS_v1.1.0.csv, the outputs of release
1.1.0, which later releases must reproduce; the fields that differ from the deposited 1.0.3 reports are
counted separately (release 1.1.0 labels contact networks that close through the boundary as
extended_periodic).

  python replay_real.py --data-root <vibe_data_update_20260911> --output <new dir> [--all-frames]
  python replay_real.py ... --write-fields <csv>   # record the fields of this run (release builder)
"""
from pathlib import Path
import argparse, csv, json, os, subprocess, time, hashlib

FIELDS = ('status', 'mode', 'component_count', 'evidence_state', 'fallback_reason', 'hard_feasible',
          'hard_certified', 'evidence_graph_connected', 'soft_observed_contacts_checked',
          'soft_observed_selected_compatible_contact_count', 'soft_observed_hypothesis_count',
          'soft_observed_compatible_hypothesis_count', 'soft_observed_ambiguous_relation_count',
          'soft_observed_no_support_relation_count', 'soft_observed_alternative_hypothesis_contact_count',
          'feasible_assignment_set_enumerated', 'feasible_assignments_within_domain', 'wrapped_coordinate_hash',
          'contact_cycle_rank', 'contact_graph_components', 'self_image_units', 'input_overlap_pairs')
# Recorded, not compared: the search route name (equal outputs can be reached by different routes) and
# the input checksum, which follows the last bits of the cell matrix; the coordinates are compared.
RECORDED = ('mode', 'wrapped_coordinate_hash')


def default_binary():
    """titan-vibe.exe in this package on Windows; source/titan-vibe (the Linux build) elsewhere."""
    root = Path(__file__).resolve().parent
    names = ['titan-vibe.exe'] if os.name == 'nt' else ['titan-vibe']
    candidates = [root / n for n in names] + [root / 'source' / n for n in names]
    if os.name != 'nt':
        candidates.reverse()
    return next((c for c in candidates if c.is_file()), candidates[0])


def coordinates(path, first_only):
    rows = []
    for line in path.read_text().splitlines():
        if line.startswith(('ATOM  ', 'HETATM')):
            rows.append(line[:54])
        elif line.startswith('ENDMDL') and rows and first_only:
            break
    return rows


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--data-root', required=True, type=Path,
                   help='Extracted vibe_data_update_20260911 directory containing INPUT_PATH_MAP.json')
    p.add_argument('--output', required=True, type=Path, help='New directory; existing directories are refused')
    p.add_argument('--all-frames', action='store_true', help='All 367 comparison frames; default first frame of each system')
    p.add_argument('--only-extras', action='store_true', help='Only the two additional Figure 1 frames')
    p.add_argument('--binary', type=Path, default=default_binary(),
                   help='Default: titan-vibe.exe (Windows) or source/titan-vibe (Linux) in this package')
    p.add_argument('--write-fields', type=Path, help='Write the report fields of this run to this CSV')
    a = p.parse_args()
    here = Path(__file__).resolve().parent
    root = a.data_root.resolve()
    out = a.output.resolve()
    out.mkdir(parents=True, exist_ok=False)
    reference = {}
    ref_csv = here / 'validation' / 'REAL367_FIELDS_v1.1.0.csv'
    if ref_csv.is_file() and not a.write_fields:
        for r in csv.DictReader(ref_csv.open(newline='')):
            reference[(r['system'], r['tag'], int(r['frame_index']))] = r
    bindings = json.loads((root / 'INPUT_PATH_MAP.json').read_text())
    records, written = {}, []
    data = root / 'jcim_revision_controls_20260911'
    jobs = []
    for sid in sorted(bindings):
        if not a.only_extras:
            jobs.append((sid, 'A'))
        if a.all_frames or a.only_extras:
            jobs.extend((sid, d.name) for d in sorted((data / 'A_vibe' / sid).glob('t*')) if d.is_dir())
    for sid, tag in jobs:
        b = bindings[sid]
        archived = data / 'A_vibe' / sid / tag
        target = out / sid / tag
        target.mkdir(parents=True)
        first_only = tag == 'A' and not a.all_frames
        values = {'-s': root / b['pdb'], '-n': root / b['ndx'], '-common_topo': root / b['top'],
                  '-traj': data / 'A_inputs' / sid / f'full_{tag}.xtc'}
        for f in values.values():
            if not f.is_file():
                raise FileNotFoundError(f)
        lines = []
        for line in (archived / 'arg.in').read_text().splitlines():
            k, sep, v = line.partition('=')
            if k in values:
                v = values[k].as_posix()
            elif k == '-e' and first_only:
                v = '0'
            lines.append(k + sep + v)
        (target / 'arg.in').write_text('\n'.join(lines) + '\n')
        t = time.perf_counter()
        with (target / 'stdout.txt').open('wb') as so, (target / 'stderr.txt').open('wb') as se:
            run = subprocess.run([str(a.binary.resolve()), '-a', 'arg.in'], cwd=target, stdout=so, stderr=se,
                                 timeout=1800)
        result = {'exit_code': run.returncode, 'elapsed_s': time.perf_counter() - t}
        name = next(s.split('=', 1)[1] for s in lines if s.startswith('-opdb='))
        if run.returncode == 0 and (target / name).exists():
            x = coordinates(target / name, False)
            y = coordinates(archived / name, first_only)
            result.update({'atom_records': len(x), 'coordinates_identical_at_PDB_precision': bool(x) and x == y})
            current = list(csv.DictReader((target / 'pbctopo_frame_report.csv').open()))
            deposited = list(csv.DictReader((archived / 'pbctopo_frame_report.csv').open()))
            if first_only:
                deposited = deposited[:1]
            result['frames'] = len(current)
            for i, row in enumerate(current):
                written.append({'system': sid, 'tag': tag, 'frame_index': i, **{k: row.get(k, '') for k in FIELDS}})
            changed = {}
            for x_row, y_row in zip(current, deposited):
                for k in FIELDS:
                    if k not in RECORDED and k in y_row and x_row.get(k) != y_row.get(k):
                        changed[k] = changed.get(k, 0) + 1
            result['fields_changed_from_1_0_3'] = changed
            if reference:
                mism = []
                for i, row in enumerate(current):
                    want = reference.get((sid, tag, i))
                    if want is None:
                        mism.append(f'frame {i}: no reference row')
                        continue
                    mism += [f'frame {i}: {k} {row.get(k)!r} != {want[k]!r}' for k in FIELDS
                             if k not in RECORDED and row.get(k, '') != want[k]]
                result['field_mismatches'] = mism[:20]
                result['scientific_fields_identical'] = not mism
            result['route_changes'] = sum(x_row.get('mode') != y_row.get('mode') for x_row, y_row in zip(current, deposited))
            result['input_checksum_changes'] = sum(x_row.get('wrapped_coordinate_hash') != y_row.get('wrapped_coordinate_hash')
                                                   for x_row, y_row in zip(current, deposited))
        result['pass'] = (run.returncode == 0 and result.get('coordinates_identical_at_PDB_precision', False)
                          and (a.write_fields is not None or result.get('scientific_fields_identical', False)))
        records[sid + '/' + tag] = result
        (out / 'REAL_REPLAY_RESULT.json').write_text(json.dumps(
            {'binary_sha256': hashlib.sha256(a.binary.read_bytes()).hexdigest(), 'all_frames_requested': a.all_frames,
             'all_pass': all(v['pass'] for v in records.values()), 'systems_completed': len(records),
             'results': records}, indent=2))
        print(sid, tag, {k: v for k, v in result.items() if k != 'field_mismatches'}, flush=True)
    if a.write_fields:
        with a.write_fields.open('w', newline='') as fh:
            w = csv.DictWriter(fh, fieldnames=['system', 'tag', 'frame_index', *FIELDS], lineterminator='\n')
            w.writeheader()
            w.writerows(written)
    return 0 if all(x['pass'] for x in records.values()) else 1


if __name__ == '__main__':
    raise SystemExit(main())
