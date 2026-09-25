#!/usr/bin/env python3
"""Does the rQTY substitution set survive a change of secondary-structure source?

The design rule converts a glutamine, threonine or tyrosine if and only if its
DSSP assignment is helical.  In the manuscript that assignment comes from the
top-ranked AlphaFold3 model of full-length mature ApoE (pTM 0.53).  This script
repeats the assignment on two experimental structures of the same protein and
reports, position by position and helix by helix, whether the rule would have
made the same decisions.

  1LPE  Wilson et al. 1991, X-ray 2.25 A, wild-type ApoE3 N-terminal domain
  2L7B  Chen, Li & Wang 2011, solution NMR, full-length monomerising ApoE3
        variant, 20 conformers

Written for Supplementary Tables S4 and S5.  Output: _analysis/dssp_xtal.json.

Two guards run on every file before anything is compared.  The mature-numbering
offset is derived by exact substring match of the observed sequence against the
native design segment, never assumed, and the segment itself is then required to
match that sequence residue for residue; a file that fails either check is
rejected rather than silently misaligned.  DSSP is computed with the same
implementation for all three sources, so the comparison is internally
consistent even though it is not the reference Kabsch-Sander program.
"""
import json, os
import numpy as np
import mdtraj as md

HERE = os.path.dirname(os.path.abspath(__file__))

def _pick(*cands):
    """First candidate that exists; the last is the documented default.
    Lets the scripts run unchanged from the project folder or from the
    published repository, and any of them can still be overridden by the
    matching environment variable."""
    for c in cands:
        if c and os.path.isdir(c):
            return c
    return cands[-1]

ROOT     = os.environ.get('APOE_ROOT') or os.path.normpath(os.path.join(HERE, '..', '..'))
ANALYSIS = _pick(os.environ.get('ANALYSIS_DIR'),
                 ANALYSIS,
                 os.path.join(ROOT, 'data', 'analysis_json'))
AF3DIR   = _pick(os.environ.get('AF3_DIR'),
                 os.path.join(ROOT, 'data', 'AF3', 'AlphaFold3'),
                 os.path.join(ROOT, 'data', 'AF3'))
STRUCTD  = _pick(os.environ.get('STRUCT_DIR'),
                 os.path.join(ROOT, 'data', 'Structures'))
MDDIR    = _pick(os.environ.get('MD_DIR'),
                 os.path.join(ROOT, 'data', 'MD'))
OUT = ANALYSIS; os.makedirs(OUT, exist_ok=True)
STRUCT = STRUCTD
FIRST, LAST = 11, 167
HELICAL = set('HGI')
NMR_THRESH = 0.5           # a residue is helical when most conformers say so
MIN_SEG = 4                # helical runs shorter than this are not reported

# the reference assignment: the six DSSP helices of the design model
REF_HEL = {'H0': (11, 19), 'H1': (25, 41), 'H1b': (45, 52),
           'H2': (55, 81), 'H3': (87, 124), 'H4': (131, 167)}

T = json.load(open(os.path.join(HERE, 'si_tables.json')))
NAT = T['NAT']                                # native segment, mature 11-167
DESIGN = {}                                   # straight out of Supplementary Table S1
for row in T['S1'][1:]:
    DESIGN[int(row[0])] = dict(res=row[1][0], dssp_af3=row[2], helix=row[3],
                               converted=(row[4] == 'yes'))

AA3 = {'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C', 'GLN': 'Q',
       'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
       'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S', 'THR': 'T', 'TRP': 'W',
       'TYR': 'Y', 'VAL': 'V', 'MSE': 'M', 'SEC': 'U'}
show = lambda c: ('C' if c == ' ' else str(c))


def chain_seq(chain):
    seq, nums, res = '', [], []
    for r in chain.residues:
        c = AA3.get(r.name.upper())
        if c is None:
            continue
        seq += c; nums.append(r.resSeq); res.append(r)
    return seq, nums, res


def find_offset(seq, nums):
    """mature = author + offset, referenced to the author number of the match.

    Some depositions already number in mature ApoE coordinates, others start at
    1, and some carry an expression tag in front of the native sequence.
    """
    for lead in range(0, min(40, len(seq))):
        probe = seq[lead:lead + 60]
        if len(probe) < 30:
            break
        k = NAT.find(probe)
        if k != -1:
            return (FIRST + k) - nums[lead], lead
    return None, None


def analyse(path, key, label):
    t = md.load(path)
    t = t.atom_slice(t.topology.select('protein'))
    for ch in t.topology.chains:
        seq, nums, res = chain_seq(ch)
        if len(seq) < 30:
            continue
        off, lead = find_offset(seq, nums)
        if off is not None:
            break
    else:
        raise SystemExit(f'{key}: no chain aligns to the design segment')

    # guard: the segment in this file must BE the native design sequence
    inseg = [(c, n + off) for c, n in zip(seq, nums) if FIRST <= n + off <= LAST]
    lo, hi = inseg[0][1], inseg[-1][1]
    got = ''.join(c for c, _ in inseg)
    want = NAT[lo - FIRST:hi - FIRST + 1]
    if got != want:
        bad = [lo + i for i, (a, b) in enumerate(zip(got, want)) if a != b]
        raise SystemExit(f'{key}: segment {lo}-{hi} is not the native sequence '
                         f'(differs at mature {bad[:12]})')

    ss = md.compute_dssp(t, simplified=False)
    nf = t.n_frames
    cols = {r.resSeq + off: r.index for r in res}
    name = {r.resSeq + off: AA3[r.name.upper()] for r in res}

    # per-residue helical fraction over the conformers, for helix boundaries
    frac, code, spread = {}, {}, {}
    for mat, col in cols.items():
        c = ss[:, col]
        frac[mat] = float(np.mean(np.isin(c, list(HELICAL))))
        u, n = np.unique(c, return_counts=True)
        code[mat] = show(u[n.argmax()])
        spread[mat] = {show(a): int(b) for a, b in zip(u, n)}

    # helical segments: contiguous runs of residues above threshold
    segs, ordered = [], sorted(m for m in cols if FIRST <= m <= LAST)
    i = 0
    while i < len(ordered):
        if frac[ordered[i]] >= NMR_THRESH:
            j = i
            while (j + 1 < len(ordered) and frac[ordered[j + 1]] >= NMR_THRESH
                   and ordered[j + 1] == ordered[j] + 1):
                j += 1
            if ordered[j] - ordered[i] + 1 >= MIN_SEG:
                segs.append((ordered[i], ordered[j]))
            i = j + 1
        else:
            i += 1

    print(f'{key}: {nf} frame(s), chain covers mature {lo}-{hi} '
          f'(author {off:+d}), segment sequence verified')
    print(f'       helical segments: {segs}')
    return dict(key=key, label=label, n_frames=nf, covers=[lo, hi], offset=off,
                segments=segs, code=code, frac=frac, spread=spread, name=name)


def main():
    A = analyse(os.path.join(STRUCT, '1LPE.cif'), 'xray',
                'X-ray, wild-type ApoE3 N-terminal domain')
    B = analyse(os.path.join(STRUCT, '2L7B.cif'), 'nmr',
                'Solution NMR, full-length ApoE3, 20 conformers')

    # ---- helix boundaries, one row per reference segment ----
    def match(src, a, b):
        """the source's helical run overlapping the reference segment"""
        if not any(a <= m <= b for m in src['code']):
            return None
        hits = [s for s in src['segments'] if s[0] <= b and s[1] >= a]
        if not hits:
            return 'none'
        return max(hits, key=lambda s: min(s[1], b) - max(s[0], a))

    helix_rows = []
    for nm, (a, b) in REF_HEL.items():
        row = {'segment': nm, 'af3': f'{a}–{b}'}
        for src in (A, B):
            m = match(src, a, b)
            row[src['key']] = ('–' if m is None else
                               'none' if m == 'none' else f'{m[0]}–{m[1]}')
        helix_rows.append(row)

    # ---- one row per glutamine, threonine and tyrosine in the segment ----
    pos_rows, unanimous, changed = [], 0, []
    for pos in sorted(DESIGN):
        d = DESIGN[pos]
        row = dict(pos=pos, res=d['res'], helix=d['helix'],
                   converted=d['converted'], af3=d['dssp_af3'])
        dissent = []
        for src in (A, B):
            k = src['key']
            if pos not in src['code']:
                row[k] = '–'
                row[k + '_conv'] = None
                continue
            assert src['name'][pos] == d['res'], f'{k} residue mismatch at {pos}'
            c = src['code'][pos]
            n = src['spread'][pos][c]
            row[k] = c if src['n_frames'] == 1 else f'{c} ({n}/{src["n_frames"]})'
            wc = c in HELICAL
            row[k + '_conv'] = wc
            if wc != d['converted']:
                dissent.append('1LPE' if k == 'xray' else '2L7B')
        row['dissent'] = dissent
        if not dissent:
            unanimous += 1
        else:
            changed.append(pos)
        pos_rows.append(row)

    seen_x = [r for r in pos_rows if r['xray'] != '–']
    agree_x = sum(1 for r in seen_x if '1LPE' not in r['dissent'])
    agree_n = sum(1 for r in pos_rows if '2L7B' not in r['dissent'])
    both = [r['pos'] for r in pos_rows if len(r['dissent']) == 2]

    print()
    print(f'  positions                 : {len(pos_rows)}')
    print(f'  unanimous across sources  : {unanimous}')
    print(f'  1LPE agrees               : {agree_x} of {len(seen_x)} resolved')
    print(f'  2L7B agrees               : {agree_n} of {len(pos_rows)}')
    print(f'  both experimental dissent : {both}')
    print(f'  any dissent               : {changed}')
    for r in pos_rows:
        if r['dissent']:
            print(f"    {r['res']}{r['pos']:<4} converted={str(r['converted']):<5} "
                  f"AF3={r['af3']:<2} 1LPE={r['xray']:<10} 2L7B={r['nmr']:<10} "
                  f"dissent={','.join(r['dissent'])}")

    json.dump(dict(helix_rows=helix_rows, pos_rows=pos_rows,
                   summary=dict(n_positions=len(pos_rows), unanimous=unanimous,
                                xray_agree=agree_x, xray_resolved=len(seen_x),
                                nmr_agree=agree_n, both_dissent=both,
                                any_dissent=changed),
                   sources={s['key']: {k: s[k] for k in
                                       ('label', 'n_frames', 'covers', 'offset', 'segments')}
                            for s in (A, B)}),
              open(os.path.join(OUT, 'dssp_xtal.json'), 'w'), indent=1)
    print(f'\n  wrote {os.path.join(OUT, "dssp_xtal.json")}')


if __name__ == '__main__':
    main()
