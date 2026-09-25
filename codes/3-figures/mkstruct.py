#!/usr/bin/env python3
"""Ribbon renders for Figure 2, drawn with PyMOL.

  struct_subs.png   native bundle, the 22 converted positions as sticks
  struct_super.png  native and rQTY top-ranked models superposed on the bundle
  struct_xtal.png   the 1LPE crystal structure superposed on the native model

NUMBERING.  The AlphaFold3 segment models are numbered 1-157 in their own
coordinates, while 1LPE is numbered in mature ApoE coordinates (23-166).  Every
selection here is written in mature numbering and converted per file through an
offset that is derived by exact sequence match against the design segment,
never assumed, and asserted over the whole overlap before anything is drawn.
Getting this wrong shifts the substitution sticks and the receptor-binding
block by ten residues with no visible error.

Every superposition RMSD PyMOL reports is re-derived here by an independent
Kabsch fit over explicitly matched mature residue numbers, and the script stops
if the two disagree, because cmd.align does its own sequence alignment and can
silently pair residues that were not intended.

All panels are ray traced on white with an orthoscopic camera and the same
orientation.  Colours come from figstyle, so the ribbons match the plots.  The
superposition is over the bundle (25-167) only, which is the comparison the
manuscript reports; the N-terminal arm is left free, because where it sits is
exactly what differs between the two ensembles.
"""
import json, os
import numpy as np
import figstyle as F

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
AF3 = AF3DIR
STRUCT = STRUCTD
OUT = os.environ.get('STRUCTOUT') or os.path.join(ROOT, 'figures', '_struct')
os.makedirs(OUT, exist_ok=True)

NAT_CIF = os.path.join(AF3, 'fold_apoe_nat/fold_apoe_nat_model_0.cif')
RQ_CIF = os.path.join(AF3, 'fold_apoe_rqty/fold_apoe_rqty_model_0.cif')
XTAL = os.path.join(STRUCT, '1LPE.cif')

TB = json.load(open(os.path.join(HERE, 'si_tables.json')))
SEQS = {'native': TB['NAT'], 'rQTY': TB['RQ']}        # mature 11-167
FIRST = 11
BUNDLE = (25, 167)
LDLR = (136, 150)
SUBS = {'Q': [16, 17, 41, 46, 48, 55, 58, 81, 98, 101, 117, 123, 156, 163],
        'T': [18, 57, 67, 89],
        'Y': [36, 74, 118, 162]}
ALLSUBS = sorted(sum(SUBS.values(), []))
assert len(ALLSUBS) == 22, len(ALLSUBS)
W, H = 1400, 1500

AA3 = {'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C', 'GLN': 'Q',
       'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
       'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S', 'THR': 'T', 'TRP': 'W',
       'TYR': 'Y', 'VAL': 'V', 'MSE': 'M'}

import pymol
pymol.finish_launching(['pymol', '-qc'])
from pymol import cmd

PALETTE = {'natgreen': F.NAT, 'rqgold': F.RQ, 'subq': F.SUB_COLOR['Q'],
           'subt': F.SUB_COLOR['T'], 'suby': F.SUB_COLOR['Y'],
           'ldlr': F.LDLR, 'xgrey': '#8A8A86', 'armgreen': '#7FA98A'}
NATC, RQC, QC, TC, YC, LDC, GREY, ARMC = ('natgreen', 'rqgold', 'subq', 'subt',
                                          'suby', 'ldlr', 'xgrey', 'armgreen')


def setup():
    cmd.reinitialize()
    for name, h in PALETTE.items():      # reinitialize discards custom colours
        cmd.set_color(name, [int(h[i:i + 2], 16) / 255 for i in (1, 3, 5)])
    for k, v in [('bg_rgb', [1, 1, 1]), ('opaque_background', 1),
                 ('ray_opaque_background', 1), ('orthoscopic', 1),
                 ('antialias', 2), ('cartoon_fancy_helices', 1),
                 ('cartoon_smooth_loops', 1), ('specular', 0.15),
                 ('ambient', 0.22), ('direct', 0.55), ('reflect', 0.35),
                 ('depth_cue', 0), ('ray_shadow', 0), ('stick_radius', 0.28),
                 ('cartoon_oval_length', 1.0)]:
        cmd.set(k, v)


def residues(name):
    got = []
    cmd.iterate(f'{name} and polymer and name CA',
                'got.append((int(resi), resn))', space={'got': got})
    got.sort()
    return [n for n, _ in got], ''.join(AA3.get(r.upper(), 'X') for _, r in got)


def offset_of(name):
    """mature = file_resi + offset, derived and then verified in full."""
    nums, seq = residues(name)
    for which, ref in SEQS.items():
        k = ref.find(seq[:45])
        if k == -1:
            continue
        off = (FIRST + k) - nums[0]
        lo, hi = nums[0] + off, nums[-1] + off
        assert seq == ref[lo - FIRST:hi - FIRST + 1], \
            f'{name}: sequence does not match {which} over {lo}-{hi}'
        print(f'  {name:5s} {len(nums):3d} res, file {nums[0]}-{nums[-1]}, '
              f'mature {lo}-{hi} (offset {off:+d}), verified against {which}')
        return off, lo, hi
    raise SystemExit(f'{name}: matches neither design sequence')


def load(path, name):
    cmd.load(path, name)
    return offset_of(name)


def rng(name, off, lo, hi):
    return f'{name} and resi {lo - off}-{hi - off}'


def lst(name, off, poss):
    return f'{name} and resi ' + '+'.join(str(p - off) for p in poss)


def ca_of(name, off, lo, hi):
    got = []
    cmd.iterate_state(1, f'{name} and polymer and name CA',
                      'got.append((int(resi), x, y, z))', space={'got': got})
    return {r + off: np.array([x, y, z]) for r, x, y, z in got if lo <= r + off <= hi}


def kabsch(P, Q):
    Pc, Qc = P - P.mean(0), Q - Q.mean(0)
    V, _, Wt = np.linalg.svd(Pc.T @ Qc)
    R = V @ np.diag([1, 1, np.sign(np.linalg.det(V @ Wt))]) @ Wt
    return float(np.sqrt(((Pc @ R - Qc) ** 2).sum(1).mean()))


def verify(a, offa, b, offb, lo, hi, label, reported):
    ca, cb = ca_of(a, offa, lo, hi), ca_of(b, offb, lo, hi)
    common = sorted(set(ca) & set(cb))
    ind = kabsch(np.array([ca[i] for i in common]), np.array([cb[i] for i in common]))
    agree = abs(ind - reported) < 0.05
    print(f'  {label}: PyMOL {reported:.2f} A, independent Kabsch over '
          f'{len(common)} matched CA {ind:.2f} A -> {"agree" if agree else "DISAGREE"}')
    assert agree, f'{label}: the superposition does not reproduce'
    return round(ind, 2), len(common)


def orient_on(sel):
    cmd.orient(sel); cmd.turn('z', 90); cmd.zoom('polymer', 8.0)


def render(name):
    from PIL import Image
    path = os.path.join(OUT, name + '.png')
    cmd.ray(W, H); cmd.png(path, dpi=600)
    a = np.array(Image.open(path).convert('RGB'))
    nz = np.where((a < 245).any(2))
    assert not (nz[0].min() == 0 or nz[1].min() == 0
                or nz[0].max() == H - 1 or nz[1].max() == W - 1), \
        f'{name}: ribbon touches the frame edge'
    print(f'  {name:14s} rendered, clear of frame')


res = {}

print('panel a: the 22 substitutions on the native fold')
setup()
offn, _, _ = load(NAT_CIF, 'nat')
cmd.hide('everything'); cmd.show('cartoon', 'nat and polymer')
cmd.color(ARMC, 'nat')
cmd.color(NATC, rng('nat', offn, *BUNDLE))
cmd.color(LDC, rng('nat', offn, *LDLR))
drawn = 0
for col, poss in [(QC, SUBS['Q']), (TC, SUBS['T']), (YC, SUBS['Y'])]:
    sel = lst('nat', offn, poss) + ' and sidechain'
    cmd.show('sticks', sel); cmd.color(col, sel)
    cmd.set('stick_radius', 0.30, sel)
    drawn += cmd.count_atoms(lst('nat', offn, poss) + ' and name CB')
print(f'  {drawn} of 22 substituted residues drawn as side-chain sticks')
assert drawn == 22, drawn
# the residues under the sticks must be the ones the design converted
_, natseq = residues('nat')
ident = ''.join(natseq[p - offn - 1] for p in ALLSUBS)
expect = ''.join('Q' * 14 + 'T' * 4 + 'Y' * 4)
assert sorted(ident) == sorted(expect), f'stick residues are {ident}, not Q/T/Y'
print(f'  residue identities under the sticks: {ident} (14 Q, 4 T, 4 Y)')
orient_on('nat and name CA')
render('struct_subs')

print('\npanel b: native and rQTY superposed on the bundle')
setup()
offn, _, _ = load(NAT_CIF, 'nat')
offr, _, _ = load(RQ_CIF, 'rq')
rms = cmd.align(rng('rq', offr, *BUNDLE) + ' and name CA',
                rng('nat', offn, *BUNDLE) + ' and name CA', cycles=0)[0]
res['super'] = verify('rq', offr, 'nat', offn, *BUNDLE,
                      'bundle 25-167, native vs rQTY', rms)
cmd.hide('everything'); cmd.show('cartoon', 'polymer')
cmd.color(NATC, 'nat'); cmd.color(RQC, 'rq')
orient_on(rng('nat', offn, *BUNDLE) + ' and name CA')
render('struct_super')

print('\npanel c: crystal structure on the native model')
setup()
offn, _, _ = load(NAT_CIF, 'nat')
cmd.load(XTAL, 'xtal'); cmd.remove('xtal and not polymer')
offx, xlo, xhi = offset_of('xtal')
lo, hi = max(BUNDLE[0], xlo), min(BUNDLE[1], xhi)
rmsx = cmd.align(rng('xtal', offx, lo, hi) + ' and name CA',
                 rng('nat', offn, lo, hi) + ' and name CA', cycles=0)[0]
res['xtal'] = verify('xtal', offx, 'nat', offn, lo, hi,
                     f'bundle {lo}-{hi}, crystal vs native model', rmsx)
cmd.hide('everything'); cmd.show('cartoon', 'polymer')
cmd.color(NATC, 'nat'); cmd.color(GREY, 'xtal')
orient_on(rng('nat', offn, lo, hi) + ' and name CA')
render('struct_xtal')

json.dump({'super_rmsd': res['super'][0], 'super_n': res['super'][1],
           'xtal_rmsd': res['xtal'][0], 'xtal_n': res['xtal'][1],
           'xtal_range': [lo, hi]},
          open(os.path.join(OUT, 'struct.json'), 'w'), indent=1)
print(f"\n  superposition {res['super'][0]:.2f} A over {res['super'][1]} CA")
print(f"  crystal       {res['xtal'][0]:.2f} A over {res['xtal'][1]} CA ({lo}-{hi})")
