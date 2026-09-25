#!/usr/bin/env python3
"""Supplementary Figures S1 and S2.

S1  complete pairwise Ca RMSD matrices for the ten AlphaFold3 models
S2  per-residue relative solvent accessibility, and membrane depth
"""
import json, os
import numpy as np
import matplotlib.pyplot as plt
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
A = json.load(open(os.path.join(ANALYSIS, 'af3.json')))
S = json.load(open(os.path.join(ANALYSIS, 'rsa.json')))
M = {k: json.load(open(os.path.join(ANALYSIS, f'{k}.json')))
     for k in ['rqty_bilayer', 'rqty_hmmm']}
FIRST = 11
LAB = [f'N{i}' for i in range(1, 6)] + [f'R{i}' for i in range(1, 6)]
F.style()

# ---------------- Figure S1 ----------------
fig, axes = plt.subplots(1, 2, figsize=(8.6, 3.9), gridspec_kw={'wspace': 0.34})
for ax, key, title, vmax in [(axes[0], 'full_mat', 'All 157 residues', 7.0),
                             (axes[1], 'bundle_mat', 'Bundle, residues 25–167', 1.1)]:
    Mx = np.array(A[key])
    im = ax.imshow(Mx, cmap=F.SEQ_GREEN, vmin=0, vmax=vmax)
    ax.set_xticks(range(10)); ax.set_xticklabels(LAB, fontsize=6.5)
    ax.set_yticks(range(10)); ax.set_yticklabels(LAB, fontsize=6.5)
    ax.set_title(title, fontsize=8, pad=6)
    for i in range(10):
        for j in range(10):
            ax.text(j, i, f'{Mx[i, j]:.2f}', ha='center', va='center', fontsize=6.5,
                    color='white' if Mx[i, j] > vmax * 0.55 else F.INK)
    for s in ax.spines.values():
        s.set_visible(True)
    ax.axhline(4.5, color=F.INK, lw=1.1); ax.axvline(4.5, color=F.INK, lw=1.1)
    cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03)
    cb.set_label('C$\\alpha$ RMSD (Å)', fontsize=7)
    cb.ax.tick_params(labelsize=6.5)
for ax, tag in zip(axes, 'ab'):
    ax.text(-0.30, 1.16, tag, transform=ax.transAxes, fontsize=10, fontweight='bold',
            va='top', color=F.INK)
fig.text(0.5, -0.02, 'N1–N5, the five AlphaFold3 models of the native segment; '
                     'R1–R5, the five models of the rQTY variant.',
         ha='center', fontsize=6.6, color=F.MUTED)
paths = F.save(fig, 'FigureS1')
plt.close(fig)
for _p in paths:
    print('FigureS1 ->', os.path.normpath(_p))

# ---------------- Figure S2 ----------------
fig, axes = plt.subplots(2, 1, figsize=(7.2, 5.7), gridspec_kw={'hspace': 0.42})
r = np.arange(FIRST, FIRST + 157)

ax = axes[0]
for p in S['subs']:
    ax.axvline(p, color='#C9C8C3', lw=0.7, zorder=0)
ax.axhline(25, color=F.RULE, lw=0.8, ls='--', zorder=1)
ax.plot(r, S['rsa_nat'], color=F.NAT, lw=1.2, zorder=3, label='native')
ax.plot(r, S['rsa_rqty'], color=F.RQ, lw=1.2, zorder=4, label='rQTY')
ax.set_xlim(FIRST, FIRST + 156)
ax.set_xlabel('Residue', fontsize=7.5)
ax.set_ylabel('Relative solvent\naccessibility (%)', fontsize=7.5)
ax.legend(fontsize=7, loc='upper right', handlelength=1.3)
ax.text(-0.075, 1.13, 'a', transform=ax.transAxes, fontsize=10, fontweight='bold',
        va='top', color=F.INK)

ax = axes[1]
for p in S['subs']:
    ax.axvline(p, color='#F0EFEC', lw=0.7, zorder=0)
for key, ls, lab in [('rqty_bilayer', '-', 'neuronal bilayer'),
                     ('rqty_hmmm', '--', 'HMMM')]:
    ax.plot(r, M[key]['perz'], color=F.RQ, ls=ls, lw=1.2, zorder=3, label=lab)
ax.set_xlim(FIRST, FIRST + 156)
ax.set_xlabel('Residue', fontsize=7.5)
ax.set_ylabel('Mean distance from\nbilayer midplane (Å)', fontsize=7.5)
ax.legend(fontsize=7, loc='upper right', handlelength=1.8)
ax.text(-0.075, 1.13, 'b', transform=ax.transAxes, fontsize=10, fontweight='bold',
        va='top', color=F.INK)

paths = F.save(fig, 'FigureS2')
for _p in paths:
    print('FigureS2 ->', os.path.normpath(_p))
