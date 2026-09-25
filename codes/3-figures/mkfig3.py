#!/usr/bin/env python3
"""Figure 3 -- AlphaFold3 comparison of the native and rQTY segments.

(a) per-residue pLDDT of the top-ranked model of each sequence
(b) per-residue deviation of each rQTY model from the top native model
(c) pairwise Ca RMSD within and between the two ensembles
(d) total and apolar solvent-accessible surface area
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
FIRST = 11
F.style()

fig = plt.figure(figsize=(7.2, 5.5))
gs = fig.add_gridspec(2, 2, hspace=0.55, wspace=0.34)

# ---------------- (a) pLDDT ----------------
ax = fig.add_subplot(gs[0, 0])
pn = np.array(A['plddt_nat']); pr = np.array(A['plddt_rqty'])
r = np.arange(FIRST, FIRST + len(pn))
ax.axvspan(F.ARM[0], F.ARM[1], color=F.BAND, lw=0, zorder=0)
for y in (70, 90):
    ax.axhline(y, color=F.RULE, lw=0.7, ls=':', zorder=1)
ax.plot(r, pn, color=F.NAT, lw=1.3, zorder=3, label='native')
ax.plot(r, pr, color=F.RQ, lw=1.3, zorder=4, label='rQTY')
ax.set_xlim(FIRST, FIRST + len(pn) - 1); ax.set_ylim(20, 100)
ax.set_xlabel('Residue', fontsize=7.5); ax.set_ylabel('pLDDT', fontsize=7.5)
ax.legend(fontsize=6.8, loc='lower right', handlelength=1.3)
ax.text(17.5, 26, 'arm', fontsize=6.3, color=F.MUTED, ha='center')
F.panel(ax, 'a')

# ---------------- (b) per-residue deviation ----------------
ax = fig.add_subplot(gs[0, 1])
D = np.array(A['dev'])
ax.axvspan(F.ARM[0], F.ARM[1], color=F.BAND, lw=0, zorder=0)
for row in D:
    ax.plot(r, row, color=F.RQ, lw=0.6, alpha=0.45, zorder=2)
ax.plot(r, D.mean(0), color=F.RQ, lw=1.6, zorder=3, label='mean of 5 models')
ax.set_xlim(FIRST, FIRST + len(pn) - 1)
ax.set_xlabel('Residue', fontsize=7.5)
ax.set_ylabel('C$\\alpha$ deviation from\ntop native model (Å)', fontsize=7.5)
ax.legend(fontsize=6.8, loc='upper right', handlelength=1.3)
F.panel(ax, 'b')

# ---------------- (c) pairwise RMSD ----------------
ax = fig.add_subplot(gs[1, 0])
groups = [('within\nnative', 'nat'), ('within\nrQTY', 'rqty'), ('native vs\nrQTY', 'cross')]
rng = np.random.default_rng(0)
for i, (lab, key) in enumerate(groups):
    for tag, mk, fc, off in [('full', 'o', None, -0.16), ('bundle', 'D', 'white', 0.16)]:
        v = np.array(A[tag][key])
        col = F.NAT if key == 'nat' else (F.RQ if key == 'rqty' else F.MUTED)
        ax.scatter(i + off + rng.uniform(-0.05, 0.05, len(v)), v, s=11, marker=mk,
                   facecolor=(fc if fc else col), edgecolor=col, linewidth=0.8, zorder=3)
        ax.plot([i + off - 0.11, i + off + 0.11], [v.mean()] * 2, color=col, lw=1.6, zorder=4)
ax.set_xticks(range(3)); ax.set_xticklabels([g[0] for g in groups], fontsize=7)
ax.set_ylabel('Pairwise C$\\alpha$ RMSD (Å)', fontsize=7.5)
ax.set_yscale('log'); ax.set_ylim(0.1, 12)
ax.set_yticks([0.1, 0.3, 1, 3, 10]); ax.set_yticklabels(['0.1', '0.3', '1', '3', '10'], fontsize=7)
hs = [plt.Line2D([], [], ls='none', marker='o', ms=4, color=F.MUTED, label='all 157 residues'),
      plt.Line2D([], [], ls='none', marker='D', ms=4, mfc='white', mec=F.MUTED, label='bundle 25–167')]
ax.legend(handles=hs, fontsize=6.6, loc='upper left', handletextpad=0.3)
F.panel(ax, 'c')

# ---------------- (d) SASA ----------------
ax = fig.add_subplot(gs[1, 1])
lab = ['total', 'apolar']
xi = np.arange(2)
for j, (key, col, name) in enumerate([('sasa_nat', F.NAT, 'native'), ('sasa_rqty', F.RQ, 'rQTY')]):
    for i, k in enumerate(lab):
        v = np.array(A[key][k]); off = -0.19 + 0.38 * j
        ax.bar(i + off, v.mean(), 0.36, color=col, lw=0, zorder=2,
               label=name if i == 0 else None)
        ax.scatter([i + off] * len(v), v, s=8, color='white', edgecolor=col,
                   linewidth=0.7, zorder=3)
ax.set_xticks(xi); ax.set_xticklabels(['Total SASA', 'Apolar SASA'], fontsize=7.5)
ax.set_ylabel('Solvent-accessible surface (Å$^2$)', fontsize=7.5)
ax.legend(fontsize=7, loc='upper right', handlelength=1.2)
an = np.mean(A['sasa_nat']['apolar']); ar = np.mean(A['sasa_rqty']['apolar'])
ax.annotate(f'+{100*(ar-an)/an:.0f}%', xy=(1.19, ar), xytext=(1.19, ar + 900),
            fontsize=7, color=F.RQ, ha='center',
            arrowprops=dict(arrowstyle='-|>', color=F.RQ, lw=0.9, shrinkA=1, shrinkB=2))
F.panel(ax, 'd')

paths = F.save(fig, 'Figure3')
for _p in paths:
    print('Figure3 ->', os.path.normpath(_p))
