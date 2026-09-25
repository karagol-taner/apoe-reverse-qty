#!/usr/bin/env python3
"""Figure 5 -- membrane engagement of the embedded rQTY variant.

(a) rendered snapshot in the neuronal bilayer
(b) leaflet phosphate planes against lateral distance from the protein
(c) residues within 4 A of a lipid heavy atom against time
(d) first-shell lipid enrichment
"""
import json, os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
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
D = {k: json.load(open(os.path.join(ANALYSIS, f'{k}.json')))
     for k in ['rqty_bilayer', 'rqty_hmmm']}
SNAP = os.environ.get('SNAPSHOT') or os.path.join(HERE, '..', 'apoe_rqty.png')
BINS = ['0-12', '12-18', '18-24', '24-30', '30-40', '40-60']
CENT = [6, 15, 21, 27, 35, 50]
F.style()


def broken(t, y, gap=2.0):
    t = np.asarray(t, float); y = np.asarray(y, float)
    ot, oy = [t[0]], [y[0]]
    for i in range(1, len(t)):
        if t[i] - t[i - 1] > gap:
            ot.append(np.nan); oy.append(np.nan)
        ot.append(t[i]); oy.append(y[i])
    return np.array(ot), np.array(oy)


fig = plt.figure(figsize=(7.2, 5.6))
gs = fig.add_gridspec(2, 3, hspace=0.46, wspace=0.46, width_ratios=[1.25, 1, 1])

# ---------------- (a) snapshot ----------------
ax = fig.add_subplot(gs[0, 0])
im = mpimg.imread(SNAP)
h, w = im.shape[:2]
ax.imshow(im[int(0.02 * h):int(0.93 * h), int(0.03 * w):int(0.97 * w)])
ax.axis('off')
ax.text(-0.06, 1.05, 'a', transform=ax.transAxes, fontsize=10, fontweight='bold',
        va='top', color=F.INK)
ax.text(0.5, -0.045, 'rQTY variant in the neuronal bilayer, 50 ns', transform=ax.transAxes,
        ha='center', va='top', fontsize=6.6, color=F.MUTED)

# ---------------- (b) phosphate-plane profile ----------------
ax = fig.add_subplot(gs[0, 1:])
for key, lab, ls, mk in [('rqty_bilayer', 'neuronal bilayer', '-', 'o'),
                         ('rqty_hmmm', 'HMMM', '--', 's')]:
    pr = D[key]['prof']
    up, lo, xs = [], [], []
    for b, c in zip(BINS, CENT):
        v = np.array(pr[b], float)
        if np.all(np.isnan(v)):
            continue
        half = np.nanmean(v) / 2.0
        up.append(half); lo.append(-half); xs.append(c)
    ax.plot(xs, up, color=F.RQ, ls=ls, lw=1.4, marker=mk, ms=3.4, label=lab)
    ax.plot(xs, lo, color=F.RQ, ls=ls, lw=1.4, marker=mk, ms=3.4)
ax.axhline(0, color=F.RULE, lw=0.8)
ax.set_xlabel('Lateral distance from protein (Å)', fontsize=7.5)
ax.set_ylabel('Phosphate plane position\nrelative to midplane (Å)', fontsize=7.5)
ax.legend(fontsize=6.6, loc='center left', bbox_to_anchor=(0.02, 0.42), handlelength=2.0)
pb = D['rqty_bilayer']['prof']
inner = np.nanmean(np.array(pb['12-18'], float))
far = np.nanmean(np.array(pb['40-60'], float))
ax.annotate('', xy=(15, inner / 2), xytext=(15, -inner / 2),
            arrowprops=dict(arrowstyle='<->', color=F.RQ, lw=0.9))
ax.text(16.5, 6.5, f'{inner:.1f} Å', fontsize=6.3, color=F.RQ, va='center')
ax.annotate('', xy=(50, far / 2), xytext=(50, -far / 2),
            arrowprops=dict(arrowstyle='<->', color=F.RQ, lw=0.9))
ax.text(48.0, 6.5, f'{far:.1f} Å', fontsize=6.3, color=F.RQ, va='center', ha='right')
F.panel(ax, 'b', dx=-0.115, dy=1.12)

# ---------------- (c) residues in contact ----------------
ax = fig.add_subplot(gs[1, 0])
for key, ls in [('rqty_bilayer', '-'), ('rqty_hmmm', '--')]:
    t, y = broken(D[key]['time'], D[key]['nresc'])
    ax.plot(t, y, color=F.RQ, ls=ls, lw=1.3)
ax.set_xlabel('Time (ns)', fontsize=7.5)
ax.set_ylabel('Residues in lipid contact\n(< 4 Å, of 157)', fontsize=7.5)
hs = [plt.Line2D([], [], color=F.RQ, ls='-', lw=1.3, label='neuronal bilayer'),
      plt.Line2D([], [], color=F.RQ, ls='--', lw=1.3, label='HMMM')]
ax.legend(handles=hs, fontsize=6.3, loc='lower right', handlelength=1.8)
F.panel(ax, 'c')

# ---------------- (d) first-shell enrichment ----------------
ax = fig.add_subplot(gs[1, 1:])
NAMES = {'POPS': 'PS', 'POPE': 'PE', 'POPC': 'PC', 'POPI': 'PI', 'CHL1': 'Chol', 'PSM': 'SM'}


def enrich(d):
    half = lambda v: np.asarray(v, float)[len(v) // 2:].mean()
    tot = sum(half(v) for v in d['enr'].values())
    tb = sum(d['bulk'].values())
    return {k: (half(v) / tot) / (d['bulk'][k] / tb) for k, v in d['enr'].items()}


eb, eh = enrich(D['rqty_bilayer']), enrich(D['rqty_hmmm'])
order = sorted(eb, key=lambda k: -eb[k])
y = np.arange(len(order))
ax.barh(y + 0.19, [eb[k] for k in order], 0.34, color=F.RQ, lw=0, label='neuronal bilayer')
ax.barh(y - 0.19, [eh[k] for k in order], 0.34, color='white', edgecolor=F.RQ, lw=0.9,
        hatch='////', label='HMMM')
ax.axvline(1, color=F.RULE, lw=0.9)
ax.set_yticks(y); ax.set_yticklabels([NAMES[k] for k in order], fontsize=7.5)
ax.invert_yaxis()
ax.set_xlabel('First-shell enrichment (shell / bulk mole fraction)', fontsize=7.5)
for i, k in enumerate(order):
    if eb[k] < 0.02:
        ax.text(0.035, i + 0.19, '0.00', fontsize=6.2, color=F.MUTED, va='center')
    if eh[k] < 0.02:
        ax.text(0.035, i - 0.19, '0.00', fontsize=6.2, color=F.MUTED, va='center')
ax.legend(fontsize=6.5, loc='lower right', handlelength=1.5)
F.panel(ax, 'd', dx=-0.115, dy=1.12)

paths = F.save(fig, 'Figure5')
for _p in paths:
    print('Figure5 ->', os.path.normpath(_p))
