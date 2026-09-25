#!/usr/bin/env python3
"""Figure 1 -- design of the reverse-QTY variant.

(a) substitution map along mature residues 11-167
(b) Kyte-Doolittle hydropathy in a nine-residue sliding window
(c) mean GRAVY per helical segment
(d) mean hydrophobic moment per helical segment
"""
import json, os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
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
H = json.load(open(os.path.join(ANALYSIS, 'hydro.json')))
S = json.load(open(os.path.join(ANALYSIS, 'rsa.json')))
NAT, RQ = S['nat'], S['rq']
FIRST = 11
F.style()

fig = plt.figure(figsize=(7.2, 6.6))
gs = fig.add_gridspec(3, 2, height_ratios=[0.78, 1, 0.95], hspace=0.70, wspace=0.34)

# ---------------- (a) substitution map ----------------
ax = fig.add_subplot(gs[0, :])
ax.set_xlim(6, 172); ax.set_ylim(-2.35, 1.15); ax.axis('off')
ax.plot([11, 167], [0, 0], color=F.HELIX_EDGE, lw=2.0, solid_capstyle='butt', zorder=1)
for nm, a, b in F.HEL:
    ax.add_patch(Rectangle((a, -0.26), b - a, 0.52, facecolor=F.HELIX_FILL,
                           edgecolor=F.HELIX_EDGE, lw=0.6, zorder=2))
    tx = 160 if nm == 'H4' else (a + b) / 2      # keep H4 clear of the LDLR block
    ax.text(tx, 0, nm, ha='center', va='center', fontsize=6.5, color=F.MUTED, zorder=3)
ax.add_patch(Rectangle((136, -0.26), 14, 0.52, facecolor=F.LDLR, alpha=0.16,
                       edgecolor=F.LDLR, lw=0.9, zorder=3))
ax.plot([136, 150], [-1.05, -1.05], color=F.LDLR, lw=1.0, solid_capstyle='butt', zorder=4)
ax.text(143, -1.16, 'LDLR region 136–150\n(no substitution)', ha='center', va='top',
        fontsize=6.0, color=F.LDLR, linespacing=1.25)

for p in F.SUBS:
    k = NAT[p - FIRST]
    ax.plot([p, p], [0.28, 0.55], color=F.SUB_COLOR[k], lw=0.7, zorder=4)
    ax.plot([p], [0.62], marker=F.SUB_MARK[k], ms=4.0, color=F.SUB_COLOR[k],
            mec='white', mew=0.5, zorder=5, clip_on=False)
for p in F.NOSUB:
    ax.plot([p], [-0.58], marker='x', ms=4.0, color=F.MUTED, mew=1.1, zorder=5)

ax.plot([11, 167], [-1.78, -1.78], color=F.HELIX_EDGE, lw=0.8, solid_capstyle='butt')
for t in range(20, 161, 20):
    ax.plot([t, t], [-1.78, -1.86], color=F.HELIX_EDGE, lw=0.8)
    ax.text(t, -1.96, str(t), ha='center', va='top', fontsize=6.0, color=F.MUTED)

counts = {k: sum(1 for p in F.SUBS if NAT[p - FIRST] == k) for k in 'QTY'}
hs = [plt.Line2D([], [], ls='none', marker=F.SUB_MARK[k], ms=4.5, color=F.SUB_COLOR[k],
                 mec='white', mew=0.5,
                 label=f'{k}→{"LVF"["QTY".index(k)]} ({counts[k]})') for k in 'QTY']
hs.append(plt.Line2D([], [], ls='none', marker='x', ms=4.5, color=F.MUTED, mew=1.1,
                     label=f'Q/T/Y in loops, not converted ({len(F.NOSUB)})'))
ax.legend(handles=hs, loc='upper center', bbox_to_anchor=(0.5, 1.30), ncol=4,
          fontsize=6.6, handletextpad=0.4, columnspacing=1.6)
F.panel(ax, 'a', dx=-0.045, dy=1.06)

# ---------------- (b) sliding-window hydropathy ----------------
ax = fig.add_subplot(gs[1, :])
x = np.array(H['win_x']); yn = np.array(H['win_nat']); yr = np.array(H['win_rq'])
for nm, a, b in F.HEL:
    ax.axvspan(a, b, color=F.BAND, lw=0, zorder=0)
ax.fill_between(x, yn, yr, where=yr > yn, color=F.RQ, alpha=0.16, lw=0, zorder=1)
ax.axhline(0, color=F.RULE, lw=0.8, zorder=2)
ax.plot(x, yn, color=F.NAT, lw=1.5, zorder=3, label='native')
ax.plot(x, yr, color=F.RQ, lw=1.5, zorder=4, label='rQTY')
ax.set_xlim(x.min() - 1, x.max() + 1)
ax.set_ylim(min(yn.min(), yr.min()) - 0.85, max(yn.max(), yr.max()) + 0.15)
ax.set_xlabel('Residue (mature ApoE numbering)', fontsize=8)
ax.set_ylabel('Kyte–Doolittle GRAVY\n(9-residue window)', fontsize=7.5)
ax.legend(fontsize=7.2, loc='upper left', handlelength=1.4)
nn = int((yn > 0).sum()); nr = int((yr > 0).sum())
ax.text(0.5, 0.015, f'net-hydrophobic windows: {nn}/{len(yn)} native, {nr}/{len(yr)} rQTY',
        transform=ax.transAxes, ha='center', va='bottom', fontsize=6.6, color=F.MUTED)
F.panel(ax, 'b', dx=-0.075, dy=1.13)

# ---------------- (c) mean GRAVY and (d) mean hydrophobic moment ----------------
names = [s['name'] for s in H['seg']]
xi = np.arange(len(names))
for k, (key_n, key_r, lab, tag) in enumerate([
        ('gravy_nat', 'gravy_rq', 'Mean GRAVY', 'c'),
        ('uH_nat', 'uH_rq', 'Mean hydrophobic moment $\\langle\\mu H\\rangle$', 'd')]):
    ax = fig.add_subplot(gs[2, k])
    vn = [s[key_n] for s in H['seg']]; vr = [s[key_r] for s in H['seg']]
    ax.bar(xi - 0.19, vn, 0.36, color=F.NAT, lw=0, label='native')
    ax.bar(xi + 0.19, vr, 0.36, color=F.RQ, lw=0, label='rQTY')
    ax.axhline(0, color=F.RULE, lw=0.8)
    ax.set_xticks(xi); ax.set_xticklabels(names, fontsize=7)
    ax.set_ylabel(lab, fontsize=7.5, labelpad=2)
    if k == 0:
        ax.legend(fontsize=7, loc='upper right', handlelength=1.2)
    F.panel(ax, tag, dx=-0.235, dy=1.20)

paths = F.save(fig, 'Figure1')
for _p in paths:
    print('Figure1 ->', os.path.normpath(_p))
