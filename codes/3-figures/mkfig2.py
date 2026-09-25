#!/usr/bin/env python3
"""Figure 2 -- the designed fold, drawn from structures.

(a) the 22 converted positions as side-chain sticks on the native bundle
(b) native and rQTY top-ranked AlphaFold3 models superposed on the bundle
(c) the 1LPE crystal structure superposed on the native model

The three ribbons are ray traced by mkstruct.py, which derives and verifies the
mature-numbering offset of every file and cross-checks each superposition RMSD
against an independent Kabsch fit.  This script only crops them to their ink,
lays them out and adds the keys, so the numbers quoted in the panel labels come
from struct.json rather than from anything typed here.
"""
import json, os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
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
SDIR = os.environ.get('STRUCTOUT') or os.path.join(ROOT, 'figures', '_struct')
S = json.load(open(os.path.join(SDIR, 'struct.json')))
ARMC = '#7FA98A'
GREY = '#8A8A86'
F.style()


def ink(name):
    """Load a render and crop it to its own ink, so the three panels scale alike."""
    a = mpimg.imread(os.path.join(SDIR, name + '.png'))
    rgb = a[..., :3] if a.ndim == 3 else a
    mask = (rgb < 0.96).any(-1)
    r, c = np.where(mask)
    pad = 6
    return a[max(0, r.min() - pad):r.max() + pad, max(0, c.min() - pad):c.max() + pad]


fig = plt.figure(figsize=(7.2, 5.0))
# the panel band stops at 0.27 so the two keys have a clear strip beneath it;
# at a smaller bottom margin the 'ribbons' title collided with panel b's caption
gs = fig.add_gridspec(1, 3, wspace=0.06, left=0.02, right=0.98,
                      top=0.985, bottom=0.27)

PANELS = [
    ('struct_subs', 'a',
     'native fold, 22 converted positions'),
    ('struct_super', 'b',
     f'native and rQTY models, {S["super_rmsd"]:.2f} ' + 'Å over the bundle'),
    ('struct_xtal', 'c',
     f'crystal structure and native model, {S["xtal_rmsd"]:.2f} ' + 'Å'),
]
for k, (name, tag, sub) in enumerate(PANELS):
    ax = fig.add_subplot(gs[0, k])
    ax.imshow(ink(name))
    ax.axis('off')
    ax.text(-0.02, 1.0, tag, transform=ax.transAxes, fontsize=10,
            fontweight='bold', va='top', ha='left', color=F.INK)
    ax.text(0.5, -0.035, sub, transform=ax.transAxes, ha='center', va='top',
            fontsize=6.5, color=F.MUTED)

# ---- keys, one row under the panels ----
sub_keys = [Line2D([], [], color=F.SUB_COLOR[k], lw=2.6,
                   label=f'{k}→{v}  ({n})')
            for k, v, n in [('Q', 'L', 14), ('T', 'V', 4), ('Y', 'F', 4)]]
struct_keys = [Patch(facecolor=F.NAT, edgecolor='none', label='native'),
               Patch(facecolor=F.RQ, edgecolor='none', label='rQTY variant'),
               Patch(facecolor=ARMC, edgecolor='none', label='N-terminal arm (11–24)'),
               Patch(facecolor=F.LDLR, edgecolor='none', label='LDLR region (136–150)'),
               Patch(facecolor=GREY, edgecolor='none', label='crystal structure (1LPE)')]

l1 = fig.legend(handles=sub_keys, loc='upper left', bbox_to_anchor=(0.045, 0.155),
                ncol=1, fontsize=6.6, handlelength=1.5, columnspacing=1.4,
                frameon=False, title='substitutions (sticks)',
                title_fontsize=6.6, alignment='left')
fig.add_artist(l1)
fig.legend(handles=struct_keys, loc='upper left', bbox_to_anchor=(0.30, 0.155),
           ncol=2, fontsize=6.6, handlelength=1.1, columnspacing=1.6,
           frameon=False, title='ribbons', title_fontsize=6.6, alignment='left')

for p in F.save(fig, 'Figure2'):
    print('Figure2 ->', os.path.normpath(p))
