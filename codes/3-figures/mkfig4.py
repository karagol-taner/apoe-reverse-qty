#!/usr/bin/env python3
"""Figure 4 -- conformational stability in all-atom molecular dynamics.

(a) Ca RMSD of the bundle, (b) radius of gyration, (c) DSSP helix content
against simulation time; (d) per-residue Ca RMSF over the second half.

The bilayer trace is broken across 10-20 ns, where the saved coordinates are
unavailable because that segment was overwritten by an interrupted re-run
(see Supplementary Note S2).  Gaps are inserted rather than interpolated.
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
D = {k: json.load(open(os.path.join(ANALYSIS, f'{k}.json')))
     for k in ['nat_water', 'rqty_bilayer', 'rqty_hmmm']}
FIRST = 11
F.style()

SERIES = [('nat_water',    F.NAT, '-',  'native, water (100 ns)'),
          ('rqty_bilayer', F.RQ,  '-',  'rQTY, neuronal bilayer (50 ns)'),
          ('rqty_hmmm',    F.RQ,  '--', 'rQTY, HMMM (90 ns)')]


def broken(t, y, gap=2.0):
    """Insert NaN wherever the time axis jumps, so the line is not interpolated."""
    t = np.asarray(t, float); y = np.asarray(y, float)
    out_t, out_y = [t[0]], [y[0]]
    for i in range(1, len(t)):
        if t[i] - t[i - 1] > gap:
            out_t.append(np.nan); out_y.append(np.nan)
        out_t.append(t[i]); out_y.append(y[i])
    return np.array(out_t), np.array(out_y)


fig = plt.figure(figsize=(7.2, 5.7))
gs = fig.add_gridspec(2, 2, hspace=0.50, wspace=0.34)

PANELS = [('rmsd_bundle', 'C$\\alpha$ RMSD, bundle (Å)', 'a'),
          ('rg',          'Radius of gyration (Å)',      'b'),
          ('helix',       'Helix content (%)',                'c')]
for k, (key, lab, tag) in enumerate(PANELS):
    ax = fig.add_subplot(gs[k // 2, k % 2])
    for name, col, ls, leg in SERIES:
        t, y = broken(D[name]['time'], D[name][key])
        ax.plot(t, y, color=col, ls=ls, lw=1.3, label=leg, zorder=3)
    ax.set_xlabel('Time (ns)', fontsize=7.5)
    ax.set_ylabel(lab, fontsize=7.5)
    ax.set_xlim(0, 100)
    if k == 0:
        ax.legend(fontsize=6.4, loc='lower right', handlelength=1.8)
    F.panel(ax, tag)

# ---------------- (d) RMSF ----------------
ax = fig.add_subplot(gs[1, 1])
ax.axvspan(F.ARM[0], F.ARM[1], color=F.BAND, lw=0, zorder=0)
for name, col, ls, leg in SERIES:
    y = np.array(D[name]['rmsf'])
    ax.plot(np.arange(FIRST, FIRST + len(y)), y, color=col, ls=ls, lw=1.2, zorder=3)
ax.set_xlim(FIRST, FIRST + len(D['nat_water']['rmsf']) - 1)
ax.set_xlabel('Residue', fontsize=7.5)
ax.set_ylabel('C$\\alpha$ RMSF (Å)', fontsize=7.5)
ax.set_ylim(0, 3.2)
ax.text(17.5, 3.0, 'arm', fontsize=6.3, color=F.MUTED, ha='center')
_cmax = max(D[n]['rmsf'][-1] for n, _, _, _ in SERIES)
ax.annotate(f'native C-terminus\n{_cmax:.1f} \u00c5, off scale', xy=(165, 3.05),
            xytext=(0.50, 0.62), textcoords='axes fraction',
            fontsize=6.2, color=F.MUTED, ha='left', va='center',
            bbox=dict(boxstyle='round,pad=0.25', facecolor='white', edgecolor='none', alpha=0.92),
            arrowprops=dict(arrowstyle='-|>', color=F.MUTED, lw=0.7, shrinkA=3, shrinkB=2))
F.panel(ax, 'd')

paths = F.save(fig, 'Figure4')
for _p in paths:
    print('Figure4 ->', os.path.normpath(_p))
