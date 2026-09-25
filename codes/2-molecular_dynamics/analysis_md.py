#!/usr/bin/env python3
"""
Trajectory analysis for the ApoE reverse-QTY manuscript.

Reads the GROMACS production trajectories under data/MD/ and writes the
per-frame arrays that the figure scripts plot, as JSON.

Run one system at a time:   python3 analysis_md.py nat_water
                            python3 analysis_md.py rqty_bilayer
                            python3 analysis_md.py rqty_hmmm

IMPORTANT -- bilayer segment 2.
step7_2.xtc was re-run on 2025-09-10 and interrupted at 4.18 ns, overwriting the
original complete 10-20 ns segment.  The saved coordinates therefore cover
0-14 ns and 20-50 ns.  The 5 frames of that abandoned branch are excluded here;
the remaining 44 frames carry their true physical time, with a gap at 14-20 ns.
The water (10 x 10 ns) and mimetic (9 x 10 ns) runs are complete.
"""
import sys, os, json, warnings
import numpy as np
import MDAnalysis as mda
from MDAnalysis.lib.distances import distance_array
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
warnings.filterwarnings('ignore')

ROOT = MDDIR
OUT  = ANALYSIS
os.makedirs(OUT, exist_ok=True)

# segment file, and the physical start time (ns) of each segment
SYSTEMS = {
    'nat_water': dict(
        d='Apoenat_water', top='gromacs/step3_input.gro',
        segs=[(f'gromacs/step5_{i}.xtc', (i - 1) * 10.0) for i in range(1, 11)],
        membrane=False),
    'rqty_bilayer': dict(
        d='Apoerqtybilayer', top='gromacs/step5_input.gro',
        # segment 2 omitted: overwritten by an interrupted re-run (see module docstring)
        segs=[('gromacs/step7_1.xtc', 0.0),
              ('gromacs/step7_3.xtc', 20.0),
              ('gromacs/step7_4.xtc', 30.0),
              ('gromacs/step7_5.xtc', 40.0)],
        membrane=True),
    'rqty_hmmm': dict(
        d='Apoerqty_membrane_hmm', top='gromacs/step5_input.gro',
        segs=[(f'gromacs/step7_{i}.xtc', (i - 1) * 10.0) for i in range(1, 10)],
        membrane=True),
}

BUNDLE = (25, 167)          # bundle residues, mature ApoE numbering
FIRST  = 11                 # first residue of the construct
APOLAR = {'ALA','VAL','LEU','ILE','PHE','MET','TRP','PRO','CYS'}
BINS   = [(0,12),(12,18),(18,24),(24,30),(30,40),(40,60)]
LIPIDS = ['POPC','POPE','POPS','POPI','PSM','CHL1']


def unwrap(pos, box):
    """Make the protein whole by a sequential minimum-image walk along atom order."""
    out = pos.copy()
    for i in range(1, len(out)):
        d = out[i] - out[i - 1]
        out[i] -= box[:3] * np.round(d / box[:3])
    return out


def kabsch_rmsd(P, Q):
    Pc, Qc = P - P.mean(0), Q - Q.mean(0)
    U, S, Vt = np.linalg.svd(Pc.T @ Qc)
    d = np.sign(np.linalg.det(U @ Vt))
    R = U @ np.diag([1, 1, d]) @ Vt
    return np.sqrt(((Pc @ R - Qc) ** 2).sum(1).mean()), Pc @ R, Qc


def run(key):
    cfg = SYSTEMS[key]
    base = os.path.join(ROOT, cfg['d'])
    traj = [os.path.join(base, s) for s, _ in cfg['segs']]
    u = mda.Universe(os.path.join(base, cfg['top']), traj)

    prot = u.select_atoms('protein and not name H*')
    ca   = u.select_atoms('protein and name CA')
    resids = np.array([a.resid for a in ca]) + FIRST - 1
    bmask  = (resids >= BUNDLE[0]) & (resids <= BUNDLE[1])
    apol_idx = np.array([i for i, a in enumerate(prot) if a.resname in APOLAR])

    # physical time of every frame, from each segment's start offset
    times = []
    for (fn, t0), n in zip(cfg['segs'], seg_frame_counts(traj)):
        times.extend(t0 + np.linspace(0, 10.0, n))
    times = np.array(times)

    if cfg['membrane']:
        phos  = u.select_atoms('name P')
        basic = u.select_atoms('(resname ARG and name NE NH1 NH2) or (resname LYS and name NZ)')
        allPO = u.select_atoms('(resname POPC POPE POPS POPI PSM) and name O11 O12 O13 O14')
        aniPO = u.select_atoms('(resname POPS POPI) and name O11 O12 O13 O14')
        wat   = u.select_atoms('resname TIP3 and name OH2')
        lipall = u.select_atoms('resname ' + ' '.join(LIPIDS) + ' and not name H*')
        lipres = {rn: u.select_atoms(f'resname {rn} and not name H*') for rn in LIPIDS
                  if len(u.select_atoms(f'resname {rn}'))}
        bulk = {rn: lipres[rn].n_residues for rn in lipres}
        # heavy-atom slices per lipid residue, for first-shell counting
        lslices = {}
        for rn, g in lipres.items():
            j, sl = 0, []
            for r in g.residues:
                m = len(r.atoms.select_atoms('not name H*')); sl.append((j, j + m)); j += m
            lslices[rn] = sl
        pslices, j = [], 0
        for r in prot.residues:
            m = len(r.atoms.select_atoms('not name H*')); pslices.append((j, j + m)); j += m

    ref = None
    rec = dict(time=times.tolist(), rmsd_all=[], rmsd_bundle=[], rg=[], helix=[],
               sasa=[], hsasa=[], perhelix={}, frac=[], zcom=[], ncont=[], nresc=[],
               wcore=[], ec_all=[], ec_ani=[], ec_n=[], prof={f'{a}-{b}': [] for a, b in BINS},
               enr={rn: [] for rn in (lipres if cfg['membrane'] else {})})
    permin_acc, perz_acc, ss_acc, xyz_acc = [], [], [], []

    for ts in u.trajectory:
        box = ts.dimensions
        P = unwrap(prot.positions, box)
        C = P[[i for i, a in enumerate(prot) if a.name == 'CA']]
        if ref is None:
            ref = C.copy(); refb = C[bmask].copy()
        rec['rmsd_all'].append(kabsch_rmsd(C, ref)[0])
        rec['rmsd_bundle'].append(kabsch_rmsd(C[bmask], refb)[0])
        rec['rg'].append(float(np.sqrt(((C - C.mean(0)) ** 2).sum(1).mean())))
        xyz_acc.append(C)

        if cfg['membrane']:
            z = phos.positions[:, 2]; mid = z.mean()
            up, lo = z[z > mid].mean(), z[z <= mid].mean()
            pz = P[:, 2]
            rec['frac'].append(float(((pz > lo) & (pz < up)).mean()))
            pc = P.mean(0)
            rec['zcom'].append(float(pc[2] - mid))
            D = distance_array(P, lipall.positions, box=box)
            mins = np.array([D[a:b].min() for a, b in pslices])
            permin_acc.append(mins)
            perz_acc.append(np.abs(C[:, 2] - mid))
            rec['ncont'].append(int((D < 4.0).sum()))
            rec['nresc'].append(int((mins < 4.0).sum()))
            # annular phosphate-plane profile
            dxy = np.linalg.norm(phos.positions[:, :2] - pc[:2], axis=1)
            for a, b in BINS:
                m = (dxy >= a) & (dxy < b)
                zu, zl = z[m & (z > mid)], z[m & (z <= mid)]
                rec['prof'][f'{a}-{b}'].append(
                    float(zu.mean() - zl.mean()) if len(zu) >= 3 and len(zl) >= 3 else np.nan)
            # core hydration
            wp = wat.positions
            rec['wcore'].append(int((( np.abs(wp[:, 2] - mid) < 10) &
                (np.linalg.norm(wp[:, :2] - pc[:2], axis=1) < 25)).sum()))
            # basic-N contacts, nested atom sets
            rec['ec_all'].append(int((distance_array(basic.positions, allPO.positions, box=box) < 4).sum()))
            rec['ec_ani'].append(int((distance_array(basic.positions, aniPO.positions, box=box) < 4).sum()))
            rec['ec_n'].append(int((distance_array(basic.positions, allPO.positions, box=box) < 4).any(1).sum()))
            # first-shell lipid counts
            for rn, g in lipres.items():
                Dl = distance_array(P, g.positions, box=box)
                rec['enr'][rn].append(int(sum(1 for a, b in lslices[rn] if Dl[:, a:b].min() < 5.0)))

    # secondary structure and SASA in one mdtraj pass over the same frames
    t = md.load(traj, top=os.path.join(base, cfg['top']), atom_indices=prot.indices)
    ss = md.compute_dssp(t, simplified=False)
    hel = np.isin(ss, ['H', 'G', 'I'])
    rec['helix'] = (100 * hel.mean(1)).tolist()
    sasa = md.shrake_rupley(t, probe_radius=0.14, n_sphere_points=960, mode='atom') * 100  # A^2
    rec['sasa'] = sasa.sum(1).tolist()
    rec['hsasa'] = sasa[:, apol_idx].sum(1).tolist()
    HEL = {'H0': (11,19), 'H1': (25,41), 'H1b': (45,52), 'H2': (55,81), 'H3': (87,124), 'H4': (131,167)}
    h = slice(len(hel) // 2, None)
    rr = np.arange(FIRST, FIRST + ss.shape[1])
    inh = np.zeros(ss.shape[1], bool)
    for nm, (a, b) in HEL.items():
        m = (rr >= a) & (rr <= b); inh |= m
        rec['perhelix'][nm] = float(100 * hel[h][:, m].mean())
    rec['perhelix']['loops'] = float(100 * hel[h][:, ~inh].mean())

    # RMSF over the second half, superposed on the bundle
    X = np.array(xyz_acc); hh = slice(len(X) // 2, None)
    Xs = []
    for fr in X[hh]:
        _, a, b = kabsch_rmsd(fr[bmask], X[hh][0][bmask])
        c = fr - fr[bmask].mean(0)
        U, S, Vt = np.linalg.svd((fr[bmask] - fr[bmask].mean(0)).T @
                                 (X[hh][0][bmask] - X[hh][0][bmask].mean(0)))
        d = np.sign(np.linalg.det(U @ Vt)); R = U @ np.diag([1, 1, d]) @ Vt
        Xs.append(c @ R)
    Xs = np.array(Xs)
    rec['rmsf'] = np.sqrt(((Xs - Xs.mean(0)) ** 2).sum(-1).mean(0)).tolist()

    if cfg['membrane']:
        rec['permin'] = np.array(permin_acc).mean(0).tolist()
        rec['perz']   = np.array(perz_acc)[hh].mean(0).tolist()
        rec['bulk']   = bulk
        rec['nbasic'] = len(basic)

    json.dump(rec, open(os.path.join(OUT, key + '.json'), 'w'))
    print(f'{key}: {len(times)} frames, {times[0]:.0f}-{times[-1]:.0f} ns -> {key}.json')


def seg_frame_counts(paths):
    import struct
    counts = []
    for p in paths:
        data = open(p, 'rb').read(); off = 0; n = 0
        while off + 16 <= len(data):
            magic, natoms, step = struct.unpack_from('>iii', data, off)
            if magic != 1995: break
            n += 1; q = off + 52
            na = struct.unpack_from('>i', data, q)[0]; q += 4
            if na <= 9: q += na * 3 * 4
            else:
                q += 4 + 6 * 4 + 4
                nb = struct.unpack_from('>i', data, q)[0]; q += 4 + nb; q = (q + 3) // 4 * 4
            off = q
        counts.append(n)
    return counts


if __name__ == '__main__':
    run(sys.argv[1])
