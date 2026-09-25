#!/usr/bin/env python3
"""
Sequence, hydropathy and AlphaFold3 analysis for the ApoE reverse-QTY manuscript.

Writes the JSON that mkfig1.py, mkfig2.py and mkfigS.py plot:
    _analysis/hydro.json   hydropathy and amphipathicity, whole-segment and per-helix
    _analysis/af3.json     pLDDT, pairwise RMSD matrices, per-residue deviation, SASA
    _analysis/rsa.json     per-residue relative solvent accessibility, both ensembles

Inputs are the five AlphaFold3 models of each sequence and the two sequence files.
Every derived number is printed against its published value so the run self-checks.
"""
import os, json, itertools
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
OUT  = ANALYSIS; os.makedirs(OUT, exist_ok=True)
ROOT = os.path.normpath(os.path.join(HERE, '..', '..'))          # the APOE-rQTY folder
AF3  = AF3DIR
SEQD = os.environ.get('SEQ_DIR') or os.path.join(ROOT, 'data')

FIRST, LAST = 11, 167
BUNDLE = (25, 167)
ARM    = (11, 24)
HEL = [('H0', 11, 19), ('H1', 25, 41), ('H1b', 45, 52),
       ('H2', 55, 81), ('H3', 87, 124), ('H4', 131, 167)]
SUBS = [16, 17, 18, 36, 41, 46, 48, 55, 57, 58, 67, 74, 81, 89, 98, 101, 117, 118, 123, 156, 162, 163]
NOSUB = [21, 24, 42, 83, 128, 130]

KD = dict(A=1.8, R=-4.5, N=-3.5, D=-3.5, C=2.5, Q=-3.5, E=-3.5, G=-0.4, H=-3.2, I=4.5,
          L=3.8, K=-3.9, M=1.9, F=2.8, P=-1.6, S=-0.8, T=-0.7, W=-0.9, Y=-1.3, V=4.2)
EIS = dict(A=0.62, R=-2.53, N=-0.78, D=-0.90, C=0.29, Q=-0.85, E=-0.74, G=0.48, H=-0.40,
           I=1.38, L=1.06, K=-1.50, M=0.64, F=1.19, P=0.12, S=-0.18, T=-0.05, W=0.81,
           Y=0.26, V=1.08)
MAXASA = dict(A=129, R=274, N=195, D=193, C=167, Q=225, E=223, G=104, H=224, I=197, L=201,
              K=236, M=224, F=240, P=159, S=155, T=172, W=285, Y=263, V=174)
APOLAR = set('AVLIFMWPC')


def read_seq(path):
    s = ''.join(l.strip() for l in open(path) if not l.startswith('>'))
    return ''.join(c for c in s if c.isalpha())


def gravy(seq):
    return sum(KD[c] for c in seq) / len(seq)


def moment(seq, delta=100.0):
    """Eisenberg mean hydrophobicity and mean hydrophobic moment."""
    h = np.array([EIS[c] for c in seq])
    ang = np.deg2rad(delta) * np.arange(len(seq))
    mu = np.hypot((h * np.sin(ang)).sum(), (h * np.cos(ang)).sum())
    return h.mean(), mu / len(seq)


def ca(t, lo=None, hi=None):
    return np.array([a.index for a in t.topology.atoms if a.name == 'CA'
                     and (lo is None or lo <= a.residue.index + FIRST <= hi)])


def superpose(P, Q):
    """Kabsch: rotate P onto Q. Returns rmsd and the rotated, centred P."""
    cp, cq = P.mean(0), Q.mean(0)
    U, S, Vt = np.linalg.svd((P - cp).T @ (Q - cq))
    d = np.sign(np.linalg.det(U @ Vt))
    R = U @ np.diag([1, 1, d]) @ Vt
    Pr = (P - cp) @ R
    return float(np.sqrt(((Pr - (Q - cq)) ** 2).sum(1).mean())), Pr, Q - cq, R, cp


def main():
    nat = read_seq(os.path.join(SEQD, 'APOE_N'))
    rq  = read_seq(os.path.join(SEQD, 'APOE_N_rqty'))
    assert len(nat) == len(rq) == 157
    diffs = [i + FIRST for i, (a, b) in enumerate(zip(nat, rq)) if a != b]
    assert diffs == SUBS, diffs

    # ---------------- hydropathy ----------------
    H = dict(nat_gravy=gravy(nat), rq_gravy=gravy(rq))
    w = 9
    H['win_nat'] = [gravy(nat[i:i + w]) for i in range(len(nat) - w + 1)]
    H['win_rq']  = [gravy(rq[i:i + w])  for i in range(len(rq) - w + 1)]
    H['win_x']   = [FIRST + i + w // 2 for i in range(len(nat) - w + 1)]
    H['seg'] = []
    for nm, a, b in HEL:
        sn, sr = nat[a - FIRST:b - FIRST + 1], rq[a - FIRST:b - FIRST + 1]
        hn, mn = moment(sn); hr, mr = moment(sr)
        H['seg'].append(dict(name=nm, a=a, b=b, n=len(sn),
                             subs=sum(1 for p in SUBS if a <= p <= b),
                             gravy_nat=gravy(sn), gravy_rq=gravy(sr),
                             H_nat=hn, H_rq=hr, uH_nat=mn, uH_rq=mr))
    json.dump(H, open(os.path.join(OUT, 'hydro.json'), 'w'))

    nnat = sum(1 for v in H['win_nat'] if v > 0)
    nrq  = sum(1 for v in H['win_rq'] if v > 0)
    print('HYDROPATHY')
    print(f'  GRAVY  native {H["nat_gravy"]:+.4f} (published -0.757)   '
          f'rQTY {H["rq_gravy"]:+.4f} (published +0.123)')
    print(f'  net-hydrophobic windows  {nnat}/{len(H["win_nat"])} native (published 17/149), '
          f'{nrq}/{len(H["win_rq"])} rQTY (published 77/149)')
    print('  per-segment (published Table 2 in brackets)')
    pub = {'H0': (-2.278, -0.111, -0.591, -0.041, 0.306, 0.445),
           'H1': (-0.482, 0.188, -0.049, 0.118, 0.392, 0.462),
           'H1b': (-0.713, 1.113, -0.090, 0.388, 0.390, 0.331),
           'H2': (-1.033, 0.293, -0.263, 0.067, 0.366, 0.404),
           'H3': (-0.682, 0.324, -0.249, 0.007, 0.183, 0.225),
           'H4': (-0.814, -0.308, -0.452, -0.324, 0.330, 0.392)}
    for s in H['seg']:
        p = pub[s['name']]
        print(f'    {s["name"]:4s} GRAVY {s["gravy_nat"]:+.3f}/{s["gravy_rq"]:+.3f} '
              f'[{p[0]:+.3f}/{p[1]:+.3f}]  <H> {s["H_nat"]:+.3f}/{s["H_rq"]:+.3f} '
              f'[{p[2]:+.3f}/{p[3]:+.3f}]  <uH> {s["uH_nat"]:.3f}/{s["uH_rq"]:.3f} '
              f'[{p[4]:.3f}/{p[5]:.3f}]')

    # ---------------- AlphaFold3 ----------------
    N = [md.load(f'{AF3}/fold_apoe_nat/fold_apoe_nat_model_{i}.cif') for i in range(5)]
    R = [md.load(f'{AF3}/fold_apoe_rqty/fold_apoe_rqty_model_{i}.cif') for i in range(5)]
    A = {}
    for tag, reg in [('full', (None, None)), ('bundle', BUNDLE)]:
        def X(t):
            return t.xyz[0][ca(t, *reg)] * 10
        A[tag] = dict(
            nat=[superpose(X(N[i]), X(N[j]))[0] for i, j in itertools.combinations(range(5), 2)],
            rqty=[superpose(X(R[i]), X(R[j]))[0] for i, j in itertools.combinations(range(5), 2)],
            cross=[superpose(X(R[i]), X(N[j]))[0] for i in range(5) for j in range(5)])
        # full 10x10 matrix for Supplementary Figure S1
        mods = N + R
        A[tag + '_mat'] = [[superpose(X(a), X(b))[0] for b in mods] for a in mods]

    # per-residue deviation of each rQTY model from the top native model, on the bundle
    ib = ca(N[0], *BUNDLE); ia = ca(N[0])
    Qb = N[0].xyz[0][ib] * 10
    dev = []
    for t in R:
        Pb = t.xyz[0][ca(t, *BUNDLE)] * 10
        _, _, _, Rot, cp = superpose(Pb, Qb)
        P = (t.xyz[0][ca(t)] * 10 - cp) @ Rot
        Q = N[0].xyz[0][ia] * 10 - Qb.mean(0)
        dev.append(np.linalg.norm(P - Q, axis=1).tolist())
    A['dev'] = dev

    def atom_plddt(path):
        return np.array([float(l.split()[14]) for l in open(path) if l.startswith('ATOM')])

    def res_plddt(t, path):
        v = atom_plddt(path); out = []
        for r in t.topology.residues:
            idx = [a.index for a in r.atoms]
            out.append(float(v[idx].mean()))
        return out

    A['plddt_nat'] = res_plddt(N[0], f'{AF3}/fold_apoe_nat/fold_apoe_nat_model_0.cif')
    A['plddt_rqty'] = res_plddt(R[0], f'{AF3}/fold_apoe_rqty/fold_apoe_rqty_model_0.cif')

    def sasa(models, seq):
        tot, apo = [], []
        for t in models:
            s = md.shrake_rupley(t, probe_radius=0.14, n_sphere_points=960, mode='atom')[0] * 100
            tot.append(float(s.sum()))
            ai = [a.index for a in t.topology.atoms if a.residue.code in APOLAR]
            apo.append(float(s[ai].sum()))
        return tot, apo

    A['sasa_nat'] = dict(zip(('total', 'apolar'), sasa(N, nat)))
    A['sasa_rqty'] = dict(zip(('total', 'apolar'), sasa(R, rq)))
    json.dump(A, open(os.path.join(OUT, 'af3.json'), 'w'))

    # ---------------- per-residue RSA ----------------
    def rsa(models, seq):
        acc = []
        for t in models:
            s = md.shrake_rupley(t, probe_radius=0.14, n_sphere_points=960, mode='residue')[0] * 100
            acc.append([s[i] / MAXASA[seq[i]] * 100 for i in range(len(seq))])
        return np.array(acc).mean(0).tolist()

    S = dict(rsa_nat=rsa(N, nat), rsa_rqty=rsa(R, rq), subs=SUBS, nosub=NOSUB,
             nat=nat, rq=rq)
    json.dump(S, open(os.path.join(OUT, 'rsa.json'), 'w'))

    print('\nALPHAFOLD3')
    for tag, p in [('full', '0.52 / 4.18 / 4.34'), ('bundle', '0.21 / 0.37 / 0.81')]:
        a = A[tag]
        print(f'  {tag:7s} within-nat {np.mean(a["nat"]):.2f} ({min(a["nat"]):.2f}-{max(a["nat"]):.2f})  '
              f'within-rQTY {np.mean(a["rqty"]):.2f} ({min(a["rqty"]):.2f}-{max(a["rqty"]):.2f})  '
              f'cross {np.mean(a["cross"]):.2f} ({min(a["cross"]):.2f}-{max(a["cross"]):.2f})   published {p}')
    for tag, lab in [('plddt_nat', 'native'), ('plddt_rqty', 'rQTY')]:
        v = np.array(A[tag]); rr = np.arange(FIRST, FIRST + len(v))
        print(f'  pLDDT {lab:6s} arm {v[(rr>=ARM[0])&(rr<=ARM[1])].mean():.1f}  '
              f'bundle {v[rr>=BUNDLE[0]].mean():.1f}   published '
              f'{"57.2 / 90.7" if lab=="native" else "52.2 / 84.4"}')
    for tag, p in [('sasa_nat', '9776±82 / 2079±32'), ('sasa_rqty', '9646±236 / 3313±174')]:
        t_, a_ = np.array(A[tag]['total']), np.array(A[tag]['apolar'])
        print(f'  SASA {tag[5:]:5s} total {t_.mean():.0f}±{t_.std():.0f}  '
              f'apolar {a_.mean():.0f}±{a_.std():.0f}   published {p}')
    an = np.mean(A['sasa_nat']['apolar']); ar = np.mean(A['sasa_rqty']['apolar'])
    print(f'  apolar increase {100*(ar-an)/an:.1f}%   published 59%')
    rn = np.array(S['rsa_nat']); rr2 = np.array(S['rsa_rqty'])
    si = [p - FIRST for p in SUBS]
    print(f'  RSA at substituted positions {rn[si].mean():.1f}% / {rr2[si].mean():.1f}%'
          f'   published 27.2 / 27.9')


if __name__ == '__main__':
    main()
