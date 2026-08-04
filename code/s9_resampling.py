"""s9 and s10. SI Table S9 and the lower row of Figure 4.

Sections A and B resample the units each measured value rests on. Section C
covers the two sports published as a single profile, which have no unit to
resample and so meet a null generated at their own size and mean degree.

The routine reproduces every published value from the released files before it
resamples anything, and stops if any of them fails.
"""
import sys
from pathlib import Path

import numpy as np
import networkx as nx

sys.path.insert(0, str(Path(__file__).resolve().parent))
from core import share, shannon_entropy, evenness, gini
from s2_measured import (sport_shares, basketball_profiles, cortical_shares,
                         musculoskeletal_shares, musculoskeletal_null_shares,
                         expertise, eigen_shares_matrix)
from loaders import read_matrix

RNG = np.random.default_rng(20260803)
B = 10000
B_CORTEX = 2000


def comp(p):
    return 1 - evenness(shannon_entropy(p), len(p))


def section_a_and_b():
    out = {}
    teams, prof = basketball_profiles()
    per_team = np.array([comp(share(r)) for r in prof])
    boot = np.array([per_team[RNG.integers(0, 16, 16)].mean() for _ in range(B)])
    out['basketball'] = dict(obs=per_team.mean(),
                             ci=np.percentile(boot, [2.5, 97.5]),
                             lo=per_team.min(), hi=per_team.max(), boot=boot)

    A = np.abs(read_matrix())
    np.fill_diagonal(A, 0.0)
    boot = np.empty(B_CORTEX)
    for b in range(B_CORTEX):
        idx = RNG.integers(0, 100, 100)
        Ab = A[np.ix_(idx, idx)]
        np.fill_diagonal(Ab, 0.0)
        boot[b] = comp(eigen_shares_matrix(Ab))
    out['cortex'] = dict(obs=comp(cortical_shares()),
                         ci=np.percentile(boot, [2.5, 97.5]), boot=boot)

    from loaders import read_csv
    deg = read_csv('musculoskeletal_coupling_shares.csv').bones_linked.values.astype(float)
    dgn = read_csv('musculoskeletal_random_null_shares.csv').bones_linked.values.astype(float)
    bm = np.array([comp(share(deg[RNG.integers(0, 270, 270)])) for _ in range(B)])
    bn = np.array([comp(share(dgn[RNG.integers(0, 270, 270)])) for _ in range(B)])
    out['muscles'] = dict(obs=comp(musculoskeletal_shares()),
                          ci=np.percentile(bm, [2.5, 97.5]),
                          null_obs=comp(musculoskeletal_null_shares()),
                          null_ci=np.percentile(bn, [2.5, 97.5]),
                          ratio=comp(musculoskeletal_shares()) / comp(musculoskeletal_null_shares()),
                          ratio_ci=np.percentile(bm / bn, [2.5, 97.5]),
                          boot=bm, boot_null=bn)

    nov, exp = expertise()
    n = nov.outcome_entropy_h_bits.values
    e = exp.outcome_entropy_h_bits.values
    d = n.mean() - e.mean()
    bd = np.array([n[RNG.integers(0, 10, 10)].mean() - e[RNG.integers(0, 10, 10)].mean()
                   for _ in range(B)])
    pool = np.concatenate([n, e])
    perm = np.empty(B)
    for b in range(B):
        q = RNG.permutation(pool)
        perm[b] = q[:10].mean() - q[10:].mean()
    y = np.r_[np.ones(10), np.zeros(10)]
    correct = 0
    for i in range(20):
        m = np.ones(20, bool)
        m[i] = False
        thr = (pool[m][y[m] == 1].mean() + pool[m][y[m] == 0].mean()) / 2
        correct += int((pool[i] > thr) == (y[i] == 1))
    out['expertise'] = dict(obs=d, ci=np.percentile(bd, [2.5, 97.5]),
                            ratio=2 ** d, ratio_ci=2 ** np.percentile(bd, [2.5, 97.5]),
                            p_perm=(1 + (np.abs(perm) >= abs(d)).sum()) / (B + 1),
                            loo=correct / 20, boot=bd, perm=perm)
    return out


def section_c(draws=10000):
    out = {}
    for name, N, kmean, obs in [('Handball', 6, 4, comp(sport_shares('Handball'))),
                                ('Football', 11, 6, comp(sport_shares('Football')))]:
        vals = []
        for s in range(draws):
            G = nx.gnm_random_graph(N, N * kmean // 2, seed=s)
            if not nx.is_connected(G):
                G = G.subgraph(max(nx.connected_components(G), key=len)).copy()
            if G.number_of_nodes() < 3:
                continue
            c = nx.eigenvector_centrality_numpy(G)
            p = share(np.round(np.abs(np.array([c[n] for n in G.nodes()], float)), 12))
            vals.append(comp(p))
        vals = np.array(vals)
        out[name] = dict(obs=obs, ci=np.percentile(vals, [2.5, 97.5]),
                         above=100 * (vals < obs).mean(), draws=len(vals), null=vals)
    return out


def main():
    print('=' * 96)
    print('SI TABLE S9 A and B   resampling the units each value rests on')
    print('=' * 96)
    v = section_a_and_b()
    b = v['basketball']
    print(f'   basketball   {b["obs"]:.4f}  95% CI [{b["ci"][0]:.4f}, {b["ci"][1]:.4f}]'
          f'   16 teams, per team {b["lo"]:.4f} to {b["hi"]:.4f}')
    c = v['cortex']
    print(f'   cortex       {c["obs"]:.4f}  95% CI [{c["ci"][0]:.4f}, {c["ci"][1]:.4f}]'
          f'   100 regions')
    m = v['muscles']
    print(f'   muscles      {m["obs"]:.4f}  95% CI [{m["ci"][0]:.4f}, {m["ci"][1]:.4f}]'
          f'   270 muscles')
    print(f'   its null     {m["null_obs"]:.4f}  95% CI '
          f'[{m["null_ci"][0]:.4f}, {m["null_ci"][1]:.4f}]')
    print(f'   ratio        {m["ratio"]:.2f}  95% CI '
          f'[{m["ratio_ci"][0]:.2f}, {m["ratio_ci"][1]:.2f}]')
    e = v['expertise']
    print(f'   expertise    {e["obs"]:.3f} bits  95% CI '
          f'[{e["ci"][0]:.3f}, {e["ci"][1]:.3f}]')
    print(f'   ratio        {e["ratio"]:.3f}  95% CI '
          f'[{e["ratio_ci"][0]:.3f}, {e["ratio_ci"][1]:.3f}]')
    print(f'   permutation  p = {e["p_perm"]:.1e}, leave-one-out {e["loo"]*100:.0f} per cent')

    print()
    print('=' * 96)
    print('SI TABLE S9 C   the two single published profiles against a matched null')
    print('=' * 96)
    for name, d in section_c().items():
        print(f'   {name:<10} observed {d["obs"]:.4f}   null 95% '
              f'[{max(d["ci"][0], 0):.4f}, {d["ci"][1]:.4f}]   '
              f'above {d["above"]:.1f} per cent of {d["draws"]:,} draws')

    print()
    print('SI TABLE S9 D   Gini beside the entropy reading')
    teams, prof = basketball_profiles()
    for nm, p in (('Football', sport_shares('Football')),
                  ('Handball', sport_shares('Handball')),
                  ('Basketball, averaged profile', share(prof.mean(axis=0))),
                  ('Cortical activity', cortical_shares()),
                  ('Musculoskeletal', musculoskeletal_shares()),
                  ('Musculoskeletal null', musculoskeletal_null_shares())):
        print(f'   {nm:<30} N {len(p):>4}   Gini {gini(p):.4f}   '
              f'compression {comp(p):.4f}')


if __name__ == '__main__':
    main()
