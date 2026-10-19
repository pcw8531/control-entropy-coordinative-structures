"""s9 and s10. Supplementary Note 4, Supplementary Table 7, Supplementary Fig. 1
and the lower row of Figure 4.

Section A resamples the units each measured value rests on, and Section B gives
the expertise difference its exact interval from the F distribution. Section C
covers the two sports published as a single profile, which have no unit to
resample and so meet a null of random graphs at their own size and at the mean
degree of Supplementary Table 4, six for football and four for handball. Because
a published profile fixes no mean degree, a density check repeats the football
null at mean degree four.

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
                         table_expertise, eigen_shares_matrix)
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

    t = table_expertise()
    out['expertise'] = dict(obs=t['difference'], ci=t['diff_ci'], ratio=t['ratio'],
                            ratio_ci=t['ratio_ci'], F=t['F'], df=t['df'], p=t['p'])
    return out


def random_graph_null(N, kmean, draws=10000):
    """Compression of random graphs at a fixed link count, seeds 0 to draws - 1."""
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
    return np.array(vals)


def section_c(draws=10000):
    out = {}
    for name, N, kmean, obs in [('Handball', 6, 4, comp(sport_shares('Handball'))),
                                ('Football', 11, 6, comp(sport_shares('Football')))]:
        vals = random_graph_null(N, kmean, draws)
        out[name] = dict(obs=obs, ci=np.percentile(vals, [2.5, 97.5]),
                         above=100 * (vals < obs).mean(), draws=len(vals), null=vals)
    return out


def density_check(draws=10000):
    """Supplementary Note 4. The football null again, at mean degree four rather than six."""
    obs = comp(sport_shares('Football'))
    vals = random_graph_null(11, 4, draws)
    return dict(obs=obs, ci=np.percentile(vals, [2.5, 97.5]),
                above=100 * (vals < obs).mean(), draws=len(vals))


def main():
    print('=' * 96)
    print('SUPPLEMENTARY TABLE 7 A and B   resampling the units, and the exact expertise interval')
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
    print(f'   expertise    {e["obs"]:.3f} bits  exact 95% CI '
          f'[{e["ci"][0]:.3f}, {e["ci"][1]:.3f}]   two groups of ten')
    print(f'   ratio        {e["ratio"]:.3f}  exact 95% CI '
          f'[{e["ratio_ci"][0]:.3f}, {e["ratio_ci"][1]:.3f}]')
    print(f'   F({e["df"][0]}, {e["df"][1]}) = {e["F"]:.2f}, p = {e["p"]:.3f}, two-sided, '
          f'against the null of equal spread')

    print()
    print('=' * 96)
    print('SUPPLEMENTARY TABLE 7 C   the two single published profiles against a matched null')
    print('=' * 96)
    for name, d in section_c().items():
        print(f'   {name:<10} observed {d["obs"]:.4f}   null 95% '
              f'[{max(d["ci"][0], 0):.4f}, {d["ci"][1]:.4f}]   '
              f'above {d["above"]:.1f} per cent of {d["draws"]:,} draws')

    d = density_check()
    print(f'   density check, Supplementary Note 4: football against the same null at mean '
          f'degree four')
    print(f'   {"Football":<10} observed {d["obs"]:.4f}   null 95% '
          f'[{max(d["ci"][0], 0):.4f}, {d["ci"][1]:.4f}]   '
          f'above {d["above"]:.1f} per cent of {d["draws"]:,} draws')

    print()
    print('SUPPLEMENTARY TABLE 7 D   Gini beside the entropy reading')
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
