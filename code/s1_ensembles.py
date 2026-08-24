"""s1. The simulated ensembles behind Table 1, Figure 3 and Supplementary Tables 1 to 5.

Run directly to print every table. Results are written to results/ as CSV.
"""
import sys
from pathlib import Path

import numpy as np
import networkx as nx

sys.path.insert(0, str(Path(__file__).resolve().parent))
from core import (share, shannon_entropy, max_entropy, evenness, effective_dof,
                  compression, degree_heterogeneity, degree_assortativity, gini)
from generators import make, KINDS, LABEL, SEEDS

RESULTS = Path(__file__).resolve().parent.parent / 'results'


def giant(G):
    """The largest connected component.

    The eigenvector calculation uses it because the leading eigenvector is not
    unique on a disconnected graph, while the degree statistics use the full
    generated graph, each being a property of the whole degree sequence.
    """
    if nx.is_connected(G):
        return G
    return G.subgraph(max(nx.connected_components(G), key=len)).copy()


def eigen_shares(G):
    c = nx.eigenvector_centrality_numpy(G)
    return share(np.round(np.abs(np.array([c[n] for n in G.nodes()], float)), 12))


def degree_shares(G):
    return share(np.array([d for _, d in G.degree()], float))


def ensemble(kind, N=100, kmean=6, seeds=SEEDS):
    """One ensemble, every statistic the tables report."""
    out = dict(H=[], Hdeg=[], kappa=[], r=[], gini=[], k=[], links=[])
    for s in seeds:
        full = make(kind, N, kmean, s)
        big = giant(full)
        p = eigen_shares(big)
        out['H'].append(shannon_entropy(p))
        out['Hdeg'].append(shannon_entropy(degree_shares(big)))
        out['gini'].append(gini(p))
        out['kappa'].append(degree_heterogeneity(full))
        out['r'].append(degree_assortativity(full))
        out['k'].append(2 * full.number_of_edges() / full.number_of_nodes())
        out['links'].append(full.number_of_edges())
    return {k: np.array(v, float) for k, v in out.items()}


def table_1(N=100, kmean=6, seeds=SEEDS):
    ceil = max_entropy(N)
    rows = []
    for kind in KINDS:
        e = ensemble(kind, N, kmean, seeds)
        h = e['H'].mean()
        finite = e['r'][np.isfinite(e['r'])]
        rows.append(dict(topology=LABEL[kind], H=h, J=evenness(h, N),
                         D_eff=effective_dof(h), dH=compression(h, N),
                         kappa=e['kappa'].mean(), H_degree=e['Hdeg'].mean(),
                         gini=e['gini'].mean(), mean_degree=e['k'].mean(),
                         links=e['links'].mean(),
                         r=finite.mean() if finite.size else float('nan'),
                         H_sd=e['H'].std(ddof=1)))
    return rows


def _print(rows, ceil):
    head = (f'{"Topology":<20}{"H":>9}{"J":>9}{"2^H":>8}{"dH":>8}{"kappa":>8}'
            f'{"H deg":>9}{"Gini":>8}{"<k>":>7}{"links":>7}{"r":>9}')
    print(head)
    for r in rows:
        rr = 'undefined' if not np.isfinite(r['r']) else f'{r["r"]:.2f}'
        print(f'{r["topology"]:<20}{r["H"]:>9.4f}{r["J"]:>9.4f}{r["D_eff"]:>8.1f}'
              f'{r["dH"]:>8.3f}{r["kappa"]:>8.2f}{r["H_degree"]:>9.4f}'
              f'{r["gini"]:>8.3f}{r["mean_degree"]:>7.2f}{r["links"]:>7.1f}{rr:>9}')
    print(f'ceiling log2 N = {ceil:.4f} bits')


def main():
    RESULTS.mkdir(exist_ok=True)
    print('=' * 96)
    print('TABLE 1 and SUPPLEMENTARY TABLE 1   N = 100, mean degree 6, seeds 42 to 141')
    print('=' * 96)
    rows = table_1()
    _print(rows, max_entropy(100))
    import csv
    with open(RESULTS / 'table_1.csv', 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)

    print()
    print('SUPPLEMENTARY TABLE 2   control entropy across network size, mean degree 6')
    print(f'{"N":>6}' + ''.join(f'{LABEL[k]:>20}' for k in KINDS))
    for N in (50, 100, 200, 500):
        vals = [ensemble(k, N, 6)['H'].mean() for k in KINDS]
        print(f'{N:>6}' + ''.join(f'{v:>20.4f}' for v in vals))

    print()
    print('SUPPLEMENTARY TABLE 3   control entropy across mean degree, N = 100')
    print(f'{"<k>":>6}' + ''.join(f'{LABEL[k]:>20}' for k in KINDS))
    for km in (4, 6, 10):
        vals = [ensemble(k, 100, km)['H'].mean() for k in KINDS]
        print(f'{km:>6}' + ''.join(f'{v:>20.4f}' for v in vals))

    print()
    print('SUPPLEMENTARY TABLE 4   the two coupling weightings, N = 100')
    print(f'{"weighting":>12}' + ''.join(f'{LABEL[k]:>20}' for k in KINDS))
    ens = {k: ensemble(k) for k in KINDS}
    print(f'{"degree":>12}' + ''.join(f'{ens[k]["Hdeg"].mean():>20.4f}' for k in KINDS))
    print(f'{"eigenvector":>12}' + ''.join(f'{ens[k]["H"].mean():>20.4f}' for k in KINDS))

    print()
    print('SUPPLEMENTARY NOTE 3   bootstrap intervals and the separation')
    rng = np.random.default_rng(1)
    for k in KINDS:
        v = ens[k]['H']
        bs = np.array([v[rng.integers(0, len(v), len(v))].mean() for _ in range(10000)])
        lo, hi = np.percentile(bs, [2.5, 97.5])
        print(f'   {LABEL[k]:<20} mean {v.mean():.4f}  95% CI [{lo:.4f}, {hi:.4f}]')
    a, b = ens['random']['H'], ens['scale-free']['H']
    sp = np.sqrt((a.var(ddof=1) + b.var(ddof=1)) / 2)
    print(f'   random minus scale-free {a.mean()-b.mean():.4f} bits, '
          f'Cohen d {(a.mean()-b.mean())/sp:.2f}')

    print()
    print('SUPPLEMENTARY TABLE 5   five-fold accuracy from a single feature, four classes')
    try:
        table_accuracy(ens)
    except ImportError:
        print('   scikit-learn is not installed, skipped')

    print()
    print('RESULTS TEXT   the same four generators at 13 elements and 26 links')
    vals = [ensemble(k, 13, 4)['H'].mean() / max_entropy(13) for k in KINDS]
    print('   evenness ' + '  '.join(f'{v:.4f}' for v in vals))


def table_accuracy(ens=None):
    from sklearn.svm import SVC
    from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
    from sklearn.linear_model import LogisticRegression
    from sklearn.model_selection import cross_val_score, StratifiedKFold
    from sklearn.pipeline import make_pipeline
    from sklearn.preprocessing import StandardScaler
    ens = ens or {k: ensemble(k) for k in KINDS}
    X_h, X_k, y = [], [], []
    for ci, k in enumerate(KINDS):
        X_h += list(ens[k]['H'])
        X_k += list(ens[k]['kappa'])
        y += [ci] * len(ens[k]['H'])
    X_h = np.array(X_h).reshape(-1, 1)
    X_k = np.array(X_k).reshape(-1, 1)
    y = np.array(y)
    cv = StratifiedKFold(5, shuffle=True, random_state=0)
    models = [('SVM (rbf)', make_pipeline(StandardScaler(), SVC())),
              ('Random forest', RandomForestClassifier(random_state=0)),
              ('Gradient boosting', GradientBoostingClassifier(random_state=0)),
              ('Logistic regression', make_pipeline(
                  StandardScaler(), LogisticRegression(max_iter=2000)))]
    print(f'   {"model":<22}{"entropy alone":>16}{"kappa alone":>14}')
    for nm, mdl in models:
        e = 100 * cross_val_score(mdl, X_h, y, cv=cv).mean()
        kk = 100 * cross_val_score(mdl, X_k, y, cv=cv).mean()
        print(f'   {nm:<22}{e:>16.1f}{kk:>14.1f}')


if __name__ == '__main__':
    main()
