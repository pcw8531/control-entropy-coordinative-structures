"""s3. Supplementary Table 6. Each measured structure against the four generated
topologies run at its own size, sixty realisations per cell, compression 1 - J.

Two columns of results are printed, and the difference between them is on the
record rather than buried.

    matched     the generators of SI Appendix section 2, which hold the link
                budget wherever the size allows it
    bare        the superseded scale-free form, grown from a bare seed at three
                links per node at every size

The values printed in the current SUPPLEMENTARY TABLE 6 come from the bare form. The three
non scale-free columns are the same either way, which is why only that one column
moves. Where a Barabasi-Albert graph cannot reach the row's mean degree the link
count actually used is printed, because at five, six and eleven elements no
integer number of links per added node reaches it.
"""
import sys
from pathlib import Path

import numpy as np
import networkx as nx

sys.path.insert(0, str(Path(__file__).resolve().parent))
from core import share, shannon_entropy, evenness
from generators import make, make_bare_scalefree, KINDS, SEEDS_60

ROWS = [
    ('Basketball', 5, 4, 0.0390),
    ('Handball', 6, 4, 0.1168),
    ('Football', 11, 6, 0.0208),
    ('Cortical activity', 100, 6, 0.0146),
    ('Musculoskeletal', 270, 6, 0.0421),
]


def giant(G):
    if nx.is_connected(G):
        return G
    return G.subgraph(max(nx.connected_components(G), key=len)).copy()


def compression_of(G):
    big = giant(G)
    c = nx.eigenvector_centrality_numpy(big)
    p = share(np.round(np.abs(np.array([c[n] for n in big.nodes()], float)), 12))
    return 1 - evenness(shannon_entropy(p), len(p))


def cell(kind, N, kmean, seeds=SEEDS_60, bare=False):
    vals, links = [], []
    for s in seeds:
        G = make_bare_scalefree(N, 3, s) if bare else make(kind, N, kmean, s)
        if giant(G).number_of_nodes() < 3:
            continue
        vals.append(compression_of(G))
        links.append(G.number_of_edges())
    return float(np.mean(vals)), float(np.std(vals, ddof=1)), float(np.mean(links))


def main():
    print('=' * 104)
    print('SI TABLE S8   sixty realisations per cell, seeds 42 to 101')
    print('=' * 104)
    for name, N, kmean, measured in ROWS:
        budget = N * kmean // 2
        print(f'\n{name}, N = {N}, stated mean degree {kmean}, '
              f'link budget {budget}')
        line = []
        for kind in KINDS:
            m, sd, lk = cell(kind, N, kmean)
            flag = '' if abs(lk - budget) < 1e-9 else f'  [{lk:.0f} links, not {budget}]'
            line.append((kind, m, sd, flag))
            print(f'    {kind:<14} {m:>8.4f} +/- {sd:.4f}{flag}')
        bm, bsd, blk = cell('scale-free', N, kmean, bare=True)
        print(f'    {"scale-free":<14} {bm:>8.4f} +/- {bsd:.4f}'
              f'  [{blk:.0f} links]   superseded bare-seed form, '
              f'the value currently printed in Table S8')
        sf = [v for v in line if v[0] == 'scale-free'][0][1]
        print(f'    measured {measured:.4f}   '
              f'{measured/sf:.2f}x the matched scale-free ensemble, '
              f'{measured/bm:.2f}x the bare-seed one')


if __name__ == '__main__':
    main()
