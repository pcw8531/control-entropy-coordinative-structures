"""s3. Supplementary Table 4. Each measured structure against the four generated
topologies run at its own size, sixty realisations per cell, compression 1 - J.

The generators are those of Supplementary Note 1, which hold the link budget
wherever the size allows it. The scale-free ensemble grows from a ring lattice
seed of twelve elements at N = 100 and N = 270 and of seven elements at N = 11,
the sizes Supplementary Table 4 states. At five and six elements no scale-free
ensemble exists at the row's density, since ten links on five elements is the
complete graph and twelve links on six leave only three pairs unlinked, so those
two cells are empty.

The superseded scale-free form, grown from a bare seed at three links per node at
every size, is printed beside each row for the record. Supplementary Table 4 does
not use it.
"""
import sys
from pathlib import Path

import numpy as np
import networkx as nx

sys.path.insert(0, str(Path(__file__).resolve().parent))
from core import share, shannon_entropy, evenness
from generators import make, make_bare_scalefree, KINDS, SEEDS_60

# name, N, mean degree, measured compression, scale-free seed size (None where no
# scale-free ensemble exists at the row's density)
ROWS = [
    ('Basketball', 5, 4, 0.0390, None),
    ('Handball', 6, 4, 0.1168, None),
    ('Football', 11, 6, 0.0208, 7),
    ('Cortical activity', 100, 6, 0.0146, 12),
    ('Musculoskeletal', 270, 6, 0.0421, 12),
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


def cell(kind, N, kmean, seeds=SEEDS_60, bare=False, seed_size=None):
    vals, links = [], []
    for s in seeds:
        G = make_bare_scalefree(N, 3, s) if bare else make(kind, N, kmean, s, seed_size)
        if giant(G).number_of_nodes() < 3:
            continue
        vals.append(compression_of(G))
        links.append(G.number_of_edges())
    return float(np.mean(vals)), float(np.std(vals, ddof=1)), float(np.mean(links))


def main():
    print('=' * 104)
    print('SUPPLEMENTARY TABLE 4   sixty realisations per cell, seeds 42 to 101')
    print('=' * 104)
    for name, N, kmean, measured, sf_seed in ROWS:
        budget = N * kmean // 2
        print(f'\n{name}, N = {N}, stated mean degree {kmean}, '
              f'link budget {budget}')
        line = []
        for kind in KINDS:
            if kind == 'scale-free' and sf_seed is None:
                print(f'    {kind:<14} {"--":>8}   no scale-free ensemble at this density, '
                      f'cell left empty')
                continue
            m, sd, lk = cell(kind, N, kmean,
                             seed_size=sf_seed if kind == 'scale-free' else None)
            flag = '' if abs(lk - budget) < 1e-9 else f'  [{lk:.0f} links, not {budget}]'
            if kind == 'scale-free':
                flag += f'  [{sf_seed}-element ring seed]'
            line.append((kind, m, sd, flag))
            print(f'    {kind:<14} {m:>8.4f} +/- {sd:.4f}{flag}')
        bm, bsd, blk = cell('scale-free', N, kmean, bare=True)
        print(f'    {"scale-free":<14} {bm:>8.4f} +/- {bsd:.4f}'
              f'  [{blk:.0f} links]   superseded bare-seed form, for the record, '
              f'not used in Supplementary Table 4')
        sf = [v for v in line if v[0] == 'scale-free']
        if sf:
            # The ratios use the values as Supplementary Table 4 prints them, to
            # four decimals, so that they can be checked against the table.
            sf4, bm4 = round(sf[0][1], 4), round(bm, 4)
            print(f'    measured {measured:.4f}   '
                  f'{measured/sf4:.2f}x the matched scale-free ensemble, '
                  f'{measured/bm4:.2f}x the bare-seed one')
        else:
            print(f'    measured {measured:.4f}   no matched scale-free ensemble')


if __name__ == '__main__':
    main()
