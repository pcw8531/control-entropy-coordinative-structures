"""The four matched topology ensembles, exactly as Materials and Methods and
Supplementary Note 2 describes them.

Every realisation of every ensemble holds N * kmean / 2 links, so the mean degree
is identical in all four and the ensembles differ only in how those links are
placed. At N = 100 and mean degree 6 that is 300 links in every realisation.

    regular       ring lattice, rewiring probability 0
    small-world   the same generator at rewiring probability 0.1
    random        Erdos-Renyi at a fixed link count, not a fixed edge probability
    scale-free    Barabasi-Albert adding m = kmean / 2 links per node, grown from
                  a ring lattice seed of 4m elements at the same mean degree

The scale-free seed is what holds the link budget. Growth from a bare seed stops
at 291 links and a mean degree of 5.82 at N = 100, which is why the seed is
there. Where 4m is not smaller than N the seed cannot be built and the bare form
is used; `link_budget_met` reports whether the row reached its stated mean degree.
"""
import networkx as nx


def ring_seed(m, kmean, seed):
    """The 4m element ring lattice the scale-free ensemble grows from."""
    return nx.watts_strogatz_graph(4 * m, kmean, 0.0, seed=seed)


def make(kind, N, kmean, seed):
    """One realisation of one ensemble."""
    m = max(1, kmean // 2)
    if kind == 'regular':
        return nx.watts_strogatz_graph(N, kmean, 0.0, seed=seed)
    if kind == 'small-world':
        return nx.watts_strogatz_graph(N, kmean, 0.1, seed=seed)
    if kind == 'random':
        return nx.gnm_random_graph(N, N * kmean // 2, seed=seed)
    if kind == 'scale-free':
        if 4 * m < N:
            return nx.barabasi_albert_graph(
                N, m, seed=seed, initial_graph=ring_seed(m, kmean, seed))
        return nx.barabasi_albert_graph(N, m, seed=seed)
    raise ValueError(f'unknown topology: {kind}')


def make_bare_scalefree(N, m, seed):
    """The superseded form, kept so the pre-tightening values stay reproducible.

    Grown from a bare seed, so the link count falls short of the budget the other
    three ensembles hold. Used nowhere in the current tables; retained because
    earlier drafts reported values from it and a reader may want to recover them.
    """
    return nx.barabasi_albert_graph(N, m, seed=seed)


def link_budget_met(kind, N, kmean, seed):
    """True when this realisation actually carries N * kmean / 2 links."""
    return make(kind, N, kmean, seed).number_of_edges() == N * kmean // 2


KINDS = ['regular', 'small-world', 'random', 'scale-free']
LABEL = {'regular': 'Regular lattice', 'small-world': 'Small-world (WS)',
         'random': 'Random (ER)', 'scale-free': 'Scale-free (BA)'}
SEEDS = range(42, 142)
SEEDS_60 = range(42, 102)
