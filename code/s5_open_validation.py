"""s5. Football passing networks, Supplementary Note 5.

Specified in advance. No value in the paper comes from this script, and the
results will be reported once the release has been analysed.

Why this script exists. Football and handball are each published as a single
profile, so each gives one number with no interval (Supplementary Note 4). Open
event data remove that limit. Control entropy per team and match turns the value
into a distribution with a confidence interval, and lets levels of play be
contrasted as an axis of expertise.

Source. The StatsBomb open data release (StatsBomb, n.d.),

    https://github.com/statsbomb/open-data

free for research with attribution. For each team in each match, completed passes
define a weighted passer-to-receiver network with one node per player (Korte et
al., 2019). Eigenvector centrality gives the shares by Eq 7, and control entropy,
evenness and the effective count follow by Eqs 4, 8 and 9. Squad sizes differ, so
matches are compared by evenness J and not by raw H.

Contrast. Better organised play is predicted to concentrate ball flow more, so its
J sits lower. Networks are labelled by level of play, and the two-group contrast
is a Mann-Whitney test with a rank-biserial effect size, reported alongside a bit
difference between the group medians in the manner of Eq 12, here the log2 ratio
of the median compressions 1 - J.

Acceptance. A difference is read as functional only when both conditions of
Supplementary Note 5 hold together. Control entropy falls with task error inside
the group in which the change appears, rho(H, E) > 0 (Eq S3), and the fall is
larger across the task-relevant elements than across the task-irrelevant ones,
dH_rel < dH_irrel (Eq S4). A difference satisfying neither is recorded as
unstructured, and no functional reading is placed on it. The outcome variable,
the rule that splits the players and the two groups are fixed in
FIXED_BEFORE_RUN before the routine runs.

load_passing_networks() is completed when the release is analysed, and raises
NotImplementedError until then, which is intended. The analysis code is checked
on simulated data with

    python code/s5_open_validation.py --selftest

whose numbers are NOT a result. No files are written.
"""
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import networkx as nx
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
from core import share, shannon_entropy, evenness, effective_dof

# Fixed before the routine runs and recorded with it (Supplementary Note 5). The
# run on the release stops while any entry is None.
FIXED_BEFORE_RUN = {
    'outcome_variable': None,   # the task error E recorded for each team and match
    'relevance_rule': None,     # how a player is assigned to the task-relevant set
    'reference_group': None,    # the level of play the change is read against
    'change_group': None,       # the level of play in which the change is read
}


@dataclass
class PassingNetwork:
    """One team in one match."""
    match_id: str
    team: str
    level: str                     # e.g. 'international' or 'domestic'
    centrality: np.ndarray = None  # per-player eigenvector centrality, if precomputed
    adjacency: np.ndarray = None   # weighted pass matrix, used when centrality is None
    outcome: float = None          # the task error E, as fixed in FIXED_BEFORE_RUN
    relevant: np.ndarray = None    # True for the task-relevant players, by the fixed rule


# ---------------------------------------------------------------------------
# The data loader, to be completed when the release is analysed
# ---------------------------------------------------------------------------
def load_passing_networks(source_dir):
    """Return a list of PassingNetwork records from the open release.

    To be completed. For StatsBomb, read the event JSON per match, keep completed
    passes, build the weighted passer -> receiver matrix per team, and return one
    record per team per match with its level label, its outcome E and its
    relevant-player mask, both as fixed in FIXED_BEFORE_RUN.

    Nothing here is fabricated. Until this is completed the script produces no
    validation number.
    """
    raise NotImplementedError(
        'Load the open release here. Until then, run with --selftest, whose numbers '
        'are NOT a result.')


# ---------------------------------------------------------------------------
# Analysis code
# ---------------------------------------------------------------------------
def eigenvector_shares(adjacency):
    """Eq 7, from a weighted adjacency matrix."""
    A = np.asarray(adjacency, dtype=float)
    G = nx.from_numpy_array(np.maximum(A, A.T))   # symmetrise for centrality
    c = nx.eigenvector_centrality_numpy(G, weight='weight')
    return share(np.abs(np.array([c[n] for n in G.nodes()], dtype=float)))


def entropy_of(net):
    p = share(net.centrality) if net.centrality is not None else eigenvector_shares(net.adjacency)
    n = len(p)
    h = shannon_entropy(p)
    return {'N': n, 'H': h, 'J': evenness(h, n), 'D_eff': effective_dof(h), 'p': p}


def bootstrap_ci(values, n_boot=10000, alpha=0.05, seed=0):
    """Percentile bootstrap interval for the mean."""
    v = np.asarray(values, dtype=float)
    rng = np.random.default_rng(seed)
    m = np.array([rng.choice(v, size=v.size, replace=True).mean() for _ in range(n_boot)])
    lo, hi = np.percentile(m, [100 * alpha / 2, 100 * (1 - alpha / 2)])
    return float(v.mean()), float(lo), float(hi)


def rank_contrast(a, b):
    """Mann-Whitney contrast and the rank-biserial effect size between two groups."""
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    u, p = stats.mannwhitneyu(a, b, alternative='two-sided')
    return {'U': float(u), 'p': float(p), 'rank_biserial': float(1.0 - 2.0 * u / (a.size * b.size))}


def summarise(nets):
    """Evenness across networks, and the contrast between two levels of play."""
    per = [entropy_of(n) for n in nets]
    J = np.array([d['J'] for d in per], dtype=float)
    mean, lo, hi = bootstrap_ci(J)
    print('Football passing networks, control entropy per team and match\n')
    print(f'networks analysed : {len(nets)}')
    print(f'J mean            : {mean:.4f}   95% CI [{lo:.4f}, {hi:.4f}]')
    print(f'J range           : {J.min():.4f} to {J.max():.4f}')
    print(f"median share of positions retained (2^H / N): "
          f"{np.median([d['D_eff'] / d['N'] for d in per]):.3f}")
    levels = sorted({n.level for n in nets})
    if len(levels) == 2:
        g0 = [d['J'] for d, n in zip(per, nets) if n.level == levels[0]]
        g1 = [d['J'] for d, n in zip(per, nets) if n.level == levels[1]]
        if len(g0) >= 3 and len(g1) >= 3:
            res = rank_contrast(g0, g1)
            m0, m1 = np.median(g0), np.median(g1)
            print(f'\nLevel contrast: {levels[0]} (n={len(g0)}) vs {levels[1]} (n={len(g1)})')
            print(f'  median J: {m0:.4f} vs {m1:.4f}')
            print(f"  Mann-Whitney U = {res['U']:.0f}, p = {res['p']:.3g}, "
                  f"rank-biserial = {res['rank_biserial']:+.2f}")
            if 0 < m0 < 1 and 0 < m1 < 1:
                print(f'  bit difference of the medians, log2((1 - J) {levels[1]} / '
                      f'(1 - J) {levels[0]}) = {np.log2((1 - m1) / (1 - m0)):+.3f} bits')
            print('  A lower J means ball flow is more concentrated on a hub.')
    return per


def subset_entropy(p, mask):
    """Control entropy, in bits, of the shares of one set of elements, renormalised
    over that set. Used for the task-relevant and task-irrelevant sets of Eq S4."""
    q = np.asarray(p, dtype=float)[np.asarray(mask, dtype=bool)]
    return shannon_entropy(q) if q.size and q.sum() > 0 else float('nan')


def acceptance(records, per, reference, changed):
    """Eqs S3 and S4 of Supplementary Note 5, read together.

    S3. Inside the group in which the change appears, control entropy falls with
        task error, rho(H, E) > 0, as a rank correlation over the records that
        carry E.
    S4. The change in median control entropy from the reference group to the
        changed group is more negative across the task-relevant elements than
        across the task-irrelevant ones, dH_rel < dH_irrel.
    Both together: functional. Neither: unstructured, no functional reading.
    """
    ch = [(d, r) for d, r in zip(per, records) if r.level == changed]
    ref = [(d, r) for d, r in zip(per, records) if r.level == reference]
    he = [(d['H'], r.outcome) for d, r in ch if r.outcome is not None]
    rho, p_rho = (stats.spearmanr([h for h, _ in he], [e for _, e in he])
                  if len(he) >= 3 else (float('nan'), float('nan')))

    def median_subset(group, relevant_side):
        vals = []
        for d, r in group:
            if r.relevant is None:
                continue
            mask = np.asarray(r.relevant, dtype=bool)
            vals.append(subset_entropy(d['p'], mask if relevant_side else ~mask))
        vals = [v for v in vals if np.isfinite(v)]
        return float(np.median(vals)) if vals else float('nan')

    d_rel = median_subset(ch, True) - median_subset(ref, True)
    d_irr = median_subset(ch, False) - median_subset(ref, False)
    if not (np.isfinite(rho) and np.isfinite(d_rel) and np.isfinite(d_irr)):
        verdict = 'undetermined, E or the relevant set is missing'
        s3 = s4 = None
    else:
        s3, s4 = bool(rho > 0), bool(d_rel < d_irr)
        verdict = ('functional' if s3 and s4 else
                   'unstructured, no functional reading' if not (s3 or s4) else
                   'one condition only, not read as functional')
    return {'rho': float(rho), 'p_rho': float(p_rho), 'n_E': len(he), 'd_rel': d_rel,
            'd_irr': d_irr, 'S3': s3, 'S4': s4, 'verdict': verdict}


def print_acceptance(a):
    def word(x):
        return 'not computed' if x is None else ('holds' if x else 'fails')
    print('\nAcceptance criterion, Supplementary Note 5')
    print(f"  Eq S3  rho(H, E) = {a['rho']:+.2f} over {a['n_E']} records "
          f"(p = {a['p_rho']:.3g})  -> {word(a['S3'])}")
    print(f"  Eq S4  dH_rel = {a['d_rel']:+.4f} bits, dH_irrel = {a['d_irr']:+.4f} bits"
          f"  -> {word(a['S4'])}")
    print(f"  reading: {a['verdict']}")


def plot_distribution(nets, title='Passing-network evenness J'):
    import matplotlib.pyplot as plt
    per = [entropy_of(n) for n in nets]
    fig, ax = plt.subplots(figsize=(6.4, 4))
    for lv in sorted({n.level for n in nets}):
        Js = [d['J'] for d, n in zip(per, nets) if n.level == lv]
        ax.hist(Js, bins=15, alpha=0.6, label=f'{lv} (n={len(Js)})')
    ax.set_xlabel('evenness J')
    ax.set_ylabel('matches')
    ax.set_title(title)
    ax.legend()
    fig.tight_layout()
    plt.show()


# ---------------------------------------------------------------------------
# Check on simulated data. NOT a result.
# ---------------------------------------------------------------------------
def _synthetic_networks(n_per_group=40, seed=0):
    """Simulated passing networks for two levels. NOT a result."""
    rng = np.random.default_rng(seed)
    nets = []
    for level, m in [('domestic', 2), ('international', 4)]:
        for i in range(n_per_group):
            n = int(rng.integers(10, 12))
            G = nx.barabasi_albert_graph(n, m, seed=int(rng.integers(0, 10 ** 6)))
            A = nx.to_numpy_array(G) * rng.uniform(1, 5, size=(n, n))
            nets.append(PassingNetwork(f'm{level}{i}', f'team{i}', level, adjacency=A))
    return nets


def _acceptance_case(kind, n_per_group=30, n=11, seed=1):
    """Two simulated groups built so the criterion has a known answer. NOT a result.

    In the functional case the changed group concentrates on a task-relevant
    player and E rises with H. In the unstructured case it concentrates on a
    task-irrelevant player and E falls with H.
    """
    rng = np.random.default_rng(seed)
    relevant = np.arange(n) < 5
    target = 0 if kind == 'functional' else n - 1
    nets = []
    for i in range(n_per_group):
        c = 1.0 + 0.2 * rng.uniform(size=n)
        nets.append(PassingNetwork(f'r{i}', f'team{i}', 'reference', centrality=c,
                                   relevant=relevant))
    for i in range(n_per_group):
        c = 1.0 + 0.2 * rng.uniform(size=n)
        c[target] += rng.uniform(0.5, 2.0)
        h = shannon_entropy(c)
        e = (h if kind == 'functional' else -h) + 0.01 * rng.standard_normal()
        nets.append(PassingNetwork(f'c{i}', f'team{i}', 'changed', centrality=c,
                                   outcome=e, relevant=relevant))
    return nets


def selftest(plot=False):
    print('>>> Check on SIMULATED data. These numbers are NOT a result. <<<\n')
    nets = _synthetic_networks()
    summarise(nets)
    if plot:
        plot_distribution(nets, 'Passing-network evenness J, simulated')
    print('\nThe acceptance criterion reads back the case it was given:')
    ok_all = True
    for kind, expect in (('functional', 'functional'),
                         ('unstructured', 'unstructured, no functional reading')):
        cases = _acceptance_case(kind)
        a = acceptance(cases, [entropy_of(c) for c in cases], 'reference', 'changed')
        ok = a['verdict'] == expect
        ok_all &= ok
        print(f"  built as {kind:<12} -> rho {a['rho']:+.2f}, dH_rel {a['d_rel']:+.3f}, "
              f"dH_irrel {a['d_irr']:+.3f} -> {a['verdict']:<36} [{'PASS' if ok else 'FAIL'}]")
    return ok_all


def main():
    if '--selftest' in sys.argv:
        return 0 if selftest(plot='--plot' in sys.argv) else 1
    unset = [k for k, v in FIXED_BEFORE_RUN.items() if v is None]
    if unset:
        raise ValueError('Fix these before the routine runs (Supplementary Note 5): '
                         + ', '.join(unset))
    nets = load_passing_networks(source_dir='statsbomb/open-data')
    per = summarise(nets)
    print_acceptance(acceptance(nets, per, FIXED_BEFORE_RUN['reference_group'],
                                FIXED_BEFORE_RUN['change_group']))
    if '--plot' in sys.argv:
        plot_distribution(nets)
    return 0


if __name__ == '__main__':
    sys.exit(main())
