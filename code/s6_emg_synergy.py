"""s6. Muscle coupling networks and the synergy test, Supplementary Note 5.

Specified in advance. No value in the paper comes from this script, and the
results will be reported once the release has been analysed.

Why this script exists. s5 tests control entropy on ball flow. This one tests it
on the body, which is where the coordinative-structure claim lives. From the open
multi-muscle electromyography release named in the Methods (Atzori et al., 2014),
each recording is band-passed, rectified and reduced to a linear envelope, and
the coupling between two envelopes, by correlation or by coherence in a
task-relevant band, defines a weighted network with one node per muscle. The
shares follow by Eq 7, with control entropy, evenness and the effective count by
Eqs 4, 8 and 9, compared across conditions by evenness.

Synergy test. Non-negative matrix factorisation of the envelope matrix yields
synergies (d'Avella et al., 2003). The effective count is set against the number
of synergies extracted, and the leading muscle of the first synergy is compared
with the centrality hub, which tests the prediction of the earlier analysis
(Park, 2026b) that the hub muscles are the synergy leaders.

Acceptance. A difference is read as functional only when both conditions of
Supplementary Note 5 hold together. Control entropy falls with task error inside
the group in which the change appears, rho(H, E) > 0 (Eq S3), and the fall is
larger across the task-relevant muscles than across the task-irrelevant ones,
dH_rel < dH_irrel (Eq S4). A difference satisfying neither is recorded as
unstructured, and no functional reading is placed on it.

The coupling measure, the synergy extraction, the outcome variable, the rule that
splits the muscles and the two conditions are fixed in FIXED_BEFORE_RUN before
the routine runs, and are then applied uniformly.

load_emg_coupling() is completed when the release is analysed, and raises
NotImplementedError until then, which is intended. The analysis code is checked
on simulated data with

    python code/s6_emg_synergy.py --selftest

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
    'coupling_measure': None,    # correlation, or coherence in a named band
    'synergy_extraction': None,  # the NMF settings and how many synergies are kept
    'outcome_variable': None,    # the task error E recorded for each performer
    'relevance_rule': None,      # how a muscle is assigned to the task-relevant set
    'reference_group': None,     # the condition the change is read against
    'change_group': None,        # the condition in which the change is read
}


@dataclass
class MuscleCoupling:
    """One subject in one condition."""
    subject: str
    level: str                      # e.g. 'expert' or 'novice', or a task label
    coupling: np.ndarray            # M x M muscle coupling matrix, non-negative
    synergy_leader: int = None      # leading muscle of the first NMF synergy
    n_synergies: int = None         # number of synergies extracted
    outcome: float = None           # the task error E, as fixed in FIXED_BEFORE_RUN
    relevant: np.ndarray = None     # True for the task-relevant muscles, by the fixed rule


# ---------------------------------------------------------------------------
# The data loader, to be completed when the release is analysed
# ---------------------------------------------------------------------------
def load_emg_coupling(source):
    """Return a list of MuscleCoupling records from the open release.

    To be completed. Per subject and condition: band-pass, rectify, take the
    linear envelope, compute the M x M coupling between envelopes by the fixed
    measure, and return one record. Run NMF on the envelope matrix by the fixed
    settings, and store the leading muscle of the first synergy as synergy_leader
    and the number of synergies extracted as n_synergies. Attach the outcome E and
    the relevant-muscle mask, both as fixed in FIXED_BEFORE_RUN.

    Nothing here is fabricated. Until this is completed the script produces no
    validation number.
    """
    raise NotImplementedError(
        'Load the open release here. Until then, run with --selftest, whose numbers '
        'are NOT a result.')


# ---------------------------------------------------------------------------
# Analysis code
# ---------------------------------------------------------------------------
def eigenvector_shares(coupling):
    """Eq 7, from a non-negative coupling matrix."""
    A = np.asarray(coupling, dtype=float)
    A = np.maximum(A, A.T)
    np.fill_diagonal(A, 0.0)
    G = nx.from_numpy_array(A)
    c = nx.eigenvector_centrality_numpy(G, weight='weight')
    return share(np.abs(np.array([c[n] for n in G.nodes()], dtype=float)))


def entropy_of(rec):
    p = eigenvector_shares(rec.coupling)
    n = len(p)
    h = shannon_entropy(p)
    return {'N': n, 'H': h, 'J': evenness(h, n), 'D_eff': effective_dof(h),
            'hub': int(np.argmax(p)), 'p': p}


def bootstrap_ci(values, n_boot=10000, alpha=0.05, seed=0):
    """Percentile bootstrap interval for the mean."""
    v = np.asarray(values, dtype=float)
    rng = np.random.default_rng(seed)
    m = np.array([rng.choice(v, size=v.size, replace=True).mean() for _ in range(n_boot)])
    lo, hi = np.percentile(m, [100 * alpha / 2, 100 * (1 - alpha / 2)])
    return float(v.mean()), float(lo), float(hi)


def summarise(records):
    """Evenness across records, the contrast between two conditions, and the
    synergy test."""
    per = [entropy_of(r) for r in records]
    J = np.array([d['J'] for d in per], dtype=float)
    mean, lo, hi = bootstrap_ci(J)
    print('Muscle coupling networks, control entropy per subject and condition\n')
    print(f'muscle networks analysed : {len(records)}')
    print(f'J mean                   : {mean:.4f}   95% CI [{lo:.4f}, {hi:.4f}]')
    print(f'J range                  : {J.min():.4f} to {J.max():.4f}')
    levels = sorted({r.level for r in records})
    if len(levels) == 2:
        g0 = [d['J'] for d, r in zip(per, records) if r.level == levels[0]]
        g1 = [d['J'] for d, r in zip(per, records) if r.level == levels[1]]
        if len(g0) >= 3 and len(g1) >= 3:
            u, p = stats.mannwhitneyu(g0, g1, alternative='two-sided')
            rbc = 1.0 - 2.0 * u / (len(g0) * len(g1))
            print(f'\nCondition contrast: {levels[0]} (n={len(g0)}) vs {levels[1]} (n={len(g1)})')
            print(f'  median J: {np.median(g0):.4f} vs {np.median(g1):.4f}')
            print(f'  Mann-Whitney U = {u:.0f}, p = {p:.3g}, rank-biserial = {rbc:+.2f}')
    pairs = [(d['hub'], r.synergy_leader) for d, r in zip(per, records)
             if r.synergy_leader is not None]
    if pairs:
        agree = float(np.mean([h == s for h, s in pairs]))
        print(f'\nSynergy test: the centrality hub is the leading muscle of the first '
              f'synergy in {agree * 100:.0f}% of {len(pairs)} records.')
    syn = [(d['D_eff'], r.n_synergies) for d, r in zip(per, records)
           if r.n_synergies is not None]
    if len(syn) >= 3:
        de = np.array([a for a, _ in syn], dtype=float)
        ns = np.array([b for _, b in syn], dtype=float)
        rho, p = stats.spearmanr(de, ns)
        print(f'Synergy test: effective count against the number of synergies, '
              f'{len(syn)} records')
        print(f'  mean effective count {de.mean():.2f}, mean number of synergies '
              f'{ns.mean():.2f}, rank correlation {rho:+.2f} (p = {p:.3g})')
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
        changed group is more negative across the task-relevant muscles than
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


def plot_distribution(records, title='Muscle-network evenness J'):
    import matplotlib.pyplot as plt
    per = [entropy_of(r) for r in records]
    fig, ax = plt.subplots(figsize=(6.4, 4))
    for lv in sorted({r.level for r in records}):
        Js = [d['J'] for d, r in zip(per, records) if r.level == lv]
        ax.hist(Js, bins=12, alpha=0.6, label=f'{lv} (n={len(Js)})')
    ax.set_xlabel('evenness J')
    ax.set_ylabel('subjects')
    ax.set_title(title)
    ax.legend()
    fig.tight_layout()
    plt.show()


# ---------------------------------------------------------------------------
# Check on simulated data. NOT a result.
# ---------------------------------------------------------------------------
def _synthetic_records(n_per_group=30, n_muscles=12, seed=0):
    """Simulated muscle couplings for two conditions. NOT a result."""
    rng = np.random.default_rng(seed)
    out = []
    for level, conc in [('novice', 0.2), ('expert', 0.6)]:
        for i in range(n_per_group):
            base = rng.uniform(0.1, 0.4, size=(n_muscles, n_muscles))
            leader = int(rng.integers(0, n_muscles))
            base[leader, :] += conc
            base[:, leader] += conc
            A = (base + base.T) / 2
            np.fill_diagonal(A, 0.0)
            out.append(MuscleCoupling(f's{level}{i}', level, A, synergy_leader=leader,
                                      n_synergies=int(rng.integers(3, 7))))
    return out


def _acceptance_case(kind, n_per_group=30, n_muscles=12, seed=1):
    """Two simulated conditions built so the criterion has a known answer. NOT a
    result. In the functional case coupling gathers on a task-relevant muscle and
    E rises with H. In the unstructured case it gathers on a task-irrelevant
    muscle and E falls with H."""
    rng = np.random.default_rng(seed)
    relevant = np.arange(n_muscles) < n_muscles // 2
    target = 0 if kind == 'functional' else n_muscles - 1

    def coupling(boost):
        base = rng.uniform(0.1, 0.4, size=(n_muscles, n_muscles))
        base[target, :] += boost
        base[:, target] += boost
        A = (base + base.T) / 2
        np.fill_diagonal(A, 0.0)
        return A

    out = [MuscleCoupling(f'r{i}', 'reference', coupling(0.0), relevant=relevant)
           for i in range(n_per_group)]
    for i in range(n_per_group):
        rec = MuscleCoupling(f'c{i}', 'changed', coupling(rng.uniform(0.3, 1.0)),
                             relevant=relevant)
        h = entropy_of(rec)['H']
        rec.outcome = (h if kind == 'functional' else -h) + 0.001 * rng.standard_normal()
        out.append(rec)
    return out


def selftest(plot=False):
    print('>>> Check on SIMULATED data. These numbers are NOT a result. <<<\n')
    recs = _synthetic_records()
    summarise(recs)
    if plot:
        plot_distribution(recs, 'Muscle-network evenness J, simulated')
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
    records = load_emg_coupling(source='the release named in the Methods')
    per = summarise(records)
    print_acceptance(acceptance(records, per, FIXED_BEFORE_RUN['reference_group'],
                                FIXED_BEFORE_RUN['change_group']))
    if '--plot' in sys.argv:
        plot_distribution(records)
    return 0


if __name__ == '__main__':
    sys.exit(main())
