"""s4. Learning dynamics and the order of change, Supplementary Note 5.

Specified in advance. No value in the paper comes from this script, and the
trajectory will be reported once it has been run with the earlier model's code.

Why a fixed node set can still grow a hub. A regular, random or small-world
construction cannot produce a hub at all. Every node is there from the start and
every edge is placed by the same rule, so no element can get ahead of another.
Growth with preferential attachment answers this for networks that add nodes. A
body does not add joints, so the growth here is in the coupling, not in the node
set, and the advantage that breaks the equality is functional, co-activation in
the task rather than seniority.

What it computes. Three curves along one practice trajectory of the earlier
model (Park, 2026b):

  1. H(t), control entropy on the strength shares, expected to fall as weight
     concentrates on hubs.
  2. r(t), the strength assortativity of the weighted structure. Degree
     assortativity is the wrong curve here, because the edge set is fixed during
     learning and only the weights move, so it does not change. Strength
     assortativity does.
  3. The order-of-change test, which compares the onset of coupling
     concentration with the onset of activation similarity, onset being the
     first step at which a signal reaches half its baseline-to-final rise. If
     coupling leads the process is influence-like, and if similarity leads it is
     selection-like.

Why one function is left to the earlier model's code. Supplementary Eq S2,

    W_ij(t+1) = W_ij(t) + eta A_i(t) A_j(t) - delta W_ij(t),
    eta = 0.01, delta = 0.001,

depends on the activation state A_i(t) of the earlier model's cascade, which is
not a library routine. A fresh re-implementation would give traces that differ
from the ones behind the published H_W = 5.5, and every curve downstream would
then be wrong in a way that is hard to see. So learning_trajectory() is left to
call the released routine, in simulations/04_coupled_learning.py of

    https://github.com/pcw8531/Coordinative-structures-scale-free-networks
    (Zenodo 10.5281/zenodo.20694466)

and raises NotImplementedError until then, which is intended.

The analysis code is checked on simulated data with

    python code/s4_learning_entropy.py --selftest

whose numbers are NOT a result. Add --plot to draw C(t) and S(t), since both
curves are to be inspected before the verdict is read. No files are written.
"""
import sys
from pathlib import Path

import numpy as np
import networkx as nx

sys.path.insert(0, str(Path(__file__).resolve().parent))
from core import shannon_entropy, effective_dof, evenness, strength_shares, max_entropy

# The earlier model's skeleton and schedule
N = 100
M_BA = 3
SEED = 42
TIMEPOINTS = [0, 20, 56, 100, 300, 1000]   # t_0.5 = 56 is the earlier paper's value
ETA = 0.01
DELTA = 0.001
EPS = 1e-9   # weight above which a pair counts as coupled


# ---------------------------------------------------------------------------
# To be completed with the earlier model's released routine
# ---------------------------------------------------------------------------
def learning_trajectory(G, timepoints, n_steps, seed=SEED):
    """Run Supplementary Eq S2 on the skeleton G and return the trajectory.

    To be completed with a call into the earlier model's released routine, for
    instance the learning loop of simulations/04_coupled_learning.py, recording
    the weights, the node strengths and the activations as it runs. The
    requirement is that A_i(t) comes from the released cascade code, not from a
    fresh re-implementation.

    Must return a dict with:
        'W_at'            : {t: (N, N) weight matrix} for every t in timepoints
        'strength_time'   : (n_steps, N) node strengths at each fine step
        'activation_time' : (n_steps, N) activations A_i at each fine step
    """
    raise NotImplementedError(
        "Link the earlier model's Supplementary Eq S2 routine here. Until then, run "
        "with --selftest, whose numbers are NOT a result.")


# ---------------------------------------------------------------------------
# Analysis code
# ---------------------------------------------------------------------------
def strength_assortativity(W, eps=EPS):
    """Weighted analogue of Eq 11: correlation of node strength across edges.

    Built on the support of W, the pairs with weight above eps. Returns nan if
    the support has no usable variance.
    """
    W = np.asarray(W, dtype=float)
    n = W.shape[0]
    s = W.sum(axis=1)
    G = nx.Graph()
    G.add_nodes_from(range(n))
    iu = np.triu_indices(n, k=1)
    for i, j in zip(*iu):
        if W[i, j] > eps or W[j, i] > eps:
            G.add_edge(int(i), int(j))
    if G.number_of_edges() == 0:
        return float('nan')
    nx.set_node_attributes(G, {i: float(s[i]) for i in range(n)}, 'strength')
    try:
        return float(nx.numeric_assortativity_coefficient(G, 'strength'))
    except (ValueError, ZeroDivisionError):
        return float('nan')


def concentration_signal(strength_time):
    """C(t) = 1 - J(t). Rises as coupling gathers on a few nodes."""
    S = np.asarray(strength_time, dtype=float)
    logn = max_entropy(S.shape[1])
    out = np.zeros(S.shape[0])
    for k in range(S.shape[0]):
        out[k] = 0.0 if S[k].sum() <= 0 else 1.0 - shannon_entropy(S[k]) / logn
    return out


def similarity_signal(activation_time, window=25):
    """S(t), the mean pairwise activation correlation in a sliding window."""
    A = np.asarray(activation_time, dtype=float)
    T, n = A.shape
    out = np.zeros(T)
    iu = np.triu_indices(n, k=1)
    for k in range(T):
        block = A[max(0, k - window + 1):k + 1]
        if block.shape[0] < 3:
            out[k] = np.nan
            continue
        with np.errstate(invalid='ignore', divide='ignore'):
            out[k] = np.nanmean(np.corrcoef(block.T)[iu])
    return out


def onset_time(signal, frac=0.5):
    """First index at which a signal reaches frac of its baseline-to-final rise."""
    x = np.asarray(signal, dtype=float)
    good = np.isfinite(x)
    if good.sum() < 2:
        return None
    base = np.nanmin(x)
    final = np.nanmedian(x[np.where(good)[0][-5:]])
    if final <= base:
        return None
    idx = np.where(good & (x >= base + frac * (final - base)))[0]
    return int(idx[0]) if idx.size else None


def precedence_verdict(strength_time, activation_time):
    """Does coupling concentration or activation similarity move first?"""
    tC = onset_time(concentration_signal(strength_time))
    tS = onset_time(similarity_signal(activation_time))
    if tC is None or tS is None:
        v = 'undetermined (a signal did not rise)'
    elif tC < tS:
        v = 'coupling leads -> influence-like'
    elif tS < tC:
        v = 'similarity leads -> selection-like'
    else:
        v = 'simultaneous within resolution'
    return {'onset_coupling': tC, 'onset_similarity': tS, 'verdict': v}


def report_trajectory(traj):
    """Print H(t), J(t), D_eff(t), r(t) and the order-of-change verdict."""
    print('H(t), J(t), D_eff(t) and strength assortativity r(t) along the trajectory\n')
    hdr = f"{'t':>6}{'H(t) bits':>12}{'J':>9}{'D_eff':>9}{'r(t)':>10}{'dH from t=0':>14}"
    print(hdr)
    print('-' * len(hdr))
    h0 = None
    for t in TIMEPOINTS:
        W = traj['W_at'][t]
        h = shannon_entropy(strength_shares(W))
        r = strength_assortativity(W)
        h0 = h if h0 is None else h0
        rt = '  nan' if np.isnan(r) else f'{r:+.3f}'
        print(f'{t:>6}{h:>12.4f}{evenness(h, N):>9.4f}{effective_dof(h):>9.1f}'
              f'{rt:>10}{h - h0:>14.4f}')
    pv = precedence_verdict(traj['strength_time'], traj['activation_time'])
    print('\nOrder-of-change test:')
    print(f"  onset of coupling concentration : step {pv['onset_coupling']}")
    print(f"  onset of activation similarity  : step {pv['onset_similarity']}")
    print(f"  reading: {pv['verdict']}")
    print('  The onset compares each signal against its own range, so it is sensitive')
    print('  to signal shape as well as timing. Plot C(t) and S(t), with --plot, and')
    print('  inspect both before the verdict is read.')


def plot_signals(traj, title=''):
    """Draw C(t) and S(t), each scaled to its own range, with both onsets marked."""
    import matplotlib.pyplot as plt
    C = concentration_signal(traj['strength_time'])
    S = similarity_signal(traj['activation_time'])
    pv = precedence_verdict(traj['strength_time'], traj['activation_time'])

    def nrm(v):
        return (v - np.nanmin(v)) / (np.nanmax(v) - np.nanmin(v))

    fig, ax = plt.subplots(figsize=(6.6, 4))
    ax.plot(nrm(C), color='#4C72B0', label='coupling concentration C(t)')
    ax.plot(nrm(S), color='#C44E52', label='activation similarity S(t)')
    for t, c in ((pv['onset_coupling'], '#4C72B0'), (pv['onset_similarity'], '#C44E52')):
        if t is not None:
            ax.axvline(t, ls=':', color=c)
    ax.set_xlabel('step')
    ax.set_ylabel('rise, scaled to its own range')
    ax.set_title(title)
    ax.legend(fontsize=7)
    fig.tight_layout()
    plt.show()


# ---------------------------------------------------------------------------
# Check on simulated data. NOT a result.
# ---------------------------------------------------------------------------
def _synthetic_trajectory(lead='coupling', n_steps=400, seed=0):
    """A simulated trajectory with a known lead direction. NOT a result.

    Both signals are built as smooth ramps of the same shape, one delayed against
    the other, so the only thing that differs between the two cases is which
    moves first. That makes this a test of the detector, since the onset compares
    each signal against its own range and two signals of different shape could
    otherwise give a verdict that reflects shape rather than order.
    """
    rng = np.random.default_rng(seed)
    A0 = nx.to_numpy_array(nx.barabasi_albert_graph(N, M_BA, seed=SEED))
    deg = A0.sum(axis=1)
    w = deg / deg.max()
    x = np.arange(n_steps) / n_steps

    def ramp(s):
        return np.clip((x - s) / 0.35, 0.0, 1.0)

    s_c, s_s = (0.05, 0.45) if lead == 'coupling' else (0.45, 0.05)
    conc, simi = ramp(s_c), ramp(s_s)

    def W_of(c):
        return A0 * (1.0 + 12.0 * c * np.outer(deg, deg) / deg.max() ** 2)

    st = np.zeros((n_steps, N))
    at = np.zeros((n_steps, N))
    u = rng.standard_normal(n_steps)          # one shared driver
    for k in range(n_steps):
        st[k] = W_of(conc[k]).sum(axis=1)
        # as simi rises the nodes follow the common driver, so their windowed
        # correlation rises smoothly from noise toward one
        at[k] = simi[k] * u[k] * w + (1.0 - simi[k]) * rng.standard_normal(N)
    W_at = {t: W_of(conc[min(int(t / 1000 * (n_steps - 1)), n_steps - 1)])
            for t in TIMEPOINTS}
    return {'W_at': W_at, 'strength_time': st, 'activation_time': at}


def selftest(plot=False):
    print('>>> Check on SIMULATED data. These numbers are NOT a result. <<<')
    print('The only claim is that the detector reads back the direction it was given.\n')
    ok_all = True
    for lead in ('coupling', 'similarity'):
        tr = _synthetic_trajectory(lead)
        pv = precedence_verdict(tr['strength_time'], tr['activation_time'])
        ok = ('coupling leads' in pv['verdict']) == (lead == 'coupling')
        ok_all &= ok
        print(f"built with {lead:>10} leading -> onset C = {pv['onset_coupling']:>4}, "
              f"onset S = {pv['onset_similarity']:>4} -> {pv['verdict']:<38} "
              f"[{'PASS' if ok else 'FAIL'}]")
        if plot:
            plot_signals(tr, f'Simulated, built so that {lead} leads')
    print()
    report_trajectory(_synthetic_trajectory('coupling'))
    return ok_all


def main():
    if '--selftest' in sys.argv:
        return 0 if selftest(plot='--plot' in sys.argv) else 1
    G0 = nx.barabasi_albert_graph(N, M_BA, seed=SEED)
    print(f'BA skeleton of the earlier model: N = {N}, m = {M_BA}, seed = {SEED}')
    print(f'Supplementary Eq S2: eta = {ETA}, delta = {DELTA}\n')
    traj = learning_trajectory(G0, TIMEPOINTS, n_steps=1000, seed=SEED)
    report_trajectory(traj)
    if '--plot' in sys.argv:
        plot_signals(traj, 'Coupling concentration and activation similarity over practice')
    return 0


if __name__ == '__main__':
    sys.exit(main())
