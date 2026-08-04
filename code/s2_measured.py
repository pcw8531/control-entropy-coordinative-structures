"""s2. The five measured coordinative structures and the expertise data.

Behind SI Appendix section 4, section 5, Table S6 and Table S7, and behind the
upper row of Figure 4. Everything here is read from a released file or from a
published table transcribed into one; nothing was collected for this paper.
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from core import (share, shannon_entropy, max_entropy, evenness, effective_dof,
                  gini, outcome_entropy_gaussian, delta_outcome_entropy,
                  predictability_ratio)
from loaders import read_csv, read_matrix

POSITIONS = ['PG', 'SG', 'SF', 'PF', 'CN']


def eigen_shares_matrix(A):
    """Eq 3 on a weighted adjacency matrix."""
    w, v = np.linalg.eigh(A)
    return share(np.abs(v[:, np.argmax(w)]))


def sport_shares(sport):
    cent = read_csv('team_sport_position_centrality.csv')
    return share(cent[cent.sport == sport].centrality.values)


def basketball_profiles():
    """Sixteen teams, each its own profile. Published team by team."""
    cent = read_csv('team_sport_position_centrality.csv')
    q = cent[cent.sport == 'Basketball']
    teams = sorted(q.unit.unique())
    prof = np.array([[q[(q.unit == t) & (q.position == p)].centrality.values[0]
                      for p in POSITIONS] for t in teams])
    return teams, prof


def cortical_shares(from_matrix=True):
    """Absolute correlation with a zero diagonal, no threshold imposed."""
    if from_matrix:
        A = np.abs(read_matrix())
        np.fill_diagonal(A, 0.0)
        return eigen_shares_matrix(A)
    return read_csv('cortical_coupling_shares.csv').share_p_i_eigenvector.values


def musculoskeletal_shares():
    return share(read_csv('musculoskeletal_coupling_shares.csv').bones_linked.values)


def musculoskeletal_null_shares():
    return share(read_csv('musculoskeletal_random_null_shares.csv').bones_linked.values)


def expertise():
    e = read_csv('expertise_outcome_entropy.csv')
    return (e[e.group == 'Novice'], e[e.group == 'Expert'])


def table_s6():
    """Control entropy of the five measured structures."""
    teams, prof = basketball_profiles()
    per_team = np.array([shannon_entropy(share(r)) for r in prof])
    rows = []

    def row(name, p, H=None, note=''):
        N = len(p)
        h = shannon_entropy(p) if H is None else H
        rows.append(dict(system=name, N=N, H=h, ceiling=max_entropy(N),
                         J=evenness(h, N), D_eff=effective_dof(h),
                         hub_share=100 * max(p), uniform_share=100 / N, note=note))

    row('Football', sport_shares('Football'))
    row('Basketball', share(prof.mean(axis=0)), H=per_team.mean(),
        note=f'mean of sixteen team entropies, SD {per_team.std(ddof=1):.4f}; '
             f'the profile drawn is the averaged one at '
             f'{shannon_entropy(share(prof.mean(axis=0))):.4f} bits')
    row('Handball', sport_shares('Handball'))
    row('Cortical activity', cortical_shares())
    row('Musculoskeletal', musculoskeletal_shares())
    return rows


def table_s7():
    nov, exp = expertise()
    n_h = nov.outcome_entropy_h_bits.values
    e_h = exp.outcome_entropy_h_bits.values
    return dict(
        novice_mean_sd=float(nov.outcome_sd_cm.mean()),
        expert_mean_sd=float(exp.outcome_sd_cm.mean()),
        novice_h=float(n_h.mean()), novice_h_sd=float(n_h.std(ddof=1)),
        expert_h=float(e_h.mean()), expert_h_sd=float(e_h.std(ddof=1)),
        novice_range=(float(n_h.min()), float(n_h.max())),
        expert_range=(float(e_h.min()), float(e_h.max())),
        difference_group_sd=delta_outcome_entropy(
            nov.outcome_sd_cm.mean(), exp.outcome_sd_cm.mean()),
        difference_per_participant=float(n_h.mean() - e_h.mean()),
        ratio=predictability_ratio(float(n_h.mean() - e_h.mean())),
        non_overlap=float(n_h.min() - e_h.max()))


def main():
    print('=' * 96)
    print('SI TABLE S6   five measured coordinative structures')
    print('=' * 96)
    print(f'{"system":<20}{"N":>5}{"H":>10}{"log2 N":>10}{"J":>9}{"2^H":>10}'
          f'{"hub %":>8}{"uniform %":>11}')
    for r in table_s6():
        print(f'{r["system"]:<20}{r["N"]:>5}{r["H"]:>10.4f}{r["ceiling"]:>10.4f}'
              f'{r["J"]:>9.4f}{r["D_eff"]:>10.2f}{r["hub_share"]:>8.1f}'
              f'{r["uniform_share"]:>11.1f}')
        if r['note']:
            print(f'    note: {r["note"]}')
    print()
    print('Volleyball is uniform by the rotation rule, so its control entropy is '
          f'log2 6 = {max_entropy(6):.4f} bits at evenness exactly one.')

    print()
    print('=' * 96)
    print('SI TABLE S7   outcome entropy in the expertise data')
    print('=' * 96)
    t = table_s7()
    print(f'   novice  mean SD {t["novice_mean_sd"]:.3f} cm   '
          f'h {t["novice_h"]:.3f} +/- {t["novice_h_sd"]:.3f} bits   '
          f'range {t["novice_range"][0]:.2f} to {t["novice_range"][1]:.2f}')
    print(f'   expert  mean SD {t["expert_mean_sd"]:.3f} cm   '
          f'h {t["expert_h"]:.3f} +/- {t["expert_h_sd"]:.3f} bits   '
          f'range {t["expert_range"][0]:.2f} to {t["expert_range"][1]:.2f}')
    print(f'   difference from the group mean SDs      '
          f'{t["difference_group_sd"]:.3f} bits')
    print(f'   difference of the per-participant means '
          f'{t["difference_per_participant"]:.3f} bits')
    print(f'   predictability ratio {t["ratio"]:.3f}, '
          f'groups apart by {t["non_overlap"]:.3f} bits')
    print()
    print('   The mean absolute error column of Table S7 is not derived here.')
    print('   It is not present in expertise_outcome_entropy.csv or in any other')
    print('   released file, and the source publication for this dataset is')
    print('   recorded in SI section 5 as still to be supplied.')

    print()
    print('=' * 96)
    print('SI SECTION 4   the cortical threshold sweep')
    print('=' * 96)
    sweep = read_csv('cortical_threshold_sweep.csv')
    print(sweep.to_string(index=False))
    print()
    print('   Note on the compression column. It is computed against log2 of the')
    print('   number of regions that remain connected at that threshold, not')
    print('   against log2 of 100. The two agree from 1200 links upward, where')
    print('   every region stays connected. At 300 links a fixed ceiling of')
    print('   log2 100 would give 0.1762 rather than 0.1438, and at 600 links')
    print('   0.1247 rather than 0.1189.')


if __name__ == '__main__':
    main()
