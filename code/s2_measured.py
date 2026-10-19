"""s2. The five measured coordinative structures and the expertise data.

Behind Supplementary Tables 5 and 6, Supplementary Notes 2 and 3, and the upper
row of Figure 4. Everything here is read from a released file or from a
published table transcribed into one; nothing was collected for this paper.
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from core import (share, shannon_entropy, max_entropy, evenness, effective_dof,
                  gini, outcome_entropy_gaussian, delta_outcome_entropy,
                  predictability_ratio, two_level_entropy)
from loaders import read_csv, read_matrix

POSITIONS = ['PG', 'SG', 'SF', 'PF', 'CN']


def eigen_shares_matrix(A):
    """Eq 7 on a weighted adjacency matrix."""
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


def table_measured():
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
    hubs = [share(r).max() for r in prof]
    pg_leads = sum(int(np.argmax(r)) == POSITIONS.index('PG') for r in prof)
    row('Basketball', share(prof.mean(axis=0)), H=per_team.mean(),
        note=f'mean of sixteen team entropies, SD {per_team.std(ddof=1):.4f}; '
             f'the profile drawn is the averaged one at '
             f'{shannon_entropy(share(prof.mean(axis=0))):.4f} bits; '
             f'the sixteen hub shares average {100 * np.mean(hubs):.1f} per cent, '
             f'and the point guard leads {pg_leads} of the {len(prof)} teams')
    row('Handball', sport_shares('Handball'))
    row('Cortical activity', cortical_shares())
    row('Musculoskeletal', musculoskeletal_shares())
    return rows


def leader_scan():
    """Supplementary Note 2. Eq 15 of the paper (Eq 9 in core.py) over the number
    of leaders h, as the predicted control entropy minus the measured one, in
    bits. Basketball is formed team by team and then averaged."""
    out = {}
    for sport in ('Football', 'Handball'):
        p = sport_shares(sport)
        q = np.sort(p)[::-1]
        H = shannon_entropy(p)
        out[sport] = [two_level_entropy(len(q), h, q[:h].sum()) - H
                      for h in range(1, len(q))]
    teams, prof = basketball_profiles()
    H = np.mean([shannon_entropy(share(r)) for r in prof])
    out['Basketball'] = [float(np.mean([two_level_entropy(5, h, np.sort(share(r))[::-1][:h].sum())
                                        for r in prof])) - H for h in range(1, 5)]
    out['handball_three'] = float(np.sort(sport_shares('Handball'))[::-1][:3].sum())
    return out


def cortical_file():
    """Supplementary Note 2. What the released cortical matrix holds, and how
    closely the eigenvector and strength weightings order the hundred regions."""
    from scipy.stats import spearmanr
    M = read_matrix()
    off = M[~np.eye(len(M), dtype=bool)]
    A = np.abs(M)
    np.fill_diagonal(A, 0.0)
    rho = spearmanr(cortical_shares(), share(A.sum(axis=1)))[0]
    return dict(diag_min=float(np.diag(M).min()), diag_max=float(np.diag(M).max()),
                off_min=float(off.min()), off_max=float(off.max()),
                negative=float(100 * (off < 0).mean()), rho=float(rho))


def muscle_columns(draws=20000, seed=1):
    """Supplementary Note 2. The two released columns of fig2e.csv, read with the
    first bin at degree two, and nulls rebuilt at the measured attachment count.

    Each rebuilt null gives every muscle two attachments, an origin and an
    insertion, and spreads the remaining attachments over the muscles at random.
    """
    md = read_csv('muscle_degree_counts.csv')
    n = int(md.muscles_measured.sum())
    att_measured = int((md.bones_linked * md.muscles_measured).sum())
    att_null = int((md.bones_linked * md.muscles_random_null).sum())
    extra = att_measured - 2 * n
    rng = np.random.default_rng(seed)
    rebuilt = np.empty(draws)
    for b in range(draws):
        p = share(2 + np.bincount(rng.integers(0, n, extra), minlength=n))
        rebuilt[b] = 1 - evenness(shannon_entropy(p), n)
    measured = 1 - evenness(shannon_entropy(musculoskeletal_shares()), n)
    released = 1 - evenness(shannon_entropy(musculoskeletal_null_shares()), n)
    return dict(muscles=n, att_measured=att_measured, att_null=att_null, extra=extra,
                draws=draws, rebuilt=float(rebuilt.mean()), measured=measured,
                released=released, factor_released=measured / released,
                factor_rebuilt=measured / float(rebuilt.mean()))


def table_expertise():
    """Eq 12 on the twenty published participant values, with the exact F test."""
    from scipy.stats import f as f_dist
    nov, exp = expertise()
    x_n = nov.radial_error_cm.values
    x_e = exp.radial_error_cm.values
    sd_n, sd_e = float(x_n.std(ddof=1)), float(x_e.std(ddof=1))
    diff = delta_outcome_entropy(sd_n, sd_e)
    dfn, dfd = len(x_n) - 1, len(x_e) - 1
    F = (sd_n / sd_e) ** 2
    p = 2 * min(f_dist.sf(F, dfn, dfd), f_dist.cdf(F, dfn, dfd))
    lo, hi = F / f_dist.ppf(0.975, dfn, dfd), F / f_dist.ppf(0.025, dfn, dfd)
    return dict(
        novice_mean=float(x_n.mean()), expert_mean=float(x_e.mean()),
        novice_sd=sd_n, expert_sd=sd_e,
        novice_h=float(outcome_entropy_gaussian(x_n.std(ddof=0))),
        expert_h=float(outcome_entropy_gaussian(x_e.std(ddof=0))),
        novice_h_column=float(nov.outcome_information_bits.mean()),
        expert_h_column=float(exp.outcome_information_bits.mean()),
        difference=diff, ratio=predictability_ratio(diff),
        F=float(F), df=(dfn, dfd), p=float(p),
        diff_ci=(0.5 * float(np.log2(lo)), 0.5 * float(np.log2(hi))),
        ratio_ci=(float(np.sqrt(lo)), float(np.sqrt(hi))))


def main():
    print('=' * 96)
    print('SUPPLEMENTARY TABLE 5   five measured coordinative structures')
    print('=' * 96)
    print(f'{"system":<20}{"N":>5}{"H":>10}{"log2 N":>10}{"J":>9}{"2^H":>10}'
          f'{"hub %":>8}{"uniform %":>11}')
    for r in table_measured():
        print(f'{r["system"]:<20}{r["N"]:>5}{r["H"]:>10.4f}{r["ceiling"]:>10.4f}'
              f'{r["J"]:>9.4f}{r["D_eff"]:>10.2f}{r["hub_share"]:>8.1f}'
              f'{r["uniform_share"]:>11.1f}')
        if r['note']:
            print(f'    note: {r["note"]}')

    print()
    print('=' * 96)
    print('SUPPLEMENTARY NOTE 2   the leader scan and the cortical file')
    print('=' * 96)
    s = leader_scan()
    print('   Eq 15 over the number of leaders h, predicted minus measured control entropy, bits')
    for sport in ('Football', 'Handball', 'Basketball'):
        g = s[sport]
        tail = ''
        if sport == 'Handball':
            tail = f'   three leaders hold {s["handball_three"]:.3f}'
        elif sport == 'Basketball':
            tail = '   formed team by team, then averaged'
        print(f'   {sport:<12} h = 1 to {len(g):<3}' + ' '.join(f'{x:.3f}' for x in g) + tail)
    c = cortical_file()
    diag = (f'{c["diag_min"]:.5f} in every entry' if c["diag_max"] - c["diag_min"] < 5e-6
            else f'{c["diag_min"]:.5f} to {c["diag_max"]:.5f}')
    print(f'   cortical file: diagonal {diag}, '
          f'off-diagonal {c["off_min"]:+.2f} to {c["off_max"]:+.2f}, '
          f'{c["negative"]:.0f} per cent negative')
    print(f'   eigenvector and strength shares order the regions at a rank correlation '
          f'of {c["rho"]:.3f}')

    print()
    print('=' * 96)
    print('SUPPLEMENTARY TABLE 6   outcome entropy in the expertise data')
    print('=' * 96)
    t = table_expertise()
    print(f'   novice  mean RE {t["novice_mean"]:.2f} cm   SD {t["novice_sd"]:.3f} cm   '
          f'h {t["novice_h"]:.3f} bits')
    print(f'   expert  mean RE {t["expert_mean"]:.2f} cm   SD {t["expert_sd"]:.3f} cm   '
          f'h {t["expert_h"]:.3f} bits')
    print(f'   difference by Eq 12 {t["difference"]:.3f} bits, '
          f'predictability ratio {t["ratio"]:.3f}')
    print(f'   h is the mean information of the ten values under the Gaussian fitted by '
          f'maximum likelihood, {t["novice_h_column"]:.3f} and {t["expert_h_column"]:.3f} '
          f'bits from the released column')
    print(f'   two-sided F test on the ratio of the two variances, Supplementary Note 3: '
          f'F({t["df"][0]}, {t["df"][1]}) = {t["F"]:.2f}, p = {t["p"]:.3f}')
    print()
    print('   The participant values are those of Appendix B of the supplementary material')
    print('   of Park (2027), as Supplementary Note 3 records.')

    print()
    print('=' * 96)
    print('SUPPLEMENTARY NOTE 2 and SUPPLEMENTARY FIG. 1D   the cortical threshold sweep')
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

    print()
    print('=' * 96)
    print('SUPPLEMENTARY NOTE 2   the two released muscle columns')
    print('=' * 96)
    m = muscle_columns()
    print(f'   {m["muscles"]} muscles. Attachments, first bin read as degree two: '
          f'measured column {m["att_measured"]:,}, released null column {m["att_null"]:,}')
    print(f'   released null   1 - J {m["released"]:.4f}, measured / null {m["factor_released"]:.2f}')
    print(f'   nulls rebuilt at {m["att_measured"]} attachments, each muscle given two and the '
          f'other {m["extra"]} spread at random,')
    print(f'   {m["draws"]:,} draws: mean 1 - J {m["rebuilt"]:.4f}, measured / null '
          f'{m["factor_rebuilt"]:.2f}')


if __name__ == '__main__':
    main()
