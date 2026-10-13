"""The reproduction gate. Every value the paper prints, recomputed from the
released files and from the generators the Methods describes, and compared
against what is printed.

Run this before anything else. It exits non-zero if a value fails to reproduce.
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from core import (share, shannon_entropy, max_entropy, evenness, effective_dof,
                  compression, gini, two_level_entropy, delta_outcome_entropy,
                  predictability_ratio)
from loaders import check_provenance, read_csv, read_matrix
from s2_measured import (sport_shares, basketball_profiles, cortical_shares,
                         musculoskeletal_shares, musculoskeletal_null_shares,
                         expertise, eigen_shares_matrix)

CHECKS = []


def chk(label, got, printed, tol=5e-4, unit=''):
    ok = abs(got - printed) <= tol
    CHECKS.append((label, got, printed, ok))
    print(f'  {"OK " if ok else ">>>"} {label:<50} {got:>12.4f}   '
          f'printed {printed:>10.4f} {unit}')


def comp(p):
    return 1 - evenness(shannon_entropy(p), len(p))


def main():
    print('=' * 92)
    print('PROVENANCE')
    print('=' * 92)
    provenance_ok = check_provenance()

    print()
    print('=' * 92)
    print('MEASURED STRUCTURES, from the released files')
    print('=' * 92)
    p_fb, p_hb = sport_shares('Football'), sport_shares('Handball')
    chk('football H', shannon_entropy(p_fb), 3.3875, unit='bits')
    chk('football evenness', evenness(shannon_entropy(p_fb), 11), 0.9792)
    chk('football effective count', effective_dof(shannon_entropy(p_fb)), 10.47, tol=6e-3)
    chk('football compression', comp(p_fb), 0.0208)
    chk('football Gini', gini(p_fb), 0.1701)
    chk('handball H', shannon_entropy(p_hb), 2.2829, unit='bits')
    chk('handball evenness', evenness(shannon_entropy(p_hb), 6), 0.8832)
    chk('handball effective count', effective_dof(shannon_entropy(p_hb)), 4.87, tol=6e-3)
    chk('handball compression', comp(p_hb), 0.1168)
    chk('handball Gini', gini(p_hb), 0.3303)

    teams, prof = basketball_profiles()
    per = np.array([shannon_entropy(share(r)) for r in prof])
    chk('basketball H, mean of sixteen teams', per.mean(), 2.2314, unit='bits')
    chk('basketball SD across teams', per.std(ddof=1), 0.0477, unit='bits')
    chk('basketball H of the averaged profile',
        shannon_entropy(share(prof.mean(axis=0))), 2.2620, unit='bits')
    chk('basketball evenness', evenness(per.mean(), 5), 0.9610)
    ct = np.array([comp(share(r)) for r in prof])
    chk('basketball compression', ct.mean(), 0.0390)
    chk('basketball per-team lowest', ct.min(), 0.0083)
    chk('basketball per-team highest', ct.max(), 0.0971)

    p_cx = cortical_shares()
    chk('cortical H from the released matrix', shannon_entropy(p_cx), 6.5468, unit='bits')
    chk('cortical compression', comp(p_cx), 0.0146)
    chk('cortical effective count', effective_dof(shannon_entropy(p_cx)), 93.49, tol=6e-3)
    chk('cortical Gini', gini(p_cx), 0.2093)
    chk('cortical H from the shares file',
        shannon_entropy(read_csv('cortical_coupling_shares.csv')
                        .share_p_i_eigenvector.values), 6.5468, unit='bits')

    p_ms, p_mn = musculoskeletal_shares(), musculoskeletal_null_shares()
    chk('musculoskeletal H', shannon_entropy(p_ms), 7.7369, unit='bits')
    chk('musculoskeletal null H', shannon_entropy(p_mn), 7.9876, unit='bits')
    chk('musculoskeletal effective count', effective_dof(shannon_entropy(p_ms)),
        213.33, tol=6e-2)
    chk('musculoskeletal compression', comp(p_ms), 0.0421)
    chk('musculoskeletal null compression', comp(p_mn), 0.0111)
    chk('ratio to its own null', comp(p_ms) / comp(p_mn), 3.81, tol=6e-3)
    chk('musculoskeletal Gini', gini(p_ms), 0.3417)
    chk('musculoskeletal null Gini', gini(p_mn), 0.1940)
    md = read_csv('muscle_degree_counts.csv')
    chk('muscle count', float(md.muscles_measured.sum()), 270.0, tol=0)
    chk('attachments', float((md.bones_linked * md.muscles_measured).sum()), 977.0, tol=0)
    chk('attachments per bone', float((md.bones_linked * md.muscles_measured).sum() / 173),
        5.65, tol=6e-3)

    print()
    print('=' * 92)
    print('EXPERTISE')
    print('=' * 92)
    nov, exp = expertise()
    n = nov.outcome_entropy_h_bits.values
    e = exp.outcome_entropy_h_bits.values
    chk('novice outcome entropy', n.mean(), 3.316, unit='bits')
    chk('expert outcome entropy', e.mean(), 2.021, unit='bits')
    chk('difference of per-participant means', n.mean() - e.mean(), 1.295, unit='bits')
    chk('difference from the group mean SDs',
        delta_outcome_entropy(nov.outcome_sd_cm.mean(), exp.outcome_sd_cm.mean()),
        1.293, unit='bits')
    chk('predictability ratio', predictability_ratio(n.mean() - e.mean()), 2.454, tol=1e-3)
    chk('non-overlap', n.min() - e.max(), 0.836, unit='bits')
    from scipy import stats
    t, _ = stats.ttest_ind(n, e)
    chk('t statistic', t, 20.55, tol=6e-3)
    sp = np.sqrt((n.var(ddof=1) + e.var(ddof=1)) / 2)
    chk('Cohen d', (n.mean() - e.mean()) / sp, 9.19, tol=6e-3)

    print()
    print('=' * 92)
    print('THE LEADER-FOLLOWER CLOSED FORM, Eq 9')
    print('=' * 92)
    fb = np.sort(p_fb)[::-1]
    hb = np.sort(p_hb)[::-1]
    chk('football, one leader', two_level_entropy(11, 1, fb[0]), 3.4265, unit='bits')
    chk('handball, one leader', two_level_entropy(6, 1, hb[0]), 2.5195, unit='bits')
    chk('handball, three leaders, combined share', hb[:3].sum(), 0.806, tol=6e-4)
    chk('handball, three leaders', two_level_entropy(6, 3, hb[:3].sum()), 2.2946, unit='bits')
    pt = np.array([share(r) for r in prof])
    chk('basketball, one leader, averaged over teams',
        float(np.mean([two_level_entropy(5, 1, np.sort(p)[::-1][0]) for p in pt])),
        2.2580, unit='bits')
    chk('volleyball, uniform by rule', max_entropy(6), 2.5850, unit='bits')

    print()
    print('=' * 92)
    print('THE SIMULATED ENSEMBLES, Table 1')
    print('=' * 92)
    from s1_ensembles import table_1
    printed = {'Regular lattice': (6.6439, 1.0000, 100.0, 0.000, 6.00, 6.6439),
               'Small-world (WS)': (6.5832, 0.9909, 95.9, 0.061, 6.09, 6.6324),
               'Random (ER)': (6.4820, 0.9756, 89.4, 0.162, 6.90, 6.5313),
               'Scale-free (BA)': (6.3176, 0.9509, 79.8, 0.326, 10.10, 6.2865)}
    for r in table_1():
        pv = printed[r['topology']]
        chk(f'{r["topology"]} H', r['H'], pv[0], unit='bits')
        chk(f'{r["topology"]} evenness', r['J'], pv[1])
        chk(f'{r["topology"]} effective count', r['D_eff'], pv[2], tol=6e-2)
        chk(f'{r["topology"]} compression', r['dH'], pv[3], tol=6e-4, unit='bits')
        chk(f'{r["topology"]} heterogeneity', r['kappa'], pv[4], tol=6e-3)
        chk(f'{r["topology"]} H degree-weighted', r['H_degree'], pv[5], unit='bits')
        chk(f'{r["topology"]} links per realisation', r['links'], 300.0, tol=0)

    print()
    print('=' * 92)
    bad = [c for c in CHECKS if not c[3]]
    print(f'{len(CHECKS) - len(bad)} of {len(CHECKS)} values reproduce.')
    if not provenance_ok:
        print('PROVENANCE FAILED: a released file is not the one the paper used.')
    for label, got, printed_v, _ in bad:
        print(f'   MISMATCH  {label}: computed {got:.6f}, printed {printed_v}')
    return 0 if (not bad and provenance_ok) else 1


if __name__ == '__main__':
    sys.exit(main())
