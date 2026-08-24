"""Run everything, in order, and report.

    python reproduce.py            the gate plus every table
    python reproduce.py --gate     the gate alone, about a minute
"""
import hashlib
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
CODE = ROOT / 'code'
CORE_MARK = '# CORE  -  Equations 1 to 13'


def module_digest():
    """One definition stands behind every number in the paper. This prints the
    hash of the equation module so a reader can confirm the copy they hold is
    the copy the values were computed from.

    The notebooks that draw the figures and the movie are not part of this
    deposit. They arrange values these scripts compute; nothing is calculated
    in them that is not calculated here.
    """
    def digest(text):
        return hashlib.sha256(text[text.index(CORE_MARK):].rstrip().encode()
                              ).hexdigest()[:16]

    return digest((CODE / 'core.py').read_text())


def local_definitions_agree():
    """An independent re-implementation of share, entropy, evenness and Gini,
    written inline here, must return what the module returns. This catches a
    module that has been edited without the values being recomputed."""
    import numpy as np
    sys.path.insert(0, str(CODE))
    from core import share, shannon_entropy, evenness, gini

    def local_share(x):
        x = np.asarray(x, float).ravel()
        return x / x.sum()

    def local_H(p):
        p = np.asarray(p, float)
        p = p[p > 0]
        return float(-(p * np.log2(p)).sum())

    def local_gini(x):
        x = np.sort(np.asarray(x, float))
        n = len(x)
        i = np.arange(1, n + 1)
        return float((2 * (i * x).sum()) / (n * x.sum()) - (n + 1) / n)

    rng = np.random.default_rng(0)
    worst = 0.0
    for _ in range(500):
        v = rng.gamma(1.4, 2.0, rng.integers(5, 300))
        p, q = share(v), local_share(v)
        worst = max(worst, float(np.abs(p - q).max()),
                    abs(shannon_entropy(p) - local_H(q)),
                    abs(1 - evenness(shannon_entropy(p), len(p))
                        - (1 - local_H(q) / np.log2(len(q)))),
                    abs(gini(p) - local_gini(q)))
    return worst


def run(script):
    print()
    print('#' * 100)
    print(f'# {script}')
    print('#' * 100)
    return subprocess.call([sys.executable, str(CODE / script)])


def main():
    print('=' * 100)
    print('THE EQUATION MODULE, one definition behind every number')
    print('=' * 100)
    print(f"  core.py                                      {module_digest()}")
    worst = local_definitions_agree()
    print(f'  largest disagreement between an independent inline definition and')
    print(f'  the module, over 500 random distributions: {worst:.2e}')

    code = run('verify_published_values.py')
    if code != 0:
        print('\nThe gate failed. Nothing further is run.')
        return code
    if '--gate' in sys.argv:
        return 0
    for s in ('s1_ensembles.py', 's2_measured.py', 's3_matched_nulls.py',
              's9_resampling.py'):
        code |= run(s)
    return code


if __name__ == '__main__':
    sys.exit(main())
