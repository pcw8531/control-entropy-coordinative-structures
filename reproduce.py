"""Run everything, in order, and report.

    python reproduce.py            the gate plus every table
    python reproduce.py --gate     the gate alone, about a minute
"""
import hashlib
import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
CODE = ROOT / 'code'
CORE_MARK = '# CORE  -  Equations 1 to 13'


def core_cell_hashes():
    """The Methods states that the equation module is byte-identical across the
    analysis notebooks and verified by hash. This checks it.

    Three notebooks do not carry the module. fig4_bottom_layer, fig5_publication
    and fig_S3_checks define the measures they need locally, so that each runs
    standalone. Those local definitions are checked against the module below
    rather than assumed to agree.
    """
    def digest(text):
        return hashlib.sha256(text[text.index(CORE_MARK):].rstrip().encode()
                              ).hexdigest()[:16]

    carries, standalone = {}, []
    carries['core.py'] = digest((CODE / 'core.py').read_text())
    for nb_path in sorted((CODE / 'notebooks').glob('*.ipynb')):
        nb = json.loads(nb_path.read_text())
        found = False
        for c in nb['cells']:
            src = ''.join(c['source'])
            if CORE_MARK in src:
                carries[nb_path.name] = digest(src)
                found = True
                break
        if not found:
            standalone.append(nb_path.name)
    return carries, standalone


def local_definitions_agree():
    """The standalone notebooks define share, entropy, evenness and Gini inline.
    Check those definitions return what the module returns."""
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
    carries, standalone = core_cell_hashes()
    for name, digest in carries.items():
        print(f'  {name:<44} {digest}')
    identical = len(set(carries.values())) == 1
    print(f'  identical across {len(carries)} files: {identical}')
    if not identical:
        print('  >>> the module has drifted between files. Fix before trusting anything.')
    if standalone:
        print()
        print('  These notebooks define their measures locally rather than carrying')
        print('  the module, so that each runs standalone:')
        for name in standalone:
            print(f'    {name}')
        worst = local_definitions_agree()
        print(f'  largest disagreement between the local definitions and the '
              f'module, over 500 random distributions: {worst:.2e}')

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
