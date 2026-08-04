"""Fetch the one external file that is not redistributed here, and check it.

cov_MOTOR comes from the task_covs release of Weaver et al. (2026). That
repository declares no licence, so the file is not included in this deposit.
This script downloads it and confirms it is the one the paper was computed from:
80,162 bytes with sha256 beginning e00c39db9c61fada, which is what Materials and
Methods records.

fig2e.csv is already here. It is released under CC BY-SA 4.0 by Murphy et al.
(2018) at Zenodo doi:10.5281/zenodo.1069104 and is redistributed with attribution.
"""
import hashlib
import sys
import urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parent

TARGETS = {
    'cov_MOTOR': dict(
        url='https://raw.githubusercontent.com/NicholasJWeaver/'
            'BrainCompressibility2025/main/task_covs/cov_MOTOR',
        size=80162,
        sha256='e00c39db9c61fada',
        source='Weaver NJ, Faskowitz J, Betzel RF, Lynn CW. Quantifying the '
               'compressibility of the human brain. PNAS 2026;123(4):e2531115123.',
    ),
}


def sha16(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()[:16]


def main():
    failed = False
    for name, spec in TARGETS.items():
        dest = HERE / name
        if dest.exists():
            print(f'{name}: already here')
        else:
            print(f'{name}: downloading from {spec["url"]}')
            try:
                urllib.request.urlretrieve(spec['url'], dest)
            except Exception as exc:
                print(f'  download failed: {exc}')
                print(f'  fetch it by hand from the repository and place it at {dest}')
                failed = True
                continue
        size, digest = dest.stat().st_size, sha16(dest)
        ok = size == spec['size'] and digest == spec['sha256']
        print(f'  {size} bytes, sha256 {digest}   '
              f'{"matches the paper" if ok else "DOES NOT MATCH THE PAPER"}')
        print(f'  source: {spec["source"]}')
        failed |= not ok
    return 1 if failed else 0


if __name__ == '__main__':
    sys.exit(main())
