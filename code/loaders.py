"""Readers for the released source files, and the provenance record.

Every file the analysis touches is listed here with its size and hash, so a
reader can confirm they hold what the paper was computed from before running
anything.

Two of the inputs belong to other groups and are not redistributed here. Run
`data/external/fetch_external.py` to obtain them; the hashes below are checked on
arrival.
"""
import hashlib
import io
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
DERIVED = ROOT / 'data' / 'derived'
EXTERNAL = ROOT / 'data' / 'external'

# name -> (bytes, sha256 prefix). Recorded in Materials and Methods for the two
# external files, so these are the values the paper itself commits to.
PROVENANCE = {
    'team_sport_position_centrality.csv': (6075, '950150735ce4fddd'),
    'team_sport_coupling_shares.csv': (1093, '2d501f8b79bcc26d'),
    'cortical_coupling_shares.csv': (4617, 'ad11ced52a2fbb42'),
    'cortical_threshold_sweep.csv': (610, 'b19446c8c9ccd523'),
    'cov_MOTOR_matrix.csv': (158979, 'bded170235212706'),
    'muscle_degree_counts.csv': (499, 'f331e4fd5fcebd7f'),
    'musculoskeletal_coupling_shares.csv': (8119, 'f263f68a371ab8a7'),
    'musculoskeletal_random_null_shares.csv': (5440, 'bf040ad84443ac00'),
    'expertise_outcome_entropy.csv': (906, '6125d71889dcf310'),
}

EXTERNAL_PROVENANCE = {
    'cov_MOTOR': (80162, 'e00c39db9c61fada'),
    'fig2e.csv': (135, '516bfd49deeeae4e'),
}


def sha16(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()[:16]


def _strip_comments(text):
    return '\n'.join(ln for ln in text.splitlines()
                     if not ln.lstrip().lstrip('"').startswith('#'))


def read_csv(name, folder=None):
    """A released CSV, with its comment header stripped."""
    path = (Path(folder) if folder else DERIVED) / name
    return pd.read_csv(io.StringIO(
        _strip_comments(path.read_text(encoding='utf-8-sig'))))


def read_matrix(name='cov_MOTOR_matrix.csv', folder=None):
    """The cortical covariance matrix as a plain array."""
    path = (Path(folder) if folder else DERIVED) / name
    text = _strip_comments(path.read_text(encoding='utf-8-sig'))
    return np.array([[float(v) for v in ln.split(',')]
                     for ln in text.splitlines() if ln.strip()])


def check_provenance(verbose=True):
    """Confirm every released file is the one the paper was computed from."""
    ok = True
    for name, (size, digest) in PROVENANCE.items():
        path = DERIVED / name
        if not path.exists():
            ok = False
            if verbose:
                print(f'  MISSING  {name}')
            continue
        got_size, got_hash = path.stat().st_size, sha16(path)
        good = got_size == size and got_hash == digest
        ok &= good
        if verbose:
            print(f'  {"OK " if good else ">>>"} {name:<42} '
                  f'{got_size:>8} bytes  sha256 {got_hash}')
    for name, (size, digest) in EXTERNAL_PROVENANCE.items():
        path = EXTERNAL / name
        if not path.exists():
            if verbose:
                print(f'  --  {name:<42} not fetched, skipped')
            continue
        got_size, got_hash = path.stat().st_size, sha16(path)
        good = got_size == size and got_hash == digest
        ok &= good
        if verbose:
            print(f'  {"OK " if good else ">>>"} {name:<42} '
                  f'{got_size:>8} bytes  sha256 {got_hash}')
    return ok
