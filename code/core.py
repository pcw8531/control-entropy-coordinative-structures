"""Control entropy of coordinative structures: the single equation module.

This file is byte-identical, line for line, to the CORE cell carried in every
analysis notebook, so every number in the paper traces to one definition.
Verified by hash; see reproduce.py.
"""
# ============================================================================
# CORE  -  Equations 1 to 13 of the PNAS manuscript (Methods numbering).
# Byte-identical cell in every notebook, verified by hash (roadmap step 6.4).
# ============================================================================
# Eq 1   p_i = c_i / sum_j c_j                       coupling share
# Eq 2   H = - sum_i p_i log2 p_i                     control entropy (avg information)
# Eq 3   p_i = C_i / sum_j C_j,  A C = lambda C       eigenvector weighting
# Eq 4   J = H / log2 N                               evenness
# Eq 5   D_eff = 2^H                                  effective count of DOF
# Eq 6   D_eff / N = 2^(-dH),  dH = log2 N - H        retained fraction / compression
# Eq 7   H <= log2 N                                  the bound (Gibbs / Jensen)
# Eq 8   dH/deps = log2(p_j / p_i)                    transfer gradient, concentration lowers H
# Eq 9   H = H2(phi) + phi log2 h + (1-phi) log2(N-h) leader-follower closed form
# Eq 10  kappa = <k^2> / <k>                          degree heterogeneity
# Eq 11  r                                            degree assortativity
# Eq 12  dh = log2(sigma_A / sigma_B)                 outcome entropy difference
# Eq 13 (Hebbian) lives in s4, since it needs the parent model's activations.
#
# I(p) = log2(1/p) is the summand of Eq 2 (Introduction Eq 2).
# Compression dH = log2 N - H is the bits saved (Introduction Eq 5), used by Eq 6.
# ============================================================================
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt

CORE_VERSION = "v5 / 24 July 2026 (Methods Eq 1-13)"

def share(c):
    """Eq 1. Normalise any non-negative vector into coupling shares."""
    c = np.asarray(c, dtype=float).ravel()
    if np.any(c < 0):
        raise ValueError("coupling values must be non-negative")
    t = c.sum()
    if t <= 0:
        raise ValueError("total coupling is zero")
    return c / t

def information(p):
    """Information of a share, in bits, the summand of Eq 2. Small share -> large information."""
    p = np.asarray(p, dtype=float)
    out = np.full(p.shape, np.inf)
    np.log2(1.0 / p, out=out, where=(p > 0))
    return out

def shannon_entropy(c):
    """Eq 2. H as the average of the information over the distribution. Input may be raw."""
    p = share(c)
    p = p[p > 0]
    return float((p * information(p)).sum())

def max_entropy(n):
    """Eq 7 (right side). The ceiling log2 N, attained only by the uniform distribution."""
    if n < 1:
        raise ValueError("n must be at least 1")
    return float(np.log2(n))

def evenness(h, n):
    """Eq 4. J = H / log2 N. Compare systems of different size with this."""
    m = max_entropy(n)
    return 1.0 if m == 0 else float(h / m)

def effective_dof(h):
    """Eq 5. D_eff = 2^H, a number of equally coupled degrees of freedom."""
    return float(2.0 ** h)

def compression(H, N):
    """dH = log2 N - H (Introduction Eq 5). Bits saved against even coupling.
    Regular gives 0 exactly; the scale-free ensemble gives 0.326 at N = 100.
    """
    return float(np.log2(N) - H)

def retained_fraction(H, N):
    """Eq 6. 2^(-dH) = D_eff / N, the fraction of degrees of freedom still to be
    specified. Equals 2^H / N. Scale-free gives 0.798, the 79.8 out of 100.
    """
    return float(2.0 ** (-(np.log2(N) - H)))

def coupling_shares(G, weight="degree"):
    """Eq 1 with c = degree, or Eq 3 with c = eigenvector centrality."""
    nodes = list(G.nodes())
    if weight == "degree":
        v = np.array([G.degree(n) for n in nodes], dtype=float)
    elif weight == "eigenvector":
        c = nx.eigenvector_centrality_numpy(G)
        v = np.abs(np.array([c[n] for n in nodes], dtype=float))
    else:
        raise ValueError("weight must be 'degree' or 'eigenvector'")
    return share(v)

def control_entropy(G, weight="degree"):
    """Eq 2 applied to Eq 1 or Eq 3."""
    return shannon_entropy(coupling_shares(G, weight=weight))

def strength_shares(W):
    """Eq 1 for a weighted structure: c_i = node strength = sum_j W_ij."""
    W = np.asarray(W, dtype=float)
    if W.ndim != 2 or W.shape[0] != W.shape[1]:
        raise ValueError("W must be square")
    return share(W.sum(axis=1))

def transfer_gradient(p, i, j):
    """Eq 8. dH/deps = log2(p_j / p_i) when eps moves coupling FROM j TO i.
    Negative whenever p_i > p_j, so any transfer toward an already-strong
    element strictly lowers H. i is the element that receives coupling.
    """
    p = np.asarray(p, dtype=float)
    return float(np.log2(p[j] / p[i]))

def _binary_entropy(phi):
    """H2(phi) = -phi log2 phi - (1-phi) log2(1-phi), with the 0 log 0 = 0 rule."""
    phi = float(phi)
    if phi <= 0.0 or phi >= 1.0:
        return 0.0
    return float(-phi * np.log2(phi) - (1.0 - phi) * np.log2(1.0 - phi))

def two_level_entropy(N, h, phi):
    """Eq 9. H = H2(phi) + phi log2 h + (1-phi) log2(N-h).
    h leaders share phi equally and N-h followers share (1-phi) equally.
    Returns log2 N exactly when phi = h/N, the built-in unit test.
    """
    if not (1 <= h < N):
        raise ValueError("need 1 <= h < N")
    if not (0.0 <= phi <= 1.0):
        raise ValueError("phi must be in [0, 1]")
    lead = phi * np.log2(h) if h > 0 else 0.0
    return float(_binary_entropy(phi) + lead + (1.0 - phi) * np.log2(N - h))

def degree_heterogeneity(G):
    """Eq 10. kappa = <k^2> / <k>, the parent paper's normalisation.

    NOTE ON NORMALISATION. The handwritten note of 23 July writes the ratio as
    <k^2>/<k>^2. That differs from this by a factor <k>, and the parent paper's
    published values (6.84, 6.09, 9.67) are on THIS definition. Table 1 of the
    manuscript uses this one so the two papers stay comparable. Convert with
    kappa_alt = kappa / <k> if a source uses the other.
    """
    k = np.array([d for _, d in G.degree()], dtype=float)
    if k.mean() <= 0:
        return float("nan")
    return float((k ** 2).mean() / k.mean())

def degree_assortativity(G):
    """Eq 11. Pearson correlation of degrees at the two ends of an edge.

    Undefined for a regular graph, where the degrees carry no variance. Returns
    nan there rather than a fabricated zero.
    """
    k = np.array([d for _, d in G.degree()], dtype=float)
    if k.std() == 0 or G.number_of_edges() == 0:
        return float("nan")
    try:
        r = nx.degree_assortativity_coefficient(G)
    except (ValueError, ZeroDivisionError):
        return float("nan")
    return float(r) if np.isfinite(r) else float("nan")

def gini(x):
    """Gini coefficient, matching the parent paper's definition."""
    x = np.sort(np.asarray(x, dtype=float))
    n = len(x)
    if x.sum() == 0:
        return 0.0
    idx = np.arange(1, n + 1)
    return float((2.0 * (idx * x).sum()) / (n * x.sum()) - (n + 1.0) / n)

GAUSS_CONST = 0.5 * np.log2(2.0 * np.pi * np.e)

def outcome_entropy_gaussian(sigma):
    """h = 0.5 log2(2 pi e sigma^2). A DIFFERENTIAL entropy, not Eq 2.

    It shares the unit and the reading of predictability with control entropy
    and is a different quantity on a different object. Only differences on a
    common scale are interpreted, which is why Eq 12 exists.
    """
    s = np.asarray(sigma, dtype=float)
    if np.any(s <= 0):
        raise ValueError("sigma must be positive")
    return GAUSS_CONST + np.log2(s)

def delta_outcome_entropy(sigma_a, sigma_b):
    """Eq 12. The additive constant and the unit both cancel."""
    return float(np.log2(float(sigma_a) / float(sigma_b)))

def predictability_ratio(dh):
    """R = 2^dh, the predictability ratio that follows Eq 12."""
    return float(2.0 ** float(dh))

print("core loaded:", CORE_VERSION)

