import matplotlib.pyplot as plt
import numpy as np


##############################################
def _compute_perf_profile_T(H, gate):
    """
    Pure (non-plotting) computation of the performance-profile T matrix.

      H contains a three dimensional array of function values.
        H[f,p,s] = function value # f for problem p and solver s.
      gate is a positive constant reflecting the convergence tolerance.

    Returns the unsorted, pre-NaN-replacement T matrix (np-by-ns), where
    T[p, s] is the number of evaluations required for solver s to reach the
    cutoff value on problem p (NaN if it never does).
    """

    H = H.copy()  # Avoid mutating the caller's array

    nf, nprob, ns = H.shape  # Grab the dimensions

    # Produce a suitable history array with sorted entries:
    for j in range(ns):
        for i in range(1, nf):
            H[i, :, j] = np.minimum(H[i, :, j], H[i - 1, :, j])

    prob_min = np.nanmin(H, axis=(0, 2))  # The minimum value seen for each problem
    prob_max = H[0, :, 0]  # The starting value for each problem

    # For each problem and solver, determine the number of evaluations
    # required to reach the cutoff value
    T = np.zeros((nprob, ns))
    for p in range(nprob):
        cutoff = prob_min[p] + gate * (prob_max[p] - prob_min[p])
        for s in range(ns):
            nfevs = np.flatnonzero(H[:, p, s] <= cutoff)
            if nfevs.size == 0:
                T[p, s] = np.nan
            else:
                T[p, s] = nfevs[0] + 1

    return T


def _ratios_from_T(T):
    """Compute performance ratios: T divided by the smallest element in each row."""
    row_min = np.nanmin(T, axis=1, keepdims=True)
    return T / row_min


##############################################
def plot_perf_profile(H, gate, logplot=False, legendstr=None):
    """
    This subroutine produces a performance profile as described in:

    Benchmarking Derivative-Free Optimization Algorithms
    Jorge J. More' and Stefan M. Wild
    SIAM J. Optimization, Vol. 20 (1), pp.172-191, 2009.

    Performance profiles were originally introduced in
    Benchmarking optimization software with performance profiles,
    E.D. Dolan and J.J. More',
    Mathematical Programming, 91 (2002), 201--213.

    The subroutine returns a handle to lines in a performance profile.

      H contains a three dimensional array of function values.
        H[f,p,s] = function value # f for problem p and solver s.
      gate is a positive constant reflecting the convergence tolerance.
      logplot=True is used to indicate that a log (base 2) plot is desired.

    Argonne National Laboratory
    Jorge More' and Stefan Wild. January 2008.
    """

    nf, nprob, ns = H.shape  # Grab the dimensions

    if legendstr is None:
        legendstr = [f"solver {s}" for s in range(ns)]

    T = _compute_perf_profile_T(H, gate)

    # Other colors, lines, and markers are easily possible:
    lines = ["-", "-.", "--"]
    colors = ["b", "r", "k", "m", "c", "g", "y"]
    markers = ["s", "o", "^", "v", "p", "<", "x", "h", "+", "d", "*", "<"]

    # Compute ratios and divide by the smallest element in each row.
    r = _ratios_from_T(T)

    # Replace all NaN's with twice the max_ratio and sort.
    max_ratio = np.nanmax(r)
    r[np.isnan(r)] = 2 * max_ratio
    r = np.sort(r, axis=0)

    # Plot stair graphs with markers.
    hl = [None] * ns
    for s in range(ns):
        col = r[:, s]
        xs = np.append(col, col[-1])
        ys = np.arange(1, nprob + 2) / nprob

        # Only plot one marker at the intercept
        if xs[0] == 1:
            vv = np.nonzero(xs == 1)[0][-1]
            xs = xs[vv:]
            ys = ys[vv:]

        sl = s % len(lines)
        sc = s % len(colors)
        sm = s % len(markers)
        option1 = f"{lines[sl]}{colors[sc]}{markers[sm]}"

        (hl[s],) = plt.step(xs, ys, option1, where="post", label=legendstr[s])

    # Axis properties are set so that failures are not shown, but with the
    # max_ratio data points shown. This highlights the "flatline" effect.
    if logplot:
        plt.xscale("log")
        twop = int(np.floor(np.log2(1.1 * max_ratio)))
        plt.xticks(2.0 ** np.arange(0, twop + 1))
    plt.axis([1, 1.1 * max_ratio, 0, 1])
    plt.xlabel("Performance Ratio")
    plt.ylabel("Proportion of Solved Problems")
    plt.legend()

    return hl
