import matplotlib.pyplot as plt
import numpy as np


##############################################
def _compute_data_profile_T(HIST, N, gate, optimality_type="value"):
    """
    Pure (non-plotting) computation of the data-profile T matrix.

      HIST contains a three dimensional array of function values.
        H[f,p,s] = function value # f for problem p and solver s.
      N is an np-vector of (positive) budget units. If simplex
        gradients are desired, then N(p) would be n(p)+1, where n(p) is
        the number of variables for problem p.
      gate is a positive constant reflecting the convergence tolerance.
           meaning is interpreted by the optimality type
      optimality_type is a string indicating the measure of optimality to use
           "value"   -> use smallest function value seen by all solvers
           "absgrad" -> solved when absolute gradient is small, ||nabla f_k|| <= gate
           "relgrad" -> solved when relative gradient is small, ||nabla f_k|| <= gate*||nabla f_0||

    Returns the unsorted, pre-NaN-replacement T matrix (np-by-ns), where
    T[p, s] is the number of N-function bundles required for solver s to
    reach the cutoff value on problem p (NaN if it never does).
    """

    HIST = HIST.copy()  # Avoid mutating the caller's array

    nf, nprob, ns = HIST.shape  # Grab the dimensions

    # Produce a suitable history array with sorted entries:
    for j in range(ns):
        for i in range(1, nf):
            HIST[i, :, j] = np.minimum(HIST[i, :, j], HIST[i - 1, :, j])

    # build T with T[p,s] = number of iters to solve problem p with solver s
    match optimality_type:

        case "value":

            prob_min = np.nanmin(HIST, axis=(0, 2))  # The minimum value seen for each problem
            prob_max = HIST[0, :, 0]  # The starting value for each problem

            # For each problem and solver, determine the number of
            # N-function bundles (e.g.- gradients) required to reach the cutoff value
            T = np.zeros((nprob, ns))
            for p in range(nprob):
                cutoff = prob_min[p] + gate * (prob_max[p] - prob_min[p])
                for s in range(ns):
                    nfevs = np.argmax(HIST[:, p, s] <= cutoff) + 1  # use argmax to find first occurrence; +1 for zero index
                    if nfevs == 1 and not (HIST[0, p, s] <= cutoff):  # all HIST[:,p,s] values are above the cutoff; argmax returns first index
                        T[p, s] = np.nan
                    else:
                        T[p, s] = nfevs / N[p]

        case "absgrad":

            # For each problem and solver, determine the number of
            # N-function bundles (e.g.- gradients) required to reach the cutoff value
            T = np.zeros((nprob, ns))
            for p in range(nprob):
                cutoff = gate
                for s in range(ns):
                    nfevs = np.argmax(HIST[:, p, s] <= cutoff) + 1  # use argmax to find first occurrence; +1 for zero index
                    if nfevs == 1 and not (HIST[0, p, s] <= cutoff):
                        T[p, s] = np.nan
                    else:
                        T[p, s] = nfevs / N[p]

        case "relgrad":

            prob_max = HIST[0, :, 0]  # The starting grad norm for each problem

            # For each problem and solver, determine the number of
            # N-function bundles (e.g.- gradients) required to reach the cutoff value
            T = np.zeros((nprob, ns))
            for p in range(nprob):
                cutoff = gate * prob_max[p]
                for s in range(ns):
                    nfevs = np.argmax(HIST[:, p, s] <= cutoff) + 1  # use argmax to find first occurrence; +1 for zero index
                    if nfevs == 1 and not (HIST[0, p, s] <= cutoff):
                        T[p, s] = np.nan
                    else:
                        T[p, s] = nfevs / N[p]

    return T


##############################################
def plot_data_profile(HIST, N, gate, optimality_type="value", legendstr=None):
    """
    This subroutine produces a data profile as described in:

    Benchmarking Derivative-Free Optimization Algorithms
    Jorge J. More' and Stefan M. Wild
    SIAM J. Optimization, Vol. 20 (1), pp.172-191, 2009.

    The subroutine returns a handle to lines in a data profile.

      HIST contains a three dimensional array of function values.
        H[f,p,s] = function value # f for problem p and solver s.
      N is an np-vector of (positive) budget units. If simplex
        gradients are desired, then N(p) would be n(p)+1, where n(p) is
        the number of variables for problem p.
      gate is a positive constant reflecting the convergence tolerance.
           meaning is interpreted by the optimality type
      optimality_type is a string indicating the measure of optimality to use
           "value"   -> use smallest function value seen by all solvers
           "absgrad" -> solved when absolute gradient is small, ||nabla f_k|| <= gate
           "relgrad" -> solved when relative gradient is small, ||nabla f_k|| <= gate*||nabla f_0||

    Argonne National Laboratory
    Jorge More' and Stefan Wild. January 2008.
    """

    nf, nprob, ns = HIST.shape  # Grab the dimensions

    if legendstr is None:
        legendstr = [f"solver {s}" for s in range(ns)]

    T = _compute_data_profile_T(HIST, N, gate, optimality_type)

    ##############################################################
    # plot
    ##############################################################

    # Other colors, lines, and markers are easily possible:
    lines = ["-", "-.", "--"]
    colors = ["b", "r", "k", "m", "c", "g", "y"]
    markers = ["s", "o", "^", "v", "p", "<", "x", "h", "+", "d", "*", "<"]

    # Replace all NaN's with twice the max_ratio and sort.
    max_data = np.nanmax(T)
    T[np.isnan(T)] = 2 * max_data
    T = np.sort(T, axis=0)

    # For each solver, plot stair graphs with markers.
    hl = [None] * ns
    for s in range(ns):
        xs, ys = np.append(T[:, s], T[-1, s]), np.arange(1, nprob + 2) / nprob
        sl = s % len(lines)
        sc = s % len(colors)
        sm = s % len(markers)
        fstring = f"{lines[sl]}{colors[sc]}{markers[sm]}"

        (hl[s],) = plt.step(xs, ys, fstring, where="post", label=legendstr[s])

    plt.axis([0, 1.1 * max_data, 0, 1])
    plt.xlabel("Normalized Iterations")
    plt.ylabel("Proportion of Solved Problems")
    plt.legend()

    return hl
