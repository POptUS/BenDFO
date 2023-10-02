import numpy as np
from g_dfovec_1d import g_dfovec_1d


def jacobian(m, n, x, nprob):
    """
    This subroutine computes the Jacobian of the nonlinear equations
    defining the benchmark problems in
    
    Benchmarking Derivative-Free Optimization Algorithms
    Jorge J. More' and Stefan M. Wild
    SIAM J. Optimization, Vol. 20 (1), pp.172-191, 2009.
    
    The dependencies of this function are based on a Python translation of
    the Matlab AD software "adimat" on a modified from of the dfovec function
    from
    http://www.mcs.anl.gov/~more/dfo/
    See the instructions of dfovec.m for additional details on these
    nonlinear benchmark problems (and appropriate values of m and n).
    
      J is an output array of size m-by-n, with J(i,j) denoting the
        derivative (evaluated at x) of the ith equation with respect to
        the jth variable.
      fvec returns the usual dfovec
      m and n are positive integer input variables. n must not
        exceed m.
      x is an input array of length n.
      nprob is a positive integer input variable which defines the
        number of the problem. nprob must not exceed 22.
    """

    t = 0
    g_t = 1

    J = np.zeros((m, n))

    for ind in range(n):
        g_fvec, fvec = g_dfovec_1d(g_t, t, ind, m, n, np.zeros(x.shape), x, nprob)
        J[:, ind] = g_fvec

    return J, fvec
