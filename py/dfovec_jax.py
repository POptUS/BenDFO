import jax
import jax.numpy as np


def set_constants():
    c13 = 1.3e1
    c14 = 1.4e1
    c29 = 2.9e1
    c45 = 4.5e1
    v = np.array([4.0, 2.0, 1.0, 0.5, 0.25, 0.167, 0.125, 0.1, 0.0833, 0.0714, 0.0625])
    y1 = np.array([0.14, 0.18, 0.22, 0.25, 0.29, 0.32, 0.35, 0.39, 0.37, 0.58, 0.73, 0.96, 1.34, 2.1, 4.39])
    y2 = np.array([0.1957, 0.1947, 0.1735, 0.16, 0.0844, 0.0627, 0.0456, 0.0342, 0.0323, 0.0235, 0.0246])
    y3 = np.array([34780.0, 28610.0, 23650.0, 19630.0, 16370.0, 13720.0, 11540.0, 9744.0, 8261.0, 7030.0, 6005.0, 5147.0, 4427.0, 3820.0, 3307.0, 2872.0])
    y4 = np.array(
        [
            0.844,
            0.908,
            0.932,
            0.936,
            0.925,
            0.908,
            0.881,
            0.85,
            0.818,
            0.784,
            0.751,
            0.718,
            0.685,
            0.658,
            0.628,
            0.603,
            0.58,
            0.558,
            0.538,
            0.522,
            0.506,
            0.49,
            0.478,
            0.467,
            0.457,
            0.448,
            0.438,
            0.431,
            0.424,
            0.42,
            0.414,
            0.411,
            0.406,
        ]
    )
    y5 = np.array(
        [
            1.366,
            1.191,
            1.112,
            1.013,
            0.991,
            0.885,
            0.831,
            0.847,
            0.786,
            0.725,
            0.746,
            0.679,
            0.608,
            0.655,
            0.616,
            0.606,
            0.602,
            0.626,
            0.651,
            0.724,
            0.649,
            0.649,
            0.694,
            0.644,
            0.624,
            0.661,
            0.612,
            0.558,
            0.533,
            0.495,
            0.5,
            0.423,
            0.395,
            0.375,
            0.372,
            0.391,
            0.396,
            0.405,
            0.428,
            0.429,
            0.523,
            0.562,
            0.607,
            0.653,
            0.672,
            0.708,
            0.633,
            0.668,
            0.645,
            0.632,
            0.591,
            0.559,
            0.597,
            0.625,
            0.739,
            0.71,
            0.729,
            0.72,
            0.636,
            0.581,
            0.428,
            0.292,
            0.162,
            0.098,
            0.054,
        ]
    )
    return c13, c14, c29, c45, v, y1, y2, y3, y4, y5


def dfovec_jax(m, n, x, nprob):
    """
    This is a Python translation of the Matlab version of the subroutine dfovec.f
    This subroutine specifies the nonlinear benchmark problems in

    Benchmarking Derivative-Free Optimization Algorithms
    Jorge J. More' and Stefan M. Wild
    SIAM J. Optimization, Vol. 20 (1), pp.172-191, 2009.

    The latest version of this subroutine is always available at
          https://github.com/POptUS/BenDFO/
    The authors would appreciate feedback and experiences from numerical
    studies conducted using this subroutine.

    The data file dfo.dat defines suitable values of m and n
    for each problem number nprob.

    This subroutine defines the functions of 22 nonlinear
    least squares problems. The allowable values of (m,n) for
    functions 1,2 and 3 are variable but with m .ge. n.
    For functions 4,5,6,7,8,9 and 10 the values of (m,n) are
    (2,2),(3,3),(4,4),(2,2),(15,3),(11,4) and (16,3), respectively.
    Function 11 (Watson) has m = 31 with n usually 6 or 9.
    However, any n, n = 2,...,31, is permitted.
    Functions 12,13 and 14 have n = 3,2 and 4, respectively, but
    allow any m .ge. n, with the usual choices being 10,10 and 20.
    Function 15 (Chebyquad) allows m and n variable with m .ge. n.
    Function 16 (Brown) allows n variable with m = n.
    For functions 17 and 18, the values of (m,n) are
    (33,5) and (65,11), respectively.

    fvec = dfovec(m, n, x, nprob)
        fvec is an output array of length m which contains the nprob
          function evaluated at x.
        m and n are positive integer input variables. n must not
          exceed m.
        x is an input array of length n.
        nprob is a positive integer input variable which defines the
          number of the problem. nprob must not exceed 22.

    Argonne National Laboratory
    Jorge More' and Stefan Wild. January 2008.
    """

    c13, c14, c29, c45, v, y1, y2, y3, y4, y5 = set_constants()

    # Initialize things
    fvec = np.zeros(m)
    total = 0

    if nprob == 1:  # Linear function - full rank.
        total = np.sum(x[:n])
        temp = 2 * total / m + 1
        fvec = -temp * np.ones(m)
        fvec = fvec.at[:n].add(x[:n])

    elif nprob == 2:  # Linear function - rank 1.
        weights = np.arange(1, n + 1)
        total = np.sum(weights * x)
        fvec = (np.arange(1, m + 1) * total) - 1

    elif nprob == 3:  # Linear function - rank 1 with zero columns and rows.
        weights = np.arange(2, n)
        total = np.sum(weights * x[1:-1])
        fvec = np.zeros(m).at[: m - 1].set(np.arange(m - 1) * total - 1).at[m - 1].set(-1.0)

    elif nprob == 4:  # Rosenbrock function.
        fvec = np.array([10 * (x[1] - x[0] ** 2), 1 - x[0]])

    elif nprob == 5:  # Helical valley function.
        th = np.arctan2(x[1], x[0]) / (2.0 * np.pi)
        r = np.sqrt(x[0] ** 2 + x[1] ** 2)

        fvec = fvec.at[0].set(10.0 * (x[2] - 10.0 * th))
        fvec = fvec.at[1].set(10.0 * (r - 1.0))
        fvec = fvec.at[2].set(x[2])

    elif nprob == 6:  # Powell singular function.
        fvec = np.array([x[0] + 10 * x[1], np.sqrt(5.0) * (x[2] - x[3]), (x[1] - 2 * x[2]) ** 2, np.sqrt(10.0) * (x[0] - x[3]) ** 2])

    elif nprob == 7:  # Freudenstein and Roth function.
        fvec = np.array([-c13 + x[0] + ((5 - x[1]) * x[1] - 2) * x[1], -c29 + x[0] + ((1 + x[1]) * x[1] - c14) * x[1]])

    elif nprob == 8:  # Bard function.
        i = np.arange(15)
        tmp1 = i + 1
        tmp2 = 15 - i
        tmp3 = np.where(i > 7, tmp2, tmp1)
        denom = x[1] * tmp2 + x[2] * tmp3
        fvec = y1 - (x[0] + tmp1 / denom)

    elif nprob == 9:  # Kowalik and Osborne function.
        tmp1 = v * (v + x[1])
        tmp2 = v * (v + x[2]) + x[3]
        fvec = y2 - x[0] * tmp1 / tmp2

    elif nprob == 10:  # Meyer function.
        i = np.arange(16)
        temp = 5 * (i + 1) + c45 + x[2]
        tmp1 = x[1] / temp
        tmp2 = np.exp(tmp1)
        fvec = x[0] * tmp2 - y3

    elif nprob == 11:  # Watson function.
        fvec_body = []
        for i in range(29):
            div = (i + 1.0) / c29
            s1 = 0.0
            dx = 1.0
            for j in range(1, n):
                s1 += j * dx * x[j]
                dx *= div
            s2 = 0.0
            dx = 1.0
            for j in range(n):
                s2 += dx * x[j]
                dx *= div
            fvec_body.append(s1 - s2 * s2 - 1.0)
        fvec_tail = [x[0], x[1] - x[0] ** 2 - 1.0]
        fvec = np.array(fvec_body + fvec_tail)

    elif nprob == 12:  # Box 3-dimensional function.
        i = np.arange(1, m + 1)
        tmp1 = i / 10.0
        term = np.exp(-tmp1[:, None] * x[:2])
        const = np.exp(-i) - np.exp(-tmp1)
        fvec = term[:, 0] - term[:, 1] + const * x[2]

    elif nprob == 13:  # Jennrich and Sampson function.
        i = np.arange(1, m + 1)
        fvec = 2 + 2 * i - np.exp(i * x[0]) - np.exp(i * x[1])

    elif nprob == 14:  # Brown and Dennis function.
        i = np.arange(1, m + 1)
        temp = i / 5.0
        tmp1 = x[0] + temp * x[1] - np.exp(temp)
        tmp2 = x[2] + np.sin(temp) * x[3] - np.cos(temp)
        fvec = tmp1**2 + tmp2**2

    elif nprob == 15:  # Chebyquad function.
        for j in range(n):
            t1 = 1.0
            t2 = 2.0 * x[j] - 1.0
            t = 2.0 * t2
            for i in range(m):
                fvec = fvec.at[i].add(t2)
                th = t * t2 - t1
                t1 = t2
                t2 = th

        iev = -1
        for i in range(m):
            val = fvec[i] / n
            if iev > 0:
                val += 1.0 / ((i + 1) ** 2 - 1.0)
            fvec = fvec.at[i].set(val)
            iev = -iev

    elif nprob == 16:  # Brown almost-linear function.
        total = np.sum(x) - (n + 1)
        prod = np.prod(x)
        fvec = np.concatenate([x[:-1] + total, np.array([prod - 1])])

    elif nprob == 17:  # Osborne 1 function.
        i = np.arange(33)
        tmp1 = np.exp(-x[3] * 10 * i)
        tmp2 = np.exp(-x[4] * 10 * i)
        fvec = y4 - (x[0] + x[1] * tmp1 + x[2] * tmp2)

    elif nprob == 18:  # Osborne 2 function.
        i = np.arange(65)
        t = i / 10.0
        tmp1 = np.exp(-x[4] * t)
        tmp2 = np.exp(-x[5] * (t - x[8]) ** 2)
        tmp3 = np.exp(-x[6] * (t - x[9]) ** 2)
        tmp4 = np.exp(-x[7] * (t - x[10]) ** 2)
        fvec = y5 - (x[0] * tmp1 + x[1] * tmp2 + x[2] * tmp3 + x[3] * tmp4)

    elif nprob == 19:  # Bdqrtic
        for i in range(n - 4):
            fvec = fvec.at[i].set(-4.0 * x[i] + 3.0)
            quad = x[i] ** 2 + 2.0 * x[i + 1] ** 2 + 3.0 * x[i + 2] ** 2 + 4.0 * x[i + 3] ** 2 + 5.0 * x[n - 1] ** 2
            fvec = fvec.at[n - 4 + i].set(quad)

    elif nprob == 20:  # Cube
        fvec = np.zeros(n)
        fvec = fvec.at[0].set(x[0] - 1.0)
        fvec = fvec.at[1:].set(10 * (x[1:] - x[:-1] ** 3))

    elif nprob == 21:  # Mancino

        def mancino_term(i):
            j = np.arange(n)
            v2 = np.sqrt(x[i] ** 2 + (i + 1) / (j + 1))
            return 1400 * x[i] + (i - 49) ** 3 + np.sum(v2 * (np.sin(np.log(v2)) ** 5 + np.cos(np.log(v2)) ** 5))

        fvec = np.array([mancino_term(i) for i in range(n)])

    elif nprob == 22:  # Heart8ls
        fvec = np.zeros(8)
        fvec = fvec.at[0].set(x[0] + x[1] + 0.69)
        fvec = fvec.at[1].set(x[2] + x[3] + 0.044)
        fvec = fvec.at[2].set(x[4] * x[0] + x[5] * x[1] - x[6] * x[2] - x[7] * x[3] + 1.57)
        fvec = fvec.at[3].set(x[6] * x[0] + x[7] * x[1] + x[4] * x[2] + x[5] * x[3] + 1.31)
        fvec = fvec.at[4].set(x[0] * (x[4] ** 2 - x[6] ** 2) - 2 * x[2] * x[4] * x[6] + x[1] * (x[5] ** 2 - x[7] ** 2) - 2 * x[3] * x[5] * x[7] + 2.65)
        fvec = fvec.at[5].set(x[2] * (x[4] ** 2 - x[6] ** 2) + 2 * x[0] * x[4] * x[6] + x[3] * (x[5] ** 2 - x[7] ** 2) + 2 * x[1] * x[5] * x[7] - 2.0)
        fvec = fvec.at[6].set(
            x[0] * x[4] * (x[4] ** 2 - 3 * x[6] ** 2) + x[2] * x[6] * (x[6] ** 2 - 3 * x[4] ** 2) + x[1] * x[5] * (x[5] ** 2 - 3 * x[7] ** 2) + x[3] * x[7] * (x[7] ** 2 - 3 * x[5] ** 2) + 12.6
        )
        fvec = fvec.at[7].set(
            x[2] * x[4] * (x[4] ** 2 - 3 * x[6] ** 2) - x[0] * x[6] * (x[6] ** 2 - 3 * x[4] ** 2) + x[3] * x[5] * (x[5] ** 2 - 3 * x[7] ** 2) - x[1] * x[7] * (x[7] ** 2 - 3 * x[5] ** 2) - 9.48
        )

    else:
        raise NotImplementedError(f"nprob={nprob} not implemented")

    return fvec
