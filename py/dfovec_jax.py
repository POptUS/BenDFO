import jax.numpy as jnp


def set_constants():
    c13 = 1.3e1
    c14 = 1.4e1
    c29 = 2.9e1
    c45 = 4.5e1
    v = jnp.array([4.0, 2.0, 1.0, 0.5, 0.25, 0.167, 0.125, 0.1, 0.0833, 0.0714, 0.0625])
    y1 = jnp.array([0.14, 0.18, 0.22, 0.25, 0.29, 0.32, 0.35, 0.39, 0.37, 0.58, 0.73, 0.96, 1.34, 2.1, 4.39])
    y2 = jnp.array([0.1957, 0.1947, 0.1735, 0.16, 0.0844, 0.0627, 0.0456, 0.0342, 0.0323, 0.0235, 0.0246])
    y3 = jnp.array([34780.0, 28610.0, 23650.0, 19630.0, 16370.0, 13720.0, 11540.0, 9744.0, 8261.0, 7030.0, 6005.0, 5147.0, 4427.0, 3820.0, 3307.0, 2872.0])
    y4 = jnp.array(
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
    y5 = jnp.array(
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
    fvec = jnp.zeros(m)
    total = 0

    if nprob == 1:  # Linear function - full rank.
        total = jnp.sum(x[:n])
        temp = 2 * total / m + 1
        fvec = -temp * jnp.ones(m)
        fvec = fvec.at[:n].add(x[:n])

    elif nprob == 2:  # Linear function - rank 1.
        weights = jnp.arange(1, n + 1)
        total = jnp.sum(weights * x)
        fvec = (jnp.arange(1, m + 1) * total) - 1

    elif nprob == 3:  # Linear function - rank 1 with zero columns and rows.
        weights = jnp.arange(2, n)
        total = jnp.sum(weights * x[1:-1])
        fvec = jnp.zeros(m).at[: m - 1].set(jnp.arange(m - 1) * total - 1).at[m - 1].set(-1.0)

    elif nprob == 4:  # Rosenbrock function.
        fvec = jnp.array([10 * (x[1] - x[0] ** 2), 1 - x[0]])

    elif nprob == 5:  # Helical valley function.

        def theta(x0, x1):
            angle = jnp.arctan2(x1, x0) / (2 * jnp.pi)
            return jnp.where(x0 > 0, angle, jnp.where(x0 < 0, angle + 0.5, jnp.where((x0 == 0) & (x1 == 0), 0.0, 0.25)))

        th = theta(x[0], x[1])
        r = jnp.sqrt(x[0] ** 2 + x[1] ** 2)
        fvec = jnp.array([10 * (x[2] - 10 * th), 10 * (r - 1), x[2]])

    elif nprob == 6:  # Powell singular function.
        fvec = jnp.array([x[0] + 10 * x[1], jnp.sqrt(5.0) * (x[2] - x[3]), (x[1] - 2 * x[2]) ** 2, jnp.sqrt(10.0) * (x[0] - x[3]) ** 2])

    elif nprob == 7:  # Freudenstein and Roth function.
        fvec = jnp.array([-c13 + x[0] + ((5 - x[1]) * x[1] - 2) * x[1], -c29 + x[0] + ((1 + x[1]) * x[1] - c14) * x[1]])

    elif nprob == 8:  # Bard function.
        i = jnp.arange(15)
        tmp1 = i + 1
        tmp2 = 15 - i
        tmp3 = jnp.where(i > 7, tmp2, tmp1)
        denom = x[1] * tmp2 + x[2] * tmp3
        fvec = y1 - (x[0] + tmp1 / denom)

    elif nprob == 9:  # Kowalik and Osborne function.
        tmp1 = v * (v + x[1])
        tmp2 = v * (v + x[2]) + x[3]
        fvec = y2 - x[0] * tmp1 / tmp2

    elif nprob == 10:  # Meyer function.
        i = jnp.arange(16)
        temp = 5 * (i + 1) + c45 + x[2]
        tmp1 = x[1] / temp
        tmp2 = jnp.exp(tmp1)
        fvec = x[0] * tmp2 - y3

    elif nprob == 11:  # Watson function.

        def watson_term(i):
            div = (i + 1) / c29
            dx = div ** jnp.arange(n)
            dx1 = jnp.arange(1, n) * div ** jnp.arange(1, n)
            s1 = jnp.sum(dx1 * x[1:])
            s2 = jnp.sum(dx * x)
            return s1 - s2**2 - 1

        fvec = jnp.array([watson_term(i) for i in range(29)] + [x[0], x[1] - x[0] ** 2 - 1])

    elif nprob == 12:  # Box 3-dimensional function.
        i = jnp.arange(1, m + 1)
        tmp1 = i / 10.0
        term = jnp.exp(-tmp1[:, None] * x[:2])
        const = jnp.exp(-i) - jnp.exp(-tmp1)
        fvec = term[:, 0] - term[:, 1] + const * x[2]

    elif nprob == 13:  # Jennrich and Sampson function.
        i = jnp.arange(1, m + 1)
        fvec = 2 + 2 * i - jnp.exp(i * x[0]) - jnp.exp(i * x[1])

    elif nprob == 14:  # Brown and Dennis function.
        i = jnp.arange(1, m + 1)
        temp = i / 5.0
        tmp1 = x[0] + temp * x[1] - jnp.exp(temp)
        tmp2 = x[2] + jnp.sin(temp) * x[3] - jnp.cos(temp)
        fvec = tmp1**2 + tmp2**2

    elif nprob == 15:  # Chebyquad function.
        t = 2 * x - 1
        T = jnp.polynomial.chebyshev.chebvander(t, m - 1).T
        coeffs = jnp.mean(T, axis=1)
        correction = jnp.array([0.0 if i % 2 == 0 else 1 / ((i + 1) ** 2 - 1) for i in range(m)])
        fvec = coeffs + correction

    elif nprob == 16:  # Brown almost-linear function.
        total = jnp.sum(x) - (n + 1)
        prod = jnp.prod(x)
        fvec = jnp.concatenate([x[:-1] + total, jnp.array([prod - 1])])

    elif nprob == 17:  # Osborne 1 function.
        i = jnp.arange(33)
        tmp1 = jnp.exp(-x[3] * 10 * i)
        tmp2 = jnp.exp(-x[4] * 10 * i)
        fvec = y4 - (x[0] + x[1] * tmp1 + x[2] * tmp2)

    elif nprob == 18:  # Osborne 2 function.
        i = jnp.arange(65)
        t = i / 10.0
        tmp1 = jnp.exp(-x[4] * t)
        tmp2 = jnp.exp(-x[5] * (t - x[8]) ** 2)
        tmp3 = jnp.exp(-x[6] * (t - x[9]) ** 2)
        tmp4 = jnp.exp(-x[7] * (t - x[10]) ** 2)
        fvec = y5 - (x[0] * tmp1 + x[1] * tmp2 + x[2] * tmp3 + x[3] * tmp4)

    elif nprob == 19:  # Bdqrtic
        f1 = -4 * x[: n - 4] + 3
        f2 = sum((i + 1) * x[i + j] ** 2 for j, i in enumerate(range(n - 4)))
        fvec = jnp.concatenate([f1, f2.reshape(-1)])

    elif nprob == 20:  # Cube
        fvec = jnp.zeros(n)
        fvec = fvec.at[0].set(x[0] - 1.0)
        fvec = fvec.at[1:].set(10 * (x[1:] - x[:-1] ** 3))

    elif nprob == 21:  # Mancino

        def mancino_term(i):
            j = jnp.arange(n)
            v2 = jnp.sqrt(x[i] ** 2 + (i + 1) / (j + 1))
            return 1400 * x[i] + (i - 49) ** 3 + jnp.sum(v2 * (jnp.sin(jnp.log(v2)) ** 5 + jnp.cos(jnp.log(v2)) ** 5))

        fvec = jnp.array([mancino_term(i) for i in range(n)])

    elif nprob == 22:  # Heart8ls
        fvec = jnp.zeros(8)
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
