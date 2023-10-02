# This is a python implementation of dfovec.m,
# provided at https://github.com/POptUS/BenDFO
import numpy as np
from constants import set_constants

def dfovec(m, n, x, nprob):

    [c13, c14, c29, c45, v, y1, y2, y3, y4, y5] = set_constants()

    # Initialize things
    fvec = np.zeros(m)
    total = 0

    if nprob == 1:  # Linear function - full rank.
        for j in range(n):
            total = total + x[j]
        temp = 2 * total / m + 1
        for i in range(m):
            fvec[i] = -temp
            if i < n:
                fvec[i] = fvec[i] + x[i]
    elif nprob == 2:  # Linear function - rank 1.
        for j in range(n):
            total = total + (j + 1) * x[j]
        for i in range(m):
            fvec[i] = (i + 1) * total - 1
    elif nprob == 3:  # Linear function - rank 1 with zero columns and rows.
        for j in range(1, n - 1):
            total = total + (j + 1) * x[j]
        for i in range(m - 1):
            fvec[i] = i * total - 1
        fvec[m - 1] = -1
    elif nprob == 4:  # Rosenbrock function.
        fvec[0] = 10 * (x[1] - x[0] * x[0])
        fvec[1] = 1 - x[0]
    elif nprob == 5:  # Helical valley function.
        if x[0] > 0:
            th = np.arctan(x[1] / x[0]) / (2 * np.pi)
        elif x[0] < 0:
            th = np.arctan(x[1] / x[0]) / (2 * np.pi) + 0.5
        elif x[0] == x[1] and x[1] == 0:
            th = 0.0
        else:
            th = 0.25
        r = np.sqrt(x[0] * x[0] + x[1] * x[1])
        fvec[0] = 10 * (x[2] - 10 * th)
        fvec[1] = 10 * (r - 1)
        fvec[2] = x[2]
    elif nprob == 6:  # Powell singular function.
        fvec[0] = x[0] + 10 * x[1]
        fvec[1] = np.sqrt(5) * (x[2] - x[3])
        fvec[2] = (x[1] - 2 * x[2]) ** 2
        fvec[3] = np.sqrt(10) * (x[0] - x[3]) ** 2
    elif nprob == 7:  # Freudenstein and Roth function.
        fvec[0] = -c13 + x[0] + ((5 - x[1]) * x[1] - 2) * x[1]
        fvec[1] = -c29 + x[0] + ((1 + x[1]) * x[1] - c14) * x[1]
    elif nprob == 8:  # Bard function.
        for i in range(15):
            tmp1 = i + 1
            tmp2 = 15 - i
            tmp3 = tmp1
            if i > 7:
                tmp3 = tmp2
            fvec[i] = y1[i] - (x[0] + tmp1 / (x[1] * tmp2 + x[2] * tmp3))
    elif nprob == 9:  # Kowalik and Osborne function.
        for i in range(11):
            tmp1 = v[i] * (v[i] + x[1])
            tmp2 = v[i] * (v[i] + x[2]) + x[3]
            fvec[i] = y2[i] - x[0] * tmp1 / tmp2
    elif nprob == 10:  # Meyer function.
        for i in range(16):
            temp = 5 * (i + 1) + c45 + x[2]
            tmp1 = x[1] / temp
            tmp2 = np.exp(tmp1)
            fvec[i] = x[0] * tmp2 - y3[i]
    elif nprob == 11:  # Watson function.
        for i in range(29):
            div = (i + 1) / c29
            s1 = 0
            dx = 1
            for j in range(1, n):
                s1 = s1 + j * dx * x[j]
                dx = div * dx
            s2 = 0
            dx = 1
            for j in range(n):
                s2 = s2 + dx * x[j]
                dx = div * dx
            fvec[i] = s1 - s2 * s2 - 1
        fvec[29] = x[0]
        fvec[30] = x[1] - x[0] * x[0] - 1
    elif nprob == 12:  # Box 3-dimensional function.
        for i in range(m):
            temp = i + 1
            tmp1 = temp / 10
            fvec[i] = np.exp(-tmp1 * x[0]) - np.exp(-tmp1 * x[1]) + (np.exp(-temp) - np.exp(-tmp1)) * x[2]
    elif nprob == 13:  # Jennrich and Sampson function.
        for i in range(m):
            temp = i + 1
            fvec[i] = 2 + 2 * temp - np.exp(temp * x[0]) - np.exp(temp * x[1])
    elif nprob == 14:  # Brown and Dennis function.
        for i in range(m):
            temp = (i + 1) / 5
            tmp1 = x[0] + temp * x[1] - np.exp(temp)
            tmp2 = x[2] + np.sin(temp) * x[3] - np.cos(temp)
            fvec[i] = tmp1 * tmp1 + tmp2 * tmp2
    elif nprob == 15:  # Chebyquad function.
        for j in range(n):
            t1 = 1
            t2 = 2 * x[j] - 1
            t = 2 * t2
            for i in range(m):
                fvec[i] = fvec[i] + t2
                th = t * t2 - t1
                t1 = t2
                t2 = th
        iev = -1
        for i in range(m):
            fvec[i] = fvec[i] / n
            if iev > 0:
                fvec[i] = fvec[i] + 1 / ((i + 1) ** 2 - 1)
            iev = -iev
    elif nprob == 16:  # Brown almost-linear function.
        total1 = -(n + 1)
        prod1 = 1
        for j in range(n):
            total1 = total1 + x[j]
            prod1 = x[j] * prod1
        for i in range(n - 1):
            fvec[i] = x[i] + total1
        fvec[n - 1] = prod1 - 1
    elif nprob == 17:  # Osborne 1 function.
        for i in range(33):
            temp = 10 * i
            tmp1 = np.exp(-x[3] * temp)
            tmp2 = np.exp(-x[4] * temp)
            fvec[i] = y4[i] - (x[0] + x[1] * tmp1 + x[2] * tmp2)
    elif nprob == 18:  # Osborne 2 function.
        for i in range(65):
            temp = i / 10
            tmp1 = np.exp(-x[4] * temp)
            tmp2 = np.exp(-x[5] * (temp - x[8]) ** 2)
            tmp3 = np.exp(-x[6] * (temp - x[9]) ** 2)
            tmp4 = np.exp(-x[7] * (temp - x[10]) ** 2)
            fvec[i] = y5[i] - (x[0] * tmp1 + x[1] * tmp2 + x[2] * tmp3 + x[3] * tmp4)  # noqa
    elif nprob == 19:  # Bdqrtic
        # n >= 5, m = (n-4)*2
        for i in range(n - 4):
            fvec[i] = -4 * x[i] + 3.0
            fvec[n - 4 + i] = x[i] ** 2 + 2 * x[i + 1] ** 2 + 3 * x[i + 2] ** 2 + 4 * x[i + 3] ** 2 + 5 * x[n - 1] ** 2
    elif nprob == 20:  # Cube
        # n = 2, m = n
        fvec[0] = x[0] - 1.0
        for i in range(1, n):
            fvec[i] = 10 * (x[i] - x[i - 1] ** 3)
    elif nprob == 21:  # Mancino
        # n = 2, m = n
        for i in range(n):
            ss = 0
            for j in range(n):
                v2 = np.sqrt(x[i] ** 2 + (i + 1) / (j + 1))
                ss = ss + v2 * ((np.sin(np.log(v2))) ** 5 + (np.cos(np.log(v2))) ** 5)  # noqa
            fvec[i] = 1400 * x[i] + (i - 49) ** 3 + ss
    elif nprob == 22:  # Heart8ls
        # m = n = 8
        fvec[0] = x[0] + x[1] + 0.69
        fvec[1] = x[2] + x[3] + 0.044
        fvec[2] = x[4] * x[0] + x[5] * x[1] - x[6] * x[2] - x[7] * x[3] + 1.57
        fvec[3] = x[6] * x[0] + x[7] * x[1] + x[4] * x[2] + x[5] * x[3] + 1.31
        fvec[4] = x[0] * (x[4] ** 2 - x[6] ** 2) - 2.0 * x[2] * x[4] * x[6] + x[1] * (x[5] ** 2 - x[7] ** 2) - 2.0 * x[3] * x[5] * x[7] + 2.65
        fvec[5] = x[2] * (x[4] ** 2 - x[6] ** 2) + 2.0 * x[0] * x[4] * x[6] + x[3] * (x[5] ** 2 - x[7] ** 2) + 2.0 * x[1] * x[5] * x[7] - 2.0
        fvec[6] = (
            x[0] * x[4] * (x[4] ** 2 - 3.0 * x[6] ** 2) + x[2] * x[6] * (x[6] ** 2 - 3.0 * x[4] ** 2) + x[1] * x[5] * (x[5] ** 2 - 3.0 * x[7] ** 2) + x[3] * x[7] * (x[7] ** 2 - 3.0 * x[5] ** 2) + 12.6
        )
        fvec[7] = (
            x[2] * x[4] * (x[4] ** 2 - 3.0 * x[6] ** 2) - x[0] * x[6] * (x[6] ** 2 - 3.0 * x[4] ** 2) + x[3] * x[5] * (x[5] ** 2 - 3.0 * x[7] ** 2) - x[1] * x[7] * (x[7] ** 2 - 3.0 * x[5] ** 2) - 9.48
        )
    else:
        print(f"unrecognized function number {nprob}")
        return None
    return fvec
