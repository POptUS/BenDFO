import os
import sys

import numpy as np
import scipy as sp
from calfun import calfun
from dfoxs import dfoxs

# np.seterr("raise")


def doit():
    dfo = np.loadtxt("../data/dfo.dat")
    Results = {}

    for p, probtype in enumerate(["absnormal", "absuniform", "abswild", "noisy3", "nondiff", "relnormal", "reluniform", "relwild", "smooth", "wild3"]):
        for row, (nprob, n, m, factor_power) in enumerate(dfo[0:1]):
            n = int(n)
            m = int(m)

            X0 = dfoxs(n, nprob, int(10**factor_power)).T
            [y, F, G, J] = calfun(X0, m, int(nprob), probtype, gradout=True)

            Results[str(p) + "_" + str(row)] = {}
            Results[str(p) + "_" + str(row)]["y"] = y
            Results[str(p) + "_" + str(row)]["F"] = F
            Results[str(p) + "_" + str(row)]["G"] = G
            Results[str(p) + "_" + str(row)]["J"] = J


if __name__ == "__main__":
    doit()
