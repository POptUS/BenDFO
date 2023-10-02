import os
import sys

import numpy as np
import scipy as sp
from calfun import calfun
from dfoxs import dfoxs

# np.seterr("raise")


dfo = np.loadtxt("../data/dfo.dat")
probtypes = ["absnormal", "absuniform", "abswild", "noisy3", "nondiff", "relnormal", "reluniform", "relwild", "smooth", "wild3"]

Results = {}

for p, probtype in enumerate(probtypes):
    for row, (nprob, n, m, factor_power) in enumerate(dfo):
        n = int(n)
        m = int(m)

        for pt in range(3):
            if pt == 0:
                X0 = dfoxs(n, nprob, int(10**factor_power))
            elif pt == 1:
                X0 = 0.1 * np.ones(n)
            elif pt == 2:
                X0 = 0.1 * np.arange(1, n + 1)

            [y, F, G, J] = calfun(X0, m, int(nprob), probtype, num_outs = 4)

            Results["prob_" + str(p + 1) + "_" + str(row + 1) + "_" + str(pt + 1)] = {}
            Results["prob_" + str(p + 1) + "_" + str(row + 1) + "_" + str(pt + 1)]["X0"] = X0
            Results["prob_" + str(p + 1) + "_" + str(row + 1) + "_" + str(pt + 1)]["y"] = y
            Results["prob_" + str(p + 1) + "_" + str(row + 1) + "_" + str(pt + 1)]["F"] = F
            Results["prob_" + str(p + 1) + "_" + str(row + 1) + "_" + str(pt + 1)]["G"] = G
            Results["prob_" + str(p + 1) + "_" + str(row + 1) + "_" + str(pt + 1)]["J"] = J

sp.io.savemat("fvec_and_gradients_at_starting_values_python.mat", Results)
