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

        X0 = dfoxs(n, nprob, int(10**factor_power)).T
        [y, F, G, J] = calfun(X0, m, int(nprob), probtype, gradout=True)

        Results["prob_" + str(p + 1) + "_" + str(row + 1)] = {}
        Results["prob_" + str(p + 1) + "_" + str(row + 1)]["X0"] = X0
        Results["prob_" + str(p + 1) + "_" + str(row + 1)]["y"] = y
        Results["prob_" + str(p + 1) + "_" + str(row + 1)]["F"] = F
        Results["prob_" + str(p + 1) + "_" + str(row + 1)]["G"] = G
        Results["prob_" + str(p + 1) + "_" + str(row + 1)]["J"] = J

sp.io.savemat("fvec_and_gradients_at_starting_values_python.mat", Results)
