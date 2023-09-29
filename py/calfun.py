# This is a python implementation of calfun.m,
# provided at https://github.com/POptUS/BenDFO
import numpy as np
from dfovec import dfovec
from jacobian import jacobian


def norm(x, type=2):
    if type == 1:
        return np.sum(np.abs(x))
    elif type == 2:
        return np.sqrt(x**2)
    else:  # type==np.inf:
        return max(np.abs(x))


def calfun(x, m, nprob, probtype="smooth", noise_level=1e-3, vecout=False, gradout=False):
    # Note: The vecout=True outputs are independent of probtype and noise_level

    n = len(x)

    # Restrict domain for some nondiff problems
    xc = x
    if probtype == "nondiff":
        if nprob == 8 or nprob == 9 or nprob == 13 or nprob == 16 or nprob == 17 or nprob == 18:
            xc = np.maximum(x, 0)

    # Generate the vector
    fvec = dfovec(m, n, xc, nprob)

    if vecout:
        return fvec

    # Calculate the function value
    if probtype == "absnormal":
        sigma = noise_level
        z = sigma * np.random.randn(m)
        fvec = fvec + z
        y = np.sum(fvec**2)
    elif probtype == "absuniform":
        sigma = noise_level
        z = (sigma * np.sqrt(3)) * (2 * np.random.rand(m) - 1)
        fvec = fvec + z
        y = np.sum(fvec**2)
    elif probtype == "abswild":
        z = 0.9 * np.sin(100 * np.linalg.norm(x, 1)) * np.cos(100 * np.linalg.norm(x, np.inf)) + 0.1 * np.cos(np.linalg.norm(x, 2))
        z = z * (4 * z**2 - 3)
        y = np.sum(fvec**2) + z
    elif probtype == "relnormal":
        sigma = noise_level
        z = sigma * np.random.randn(m)
        fvec = fvec * (1 + z)
        y = np.sum(fvec**2)
    elif probtype == "reluniform":
        sigma = noise_level
        z = (sigma * np.sqrt(3)) * (2 * np.random.rand(m) - 1)
        fvec = fvec * (1 + z)
        y = np.sum(fvec**2)
    elif probtype == "relwild":
        sigma = noise_level
        z = 0.9 * np.sin(100 * np.linalg.norm(x, 1)) * np.cos(100 * np.linalg.norm(x, np.inf)) + 0.1 * np.cos(np.linalg.norm(x, 2))
        z = z * (4 * z**2 - 3)
        y = (1 + sigma * z) * np.sum(fvec**2)
    elif probtype == "noisy3":
        sigma = 10**-3
        u = sigma * (-1 + 2 * np.random.rand(m))
        fvec = fvec * (1 + u)
        y = np.sum(fvec**2)
    elif probtype == "wild3":
        sigma = 10**-3
        phi = 0.9 * np.sin(100 * np.linalg.norm(x, 1)) * np.cos(100 * np.linalg.norm(x, np.inf)) + 0.1 * np.cos(np.linalg.norm(x, 2))
        phi = phi * (4 * phi**2 - 3)
        y = (1 + sigma * phi) * np.sum(fvec**2)
    elif probtype == "smooth":
        y = np.sum(fvec**2)
    elif probtype == "nondiff":
        y = np.sum(np.abs(fvec))
    else:
        print(f"invalid probtype {probtype}")
        return None

    # Never return nan. Return inf instead so that
    # optimization algorithms treat it as out of bounds.
    if np.isnan(y):
        return np.inf

    if gradout:
        G = jacobian(m, n, xc, nprob)

        if probtype == "nondiff":
            G = J * np.sign(fvec)
        elif probtype in ["relnormal", "reluniform", "noisy3"]:
            G = (1 + noise_level**2) * J * np.sign(fvec)
        else:
            G = 2 * J * fvec
        return y, fvec, G, J
    else:
        return y
