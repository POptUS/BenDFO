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


def calfun(x, m, nprob, probtype="smooth", noise_level=1e-3, num_outs=1):
    """
    This is a Python translation of a Matlab version of the subroutine calfun.f
    This subroutine returns a function value as used in:

    Benchmarking Derivative-Free Optimization Algorithms
    Jorge J. More' and Stefan M. Wild
    SIAM J. Optimization, Vol. 20 (1), pp.172-191, 2009.

    The latest version of this subroutine is always available at
          https://github.com/POptUS/BenDFO/
    The authors would appreciate feedback and experiences from numerical
    studies conducted using this subroutine.

    The subroutine returns the function value f(x)

      x is an input array of length n.
      y is an output that contains the function value at x.
      fvec is an m-by-1 array containing component function values at x.
      G is a 1-by-n array containing the gradient of the function f at x.
          For stochastic problem types, this is the gradient of the expectation of f.
          For deterministically noisy problem types, this ignore the noise.
          For the nondiff problem type, this is the subgradient J * sign(fvec).
      J is an n-by-m array containing the gradients of the component
          functions at x.
          J(i,j) contains the derivative of the jth equation wrt x(i).

    # Note: The vecout=True outputs are independent of probtype and noise_level
    If reproducibility is needed, the np.random seed should be set before
    this function is called.

    Additional problem descriptors are passed through the following fields:
      m is a positive integer (length of output from dfovec).
         n must not exceed m.
      nprob is a positive integer that defines the number of the problem.
         nprob must not exceed 22.
      probtype is a string specifying the type of problem desired:
          'smooth' corresponds to smooth (noise-free) problems
          'absnormal' corresponds to stochastic Gaussian absolute noise
          'absuniform' corresponds to stochastic uniform absolute noise
          'abswild' corresponds to deterministic absolute noise
          'relnormal' corresponds to stochastic Gaussian relative noise
          'reluniform' corresponds to stochastic uniform relative noise
          'relwild' corresponds to deterministic relative noise
          'nondiff' corresponds to piecewise-smooth problems'smooth' corresponds to smooth problems
          'wild3' corresponds to deterministic relative noise with
          'noisy3' corresponds to stochastically noisy problems
      noise_level is a standard deviation; it is ignored for 'smooth', 'nondiff',
         'noisy3', 'wild3', and 'abswild' problem types
      num_outs is the desired number of outputs
            1: returns y
            2: returns y, fvec
            3: returns y, fvec, G,
            4: returns y, fvec, G, J

    Argonne National Laboratory
    Jorge More' and Stefan Wild. January 2008.
    """

    n = len(x)

    # Restrict domain for some nondiff problems
    xc = x
    if probtype == "nondiff":
        if nprob == 8 or nprob == 9 or nprob == 13 or nprob == 16 or nprob == 17 or nprob == 18:
            xc = np.maximum(x, 0)

    # Generate the vector
    fvec = dfovec(m, n, xc, nprob)

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

    # The scipy.benchmark version of calfun doesn't return Nans.
    # if np.isnan(y):
    #     return np.inf

    if num_outs == 1:
        return y
    elif num_outs == 2:
        return y, fvec
    else:
        J, dummy = jacobian(m, n, xc, nprob)
        J = J.T
        # if probtype == "smooth":
        #     assert np.all(dummy == fvec), "Why do the fvecs from jacobian and dfovec disagree?"

        if probtype == "nondiff":
            G = J @ np.sign(fvec)
        elif probtype in ["relnormal", "reluniform", "noisy3"]:
            G = (1 + noise_level**2) * J @ np.sign(fvec)
        else:
            G = 2 * J @ fvec

        if num_outs == 3:
            return y, fvec, G
        elif num_outs == 4:
            return y, fvec, G, J
        else:
            raise ValueError("Unknown value for num_outs")
