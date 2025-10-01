import numpy as np
import jax
import jax.numpy as jnp

from dfovec import dfovec
from dfovec_jax import dfovec_jax
from jacobian import jacobian
from dfoxs import dfoxs  # for generating X0s

# Load benchmark problem definitions
dfo_table = np.loadtxt("../data/dfo.dat")


def compare_jacobians(nprob, m, n, x0):
    x0 = np.array(x0)
    fvec_np = dfovec(m, n, x0, nprob)

    # Get FD-based Jacobian
    J_fd, _ = jacobian(m, n, x0, nprob)  # shape (m, n)

    # Get JAX-based Jacobian
    def f(x):
        return dfovec_jax(m, n, x, nprob)

    J_jax = jax.jacfwd(f)(jnp.array(x0))  # shape (m, n)

    # Compare
    diff = np.linalg.norm(J_fd - np.array(J_jax))
    rel_diff = diff / (np.linalg.norm(J_fd) + 1e-12)

    print(f"nprob={nprob:2d}, m={m:3d}, n={n:2d} | abs diff: {diff:.2e}, rel diff: {rel_diff:.2e}")
    if diff > 6e-7 and rel_diff > 6e-7:
        raise Exception("Difference is too large")


# Loop over all benchmark problems and 3 starting points
for nprob, n, m, factor_power in dfo_table:
    n = int(n)
    m = int(m)
    nprob = int(nprob)
    scale = int(10**factor_power)

    for pt in range(3):
        if pt == 0:
            x0 = dfoxs(n, nprob, scale)
        elif pt == 1:
            x0 = 0.1 * np.ones(n)
        elif pt == 2:
            x0 = 0.1 * np.arange(1, n + 1)

        try:
            compare_jacobians(nprob, m, n, x0)
        except Exception as e:
            print(f"nprob={nprob:2d} failed at pt={pt}: {e}")
