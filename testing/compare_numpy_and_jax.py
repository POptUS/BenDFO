import sys
from pathlib import Path

import jax
import jax.numpy as jnp
import numpy as np
jax.config.update("jax_enable_x64", True)

sys.path.append("../py")
from calfun import calfun
from dfovec_jax import dfovec_jax
from dfoxs import dfoxs

# Load problem definitions: nprob, n, m, scale_power
dfo_table = np.loadtxt("../data/dfo.dat")  # Adjust path as needed


def compare_outputs(nprob, m, n, x0):
    x0 = np.array(x0)

    # --- Reference values from calfun ---
    fval_ref, fvec_ref, grad_ref, J_ref = calfun(x0, m, nprob, probtype="smooth", num_outs=4)
    J_ref = J_ref.T  # Transpose back: calfun gives (n x m), we want (m x n)

    # --- JAX values ---
    x_jax = jnp.array(x0)

    def fvec_fun(x):
        return dfovec_jax(m, n, x, nprob)

    def scalar_fun(x):
        return jnp.sum(fvec_fun(x) ** 2)

    fvec_jax = fvec_fun(x_jax)
    fval_jax = scalar_fun(x_jax)
    grad_jax = jax.grad(scalar_fun)(x_jax)
    J_jax = jax.jacfwd(fvec_fun)(x_jax)

    # --- Differences ---
    fval_diff = np.abs(fval_ref - float(fval_jax))
    fvec_diff = np.linalg.norm(fvec_ref - np.array(fvec_jax))
    grad_diff = np.linalg.norm(grad_ref - np.array(grad_jax))
    jac_diff = np.linalg.norm(J_ref - np.array(J_jax))

    fval_rel = fval_diff / (np.abs(fval_ref) + 1e-14)
    fvec_rel = fvec_diff / (np.linalg.norm(fvec_ref) + 1e-14)
    grad_rel = grad_diff / (np.linalg.norm(grad_ref) + 1e-14)
    jac_rel = jac_diff / (np.linalg.norm(J_ref) + 1e-14)

    # --- Print Summary ---
    print(
        f"nprob={nprob:2d}, m={m:3d}, n={n:2d} | "
        f"fval Δ={fval_diff:.2e} ({fval_rel:.2e}), "
        f"fvec Δ={fvec_diff:.2e} ({fvec_rel:.2e}), "
        f"grad Δ={grad_diff:.2e} ({grad_rel:.2e}), "
        f"J Δ={jac_diff:.2e} ({jac_rel:.2e})"
    )

    # --- Fail if any relative diff is large ---
    if any(val > 6e-7 for val in [fval_rel, fvec_rel, grad_rel, jac_rel]):
        raise Exception("Large differences detected")


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
            compare_outputs(nprob, m, n, x0)
        except Exception as e:
            print(f"nprob={nprob:2d}, pt={pt} failed: {e}")
