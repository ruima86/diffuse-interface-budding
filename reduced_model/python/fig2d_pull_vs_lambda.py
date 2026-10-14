"""Interfacial pull of the reduced model for decreasing Lambda~ at fixed planar line tension
gamma_pl R0/kappa = 6.5 (Fig. 2(d) of the paper and its inset).

Writes output/fig2d_pull_vs_Lambda.txt with columns
  Lambda~  chi~  r_n~  r_n~ (chi~/Lambda~)^(1/2)  Pi_p [kBT/R0]  Pi_p/Pi_p(r_n~ = 0.3).
Runtime: a few minutes.
"""
import os
import numpy as np
import reduced_model as M
from line_tension import chi_at

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, 'output'); os.makedirs(OUT, exist_ok=True)
a_grid = np.geomspace(1e-4, 0.4, 110)
with open(os.path.join(OUT, 'fig2d_pull_vs_Lambda.txt'), 'w') as f:
    f.write("# Reduced model: interfacial pull along gamma_pl R0/kappa = 6.5 (Fig. 2(d))\n"
            "# columns: Lambda~  chi~  r_n~  r_n~*sqrt(chi~/Lambda~)  Pi_p [kBT/R0]  Pi_p/Pi_p(r_n~=0.3)\n")
    for lam in [0.45, 0.30, 0.15, 0.10, 0.05]:
        chi = chi_at(lam, 6.5, lo=3.0, hi=200.0)
        L = M.landscape(lam, chi, a_grid)
        a, Pb, Pp = M.forces(L)
        P0 = np.interp(0.3, a, Pp)
        sel = (a > 1e-3) & (a < 0.35)
        i = np.where(sel)[0][np.argmin(Pp[sel])]
        print(f"Lambda~ = {lam:4.2f}  chi~ = {chi:7.3f}:  weak point at r_n = {a[i]:.4f} "
              f"= {a[i]/np.sqrt(lam/chi):.3f} (Lambda~/chi~)^1/2, depth {1 - Pp[i]/P0:.3f}", flush=True)
        for x, p in zip(a, Pp):
            if 2e-4 < x < 0.32:
                f.write(f"{lam:.2f} {chi:.3f} {x:.6e} {x*np.sqrt(chi/lam):.6e} {p:.4f} {p/P0:.5f}\n")
print("written output/fig2d_pull_vs_Lambda.txt")
