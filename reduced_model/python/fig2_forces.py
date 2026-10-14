"""Force decomposition of the reduced-model landscapes (Fig. 2(c) of the paper and the numbers quoted with it).

Reads output/landscapes.npz from fig2_landscapes.py. For each landscape:
  - checks Pi_p - Pi_b = dF_min/dr_n, with Pi_b = -d(bending)/dr_n and Pi_p = d(gradient + bulk)/dr_n;
  - locates the minimum of the pull Pi_p (the weak point), its position in units of (Lambda~/chi~)^(1/2)
    and the optimal width xi* there in units of r_n;
  - compares Pi_p at r_n = 0.3 with the line-tension force 2 pi gamma, gamma evaluated on the tanh profile.
Also prints the bending resistance as r_n -> 0, to be compared with 2 pi gamma_ideal.
Lengths in R0, energies in kBT, forces in kBT/R0.
"""
import os
import numpy as np
import reduced_model as M

_trapz = M._trapz
HERE = os.path.dirname(os.path.abspath(__file__))
D = np.load(os.path.join(HERE, 'output', 'landscapes.npz'))


def profile_line_tension(a, lam, xi, Rm, Rp):
    """gamma = 2 Lambda~ (R0^2/a_p) int (dphi/ds)^2 ds on the reduced-model profile, in kBT/R0."""
    (rm, Jm, Hm, phim, dpm), (rp, Jp, Hp, phip, dpp), _, _ = M.integrals(a, Rm, Rp, xi)
    return 2*lam*M.INV_AP*(_trapz(dpm**2*Jm, M.THETA) + _trapz(dpp**2*Jp, M.THETA))


print(f"{'Lam~':>5} {'chi~':>5} | {'split':>7} | {'r_dip':>8} {'r_dip/sqrt(L/chi)':>18} {'xi*/r_n':>8} | "
      f"{'Pi_p(0.3)':>9} {'2pi*gamma':>9} {'ratio':>6}")
for key in D.files:
    L = D[key]
    lam = float(key.split('_')[0][3:]); chi = float(key.split('_')[1][3:])
    a, Pb, Pp = M.forces(L)
    dF = np.gradient(L[:, 1], a)
    split = np.abs((Pp - Pb) - dF).max()/np.abs(dF).max()
    sel = (a > 1e-3) & (a < 0.35)                       # interior minimum, away from the grid ends
    i = np.where(sel)[0][np.argmin(Pp[sel])]
    j = np.argmin(abs(a - 0.3))
    g = profile_line_tension(a[j], lam, L[j, 5], L[j, 6], L[j, 7])
    print(f"{lam:5.2f} {chi:5.1f} | {split:7.1e} | {a[i]:8.4f} {a[i]/np.sqrt(lam/chi):18.3f} {L[i, 5]/a[i]:8.2f} | "
          f"{Pp[j]:9.1f} {2*np.pi*g:9.1f} {Pp[j]/(2*np.pi*g):6.3f}")

L = D['lam0.15_chi15.0']; a, Pb, Pp = M.forces(L)
R1, R2 = np.sqrt(1 - M.PHI0), np.sqrt(M.PHI0)
g_ideal = 2*M.KAPPA*(1/R1 + 1/R2)
print(f"\nPi_b at the smallest neck radius (r_n = {a[1]:.1e}): {Pb[1]:.1f} kBT/R0;  "
      f"2 pi gamma_ideal = {2*np.pi*g_ideal:.1f} kBT/R0  (gamma_ideal = 2 kappa (1/R1 + 1/R2) = {g_ideal:.2f} kBT/R0)")
