"""Dependence on the area a_p of a phase-separating unit (SM Sec. S10), in the reduced model with R0 = 100 nm.

(1) Table of SM Sec. S10: a_p = 100 nm^2 (Lambda~ = 0.15, chi~ = 15) and a_p = 25 nm^2, with Lambda~ and chi~
    chosen so that the planar interface has the same width (22.2 nm) and line tension (4.32 pN).
(2) Scan of chi~ at a_p = 25 nm^2: line tension at which the open-neck minimum disappears, and the barrier
    versus the distance from it (to be compared with the a_p = 100 nm^2 barriers of fig2_landscapes.py).
(3) L_CB estimate of the reduced model: closed-neck energy, corrected for the residual gradient energy
    93 Lambda~ kBT of the prescribed profile (SM Sec. S5.1), compared with the open-neck minimum.
Runtime: several minutes.
"""
import numpy as np
import reduced_model as M
from line_tension import gamma_pl

R0 = 100.0
a_grid = np.geomspace(1e-4, 0.4, 110)


def set_ap(ap):
    M.INV_AP = R0**2/ap
    return M.INV_AP


print("(1) barrier at fixed interface width and line tension")
for ap, lam, chi in [(100.0, 0.15, 15.0), (25.0, 0.03932, 5.932)]:
    n = set_ap(ap)
    mn, mx, b = M.extrema(M.landscape(lam, chi, a_grid))
    top = f"max = {mx*R0:.3f} nm  barrier = {b:.2f} kBT" if b is not None else "barrier top below the grid"
    print(f"  a_p = {ap:5.1f} nm^2  eps = {M.KAPPA/n:.4f}  Lambda~ = {lam:.5f}  chi~ = {chi:.3f}  "
          f"gamma_pl R0/kappa = {gamma_pl(lam, chi, inv_ap=n):.3f}  open-neck minimum at r_n = {mn*R0:.2f} nm  {top}",
          flush=True)

print("(2) chi~ scan at a_p = 25 nm^2, Lambda~ = 0.03932")
n = set_ap(25.0)
for chi in [6.5, 7.0, 7.5, 8.0, 9.0]:
    L = M.landscape(0.03932, chi, a_grid); mn, mx, b = M.extrema(L)
    g = gamma_pl(0.03932, chi, inv_ap=n)
    if b is not None:
        print(f"  chi~ = {chi:4.1f}  gamma_pl R0/kappa = {g:6.3f}  min = {mn*R0:6.2f} nm  max = {mx*R0:6.3f} nm  "
              f"barrier = {b:5.2f} kBT", flush=True)
    else:
        print(f"  chi~ = {chi:4.1f}  gamma_pl R0/kappa = {g:6.3f}  "
              f"{'no open-neck minimum (neck closes)' if mn is None else 'minimum without barrier'}", flush=True)

print("(3) reduced-model L_CB at a_p = 100 nm^2, Lambda~ = 0.15")
n = set_ap(100.0); lam = 0.15; res = 93.0*lam
for chi in [10.0, 11.5, 13.0]:
    L = M.landscape(lam, chi, a_grid); a, F = L[:, 0], L[:, 1]; mn, mx, b = M.extrema(L)
    if mn is None:
        print(f"  chi~ = {chi:5.2f}  no open-neck minimum"); continue
    d = F[0] - np.interp(np.log(mn), np.log(a), F)
    print(f"  chi~ = {chi:5.2f}  F(closed) - F(min) = {d:7.2f} kBT, corrected by the residue {res:.2f}: {d - res:+7.2f} kBT",
          flush=True)
