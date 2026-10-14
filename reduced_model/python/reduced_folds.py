"""Folds of the reduced model (L_fold^re in Fig. 3; reduced-model columns of the fold table in SM Sec. S9;
the fold at Lambda~ = 0.10 quoted in the main text).

For each Lambda~, chi~ is bisected between a value with an open-neck minimum and one without; the fold is
where the minimum and the barrier top merge. Reports chi~_c, the neck radius at the fold and gamma_pl R0/kappa
there, and the barrier at a few chi~ below the fold. Runtime: tens of minutes (one process per Lambda~ is
started if the argument --parallel is given).
"""
import sys
import numpy as np
from multiprocessing import Pool
import reduced_model as M
from line_tension import gamma_pl

CASES = {0.10: (27.0, 37.0, np.geomspace(5e-5, 0.30, 150)),
         0.15: (20.0, 22.0, np.geomspace(2e-4, 0.30, 130)),
         0.30: (11.0, 13.0, np.geomspace(2e-4, 0.30, 130)),
         0.45: (8.6, 10.2, np.geomspace(2e-4, 0.30, 130))}
NBISECT = 7


def fold(lam):
    lo, hi, grid = CASES[lam]
    best = None; log = []
    for _ in range(NBISECT):
        mid = 0.5*(lo + hi)
        mn, mx, b = M.extrema(M.landscape(lam, mid, grid))
        log.append(f"  Lambda~ = {lam:.2f}  chi~ = {mid:7.4f}  min = {mn}  max = {mx}  barrier = {b}")
        if mn is not None and mx is not None:
            lo = mid; best = (mid, mn, mx)
        else:
            hi = mid
    chi_c = 0.5*(lo + hi)
    r = np.sqrt(best[1]*best[2]) if best else np.nan
    return lam, chi_c, lo, hi, r, gamma_pl(lam, chi_c), log


if __name__ == '__main__':
    lams = list(CASES)
    res = Pool(len(lams)).map(fold, lams) if '--parallel' in sys.argv else [fold(l) for l in lams]
    for lam, chi_c, lo, hi, r, g, log in res:
        print("\n".join(log))
    print(f"\n{'Lambda~':>8} {'chi~_c':>8} {'bracket':>18} {'r_n at fold':>12} {'gamma_pl R0/kappa':>18}")
    for lam, chi_c, lo, hi, r, g, log in res:
        print(f"{lam:8.2f} {chi_c:8.3f}   ({lo:7.3f},{hi:7.3f}) {r:12.4f} {g:18.3f}")
