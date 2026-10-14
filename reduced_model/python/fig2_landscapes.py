"""Constrained landscapes F_min(r_n) of the reduced model (Fig. 2(a) and 2(c) of the paper).

Lambda~ = 0.15 with chi~ = 15, 18, 21, 25 (the four curves of Fig. 2(a),(c)), and landscapes around the
reduced-model folds at Lambda~ = 0.3 and 0.45. Saves output/landscapes.npz; prints the stationary points.
Runtime: a few minutes.
"""
import os, time
import numpy as np
import reduced_model as M

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, 'output'); os.makedirs(OUT, exist_ok=True)
a_grid = np.geomspace(1e-4, 0.4, 110)
plan = {0.15: [15.0, 18.0, 21.0, 25.0], 0.3: [10.0, 11.0, 12.5], 0.45: [8.0, 8.6, 9.8]}

data = {}
for lam, chis in plan.items():
    for chi in chis:
        t = time.time()
        L = M.landscape(lam, chi, a_grid)
        data[f'lam{lam}_chi{chi}'] = L
        mn, mx, b = M.extrema(L)
        desc = (f"min r_n = {mn:.3e}, max r_n = {mx:.3e}, barrier = {b:.2f} kBT" if b is not None
                else "no open-neck minimum (neck closes)" if mn is None else f"min r_n = {mn:.3e}, no maximum")
        print(f"Lambda~ = {lam:4.2f}  chi~ = {chi:5.2f}:  {desc}   [{time.time()-t:.0f} s]", flush=True)
np.savez(os.path.join(OUT, 'landscapes.npz'), **data)
print("saved output/landscapes.npz  (rows: a, F, bend, grad, bulk, xi*, R_-, R_+)")
