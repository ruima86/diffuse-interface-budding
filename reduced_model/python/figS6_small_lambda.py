"""Bending resistance and interfacial pull of the reduced model for decreasing Lambda~ (Fig. S6 of the
Supplemental Material). For Lambda~ = 0.15, 0.05, 0.02, chi~ is chosen to give the planar line tensions of
the four curves of Fig. 2(c) (Lambda~ = 0.15, chi~ = 15, 18, 21, 25).

Saves output/small_lambda.npz and output/figS6_small_lambda.png. Runtime: about ten minutes.
Uses LaTeX for the labels if available.
"""
import os, shutil
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import reduced_model as M
from line_tension import gamma_pl, chi_at

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, 'output'); os.makedirs(OUT, exist_ok=True)
a_grid = np.geomspace(3e-5, 0.4, 130)
lams = [0.15, 0.05, 0.02]
gams = [gamma_pl(0.15, c) for c in (15.0, 18.0, 21.0, 25.0)]
print("gamma_pl R0/kappa of the Fig. 2(c) curves:", [round(g, 3) for g in gams], flush=True)

data = {}
for lam in lams:
    for k, g in enumerate(gams):
        chi = chi_at(lam, g, lo=3.0, hi=600.0)
        L = M.landscape(lam, chi, a_grid); data[f'lam{lam}_g{k}'] = L; data[f'chi_lam{lam}_g{k}'] = chi
        mn, mx, b = M.extrema(L)
        print(f"Lambda~ = {lam:4.2f}  gamma_pl = {g:5.3f}  chi~ = {chi:7.2f}  min = {mn}  max = {mx}", flush=True)
np.savez(os.path.join(OUT, 'small_lambda.npz'), **data)

def make_figure(usetex):
    plt.rcParams.update({'text.usetex': usetex, 'font.family': 'serif', 'font.size': 11,
                         'axes.linewidth': 0.8, 'xtick.direction': 'in', 'ytick.direction': 'in',
                         'xtick.top': True, 'ytick.right': True})
    cols = ['#3b4cc0', '#1a8a99', '#8c6d1f', '#b2182b']
    fig, ax = plt.subplots(1, 3, figsize=(7.2, 2.75), sharey=True,
                           gridspec_kw=dict(wspace=0.07, left=0.095, right=0.995, bottom=0.17, top=0.90))
    for k, (axi, lam) in enumerate(zip(ax, lams)):
        for j, g in enumerate(gams):
            a, Pb, Pp = M.forces(data[f'lam{lam}_g{j}']); m = (a >= 1.5e-4) & (a <= 0.35)
            if j == 0:
                axi.plot(a[m], Pb[m], 'k', lw=1.6, label=r'$\Pi_b$')
            axi.plot(a[m], Pp[m], color=cols[j], lw=1.4, label=rf'$\Pi_p$, $\widetilde\gamma_{{\rm pl}}={g:.2f}$')
            d = (Pp - Pb)[m]; am = a[m]; Pbm = Pb[m]
            for i in np.where(np.diff(np.sign(d)) != 0)[0]:              # stationary states: v stable, ^ barrier top
                x0, x1 = np.log(am[i]), np.log(am[i+1]); xr = np.exp(x0 - d[i]*(x1 - x0)/(d[i+1] - d[i]))
                yr = np.interp(np.log(xr), np.log(am), Pbm)
                axi.plot(xr, yr, 'v' if d[i] < 0 else '^', color=cols[j], ms=6.5, mec='k', mew=.5, zorder=5)
        axi.axhline(2*np.pi*2*M.KAPPA*(1/np.sqrt(1 - M.PHI0) + 1/np.sqrt(M.PHI0)), color='0.5', ls='--', lw=0.9,
                    zorder=0, label=r'$2\pi\gamma_{\rm ideal}$')
        axi.set_xscale('log'); axi.set_xlim(1.5e-4, 0.35); axi.set_ylim(500, 1420); axi.set_xlabel(r'$\widetilde r_n$')
        axi.set_title(rf'({"abc"[k]})$\quad\widetilde\Lambda={lam}$', fontsize=11, loc='left')
    ax[0].set_ylabel(r'$\Pi_b,\ \Pi_p$')
    ax[2].legend(frameon=False, fontsize=8.5, loc='upper right', handlelength=1.6, borderaxespad=0.3)
    fig.savefig(os.path.join(OUT, 'figS6_small_lambda.png'), dpi=300)


try:
    make_figure(shutil.which('latex') is not None)
except RuntimeError:                       # incomplete LaTeX installation
    plt.close('all'); make_figure(False)
print("saved output/figS6_small_lambda.png")
