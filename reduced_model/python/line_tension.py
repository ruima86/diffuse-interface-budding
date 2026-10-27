"""Planar line tension gamma_pl R0/kappa of the flat interface between the coexisting phases
(Supplemental Material, Sec. S9.2):

    gamma_pl R0/kappa = (2/eps) Lambda~^(1/2) int_{phi_1}^{phi_2} sqrt(Delta g~) dphi,   eps = kappa a_p/(kBT R0^2),

with Delta g~ the Flory-Huggins free energy measured from the common tangent at the binodal compositions.

Run as a script, it writes the curve gamma_pl R0/kappa = 6.708 (the classical ideal-neck value) in the
(Lambda~, chi~) plane to output/isogamma_6p708.txt.
"""
import os
import numpy as np
from scipy.optimize import brentq

_trapz = getattr(np, 'trapezoid', None) or np.trapz
INV_AP, KAPPA = 100.0, 20.0          # R0^2/a_p, kappa/kBT
GAMMA_IDEAL = 6.708                  # 2 kappa (1/R1 + 1/R2) in units of kappa/R0 at phi0 = 0.2


def xlogx(p):
    return np.where(p > 1e-300, p*np.log(np.where(p > 1e-300, p, 1.0)), 0.0)


def binodal(chi):
    """Protein-poor coexisting fraction of the symmetric Flory-Huggins model (root found in ln phi)."""
    f = lambda u: u - np.log1p(-np.exp(u)) - chi*(2*np.exp(u) - 1)
    return float(np.exp(brentq(f, -700.0, np.log(0.5 - 1e-6))))


def gamma_pl(lam, chi, inv_ap=INV_AP, kappa=KAPPA, n=40001):
    """gamma_pl R0/kappa for Lambda~ = lam, chi~ = chi."""
    pa = binodal(chi); pb = 1 - pa
    g = lambda p: xlogx(p) + xlogx(1-p) + chi*p*(1-p)
    mu = (g(pb) - g(pa))/(pb - pa)
    # integrate over [pa, pb] only, with p = pa + (pb - pa)(1 + tanh s)/2; outside the binodal
    # Delta g is positive again and must not contribute (matters near the critical point chi~ = 2)
    S = max(25.0, 0.7*chi + 15.0)
    s = np.linspace(-S, S, n)
    p = pa + (pb - pa)*0.5*(1 + np.tanh(s))
    with np.errstate(over='ignore'):
        jac = 0.5*(pb - pa)/np.cosh(s)**2              # -> 0 in the tails
    dg = np.clip(g(p) - (g(pa) + mu*(p - pa)), 0.0, None)
    I = _trapz(np.sqrt(dg)*jac, s)
    return 2*inv_ap*np.sqrt(lam)*I/kappa


def chi_at(lam, target, lo=2.6, hi=400.0, **kw):
    """chi~ at which gamma_pl R0/kappa equals target, at fixed Lambda~."""
    return brentq(lambda c: gamma_pl(lam, c, **kw) - target, lo, hi)


if __name__ == '__main__':
    out = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'output')
    os.makedirs(out, exist_ok=True)
    lams = np.geomspace(0.08, 0.56, 25)
    rows = [(lam, chi_at(lam, GAMMA_IDEAL)) for lam in lams]
    np.savetxt(os.path.join(out, 'isogamma_6p708.txt'), np.array(rows), fmt=['%.5f', '%.4f'],
               header='classical ideal-neck contour: planar gamma*R0/kappa = 6.708, phi0 = 0.2, R0^2/a_p = 100\n'
                      'Lambda~   chi~')
    for lam, chi in rows[::4]:
        print(f"Lambda~ = {lam:.3f}   chi~ = {chi:.3f}")
    print("written output/isogamma_6p708.txt")
