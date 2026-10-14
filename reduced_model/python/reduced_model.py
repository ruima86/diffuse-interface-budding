"""Reduced (four-region) budding model, Supplemental Material Sec. S5.

Python implementation of the same model as ../matlab/reduced_landscape.m.

Units: lengths in R0 (R0 = 1), energies in kBT. kappa = 20 kBT, R0^2/a_p = 100, A0 = 4 pi, phi0 = 0.2.
In the code lam = Lambda~ and chi = chi~.
    f = (kappa/2)(2H)^2 + (R0^2/a_p) [ lam (dphi/ds)^2 + chi phi(1-phi) + phi ln phi + (1-phi) ln(1-phi) ]
    F(a, xi) = 8 pi kappa + F_II + F_III ,   F_II/III = 2 pi int r J f dtheta
Caps I and IV are spheres (pure phases, bending energy 4 pi kappa each). The neck curvature is prescribed,
|2H(0)| = C(a) = 2 - 0.68 a, which fixes omega_pm = sqrt((R-a)/((R+a)(1-aC))).
R_- and R_+ are fixed by area and protein conservation for each (a, xi); F_min(a) = min over xi.
The neck radius a is r_n~ of the paper.

To change a_p, set INV_AP = (R0^2/a_p) before calling landscape().
"""
import numpy as np
from scipy.optimize import root, minimize_scalar
from scipy.integrate import cumulative_trapezoid

_trapz = getattr(np, 'trapezoid', None) or np.trapz

KAPPA, INV_AP, A0, PHI0 = 20.0, 100.0, 4*np.pi, 0.2
NTH = 6000
_u = np.linspace(0.0, 1.0, NTH)
THETA = (np.pi/2)*_u**4                     # clustered at the neck (theta = 0)


def C_neck(a):
    return 2.0 - 0.68*a


def region(a, R):
    """Geometry of one connecting region on THETA: r, J = ds/dtheta, 2H, ell (unsigned arclength from neck)."""
    th = THETA; s, c = np.sin(th), np.cos(th)
    D = R*R - a*a
    r = np.sqrt(a*a + D*s*s)
    r_t = D*s*c/r
    r_tt = D*(c*c - s*s)/r - (D*s*c)**2/r**3
    om = np.sqrt((R - a)/((R + a)*(1.0 - a*C_neck(a))))
    z_t = om*(r + a*R/r)                      # region III sign; region II is the mirror image (same |H|, J)
    z_tt = om*r_t*(1.0 - a*R/r**2)
    J = np.sqrt(r_t**2 + z_t**2)
    twoH = (r_t*z_tt - z_t*r_tt)/J**3 + z_t/(r*J)
    ell = cumulative_trapezoid(J, th, initial=0.0)
    return r, J, twoH, ell


def xlogx(p):
    return np.where(p > 0, p*np.log(np.where(p > 0, p, 1.0)), 0.0)


def profile(ell_signed, lm, lp, xi):
    """tanh profile of width xi, normalized to phi = 0 and 1 at the ends of regions II and III."""
    den = np.tanh(lp/xi) - np.tanh(lm/xi)
    with np.errstate(over='ignore'):
        dphi = (1.0/xi)/np.cosh(ell_signed/xi)**2/den      # d phi / d ell ; |d phi/ds| = |d phi/d ell|
    phi = (np.tanh(ell_signed/xi) - np.tanh(lm/xi))/den
    return np.clip(phi, 0.0, 1.0), dphi


def integrals(a, Rm, Rp, xi):
    rm, Jm, Hm, lm_ = region(a, Rm)
    rp, Jp, Hp, lp_ = region(a, Rp)
    lm, lp = -lm_[-1], lp_[-1]                               # ell_-^b < 0 < ell_+^b
    phim, dpm = profile(-lm_, lm, lp, xi)
    phip, dpp = profile(lp_, lm, lp, xi)
    T = THETA
    area = 2*np.pi*(Rm**2 + Rp**2 + _trapz(rm*Jm, T) + _trapz(rp*Jp, T))
    aphi = 2*np.pi*(Rp**2 + _trapz(rm*phim*Jm, T) + _trapz(rp*phip*Jp, T))
    return (rm, Jm, Hm, phim, dpm), (rp, Jp, Hp, phip, dpp), area, aphi


def solve_R(a, xi, guess):
    def eqs(x):
        Rm, Rp = x
        if Rm <= a*1.0001 or Rp <= a*1.0001:
            return [1e3, 1e3]
        _, _, area, aphi = integrals(a, Rm, Rp, xi)
        return [area - A0, aphi - PHI0*A0]
    sol = root(eqs, guess, method='hybr', options=dict(xtol=1e-12))
    return sol.x, sol.success


def energy(a, xi, lam, chi, guess=(np.sqrt(1-PHI0), np.sqrt(PHI0))):
    """Total F and its parts (bending incl. caps, gradient, bulk) in kBT, plus (Rm, Rp)."""
    (Rm, Rp), ok = solve_R(a, xi, guess)
    if not ok:
        return None
    reg = integrals(a, Rm, Rp, xi)
    Fb = 8*np.pi*KAPPA; Fg = 0.0; Fbulk = 0.0
    for (r, J, twoH, phi, dphi) in reg[:2]:
        w = 2*np.pi*r*J
        Fb += _trapz(w*(KAPPA/2)*twoH**2, THETA)
        Fg += _trapz(w*INV_AP*lam*dphi**2, THETA)
        Fbulk += _trapz(w*INV_AP*(chi*phi*(1-phi) + xlogx(phi) + xlogx(1-phi)), THETA)
    return dict(F=Fb+Fg+Fbulk, bend=Fb, grad=Fg, bulk=Fbulk, Rm=Rm, Rp=Rp)


def e_min(a, lam, chi, guess, xi_bounds=(-7.0, 0.5)):
    """Minimize over xi (bounded in ln xi). Returns the energy dict at xi*, with 'xi' added."""
    def obj(lx):
        e = energy(a, np.exp(lx), lam, chi, guess)
        return np.inf if e is None else e['F']
    res = minimize_scalar(obj, bounds=xi_bounds, method='bounded', options=dict(xatol=1e-6))
    e = energy(a, np.exp(res.x), lam, chi, guess)
    e['xi'] = np.exp(res.x)
    return e


def landscape(lam, chi, a_grid):
    """Constrained landscape on a_grid. Rows: (a, F, bend, grad, bulk, xi*, R_-, R_+)."""
    rows, guess = [], (np.sqrt(1-PHI0), np.sqrt(PHI0))
    for a in a_grid:
        e = e_min(a, lam, chi, guess)
        guess = (e['Rm'], e['Rp'])
        rows.append((a, e['F'], e['bend'], e['grad'], e['bulk'], e['xi'], e['Rm'], e['Rp']))
    return np.array(rows)


def forces(L):
    """Bending resistance Pi_b = -d(bend)/da and interfacial pull Pi_p = d(grad + bulk)/da (kBT/R0)
    along the optimal width xi*(a); Pi_p - Pi_b = dF_min/da."""
    a = L[:, 0]
    return a, -np.gradient(L[:, 2], a), np.gradient(L[:, 3] + L[:, 4], a)


def extrema(L):
    """Neck radii of the local minimum (stable open-neck state) and maximum (barrier top) of F_min,
    from log-linear interpolation of the zeros of dF/da, and the barrier F(max) - F(min).
    Entries are None when absent. If several sign changes occur, the last of each kind is kept."""
    a, F = L[:, 0], L[:, 1]
    dF = np.gradient(F, a); s = np.sign(dF); mn = mx = None
    for i in np.where(np.diff(s) != 0)[0]:
        x0, x1 = np.log(a[i]), np.log(a[i+1]); y0, y1 = dF[i], dF[i+1]
        ar = np.exp(x0 - y0*(x1 - x0)/(y1 - y0))
        if y0 < 0 < y1:
            mn = ar
        elif y0 > 0 > y1:
            mx = ar
    b = None
    if mn is not None and mx is not None:
        b = np.interp(np.log(mx), np.log(a), F) - np.interp(np.log(mn), np.log(a), F)
    return mn, mx, b
