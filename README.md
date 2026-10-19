# Diffuse interfaces turn vesicle neck closure into a fold bifurcation

Code accompanying

> W. Xie and R. Ma, *Diffuse interfaces turn vesicle neck closure into a fold bifurcation* (2026).

The paper studies budding of a closed, axisymmetric vesicle carrying a phase-separating protein field
(Helfrich bending energy + Flory–Huggins free energy + Cahn–Hilliard gradient energy). The code computes

* stationary branches of the **full model** by numerical continuation of the shape equations (MATLAB);
* constrained energy landscapes, forces and folds of the **reduced model** (MATLAB and Python);
* the linear stability of the uniform sphere and the energy of the complete-budding reference state
  (Mathematica).

Section numbers below refer to the Supplemental Material (SM) of the paper.

## Parameters and units

The stationary states depend on four dimensionless parameters, with R0 defined by the total area,
A = 4πR0²:

| paper | meaning | value(s) in the paper |
|---|---|---|
| ε = κ a_p/(k_BT R0²) | bending rigidity relative to the entropy of mixing | 0.2 |
| χ̃ = χ/k_BT | Flory–Huggins interaction | varied |
| Λ̃ = Λ/(R0² k_BT) | gradient-energy coefficient (sets the interface width) | varied |
| φ0 | mean protein area fraction | 0.2 |

with κ = 20 k_BT and R0²/a_p = 100 (R0 = 100 nm, a_p = 100 nm²).

**The MATLAB code works in its own units.** The conversion is

| MATLAB variable | value | paper |
|---|---|---|
| `R0` | 1 | length unit |
| `area0` | 2 | A/(2π), i.e. A = 4πR0² |
| `ka` | 0.4 | κ |
| `kT` | 0.02 | k_BT, so κ = 20 k_BT |
| `ap` (full model), `p.a0` (reduced model) | 0.01 | a_p/R0², so ε = 0.2 |
| `Lambda`, `p.Lambda` | Λ̃ × `kT` | Λ̃ = `Lambda`/`kT` (e.g. 0.006 → 0.3) |
| `chi`, `p.chi` | χ̃ × `kT` | χ̃ = `chi`/`kT` |
| `rho0`, `p.phi0` | 0.2 | φ0 |
| energies | — | divide by `kT` to obtain F/k_BT |

The Python code uses the paper's units directly: lengths in R0, energies in k_BT, `lam` = Λ̃, `chi` = χ̃.

## Contents

```
full_model/                 MATLAB, full model (SM Sec. S1)
  shapef.m, twobcf.m        shape equations and boundary conditions, chi fixed
  shapef2.m, twobcf2.m      same, with chi unknown and the pole-to-pole height fixed
  iterate_chi.m             continuation in chi
  iterate_Hchi.m            continuation in the pole-to-pole height (passes folds in chi)
  run_chi_continuation.m    trace a branch in chi from a stored solution
  run_lambda_then_chi.m     move a stored solution to a new Lambda~, then trace the branch in chi
  entropy_cutoff.m          regularization threshold of the mixing entropy (1e-4; SM Sec. S1.5)
  fig1ab_energy_curves.m    energies of a branch, the uniform sphere and Division (Fig. 1(a,b))
  fig1cde_neck_width_tension.m  neck radius, interface width, line tension along a branch (Fig. 1(c-e))
  data/branch_Lambda0p30.mat  751 stationary solutions at Lambda~ = 0.3, chi~ = 3.7-11.6
  data/branch_Lambda1p00_sphere.mat  uniform sphere, chi~ = 0-15.95 (reference curve "Sphere")
reduced_model/
  matlab/reduced_landscape.m  reduced model: F_min(r_n) at given (Lambda~, chi~) (SM Sec. S5)
  matlab/fig2b_barriers.m   barriers of the reduced and the full model versus chi (Fig. 2(b))
  python/                   Python implementation of the same model and the analyses below
mathematica/
  stability_incompressible.nb, Matrices_phi.nb   linear stability of the uniform sphere (SM Sec. S2)
  energy_division.nb        energy of the complete-budding reference state "Division" (SM Sec. S3)
fig3_boundaries.txt         transition points of Fig. 3 (D_bud, L_CB, D_co, L_fold)
data_all/                   full data set, downloaded separately (see Data below)
```

## Requirements

* MATLAB (R2019b or later; `bvp5c`, `fminbnd`, `fminsearch` from base MATLAB; scripts with local functions).
* Python ≥ 3.8 with NumPy, SciPy (≥ 1.6) and Matplotlib; LaTeX is used for figure labels if available.
* Mathematica ≥ 12.1 for the notebooks.

## Running

### Full model

From `full_model/`:

```matlab
run_chi_continuation     % branch at Lambda~ = 0.3, from chi~ = 9.795 through the fold of L_fold
run_lambda_then_chi      % move to Lambda~ = 0.4 at fixed chi~ = 5, then trace the branch in chi
```

Each script saves the computed branch (`cchi`, `hheight`, `ssol`: chi, pole-to-pole height and the
`bvp5c` solution structures) to a `.mat` file. The starting index, the direction of the first stage
and the number of alternating chi/height stages are set at the top of each script. The boundary-value
problem is solved with `bvp5c` (RelTol = 1e-3, AbsTol = 1e-8); the poles are placed at r = 1e-5 R0.
Continuation in chi stops near a fold in chi, where continuation in the height takes over; only
warnings issued by `bvp5c` stop a stage. Continuation in Lambda at fixed chi is a natural-parameter
continuation as well and cannot pass a fold; it is therefore started from a prolate state. In practice a branch at a
new Lambda~ is obtained by choosing a (chi~, Lambda~) pair on an existing branch and stepping Lambda~;
the starting state has to be chosen so that no fold lies in between. At small Lambda~ the prolate part
and the open-neck part of a branch are best obtained separately, each from a starting shape of its own
type, or with a mesh refined near the neck.

The entropy cutoff is set in `entropy_cutoff.m`; the convergence test of SM Sec. S1.5 changes it
there.

From a computed branch, `fig1ab_energy_curves` gives the energy curves of Fig. 1(a,b) and
`fig1cde_neck_width_tension` the neck radius r_n (interior minimum of r), interface width
w = (phi_max - phi_min)/max|dphi/ds| and line tension gamma = (2 Lambda/a_p) int (dphi/ds)^2 ds of
Fig. 1(c-e). Both read `data/branch_Lambda0p30.mat` as a test case; Fig. 1 uses the branches at
Lambda~ = 0.1, 0.4, 0.8 and 1.75 (`data_all/full_model/branch_Lambda0p10.mat`, `..._0p40`, `..._0p80`,
`..._1p75`), set as the input file at the top of each script.

The boundaries of Fig. 3 are read off the computed branches: D_bud, L_CB and D_co are crossings of
energy curves (located by interpolation between stored solutions), C_an is analytical, and L_fold is
the fold (largest chi~) of the stable open-neck part. The points are listed in `fig3_boundaries.txt`.

Expected results: `run_chi_continuation` passes the fold of the open-neck branch at chi~ = 11.574
(L_fold at Lambda~ = 0.3); `run_lambda_then_chi` traces the Lambda~ = 0.4 branch through the end of the
prolate part at chi~ = 8.95 to the start of the stable open-neck part at chi~ = 8.34 (Fig. 1(d)).

### Reduced model

`reduced_model/matlab/reduced_landscape.m` computes F_min(r_n) at the (Lambda~, chi~) set at its top
(Fig. 2(a)). `fig2b_barriers(0.15)`, `fig2b_barriers(0.30)` and `fig2b_barriers(0.45)` in
`reduced_model/matlab/` read such landscapes and the two full-model branches that bound the barrier
from `data_all/reduced_model/`, and plot both barriers versus chi~ (Fig. 2(b)).
The Python scripts, run from `reduced_model/python/`, write to `reduced_model/python/output/`:

| script | result | paper |
|---|---|---|
| `fig2_landscapes.py` | F_min(r_n) at Λ̃ = 0.15 (χ̃ = 15, 18, 21, 25), 0.3 and 0.45 | Fig. 2(c); cross-check of Fig. 2(a) |
| `fig2_forces.py` | Π_b, Π_p; position and depth of the weak point; Π_b(r_n → 0) vs 2πγ_ideal | Fig. 2(c), main text |
| `fig2d_pull_vs_lambda.py` | Π_p/Π_p⁰ for Λ̃ = 0.05–0.45 at γ̃_pl = 6.5 | Fig. 2(d) and inset |
| `figS6_small_lambda.py` | Π_b, Π_p for Λ̃ = 0.15, 0.05, 0.02 | SM Fig. S6 |
| `reduced_folds.py` | folds of the reduced model, r_n and γ̃_pl at the fold | Fig. 3 (L_fold^re), SM Sec. S9 |
| `ap_dependence.py` | dependence on a_p; reduced-model estimate of L_CB | SM Sec. S10 |
| `line_tension.py` | planar line tension γ̃_pl = γ_pl R0/κ; contour γ̃_pl = 6.708 | Fig. 3, SM Sec. S9.2 |

`fig2_forces.py` reads the output of `fig2_landscapes.py`. The other scripts are independent. Run times
range from seconds (`line_tension.py`) to tens of minutes (`reduced_folds.py`; use `--parallel` to run
the four Λ̃ at once).

## Data

The stationary branches and landscapes behind the figures (6.7 GB uncompressed) are too large for GitHub;
they can be regenerated with the code here and are available from the authors upon reasonable request.
The scripts expect the data set in the repository root as `data_all/`:

```
data_all/full_model/branch_LambdaXpXX.mat     stationary branches at 18 values of Lambda~ from 0.075 to 2.0
                                              (cchi, hheight, ssol, Lambda or La; MATLAB code units)
data_all/reduced_model/landscape_LambdaXpXX_chiYYpYY.mat
                                              reduced-model landscapes (output of reduced_landscape.m)
                                              at Lambda~ = 0.15, 0.30, 0.45
data_all/reduced_model/branch_LambdaXpXX_inner.mat, ..._outer.mat
                                              the two full-model branches between which the
                                              full-model barrier of Fig. 2(b) is taken
```

`full_model/data/` holds the two small files used by the examples: the branch at Lambda~ = 0.3 and the
uniform sphere.

## License

MIT; see `LICENSE`.

## Citation

If you use this code, please cite the paper above.
