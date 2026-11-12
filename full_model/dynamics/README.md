# Dynamics of the full model (SM Sec. IX, Fig. S7)

Time integration of the full model at Λ̃ = 0.4 (`Lambda` = 0.008) and φ0 = 0.2, used in SM Sec. IX to
test the stability of the open-neck states and the closure of the neck above the fold.

The dynamics follows from Onsager's variational principle: the free energy of the stationary problem
plus a dissipation from the membrane shear viscosity η_s (`etas` = 1), tangential friction Γ (`Ga2` =
0.01; normal friction `Ga1` = 0) and protein diffusion D (`Di` = 0.00312). Membrane area is conserved
locally and protein number globally; the pressure difference is zero and the shapes stay axisymmetric.
In the paper's units η_s D/κ = 7.8e-3 and Γ R0²/η_s = 0.01. Time in the code is t; the paper uses
t/τ_D = `Di`·t, with τ_D = R0²/D. The other code units are those of the main README.

## Files

| file | content |
|---|---|
| `shapef.m`, `twobcf.m` | equations of one implicit time step and their boundary conditions (solved with `bvp5c`) |
| `sol_La_0.008_chi_0.18.mat` | stationary stable open-neck state at χ̃ = 9 (r_n = 0.091 R0) |
| `sol_La_0.008_chi_0.18_2.mat` | stationary unstable open-neck state at χ̃ = 9 (r_n = 0.039 R0) |
| `run_trajectory_1.m` | stable state, neck widened by 5% plus a local protein redistribution adding 0.2 k_BT; returns (trajectory 1) |
| `run_trajectory_2_opening.m` | unstable state, neck widened by 5%; relaxes to the stable state (trajectory 2) |
| `run_trajectory_2_closing.m` | unstable state, neck narrowed by 5%; the neck narrows to r_n < 0.003 R0 (trajectory 3) |
| `run_trajectory_3.m` | stable state at χ̃ = 9 evolved at χ̃ = 9.5, above the fold χ̃_c = 9.25; the neck narrows to r_n = 0.004 R0 (Fig. S7(c)–(e)) |

Note that `shapef.m` and `twobcf.m` here differ from the stationary-problem files of the same name in
`full_model/`. Each script adds its own folder to the MATLAB path; run the dynamics in a fresh MATLAB
session, or remove this folder from the path before returning to the stationary code.

## Running

Open a script and run it; no inputs are needed. The perturbation (`neckChange`, the relative change of
the neck radius, and `proteinEnergy`, in k_BT) is applied by the local function at the end of each
script, which checks that membrane area and protein number are unchanged. Each step solves the
boundary-value problem with `bvp5c` (RelTol = 1e-3, AbsTol = 1e-8). The output is written to `data/`
next to the script as `sol_all`, a cell array with one row per step: time t, the `bvp5c` solution
(rows of `y`: 1 ψ, 2 dψ/du, 3 r, 4 z, 12 φ, 13 dφ/du; `parameters(1)` is the arclength h), the normal
velocity, and dF/dt − dW/dt, which must be negative. The first row is the state after one step.

Each time step is implicit: the shape and composition at t + dt solve a boundary-value problem in
which the rates are the differences from the configuration at t divided by dt, solved with `bvp5c` on
an adaptive mesh of at most `NMax` = 5000 points. The scripts run with a fixed step `dt` until t = `T`
or until the `bvp5c` error exceeds RelTol.

The trajectories of the paper were obtained in the same way, but continued from the last saved state
with a different step, without a new perturbation (code time t; t/τ_D = 0.00312 t):

| script | `dt` schedule | end |
|---|---|---|
| `run_trajectory_1.m` | 1 (t ≤ 703), then 10 | t = 2593 (8.1 τ_D) |
| `run_trajectory_2_opening.m` | 1 (t ≤ 752), then 10 | t = 5252 (16.4 τ_D) |
| `run_trajectory_2_closing.m` | 0.1 (t ≤ 9.2), 0.01 (t ≤ 9.35), 0.001 | t = 9.36 (0.029 τ_D), r_n = 0.0029 R0 |
| `run_trajectory_3.m` | 1 (t ≤ 20), 0.1 (t ≤ 21.3), 0.01 (t ≤ 21.44), 0.001 | t = 21.456 (0.067 τ_D), r_n = 0.0038 R0 |

The two closing runs were stopped by hand once r_n < 0.005 R0. Across each reduction of the step the
closing rate dr_n/dt changes by 3–6%, comparable to its change over the following step, so the closing
trajectories show no visible dependence on the step. The free energy, recomputed from the solutions,
decreases monotonically to within 0.002 k_BT.

The neck radius r_n is the interior minimum of r; the interface width is
w = (φ_max − φ_min)/max|dφ/ds|, as in `fig1cde_neck_width_tension.m`.
