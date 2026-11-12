function delta = entropy_cutoff()
% ENTROPY_CUTOFF  Regularization threshold of the mixing entropy (SM Sec. I E).
% For phi < delta and phi > 1 - delta, the entropic term phi*log(phi) +
% (1-phi)*log(1-phi) is replaced by its fourth-order Taylor expansion about
% phi = delta and phi = 1 - delta. Used by SHAPEF, SHAPEF2, TWOBCF, TWOBCF2 and
% FIG1AB_ENERGY_CURVES. The convergence test of SM Sec. I E uses 1e-4, 1e-5, 1e-6.
delta = 1e-4;
end
