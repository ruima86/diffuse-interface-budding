% RUN_CHI_CONTINUATION  Trace a stationary branch of the full model in chi at fixed
% Lambda, starting from a stored solution.
%
% The branch is followed by alternating continuation in chi (ITERATE_CHI) and in the
% pole-to-pole height H (ITERATE_HCHI); the latter passes the folds of the branch in
% chi. Each stage runs until bvp5c stops, then the other stage takes over.
%
% Code units and the dimensionless parameters of the paper:
%   R0 = 1, total area A = 4*pi*R0^2, area0 = A/(2*pi) = 2
%   ka = kappa = 0.4, kT = 0.02           ->  kappa = 20 kBT
%   ap = 0.01                             ->  R0^2/a_p = 100, eps = ka*ap/(kT*R0^2) = 0.2
%   Lambda (code) = Lambda~ * kT          ->  Lambda~ = Lambda/kT   (0.006 -> 0.3)
%   chi    (code) = chi~ * kT             ->  chi~    = chi/kT
%   rho0 = phi0 = 0.2, C0 = 0, pressure difference pre = 0
%   Energies in code units divided by kT give F/kBT.

clear; clc;
S = load(fullfile('data','branch_Lambda0p30.mat'));   % stored branch at Lambda~ = 0.3: cchi, hheight, ssol

% Geometry
R0 = 1;
area0 = 2*R0^2;          % A/(2*pi)

% Membrane
pre = 0;                 % pressure difference
ka = 0.4;                % bending rigidity kappa

% Proteins
rho0 = 0.2;              % mean protein area fraction phi0
kT = 0.02;
ap = 0.01;               % area of a phase-separating unit, a_p/R0^2
C0 = 0;                  % spontaneous curvature

% Forces not used in this work
fext = 0;                % axial point force
V0 = 0;                  % soft repulsion near the axis (see SHAPEF)

% Phase separation
Lambda = 0.006;          % Lambda~ = Lambda/kT = 0.3 (must match the stored branch)

% Numerics
ri = 1e-5;               % poles placed at r = ri
RTol = 1e-3;             % bvp5c relative tolerance (AbsTol = 1e-8)
Nmax = 50000;            % bvp5c maximum number of mesh points

% Starting point: the stored solution closest to chiStart
chiStart = 9.795;
[~, kstart] = min(abs(S.cchi/kT - chiStart));
dir0 = 1;                % first stage: 1 = increasing chi, 2 = decreasing chi
nstage = 4;              % number of alternating stages (chi, H, chi, ...)
increment = [0.001, -0.001];   % step in chi (code units, i.e. 0.05 in chi~) or in H

sol = S.ssol{kstart};
chi = S.cchi(kstart);
solinit.parameters(1) = sol.parameters(1);
solinit.parameters(2) = sol.parameters(2);
solinit.parameters(3) = sol.parameters(3);
solinit.x = sol.x;
solinit.y = sol.y;
fprintf('Start: Lambda~ = %.3f, chi~ = %.4f\n', Lambda/kT, chi/kT);

yeq = @(u,y,para)shapef(u,y,para,ka,C0,kT,Lambda,chi,pre,fext,ap,V0);
ybc = @(ya,yb,para)twobcf(ya,yb,para,ka,C0,kT,ri,Lambda,chi,ap,area0,rho0);
opts = bvpset('RelTol',RTol,'AbsTol',1e-8,'NMax',Nmax);
sol = bvp5c(yeq,ybc,solinit,opts);
para_end = sol.parameters;
sol_end = sol;

cchi = [];
hheight = [];
ssol = {};
flag = true;
for k = 1:nstage
    fprintf('Stage %d: continuation in %s\n', k, ternary(flag, 'chi', 'H'));
    if flag
        [sol_end,para_end,cchi,hheight,ssol] = iterate_chi(sol_end,para_end,increment,cchi,hheight,...
            ka,pre,ri,area0,fext,rho0,C0,kT,ssol,Lambda,ap,RTol,Nmax,V0,chi,dir0);
    else
        [sol_end,para_end,cchi,hheight,ssol] = iterate_Hchi(sol_end,para_end,increment,cchi,hheight,...
            ka,pre,ri,area0,fext,rho0,C0,kT,ssol,Lambda,ap,RTol,Nmax,V0);
    end
    flag = ~flag;
end

outfile = sprintf('branch_Lambda%.3f_from_chi%.3f.mat', Lambda/kT, chi/kT);
save(outfile, 'cchi', 'hheight', 'ssol', 'Lambda', 'kT', 'ka', 'ap', 'rho0');
fprintf('Saved %d solutions to %s\n', numel(cchi), outfile);

figure(2);
plot(cchi/kT, hheight);
xlabel('\chi/k_BT'); ylabel('pole-to-pole height / R_0');

function out = ternary(c, a, b)
if c, out = a; else, out = b; end
end
