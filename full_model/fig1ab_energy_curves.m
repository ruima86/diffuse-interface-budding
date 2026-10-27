% FIG1AB_ENERGY_CURVES  Energy curves of Fig. 1(a,b): deformed branch, uniform sphere and the
% complete-budding reference (Complete bud), as E~ = F/(8*pi*kappa) and as (F - F_CB)/kT, the
% quantity plotted in Fig. 1(a,b).
% Run this script from MATLAB. Change only the two input files below.
% The included Lambda~ = 0.30 branch is a test case, NOT a Fig. 1 branch;
% Fig. 1(a,b) require branch data at Lambda~ = 0.075, 0.40, 0.80, 1.75.
clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
sphere = load(fullfile(scriptDir,'data','branch_Lambda1p00_sphere.mat'));
S = load(fullfile(scriptDir,'data','branch_Lambda0p30.mat'));

% README code units: R0=1, kT=0.02, kappa=0.4, ap=0.01, phi0=0.2.
% Stored chi and Lambda are in code units. Plot/output axes use the paper's
% dimensionless chi~ = chi/kT and E~ = F/(8*pi*kappa) for Fig. 1.
R0 = 1;
kT = 0.02;
kappa = 0.4;
ap = 0.01;
phi0 = 0.2;
epsilon = kappa*ap/(kT*R0^2);

lambdaBranch = storedLambda(S);
lambdaSphere = storedLambda(sphere);
lambdaTilde = lambdaBranch/(R0^2*kT);

[chiBranch,Ebranch,Ebend,Eprotein] = branchEnergy(S,lambdaBranch,kT,kappa,ap);
[chiSphere,EsphereNumerical] = branchEnergy(sphere,lambdaSphere,kT,kappa,ap);

% Keep the original Sphere curve on the chi range actually stored in its
% input file; do not extend it by extrapolation.
Esphere = EsphereNumerical;

% Independent exact check of the sphere reference (SM Sec. S3).
EsphereExact = 1 + (chiSphere*phi0*(1-phi0)+entropyIdeal(phi0))/(2*epsilon);
sphereError = max(abs(EsphereNumerical-EsphereExact));
if sphereError > 1e-4
    warning('Sphere data differ from the uniform-sphere formula by %.3g.',sphereError);
end

% Complete bud is the minimized energy of two spheres joined at a closed neck
% (SM Sec. S3); it is independent of Lambda~. Its free-energy density uses the
% same regularized entropy as the stationary solutions (entropy_cutoff.m), so
% that energy differences and the crossing L_CB are consistent: for
% chi~ >~ 9.2 the protein-poor phase lies inside the regularized window.
chiDivision = (0:0.05:90)';
Edivision = zeros(size(chiDivision));
gReg = @(x,chi) chi*x.*(1-x)+entropyRegularized(x);
for i = 1:numel(chiDivision)
    chi = chiDivision(i);
    phiLow = fminbnd(@(x) gReg(x,chi),-0.02,0.5,optimset('TolX',1e-14));
    gDivision = gReg(phi0,chi);
    if phiLow < phi0, gDivision = min(gDivision,gReg(phiLow,chi)); end
    Edivision(i) = 2+gDivision/(2*epsilon);
end

% All numerical output below is dimensionless; the original .mat files
% remain in the solver's code units and are not modified.
curves = struct();
curves.lambda_tilde = lambdaTilde;
curves.deformed.chi_tilde = chiBranch;
curves.deformed.E_tilde = Ebranch;
curves.deformed.E_bending_tilde = Ebend;
curves.deformed.E_protein_tilde = Eprotein;
curves.sphere.chi_tilde = chiSphere;
curves.sphere.E_tilde = Esphere;
curves.complete_bud.chi_tilde = chiDivision;
curves.complete_bud.E_tilde = Edivision;
% Fig. 1(a,b): energies relative to complete budding, (F - F_CB)/kT.
toKT = 8*pi*kappa/kT;
curves.deformed.dF_kT = (Ebranch-interp1(chiDivision,Edivision,chiBranch))*toKT;
curves.sphere.dF_kT = (Esphere-interp1(chiDivision,Edivision,chiSphere))*toKT;

figure;
scatter(chiBranch,Ebranch,4,'filled'); hold on;
plot(chiSphere,Esphere,'-.','LineWidth',1.2);
plot(chiDivision,Edivision,'-.','LineWidth',1.2);
xlabel('$\widetilde\chi$','Interpreter','latex');
ylabel('$\widetilde E=F/(8\pi\kappa)$','Interpreter','latex');
legend('Deformed','Sphere','Complete bud','Location','best');
title(['$\widetilde\Lambda=',num2str(lambdaTilde,'%.3g'),'$'], ...
    'Interpreter','latex');
box on;
fprintf('Lambda~ = %.4g; %d deformed states; %d sphere states.\n', ...
    lambdaTilde,numel(chiBranch),numel(chiSphere));
fprintf('Maximum sphere numerical/analytic difference: %.3g in E~.\n',sphereError);

function lambdaCode = storedLambda(data)
% Older branch files use La; new continuation files use Lambda.
if isfield(data,'Lambda')
    lambdaCode = data.Lambda;
elseif isfield(data,'La')
    lambdaCode = data.La;
else
    error('Each input .mat file must contain its Lambda or La value.');
end
if ~isscalar(lambdaCode) || ~isfinite(lambdaCode) || lambdaCode < 0
    error('Stored Lambda/La must be one nonnegative finite number.');
end
end

function [chiTilde,Etilde,EbendTilde,EproteinTilde] = ...
        branchEnergy(data,lambdaCode,kT,kappa,ap)
if ~isfield(data,'cchi') || ~isfield(data,'ssol') || ...
        numel(data.cchi) ~= numel(data.ssol)
    error('Each input .mat file must contain matching cchi and ssol arrays.');
end
n = numel(data.cchi);
chiTilde = data.cchi(:)/kT;
Etilde = zeros(n,1);
EbendTilde = zeros(n,1);
EproteinTilde = zeros(n,1);
energyScale = 8*pi*kappa;
for i = 1:n
    sol = data.ssol{i};
    u = linspace(sol.x(1),sol.x(end),10000);
    y = deval(sol,u);
    h = sol.parameters(1);
    psi = y(1,:);
    dpsi = y(2,:);
    r = y(3,:);
    phi = y(7,:);
    dphi = y(8,:);
    if h <= 0 || any(r <= 0)
        error('Nonpositive meridian length or radius in solution %d.',i);
    end

    % dA=2*pi*h*r*du. C0, pressure and external force are zero here.
    twoH = dpsi/h + sin(psi)./r;
    bendDensity = 0.5*kappa*twoH.^2;
    proteinDensity = (lambdaCode*(dphi/h).^2 + ...
        data.cchi(i)*phi.*(1-phi) + kT*entropyRegularized(phi))/ap;
    EbendTilde(i) = 2*pi*trapz(u,h*r.*bendDensity)/energyScale;
    EproteinTilde(i) = 2*pi*trapz(u,h*r.*proteinDensity)/energyScale;
    Etilde(i) = EbendTilde(i)+EproteinTilde(i);
end
end

function value = entropyIdeal(phi)
% x*log(x) is zero by continuity at x=0.
value = zeros(size(phi));
middle = phi > 0 & phi < 1;
value(middle) = phi(middle).*log(phi(middle)) + ...
    (1-phi(middle)).*log1p(-phi(middle));
value(phi < 0 | phi > 1) = NaN;
end

function value = entropyRegularized(phi)
% Fourth-order Taylor continuation for phi outside [delta,1-delta],
% matching SHAPEF and SM Sec. S1.5.
delta = entropy_cutoff();
q = min(phi,1-phi);
value = zeros(size(q));
middle = q >= delta;
value(middle) = entropyIdeal(q(middle));
d = q(~middle)-delta;
h0 = entropyIdeal(delta);
h1 = log(delta/(1-delta));
h2 = 1/delta + 1/(1-delta);
h3 = -1/delta^2 + 1/(1-delta)^2;
h4 = 2/delta^3 + 2/(1-delta)^3;
value(~middle) = h0+h1*d+h2*d.^2/2+h3*d.^3/6+h4*d.^4/24;
end
