% REDUCED_LANDSCAPE  Constrained energy landscape F_min(r_n) of the reduced budding model
% (Supplemental Material, Sec. S5).
%
% The geometry is prescribed: two spherical caps (regions I and IV, radii Ra = R_- and
% Rb = R_+) joined through a neck of radius aNeck = r_n by two connecting regions
% (II and III) of the Juelicher-Lipowsky four-region type. The axial stretch factors
% gammaA, gammaB (omega_- and omega_+ in the SM) follow from the prescribed neck curvature
% |2H(0)| = C(a) = 2 - 0.68 a, a least-squares fit to the full stationary branches.
% The protein fraction is 0 on cap I, 1 on cap IV and a tanh profile of width xi along
% the arclength of regions II and III. For each neck radius a, Ra and Rb are fixed by
% conservation of the total area and protein number, and F is minimized over xi.
%
% Code units: R0 = 1, A0 = 4*pi, kappa = 0.4, kBT = 0.02, a0 = a_p/R0^2 = 0.01, so that
% Lambda~ = p.Lambda/p.kBT and chi~ = p.chi/p.kBT; the defaults below are Lambda~ = 0.15,
% chi~ = 25. Energies in the output divided by p.kBT give F/kBT.
%
% p.enforceGammaMax rejects configurations with omega > 1; for the parameters of the
% paper omega_pm < 1 throughout, so this safeguard is never active.
% Output: one .mat file per (Lambda, chi), updated after every neck radius.
clear; clc; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 1. Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p.Lambda = 0.003;
p.chi    = 0.5;
xiMinAbs = 0.0001;
xiMaxAbs = 0.4;
p.xi0 = 0;
aVec = logspace(log10(0.0001), log10(0.4), 120);
% Prescribed neck curvature C(a) = k*a + b (SM Sec. S5)
p.C2H_k =-0.68;
p.C2H_b =2 ;

p.A0   = 4*pi;
p.phi0 = 0.2;
p.kappaAlpha = 0.4;
p.kappaBeta  = 0.4;
p.kBT = 0.02;
p.a0  = 0.01;
 
p.enforceGammaMax = true;
p.gammaMax = 1;

p.NthetaScan = 6000;
p.radiiTol = 1e-6;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 2. Allocate arrays %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nA = length(aVec);

xiStar     = NaN(size(aVec));
wEffStar   = NaN(size(aVec));
Fmin       = NaN(size(aVec));
RaStar     = NaN(size(aVec));
RbStar     = NaN(size(aVec));
gammaAStar = NaN(size(aVec));
gammaBStar = NaN(size(aVec));
C2HStar    = NaN(size(aVec));

areaErr    = NaN(size(aVec));
proteinErr = NaN(size(aVec));

aOverXi    = NaN(size(aVec));
aOver2Xi   = NaN(size(aVec));
aOverWEff  = NaN(size(aVec));
aOver2WEff = NaN(size(aVec));

l1Arr = NaN(size(aVec));
l2Arr = NaN(size(aVec));
DwArr = NaN(size(aVec));

AIarr     = NaN(size(aVec));
AIIarr    = NaN(size(aVec));
AIIIarr   = NaN(size(aVec));
AIVarr    = NaN(size(aVec));
AtotalArr = NaN(size(aVec));

AproteinIIArr    = NaN(size(aVec));
AproteinIIIArr   = NaN(size(aVec));
AproteinIVArr    = NaN(size(aVec));
AproteinTotalArr = NaN(size(aVec));

FbendIArr   = NaN(size(aVec));
FbendIIArr  = NaN(size(aVec));
FbendIIIArr = NaN(size(aVec));
FbendIVArr  = NaN(size(aVec));

FbendArr    = NaN(size(aVec));
FcompArr    = NaN(size(aVec));
FentropyArr = NaN(size(aVec));
FtotArr     = NaN(size(aVec));

geomStar = cell(size(aVec));

phiTotalArr = NaN(size(aVec));

matFileRunning = sprintf('La_%g_chi_%g.mat', p.Lambda, p.chi);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 3. Optimize xi for each a using scan Ntheta %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

optionsScan = optimset('TolX', 1e-6, 'Display', 'off');

pScan = p;
pScan.Ntheta = p.NthetaScan;

for i = 1:nA

    pLocal = pScan;
    aNeck = aVec(i);

    xiLow  = xiMinAbs;
    xiHigh = xiMaxAbs;

    objective = @(xi) FtotToy(aNeck, xi, pLocal);

    [xiOptLocal, ~] = fminbnd(objective, xiLow, xiHigh, optionsScan);

    [fCheck, geom] = FtotToy(aNeck, xiOptLocal, pLocal);

    if isfinite(fCheck) && fCheck < 1e20 && ~isempty(geom)

        xiStar(i)     = xiOptLocal;
        wEffStar(i)   = geom.wEff;
        Fmin(i)       = fCheck;
        FtotArr(i)    = geom.Ftot;

        RaStar(i)     = geom.Ra;
        RbStar(i)     = geom.Rb;
        gammaAStar(i) = geom.gammaAlpha;
        gammaBStar(i) = geom.gammaBeta;
        C2HStar(i)    = geom.C2Htarget;

        aOverXi(i)    = aNeck / xiOptLocal;
        aOver2Xi(i)   = aNeck / (2*xiOptLocal);
        aOverWEff(i)  = aNeck / geom.wEff;
        aOver2WEff(i) = aNeck / (2*geom.wEff);

        areaErr(i)    = (geom.Atotal - pLocal.A0) / pLocal.A0;
        proteinErr(i) = (geom.Aprotein - pLocal.A0*pLocal.phi0) / ...
                        (pLocal.A0*pLocal.phi0);

        l1Arr(i) = geom.l1;
        l2Arr(i) = geom.l2;
        DwArr(i) = geom.Dw;

        AIarr(i)     = geom.AI;
        AIIarr(i)    = geom.AII;
        AIIIarr(i)   = geom.AIII;
        AIVarr(i)    = geom.AIV;
        AtotalArr(i) = geom.Atotal;

        AproteinIIArr(i)    = geom.AproteinII;
        AproteinIIIArr(i)   = geom.AproteinIII;
        AproteinIVArr(i)    = geom.AproteinIV;
        AproteinTotalArr(i) = geom.Aprotein;

        % Total protein concentration on the whole membrane
        phiTotalArr(i) = geom.Aprotein / geom.Atotal;

        FbendIArr(i)   = geom.FbendI;
        FbendIIArr(i)  = geom.FbendII;
        FbendIIIArr(i) = geom.FbendIII;
        FbendIVArr(i)  = geom.FbendIV;

        FbendArr(i)    = geom.Fbend;
        FcompArr(i)    = geom.Fcomp;
        FentropyArr(i) = geom.Fentropy;

        geomStar{i} = geom;
    end

    fprintf(['a = %.5g, xi*_scan = %.5g, wEff*_scan = %.5g, ', ...
             'a/xi* = %.5g, a/(2wEff*) = %.5g, ', ...
             'C(a) = %.5g, Ralpha = %.5g, Rbeta = %.5g, ', ...
             'gammaA = %.5g, gammaB = %.5g, ', ...
             'Atotal = %.10g, Aprotein = %.10g, phiTotal = %.10g, ', ...
             'targetPhi = %.10g, areaErr = %.3e, proteinErr = %.3e, ', ...
             'F_scan = %.8g\n'], ...
             aNeck, xiStar(i), wEffStar(i), aOverXi(i), aOver2WEff(i), ...
             C2HStar(i), RaStar(i), RbStar(i), ...
             gammaAStar(i), gammaBStar(i), ...
             AtotalArr(i), AproteinTotalArr(i), phiTotalArr(i), ...
             pLocal.phi0, areaErr(i), proteinErr(i), Fmin(i));

    save(matFileRunning, ...
        'p', 'pScan', ...
        'aVec', 'xiMinAbs', 'xiMaxAbs', ...
        'xiStar', 'wEffStar', 'Fmin', 'FtotArr', ...
        'RaStar', 'RbStar', ...
        'gammaAStar', 'gammaBStar', 'C2HStar', ...
        'aOverXi', 'aOver2Xi', 'aOverWEff', 'aOver2WEff', ...
        'areaErr', 'proteinErr', ...
        'l1Arr', 'l2Arr', 'DwArr', ...
        'AIarr', 'AIIarr', 'AIIIarr', 'AIVarr', 'AtotalArr', ...
        'AproteinIIArr', 'AproteinIIIArr', 'AproteinIVArr', 'AproteinTotalArr', ...
        'phiTotalArr', ...
        'FbendIArr', 'FbendIIArr', 'FbendIIIArr', 'FbendIVArr', ...
        'FbendArr', 'FcompArr', 'FentropyArr', ...
        'geomStar', 'i');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Local functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [F, geom] = FtotToy(aNeck, xi, p)
geom = [];
F = 1e30;
if aNeck <= 0 || xi <= 0 || ~isfinite(aNeck) || ~isfinite(xi)
    return;
end
[Ra, Rb, ~, ok] = solveRadii(aNeck, xi, p);
if ~ok
    return;
end
geom = computeToyQuantities(Ra, Rb, aNeck, xi, p);
if isempty(geom) || ~isfinite(geom.Ftot)
    F = 1e30;
    geom = [];
    return;
end
F = geom.Ftot;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [RaBest, RbBest, errBest, ok] = solveRadii(aNeck, xi, p)

ok = false;
RaBest = NaN;
RbBest = NaN;
errBest = Inf;

Rscale = sqrt(p.A0/(4*pi));

guessList = [
    max(aNeck + 1e-3, sqrt((1 - p.phi0)*p.A0/(4*pi))), max(aNeck + 1e-3, sqrt(p.phi0*p.A0/(4*pi)));
    max(aNeck + 1e-3, sqrt(p.A0/(8*pi))),              max(aNeck + 1e-3, sqrt(p.A0/(8*pi)));
    max(aNeck + 1e-3, 0.8*Rscale),                     max(aNeck + 1e-3, 0.3*Rscale);
    max(aNeck + 1e-3, 0.3*Rscale),                     max(aNeck + 1e-3, 0.8*Rscale);
    max(aNeck + 1e-3, 1.0*Rscale),                     max(aNeck + 1e-3, 0.2*Rscale);
    max(aNeck + 1e-3, 0.2*Rscale),                     max(aNeck + 1e-3, 1.0*Rscale)
];

opt = optimset('Display', 'off', ...
               'TolX', 1e-9, ...
               'TolFun', 1e-12, ...
               'MaxIter', 1200, ...
               'MaxFunEvals', 4000);

for k = 1:size(guessList,1)

    Ra0 = guessList(k,1);
    Rb0 = guessList(k,2);

    y0 = log([Ra0 - aNeck, Rb0 - aNeck]);

    [ySol, objVal] = fminsearch(@(y) conservationObjective(y, aNeck, xi, p), y0, opt);

    if objVal < errBest^2

        RaTry = aNeck + exp(ySol(1));
        RbTry = aNeck + exp(ySol(2));

        [r, isValid] = conservationResiduals(RaTry, RbTry, aNeck, xi, p);

        if isValid
            errTry = norm(r);

            if errTry < errBest
                errBest = errTry;
                RaBest = RaTry;
                RbBest = RbTry;
            end
        end
    end
end

if errBest < p.radiiTol
    ok = true;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function val = conservationObjective(y, aNeck, xi, p)

if any(~isfinite(y))
    val = 1e30;
    return;
end
d = exp(y);

if any(~isfinite(d)) || any(d <= 0)
    val = 1e30;
    return;
end

Ra = aNeck + d(1);
Rb = aNeck + d(2);
if Ra <= aNeck || Rb <= aNeck
    val = 1e30;
    return;
end
[r, isValid] = conservationResiduals(Ra, Rb, aNeck, xi, p);
if ~isValid
    val = 1e30;
    return;
end
val = sum(r.^2);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [r, isValid] = conservationResiduals(Ra, Rb, aNeck, xi, p)

isValid = false;
r = [1e10, 1e10];

try
    [Atotal, Aprotein] = computeConservedAmounts(Ra, Rb, aNeck, xi, p);
catch
    return;
end

if ~isfinite(Atotal) || ~isfinite(Aprotein) || Atotal <= 0 || Aprotein < 0
    return;
end

targetProtein = p.A0 * p.phi0;

rArea = (Atotal - p.A0) / p.A0;
rProt = (Aprotein - targetProtein) / targetProtein;

r = [rArea, rProt];
isValid = true;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Atotal, Aprotein] = computeConservedAmounts(Ra, Rb, aNeck, xi, p)

if Ra <= aNeck || Rb <= aNeck
    Atotal = NaN;
    Aprotein = NaN;
    return;
end

gammaA = gammaFromR(Ra, aNeck, p);
gammaB = gammaFromR(Rb, aNeck, p);

if ~isfinite(gammaA) || ~isfinite(gammaB)
    Atotal = NaN;
    Aprotein = NaN;
    return;
end

sA = sideGeometryGrid(Ra, aNeck, gammaA, -1, p.Ntheta);
sB = sideGeometryGrid(Rb, aNeck, gammaB, +1, p.Ntheta);

Lalpha = sA.ellAbs(end);
Lbeta  = sB.ellAbs(end);

l1 = -Lalpha;
l2 =  Lbeta;

wEff = effectiveWidth(xi, p);
Dw = tanh(l2/wEff) - tanh(l1/wEff);

if abs(Dw) < 1e-14 || ~isfinite(Dw)
    Atotal = NaN;
    Aprotein = NaN;
    return;
end

AI  = 2*pi*Ra^2;
AIV = 2*pi*Rb^2;

[AII, PII]   = sideAreaProteinFromGrid(sA, xi, l1, l2, Dw, p);
[AIII, PIII] = sideAreaProteinFromGrid(sB, xi, l1, l2, Dw, p);

Atotal = AI + AII + AIII + AIV;
Aprotein = PII + PIII + AIV;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [A, P] = sideAreaProteinFromGrid(s, xi, l1, l2, Dw, p)

psiVal = psiOfL(s.ellSigned, l1, l2, xi, Dw, p);

A = trapz(s.theta, s.dA_dtheta);
P = trapz(s.theta, psiVal .* s.dA_dtheta);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function C2H = targetC2H(aNeck, p)

% Linear form: C(a) = k*a + b.
C2H = p.C2H_k * aNeck + p.C2H_b;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function gamma = gammaFromR(R, aNeck, p)

gamma = NaN;

if R <= aNeck
    return;
end

C2H = targetC2H(aNeck, p);

% Because C(a) is used as target |2H(0)|, reject nonpositive values.
if C2H <= 0 || ~isfinite(C2H)
    return;
end

den = 1 - aNeck*C2H;

if den <= 0 || ~isfinite(den)
    return;
end

gamma2 = (R - aNeck) / ((R + aNeck) * den);

if gamma2 <= 0 || ~isfinite(gamma2)
    return;
end

gamma = sqrt(gamma2);

% Dynamic-flattening condition: 0 < gamma <= 1.
if isfield(p, 'enforceGammaMax') && p.enforceGammaMax
    tolGamma = 1e-10;
    if gamma > p.gammaMax + tolGamma
        gamma = NaN;
        return;
    end
    gamma = min(gamma, p.gammaMax);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function geom = computeToyQuantities(Ra, Rb, aNeck, xi, p)

geom = [];

if Ra <= aNeck || Rb <= aNeck
    return;
end

C2Htarget = targetC2H(aNeck, p);

gammaA = gammaFromR(Ra, aNeck, p);
gammaB = gammaFromR(Rb, aNeck, p);

if ~isfinite(gammaA) || ~isfinite(gammaB)
    return;
end

sA = sideGeometryGrid(Ra, aNeck, gammaA, -1, p.Ntheta);
sB = sideGeometryGrid(Rb, aNeck, gammaB, +1, p.Ntheta);

Lalpha = sA.ellAbs(end);
Lbeta  = sB.ellAbs(end);

l1 = -Lalpha;
l2 =  Lbeta;

wEff = effectiveWidth(xi, p);
Dw = tanh(l2/wEff) - tanh(l1/wEff);

if abs(Dw) < 1e-14 || ~isfinite(Dw)
    return;
end

AI  = 2*pi*Ra^2;
AIV = 2*pi*Rb^2;

[AII, PII, FbII, FcII, FeII]      = sideAllFromGrid(sA, xi, l1, l2, Dw, p.kappaAlpha, p);
[AIII, PIII, FbIII, FcIII, FeIII] = sideAllFromGrid(sB, xi, l1, l2, Dw, p.kappaBeta, p);

FbI  = 4*pi*p.kappaAlpha;
FbIV = 4*pi*p.kappaBeta;

Fbend    = FbI + FbII + FbIII + FbIV;
Fcomp    = FcII + FcIII;
Fentropy = FeII + FeIII;
Ftot     = Fbend + Fcomp + Fentropy;

Atotal   = AI + AII + AIII + AIV;
Aprotein = PII + PIII + AIV;

geom.Ra = Ra;
geom.Rb = Rb;
geom.aNeck = aNeck;
geom.xi = xi;
geom.wEff = wEff;
geom.xi0 = p.xi0;

geom.C2Htarget = C2Htarget;
geom.gammaAlpha = gammaA;
geom.gammaBeta  = gammaB;

geom.l1 = l1;
geom.l2 = l2;
geom.Dw = Dw;

geom.AI = AI;
geom.AII = AII;
geom.AIII = AIII;
geom.AIV = AIV;
geom.Atotal = Atotal;

geom.AproteinII = PII;
geom.AproteinIII = PIII;
geom.AproteinIV = AIV;
geom.Aprotein = Aprotein;

geom.FbendI = FbI;
geom.FbendII = FbII;
geom.FbendIII = FbIII;
geom.FbendIV = FbIV;

geom.Fbend = Fbend;
geom.Fcomp = Fcomp;
geom.Fentropy = Fentropy;
geom.Ftot = Ftot;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [A, P, Fb, Fc, Fe] = sideAllFromGrid(s, xi, l1, l2, Dw, kappa, p)

psiVal     = psiOfL(s.ellSigned, l1, l2, xi, Dw, p);
grad2Val   = grad2OfL(s.ellSigned, l1, l2, xi, Dw, p);
entropyVal = entropy_safe(psiVal);

A = trapz(s.theta, s.dA_dtheta);
P = trapz(s.theta, psiVal .* s.dA_dtheta);

Fb = trapz(s.theta, 0.5 .* kappa .* s.twoH.^2 .* s.dA_dtheta);

Fc = trapz(s.theta, ...
           (p.Lambda .* grad2Val + p.chi .* psiVal .* (1 - psiVal)) ...
           .* s.dA_dtheta ./ p.a0);

Fe = trapz(s.theta, ...
           p.kBT .* entropyVal .* s.dA_dtheta ./ p.a0);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function s = sideGeometryGrid(R, aNeck, gamma, sgn, N)

persistent cacheN theta sinTheta cosTheta cos2Theta

if isempty(cacheN) || cacheN ~= N
    u = linspace(0, 1, N);
    theta = (pi/2) .* u.^4;
    sinTheta = sin(theta);
    cosTheta = cos(theta);
    cos2Theta = cos(2.*theta);
    cacheN = N;
end
C = R^2 - aNeck^2;
rho = sqrt(aNeck^2 + C .* sinTheta.^2);
rho_t = C .* sinTheta .* cosTheta ./ rho;
rho_tt = C .* cos2Theta ./ rho - (rho_t.^2) ./ rho;
Aterm = rho + R*aNeck ./ rho;
Aterm_t = rho_t .* (1 - R*aNeck ./ rho.^2);
z_t = sgn .* gamma .* Aterm;
z_tt = sgn .* gamma .* Aterm_t;
ds_dtheta = sqrt(rho_t.^2 + z_t.^2);
dA_dtheta = 2*pi .* rho .* ds_dtheta;
ellAbs = cumtrapz(theta, ds_dtheta);
ellSigned = sgn .* ellAbs;
z = cumtrapz(theta, z_t);
k_mer = (rho_t .* z_tt - z_t .* rho_tt) ./ (ds_dtheta.^3);
k_azi = z_t ./ (rho .* ds_dtheta);
twoH = k_mer + k_azi;
s.theta = theta;
s.rho = rho;
s.z = z;
s.ellAbs = ellAbs;
s.ellSigned = ellSigned;
s.ds_dtheta = ds_dtheta;
s.dA_dtheta = dA_dtheta;
s.twoH = twoH;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function psi = psiOfL(ell, l1, l2, xi, Dw, p)
wEff = effectiveWidth(xi, p);
psi = zeros(size(ell));
idxLeft  = ell <= l1;
idxRight = ell >= l2;
idxMid   = ~idxLeft & ~idxRight;
psi(idxLeft)  = 0;
psi(idxRight) = 1;
psi(idxMid) = (tanh(ell(idxMid)./wEff) - tanh(l1./wEff)) ./ Dw;
psi = min(max(psi, 0), 1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function g2 = grad2OfL(ell, l1, l2, xi, Dw, p)
wEff = effectiveWidth(xi, p);
g2 = zeros(size(ell));
idxMid = (ell > l1) & (ell < l2);
g2(idxMid) = sech_local(ell(idxMid)./wEff).^4 ./ (wEff.^2 .* Dw.^2);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function wEff = effectiveWidth(xi, p)
wEff = sqrt(xi.^2 + p.xi0.^2);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = sech_local(x)
y = 1 ./ cosh(x);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = entropy_safe(phi)
epsPhi = 1e-12;
pp = min(max(phi, epsPhi), 1 - epsPhi);
y = pp .* log(pp) + (1 - pp) .* log(1 - pp);
end