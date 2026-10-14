% RUN_LAMBDA_THEN_CHI  Obtain a branch of the full model at a new Lambda: starting
% from a stored solution at Lambda~ = 0.3, continue in Lambda at fixed chi to the
% target value (Lambda~ = 0.4 here), then trace the branch in chi as in RUN_CHI_CONTINUATION.
%
% Code units are those of RUN_CHI_CONTINUATION (see its header and the README):
% Lambda~ = Lambda/kT, chi~ = chi/kT, kT = 0.02.

clear; clc;
S = load(fullfile('data','branch_Lambda0p30.mat'));   % stored branch at Lambda~ = 0.3

% Geometry
R0 = 1;
area0 = 2*R0^2;          % A/(2*pi)

% Membrane
pre = 0;
ka = 0.4;

% Proteins
rho0 = 0.2;
kT = 0.02;
ap = 0.01;
C0 = 0;

% Forces not used in this work
fext = 0;
V0 = 0;

% Numerics
ri = 1e-5;
RTol = 1e-3;
Nmax = 50000;

% Continuation in Lambda at fixed chi. Natural-parameter continuation cannot pass a fold,
% so it is started from a prolate state (chi~ = 5), which exists over the whole range;
% an open-neck state exists only in a window of Lambda~ bounded by folds.
LambdaStart  = 0.006;    % Lambda~ = 0.3, the stored branch
LambdaTarget = 0.008;    % Lambda~ = 0.4
dLambda = 1e-4;          % largest step in Lambda (code units)
chiStart = 5.0;          % chi~ of the stored solution used as the starting point

% Continuation in chi afterwards
[~, kstart] = min(abs(S.cchi/kT - chiStart));
dir0 = 1;                % first chi stage: 1 = increasing chi, 2 = decreasing chi
nstage = 6;
increment = [0.001, -0.001];

sol = S.ssol{kstart};
chi = S.cchi(kstart);
solinit.parameters(1) = sol.parameters(1);
solinit.parameters(2) = sol.parameters(2);
solinit.parameters(3) = sol.parameters(3);
solinit.x = sol.x;
solinit.y = sol.y;

Lambda = LambdaStart;
yeq = @(u,y,para)shapef(u,y,para,ka,C0,kT,Lambda,chi,pre,fext,ap,V0);
ybc = @(ya,yb,para)twobcf(ya,yb,para,ka,C0,kT,ri,Lambda,chi,ap,area0,rho0);
opts = bvpset('RelTol',RTol,'AbsTol',1e-8,'NMax',Nmax);
sol = bvp5c(yeq,ybc,solinit,opts);

% Step Lambda from LambdaStart to LambdaTarget in n equal steps
n = ceil(abs(LambdaTarget - LambdaStart)/dLambda);
for j = 1:n
    Lambda = LambdaStart + (LambdaTarget - LambdaStart)*j/n;
    yeq = @(u,y,para)shapef(u,y,para,ka,C0,kT,Lambda,chi,pre,fext,ap,V0);
    ybc = @(ya,yb,para)twobcf(ya,yb,para,ka,C0,kT,ri,Lambda,chi,ap,area0,rho0);
    lastwarn('');
    sol = bvp5c(yeq,ybc,sol,opts);
    [warnMsg, warnId] = lastwarn;
    if contains(lower(warnId), 'bvp')
        error('Continuation in Lambda stopped at Lambda~ = %g (%s): %s', Lambda/kT, warnId, warnMsg);
    end
    fprintf('Lambda~ = %.4f\n', Lambda/kT);
end
Lambda = LambdaTarget;
para_end = sol.parameters;
sol_end = sol;

cchi = [];
hheight = [];
ssol = {};
flag = true;
for k = 1:nstage
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
