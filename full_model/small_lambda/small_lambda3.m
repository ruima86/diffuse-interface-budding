function small_lambda3(s)
% SMALL_LAMBDA3  Prolate (s = 1) or open-neck (s = 2) branch near D_bud at Lambda~ = 0.06, 0.05, 0.04.
% The path in (Lambda, chi) avoids the folds of each branch:
%   prolate:   start at Lambda~ = 0.075, chi~ ~ 27 (well below its upper fold); lower Lambda at fixed
%              chi (the fold moves to larger chi), then scan chi upward to the fold at each target.
%   open neck: start at Lambda~ = 0.075, chi~ = 28.58; follow chi ~ Lambda^(-0.87), the scaling of the
%              lower fold, so that the state stays inside the stable open-neck window; at each target
%              scan chi in both directions.
kT = 0.02; ka = 0.4; ap = 0.01; delta = entropy_cutoff();
S = load('branch_Lambda0p075.mat'); if isfield(S,'Lambda'), Lam0 = S.Lambda; else, Lam0 = S.La; end
names = {'prolate','openneck'};
if s == 1
    [~, k0] = min(abs(S.cchi/kT - 27.0) + 10*(S.hheight > 2.8));   % prolate state near chi~ = 27
else
    k0 = 17;                                                       % open-neck state at chi~ = 28.58
end
sol = trim(S.ssol{k0}); chi0 = S.cchi(k0); Lam = Lam0; chi = chi0;
fprintf('%s start: Lambda~ %.3f chi~ %.3f H %.3f\n', names{s}, Lam/kT, chi/kT, S.hheight(k0));
[sol, ~] = solve(sol, Lam, chi);
for t = [0.06 0.05 0.04]
    LamT = t*kT; n = ceil(log(Lam/LamT)/0.005);
    Ls = exp(linspace(log(Lam), log(LamT), n+1)); Ls = Ls(2:end); ok = true;
    for L = Ls
        if s == 1, c = chi; else, c = chi0*(L/Lam0)^(-0.87); end
        try
            [sol, w] = solve(sol, L, c); if w, error('bvp5c warning'); end
        catch err
            fprintf('%s: path failed at Lambda~ = %.4f, chi~ = %.3f (%s)\n', names{s}, L/kT, c/kT, err.message); ok = false; break;
        end
        chi = c;
    end
    if ~ok, break; end
    Lam = LamT; fprintf('%s: reached Lambda~ = %.3f at chi~ = %.3f\n', names{s}, Lam/kT, chi/kT);
    tag = strrep(sprintf('%.2f', t), '.', 'p');
    save(sprintf('state3_%s_%s.mat', tag, names{s}), 'sol', 'Lam', 'chi');
    rec = []; [Er, Eb] = energy(sol, Lam, chi, kT, ka, ap, delta); rec(end+1,:) = [chi/kT, sol.y(4,1)-sol.y(4,end), Er, Eb, min(sol.y(7,:))];
    if s == 1, dirs = 1; nk = [400 0]; else, dirs = [1 -1]; nk = [80 40]; end
    for di = 1:numel(dirs)
        sl = sol; c = chi;
        for k = 1:nk(di)
            c = c + dirs(di)*0.0015*chi;
            try
                [sl, w] = solve(sl, Lam, c); if w, break; end
            catch
                break;
            end
            [Er, Eb] = energy(sl, Lam, c, kT, ka, ap, delta);
            rec(end+1,:) = [c/kT, sl.y(4,1)-sl.y(4,end), Er, Eb, min(sl.y(7,:))]; %#ok<AGROW>
        end
    end
    save(sprintf('small_lambda_%s_%s.mat', tag, names{s}), 'rec', 't', '-v7');
    fprintf('   %s Lambda~=%.2f: %d points, chi~ %.3f..%.3f\n', names{s}, t, size(rec,1), min(rec(:,1)), max(rec(:,1)));
    if s == 1, chi = chi; end   % prolate continues the Lambda path from the state reached at fixed chi
end
end

function sol = trim(sol)
if isfield(sol,'parameters'), p = sol.parameters(1:3); sol = rmfield(sol,'parameters'); sol.parameters = p; end
end

function [sol, warned] = solve(sol, Lam, chi)
kT = 0.02; ka = 0.4; ap = 0.01; rho0 = 0.2; C0 = 0; pre = 0; fext = 0; V0 = 0; ri = 1e-5; area0 = 2;
opts = bvpset('RelTol',1e-3,'AbsTol',1e-8,'NMax',200000);
lastwarn('');
sol = bvp5c(@(u,y,p)shapef(u,y,p,ka,C0,kT,Lam,chi,pre,fext,ap,V0), ...
            @(ya,yb,p)twobcf(ya,yb,p,ka,C0,kT,ri,Lam,chi,ap,area0,rho0), sol, opts);
[~, id] = lastwarn; warned = contains(lower(id), 'bvp');
end

function [Er, Eb] = energy(sol, Lam, chi, kT, kappa, ap, delta)
u = linspace(sol.x(1), sol.x(end), 10000); y = deval(sol, u); h = sol.parameters(1);
psi=y(1,:); dpsi=y(2,:); r=y(3,:); phi=y(7,:); dphi=y(8,:);
twoH = dpsi/h + sin(psi)./r;
Eb = 2*pi*trapz(u, h*r.*0.5*kappa.*twoH.^2)/(8*pi*kappa);
q = min(phi,1-phi); s = zeros(size(q)); mid = q >= delta;
s(mid) = q(mid).*log(q(mid)) + (1-q(mid)).*log1p(-q(mid));
d = q(~mid)-delta; h0 = delta*log(delta)+(1-delta)*log(1-delta); h1 = log(delta/(1-delta));
h2 = 1/delta+1/(1-delta); h3 = -1/delta^2+1/(1-delta)^2; h4 = 2/delta^3+2/(1-delta)^3;
s(~mid) = h0+h1*d+h2*d.^2/2+h3*d.^3/6+h4*d.^4/24;
Er = Eb + 2*pi*trapz(u, h*r.*(Lam*(dphi/h).^2 + chi*phi.*(1-phi) + kT*s)/ap)/(8*pi*kappa);
end
