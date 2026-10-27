function dco_small()
% DCO_SMALL  D_co at small Lambda~: continue the stable deformed state near onset from Lambda~ = 0.075
% (chi~ = 3.3, just above C_an) to Lambda~ = 0.065 ... 0.04 at fixed chi, then scan chi downward to the
% subcritical fold. Output dco_small.mat: per Lambda~, rec = [chi~, H, E~ (regularized), phi_min, phi_max].
kT = 0.02; ka = 0.4; ap = 0.01; delta = entropy_cutoff();
S = load('branch_Lambda0p075.mat'); if isfield(S,'Lambda'), Lam0 = S.Lambda; else, Lam0 = S.La; end
ok = find(abs(S.cchi/kT - 3.3) < 0.03 & S.hheight > 2.02 & S.hheight < 2.1);
[~, j] = min(abs(S.cchi(ok)/kT - 3.3)); k0 = ok(j);
sol = trim(S.ssol{k0}); chi = S.cchi(k0); Lam = Lam0;
fprintf('start: chi~ %.4f H %.4f\n', chi/kT, S.hheight(k0));
targets = [0.065 0.06 0.055 0.05 0.045 0.04]; out = struct('lam', {}, 'rec', {});
for t = targets
    Ls = exp(linspace(log(Lam), log(t*kT), 6)); Ls = Ls(2:end);
    for L = Ls, [sol, w] = solve(sol, L, chi); if w, error('path failed at %.4f', L/kT); end, end
    Lam = t*kT; rec = []; sl = sol; c = chi; dc = 0.005*kT;
    while dc > 2e-4*kT
        try
            [s2, w] = solve(sl, Lam, c - dc); if w, error('w'); end
            sl = s2; c = c - dc;
            [Er, ~] = energy(sl, Lam, c, kT, ka, ap, delta);
            rec(end+1,:) = [c/kT, sl.y(4,1)-sl.y(4,end), Er, min(sl.y(7,:)), max(sl.y(7,:))]; %#ok<AGROW>
            if max(sl.y(7,:)) - min(sl.y(7,:)) < 1e-3, break; end
        catch
            dc = dc/2;
        end
    end
    out(end+1) = struct('lam', t, 'rec', rec); %#ok<AGROW>
    fprintf('Lambda~ %.3f: %d points, chi~ down to %.4f\n', t, size(rec,1), min(rec(:,1)));
    save('dco_small.mat', 'out', '-v7');
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
