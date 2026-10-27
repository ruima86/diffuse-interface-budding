function small_lambda4(s)
% SMALL_LAMBDA4  Extend the Lambda~ = 0.04 scans of small_lambda3 so the two branches overlap and
% the open neck crosses Division.
%   s = 1: prolate, from state3_0p04_prolate (chi~ = 27.0), chi upward until bvp5c fails (upper fold).
%   s = 2: open neck, from state3_0p04_openneck (chi~ = 49.38), chi upward to ~58.
% Output: small_lambda4_0p04_<name>.mat, rec columns as in small_lambda3.
kT = 0.02; ka = 0.4; ap = 0.01; delta = entropy_cutoff();
names = {'prolate','openneck'};
S = load(sprintf('state3_0p04_%s.mat', names{s})); sol = S.sol; Lam = S.Lam; chi = S.chi;
if s == 1, step = 0.0015*chi; nk = 900; else, step = 0.0025*chi; nk = 70; end
fprintf('%s start: Lambda~ %.3f chi~ %.3f\n', names{s}, Lam/kT, chi/kT);
rec = []; c = chi;
for k = 1:nk
    c = c + step;
    try
        [sol, w] = solve(sol, Lam, c); if w, fprintf('stop at chi~ %.3f: %s\n', c/kT, lastwarn); break; end
    catch err
        fprintf('stop at chi~ %.3f: %s\n', c/kT, err.message); break;
    end
    [Er, Eb] = energy(sol, Lam, c, kT, ka, ap, delta);
    rec(end+1,:) = [c/kT, sol.y(4,1)-sol.y(4,end), Er, Eb, min(sol.y(7,:))]; %#ok<AGROW>
    if mod(k, 20) == 0, save(sprintf('small_lambda4_0p04_%s.mat', names{s}), 'rec', '-v7'); end
end
save(sprintf('small_lambda4_0p04_%s.mat', names{s}), 'rec', '-v7');
fprintf('   %s: %d points, chi~ %.3f..%.3f\n', names{s}, size(rec,1), min(rec(:,1)), max(rec(:,1)));
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
