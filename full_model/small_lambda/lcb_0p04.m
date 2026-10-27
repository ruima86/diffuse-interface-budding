function lcb_0p04()
% LCB_0P04  Continue the stable open-neck state at Lambda~ = 0.04 from chi~ = 49.38 (state3_0p04_openneck)
% upward past its crossing with the complete bud (expected near chi~ = 56.1), in steps of 0.074 in chi~
% (as in small_lambda3, which reached 55.31) with NMax = 6e5 and step halving on failure.
% Output lcb_0p04.mat: rec = [chi~, H, E~ (regularized), E~_bend, phi_min], saved after every step.
kT = 0.02; ka = 0.4; ap = 0.01; delta = entropy_cutoff();
S = load('state3_0p04_openneck.mat'); sol = S.sol; chi = S.chi; Lam = S.Lam;
fprintf('start Lambda~ %.3f chi~ %.4f\n', Lam/kT, chi/kT);
dmax = 0.0015*chi; dc = dmax; rec = [];
while chi/kT < 57.5 && dc > dmax/64
    try
        [s2, w] = solve(sol, Lam, chi + dc); if w, error('w'); end
        sol = s2; chi = chi + dc;
        [Er, Eb] = energy(sol, Lam, chi, kT, ka, ap, delta);
        rec(end+1,:) = [chi/kT, sol.y(4,1)-sol.y(4,end), Er, Eb, min(sol.y(7,:))]; %#ok<AGROW>
        fprintf('chi~ %.4f  E~ %.8f  mesh %d\n', chi/kT, Er, numel(sol.x));
        save('lcb_0p04.mat', 'rec', '-v7');
        dc = min(2*dc, dmax);
    catch
        dc = dc/2; fprintf('  halve step to %.4f\n', dc/kT);
    end
end
fprintf('end at chi~ %.4f\n', chi/kT);
end
function sol = trim(sol)
if isfield(sol,'parameters'), p = sol.parameters(1:3); sol = rmfield(sol,'parameters'); sol.parameters = p; end
end

function [sol, warned] = solve(sol, Lam, chi)
kT = 0.02; ka = 0.4; ap = 0.01; rho0 = 0.2; C0 = 0; pre = 0; fext = 0; V0 = 0; ri = 1e-5; area0 = 2;
opts = bvpset('RelTol',1e-3,'AbsTol',1e-8,'NMax',600000);
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
