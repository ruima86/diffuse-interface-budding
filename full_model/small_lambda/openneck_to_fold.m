function openneck_to_fold()
% OPENNECK_TO_FOLD  Lambda~ = 0.075: the stable open-neck branch between the end of the stored branch
% (chi~ = 35.39, H = 2.926) and the fold (chi~ = 43.52, H = 2.789), and a short unstable piece past the fold.
% Starts from the fold solution saved by lfold_track_gamma (out(2), Lambda~ = 0.075) and steps the height H
% (chi solved as unknown, shapef2/twobcf2). Output openneck_0p075.mat:
% rec = [chi~, H, E~ (regularized), r_n, w, gamma (kT/R0)] and the meridian (r, z, phi) at 160 points.
kT = 0.02; ka = 0.4; ap = 0.01; Lam = 0.075*kT; delta = entropy_cutoff();
S = load('lfold_track_gamma.mat'); i = find(abs([S.out.lam] - 0.075) < 1e-9); sol0 = S.out(i).sol;
H0 = sol0.y(4,1) - sol0.y(4,end);
rec = []; mer = [];
for dirn = [1 -1]                                   % up in H (stable side), then down (unstable side)
    s = sol0; H = H0; n = 0; nmax = 80; if dirn < 0, nmax = 12; end
    while n < nmax
        try
            [s2, w] = solveH(s, Lam, H); if w, error('w'); end
        catch
            fprintf('stop at H %.4f\n', H); break;
        end
        s = s2; [row, m] = observables(s, Lam, kT, ka, ap, delta);
        rec(end+1,:) = row; mer(end+1,:,:) = m; %#ok<AGROW>
        fprintf('H %.4f  chi~ %.4f  r_n %.4f  w %.4f  gamma %.2f\n', H, row(1), row(4), row(5), row(6));
        if dirn > 0 && row(1) < 35.3, break; end     % overlaps the stored branch
        H = H + dirn*0.002; n = n + 1;
    end
    save('openneck_0p075.mat', 'rec', 'mer', '-v7');
end
end

function [row, m] = observables(sol, Lam, kT, kappa, ap, delta)
u = linspace(sol.x(1), sol.x(end), 10000); y = deval(sol, u); h = sol.parameters(1);
psi = y(1,:); dpsi = y(2,:); r = y(3,:); z = y(4,:); phi = y(7,:); dphi = y(8,:); chi = sol.parameters(4);
twoH = dpsi/h + sin(psi)./r;
Eb = 2*pi*trapz(u, h*r.*0.5*kappa.*twoH.^2)/(8*pi*kappa);
q = min(phi, 1-phi); s = zeros(size(q)); mid = q >= delta;
s(mid) = q(mid).*log(q(mid)) + (1-q(mid)).*log1p(-q(mid));
d = q(~mid) - delta; h0 = delta*log(delta) + (1-delta)*log(1-delta); h1 = log(delta/(1-delta));
h2 = 1/delta + 1/(1-delta); h3 = -1/delta^2 + 1/(1-delta)^2; h4 = 2/delta^3 + 2/(1-delta)^3;
s(~mid) = h0 + h1*d + h2*d.^2/2 + h3*d.^3/6 + h4*d.^4/24;
Er = Eb + 2*pi*trapz(u, h*r.*(Lam*(dphi/h).^2 + chi*phi.*(1-phi) + kT*s)/ap)/(8*pi*kappa);
mm = r(islocalmin(r)); rn = NaN; if ~isempty(mm), rn = min(mm); end
w = (max(phi) - min(phi))*h/max(abs(dphi));
g = trapz(u, 2*Lam/ap*(dphi/h).^2*h)/kT;
row = [chi/kT, z(1) - z(end), Er, rn, w, g];
k = round(linspace(1, 10000, 160)); m = [r(k); z(k); phi(k)];
end

function [sol, warned] = solveH(sol, Lam, H)
kT = 0.02; ka = 0.4; ap = 0.01; rho0 = 0.2; C0 = 0; pre = 0; fext = 0; V0 = 0; ri = 1e-5; area0 = 2;
opts = bvpset('RelTol',1e-3,'AbsTol',1e-8,'NMax',200000);
lastwarn('');
sol = bvp5c(@(u,y,p)shapef2(u,y,p,ka,C0,kT,Lam,pre,fext,ap,V0), ...
            @(ya,yb,p)twobcf2(ya,yb,p,ka,C0,kT,ri,Lam,ap,area0,rho0,H), sol, opts);
[~, id] = lastwarn; warned = contains(lower(id), 'bvp');
end
