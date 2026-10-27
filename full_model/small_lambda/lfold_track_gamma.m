function lfold_track_gamma()
% LFOLD_TRACK_GAMMA  L_fold at small Lambda~ (Fig. 3, Lambda~ = 0.04-0.08). Starts from the converged fold
% solution of branch_Lambda0p0875.mat (r_n = 0.014, ~3000 mesh points). Lambda is lowered in 1% steps at
% fixed height H (chi solved as unknown, shapef2/twobcf2), slightly above the fold height so that the path
% stays on the stable open-neck side. At each target Lambda~, H is then stepped down past the maximum of
% chi~, the fold. Also records the line tension on the vesicle, gamma = (2 Lambda/a_p) int (dphi/ds)^2 ds
% (kappa/R0), and saves the solution at each fold.
% Output lfold_track_gamma.mat: out(i).lam, out(i).rec = [chi~, H, r_n, mesh, gamma] around the fold,
% out(i).fold, out(i).gamma, out(i).sol.
kT = 0.02; ka = 0.4; ap = 0.01; rho0 = 0.2; C0 = 0; pre = 0; fext = 0; V0 = 0; ri = 1e-5; area0 = 2;
S = load('branch_Lambda0p0875.mat'); on = find(S.hheight > 2.7); [~, j] = max(S.cchi(on)); k = on(j);
% start a little before the fold, on the stable side (larger H)
cand = on(S.hheight(on) > S.hheight(k) + 0.004 & S.hheight(on) < S.hheight(k) + 0.006 & S.cchi(on) > S.cchi(k) - 0.2*kT);
k0 = cand(1); sol = S.ssol{k0}; p = sol.parameters(1:3); sol = rmfield(sol, 'parameters');
sol.parameters = [p(:)', S.cchi(k0)]; H = S.hheight(k0); Lam = 0.0875*kT;
fprintf('start: chi~ %.4f H %.5f (fold chi~ %.4f H %.5f)\n', S.cchi(k0)/kT, H, S.cchi(k)/kT, S.hheight(k));
clear S
out = struct('lam', {}, 'rec', {}, 'fold', {}, 'gamma', {}, 'sol', {});
for t = [0.08 0.075 0.07 0.065 0.06 0.055 0.05 0.045 0.04]
    while Lam > t*kT*(1 + 1e-9)
        Ln = max(Lam*0.99, t*kT); done = false; dl = Lam - Ln;
        while ~done
            try
                [s2, w] = solveH(sol, Ln, H); if w, error('w'); end
                sol = s2; Lam = Ln; done = true;
            catch
                dl = dl/2; Ln = Lam - dl;
                if dl < 1e-5*Lam
                    % no solution at this H: move H up (away from the fold) and retry
                    H = H + 0.002; fprintf('  raise H to %.4f at Lambda~ %.4f\n', H, Lam/kT); dl = Lam*0.01; Ln = Lam - dl;
                end
            end
        end
    end
    % step H down through the fold
    rec = []; s = sol; h = H; dh = 0.0005; nfail = 0;
    while true
        try
            [s2, w] = solveH(s, Lam, h); if w, error('w'); end
            s = s2; r = s.y(3,:); m = r(islocalmin(r)); rn = NaN; if ~isempty(m), rn = min(m); end
            u = linspace(s.x(1), s.x(end), 20000); yy = deval(s, u); hh = s.parameters(1);
            gam = trapz(u*hh, 2*Lam/ap*(yy(8,:)/hh).^2)/0.4;
            rec(end+1,:) = [s.parameters(4)/kT, h, rn, numel(s.x), gam]; %#ok<AGROW>
            if isempty(rec(1:end-1,1)) || rec(end,1) >= max(rec(1:end-1,1)), sfold = s; end
            if size(rec,1) > 6 && rec(end,1) < max(rec(:,1)) - 0.02, break; end   % clearly past the maximum
            h = h - dh;
        catch
            nfail = nfail + 1; dh = dh/2; h = rec(end,2) - dh;
            if nfail > 4, break; end
        end
    end
    [cm, im] = max(rec(:,1)); f = [cm, rec(im,2), rec(im,3)];
    if im > 1 && im < size(rec,1)
        q = polyfit(rec(im-1:im+1,2), rec(im-1:im+1,1), 2); hf = -q(2)/(2*q(1)); f = [polyval(q, hf), hf, rec(im,3)];
    end
    gf = interp1(rec(:,2), rec(:,5), f(2), 'pchip');
    out(end+1) = struct('lam', t, 'rec', rec, 'fold', f, 'gamma', gf, 'sol', sfold); %#ok<AGROW>
    fprintf('Lambda~ %.3f: fold chi~ %.4f at H %.4f, r_n %.4f, gamma %.3f kappa/R0, mesh %d\n', t, f(1), f(2), f(3), gf, rec(im,4));
    save('lfold_track_gamma.mat', 'out', '-v7');
end
end

function [sol, warned] = solveH(sol, Lam, H)
kT = 0.02; ka = 0.4; ap = 0.01; rho0 = 0.2; C0 = 0; pre = 0; fext = 0; V0 = 0; ri = 1e-5; area0 = 2;
opts = bvpset('RelTol',1e-3,'AbsTol',1e-8,'NMax',200000);
lastwarn('');
sol = bvp5c(@(u,y,p)shapef2(u,y,p,ka,C0,kT,Lam,pre,fext,ap,V0), ...
            @(ya,yb,p)twobcf2(ya,yb,p,ka,C0,kT,ri,Lam,ap,area0,rho0,H), sol, opts);
[~, id] = lastwarn; warned = contains(lower(id), 'bvp');
end
