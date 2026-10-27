function gamma_fold_branch(tags, lf, hf)
% GAMMA_FOLD_BRANCH  Line tension on the vesicle at the open-neck fold of the stored branches,
% gamma = (2 Lambda/a_p) int (d phi/ds)^2 ds (as in Fig. 1(c)-(e)), in units kBT/R0 and kappa/R0.
% Also: neck radius, interface width, planar values are computed separately (line_tension.py).
kT = 0.02; ap = 0.01; kap = 0.4;
res = nan(numel(tags), 7);
for i = 1:numel(tags)
    S = load(['branch_Lambda' tags{i} '.mat']); if isfield(S,'Lambda'), Lam = S.Lambda; else, Lam = S.La; end
        [~, k] = min(abs(S.cchi/kT - lf(i))/0.01 + abs(S.hheight - hf(i))/0.001);   % stored fold point (chi~, H)
    [g, rn, w] = linetension(S.ssol{k}, Lam, ap);
    res(i,:) = [Lam/kT, S.cchi(k)/kT, S.hheight(k), rn, w, g/kT, g/kap];
    fprintf('Lambda~ %.4f  fold chi~ %.4f  H %.4f  r_n %.4f  w %.4f  gamma %.2f kT/R0 = %.3f kappa/R0\n', res(i,:));
    clear S
end
save(sprintf('gamma_fold_branch_%s.mat', [tags{:}]), 'res', '-v7');
end

function [g, rn, w] = linetension(sol, Lam, ap)
u = linspace(sol.x(1), sol.x(end), 20000); y = deval(sol, u); h = sol.parameters(1);
dphi = y(8,:)/h; g = trapz(u*h, 2*Lam/ap*dphi.^2);
r = y(3,:); m = r(islocalmin(r)); rn = NaN; if ~isempty(m), rn = min(m); end
w = (max(y(7,:)) - min(y(7,:)))/max(abs(dphi));
end
