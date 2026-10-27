function width_fold_track()
% WIDTH_FOLD_TRACK  Interface width, neck radius and line tension at the tracked folds
% (lfold_track_gamma.mat), with the same definitions as gamma_fold_branch.m.
ap = 0.01; kap = 0.4; kT = 0.02;
S = load('lfold_track_gamma.mat'); res = nan(numel(S.out), 6);
for i = 1:numel(S.out)
    sol = S.out(i).sol; Lam = S.out(i).lam*kT;
    u = linspace(sol.x(1), sol.x(end), 20000); y = deval(sol, u); h = sol.parameters(1);
    dphi = y(8,:)/h; g = trapz(u*h, 2*Lam/ap*dphi.^2)/kap;
    r = y(3,:); m = r(islocalmin(r)); rn = min(m);
    w = (max(y(7,:)) - min(y(7,:)))/max(abs(dphi));
    res(i,:) = [S.out(i).lam, sol.parameters(4)/kT, rn, w, rn/w, g];
    fprintf('Lambda~ %.3f  chi~ %.4f  r_n %.4f  w %.4f  r_n/w %.3f  gamma %.3f kappa/R0\n', res(i,:));
end
save('width_fold_track.mat', 'res', '-v7');
end
