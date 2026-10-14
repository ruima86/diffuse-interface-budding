function [sol_end,para_end,cchi,hheight,ssol] = iterate_Hchi(sol_end,para_end,increment,cchi,hheight,...
    ka,pre,ri,area0,fext,rho0,C0,kT,ssol,Lambda,ap,RTol,Nmax,V0)
% ITERATE_HCHI  Continuation of a stationary branch in the pole-to-pole height H
% at fixed Lambda, with chi solved for as an unknown parameter (SHAPEF2, TWOBCF2).
%
% Used to pass the folds of a branch in chi, where continuation in chi
% (ITERATE_CHI) stops. H is stepped by increment(1) (> 0) or increment(2) (< 0),
% in the direction of the last two recorded heights. The run stops when bvp5c
% fails or warns. Each converged solution is appended to cchi, hheight and ssol.

if numel(hheight) < 2
    error('iterate_Hchi: needs at least two recorded points to set the direction.');
end
if abs(hheight(end)-hheight(end-1)) < 1e-6
    i = 3;                  % no progress in the previous run: stop
elseif hheight(end)-hheight(end-1) > 0
    i = 1;
else
    i = 2;
end

solinit = sol_end;
solinit.parameters(1) = para_end(1);
solinit.parameters(2) = para_end(2);
solinit.parameters(3) = para_end(3);
solinit.parameters(4) = cchi(end);
H = hheight(end);
yeq = @(u,y,para)shapef2(u,y,para,ka,C0,kT,Lambda,pre,fext,ap,V0);
ybc = @(ya,yb,para)twobcf2(ya,yb,para,ka,C0,kT,ri,Lambda,ap,area0,rho0,H);
opts = bvpset('RelTol',RTol,'AbsTol',1e-8,'NMax',Nmax);
sol = bvp5c(yeq,ybc,solinit,opts);

figure(1);
while true
    try
        if i == 3
            break;
        end
        yeq = @(u,y,para)shapef2(u,y,para,ka,C0,kT,Lambda,pre,fext,ap,V0);
        ybc = @(ya,yb,para)twobcf2(ya,yb,para,ka,C0,kT,ri,Lambda,ap,area0,rho0,H);
        opts = bvpset('RelTol',RTol,'AbsTol',1e-8,'NMax',Nmax);
        lastwarn('');
        sol = bvp5c(yeq,ybc,sol,opts);
        [warnMsg, warnId] = lastwarn;
        if contains(lower(warnId), 'bvp')       % stop on solver warnings only
            warning('bvp5c stopped at H = %g (%s): %s', H, warnId, warnMsg);
            break;
        end
    catch err
        fprintf('bvp5c failed at H = %g: %s\n', H, err.message);
        break;
    end
    fprintf('H = %.5f, chi = %.5f (chi~ = %.4f)\n', H, sol.parameters(4), sol.parameters(4)/kT);
    cchi = [cchi, sol.parameters(4)]; ssol = [ssol, sol]; hheight = [hheight, H];
    H = H + increment(i);
    para_end = sol.parameters;
    sol_end = sol;

    plot(sol.y(3,:), sol.y(4,:));
    axis equal;
    drawnow;
    if H < 0
        break;
    end
end
end
