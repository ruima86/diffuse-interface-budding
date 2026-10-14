function [sol_end,para_end,cchi,hheight,ssol] = iterate_chi(sol_end,para_end,increment,cchi,hheight,...
    ka,pre,ri,area0,fext,rho0,C0,kT,ssol,Lambda,ap,RTol,Nmax,V0,chi0,dir0)
% ITERATE_CHI  Natural-parameter continuation of a stationary branch in chi at fixed Lambda.
%
% chi is stepped by increment(1) (> 0) or increment(2) (< 0) and the boundary-value
% problem is re-solved with bvp5c, each time starting from the previous solution.
% The run stops when bvp5c fails or issues a warning, which happens close to a fold
% of the branch in chi; ITERATE_HCHI then takes over.
%
% Direction: if cchi holds at least two points, the run continues in the direction of
% the last two; otherwise it starts at chi0 in direction dir0 (1: increasing chi,
% 2: decreasing chi; default 2).
% Each converged solution is appended to cchi (chi), hheight (pole-to-pole height
% z(0) - z(1)) and ssol (bvp5c solution structure).

if nargin < 21 || isempty(dir0)
    dir0 = 2;
end
if isempty(cchi)
    if nargin < 20 || isempty(chi0)
        error('iterate_chi: give the starting chi0 when cchi is empty.');
    end
    chi = chi0;
    i = dir0;
elseif numel(cchi) < 2
    chi = cchi(end);
    i = dir0;
else
    chi = cchi(end);
    d = cchi(end) - cchi(end-1);
    if abs(d) < 1e-7
        i = 3;              % no progress in the previous run: stop
    elseif d > 0
        i = 1;
    else
        i = 2;
    end
end

if isfield(sol_end, 'parameters')
    sol_end = rmfield(sol_end, 'parameters');
end
solinit = sol_end;
solinit.parameters(1) = para_end(1);
solinit.parameters(2) = para_end(2);
solinit.parameters(3) = para_end(3);
yeq = @(u,y,para)shapef(u,y,para,ka,C0,kT,Lambda,chi,pre,fext,ap,V0);
ybc = @(ya,yb,para)twobcf(ya,yb,para,ka,C0,kT,ri,Lambda,chi,ap,area0,rho0);
opts = bvpset('RelTol',RTol,'AbsTol',1e-8,'NMax',Nmax);
sol = bvp5c(yeq,ybc,solinit,opts);

figure(1);
while true
    try
        if i == 3
            break;
        end
        yeq = @(u,y,para)shapef(u,y,para,ka,C0,kT,Lambda,chi,pre,fext,ap,V0);
        ybc = @(ya,yb,para)twobcf(ya,yb,para,ka,C0,kT,ri,Lambda,chi,ap,area0,rho0);
        opts = bvpset('RelTol',RTol,'AbsTol',1e-8,'NMax',Nmax);
        lastwarn('');
        sol = bvp5c(yeq,ybc,sol,opts);
        [warnMsg, warnId] = lastwarn;
        if contains(lower(warnId), 'bvp')       % stop on solver warnings only
            warning('bvp5c stopped at chi = %g (%s): %s', chi, warnId, warnMsg);
            break;
        end
    catch err
        fprintf('bvp5c failed at chi = %g: %s\n', chi, err.message);
        break;
    end
    fprintf('chi = %.5f (chi~ = %.4f)\n', chi, chi/kT);
    cchi = [cchi, chi]; hheight = [hheight, sol.y(4,1)-sol.y(4,end)]; ssol = [ssol, sol];
    chi = chi + increment(i);
    para_end = sol.parameters;
    sol_end = sol;

    plot(sol.y(3,:), sol.y(4,:));
    axis equal;
    drawnow;
    if chi < 0
        break;
    end
end

end
