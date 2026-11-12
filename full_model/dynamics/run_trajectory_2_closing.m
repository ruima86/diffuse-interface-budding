% Trajectory 2: unstable state, closing direction.
% Start from the original steady state; fixed dt, no continuation.
clear;
clc;
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);

% Geometry of the sphere
R0 = 1;  
area0 = 2*R0^2;

% Membrane properties
pre = 0;     % pressure difference
ka = 0.4;

% Protein properties
C0 = 0;
rho0 = 0.2;
Di = 0.00312;
kT = 0.02;
ap=0.01;
areaf = 0.01*area0;

% external force
fext = 0;

% phase separation strength
Lambda=0.008;
chi = 0.18;
V0=0;
% initial guess of the pressure
h0 = R0*pi;

% friction coefficients
Ga1 = 0;       % normal friction
Ga2 = 0.01;       % tangential friction
% viscosity
etas = 1;

% Computational parameters
ri = 1e-5;
dt = 0.1;
T =50000;
RTol = 1e-3;
Nbin = 201;
Nmax = 5000;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Initial condition: original steady state (not a restart)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load(fullfile(scriptDir,'sol_La_0.008_chi_0.18_2.mat'));
ysol = sol.y;
yp = sol.idata.yp;
ymid = sol.idata.ymid;

sol.y(7,:) = zeros(1,length(sol.x));
sol.y(8,:) = sol.parameters(2)*ones(1,length(sol.x));
sol.y(9,:) = zeros(1,length(sol.x));
sol.y(10,:) = zeros(1,length(sol.x));
sol.y(11,:) = zeros(1,length(sol.x));
sol.y(12,:) = ysol(7,:);
sol.y(13,:) = ysol(8,:);
sol.y(14,:) = sol.parameters(3)*ones(1,length(sol.x));
% sol.y(14,:) = zeros(1,length(sol.x));
sol.y(15,:) = ysol(9,:);


sol.idata.ymid(7,:) = zeros(1,length(sol.x)-1);
sol.idata.ymid(8,:) = sol.parameters(2)*ones(1,length(sol.x)-1);
sol.idata.ymid(9,:) = zeros(1,length(sol.x)-1);
sol.idata.ymid(10,:) = zeros(1,length(sol.x)-1);
sol.idata.ymid(11,:) = zeros(1,length(sol.x)-1);
sol.idata.ymid(12,:) = ymid(7,:);
sol.idata.ymid(13,:) = ymid(8,:);
sol.idata.ymid(14,:) = sol.parameters(3)*ones(1,length(sol.x)-1);
% sol.idata.ymid(14,:) = zeros(1,length(sol.x)-1);
sol.idata.ymid(15,:) = ymid (9,:);

sol.idata.yp(7,:) = zeros(1,length(sol.x));
sol.idata.yp(8,:) = zeros(1,length(sol.x));
sol.idata.yp(9,:) = zeros(1,length(sol.x));
sol.idata.yp(10,:) = zeros(1,length(sol.x));
sol.idata.yp(11,:) = zeros(1,length(sol.x));
sol.idata.yp(12,:) = yp(7,:);
sol.idata.yp(13,:) = yp(8,:);
sol.idata.yp(14,:) = zeros(1,length(sol.x));
sol.idata.yp(15,:) = yp(9,:);


sol.parameters = sol.parameters(1);

% Set either target to zero to run the other perturbation on its own.
neckChange = -0.05;       % Relative change in the initial neck radius
proteinEnergy = 0;   % Protein perturbation energy in kBT at fixed shape
[sol,Inte] = perturbInitialState(sol,neckChange,proteinEnergy,Lambda,chi,kT,ap);

outputDir = fullfile(scriptDir,'data');
if ~isfolder(outputDir)
    mkdir(outputDir);
end


filename = [fullfile(outputDir,'trajectory_2_closing'),...
            '_chi_',num2str(chi),...
             '_La_',num2str(Lambda),...
             '_V0_',num2str(V0),...
            '_I_',num2str(Inte),...           
            '_fext_',num2str(fext),...
            '_p_',num2str(pre),...
            '_etas_',num2str(etas),...
            '_Ga1_',num2str(Ga1),...
            '_Ga2_',num2str(Ga2),...
            '_Di_',num2str(Di),...
            '_kT_',num2str(kT),...
            '_rho0_',num2str(rho0),...                      
            '_dt_',num2str(dt),...
            '_RTol_',num2str(RTol),...
            '_date_',date,'.mat'];


sol_all = cell(1,4);
Tr = dt:dt*1:T;
j = 0;
for t = Tr
    t
    % fext0 = fext*(1-exp(-0.1*t));
    uc = sol.x;
    y0 = deval(sol,uc);
    r0c = y0(3,:);
    z0c = y0(4,:);
    psi0c = y0(1,:);
    rho0c = y0(12,:);
    h0c = sol.parameters;
    yeq = @(u,y,para)shapef(u,y,para,ka,C0,Ga1,Ga2,Di,kT,etas,Lambda,chi,pre,fext,area0,areaf,ap,uc,r0c,z0c,psi0c,rho0c,h0c,dt,V0);
    % jach = @(u,y,para)jac(u,y,para,ka,C0,Ga1,Ga2,Di,kT,etas,Lambda,chi,pre,fext,area0,areaf,ap,uc,r0c,z0c,psi0c,rho0c,h0c,dt);
    ybc = @(ya,yb,para)twobcf(ya,yb,para,ka,C0,kT,ri,Lambda,chi,ap);
    % opts = bvpset('RelTol',RTol,'AbsTol',1e-8,'NMax',Nmax,'FJacobian',jach);
    opts = bvpset('RelTol',RTol,'AbsTol',1e-8,'NMax',Nmax);
    sol0 = sol;
    sol = bvp5c(yeq,ybc,sol,opts);
    uc = sol.x;
    yc = deval(sol,uc);
    rc = yc(3,:);
    zc = yc(4,:);
    psic = yc(1,:);
    rhoc = yc(12,:);
    hc = sol.parameters;
    y0c = deval(sol0,uc);
    r0c = y0c(3,:);
    z0c = y0c(4,:);
    dr = (rc - r0c)/dt;
    dz = (zc - z0c)/dt;
    psi = psic;
    vn = cos(psi).*dz + sin(psi).*dr;
    trapz(uc,rc.*hc.*rhoc)
    
    % The following code is to check whether dF + dW < 0
    dpsic = yc(2,:);
    drhoc = yc(13,:);
    vuc = yc(7,:);
    areac = yc(15,:);

    psi0c = y0c(1,:);
    dpsi0c = y0c(2,:);
    rho0c = y0c(12,:);
    drho0c = y0c(13,:);
    h0c = sol0.parameters;
    % calculate the free energy of current step
    f1 =(1/2).*hc.^(-2).*ka.*rc.^(-2).*(dpsic.*rc+hc.*((-1).*C0.*rc.*rhoc+sin( ...
        psic))).^2;
    f2 =ap.^(-1).*kT.*((1+(-1).*rhoc).*log(1+(-1).*rhoc)+rhoc.*log(rhoc));
    f3 =ap.^(-1).*chi.*(1+(-1).*rhoc).*rhoc;
    f4 =ap.^(-1).*drhoc.^2.*hc.^(-2).*Lambda;
    f5 =-fext.*sin(psic)./rc;
    % calculate the free energy of previous step
    f10 =(1/2).*h0c.^(-2).*ka.*r0c.^(-2).*(dpsi0c.*r0c+h0c.*((-1).*C0.*r0c.*rho0c+sin( ...
        psi0c))).^2;
    f20 =ap.^(-1).*kT.*((1+(-1).*rho0c).*log(1+(-1).*rho0c)+rho0c.*log(rho0c));
    f30 =ap.^(-1).*chi.*(1+(-1).*rho0c).*rho0c;
    f40 =ap.^(-1).*drho0c.^2.*h0c.^(-2).*Lambda;
    f50 = -fext.*sin(psi0c)./r0c;
    En = trapz(uc,(f1+f2+f3+f4+f5).*rc.*hc);
    En0 = trapz(uc,(f10+f20+f30+f40+f50).*r0c.*h0c); 
    % calculate the free energy difference
    dF = (En - En0)/dt;


    %fz = fext/2*(1-tanh(100*(areac-areaf)))-fext/2*(1-tanh(100*(area0-areaf-areac)));                                                        % fz(area)
    fz = 0;
    fr = 0;                                                        % fr(area)
    fu = fr.*cos(psic)-fz.*sin(psic);
    fn = pre + fr.*sin(psic)+fz.*cos(psic);                                % fn0 = pressure
    fphi = 0;
    % calculate the word done by external forces
    dW = trapz(uc,(fu.*vuc+fn.*vn).*rc.*hc);

    dF - dW   % this quantity should be smaller than zero
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if sol.stats.maxerr > RTol
        break;
    else
        if min(abs(Tr-t))<0.5*dt
            j = j+1;
            sol_all{j,1} = t;
            sol_all{j,2} = sol;
            sol_all{j,3} = vn;
            sol_all{j,4} = dF - dW;
            save(filename,'sol_all')
        end
    end
end

function [sol,amplitude] = perturbInitialState(sol,neckChange,proteinEnergy,Lambda,chi,kT,ap)
u = linspace(sol.x(1),sol.x(end),10001);
original = deval(sol,u);
h = sol.parameters;
r = original(3,:);
phi = original(12,:);
weight = h*r;
area = cumtrapz(u,weight);
proteinTotal = trapz(u,weight.*phi);
[rn,neckIndex] = neckMinimum(u,r);
neck = u(neckIndex);
v = u;

% Shape perturbation: displace the neck and compensate area on the far lobe.
if neckChange ~= 0
    neckWidth = 0.12;
    [neckBump,neckSlope] = compactBump(u,neck,neckWidth);
    [areaBump,areaSlope] = compactBump(u,0.72,0.18);
    neckMask = abs(u-neck) < neckWidth;
    targetRadius = rn*(1+neckChange);
    radiusError = @(a) min(r(neckMask)+a*neckBump(neckMask))-targetRadius;
    neckDisplacement = fzero(radiusError,sort([0,2*rn*neckChange]));

    tangentR = h*cos(original(1,:));
    tangentZ = -h*sin(original(1,:));
    displacedR = r+neckDisplacement*neckBump;
    displacedTangent = tangentR+neckDisplacement*neckSlope;
    areaError = @(b) trapz(u,(displacedR+b*areaBump).* ...
        hypot(displacedTangent+b*areaSlope,tangentZ))-area(end);
    areaDisplacement = fzero(areaError,[-0.05,0.05]);
    slopeChange = neckDisplacement*neckSlope+areaDisplacement*areaSlope;
    r = displacedR+areaDisplacement*areaBump;
    speed = hypot(tangentR+slopeChange,tangentZ);
    arc = cumtrapz(u,speed);
    h = arc(end);
    v = arc/h;
    weight = r.*speed;
    newArea = cumtrapz(u,weight);

    % Match cumulative membrane area to carry protein with the deformed surface.
    sourceU = interp1(area/area(end),u,newArea/newArea(end),'pchip');
    phi = deval(sol,sourceU,12);
    values = original;
    values(1,:) = original(1,:)+atan2(tangentZ.*slopeChange, ...
        sol.parameters^2+tangentR.*slopeChange);
    values(3,:) = r;
    values(15,:) = newArea;
    shapePP = spline(v,values);
    mappedMesh = interp1(u,v,sol.x,'pchip');
    neckCenter = interp1(u,v,neck,'pchip');
    backgroundMesh = linspace(v(1),v(end),81);
    neckMesh = linspace(max(v(1),neckCenter-0.10),min(v(end),neckCenter+0.10),91);
    mesh = unique([mappedMesh,backgroundMesh,neckMesh]);
    sol.x = mesh;
    sol.parameters = h;
    [sol.y,sol.idata.yp,sol.idata.ymid] = sampleInitialRows(shapePP,mesh);
    anglePP = spline(v,values(1,:));
    [sol.y(2,:),sol.idata.yp(2,:),sol.idata.ymid(2,:)] = ...
        sampleInitialRows(derivativePP(anglePP),mesh);
end

[newRn,neckIndex] = neckMinimum(v,r);
amplitude = 0;
energyAdded = 0;
% Protein perturbation: redistribute protein locally at fixed shape.
if proteinEnergy > 0
    nearNeck = abs(v-v(neckIndex)) < 0.12;
    phiPP = spline(v,[0,phi,0]);
    phiSlope = ppval(derivativePP(phiPP),v);
    candidates = find(nearNeck);
    [~,peak] = max(abs(phiSlope(nearNeck)));
    interface = v(candidates(peak));
    width = 0.06;
    window = compactBump(v,interface,width).*phi.*(1-phi);
    offset = (v-interface)/width;
    offset = offset-trapz(u,weight.*window.*offset)/trapz(u,weight.*window);
    profile = window.*offset;
    profile = profile/max(abs(profile));
    profilePP = spline(v,[0,profile,0]);
    profileSlope = ppval(derivativePP(profilePP),v);
    baseEnergy = compositionEnergy(u,weight,phi,phiSlope/h,Lambda,chi,kT,ap);
    energyError = @(a) compositionEnergy(u,weight,phi+a*profile, ...
        (phiSlope+a*profileSlope)/h,Lambda,chi,kT,ap)-baseEnergy-proteinEnergy;
    upper = 0.95*min([(1-phi(profile>0))./profile(profile>0), ...
                     -phi(profile<0)./profile(profile<0)]);
    if energyError(upper) < 0
        error('Cannot reach the requested protein energy with this local perturbation.');
    end
    amplitude = fzero(energyError,[0,upper]);
    phi = phi+amplitude*profile;
    energyAdded = energyError(amplitude)+proteinEnergy;
end

if neckChange ~= 0 || proteinEnergy > 0
    phiPP = spline(v,[0,phi,0]);
    [sol.y(12,:),sol.idata.yp(12,:),sol.idata.ymid(12,:)] = ...
        sampleInitialRows(phiPP,sol.x);
    [sol.y(13,:),sol.idata.yp(13,:),sol.idata.ymid(13,:)] = ...
        sampleInitialRows(derivativePP(phiPP),sol.x);
end

% Verify the state interpolated by the time-stepping loop.
checkU = linspace(sol.x(1),sol.x(end),10001);
check = deval(sol,checkU,[3,12]);
[checkRn,~] = neckMinimum(checkU,check(1,:));
relativeAreaError = trapz(checkU,h*check(1,:))/area(end)-1;
relativeProteinError = trapz(checkU,h*check(1,:).*check(2,:))/proteinTotal-1;
if any(check(2,:) <= 0 | check(2,:) >= 1) || ...
        abs(relativeAreaError) > 1e-6 || abs(relativeProteinError) > 1e-6 || ...
        abs(checkRn-newRn) > 1e-5 || abs(checkRn/rn-1-neckChange) > 5e-4
    error('The perturbed initial state failed its area, protein, or neck-radius check.');
end
fprintf(['Initial r_n: %.8f -> %.8f (%+.2f%%); protein perturbation: %.4f kBT\n', ...
    'Relative area change: %+.3g; relative protein change: %+.3g\n'], ...
    rn,checkRn,100*(checkRn/rn-1),energyAdded,relativeAreaError,relativeProteinError);
end

function [rn,index] = neckMinimum(u,r)
indices = find(islocalmin(r) & u > u(1) & u < u(end));
if numel(indices) ~= 1 || r(indices) <= 0
    error('The initial profile must have one positive interior neck minimum.');
end
index = indices(1);
rn = r(index);
end

function [bump,slope] = compactBump(u,center,width)
q = (u-center)/width;
inside = max(1-q.^2,0);
bump = inside.^4;
slope = -8*q.*inside.^3/width;
end

function [values,slopes,midpoints] = sampleInitialRows(pp,mesh)
values = ppval(pp,mesh);
slopes = ppval(derivativePP(pp),mesh);
midpoints = ppval(pp,(mesh(1:end-1)+mesh(2:end))/2);
end

function dp = derivativePP(pp)
[breaks,coefficients,~,order,dimension] = unmkpp(pp);
dp = mkpp(breaks,coefficients(:,1:end-1).*(order-1:-1:1),dimension);
end

function energy = compositionEnergy(u,weight,phi,phiS,Lambda,chi,kT,ap)
% Same fourth-order endpoint regularization used for the paper's energy.
cutoff = 1e-4;
q = min(phi,1-phi);
z = max(q,cutoff);
d = q-z;
entropy = z.*log(z)+(1-z).*log1p(-z)+log(z./(1-z)).*d ...
    +(1./z+1./(1-z)).*d.^2/2+(-1./z.^2+1./(1-z).^2).*d.^3/6 ...
    +(2./z.^3+2./(1-z).^3).*d.^4/24;
density = Lambda*phiS.^2+chi*phi.*(1-phi)+kT*entropy;
energy = 2*pi/(ap*kT)*trapz(u,weight.*density);
end
