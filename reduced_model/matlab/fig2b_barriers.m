function [fig, fullBarrier, reducedBarrier] = fig2b_barriers(lambdaTilde)
% FIG2B_BARRIERS  Fig. 2(b) barriers at one Lambda~ (0.15, 0.30 or 0.45; default 0.30):
% reduced model (from the landscapes of REDUCED_LANDSCAPE) and full model (energy
% difference between the stable open-neck branch and the neighboring unstable branch
% at equal chi). Reads data_all/reduced_model/ in the repository root (see README).

if nargin < 1
    lambdaTilde = 0.30;
end
scriptDir = fileparts(mfilename('fullpath'));
dataDir = fullfile(scriptDir,'..','..','data_all','reduced_model');
addpath(fullfile(scriptDir,'..','..','full_model'));   % entropy_cutoff
numNumericalChiSamples = 8;   % chi samples on the common range of the two branches

% README code units. Fig. 2(b) plots barriers in kBT, unlike Fig. 1 energy
% curves, which are divided by 8*pi*kappa.
R0 = 1;
kBT = 0.02;
kappa = 0.4;
ap = 0.01;
c0 = 0;
pressure = 0;
force = 0;
energyScale = kBT;
Lambda = lambdaTilde*R0^2*kBT;     % code units
lambdaToken = strrep(sprintf('%.2f',lambdaTilde),'.','p');

reducedBarrier = collectSimplifiedBarrier( ...
    dataDir,lambdaToken,kBT,energyScale);
fullBarrier = collectNumericalBarrier( ...
    dataDir,lambdaToken,Lambda,kBT,energyScale, ...
    kappa,c0,ap,pressure,force,numNumericalChiSamples);
if isempty(reducedBarrier.chiTilde) || isempty(fullBarrier.chiTilde)
    error(['Both reduced landscapes and the inner/outer full-model branches ', ...
        'are required in %s for Lambda=%s.'],dataDir,lambdaToken);
end

fig = figure(2);
clf(fig);
set(fig,'Color','w','Units','pixels','Position',[80,80,900,650], ...
    'Name','Fig. 2(b): energy barriers','NumberTitle','off');
ax = axes(fig);
hold(ax,'on');
curveColor = [0.5566667,0.4526667,0.1503333];
plot(ax,reducedBarrier.chiTilde,reducedBarrier.deltaFTilde,'--^', ...
    'Color',curveColor,'MarkerFaceColor','w','MarkerSize',7.8, ...
    'LineWidth',2.2,'DisplayName','$\Delta\widetilde F_{\rm re}$');
plot(ax,fullBarrier.chiTilde,fullBarrier.deltaFTilde,'-o', ...
    'Color',curveColor,'MarkerFaceColor','w','MarkerSize',7, ...
    'LineWidth',2.2,'DisplayName','$\Delta\widetilde F_{\rm full}$');
plot(ax,xlim(ax),[0,0],'-','Color',[0.3,0.3,0.3], ...
    'LineWidth',1,'HandleVisibility','off');
xlabel(ax,'$\widetilde\chi$','Interpreter','latex','FontSize',20);
ylabel(ax,'$\Delta\widetilde F$', ...
    'Interpreter','latex','FontSize',20);
title(ax,['$\widetilde\Lambda=',num2str(lambdaTilde,'%.3g'),'$'], ...
    'Interpreter','latex','FontSize',20,'FontWeight','normal');
legend(ax,'Interpreter','latex','FontSize',15,'Location','best','Box','off');
set(ax,'Box','on','TickDir','in','LineWidth',1, ...
    'FontName','Arial','FontSize',16,'TickLabelInterpreter','latex');
ylim(ax,[0,max([reducedBarrier.deltaFTilde;fullBarrier.deltaFTilde])*1.12]);
drawnow;
fprintf('Lambda~ = %.4g: %d reduced and %d full-model barrier points.\n', ...
    lambdaTilde,numel(reducedBarrier.chiTilde),numel(fullBarrier.chiTilde));
end

function result = collectSimplifiedBarrier(dataDir,lambdaToken,kBT,energyScale)
% Same discrete-extrema test and energy difference as the original script.
result = emptyBarrierResult();
fileList = dir(fullfile(dataDir, ...
    sprintf('landscape_Lambda%s_chi*.mat',lambdaToken)));
chiPhysical = [];
chiTilde = [];
deltaFTilde = [];

for fileNumber = 1:numel(fileList)
    fileName = fileList(fileNumber).name;
    chiVal = parseChiFromName(fileName,kBT);
    if ~isfinite(chiVal)
        continue;
    end
    dataFile = fullfile(dataDir,fileName);
    sourceData = load(dataFile,'aVec','Fmin','p');
    validateFields(sourceData,{'aVec','Fmin','p'},dataFile);
    aVec = sourceData.aVec(:);
    Fmin = real(sourceData.Fmin(:));
    valid = isfinite(aVec) & isfinite(Fmin);
    if ~any(valid)
        continue;
    end
    aScan = aVec(valid);
    FScan = Fmin(valid);
    [~,order] = sort(aScan);
    FScan = FScan(order);
    FAtZero = 8*pi*(sourceData.p.kappaAlpha+ ...
        sourceData.p.kappaBeta)/energyScale;
    FCurve = [FAtZero;FScan/energyScale];
    interior = 2:(numel(FCurve)-1);
    localMin = interior(FCurve(interior) < FCurve(interior-1) & ...
        FCurve(interior) < FCurve(interior+1));
    localMax = interior(FCurve(interior) > FCurve(interior-1) & ...
        FCurve(interior) > FCurve(interior+1));
    if numel(localMin) ~= 1 || numel(localMax) ~= 1
        continue;
    end
    deltaF = FCurve(localMax)-FCurve(localMin);
    if deltaF < 0
        continue;
    end
    chiPhysical(end+1,1) = chiVal; %#ok<AGROW>
    chiTilde(end+1,1) = chiVal/kBT; %#ok<AGROW>
    deltaFTilde(end+1,1) = deltaF; %#ok<AGROW>
end

if isempty(deltaFTilde)
    return;
end
[chiTilde,order] = sort(chiTilde);
result.chiPhysical = chiPhysical(order);
result.chiTilde = chiTilde;
result.deltaFTilde = deltaFTilde(order);
end

function result = collectNumericalBarrier(dataDir,lambdaToken,Lambda, ...
        kBT,energyScale,kappa,c0,ap,pressure,force,numChiSamples)
% Same branch-energy calculation, averaging of repeated chi, common-range
% restriction, pchip interpolation and absolute energy difference as before.
result = emptyBarrierResult();
[innerBranch,outerBranch] = loadNumericalBranches(dataDir,lambdaToken);
if isempty(innerBranch.cchi) || isempty(outerBranch.cchi)
    return;
end
innerEnergy = computeBranchEnergy(innerBranch,Lambda, ...
    kappa,kBT,c0,ap,pressure,force);
[innerChi,innerEnergy] = singleValuedCurve(innerBranch.cchi,innerEnergy);
outerEnergy = computeBranchEnergy(outerBranch,Lambda, ...
    kappa,kBT,c0,ap,pressure,force);
[outerChi,outerEnergy] = singleValuedCurve(outerBranch.cchi,outerEnergy);
commonRange = [max(min(innerChi),min(outerChi)), ...
    min(max(innerChi),max(outerChi))];
if commonRange(1) >= commonRange(2)
    return;
end
chiPhysical = linspace(commonRange(1),commonRange(2),numChiSamples);
innerEnergy = interp1(innerChi,innerEnergy,chiPhysical,'pchip');
outerEnergy = interp1(outerChi,outerEnergy,chiPhysical,'pchip');
result.chiPhysical = chiPhysical(:);
result.chiTilde = chiPhysical(:)/kBT;
result.deltaFTilde = abs(outerEnergy(:)-innerEnergy(:))/energyScale;
end

function [innerBranch,outerBranch] = loadNumericalBranches(dataDir,lambdaToken)
innerBranch = struct('cchi',[],'ssol',{{}});
outerBranch = struct('cchi',[],'ssol',{{}});
innerFile = fullfile(dataDir, ...
    sprintf('branch_Lambda%s_inner.mat',lambdaToken));
outerFile = fullfile(dataDir, ...
    sprintf('branch_Lambda%s_outer.mat',lambdaToken));
if isfile(innerFile) && isfile(outerFile)
    innerBranch = load(innerFile,'cchi','ssol');
    outerBranch = load(outerFile,'cchi','ssol');
end
end

function energy = computeBranchEnergy(data,Lambda,ka,kT,c0,ap,pressure,force)
% Kept identical to the original figure2_newmesh energy evaluation.
cchi = data.cchi(:).';
energy = zeros(size(cchi));
for solutionNumber = 1:numel(cchi)
    sol = data.ssol{solutionNumber};
    u = linspace(min(sol.x),max(sol.x),10000);
    dphi = interp1(sol.x,sol.y(8,:),u,'spline');
    h = sol.parameters(1);
    r = interp1(sol.x,sol.y(3,:),u,'spline');
    psi = interp1(sol.x,sol.y(1,:),u,'spline');
    phi = interp1(sol.x,sol.y(7,:),u,'spline');
    mixingEntropy = regularizedMixingEntropy(real(phi));
    dpsi = interp1(sol.x,sol.y(2,:),u,'spline');
    chi = cchi(solutionNumber);
    bending = 0.5*h*ka.*r.* ...
        (dpsi/h-c0*phi+sin(psi)./r).^2;
    protein = r/(ap*h).* ...
        (Lambda*dphi.^2-chi*h^2.*(phi-1).*phi+ ...
        h^2*kT.*mixingEntropy);
    forceEnergy = -force*h.*sin(psi);
    pressureEnergy = 0.5*h*pressure.*r.^2.*sin(psi);
    energy(solutionNumber) = real(2*pi*trapz( ...
        u,bending+protein+forceEnergy+pressureEnergy));
end
end

function entropy = regularizedMixingEntropy(phi)
% Match the fourth-order Taylor continuation used by the full-model solver.
delta = entropy_cutoff();
entropy = zeros(size(phi));
middle = phi >= delta & phi <= 1-delta;
lower = phi < delta;
upper = phi > 1-delta;

phiMiddle = phi(middle);
entropy(middle) = phiMiddle.*log(phiMiddle) + ...
    (1-phiMiddle).*log(1-phiMiddle);
entropy(lower) = entropyTaylor(phi(lower),delta);
entropy(upper) = entropyTaylor(phi(upper),1-delta);
end

function value = entropyTaylor(phi,center)
% Fourth-order Taylor polynomial of phi*log(phi)+(1-phi)*log(1-phi).
offset = phi-center;
f0 = center*log(center) + (1-center)*log(1-center);
f1 = log(center) - log(1-center);
f2 = 1/center + 1/(1-center);
f3 = -1/center^2 + 1/(1-center)^2;
f4 = 2/center^3 + 2/(1-center)^3;
value = f0 + f1.*offset + (f2/2).*offset.^2 + ...
    (f3/6).*offset.^3 + (f4/24).*offset.^4;
end

function chiVal = parseChiFromName(fileName,kBT)
token = regexp(fileName,'_chi([0-9]+p[0-9]+)\.mat$','tokens','once');
if isempty(token)
    chiVal = NaN;
else
    chiVal = str2double(strrep(token{1},'p','.'))*kBT;
end
end

function validateFields(data,fieldNames,sourceFile)
for fieldNumber = 1:numel(fieldNames)
    if ~isfield(data,fieldNames{fieldNumber})
        error('File %s is missing field: %s',sourceFile,fieldNames{fieldNumber});
    end
end
end

function result = emptyBarrierResult()
result = struct('chiPhysical',[],'chiTilde',[],'deltaFTilde',[]);
end

function [xUnique,yUnique] = singleValuedCurve(x,y)
x = x(:);
y = y(:);
valid = isfinite(x) & isfinite(y);
x = x(valid);
y = y(valid);
[x,order] = sort(x);
y = y(order);
[xUnique,~,groupIndex] = unique(x,'stable');
yUnique = accumarray(groupIndex,y,[],@mean);
end
