% FIG1CDE_NECK_WIDTH_TENSION  Neck radius, interface width and line tension along one
% full-model branch (Fig. 1(c-e)).
% Change only this input file, then run the script from MATLAB.
% The included Lambda~ = 0.30 branch is a test case, not a Fig. 1(c-e)
% branch; those panels use Lambda~ = 0.10, 0.40, and 0.80 respectively.
clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
S = load(fullfile(scriptDir,'data','branch_Lambda0p30.mat'));

% If the chosen branch's fold indices have been checked, enter them here.
% Do not reuse indices from another .mat file: its solution order may differ.
turningIndices = [];

% README code units. All quantities in metrics and on the axes below are
% dimensionless: chi~=chi/kBT, r_n~=r_n/R0, w~=w/R0,
% gamma~=gamma*R0/kBT, Lambda~=Lambda/(R0^2*kBT).
R0 = 1;
kBT = 0.02;
ap = 0.01;
if ~isfield(S,'cchi') || ~isfield(S,'ssol') || ...
        numel(S.cchi) ~= numel(S.ssol)
    error('The branch .mat file must contain matching cchi and ssol arrays.');
end
lambdaCode = storedLambda(S);
lambdaTilde = lambdaCode/(R0^2*kBT);
metrics = calculateMetrics(S,lambdaCode,ap,R0,kBT);
metrics.lambda_tilde = lambdaTilde;

fig = figure(1);
clf(fig);
set(fig,'Color','w','Units','pixels','Position',[80,80,900,650], ...
    'Name','Neck radius, interface width and line tension', ...
    'NumberTitle','off','Renderer','painters');
ax = axes(fig);
hold(ax,'on');

neckColor = [55,55,55]/255;
widthColor = [184,134,11]/255;
gammaColor = [190,50,45]/255;
markerArea = 8;

neckPlot = metrics.neck;
neckPlot(neckPlot <= 0 | ~metrics.validPhase) = NaN;
widthPlot = metrics.width;
widthPlot(~metrics.validPhase) = NaN;
gammaPlot = metrics.gamma;
gammaPlot(~metrics.validPhase) = NaN;

yyaxis(ax,'left');
hNeck = scatter(ax,metrics.chi,neckPlot,markerArea,neckColor,'filled');
hWidth = scatter(ax,metrics.chi,widthPlot,markerArea,widthColor,'filled');
plotTurningPoints(ax,metrics.chi,neckPlot,turningIndices);
plotTurningPoints(ax,metrics.chi,widthPlot,turningIndices);
ylabel(ax,'$\widetilde r_{\rm n},\ \widetilde w$', ...
    'Interpreter','latex','FontSize',18);

yyaxis(ax,'right');
hGamma = scatter(ax,metrics.chi,gammaPlot,markerArea,gammaColor,'filled');
plotTurningPoints(ax,metrics.chi,gammaPlot,turningIndices);
ylabel(ax,'$\widetilde\gamma$', ...
    'Interpreter','latex','FontSize',18);

xlabel(ax,'$\widetilde\chi$','Interpreter','latex','FontSize',18);
title(ax,['$\widetilde\Lambda=',num2str(lambdaTilde,'%.3g'),'$'], ...
    'Interpreter','latex','FontSize',20,'FontWeight','normal');
legend(ax,[hNeck,hWidth,hGamma], ...
    {'$\widetilde r_{\rm n}$ neck radius', ...
     '$\widetilde w$ interface width', ...
     '$\widetilde\gamma$ line tension'}, ...
    'Interpreter','latex','FontSize',14,'Location','best','Box','off');
set(ax,'Box','on','TickDir','in','LineWidth',1, ...
    'FontName','Arial','FontSize',15,'TickLabelInterpreter','latex');
ax.YAxis(1).Color = [0,0,0];      % left axis carries both r_n and w
ax.YAxis(2).Color = gammaColor;
drawnow;
fprintf('Lambda~ = %.4g; %d solutions; chi~ from %.4g to %.4g.\n', ...
    lambdaTilde,numel(metrics.chi),min(metrics.chi),max(metrics.chi));

function lambdaCode = storedLambda(data)
% Older branch files use La; new continuation files use Lambda.
if isfield(data,'Lambda')
    lambdaCode = data.Lambda;
elseif isfield(data,'La')
    lambdaCode = data.La;
else
    error('The branch file must contain Lambda or La in MATLAB code units.');
end
if ~isscalar(lambdaCode) || ~isfinite(lambdaCode) || lambdaCode < 0
    error('Stored Lambda/La must be one nonnegative finite number.');
end
end

function metrics = calculateMetrics(data,Lambda,ap,R0,kBT)
% The next three definitions are unchanged from the original script:
% r_n is an interior local minimum of r(u); w=(phi_max-phi_min)/max_s|dphi/ds|;
% gamma=(2*Lambda/ap)*integral (dphi/ds)^2 ds.
cchi = data.cchi;
ssol = data.ssol;
count = numel(cchi);
neck = zeros(1,count);
lineTension = zeros(1,count);
interfaceWidth = zeros(1,count);

for i = 1:count
    sol = ssol{i};
    u = linspace(min(sol.x),max(sol.x),1000);
    r = interp1(sol.x,sol.y(3,:),u,'spline');
    dphi = interp1(sol.x,sol.y(8,:),u,'spline');
    phi = interp1(sol.x,sol.y(7,:),u,'spline');
    h = sol.parameters(1);

    minValues = r(islocalmin(r));      % interior local minima of r(u)
    if isscalar(minValues)
        neck(i) = minValues;
    end

    lineDensity = 2*Lambda/ap .* (dphi/h).^2 .* h;
    lineTension(i) = trapz(u,lineDensity);
    interfaceWidth(i) = (max(phi)-min(phi))*h/max(abs(dphi));
end

metrics.chi = cchi./kBT;
metrics.neck = neck./R0;
metrics.width = interfaceWidth./R0;
metrics.gamma = lineTension.*R0./kBT;
metrics.validPhase = metrics.gamma > 1e-3;
end

function plotTurningPoints(ax,x,y,indices)
% Open circles identify verified folds of this particular stored branch.
indices = unique(indices(:).','stable');
indices = indices(indices >= 1 & indices <= numel(x));
indices = indices(isfinite(x(indices)) & isfinite(y(indices)));
if isempty(indices)
    return
end
plot(ax,x(indices),y(indices),'o','LineStyle','none', ...
    'MarkerSize',9,'MarkerFaceColor','none', ...
    'MarkerEdgeColor',[0,0,0],'LineWidth',1.5, ...
    'HandleVisibility','off');
end

