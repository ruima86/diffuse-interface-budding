function res = twobcf(ya,yb,para,ka,C0,kT,ri,Lambda,chi,ap,area0,rho0)
% TWOBCF  Boundary conditions for SHAPEF (SM Sec. S1, Table of boundary conditions).
% The poles are placed at r = ri (ri = 1e-5 in the run scripts) instead of r = 0
% to avoid the coordinate singularity; alpha(0) follows from the h-condition at
% the pole. area0 = A/(2*pi) and rho0 = phi0 fix the total area and protein number.
h0 = para(1);
la0 = para(2);
% h0 = para(1);
% la = para(2);
% mu = para(3);
% 
% psi = y(1,:);
% dpsi = y(2,:);
% r = y(3,:);
% %z = y(4,:);
% al = y(5,:);
% be = y(6,:);
% rho = y(7,:);
% drho = y(8,:);
% %area = y(9,:);
% %Num = y(10,:);
% la = y(11,:);

c0 = C0;

epsilon = entropy_cutoff();       % 1e-4 in this work

alpha_i = (1/6).*((-6).*ap.^(-1).*chi.*ri.*((-1)+ya(7)).*ya(7)+(-6).*c0.* ...
  h0.^(-1).*ka.*ya(7).*(ri.*ya(2)+(-1).*c0.*h0.*ri.*ya(7))+(-3).* ...
  h0.^(-2).*ka.*ri.^(-1).*(ri.*ya(2)+(-1).*c0.*h0.*ri.*ya(7)).^2+6.* ...
  ap.^(-1).*kT.*ri.*((1+(-1).*epsilon)>=(ya(7))&(ya(7))>=(epsilon)).*(( ...
  -1).*log(1+(-1).*ya(7)).*((-1)+ya(7))+log(ya(7)).*ya(7))+ap.^(-1) ...
  .*kT.*ri.*((ya(7))<(epsilon)).*((-6).*((-1)+epsilon).*log(1+( ...
  -1).*epsilon)+6.*epsilon.*log(epsilon)+6.*(log(1+(-1).*epsilon)+( ...
  -1).*log(epsilon)).*(epsilon+(-1).*ya(7))+(-3).*((-1)+epsilon).^( ...
  -1).*epsilon.^(-1).*(epsilon+(-1).*ya(7)).^2+(1/2).*((-1).*((-1)+ ...
  epsilon).^(-3)+epsilon.^(-3)).*(epsilon+(-1).*ya(7)).^4+((-1)+ ...
  epsilon).^(-2).*epsilon.^(-2).*((-1)+2.*epsilon).*((-1).*epsilon+ ...
  ya(7)).^3)+ap.^(-1).*kT.*ri.*((epsilon+ya(7))>(1)).*((-6).*(( ...
  -1)+epsilon).*log(1+(-1).*epsilon)+6.*epsilon.*log(epsilon)+6.*( ...
  log(1+(-1).*epsilon)+(-1).*log(epsilon)).*((-1)+epsilon+ya(7))+( ...
  -3).*((-1)+epsilon).^(-1).*epsilon.^(-1).*((-1)+epsilon+ya(7)).^2+ ...
  (-1).*((-1)+epsilon).^(-2).*epsilon.^(-2).*((-1)+2.*epsilon).*(( ...
  -1)+epsilon+ya(7)).^3+(1/2).*((-1).*((-1)+epsilon).^(-3)+ ...
  epsilon.^(-3)).*((-1)+epsilon+ya(7)).^4)+(-6).*ap.^(-1).*h0.^(-2) ...
  .*Lambda.*ri.*ya(8).^2);




res = [ ya(1)-0;...       % psi(0) = 0;
    ya(3)-ri;...          % r(0) = ri;
    ya(5)-alpha_i;...     % alpha(0) = alpha_i
    ya(8)-0;...           % rho'(0) = 0;
    ya(9)-0;...           % area(0) = 0;
    ya(10)-0;...          % Num(0) = 0;
    ya(11)-la0;...        % lambda(0) = la0;
    yb(1)-pi;...          % psi(1) = pi;
    yb(3)-ri;...          % r(1) = ri;
    yb(8)-0;...            % rho'(1) = 0
    yb(9)-area0;...        % area(1) = area0;
    yb(10)-rho0*area0;...  % Num(1) = area0*phi0;
    ya(4)+yb(4);...       % z(0) + z(1) = 0;
    ya(6)+yb(6);...       % beta(0) + beta(1) = 0;
    ];         
