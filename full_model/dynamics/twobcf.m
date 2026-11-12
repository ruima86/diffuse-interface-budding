function res = twobcf(ya,yb,para,ka,C0,kT,ri,Lambda,chi,ap)
h0 = para;

% psi = y(1,:);
% dpsi = y(2,:);
% r = y(3,:);
% z = y(4,:);
% al = y(5,:);
% be = y(6,:);
% vu = y(7,:);
% la = y(8,:);
% vphi = y(9,:);
% dvphi = y(10,:);
% ju = y(11,:);
% rho = y(12,:);
% drho = y(13,:);
% ddrho = y(14,:);
% area = y(15,:);
c0 = C0;
epsilon=0.0001;

alpha_i =(1/6).*((-6).*ap.^(-1).*chi.*ri.*((-1)+ya(12)).*ya(12)+(-6).*c0.* ...
  h0.^(-1).*ka.*ya(12).*(ri.*ya(2)+(-1).*c0.*h0.*ri.*ya(12))+(-3).* ...
  h0.^(-2).*ka.*ri.^(-1).*(ri.*ya(2)+(-1).*c0.*h0.*ri.*ya(12)).^2+ ...
  6.*ap.^(-1).*kT.*ri.*((1+(-1).*epsilon)>=(ya(12))&(ya(12))>=(epsilon)).* ...
  ((-1).*log(1+(-1).*ya(12)).*((-1)+ya(12))+log(ya(12)).*ya(12))+ ...
  ap.^(-1).*kT.*ri.*((ya(12))<(epsilon)).*((-6).*((-1)+epsilon) ...
  .*log(1+(-1).*epsilon)+6.*epsilon.*log(epsilon)+6.*(log(1+(-1).* ...
  epsilon)+(-1).*log(epsilon)).*(epsilon+(-1).*ya(12))+(-3).*((-1)+ ...
  epsilon).^(-1).*epsilon.^(-1).*(epsilon+(-1).*ya(12)).^2+(1/2).*(( ...
  -1).*((-1)+epsilon).^(-3)+epsilon.^(-3)).*(epsilon+(-1).*ya(12)) ...
  .^4+((-1)+epsilon).^(-2).*epsilon.^(-2).*((-1)+2.*epsilon).*((-1) ...
  .*epsilon+ya(12)).^3)+ap.^(-1).*kT.*ri.*((epsilon+ya(12))>(1) ...
  ).*((-6).*((-1)+epsilon).*log(1+(-1).*epsilon)+6.*epsilon.*log( ...
  epsilon)+6.*(log(1+(-1).*epsilon)+(-1).*log(epsilon)).*((-1)+ ...
  epsilon+ya(12))+(-3).*((-1)+epsilon).^(-1).*epsilon.^(-1).*((-1)+ ...
  epsilon+ya(12)).^2+(-1).*((-1)+epsilon).^(-2).*epsilon.^(-2).*(( ...
  -1)+2.*epsilon).*((-1)+epsilon+ya(12)).^3+(1/2).*((-1).*((-1)+ ...
  epsilon).^(-3)+epsilon.^(-3)).*((-1)+epsilon+ya(12)).^4)+(-6).* ...
  ap.^(-1).*h0.^(-2).*Lambda.*ri.*ya(13).^2);




res = [ ya(1)-0;...       % psi(0) = 0;
    ya(3)-ri;...          % r(0) = ri;
    ya(5)-alpha_i;...     % alpha(0) = alpha_i
    ya(6)+0;...        % beta(0) = 0;
    ya(7)-0;...           % vu(0) = 0;
    ya(9)-0;...           % vphi(0) = 0;
    ya(11)-0;...          % ju(0) = 0;
    ya(13)-0;...          % rho'(0) = 0;
    ya(15)-0;...          % area(0) = 0;
    yb(1)-pi;...          % psi(1) = pi;
    yb(3)-ri;...          % r(1) = ri;
    yb(6)+0;...        % beta(1) = 0;
    yb(7)-0;...           % vu(1) = 0;
    yb(9)-0;...           % vphi(1) = 0;
    yb(11)-0;...          % ju(1) = 0; 
    yb(13)-0];            % rho'(1) = 0