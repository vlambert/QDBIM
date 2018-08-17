function [yp] = odeRegDRaging(~,y,ss)
% odeRegDRaging describes the evolution of the ordinary
% differential equation y' = f(t,y), where the state 
% vector y is 
%
%        /        s          \
%    y = |       tau         |
%        | log(theta Vo / L) |
%        \         V         /
%
% based on the regularized form of Dieterich-Ruina rate-and-state friction 
% and using the aging law
%
% Velocity is determined by the balance of shear resistance and shear stress
% Resistance : tau = a sigma asinh( V/2Vo exp( (fo + b psi ) / a))
% Stress     : tau = tauo + f(z,t) - eta*(V-Vpl)  
% Where we use radiation damping for approximating the inertial terms
% as the quasi-dynamic approximation with f(z,t) = K(delta - Vpl*t)
%
% This is done by taking the time-derivative of both equations and equating
% the two.
%
%   a sigma  alpha/ sqrt(1 + alpha^2) ( 1/V dV/dt + b/a dPhi/dt) 
%     
%                   = K( V - Vpl) - eta dV/dt

% where 
%      alpha = V / 2Vo exp( ( fo + b phi) / a)
%
%
%
% Note this form assumes no time variation in sigma, the frictional
% parameters, or the loading plate rate

% Instead of directly integrating numerically the aging law
%
%    d theta / dt = 1 - V theta / L
%
% as is, we operate the following change of variable
%
%    phi = ln (theta Vo / L)
%
% where we obtain the evolution law for the new variable
% 
%    d phi / dt = ( Vo exp(-phi) - V ) / L
%
% Given the regularized form of Dieterich-Ruina R+S we can express,
%     dV      K (V - Vpl) - b sigma dphi / dt Q
%     --  =  ----------------------------------
%     dt           a sigma Q / V + eta
%  
% where
%
%                            1
%    Q = -----------------------------------------------
%        /                                             \(1/2)
%        |  1 + [2 Vo / V exp(-(fo + b phi) / a )]^2   |
%        \                                             /


% State variable
th=y(3:ss.dgf:end);

% Slip rate
V = y(4:ss.dgf:end);

% Initialize Time Derivative
yp=zeros(size(y));

% Slip
yp(1:ss.dgf:end)=V;

% State Variable
dth = (ss.Vo.*exp(-th)-V)./ss.L;
yp(3:ss.dgf:end)=dth;

% Slip Velocity
func = ss.K*(V-ss.Vpl);
f1=2*ss.Vo./V.*exp(-(ss.fo+ss.b.*th)./ss.a);
f2=1./sqrt(1+f1.^(2));

yp(4:ss.dgf:end) = (func - ss.b.*ss.sigma.*dth.*f2)./ ...
                    (ss.a.*ss.sigma.*f2./V + ss.eta);


% Evolution of shear stress 
yp(2:ss.dgf:end)=func - ss.eta.*yp(4:ss.dgf:end);



end
