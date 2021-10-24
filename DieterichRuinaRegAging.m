function [yp] = DieterichRuinaRegAging(~,y,ss)
% This file describes the evolution of the ordinary
% differential equation y' = f(t,y), where the state 
% vector y is 
%
%        /        s          \
%    y = |       tau         |
%        | log(theta Vo / D_rs) |
%        \    log(V / Vo)    /

% based on the regularized form of Dieterich-Ruina rate-and-state friction 
% and using the aging law
%
% Velocity is determined by the balance of shear resistance and shear stress
% Resistance : tau = a sigma asinh( V/2Vo exp( (fo + b psi ) / a))
% Stress     : tau = tauo + f(z,t) - eta*(V-Vpl)  
% Where we use radiation damping to approximate the inertial terms
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
% parameters, or the loading plate rate. Loading is done purely through
% back slip at plate rate Vpl.
%
% Instead of directly integrating numerically the aging law
%
%    d theta / dt = 1 - V theta / D_rs
%
% as is, we operate the following change of variable
%
%    phi = ln (theta Vo / D_rs)
%
% where we obtain the evolution law for the new variable
% 
%    d phi / dt = ( Vo exp(-phi) - V ) / D_rs
%
% Given the regularized form of Dieterich-Ruina R+S we can express,
%    1 dV      K (V - Vpl) - b sigma dphi / dt Q
%    - --  =  ----------------------------------
%    V dt           a sigma Q  + eta V
%  
% where
%
%                            1
%    Q = -----------------------------------------------
%        /                                             \(1/2)
%        |  1 + [2 Vo / V exp(-(fo + b phi) / a )]^2   |
%        \                                             /
% 
% Note that d/dt log(V/Vo) = 1/V dV/dt, and is much more efficient to
% integrate than dV/dt alone

% State variable
th=y(3:ss.dgf:end);

% Slip rate
V = ss.Vo.* exp(y(4:ss.dgf:end));

% Initialize Time Derivative
yp=zeros(size(y));

% Slip
yp(1:ss.dgf:end)=V;

% State Variable
dth = (ss.Vo.*exp(-th)-V)./ss.Drs;
yp(3:ss.dgf:end)=dth;

% Slip Velocity
func = ss.K*(V-ss.Vpl);
f1=2*ss.Vo./V.*exp(-(ss.fo+ss.b.*th)./ss.a);
f2=1./sqrt(1+f1.^2);

yp(4:ss.dgf:end) = (func - ss.b.*ss.sigma.*dth.*f2)./ ...
                    (ss.a.*ss.sigma.*f2 + ss.eta.*V);


% Evolution of shear stress 
yp(2:ss.dgf:end)=func - ss.eta.*V.*yp(4:ss.dgf:end);



end
