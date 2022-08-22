%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%   2D Quasi-dynamic simulation of the evolution of    %
%   slip and stress on a fault in antiplane strain     %
%   governed by rate- and state-dependent friction     %
%                                                      %
%                Valere Lambert, 2022                  %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%
clear all;close all;

addpath('include/');
%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                  Useful Functions                    %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%

% boxcar function
BC=@(x) (x+0.5>=0)-(x-0.5>=0);

% Heaviside function
HS=@(x) 0+x>=0;

% ramp function
Ramp=@(x) x.*BC(x-1/2)+HS(x-1);

%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%            Stress Interaction Functions              %
%                                                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%
% Full space 2D Antiplane problem (x2 x x3): 
% Out-of-plane displacements u1 are non-zero, u2 = u3 = 0
% u1 are assumed uniform along strike, only vary in x2 x x3 plane

% density (kg/m^3)
rho = 2670;

% shear wave speed (m/s)
Vs = 3464;

% shear modulus (MPa)
G = rho*Vs^2/1e6;

% Elastostatic Green's functions for displacements and stress due to 
% uniform slip on a rectangular patch (Okada 1985,1992)
% For 2D antiplane we only care about Sigma_{12} 
% y-coordinates represent source, x-coordinates represent receiver

% Sigma_{12}
s12h=@(x2,x3,y2,y3,W) G*( ...
    (x3 - y3 - 0.5*W)./( (x2-y2).^2 + (x3 - y3 - 0.5*W).^2) -...
    (x3 - y3 + 0.5*W)./( (x2-y2).^2 + (x3 - y3 + 0.5*W).^2)...
    )/2/pi;

%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%                     Set up Mesh                      %
%                                                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%
% We consider an infinitely long strike-slip fault in the x1 direction, 
% We will model variations in slip along a frictional domain Lambdaz 
% along the x3 axis

Lambdaz=40e3;                   % Frictional Domain Size
M=400;                          % Number of cells  
dz=Lambdaz/M;                   % Cell size
y3=(-0.5*M:(0.5*M-1))'*dz;       % Left/top edge of slip patch
W=ones(M,1)*dz;                 % Width of slip patch

%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%                       Kernels                        %
%                                                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%
ss.K=zeros(M,M);                        % stress kernels 
for k=1:M
    % Evaluate the stress at the center of the slip patches
    % Coordinate from source is top of patch
    % s1h(x2,x3,y2,y3,W)
    ss.K(:,k)=s12h(0,y3,0,y3(k),W(k));
    
end

%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%                    Fault Properties                  %
%                                                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%

% reference friction coefficient
ss.fo=0.6*ones(size(y3));

VWp = find(abs(y3) <= 10e3);

% Dieterich-Ruina R+S frictional parameters (velocity-weakening friction)
ss.a = 0.025*ones(size(y3));
ss.a(VWp) = 0.01;
ss.b=0.015*ones(size(y3));

% effective normal stress (MPa)
ss.sigma=50.0*ones(size(y3));

% characteristic weakening distance (m)
ss.Drs=8e-3*ones(size(y3));

% plate rate (m/s)
ss.Vpl=1e-9*ones(size(y3));

% reference slip rate (m/s)
ss.Vo=1e-6*ones(size(y3));

% Radiation damping coefficient
ss.eta = G./(2*Vs);

% Estimates of some key parameters
VWp = find(ss.a < ss.b); % VW region
% Critical nucleation size ( h* = pi/2 G b D_rs / (b-a)^2 / sigma )
hstar=min(pi/2*G*ss.Drs(VWp).*ss.b(VWp)./(ss.b(VWp)-ss.a(VWp)).^2./ss.sigma(VWp));

% Quasi-static cohesive zone ( coh0 = 9/32 G D_rs /(b*sigma) ) 
% Note that for this QD simulation the cohesive zone will not change,
% which would not be the case for a fully dynamic simulation
coh = min(9/32*pi*G*ss.Drs(VWp)./ss.b(VWp)./ss.sigma(VWp));

% Estimate of recurrence time ( T ~ 5(b-a)*sigma / G * R/Vpl ) 
Ti = 5*mean((ss.b(VWp)-ss.a(VWp)).*ss.sigma(VWp)).*0.5.*(y3(VWp(end))-y3(VWp(1)))./(G*mean(ss.Vpl(VWp)));

% Print information about discretization
fprintf('Grid size = %.2f (m)\n', dz);
fprintf('VW zone = %.2f (km)\n', (y3(VWp(end))-y3(VWp(1)))/1e3);
fprintf('Critical nucleation size = %.2f (m)\n',hstar);
fprintf('QS Cohesive zone = %.2f (m)\n',coh);
fprintf('Est. Recurrence time = %.2f (yr)\n', Ti/3.15e7);

%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%                    Numerical Solution                %
%                                                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%
% Use ode45 (Runge-Kutta 4th / 5th order accurate integration) to 
% solve the ODE time integration problem with adaptive time-steps
% yp = f(t,y)
% Y = [slip; stress; state variable; log10(slip rate / ref slip rate)]
% Degrees of Freedom
ss.dgf=4; 

% Initial conditions (start at steady-state with zero slip)
Y0=zeros(M*ss.dgf,1);   
Y0(1:ss.dgf:end)=zeros(M,1);   
Y0(2:ss.dgf:end)=max(ss.a).*ss.sigma.*asinh(ss.Vpl./ss.Vo/2.*exp((ss.fo+ss.b.*log(ss.Vo./ss.Vpl))./max(ss.a))) + ss.eta.*ss.Vpl;
Y0(3:ss.dgf:end)=ss.a./ss.b.*log(2*ss.Vo./ss.Vpl.*sinh((Y0(2:ss.dgf:end)-ss.eta.*ss.Vpl)./ss.a./ss.sigma))-ss.fo./ss.b;
Y0(4:ss.dgf:end)=log(ss.Vpl./ss.Vo);

% initialize the function handle with set constitutive parameters
yp=@(t,y) DieterichRuinaRegAging(t,y,ss);

% ODE45 Settings
% Initial step of 1e-5 seconds
% Relative tolerance of 3e-8
% [0 3e10] = simulate 3e10 seconds, 3.15e7 seconds / year
tic
options=odeset('Refine',1,'RelTol',1e-8,'InitialStep',1e-5);
[t,Y]=ode45(yp,[0 500*3.15e7],Y0,options);  
disp('Done solving');
toc

% Included in this directory are two functions myode23 and myode45 which
% are hacked versions of the matlab functions ode23 and ode45 (may be older
% than current Matlab versions). These hacks are just to output data
% periodically to disk in order to allow for larger spatial domains and 
% time periods. Otherwise the output vectors are just increasing arrays 
% that tend to exceed memory limits of the machine. Since the versions may
% be older you should check the dimensionality of Y and t, they may be
% flipped for later plottting

%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%                        Figures                       %
%                                                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%

V = ss.Vo.*exp(Y(:,4:ss.dgf:end)'); % Slip rate (m/s)
tau = Y(:,2:ss.dgf:end);            % Shear stress (MPa)
Vmax = zeros(length(t),1);          % Maximum slip rate (m/s)
Vcenter = V(floor(M/2),:);          % Slip rate at center of VW region
for ti = 1:length(t)
    Vmax(ti) = max(V(:,ti));
end

% GPS time series (displacement in x1 direction)
shistory = Y(:,1:ss.dgf:end)';

% Evolution of slip rate in the time domain
figure(1);clf;set(gcf,'name','Time evolution')

subplot(2,1,1);cla;
pcolor(t/3.15e7,y3/1e3,log10(V)), shading flat
set(gca,'YDir','reverse');
h=colorbar('Location','NorthOutside');
title(h,'Log10 Slip Velocity (m/s)')
xlabel('Time (yr)')
ylabel('Position along fault (km)');

subplot(2,1,2);cla;
plot(t/3.15e7,log10(Vmax)')
xlabel('Time (yr)')
ylabel('Maximum Velocity (m/s)')
title('Time series at the fault center')

% Plot in the time-step domain
figure(2);clf;set(gcf,'name','Evolution with time steps')
subplot(2,1,1);cla;
pcolor((1:length(t)),y3/1e3,log10(V)), shading flat
set(gca,'YDir','reverse');
h=colorbar('Location','NorthOutside');
title(h,'Log10 Slip Velocity (m/s)')
xlabel('Time Step')
ylabel('Position along fault (km)');

subplot(2,1,2);cla;
plot((1:length(t)),log10(Vmax))
xlabel('Time Step')
ylabel('Velocity (m/s) log10')
title('Evolution at the fault center')



