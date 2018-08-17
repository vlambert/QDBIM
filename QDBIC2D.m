%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%      2D Quasi-dynamic simulator for evaluating       %
%   the slip history on a fault in antiplane strain    %
%   governed by rate- and state-dependent friction     %
%                                                      %
%                Valere Lambert, 2018                  %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%
clear all;close all;
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

% Note that these solutions take into account a free-surface using
% the method of images

% Sigma_{12}
s12h=@(x2,x3,y2,y3,W) G*( ...
    -(x3-y3)./((x2-y2).^2+(x3-y3).^2)+(x3+y3)./((x2-y2).^2+(x3+y3).^2) ...
    +(x3-y3-W)./((x2-y2).^2+(x3-y3-W).^2)-(x3+y3+W)./((x2-y2).^2+(x3+y3+W).^2) ...
    )/2/pi;


% Displacement kernels for fault slip
u1h=@(x2,x3,y2,y3,W) ...
    (+atan((x3-y3)./(x2-y2))-atan((x3+y3)./(x2-y2)) ...
     -atan((x3-y3-W)./(x2-y2))+atan((x3+y3+W)./(x2-y2)) ...
    )/2/pi;

%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%                     Set up Mesh                      %
%                                                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%
% We consider an infinitely long vertical strike-slip fault 
% in the x1 direction, extending from the surface to a fixed  
% depth, with slip only varying in the x3 = z direction (down-dip)

Lambdaz=40e3;         % Depth extent of the frictional domain
M=400;                % Number of cells  
dz=Lambdaz/M;         % Cell size
y3=(0:M-1)'*dz;       % Top of slip patch
W=ones(M,1)*dz;       % Down-dip width of slip patch

% Surface points (virtual GPS receivers crossing the fault in x2)
nsta = 100;
eps = 1e-6;
x2GPSR =  (200e3*tan(eps+(0:nsta/2)'*pi/(2*nsta)));
x2GPSL = (-200e3*tan(eps+(0:nsta/2)'*pi/(2*nsta)));
x2GPS  = [flipud(x2GPSL);x2GPSR];

%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%                       Kernels                        %
%                                                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%
ss.K=zeros(M,M);                        % stress kernels 
ss.ku1=zeros(length(x2GPS),M);          % displacement kernels for GPS
for k=1:M
    % Evaluate the stress at the center of the slip patches
    % Coordinate from source is top of patch
    % s1h(x2,x3,y2,y3,W)
    ss.K(:,k)=s12h(0,y3+dz/2,0,y3(k),W(k));
    
    % Displacement kernels
    % u1h=@(x2,x3,y2,y3,W)
    ss.ku1(:,k)=u1h(x2GPS,0,0,y3(k),W(k));
end

%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%                    Fault Properties                  %
%                                                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%

% reference friction coefficient
ss.fo=0.6*ones(size(y3));

% Dieterich-Ruina R+S frictional parameters (velocity-weakening friction)
ss.a=1e-2+Ramp((y3-15e3)/3e3)*(0.025-0.01);
ss.b=0.015*ones(size(y3));

% effective normal stress (MPa)
ss.sigma=50.0*ones(size(y3));

% characteristic weakening distance (m)
ss.L=8e-3*ones(size(y3));

% plate rate (m/s)
ss.Vpl=1e-9*ones(size(y3));

% reference slip rate (m/s)
ss.Vo=1e-6*ones(size(y3));

% Radiation damping coefficient
ss.eta = G./(2*Vs);

% Estimates of some key parameters
VWp = find(ss.a < ss.b); % VW region
% Critical nucleation size ( h* = pi/2 GbL / (b-a)^2 / sigma )
hstar=min(pi/2*G*ss.L(VWp).*ss.b(VWp)./(ss.b(VWp)-ss.a(VWp)).^2./ss.sigma(VWp));

% Quasi-static cohesive zone ( coh0 = 9/32 GL/(b*sigma) ) 
% Note that for this QD simulation the cohesive zone will not change
% but this would not be the case for a truly dynamic simulation
coh = min(9/32*pi*G*ss.L(VWp)./ss.b(VWp)./ss.sigma(VWp));

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
% Use ode45 (Runge-Kutta 4th 5th order accurate integration) to 
% solve the ODE time integration problem with adaptive time-steps
% yp = f(t,y)
% Y = [slip; stress; state variable; slip rate]
% Degrees of Freedom
ss.dgf=4; 

% Initial conditions (start at steady-state with zero slip)
Y0=zeros(M*ss.dgf,1);   
Y0(1:ss.dgf:end)=zeros(M,1);   
Y0(2:ss.dgf:end)=ss.a.*ss.sigma.*asinh(ss.Vpl./ss.Vo/2.*exp((ss.fo+ss.b.*log(ss.Vo./ss.Vpl))./ss.a));
Y0(3:ss.dgf:end)=ss.a./ss.b.*log(2*ss.Vo./ss.Vpl.*sinh((Y0(2:ss.dgf:end))./ss.a./ss.sigma))-ss.fo./ss.b;
Y0(4:ss.dgf:end) =ss.Vpl;

% initialize the function handle with set constitutive parameters
yp=@(t,y) odeRegDRaging(t,y,ss);

% ODE45 Settings
% Initial step of 1e-5 seconds
% Relative tolerance of 3e-8
% [0 3e10] = simulate 3e10 seconds, 3.15e7 seconds / year
tic
options=odeset('Refine',1,'RelTol',3e-8,'InitialStep',1e-5);
[t,Y]=ode45(yp,[0 500*3.15e7],Y0,options);  
disp('Done solving');
toc

%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%                        Figures                       %
%                                                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%

V = Y(:,4:ss.dgf:end)';       % Slip rate (m/s)
tau = Y(:,2:ss.dgf:end)';     % Shear stress (MPa)
Vmax = zeros(length(t),1);    % Maximum slip rate (m/s)
Vcenter = V(floor(M/2),:);    % Slip rate at center of VW region
for ti = 1:length(t)
    Vmax(ti) = max(V(:,ti));
end

% GPS time series (displacement in x1 direction)
shistory = Y(:,1:ss.dgf:end)';
UGPS = ss.ku1*shistory;

% Evolution of slip rate in the time domain
figure(1);clf;set(gcf,'name','Time evolution')

subplot(2,1,1);cla;
pcolor(t/3.15e7,y3/1e3,log10(V)), shading flat
set(gca,'YDir','reverse');
h=colorbar('Location','NorthOutside');
title(h,'Log10 Slip Velocity (m/s)')
xlabel('Time (yr)')
ylabel('Depth (km)');

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
ylabel('Depth (km)');

subplot(2,1,2);cla;
plot((1:length(t)),log10(Vmax))
xlabel('Time Step')
ylabel('Velocity (m/s) log10')
title('Evolution at the fault center')
%%
% GPS time series
figure(3);clf
set(gcf,'Position',[50 50 590 625])
dt = 2.628e6;            % 1 month sampling
t0 = 150*3.15e7;         % Initial time for record
MaxRecordT = 500*3.15e7; % 100 years time series
RecordT = t0 + (0:dt:MaxRecordT)';
tindx = zeros(length(RecordT),1);
for tstep = 1:length(RecordT)
    tindx(tstep) = find(abs(t - RecordT(tstep)) == min(abs(t - RecordT(tstep))),1,'first');
end
colorspec = jet(length(RecordT));
subplot(3,2,[1 2])
for tstep=1:length(RecordT)
    plot(x2GPS/1e3,UGPS(:,tindx(tstep)) - UGPS(:,tindx(1)),'color',colorspec(tstep,:)); hold on;
end
cb=colorbar;
colormap(cb,colorspec)
caxis([1 length(RecordT)])
ylabel(cb,'Time Steps');
box on; grid on;
xlabel('Distance across fault (km)');
ylabel('Displacement (m)')
title('100 year time series')
xlim([min(x2GPS/1e3) max(x2GPS/1e3)])

% Plot time series for specific points
[dist1,GPS1] = min(abs(x2GPS - 5e3));   % GPS station around 5 km from fault
[dist2,GPS2] = min(abs(x2GPS - 150e3)); % GPS station around 100 km from fault
subplot(3,2,[3 4])
plot(t/3.15e7,UGPS(GPS1,:),'k');
box on; 
xlabel('Time (yrs)');
ylabel('Displacement (m)');
title(sprintf('GPS at x = %.2f km',x2GPS(GPS1)));
xlim([min(t) max(t)]./3.15e7)

subplot(3,2,[5 6])
plot(t/3.15e7,UGPS(GPS2,:),'k');
box on; 
xlabel('Time (yrs)');
ylabel('Displacement (m)');
title(sprintf('GPS at x = %.2f km',x2GPS(GPS2)));
xlim([min(t) max(t)]./3.15e7)



