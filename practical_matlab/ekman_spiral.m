%  ekman_spiral;  an elementary upper ocean Ekman layer
%  Solve a 1-D diffusion flow numerically.  This
%  version has rotation, and thus makes an Ekman 
%  layer.  The IC is a state of rest.  The fluid is forced
%  by an imposed stress at the top of the column (z=0)
%  The lower boundary condition is free-slip to minimize
%  the effect of finite depth. 
%
%
%  Original code by Jim Price, Oct 99
%  Modified version by Marcello Vichi, for UCT Oceanography


%% user inputs
dz = 5.;          %  [m] grid interval
L = 100.;         %  [m] water column depth  
ndays = 20.;      %  [day] days to integrate
A = 5.e-2;        %  [m2/s] the eddy diffusivity
Tauwy = 0.1;      %  [N/m2=Pa] the wind stress (y)
Tauwx = 0.;       %  [N/m2=Pa] the wind stress (x)
lat = -45.;       %  [deg] specify the latitude here
%  plot profiles every nplot steps
nplot = 25;

% constants
omega = 7.29e-5;  %  [s-1] 2pi/86400
rho0 = 1025.;     %  [kg/m3] nominal constant density of water
SECPERDAY = 86400;
w = 0.4;          % the Courant number dt*A/dz^2
                  % must be less than 0.5 for numerical stability

% grid specifications
z = 0:-dz:-L ;   %  the grid definition
nz = numel(z);
dt = w*dz^2/A;   %  the time step is derived from the Courant number
nstep = round(ndays*SECPERDAY/dt); % number of steps per day
time = zeros(nstep,1);

% derived parameters
f = 2.*omega*sin(lat*pi/180);      % the Coriolis parameter
IP = (2*pi/abs(f))/SECPERDAY;      % inertial period


%% Variables and Initial conditions
u = zeros(size(z));         %  the dependent variable, velocity
v = u;
uavg = u;
vavg = u; 
navg = 0;

% set figures
set(0,'DefaultLineLineWidth',1.4)
set(0,'DefaultTextFontSize',14)
set(0,'DefaultAxesLineWidth',1.3)
set(0,'DefaultAxesFontSize',14)

f1=figure;
ax1=subplot(2,1,1); hold on; box on
ylabel('depth [m]')
xlabel('Zonal current [m/s]')
title('Upper ocean Ekman layer, \tau_y = 0.1 N m^{-2}')
    
ax2=subplot(2,1,2);  hold on; box on
ylabel('depth [m]')
xlabel('Meridional current, [m/s]')

%%
%  begin time-stepping
ht = []; % handle for the time string in the figure
for n=1:nstep
    
    % advance time
    time(n) = (n-1)*dt;
        
    % evaluate the diffusion term
    delsqu = [0 (u(1:end-2)-2.*u(2:end-1)+u(3:end)) 0];
    delsqv = [0 (v(1:end-2)-2.*v(2:end-1)+v(3:end)) 0];
    
    % Euler forward solution with Coriolis
    u = u + w*delsqu + dt*f*v;
    v = v + w*delsqv - dt*f*u;
    
    % apply the surface BC
    u(1) = u(2)  + dz*(Tauwx/rho0/A);
    v(1) = v(2)  + dz*(Tauwy/rho0/A);
    % the bottom BC (Neumann, free-slip)
    u(nz) = u(nz-1);
    v(nz) = v(nz-1);
    
    % Diagnostics
    
    %  1) time-average (after the first day)
    if time(n)/SECPERDAY >= 1.
        uavg = uavg + u;
        vavg = vavg + v;
        navg = navg + 1;
    end
    
    
    % 2) time series of surface current
    usurf(n) = u(1);
    vsurf(n) = v(1);
    mid = round(2*nz/3) -1;
    umid(n) = u(mid);
    vmid(n) = v(mid);
    udeep(n) = u(nz-3);
    vdeep(n) = v(nz-3);

    % 3) time series of transport
    transu(n) = (sum(u) - u(1))*dz;
    transv(n) = (sum(v) - v(1))*dz;

    % 4) Energy dissipation and kinetic energy
    ds = 0.;
    kes = 0.;
    for k=2:nz-1
        shru = (u(k+1) - u(k))/dz;
        shrv = (v(k+1) - v(k))/dz;
        ds = ds - dz*A*(shru^2 + shrv^2);
        kes = kes + 0.5*dz*(u(k)^2 + v(k)^2);
    end
    diss(n) = ds;
    ke(n) = kes;
           
    % plot profile every nplot time steps
    if mod(n,nplot) == nplot - 1
        figure(f1)
        plot(ax1, u, z)
        plot(ax2, v, z)
        if (~isempty(ht)); delete(ht); end
        ht=text(-0.03,-90,['time=',num2str(time(n)/SECPERDAY,3)]);
        %dummy = input('Press any key to advance');
        drawnow
    end
    
end

% compute mean profile
uavg = uavg/navg;
vavg = vavg/navg;


%%  plot diagnostics
% Theoretical Ekman transport and Ekman depth
Ektrans = sqrt(Tauwx.^2+Tauwy.^2)/rho0/f;
Ekdepth = sqrt(2*A/abs(f));

%% time series at different depths
dtime=time/SECPERDAY;
figure
subplot(2,1,1)
plot(dtime, usurf, dtime, umid, dtime, udeep)
ylabel('East [m/s]')
grid
legend('Surf','Mid','Deep','orientation','horizontal'); legend boxoff
title('East and North currents')
subplot(2,1,2)
plot(dtime, vsurf)
hold on
plot(dtime, vmid)
plot(dtime, vdeep)
ylabel('North [m/s]')
xlabel('time [day]')
grid

%% transport
figure
plot(transu, transv,'.')
hold on
Tx = mean(transu);
Ty = mean(transv);
plot([0 Tx], [0 Ty], 'g-o')
plot([0 Ektrans], [0 0],'r-+')
xlabel('U transport [m2/s] ')
ylabel('V transport [m2/s]')
% axis([-0.1 0.15 -0.25 0.])
axis('equal')
title('Transport with time, average and theoretical value')
grid on
legend('istantaneous','average','theoretical')
%% mean velocity profiles
figure
plot(uavg, z,'r',vavg, z,'b','linewidth',2)
xlabel('East and North currents [m/s]')
ylabel('Depth [m]')
title('Mean current profiles')
hold on
h=line([-0.1 0.2],-Ekdepth*[1 1]);
set(h,'linestyle','--')
legend('E','N','Ekman depth','location','southeast'); legend boxoff

%% hodograph
figure
zeroz=zeros(size(uavg));
quiver(0,0,0,0.06,'r')
hold on
quiver(zeroz,zeroz,uavg,vavg)
xlabel('East component [m/s]')
ylabel('North component, [m/s]')
title('Hodograph, dz = 2.5 m')

%% spiral in 3D
figure
h1=quiver3(zeroz,zeroz,z,uavg,vavg,zeroz,0);
set(h1,'LineWidth',1.3)
hold on
x3 = [0 0 0]; y3 = [0 0 0.1]; z3 = [-L 0 0];
hold on
h2 =  plot3(x3, y3, z3, 'r-' );
set(h2,'LineWidth',1.5)
xlabel('East component [m/s]')
ylabel('North component, [m/s]')
zlabel('Depth [m]')
grid on
 
