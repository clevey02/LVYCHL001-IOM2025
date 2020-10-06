%%
% Linear advection in 1D
% Analytical solution
%
%
clear all

% problem parameters
c=1;			% Advection speed [m s-1]
dx=500;		% X resolution [m]
dt=200;			% time step [s]
t0=0.;          % initial time
tmax=0.1;		% duration of simulation [days]
xmax=10.e3;	% Length of the basin [m]
%
% Grid definition
%
NT = floor(tmax*86400/dt)+1; % number of steps (starts from 0)
t=zeros(NT,1);
x=-xmax:dx:xmax; 
IM=numel(x); % number of grid points
disp(['grid size : ',num2str(IM)])
disp(['number of timesteps : ',num2str(NT)])
disp(' ')

%
% Array to store the analytical solution
%
T = zeros(NT,IM);

%
% Initial conditions
%
T0 = 15; % baseline temperature [degC]
sigma = 2000.; % width [m]
A = 5.; % amplitude [degC]
Gaussian = @(x) T0+A*exp(-(x/sigma).^2);
T(1,:) = Gaussian(x);
t(1) = t0;

%
% First plot
%
figure
plot(x/1000,T(1,:),'r')
xlabel('X [km]')
ylabel('T [^oC]')
axis([-xmax/1000 xmax/1000 13 22]) 
hold on
%
% Main time loop to show the analytical solution
%
for n=2:NT    
    
    %
    % Plot analytical solution
    %
    t(n)=t0+(n-1)*dt;
    T(n,:) = Gaussian(x-c*t(n));
    h2=plot(x/1000,T(n,:),'-','color',[0.7 0.7 0.7]);
    
    % 
    % Get the frame for the movie
    % 
    Mov(n-1)=getframe;
end
hold off

%
% Play the movie
%
movie(Mov,2)
