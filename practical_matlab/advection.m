%%
% Linear advection in 1D
%
% Resolution of dT/dt =  -c dT/dx
% using an explicit forward time step 
%
%

% problem parameters
c=1;			% Advection speed [m s-1]
dx=500;		% X resolution [m]
dt=100;			% time step [s]
t0=0.;          % initial time
tmax=0.1;		% duration of simulation [days]
xmax=10.e3;	% Length of the basin [m]
%
% Test of the CFL criterion
%
disp(['DT = ',num2str(dt),' s'])
disp(['DX/DT = ',num2str(dx/dt)])
disp(['MAX SPEED = ',num2str(c),' m/s'])
disp(['COURANT NUMBER = ',num2str(c*dt/dx)])
disp(' ')
%
% Grid definition
%
NT = floor(tmax*86400/dt)+1; % number of steps (starts from 0)
x=-xmax:dx:xmax; 
IM=numel(x); % number of grid points
disp(['grid size : ',num2str(IM)])
disp(' ')

%
% Initial condition
%
T0 = 15; % baseline temperature [degC]
sigma = 2000.; % width [m]
A = 5.; % amplitude [degC]
Gaussian = @(x) 15+A*exp(-(x/sigma).^2);
T = Gaussian(x);
%
% First plot
%
figure
plot(x/1000,T,'r')
xlabel('X [km]')
ylabel('T [^oC]')
axis([-xmax/1000 xmax/1000 13 22]) 
hold on
%
% Main time loop
%
for n=1:NT    
    %
    % Spatial Advection term
    %
    %rhs=-c*(T(3:end)-T(2:end-1))/dx;	% first order "downstream" scheme
    %rhs=-c*(T(2:end-1)-T(1:end-2))/dx;	% first order upstream scheme
    rhs=-c*(T(3:end)-T(1:end-2))/(2*dx);% second order centered scheme
    
    %
    % Euler time stepping
    %
    T(2:end-1)=T(2:end-1)+dt*rhs;                 % Forward scheme
    %T(2:end-1)=0.5*(T(1:end-2)+T(3:end))+dt*rhs; % Lax scheme
    
    %
    % Boundary conditions (Dirichlet(=Clamped))
    %
    T(1)=15;
    T(end)=15;
    
    %
    % Plot numerical solution
    %
    h1=plot(x/1000,T,'m');
    
    %
    % Plot analytical solution
    %
    t=t0+n*dt;
    h2=plot(x/1000,Gaussian(x-c*t),'--','color',[0.7 0.7 0.7]);
    legend([h1 h2],{'Model','Analytical'})
    
    % 
    % Get the frame for the movie
    % 
    Mov(n)=getframe;
end
hold off

%
% Play the movie
%
movie(Mov,2)
