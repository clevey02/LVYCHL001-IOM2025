%
% Molecular diffusion in 1D
%
%
clear all
close all

% constants
SECPERDAY = 86400.;
M2KM = 1000;

% Problem parameters
kappa=1.e3;		% Diffusion parameter [m2 s-1]
dx=500;		% X resolution [m]
dt=100;			% time step [s]
t0=0.;          % initial time
tmax=1.;		% duration of simulation [days]
xmax=10.e3;	% Length of the basin [m]

%
% Test of the CFL criterion
%
disp(['DT = ',num2str(dt)])
disp(['DX/DT = ',num2str(dx/dt)])
disp(['COURANT NUMBER = ',num2str(2*kappa*dt/dx^2)])
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
N0 = 20.; % baseline concentration [mmol/m3]
sigma = 2000.; % width [m]
A = 100.; % amplitude [mmo/m3]
Gaussian = @(x) N0+A*exp(-(x/sigma).^2);
N = Gaussian(x);

%
% First plot
%
figure
plot(x/M2KM,N,'r')
xlabel('X [km]')
ylabel('DIN [mmol/m3]')
%axis([-xmax/M2KM xmax/M2KM 0.8*N0 1.1*(N0+A)]) 
hold on

%
% Main computing loop
%
h1 = [];
nstep=0;

for n=1:NT
    delete(h1)
    %
    % Diffusion term
    %
    rhs=kappa*(N(1:end-2)-2*N(2:end-1)+N(3:end))./(dx^2);
    
    %
    % Euler time stepping
    %
    N(2:end-1)=N(2:end-1)+dt*rhs;
    
    %
    % Boundary conditions (Dirichlet)
    %
    N(1)=N0;
    N(end)=N0;
    
    %
    % Figure
    %
    h1=plot(x/M2KM,N,'b');
    if (mod(n,5)==0)
        nstep = nstep+1;
        Mov(nstep) = getframe;
    end
end

% Play the movie
hold off
movie(Mov,2)
