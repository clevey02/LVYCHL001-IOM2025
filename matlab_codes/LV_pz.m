function xdot = LV_pz(t,x)
% A function m-file for use with ODE45 to solve a set of coupled ordinary
% differential equations for a simple ecosystem model
%
% This is the simplified version, MLD is kept constant.
%
% Started 2004:10:25 S. Doney, WHOI
% Modif'd 2009:11:24 D. Glover to make a simple Lotka-Volterra (LV_pz)

global MLD p
xdot=zeros(2,1);

% p(1) phytoplankton growth rate
% p(2) zooplankton grazing rate
% p(3) zooplankton assimilation efficiency
% p(4) zooplankton mortality

xdot(1) = x(1)*(p(1)-p(2)*x(2));
xdot(2) = x(2)*(p(3)*p(2)*x(1)-p(4));

%xdot(1) = p(1)*(1. -x(1)/p(2))*x(1) - p(3)*x(1)*x(2);
%xdot(2) = p(3)*x(1)*x(2) - p(4)*x(2);
