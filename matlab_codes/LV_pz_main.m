% LV_pz_main.m
% Toy 2-box ecosystem model
%
% Started 2004:10:25 S. Doney, WHOI
% Modif'd 2009:03:20 D. Glover to change labels and line colors
% Modif'd 2009:11:24 D. Glover to use in chap_TUT as basic Lotka-Volterra
% Modif'd 2009:11:25 D. Glover to finish up the annotations
% Modif'd 2009:11:27 D. Glover to correct labels

% define "global" paprameters and assign values
clear all
clf
global MLD p

% units 
% concentration (mmol/m^3)
% fluxes (mmol/m^3/d)

p(1)=.1;     % phytoplankton growth rate
p(2)=.4;     % zooplankton grazing rate
p(3)=.2;     % zooplankton assimilation efficiency
p(4)=0.05;   % zooplankton mortality
MLD=50;     % mixed layer depth

T=[0:0.1:400]';      % use a fixed time-scale to compare cases
n=length(T);

X0=[0.8 .174]';
%X0=[0.8 0.1]';    % Dave's original starting point

[T,X]=ode15s('LV_pz',T,X0);

%tobs=[0:5:1000]';
%xobs(:,1)=X(1:50:n,1);
%xobs(:,2)=X(1:50:n,2);

figure(1); % time series plot
subplot(221);
H1=plot(T,X(:,1),'k',T,X(:,2),'-.k');
set(H1(2),'linewidth',1);
%tobs,xobs(:,1),'bo',tobs,xobs(:,2),'rv')
ylabel('Plankton (mmol-N m^{-3})','fontsize',14);
xlabel('Time (days)','fontsize',14);
text(10,1.3,'{\bf a}','fontsize',14);
%title('Initial Model Run','fontsize',16);
legend('Phytoplankton','Zooplankton','location','northeast');
legend boxoff;

%figure(2); % phase plane plot
subplot(222);
plot(X(:,1),X(:,2),'k');hold on;
%plot(X(1,1),X(1,2),'.k','markersize',10);
%plot(X(50,1),X(50,2),'.k','markersize',10);
xlabel('Phytoplankton (mmol-N m^{-3})','fontsize',14);
ylabel('Zooplankton (mmol-N m^{-3})','fontsize',14);
text(0.23,0.38,'{\bf b}','fontsize',14);
%annotation('arrow',[X(1,1)/2 X(1,2)/0.7],[X(100,1)/2 X(100,2)/0.7]);
%quiver(X(1,1),X(1,2),(X(50,1)-X(1,1)),(X(50,2)-X(1,2)),'k','linewidth',3);
H2=quiver(X(1,1),X(1,2),(X(50,1)-X(1,1)),(X(50,2)-X(1,2)),1,'k','linewidth',2,'maxheadsize',1);
hold off
print -djpeg Long_Run.jpg
print -deps ftut_LVpzts1.eps
%save ecopz_out tobs xobs
figure(2);
N=length(X); nn=[N-1000:N];
subplot(221)
plot(T(nn),X(nn,1),T(nn),X(nn,2));
subplot(222)
plot(X(nn,1),X(nn,2))

