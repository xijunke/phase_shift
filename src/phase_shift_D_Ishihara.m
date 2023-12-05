%% phase_shift_D_Ishihara
% single degree of freedom mass每spring每dashpot system of the wing for the theoretical analysis.
clear all;clc;
% m_w=0.000245/2*10^-3; % kg
% m_w=Rho_s*c*(0.2*h_s+0.8*h); % m_w=0.000166g; 
m_w=0.000166*10^-3; % kg
R=1.27*10^-2; % m
% S=0.59; % cm^2
% c=S/(2*R); 
c=0.23*10^-2;    % c=0.23cm;
PHI=123; % degree
f_cranefly=45.5;     % Hz
%% the theoretical analytical solution of 1 DOF  mass每spring每dashpot system
k_s=1.9*10^-7;    % gcm^2(s^2rad)^-1
k_s_delta=4*k_s/(c^2);
k=k_s_delta;
% the natural frequency 
f_n=(1/(2*pi))*sqrt(k_s_delta/m_w); % f_n =172.3578;   f_n =148.0625; 
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f=linspace(0,2*f_n,10000);
% three representative values of the damping ratios 汎=0.5 (under damping), 1 (critical damping) and 2 (over damping).
% zeta=0.5; 
% zeta=1;
zeta=2;
% kappa=1./(sqrt((1-(f/f_n).^2).^2+(2*zeta*(f/f_n)).^2));
% b=pi/2-atan(2*zeta*(f/f_n)./(1-(f/f_n).^2)); % phase shift: b
b=pi/2-atan2(2*zeta*(f/f_n),(1-(f/f_n).^2)); % phase shift: b
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%The relationship between the frequency ratio f/fn and the phase shift b 
% for the damping ratios 汎=0.5 (red line), 1 (blue line) and 2 (black line).
figure(1)
hold on
plot(f/f_n,b*180/pi,'b-','LineWidth',2)
xlabel('\itFrequency ratio f/f_n','Fontsize',20,'FontName','Times','FontWeight','Bold')
ylabel('\itPhase shift b (deg.)','Fontsize',20,'FontName','Times','FontWeight','Bold')
% legend('b vs f/f_n')
title('\itThe relationship between the frequency ratio f/f_n and the phase shift b')
box on
axis([0,2,-60,90])
set(gca,'XTick',(0:0.25:2),'LineStyle','-','LineWidth',1.5,'FontSize',16,'FontName','Times','FontWeight','Bold')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% the theoretical analytical solution
f_pitch=f_cranefly; % Hz
% f_pitch=91; % Hz
% f_pitch=120; % Hz
% f_pitch=f_n; % Hz % f_n =148.0625; 
f_ratio=f_pitch/f_n;  % f_ratio=0.3073;
A_0=(1/2)*R*PHI*(pi/180); % A_0=1.36cm;
F_0=m_w*A_0/2*(2*pi*f_pitch).^2;
P_scalingfactor=1.5;
delta_0=P_scalingfactor*F_0/k
%%%%%%%%%%%%%%%%%%%%%%%%%
% the damping ratio
% zeta=c_f/(2*sqrt(m_w*k_s_delta)); % c_f: is the damping force coefficient due to the fluid;
% three representative values of the damping ratios 汎=0.5 (under damping), 1 (critical damping) and 2 (over damping).
% zeta=0.5; 
% zeta=1;
% zeta=2;
%%%%%%%%%%%%%%%%%%%%%%%%%
kappa=1./(sqrt((1-(f_ratio).^2).^2+(2*zeta*(f_ratio)).^2)) % amplitude ratio or magnified factor
% b=pi/2-atan(2*zeta*(f_ratio)./(1-(f_ratio).^2)); % phase shift: b
b=pi/2-atan2(2*zeta*(f_ratio),(1-(f_ratio).^2)) % phase shift: b
t=linspace(0,4/f_pitch,10000);
delta_t=kappa.*delta_0.*sin(2*pi.*f_pitch.*t+b); % the theoretical analytical solution
alpha=delta_t/(c/2);
%%%%%%%%%%%%%%%%%%%%%%%%%
figure(2)
plot(t*f_pitch,alpha*180/pi,'k-','LineWidth',2)
xlabel('\itNormalized time','Fontsize',20,'FontName','Times','FontWeight','Bold')
ylabel('\itPitch angle \alpha (deg.)','Fontsize',20,'FontName','Times','FontWeight','Bold')
title('\itTime histories of the simulated passive pitching motion')
box on
% axis([0,max(t),-100,100])
set(gca,'LineStyle','-','LineWidth',1.5,'FontSize',16,'FontName','Times','FontWeight','Bold')




