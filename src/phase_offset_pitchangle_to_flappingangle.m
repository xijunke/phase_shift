%% phase_offset_pitchangle_to_flappingangle
% 用于求解拍打角和扭转角的相对相位差
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(3)       % 图1——拍打角和扭转角
subplot(311)
plot(t/T,phi*180/pi,'r-',t/T,psi*180/pi,'b-','LineWidth',2)  %转换为ms 和 度数degree   *10^3   *180/pi
xlabel('\itNormalized time')
ylabel('\itAngle (°)')
legend('\it\phi(t)','\it\psi(t)')
title('拍打角和扭转角随时间的变化规律')   % 拍打角和扭转角随时间的变化规律
grid on
% axis([1.2444,4.2444,-105,105])            % ——雄蜂  t_0(0.0102)*f(122)=1.2444;
% set(gca,'XTick',(1.2444:0.2:4.2444))     % ——雄蜂  t_0(0.0102)*f(122)=1.2444;
% axis([3.2526,6.2526,-105,105])                % ——果蝇  t_0(0.0139)*f(234)=3.2526;
% set(gca,'XTick',(3.2526:0.2:6.2526))         %——果蝇  t_0(0.0139)*f(234)=3.2526;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hold on
t_phi_max=-epsilon/(2.*pi.*f)+T;                             % t_phi_max =-0.0109;
t_psi_0=(asin(-psi_0/psi_m)-zeta)/(2.*pi.*f)+T;        % t_psi_0 =1.7425e-004;
plot((t_phi_max+T/2)/T,-phi_m*180/pi,'rs',t_psi_0/T,0,'bs','LineWidth',3)
hold on
plot([t_psi_0/T,t_psi_0/T],[-100,100],'g-.',[(t_phi_max+T/2)/T,(t_phi_max+T/2)/T],[-100,100],'g-.','LineWidth',1.5)
% hold on
% plot((t_phi_max+T)/T,phi_m*180/pi,'rs',(t_psi_0+T/2)/T,0,'bs','LineWidth',4)
delta_t1=t_phi_max+T-(t_psi_0+T/2);  % delta_t1 = 0.0018;
delta_t=(pi-((epsilon+asin(-psi_0/psi_m)-zeta)))/(2.*pi.*f); % delta_t =0.0018;
delta=pi-(epsilon+asin(-psi_0/psi_m)-zeta); % delta =0.4495;
delta_deg=delta*180/pi;   % delta_deg =25.7563;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hold on
k=0;
t_phi_0=(pi/2+2*k*pi-epsilon)/(2.*pi.*f)+T;   % t_phi_0 =-0.0044;
psi_specific=psi_m.*sin(2.*pi.*f.*t_phi_0+zeta)+psi_0;  % psi_specific =-0.8767;
alpha_specific=(pi/2+psi_specific)*180/pi;    % alpha_specific =39.7662;
plot(t_phi_0/T,0,'rd',t_phi_0/T,psi_specific*180/pi,'bd','LineWidth',3)
hold on
plot([t_phi_0/T,t_phi_0/T],[0,psi_specific*180/pi],'k:','LineWidth',1.5)
hold on
t_phi_01=(pi/2+2*k*pi-epsilon)/(2.*pi.*f)+T/2+T;  % t_phi_01 =0.0085;
psi_specific1=psi_m.*sin(2.*pi.*f.*t_phi_01+zeta)+psi_0;  % psi_specific1 =0.8797;
alpha_specific1=(pi/2-psi_specific1)*180/pi;   % alpha_specific1 =39.5943;
plot(t_phi_01/T,0,'rd',t_phi_01/T,psi_specific1*180/pi,'bd','LineWidth',4)
hold on
plot([t_phi_01/T,t_phi_01/T],[0,psi_specific1*180/pi],'k:','LineWidth',1.5)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%