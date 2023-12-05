function Euler_Motion_Eq6_fruitfly_4itemTorq2_1
% ��4�ֻ��Ƶ���������: % M_x=M_xtrans+M_xrd+M_xam+M_hinge;��������Ҳ���
% ����ȫ��
% Euler_Motion_Eq����ŷ���˶�ѧ���̵����:  I_xx*ddpsi=M_x-I_xz*ddphi.*cos(x(1))+1/2*I_xx*dphi.^2.*sin(2*x(1));
% ������С���2014��6��18��,22:51:12
% ��ó�ʼֵ����2014��6��20��,13:01:02�����޸���ʱ����
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % ��һ����������������Ŷ
% Y_rcpnd=COP_Ycpnd2_fruitfly(alpha);  % �����ת���������ء������ú�����⾻ѹ�ĵ�������λ��Y_rcpnd
% Y_rcpnd_trans=COP_Ycpnd2_TransCirc(alpha);  % ���ú���COP_Ycpnd2_fruitfly��⾻ѹ�ĵ�������λ��Y_rcpnd; % ��������
% Y_rcpnd_rot=COP_Ycpnd2_RotCirc(alpha);        % ѹ�ķֲ�����Dickinson����
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ŀ��: ����psi; dpsi; ddpsi����������ʽ�����Ĵ��phi��һ�׵�dphi�Ͷ��׵�ddphi; ����Ťת��������M_x; 
% clear all; clc;
tic;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ���ú�����ò��������������ϵ���ĺ���
% wing_para_output=zeros(1,16);
% wing_para_output=[F_ndTrans,Coeff_liftdragF_N,M_xaercoeff,I1y,Z_rnd, M_xrdcoeff,...
%     F_ndRot,F_yrotcoeff,M_xRotcoeff, I2y,...
%     I_xzam,I_xxam, I5z,I6z,I7y,M_zrdcoeff];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
wing_para=wing_shape_fruitfly_nine();  %���ú���wing_shape_fruitfly;  % size(wing_para) ���������Ҫע��Ŷ������������ʹ��
% R_wingeff=wing_para(1,1);             % R_wingeff=3.004;  % mm 
% C_avereff=wing_para(1,2);              % C_avereff=0.8854;  % mm  
% F_nd=wing_para(1,3);                     % F_nd=0.46392;  ������, ���ٻ���λΪmm^4
% Z_rnd=wing_para(1,4);                   % Z_rnd=0.14016;  ������, ���ٻ���λΪmm
I_xzam=wing_para(1,5);                 % I_xzam = -0.001245    % ��λ�� mg.mm^2
I_xxam=wing_para(1,6);                 % I_xxam = 0.0002508    % ��λ�� mg.mm^2
% I3=wing_para(1,7);                     % I3=0.74851        % ������,���ٻ���λ��mm^4;
% I3z=wing_para(1,8);                   % I3z=0.0050945   % ��λ�� mg.mm
% I4z=wing_para(1,9);                   % I4z=-0.0005346  % ��λ�� mg.mm
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
wing_para=wing_shape_fruitfly_sixteen_1();    %���ú���wing_shape_fruitfly;  % size(wing_para) 
% wing_para=wing_shape_fruitfly_sixteen_2(); 
% wing_para=wing_shape_fruitfly_sixteen_good();
%%����Ťת��ǰ��%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% R_wing=3.004;
% C_aver=0.8854;
% xr0=0.3289;
% C_maxyaxis=0.025;
% C_maxyaxis=0.1;
% C_maxyaxis=0.2024; %�����Ťת���λ�ò����˳����òѧ�Ż����Ťת���λ��   
% C_maxyaxis=0.25; 
% C_maxyaxis=0.356737; 
% C_maxyaxis=0.36;
% wing_para=wing_shape_fruitfly_sixteen_good2(R_wing,C_aver,xr0,C_maxyaxis); % ע����õĺ�����wing_shape_fruitfly_sixteen_good2 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (1) ƽ�����������������ز���
% F_ndTrans=wing_para(1,1);             % F_ndTrans=0.46392;  ������, ���ٻ���λΪmm^4
% Coeff_liftdragF_N=wing_para(1,2);  % Coeff_liftdragF_N=0.0068201;  %��λ��mg*mm  %ƽ������������
M_xaercoeff=wing_para(1,3);              % M_xaercoeff=0.0060;   %��λ��: mg.mm^2   % ƽ�������������ز��������Ƴ�ƽ���µ�չ����
% I1z=wing_para(1,4);                         % I1y=0.0162    % ��λ�� mg.mm^2            % ƽ�������������ز��������Ƴ�ƽ���µ�������
% Z_rnd=wing_para(1,5)                     % Z_rnd=0.14016;  ������, ���ٻ���λΪmm
M_xrdcoeff=wing_para(1,6);               % M_xrdcoeff =1.5848e-004;  % ��λ��mg.mm^2 %ת�������������ز������Ƴ�ƽ���µ�չ����
% (2) ת�����������������ز���
% F_ndRot=wing_para(1,7);                % F_ndRot=0.74851;  ������, ���ٻ���λΪmm^4
% F_yrotcoeff=wing_para(1,8);           % F_yrotcoeff =0.0032433;  % ��λ�� mg.mm   % ת������������
M_xRotcoeff=wing_para(1,9);             % M_xRotcoeff=0.0029;   % ��λ�� mg.mm^2   % ת��������������ϵ�������Ƴ�ƽ���µ�չ����
% I2z=wing_para(1,10);                          % I2y=0.0069;        % ��λ�� mg.mm^2          % ת�������������ز��������Ƴ�ƽ���µ�������
% (3) �����������������ز���
%%%%%%%%%%%%%%%%%%%%%%%%
% I_xzam = -0.001245    % ��λ�� mg.mm^2
% I_xxam = 0.0002508    % ��λ�� mg.mm^2
% k_hinge=1.11*10^(-3);          % psi_max =70.4552 @ k_am=1;  % �����Ťת�Ƿ�ֵ��ʵ��Ťת�Ƿ�ֵ
%%%%%%%%%%%%%%%%%%%%%%%%
% I_xzam=wing_para(1,11);                    % I_xzam =0.002     % ��λ�� mg.mm^2  % �������������ز��������Ƴ�ƽ���µ�չ����
% I_xxam=wing_para(1,12);                    % I_xxam =0.00062892  % ��λ�� mg.mm^2  % �������������ز��������Ƴ�ƽ���µ�չ����
% I5y=wing_para(1,13);                      % I5z=0.0050945   % ��λ�� mg.mm       % ������������������������
% I6y=wing_para(1,14);                      % I6z=0.0011         % ��λ�� mg.mm       % ������������������������
% I7z=wing_para(1,15);                          % I7y=0.0109;        % ��λ�� mg.mm^2   % �������������ز��������Ƴ�ƽ���µ�������
% M_zrdcoeff=wing_para(1,16);             % M_zrdcoeff=0.0012; % ��λ�� mg.mm % ת�������������ز������Ƴ�ƽ���µ�������
% C_max_LtoT=wing_para(1,17);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% I5y=wing_para(1,13);                       % I5z=0.0050945   % ��λ�� mg.mm       % ������I3yӦ�ø�ΪI5y
% I6y=wing_para(1,14);                       % I6z=0.0011         % ��λ�� mg.mm       % ������I4yӦ�ø�ΪI6y
% I7z=wing_para(1,15);                       % I7y=0.0109;        % ��λ�� mg.mm^2   % ������I5zӦ�ø�ΪI7z
% I1z=wing_para(1,4);                         % I1y=0.0162        % ��λ�� mg.mm^2    % ������I7zӦ�ø�ΪI1z
% I2z=wing_para(1,10);                       % I2y=0.0069;        % ��λ�� mg.mm^2    % ������I6zӦ�ø�ΪI2z    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (1) ƽ��������������������ϵ��
k_xaero=1;   % ƽ���������ز���û��
M_xaercoeff=k_xaero*M_xaercoeff;    % mg.mm^2;  % M_xaercoeff =6.0385e-003
%% (2) ת��������������������ϵ��
k_xRot=1;  
C_R=1.55;    % ��ϵ��������-XXXXXXXX
M_xRotcoeff=k_xRot*C_R*M_xRotcoeff;
%% (3) ת��Ťת�������ء���ת��������������ϵ��
%  C_RD=1; 
% C_RD=2; 
C_RD=5.0;     %ת������������ϵ��:C_RD=C_rd��(3,6); ����C_rd=C_dmax=3.4;
% C_RD=6.0;
% C_RD=8;  
M_xrdcoeff=C_RD*M_xrdcoeff;  % mg.mm^2;  % M_xrdcoeff =7.9240e-003;
%% (4) ������Ťת��������ϵ��
% k_am=0; % ������
% k_am=0.1; % ������
% k_am=0.2; % ������
% k_am=0.3;      % lambda =2.1588;
% k_am=0.35;    % lambda =2.0166;
% k_am=0.4;  % ����
% k_am=0.6;  % lambda =1.5761;
k_am=1; % ������
%% (5) Ťת�����ظ����ء����ն�
L_h1=70e-006;     % L_h2=175e-006;   % length of wing hinge  \um�����Ѿ����㵽m
W_h=1.8e-003;                                      % width of wing hinge  \mm�����Ѿ����㵽m
t_h=7.6e-006;                                        % thickness of wing hinge \um�����Ѿ����㵽m
E_h=2.5e009;                                         % lmodulus of elasticity for wing hinge  \Gpa ע�ⵥλ
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% k_hinge=0.865;    % psi_max =70.5348 @ k_am=1; % �����Ťת�Ƿ�ֵ��ʵ��Ťת�Ƿ�ֵ
% k_hinge=1.03;      % psi_max =72.1409 @ k_am=1;
% k_hinge=0.6*1.0;      % psi_max =72.1409 @ k_am=1; % K_hinge=1.411uN.mm; f_n =380.5327;lambda =2.0166;
% k_hinge=7.3*1.08;      % psi_max =72.1409 @ k_am=1;
k_hinge=0.97;        % psi_max =70.6746@ k_am=1;
% k_hinge=1.1;        % psi_max =70.6746@ k_am=1;
% k_hinge=1.11;          % psi_max =70.4552 @ k_am=1;  % �����Ťת�Ƿ�ֵ��ʵ��Ťת�Ƿ�ֵ
% k_hinge=1.12;      % psi_max =70.1635 @ k_am=1;
% k_hinge=5;           % psi_max =70.1635 @ k_am=1;
% k_hinge=1.75;  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% k_h=k_hinge*E_h*t_h^3*W_h/(12*L_h1);  % ����k_hΪrotational stiffness of the passive hinge
% k_h=0.1*k_hinge*E_h*t_h^3*W_h/(12*L_h1)*10^9  % N/m^2*m^3*m/m=N.m=kg*m/s^2.m=mg.mm/s^2.mm:[10^12]=10^9uN.mm
k_h=1*k_hinge*E_h*t_h^3*W_h/(12*L_h1)*10^9 
% k_h=2.05*k_hinge*E_h*t_h^3*W_h/(12*L_h1)*10^9 % k_h =2.8925e+003; f_n =337.2088; lambda = 1.7870;
% k_h=1.5*k_hinge*E_h*t_h^3*W_h/(12*L_h1)*10^9 
% k_h =2.2811e+003;  % ��Ϊ�����k_xh�ĵ�λ����ΪuN.um/rad % f_n =334.0489; lambda =1.7703;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (6) ���ת������
% m_wing=2.4*10^-9;                                             % kg
% �밴������-2003-AMS-Sun Mao & 2007-JFM-Berman_Wang ZJ_Fruitfly_�ο�
% I_wing = 0.80 �� 10^(-8) (g cm^2)  =  0.8*�� 10^(-3) mg.mm2 �����µ�λ
% ע�⣺��������������Ťת��Ϊ����ͳ�������ʱ������������
% I_xx=0.000267*10^(-12);      % inertia moment of wing    \ mg.mm^2����*10^(-12) kg.m^2  % UG ���
% I_xz=-0.000594*10^(-12);     % inertia moment of wing    \ mg.mm^2����*10^(-12) kg.m^2  % UG ���
% I_xx=0.1*0.000267;      % inertia moment of wing    \ mg.mm^2  % ����
% I_xz=-0.1*0.000594;     % inertia moment of wing    \ mg.mm^2  % ����
% I_xx=4.0*0.000267;      % inertia moment of wing    \ mg.mm^2  % ����
% I_xz=-4.0*0.000594;     % inertia moment of wing    \ mg.mm^2  % ����
I_xx=0.000267;      % inertia moment of wing    \ mg.mm^2
I_xz=-0.000594;     % inertia moment of wing    \ mg.mm^2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % mass_properties_of_fruit_fly
% % 2015-PRE-Wing-pitch modulation in maneuvering fruit flies
% % is explained by an interplay between aerodynamics and a torsional spring-Beatus_PRE_2015
% R=2.1;    % mm
% a=0.7;     % mm
% b=a/50;  % mm
% m=0.03;  % mg
% I_xx=a^2*3*m/10;             % I_11 =0.0044; % mg.mm^2     % pitch rotation axis 
% I_xz=(5*a*R/6)*3*m/10;     % I_12 = 0.0110; % mg.mm^2    % inertia product
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% I_den=I_xx+k_am*I_xxam;     % ��λ��: mg.mm^2;  ע��������������Ϊ����������domega_x=ddpsi;  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (7) ��������������ϵ�� & ������ϵ��
m_wing=2.4*10^-3;                                 % mg  namely: m_wing=2.4ug
g=9.821*10^3;                                         % m/s^2=10^3mm/s^2
% xr=0.3289;                                            %  \mm       % x-root offset 
% x_com=xr+1.920243385;                     %  \mm       % ���ĵ�չ������
% % z_com=-0.149785466+0.636=0.486215; % \mm ��Ťת����������
z_com=0.149785466;                               %  \mm       % ��Ťת����������
d_com=z_com;                                         % mm
k_inert=0;
% k_inert=1;
%%%%%%%%%%%%%%%%%%%%%%%%%
I_den=I_xx+k_am*I_xxam-k_inert*z_com^2*m_wing;  
%% Ƶ�ʱ� %%%%%%%%%%%%%%%%%%%%%%%
f1=188.7;  T=1/f1; 
f_n=sqrt(k_h/I_den)/(2*pi)
lambda=f_n/f1   % f_n =334.0489; lambda =1.7703;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ODE������΢�緽��
% t=linspace(0.0052824335,0.0052824335+3*T,20000);  % t_steady1   %���
psi_exp0 =-6.5669e-010;       % ʵ��ֵ
% dpsi_exp0 =1.8709e+003;     % ʵ��ֵ
dpsi_exp0 =-1.8709e+003;     % ʵ��ֵ   -
initcond=[psi_exp0;dpsi_exp0];  % ��ʼֵ
t_range=[0.0052824335,0.0052824335+5*T];
% t_range=[0.0052824335,0.0052824335+5*T];
% options=odeset('RelTol',1e-6,'AbsTol',[1e-6 1e-6]);
% [t,x]=ode113(@passive_rot,tspan,initcond,options,M_xaercoeff,M_xrdcoeff,k_am,k_h,I_xx,I_xz,I_den,I_xyam); 
% [t,x]=ode45(@passive_rot,tspan,initcond,options,M_xaercoeff,M_xrdcoeff,k_am,k_h,I_xx,I_xz,I_den,I_xzam,m_wing,g,d_com); 
% [t,x]=ode45(@passive_rot,[0,3*T],initcond,options,M_xaercoeff,M_xrdcoeff,k_am,k_h,I_xx,I_xz,I_den,I_xzam); 
[t,x]=ode45(@passive_rot,t_range,initcond,[],M_xaercoeff,M_xrdcoeff,k_am,k_inert,M_xRotcoeff,k_h,I_xx,I_xz,I_den,I_xzam);
psi_max=max(x(:,1)*180/pi)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (A) Euler_Motion_Eq6_fruitfly�������ֵ���㱻��Ťת�Ǻͽ��ٶ�
% ����Ťת�ǣ����ٶȺ͹���������ѱ�������ʾ�ԱȽ��
psi_sim=x(:,1)*180/pi;   % ����Ťת��
dpsi_sim=x(:,2);             % ����Ťת���ٶ�
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ���Ǽ��㡪��������ʾ
% alpha=(pi/2*sign(x(:,1))-x(:,1))*180/pi;    % ������������������Ŷ�����˷����������ȡ
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
stroke=stroke_exp(t);
phi_exp=stroke(:,1);
dphi=stroke(:,2);
ddphi=stroke(:,3);
% % Ťת��:   x(1)=psi;  x(2)=dpsi;   
% % omega_x=x(2);                      % Ťת��:   x(:,1)=psi;  x(:,2)=dpsi;   
omega_y=dphi.*sin(x(:,1));            % Ťת��:   x(:,1)=psi;  x(:,2)=dpsi;  
omega_z=dphi.*cos(x(:,1));            %  Ťת��:   x(:,1)=psi;  x(:,2)=dpsi; 
% % omega_h=dphi;                             % �������ϵ�½���������  
% % domega_x=ddpsi;                                                         %  �Ǽ����ʡ����������������ļ���
% % domega_y=ddphi.*sin(x(:,1))+dphi.*x(:,2).*cos(x(:,1));         % �Ǽ��ٶȡ����������������ļ���
% % domega_z=ddphi.*cos(x(:,1))-dphi.*x(:,2).*sin(x(:,1));       % �Ǽ��ٶȡ����������������ļ���
% �������Ǽ���
v_y_nonr=omega_z;     % v_y=r*dphi*cos(psi)
v_z_nonr=-omega_y;    % v_z=-r*dphi*sin(psi)
% �����alpha������ȷ����ע�������ĵ�alpha=atan2(omega_z,-omega_y)*180/pi; ��ͬ
alpha=atan2(v_y_nonr,-v_z_nonr)*180/pi;   %����˩���
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A=[t,psi_sim,alpha,dpsi_sim];   % size(B)  % (12000*4)
xlswrite('PassiveRot_angle_Alpha_1.xlsx',A,'sheet1','A1:D12000');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (B) ʵ���Ĵ�����
B=[t,phi_exp,dphi,ddphi];
xlswrite('stroke_exp_angle_Alpha_1.xlsx',B,'sheet1','A1:D12000');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% �����˶�ѧ���
% figure(11) % ���Ƽ�����õĹ��ǡ����������������������ǳ���
% hold on
% plot(t/T,alpha,'r-','LineWidth',3)  
% xlabel('Normalized time','Fontsize',20,'FontName','Times','FontWeight','Bold');  
% ylabel('Angle of attack \it\alpha (deg)','Fontsize',20,'FontName','Times','FontWeight','Bold');
% legend('\alpha(t)')
% % title('����\alpha(t)�ı仯����')
% grid on
% box on
% axis([min(t/T),max(t/T),-inf,inf])
% % axis([min(t/T),max(t/T),-120,120])
% set(gca,'XTick',(min(t/T):0.5:max(t/T)),'LineStyle','-','LineWidth',1.5,'FontSize',16,'FontName','Times','FontWeight','Bold')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ����ת����\psi(t)��ʱ��ı仯
w =1185.6;           % ��Ƶ��    % f=188.7; T=1/f;     %����Ƶ�� (Hz)������ 
psi_exp=(4.2936+1.8536*cos(t*w)+59.6529*sin(t*w)+5.1852*cos(2*t*w)+1.6095*sin(2*t*w)-8.4569*cos(3*t*w)+9.8057*sin(3*t*w)+...
              3.9562*cos(4*t*w)+5.8064*sin(4*t*w)-3.0337*cos(5*t*w)-2.8749*sin(5*t*w)-2.8771*cos(6*t*w)+0.6686*sin(6*t*w)+...
               0.8649*cos(7*t*w)-0.6137*sin(7*t*w)+0.0771*cos(8*t*w)-1.0007*sin(8*t*w))*pi/180;
psi_sim=-x(:,1);
% psi_sim=1.5*x(:,1);
% psi_sim=2.0*x(:,1);
% psi_sim=2.515*x(:,1);
figure(12)  % ʵ����Ĵ�Ǻͼ�����õ�Ťת��
hold on
plot(t/T,phi_exp*180/pi,'k-',t/T,psi_exp*180/pi,'r-',t/T,psi_sim*180/pi,'b-.','LineWidth',3);      
xlabel('Normalized time','Fontsize',20,'FontName','Times','FontWeight','Bold');   
% ylabel('\it\phi_{exp} & \psi_{exp}(t) & \psi_{sim}(t)','Fontsize',20,'FontName','Times','FontWeight','Bold');    
ylabel('Flapping and pitch angle (deg.)','Fontsize',20,'FontName','Times','FontWeight','Bold');  
legend('\it\phi_{\rmexp}\rm(\itt\rm)','\it\psi_{\rmexp}\rm(\itt\rm)','\it\psi_{\rmsim}\rm(\itt\rm)');  
% title('ʵ����Ĵ��\phi_{exp}(t)�ͼ�����õ�Ťת��\psi(t)��ʱ��ı仯')
% grid on
box on
% axis([min(t/T),max(t/T),-inf,inf])
axis([min(t/T),max(t/T),-90,90]) 
% set(gca,'XTick',(min(t/T):0.5:max(t/T)),'LineStyle','-','LineWidth',1.5,'FontSize',16,'FontName','Times','FontWeight','Bold')
%%%%%%%%%%%%%%%%%%%%%%%%%
figure(13)  
plot(x(:,1),x(:,2),'LineWidth',2);
hold on;
xlabel('\it\psi \rm(\itt\rm) (rad)','Fontsize',20,'FontName','Times','FontWeight','Bold'); 
ylabel('d\it\psi \rm(\itt\rm) / d\itt \rm(rad/sec.)','Fontsize',20,'FontName','Times','FontWeight','Bold'); 
% legend('d\it\psi \rm(\itt\rm) / d\itt vs \it\psi \rm(\itt\rm)','Location','NorthWest'); 
% title('��ƽ��ͼ')
% grid on
box on
% set(gca,'XTick',(0.98:0.5:max(t0)/T+0.05),'LineStyle','-','LineWidth',1.5,'FontSize',20,'FontName','Times','FontWeight','Bold')
% set(gca,'LineStyle','-','LineWidth',1.5,'FontSize',16,'FontName','Times','FontWeight','Bold')
%%%%%%%%%%%%%%%%%%%%%%%%%
toc
% Elapsed time is 36.80381 seconds. for psi_max =70.4819; 
% Elapsed time is 1268.267093 seconds.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fun_psi=passive_rot(t,x,M_xaercoeff,M_xrdcoeff,k_am,k_inert,M_xRotcoeff,k_h,I_xx,I_xz,I_den,I_xzam)  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
stroke=stroke_exp(t);
% phi_exp=stroke(:,1);
dphi=stroke(:,2);
ddphi=stroke(:,3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Ťת��:   x(1)=psi;  x(2)=dpsi;   
% (1) �������ϵ�µĽ����ʺͽǼ����ʡ����������������������Գ�2DOF�˶�
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
omega_x=x(2);                              % Ťת��:   x(1)=psi;  x(2)=dpsi; 
omega_y=dphi.*sin(x(1));              % Ťת��:   x(1)=psi;  x(2)=dpsi; 
omega_z=dphi.*cos(x(1));             %  Ťת��:   x(1)=psi;  x(2)=dpsi; 
% domega_x=ddpsi;                                                  %  �Ǽ����ʡ����������������ļ���
% domega_y=ddphi.*sin(x(1))+dphi.*x(2).*cos(x(1)); % �Ǽ��ٶȡ����������������ļ���
domega_z=ddphi.*cos(x(1))-dphi.*x(2).*sin(x(1));      % �Ǽ��ٶȡ����������������ļ���
% (2)�������Ǽ��㡪�����ú���alpha_cal
v_y_nonr=omega_z;     % v_y=r*dphi*cos(psi)
v_z_nonr=-omega_y;    % v_z=-r*dphi*sin(psi)
alpha=alpha_cal(v_y_nonr,v_z_nonr);        % ע�������ǵ��ú���alpha_cal(v_y_nonr,v_z_nonr), ����atan2
V_nonr=sqrt(v_y_nonr.^2+v_z_nonr.^2); % ���������ٶ�V_nonr=omega_h=dphi;   
% omega_h=dphi;  % �������ϵ�½��������� % ��λ�� rad/s   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (3)  �湥�Ǳ仯��������ϵ��    
% alpha=45;     % �ٶ����Ǻ㶨����;���湫ʽ��������ϵ��Ϊ: C_L=1.8046;  C_D=1.7037;
C_L =sign(alpha).*(0.225+1.58*sin((2.13*abs(alpha)*180/pi-7.2)*pi/180));  
C_D =sign(alpha).*(1.92-1.55*cos((2.04*abs(alpha)*180/pi-9.82)*pi/180)); 
C_N=cos(abs(alpha)).*C_L+sin(abs(alpha)).*C_D; %��������ϵ���ϳɡ���2010-JFM-RJ Wood
% C_L =0.225+1.58*sin((2.13*abs(alpha)*180/pi-7.2)*pi/180);  
% C_D =sign(alpha).*(1.92-1.55*cos((2.04*abs(alpha)*180/pi-9.82)*pi/180)); 
% C_N=sqrt(C_L.^2+C_D.^2); %��������ϵ���ϳɡ���2010-JFM-RJ Wood
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ת����������; ת��������������;  ת�����������ء�����λ: (mN.mm)��(uN.m)�������ﵥλתΪN.m
%% (1) ƽ�������������������ء�����ת����������
% % M_xtrans1=(-sign(alpha).*Rou.*omega_h.^2.*C_N.*C_avereff^2*R_wingeff^3*F_nd.*Y_rcpnd/2)*(10^-6*10^-9); %ԭ�Ĺ�ʽN.m
% Y_rcpnd=COP_Ycpnd2_fruitfly(alpha);  % �����ת���������ء������ú�����⾻ѹ�ĵ�������λ��Y_rcpnd
% % M_xtrans1=-C_N.*omega_h.^2.*Y_rcpnd*M_xaercoeff; 
% % M_xtrans1=(-sign(dphi).*C_N.*omega_h.^2.*Y_rcpnd*M_xaercoeff);  % (rad/s)^2*mg.mm^2=mg.mm/s^2.mm=10^-6mN.mm
% M_xtrans=sign(alpha).*abs(C_N).*omega_h.^2.*Y_rcpnd*M_xaercoeff;   %   % ��ʼʱ��ʱ��
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N=length(t);
M_xtrans=zeros(N,1);
for i=1:1:N
    Y_rcpnd_trans=COP_Ycpnd2_fruitfly3(alpha);  % �����ת���������ء������ú�����⾻ѹ�ĵ�������λ��Y_rcpnd
    % Y_rcpnd_trans=COP_Ycpnd2_TransCirc(alpha(i,1));  % ���ú���COP_Ycpnd2_fruitfly��⾻ѹ�ĵ�������λ��Y_rcpnd; % ��������  % ��Ҫ����20���Ӽ���5������
    % Y_rcpnd_trans=COP_Ycpnd2_TransCirc3_1(alpha(i,1)); 
    % Y_rcpnd_trans=abs(Y_rcpnd_trans);  % ��  
    % M_xaercoeff=0.0060;   %��λ��: mg.mm^2 
    % ��λ: (rad/s)^2*mg.mm^2=mg.mm/s^2.mm=10^-9N.mm=10^-6mN.mm;  ƽ��������������ת����������һ��ʼ����ʱ���
    M_xtrans(i,1)=sign(alpha(i,1)).*abs(C_N(i,1)).*V_nonr(i,1).^2.*Y_rcpnd_trans*M_xaercoeff; 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%
% Y_rcpnd_trans=COP_Ycpnd2_TransCirc(alpha);  % ���ú���COP_Ycpnd2_fruitfly��⾻ѹ�ĵ�������λ��Y_rcpnd; % ��������
% Y_rcpnd_trans=abs(Y_rcpnd_trans);                % ��
% M_xtrans=sign(alpha).*abs(C_N).*V_nonr.^2.*Y_rcpnd_trans*M_xaercoeff;    %
%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (2) ת��Ťת�������ء���ת�������������� % ע���������ȡ����omega_x
% M_xrd=(-omega_x.*abs(omega_x)*Rou*C_RD*C_avereff^4*R_wingeff*Z_rnd/2)*10^(-15); %ԭ�Ĺ�ʽ N.m
M_xrd=-omega_x.*abs(omega_x)*M_xrdcoeff;   % (rad/s)^2*mg.mm^2=mg.mm/s^2.mm=10^-6mN.mm  % ��ʼʱ��ʱ��
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (3) ������Ťת��������
% M_xam=(-I_xyam*(domega_y-omega_x*omega_z)-I_xxam*domega_x)*10^(-12); %N.m-ԭʼ��ʽ����������ֵ����,��Ҫת����
% ��λ��: mg.mm^2*(rad/s)^2=mg.mm/s^2.mm=10^-6mN.mm; %ע��������������Ϊ����������domega_x=ddpsi;   
% ƽ�����������ء�������*��������=����*������*������ٶ�: 
% M_xam=k_am*(-I_xzam*ddphi.*cos(x(1)));  % ddphi_0=-2.2841e+006<0    % ��ʼʱ˳ʱ��
M_xam=k_am*(-I_xzam*(domega_z+omega_x.*omega_y));    % ƽ�����������ء�������*��������=����*������*������ٶ�
% M_xam=k_am*(I_xzam*(domega_z+omega_x.*omega_y));  % ƽ�����������ء��������ʵ��Ťת�ǲ����ӳ���
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (4) ת�������������������ء���Ťת��
N=length(t);
M_xrotcirc=zeros(N,1);
for i=1:1:N
    % ת��������Ťת���������ء������ú���COP_Ycpnd2_RotCirc��⾻ѹ�ĵ�������λ��Y_rcpnd_rot
    % ת�����������ġ�ѹ�ķֲ�����Dickinson���� or ѹ�������ҵ� or ѹ����c(r)/4����Ťת������
    Y_rcpnd_rot=COP_Ycpnd2_fruitfly3(alpha);  % �����ת���������ء������ú�����⾻ѹ�ĵ�������λ��Y_rcpnd
    % Y_rcpnd_rot=COP_Ycpnd2_RotCirc(alpha(i,1));  % ѹ�ķֲ�����Dickinson���� % ��Ҫ����20���Ӽ���5������
    % Y_rcpnd_rot=COP_Ycpnd2_RotCirc3_1(alpha(i,1));
    % Y_rcpnd_rot=abs(Y_rcpnd_rot);  % ��  
    % ��λ:(rad/s)^2*mg.mm^2=mg.mm/s^2.mm=10^-6mN.mm
    M_xrotcirc(i,1)=-omega_x(i,1).*V_nonr(i,1).*Y_rcpnd_rot*M_xRotcoeff;    % ת��������Ťת����������
end
%%%%%%%%%%%%%%%%%%%%%%%%%
% Y_rcpnd_rot=COP_Ycpnd2_RotCirc(alpha);        % ѹ�ķֲ�����Dickinson����
% Y_rcpnd_rot=abs(Y_rcpnd_rot);                         % ��
% M_xrotcirc=omega_x.*omega_h.*Y_rcpnd_rot*M_xRotcoeff;    % ת��������Ťת����������
%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (5) Ťת�����ظ����ء�������Ҳ���
M_hinge=-k_h*x(1); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (6) ��������������
m_wing=2.4*10^-3;                                 % mg
xr=0.3289;                                               %  \mm       % x-root offset 
x_com=xr+1.920243385;                        %  \mm       % ���ĵ�չ������
z_com=0.149785466;                              %  \mm       % ��Ťת����������
%%%%%%%%%%%%%%%%%%%%%%%%% ������������������
% F_inert_y=-m_wing*(domega_z*x_com-domega_x*z_com+omega_y.*(omega_x*x_com+omega_z*z_com))*10^3; % mN
% M_inert_x=-z_com*F_inert_y*10^3;                    % չ�򡪡���������������;  ��ʼ��ʱ��(-)
%%%%%%%%%%%%%%%%%%%%%%%%%��ƽ������ϵ��������������
% M_inert_x=-z_com*m_wing*(ddphi.*cos(psi)*x_com-ddpsi*z_com+dphi.^2.*(sin(2*psi)/2)*z_com); 
% M_inert_y=x_com*m_wing*(-ddphi.*sin(psi)*x_com-(dphi.*sin(psi)).^2*z_com-dpsi.^2*z_com); 
% M_inert_z=-x_com*m_wing*(ddphi.*cos(psi)*x_com-ddpsi*z_com+dphi.^2.*(sin(2*psi)/2)*z_com); 
M_inert_x=-k_inert*z_com*m_wing*(ddphi.*cos(x(1))*x_com+dphi.^2.*(sin(2*x(1))/2)*z_com);  % z_com^2*m_wing*ddpsi
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (7)  ������ϵ��: ����Ťת��������{��������Ťת��(x��)}: ��������Ҳ���
% M_x=M_xtrans+M_xrd+M_xam+M_xrotcirc;  % ��û�е��Իָ����ز�����Ŷ
% M_x=M_xtrans+M_xrd+M_xam+M_hinge;
M_x=M_xtrans+M_xrd+M_xam+M_xrotcirc+M_inert_x+M_hinge;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (rad/s)^2*mg.mm^2=mg.mm/s^2.mm=10^-6mN.mm 
% fun_psi=[x(2);  M_x/I_xx+I_xy*ddphi.*cos(x(1))/I_xx+1/2*dphi.^2.*sin(2*x(1))];% ԭʼ��ʽ����������ֵ����,��Ҫת����
% I_den=I_xx-k_am*I_xxam;     % ��λ��: mg.mm^2;  ע��������������Ϊ����������domega_x=ddpsi;  
% fun_psi=[x(2);(M_x+I_xz*ddphi.*cos(x(1))+1/2*I_xx*dphi.^2.*sin(2*x(1)))/I_den];      % 4itemTorq1_1
% fun_psi=[x(2); (M_x+I_xz*(ddphi.*cos(x(1))-2*dphi.*x(2).*sin(x(1)))-1/2*I_xx*dphi.^2.*sin(2*x(1)))/I_den]; % 4itemTorq3
fun_psi=[x(2);(M_x-I_xz*ddphi.*cos(x(1))+1/2*I_xx*dphi.^2.*sin(2*x(1)))/I_den];            % 4itemTorq2
end

function stroke=stroke_exp(t)
% % 2014-Science-ʵ���Ťת�ǡ��Ĵ�Ǻ��� (ע�⣺���ռ��������Ťת�ǣ���Ҫʱ���������ʵ��ʱ������׼)������Ҫ������˶�
% % syms t                     % t1  �Լ�����ʱ����
% % phi1=sym('3.9008+65.0445*cos(w*t)+4.2642*sin(w*t)+3.5806*cos(2*w*t)-2.9492*sin(2*w*t)+0.1319*cos(3*w*t)+0.3639*sin(3*w*t)+0.7844*cos(4*w*t)+0.2098*sin(4*w*t)');  %���ź���
% % dphi1=diff(phi1,t,1);
% % ddphi1=diff(phi1,t,2); 
% % dphi=inline(vectorize(dphi1),'w','t')            % ��ֵ����
% % ddphi=inline(vectorize(ddphi1),'w','t')  
w =1185.6;           %  ��Ƶ��     %  f=188.7; % Hz��������Ƶ��  % T=1/f;
phi=(3.9008+65.0445*cos(w*t)+4.2642*sin(w*t)+3.5806*cos(2*w*t)-2.9492*sin(2*w*t)...
         +0.1319*cos(3*w*t)+0.3639*sin(3*w*t)+0.7844*cos(4*w*t)+0.2098*sin(4*w*t))*pi/180;
dphi =(4.2642.*w.*cos(t.*w)-5.8984.*w.*cos(2*t.*w)+1.0917.*w.*cos(3*t.*w)+0.8392.*w.*cos(4*t.*w)...
           -65.0445.*w.*sin(t.*w)-7.1612.*w.*sin(2*t.*w)-0.3957.*w.*sin(3*t.*w)-3.1376.*w.*sin(4*t.*w))*pi/180;
ddphi =(11.7968.*w.^2.*sin(2*t.*w)-14.3224.*w.^2.*cos(2*t.*w)-1.1871.*w.^2.*cos(3*t.*w)-12.5504.*w.^2.*cos(4*t.*w)...
             -4.2642.*w.^2.*sin(t.*w)-65.0445.*w.^2.*cos(t.*w)-3.2751.*w.^2.*sin(3*t.*w)-3.3568.*w.^2.*sin(4*t.*w))*pi/180;         
stroke=[phi,dphi,ddphi];
end

function alpha_sim=alpha_cal(v_y_nonr,v_z_nonr)
% ����alpha=atan2(cot(psi))=atan2(cot(pi/2-alpha))=atan2(tan(alpha)); %����atan2������������ֵ��������alpha>pi/2ʱ
% alpha=pi/2+x(1).*sign(dphi);      % (1*400)���������������ι��ǡ���������   %���������ȷ����ļ��ι���
% alpha_sim=atan2(omega_z,-omega_y)+pi/2;   % �������: �������ι��ǡ����ᳫʹ��XXX���ι��ǡ�����ȡ����ֵ����ȷ
alpha_sim=atan2(v_y_nonr,v_z_nonr);   % �������: �������ι��ǡ����ᳫʹ��XXX���ι��ǡ�����ȡ����ֵ����ȷ
end
end



