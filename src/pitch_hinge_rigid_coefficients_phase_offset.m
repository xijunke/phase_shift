function hinge_rigidcoef=pitch_hinge_rigid_coefficients_phase_offset()
%% frequency_ratio_to_pitch_hinge_rigid_coefficients
clear all;clc;
% lambda=linspace(1,3,9);
lambda=linspace(1.5,3.5,13);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 绘制扭转角相位偏置与频率比之间的关系
% 针对  lambda=linspace(1.5,3.5,13);
% pitchhinge_coef =[0.6964,0.8598,1.0403,1.2381,1.4530,1.6852,1.9345,2.2011,2.4848,2.7857,3.1038,3.4392,3.7917];
delta =[8.4591,7.0027,7.2051,4.5086,3.5464,2.5429,0,0,-2.9378,-3.1879,-2.9656,-2.9414,-5.9187];
figure(16)
plot(lambda,delta,'g:',lambda,delta,'dk','LineWidth',3)
xlabel('Frequency ratio \it\lambda','Fontsize',20,'FontName','Times','FontWeight','Bold')
ylabel('Phase offset \it\delta \rm(deg.)','Fontsize',20,'FontName','Times','FontWeight','Bold')
legend('\it\delta vs \lambda=f_{flap}/f_{rot,n}')
% title('\itThe relationship between the frequency ratio and the phase offset')
box on
axis([min(lambda)-0.05,max(lambda)+0.05,min(delta)-2,max(delta)+2])
set(gca,'XTick',(min(lambda):0.25:max(lambda)),'LineStyle','-','LineWidth',1.5,'FontSize',16,'FontName','Times','FontWeight','Bold')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (1) 频率比 %%%%%%%%%%%%%%%%%%%%%%%
f1=188.7;  % T=1/f1; 
% f_n=sqrt(k_h/I_den)/(2*pi)
% lambda=f_n/f1  
% k_h=(2*pi*lambda*f1).^2.*I_den;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (2) 翅膀转动惯量
% m_wing=2.4*10^-9;                                             % kg
% 请按照下面-2003-AMS-Sun Mao & 2007-JFM-Berman_Wang ZJ_Fruitfly_参考
% I_wing = 0.80 × 10^(-8) (g cm^2)  =  0.8*× 10^(-3) mg.mm2 换算下单位
% 注意：下面的数据是针对扭转轴为翅根和翅尖的连线时翅膀的质量属性
% I_xx=0.000267*10^(-12);      % inertia moment of wing    \ mg.mm^2——*10^(-12) kg.m^2  % UG 结果
% I_xz=-0.000594*10^(-12);     % inertia moment of wing    \ mg.mm^2——*10^(-12) kg.m^2  % UG 结果
% I_xx=0.1*0.000267;      % inertia moment of wing    \ mg.mm^2  % 早期
% I_xz=-0.1*0.000594;     % inertia moment of wing    \ mg.mm^2  % 早期
% I_xx=4.0*0.000267;      % inertia moment of wing    \ mg.mm^2  % 早期
% I_xz=-4.0*0.000594;     % inertia moment of wing    \ mg.mm^2  % 早期
I_xx=0.000267;      % inertia moment of wing    \ mg.mm^2
% I_xz=-0.000594;     % inertia moment of wing    \ mg.mm^2
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
% I_den=I_xx+k_am*I_xxam;     % 单位是: mg.mm^2;  注意这样处理，是因为虚质量项中domega_x=ddpsi;  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (3) 翅膀自身惯性力矩系数 & 重力矩系数
m_wing=2.4*10^-3;                                 % mg  namely: m_wing=2.4ug
% g=9.821*10^3;                                         % m/s^2=10^3mm/s^2
% xr=0.3289;                                            %  \mm       % x-root offset 
% x_com=xr+1.920243385;                     %  \mm       % 质心的展向坐标
% % z_com=-0.149785466+0.636=0.486215; % \mm 到扭转轴的弦向距离
z_com=0.149785466;                               %  \mm       % 到扭转轴的弦向距离
% d_com=z_com;                                         % mm
k_inert=0;
% k_inert=1;
%%%%%%%%%%%%%%%%%%%%%%%%%
wing_para=wing_shape_fruitfly_nine();  %调用函数wing_shape_fruitfly; 
% I_xzam=wing_para(1,5);                 % I_xzam = -0.001245    % 单位是 mg.mm^2
I_xxam=wing_para(1,6);                 % I_xxam = 0.0002508    % 单位是 mg.mm^2
k_am=1; % 虚质量扭转气动力矩系数
I_den=I_xx+k_am*I_xxam-k_inert*z_com^2*m_wing;  
%% (4) pitch_hinge_rigid_coefficients%%%%%%%%%%%%%%%%%%%%%%%
k_h=(2*pi*lambda*f1).^2.*I_den;
% k_h =
%   1.0e+003 *
%   Columns 1 through 6
%     1.6378    2.0219    2.4465    2.9116    3.4171    3.9630
%   Columns 7 through 12
%     4.5494    5.1761    5.8434    6.5511    7.2992    8.0877
%   Column 13
%     8.9167
%%%%%%%%%%%%%%%%%%%%%%%%%
% 扭转铰链回复力矩——刚度
L_h1=70e-006;     % L_h2=175e-006;   % length of wing hinge  \um——已经换算到m
W_h=1.8e-003;                                      % width of wing hinge  \mm——已经换算到m
t_h=7.6e-006;                                        % thickness of wing hinge \um——已经换算到m
E_h=2.5e009;                                         % lmodulus of elasticity for wing hinge  \Gpa 注意单位
%%%%%%%%%%%%%%%%%%%%%%%%%
k_hinge=k_h./(E_h*t_h^3*W_h/(12*L_h1)*10^9);
% 扭转角峰值与扭转铰链刚度之间的函数关系
% k_h=1*k_hinge*E_h*t_h^3*W_h/(12*L_h1)*10^9;
% psi_max=max(x(:,1)*180/pi)
psi_peak=[76.3123,72.6324,69.0907,65.6204,62.2658,59.0059,56.0679,53.4473,51.2758,49.4162,47.8372,46.3858,45.0164];
% length(psi_peak) % 13
% stdev_psi_sim_to_exp =[0.2233,0.1915,0.1619,0.1417,0.1330,0.1346,0.1473,0.1665,0.1823,0.2039,0.2336,0.2601,0.2808]; 
%%%%%%%%%%%%%%%%%%%%%%%%%
figure(17)  % Peak_value_psi_sim_to_ stiffness_coeff_pitch_hinge
plot(k_h,psi_peak,'r:',k_h,psi_peak,'dk','LineWidth',3)
% plot(k_h,psi_peak,'k-','LineWidth',3)
xlabel('Stiffness coefficients of pitch hinge  \itk_{pitch,hinge} \rm(uN.um/rad)','Fontsize',18,'FontName','Times','FontWeight','Bold')
ylabel('Peak value of simulated pitch angle \it\psi_{sim,peak} \rm(deg.)','Fontsize',18,'FontName','Times','FontWeight','Bold')
legend('\it\psi_{sim,peak} vs k_{pitch,hinge}')
% title('\itPeak value of simulated pitch angle vs Stiffness coefficients of pitch hinge')
box on
axis([min(k_h)-50,max(k_h)+50,min(psi_peak)-2.5,max(psi_peak)+2.5])
set(gca,'XTick',(min(k_h):1e+003:max(k_h)),'LineStyle','-','LineWidth',1.5,'FontSize',16,'FontName','Times','FontWeight','Bold')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hinge_rigidcoef=k_hinge;
end