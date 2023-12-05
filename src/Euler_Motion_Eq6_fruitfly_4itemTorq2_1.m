function Euler_Motion_Eq6_fruitfly_4itemTorq2_1
% 含4种机制的气动力矩: % M_x=M_xtrans+M_xrd+M_xam+M_hinge;——针对右侧翅膀
% 三个全正
% Euler_Motion_Eq——欧拉运动学方程的求解:  I_xx*ddpsi=M_x-I_xz*ddphi.*cos(x(1))+1/2*I_xx*dphi.^2.*sin(2*x(1));
% 计算可行——2014年6月18日,22:51:12
% 获得初始值——2014年6月20日,13:01:02——修改了时间轴
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 第一方案——不够合理哦
% Y_rcpnd=COP_Ycpnd2_fruitfly(alpha);  % 针对旋转轴气动力矩——调用函数求解净压心的无量纲位置Y_rcpnd
% Y_rcpnd_trans=COP_Ycpnd2_TransCirc(alpha);  % 调用函数COP_Ycpnd2_fruitfly求解净压心的无量纲位置Y_rcpnd; % 正负交替
% Y_rcpnd_rot=COP_Ycpnd2_RotCirc(alpha);        % 压心分布符合Dickinson函数
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 目标: 求解出psi; dpsi; ddpsi——向量形式——拍打角phi的一阶导dphi和二阶导ddphi; 被动扭转气动力矩M_x; 
% clear all; clc;
tic;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 调用含翅形貌参数化和气动力系数的函数
% wing_para_output=zeros(1,16);
% wing_para_output=[F_ndTrans,Coeff_liftdragF_N,M_xaercoeff,I1y,Z_rnd, M_xrdcoeff,...
%     F_ndRot,F_yrotcoeff,M_xRotcoeff, I2y,...
%     I_xzam,I_xxam, I5z,I6z,I7y,M_zrdcoeff];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
wing_para=wing_shape_fruitfly_nine();  %调用函数wing_shape_fruitfly;  % size(wing_para) 在这里——需要注意哦，与下面的组合使用
% R_wingeff=wing_para(1,1);             % R_wingeff=3.004;  % mm 
% C_avereff=wing_para(1,2);              % C_avereff=0.8854;  % mm  
% F_nd=wing_para(1,3);                     % F_nd=0.46392;  无量纲, 量纲化单位为mm^4
% Z_rnd=wing_para(1,4);                   % Z_rnd=0.14016;  无量纲, 量纲化单位为mm
I_xzam=wing_para(1,5);                 % I_xzam = -0.001245    % 单位是 mg.mm^2
I_xxam=wing_para(1,6);                 % I_xxam = 0.0002508    % 单位是 mg.mm^2
% I3=wing_para(1,7);                     % I3=0.74851        % 无量纲,量纲化单位是mm^4;
% I3z=wing_para(1,8);                   % I3z=0.0050945   % 单位是 mg.mm
% I4z=wing_para(1,9);                   % I4z=-0.0005346  % 单位是 mg.mm
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
wing_para=wing_shape_fruitfly_sixteen_1();    %调用函数wing_shape_fruitfly;  % size(wing_para) 
% wing_para=wing_shape_fruitfly_sixteen_2(); 
% wing_para=wing_shape_fruitfly_sixteen_good();
%%考虑扭转轴前移%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% R_wing=3.004;
% C_aver=0.8854;
% xr0=0.3289;
% C_maxyaxis=0.025;
% C_maxyaxis=0.1;
% C_maxyaxis=0.2024; %这里的扭转轴的位置采用了翅膀形貌学优化获得扭转轴的位置   
% C_maxyaxis=0.25; 
% C_maxyaxis=0.356737; 
% C_maxyaxis=0.36;
% wing_para=wing_shape_fruitfly_sixteen_good2(R_wing,C_aver,xr0,C_maxyaxis); % 注意调用的函数是wing_shape_fruitfly_sixteen_good2 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (1) 平动环量气动力和力矩参数
% F_ndTrans=wing_para(1,1);             % F_ndTrans=0.46392;  无量纲, 量纲化单位为mm^4
% Coeff_liftdragF_N=wing_para(1,2);  % Coeff_liftdragF_N=0.0068201;  %单位是mg*mm  %平动环量法向力
M_xaercoeff=wing_para(1,3);              % M_xaercoeff=0.0060;   %单位是: mg.mm^2   % 平动环量气动力矩参数——绕翅平面下的展向轴
% I1z=wing_para(1,4);                         % I1y=0.0162    % 单位是 mg.mm^2            % 平动环量气动力矩参数——绕翅平面下的弦向轴
% Z_rnd=wing_para(1,5)                     % Z_rnd=0.14016;  无量纲, 量纲化单位为mm
M_xrdcoeff=wing_para(1,6);               % M_xrdcoeff =1.5848e-004;  % 单位是mg.mm^2 %转动气动阻尼力矩参数—绕翅平面下的展向轴
% (2) 转动环量气动力和力矩参数
% F_ndRot=wing_para(1,7);                % F_ndRot=0.74851;  无量纲, 量纲化单位为mm^4
% F_yrotcoeff=wing_para(1,8);           % F_yrotcoeff =0.0032433;  % 单位是 mg.mm   % 转动环量法向力
M_xRotcoeff=wing_para(1,9);             % M_xRotcoeff=0.0029;   % 单位是 mg.mm^2   % 转动环量气动力矩系数——绕翅平面下的展向轴
% I2z=wing_para(1,10);                          % I2y=0.0069;        % 单位是 mg.mm^2          % 转动环量气动力矩参数——绕翅平面下的弦向轴
% (3) 虚质量气动力和力矩参数
%%%%%%%%%%%%%%%%%%%%%%%%
% I_xzam = -0.001245    % 单位是 mg.mm^2
% I_xxam = 0.0002508    % 单位是 mg.mm^2
% k_hinge=1.11*10^(-3);          % psi_max =70.4552 @ k_am=1;  % 计算的扭转角幅值≈实测扭转角幅值
%%%%%%%%%%%%%%%%%%%%%%%%
% I_xzam=wing_para(1,11);                    % I_xzam =0.002     % 单位是 mg.mm^2  % 虚质量气动力矩参数——绕翅平面下的展向轴
% I_xxam=wing_para(1,12);                    % I_xxam =0.00062892  % 单位是 mg.mm^2  % 虚质量气动力矩参数——绕翅平面下的展向轴
% I5y=wing_para(1,13);                      % I5z=0.0050945   % 单位是 mg.mm       % 虚质量气动力参数—法向力
% I6y=wing_para(1,14);                      % I6z=0.0011         % 单位是 mg.mm       % 虚质量气动力参数—法向力
% I7z=wing_para(1,15);                          % I7y=0.0109;        % 单位是 mg.mm^2   % 虚质量气动力矩参数——绕翅平面下的弦向轴
% M_zrdcoeff=wing_para(1,16);             % M_zrdcoeff=0.0012; % 单位是 mg.mm % 转动气动阻尼力矩参数—绕翅平面下的弦向轴
% C_max_LtoT=wing_para(1,17);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% I5y=wing_para(1,13);                       % I5z=0.0050945   % 单位是 mg.mm       % 下文中I3y应该改为I5y
% I6y=wing_para(1,14);                       % I6z=0.0011         % 单位是 mg.mm       % 下文中I4y应该改为I6y
% I7z=wing_para(1,15);                       % I7y=0.0109;        % 单位是 mg.mm^2   % 下文中I5z应该改为I7z
% I1z=wing_para(1,4);                         % I1y=0.0162        % 单位是 mg.mm^2    % 下文中I7z应该改为I1z
% I2z=wing_para(1,10);                       % I2y=0.0069;        % 单位是 mg.mm^2    % 下文中I6z应该改为I2z    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (1) 平动环量产生的气动力矩系数
k_xaero=1;   % 平动环量力矩不能没有
M_xaercoeff=k_xaero*M_xaercoeff;    % mg.mm^2;  % M_xaercoeff =6.0385e-003
%% (2) 转动环量产生的气动力矩系数
k_xRot=1;  
C_R=1.55;    % 该系数有问题-XXXXXXXX
M_xRotcoeff=k_xRot*C_R*M_xRotcoeff;
%% (3) 转动扭转气动力矩——转动气动阻尼力矩系数
%  C_RD=1; 
% C_RD=2; 
C_RD=5.0;     %转动阻尼气动力系数:C_RD=C_rd∈(3,6); 或者C_rd=C_dmax=3.4;
% C_RD=6.0;
% C_RD=8;  
M_xrdcoeff=C_RD*M_xrdcoeff;  % mg.mm^2;  % M_xrdcoeff =7.9240e-003;
%% (4) 虚质量扭转气动力矩系数
% k_am=0; % 不可行
% k_am=0.1; % 不可行
% k_am=0.2; % 不可行
% k_am=0.3;      % lambda =2.1588;
% k_am=0.35;    % lambda =2.0166;
% k_am=0.4;  % 可行
% k_am=0.6;  % lambda =1.5761;
k_am=1; % 不可行
%% (5) 扭转铰链回复力矩——刚度
L_h1=70e-006;     % L_h2=175e-006;   % length of wing hinge  \um——已经换算到m
W_h=1.8e-003;                                      % width of wing hinge  \mm——已经换算到m
t_h=7.6e-006;                                        % thickness of wing hinge \um——已经换算到m
E_h=2.5e009;                                         % lmodulus of elasticity for wing hinge  \Gpa 注意单位
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% k_hinge=0.865;    % psi_max =70.5348 @ k_am=1; % 计算的扭转角幅值≈实测扭转角幅值
% k_hinge=1.03;      % psi_max =72.1409 @ k_am=1;
% k_hinge=0.6*1.0;      % psi_max =72.1409 @ k_am=1; % K_hinge=1.411uN.mm; f_n =380.5327;lambda =2.0166;
% k_hinge=7.3*1.08;      % psi_max =72.1409 @ k_am=1;
k_hinge=0.97;        % psi_max =70.6746@ k_am=1;
% k_hinge=1.1;        % psi_max =70.6746@ k_am=1;
% k_hinge=1.11;          % psi_max =70.4552 @ k_am=1;  % 计算的扭转角幅值≈实测扭转角幅值
% k_hinge=1.12;      % psi_max =70.1635 @ k_am=1;
% k_hinge=5;           % psi_max =70.1635 @ k_am=1;
% k_hinge=1.75;  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% k_h=k_hinge*E_h*t_h^3*W_h/(12*L_h1);  % 这里k_h为rotational stiffness of the passive hinge
% k_h=0.1*k_hinge*E_h*t_h^3*W_h/(12*L_h1)*10^9  % N/m^2*m^3*m/m=N.m=kg*m/s^2.m=mg.mm/s^2.mm:[10^12]=10^9uN.mm
k_h=1*k_hinge*E_h*t_h^3*W_h/(12*L_h1)*10^9 
% k_h=2.05*k_hinge*E_h*t_h^3*W_h/(12*L_h1)*10^9 % k_h =2.8925e+003; f_n =337.2088; lambda = 1.7870;
% k_h=1.5*k_hinge*E_h*t_h^3*W_h/(12*L_h1)*10^9 
% k_h =2.2811e+003;  % 因为下面的k_xh的单位必须为uN.um/rad % f_n =334.0489; lambda =1.7703;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (6) 翅膀转动惯量
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
% I_den=I_xx+k_am*I_xxam;     % 单位是: mg.mm^2;  注意这样处理，是因为虚质量项中domega_x=ddpsi;  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (7) 翅膀自身惯性力矩系数 & 重力矩系数
m_wing=2.4*10^-3;                                 % mg  namely: m_wing=2.4ug
g=9.821*10^3;                                         % m/s^2=10^3mm/s^2
% xr=0.3289;                                            %  \mm       % x-root offset 
% x_com=xr+1.920243385;                     %  \mm       % 质心的展向坐标
% % z_com=-0.149785466+0.636=0.486215; % \mm 到扭转轴的弦向距离
z_com=0.149785466;                               %  \mm       % 到扭转轴的弦向距离
d_com=z_com;                                         % mm
k_inert=0;
% k_inert=1;
%%%%%%%%%%%%%%%%%%%%%%%%%
I_den=I_xx+k_am*I_xxam-k_inert*z_com^2*m_wing;  
%% 频率比 %%%%%%%%%%%%%%%%%%%%%%%
f1=188.7;  T=1/f1; 
f_n=sqrt(k_h/I_den)/(2*pi)
lambda=f_n/f1   % f_n =334.0489; lambda =1.7703;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ODE——常微风方程
% t=linspace(0.0052824335,0.0052824335+3*T,20000);  % t_steady1   %输出
psi_exp0 =-6.5669e-010;       % 实测值
% dpsi_exp0 =1.8709e+003;     % 实测值
dpsi_exp0 =-1.8709e+003;     % 实测值   -
initcond=[psi_exp0;dpsi_exp0];  % 初始值
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
%% (A) Euler_Motion_Eq6_fruitfly的输出数值计算被动扭转角和角速度
% 被动扭转角，角速度和攻角输出，已被用于显示对比结果
psi_sim=x(:,1)*180/pi;   % 被动扭转角
dpsi_sim=x(:,2);             % 被动扭转角速度
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 攻角计算——用于显示
% alpha=(pi/2*sign(x(:,1))-x(:,1))*180/pi;    % 这里的输出可能有问题哦——此方案不建议采取
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
stroke=stroke_exp(t);
phi_exp=stroke(:,1);
dphi=stroke(:,2);
ddphi=stroke(:,3);
% % 扭转角:   x(1)=psi;  x(2)=dpsi;   
% % omega_x=x(2);                      % 扭转角:   x(:,1)=psi;  x(:,2)=dpsi;   
omega_y=dphi.*sin(x(:,1));            % 扭转角:   x(:,1)=psi;  x(:,2)=dpsi;  
omega_z=dphi.*cos(x(:,1));            %  扭转角:   x(:,1)=psi;  x(:,2)=dpsi; 
% % omega_h=dphi;                             % 翅膀坐标系下铰链角速率  
% % domega_x=ddpsi;                                                         %  角加速率——用于虚质量力的计算
% % domega_y=ddphi.*sin(x(:,1))+dphi.*x(:,2).*cos(x(:,1));         % 角加速度——用于虚质量力的计算
% % domega_z=ddphi.*cos(x(:,1))-dphi.*x(:,2).*sin(x(:,1));       % 角加速度——用于虚质量力的计算
% 气动攻角计算
v_y_nonr=omega_z;     % v_y=r*dphi*cos(psi)
v_z_nonr=-omega_y;    % v_z=-r*dphi*sin(psi)
% 下面的alpha计算正确——注意与下文的alpha=atan2(omega_z,-omega_y)*180/pi; 不同
alpha=atan2(v_y_nonr,-v_z_nonr)*180/pi;   %添加了﹣号
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A=[t,psi_sim,alpha,dpsi_sim];   % size(B)  % (12000*4)
xlswrite('PassiveRot_angle_Alpha_1.xlsx',A,'sheet1','A1:D12000');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (B) 实测拍打角输出
B=[t,phi_exp,dphi,ddphi];
xlswrite('stroke_exp_angle_Alpha_1.xlsx',B,'sheet1','A1:D12000');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 绘制运动学输出
% figure(11) % 绘制计算而得的攻角——————————老是出错
% hold on
% plot(t/T,alpha,'r-','LineWidth',3)  
% xlabel('Normalized time','Fontsize',20,'FontName','Times','FontWeight','Bold');  
% ylabel('Angle of attack \it\alpha (deg)','Fontsize',20,'FontName','Times','FontWeight','Bold');
% legend('\alpha(t)')
% % title('攻角\alpha(t)的变化规律')
% grid on
% box on
% axis([min(t/T),max(t/T),-inf,inf])
% % axis([min(t/T),max(t/T),-120,120])
% set(gca,'XTick',(min(t/T):0.5:max(t/T)),'LineStyle','-','LineWidth',1.5,'FontSize',16,'FontName','Times','FontWeight','Bold')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 被动转动角\psi(t)随时间的变化
w =1185.6;           % 角频率    % f=188.7; T=1/f;     %翅拍频率 (Hz)和周期 
psi_exp=(4.2936+1.8536*cos(t*w)+59.6529*sin(t*w)+5.1852*cos(2*t*w)+1.6095*sin(2*t*w)-8.4569*cos(3*t*w)+9.8057*sin(3*t*w)+...
              3.9562*cos(4*t*w)+5.8064*sin(4*t*w)-3.0337*cos(5*t*w)-2.8749*sin(5*t*w)-2.8771*cos(6*t*w)+0.6686*sin(6*t*w)+...
               0.8649*cos(7*t*w)-0.6137*sin(7*t*w)+0.0771*cos(8*t*w)-1.0007*sin(8*t*w))*pi/180;
psi_sim=-x(:,1);
% psi_sim=1.5*x(:,1);
% psi_sim=2.0*x(:,1);
% psi_sim=2.515*x(:,1);
figure(12)  % 实测的拍打角和计算而得的扭转角
hold on
plot(t/T,phi_exp*180/pi,'k-',t/T,psi_exp*180/pi,'r-',t/T,psi_sim*180/pi,'b-.','LineWidth',3);      
xlabel('Normalized time','Fontsize',20,'FontName','Times','FontWeight','Bold');   
% ylabel('\it\phi_{exp} & \psi_{exp}(t) & \psi_{sim}(t)','Fontsize',20,'FontName','Times','FontWeight','Bold');    
ylabel('Flapping and pitch angle (deg.)','Fontsize',20,'FontName','Times','FontWeight','Bold');  
legend('\it\phi_{\rmexp}\rm(\itt\rm)','\it\psi_{\rmexp}\rm(\itt\rm)','\it\psi_{\rmsim}\rm(\itt\rm)');  
% title('实测的拍打角\phi_{exp}(t)和计算而得的扭转角\psi(t)随时间的变化')
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
% title('相平面图')
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
%%  扭转角:   x(1)=psi;  x(2)=dpsi;   
% (1) 翅膀坐标系下的角速率和角加速率——————这组数据来自翅2DOF运动
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
omega_x=x(2);                              % 扭转角:   x(1)=psi;  x(2)=dpsi; 
omega_y=dphi.*sin(x(1));              % 扭转角:   x(1)=psi;  x(2)=dpsi; 
omega_z=dphi.*cos(x(1));             %  扭转角:   x(1)=psi;  x(2)=dpsi; 
% domega_x=ddpsi;                                                  %  角加速率——用于虚质量力的计算
% domega_y=ddphi.*sin(x(1))+dphi.*x(2).*cos(x(1)); % 角加速度——用于虚质量力的计算
domega_z=ddphi.*cos(x(1))-dphi.*x(2).*sin(x(1));      % 角加速度——用于虚质量力的计算
% (2)气动攻角计算——调用函数alpha_cal
v_y_nonr=omega_z;     % v_y=r*dphi*cos(psi)
v_z_nonr=-omega_y;    % v_z=-r*dphi*sin(psi)
alpha=alpha_cal(v_y_nonr,v_z_nonr);        % 注意这里是调用函数alpha_cal(v_y_nonr,v_z_nonr), 不是atan2
V_nonr=sqrt(v_y_nonr.^2+v_z_nonr.^2); % 当地来流速度V_nonr=omega_h=dphi;   
% omega_h=dphi;  % 翅膀坐标系下铰链角速率 % 单位是 rad/s   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (3)  随攻角变化的升阻力系数    
% alpha=45;     % 假定攻角恒定不变;下面公式的升阻力系数为: C_L=1.8046;  C_D=1.7037;
C_L =sign(alpha).*(0.225+1.58*sin((2.13*abs(alpha)*180/pi-7.2)*pi/180));  
C_D =sign(alpha).*(1.92-1.55*cos((2.04*abs(alpha)*180/pi-9.82)*pi/180)); 
C_N=cos(abs(alpha)).*C_L+sin(abs(alpha)).*C_D; %由升阻力系数合成——2010-JFM-RJ Wood
% C_L =0.225+1.58*sin((2.13*abs(alpha)*180/pi-7.2)*pi/180);  
% C_D =sign(alpha).*(1.92-1.55*cos((2.04*abs(alpha)*180/pi-9.82)*pi/180)); 
% C_N=sqrt(C_L.^2+C_D.^2); %由升阻力系数合成——2010-JFM-RJ Wood
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 转动气动力矩; 转动气动阻尼力矩;  转动虚质量力矩——单位: (mN.mm)或(uN.m)——这里单位转为N.m
%% (1) 平动环量产生的气动力矩——旋转轴气动力矩
% % M_xtrans1=(-sign(alpha).*Rou.*omega_h.^2.*C_N.*C_avereff^2*R_wingeff^3*F_nd.*Y_rcpnd/2)*(10^-6*10^-9); %原文公式N.m
% Y_rcpnd=COP_Ycpnd2_fruitfly(alpha);  % 针对旋转轴气动力矩——调用函数求解净压心的无量纲位置Y_rcpnd
% % M_xtrans1=-C_N.*omega_h.^2.*Y_rcpnd*M_xaercoeff; 
% % M_xtrans1=(-sign(dphi).*C_N.*omega_h.^2.*Y_rcpnd*M_xaercoeff);  % (rad/s)^2*mg.mm^2=mg.mm/s^2.mm=10^-6mN.mm
% M_xtrans=sign(alpha).*abs(C_N).*omega_h.^2.*Y_rcpnd*M_xaercoeff;   %   % 初始时逆时针
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N=length(t);
M_xtrans=zeros(N,1);
for i=1:1:N
    Y_rcpnd_trans=COP_Ycpnd2_fruitfly3(alpha);  % 针对旋转轴气动力矩——调用函数求解净压心的无量纲位置Y_rcpnd
    % Y_rcpnd_trans=COP_Ycpnd2_TransCirc(alpha(i,1));  % 调用函数COP_Ycpnd2_fruitfly求解净压心的无量纲位置Y_rcpnd; % 正负交替  % 需要花掉20分钟计算5个周期
    % Y_rcpnd_trans=COP_Ycpnd2_TransCirc3_1(alpha(i,1)); 
    % Y_rcpnd_trans=abs(Y_rcpnd_trans);  % 正  
    % M_xaercoeff=0.0060;   %单位是: mg.mm^2 
    % 单位: (rad/s)^2*mg.mm^2=mg.mm/s^2.mm=10^-9N.mm=10^-6mN.mm;  平动环量产生的旋转轴气动力矩一开始是逆时针的
    M_xtrans(i,1)=sign(alpha(i,1)).*abs(C_N(i,1)).*V_nonr(i,1).^2.*Y_rcpnd_trans*M_xaercoeff; 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%
% Y_rcpnd_trans=COP_Ycpnd2_TransCirc(alpha);  % 调用函数COP_Ycpnd2_fruitfly求解净压心的无量纲位置Y_rcpnd; % 正负交替
% Y_rcpnd_trans=abs(Y_rcpnd_trans);                % 正
% M_xtrans=sign(alpha).*abs(C_N).*V_nonr.^2.*Y_rcpnd_trans*M_xaercoeff;    %
%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (2) 转动扭转气动力矩——转动气动阻尼力矩 % 注意这个方向取决于omega_x
% M_xrd=(-omega_x.*abs(omega_x)*Rou*C_RD*C_avereff^4*R_wingeff*Z_rnd/2)*10^(-15); %原文公式 N.m
M_xrd=-omega_x.*abs(omega_x)*M_xrdcoeff;   % (rad/s)^2*mg.mm^2=mg.mm/s^2.mm=10^-6mN.mm  % 初始时逆时针
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (3) 虚质量扭转气动力矩
% M_xam=(-I_xyam*(domega_y-omega_x*omega_z)-I_xxam*domega_x)*10^(-12); %N.m-原始公式不可用于数值计算,需要转换下
% 单位是: mg.mm^2*(rad/s)^2=mg.mm/s^2.mm=10^-6mN.mm; %注意这样处理，是因为虚质量项中domega_x=ddpsi;   
% 平动虚质量力矩——力臂*虚质量力=力臂*虚质量*法向加速度: 
% M_xam=k_am*(-I_xzam*ddphi.*cos(x(1)));  % ddphi_0=-2.2841e+006<0    % 初始时顺时针
M_xam=k_am*(-I_xzam*(domega_z+omega_x.*omega_y));    % 平动虚质量力矩——力臂*虚质量力=力臂*虚质量*法向加速度
% M_xam=k_am*(I_xzam*(domega_z+omega_x.*omega_y));  % 平动虚质量力矩——相对于实测扭转角产生延迟相
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (4) 转动环量产生的气动力矩——扭转轴
N=length(t);
M_xrotcirc=zeros(N,1);
for i=1:1:N
    % 转动环量绕扭转轴气动力矩——调用函数COP_Ycpnd2_RotCirc求解净压心的无量纲位置Y_rcpnd_rot
    % 转动环量产生的—压心分布符合Dickinson函数 or 压心在中弦点 or 压心在c(r)/4处—扭转轴力矩
    Y_rcpnd_rot=COP_Ycpnd2_fruitfly3(alpha);  % 针对旋转轴气动力矩——调用函数求解净压心的无量纲位置Y_rcpnd
    % Y_rcpnd_rot=COP_Ycpnd2_RotCirc(alpha(i,1));  % 压心分布符合Dickinson函数 % 需要花掉20分钟计算5个周期
    % Y_rcpnd_rot=COP_Ycpnd2_RotCirc3_1(alpha(i,1));
    % Y_rcpnd_rot=abs(Y_rcpnd_rot);  % 正  
    % 单位:(rad/s)^2*mg.mm^2=mg.mm/s^2.mm=10^-6mN.mm
    M_xrotcirc(i,1)=-omega_x(i,1).*V_nonr(i,1).*Y_rcpnd_rot*M_xRotcoeff;    % 转动环量绕扭转轴气动力矩
end
%%%%%%%%%%%%%%%%%%%%%%%%%
% Y_rcpnd_rot=COP_Ycpnd2_RotCirc(alpha);        % 压心分布符合Dickinson函数
% Y_rcpnd_rot=abs(Y_rcpnd_rot);                         % 正
% M_xrotcirc=omega_x.*omega_h.*Y_rcpnd_rot*M_xRotcoeff;    % 转动环量绕扭转轴气动力矩
%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (5) 扭转铰链回复力矩——针对右侧翅膀
M_hinge=-k_h*x(1); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (6) 翅膀自身惯性力矩
m_wing=2.4*10^-3;                                 % mg
xr=0.3289;                                               %  \mm       % x-root offset 
x_com=xr+1.920243385;                        %  \mm       % 质心的展向坐标
z_com=0.149785466;                              %  \mm       % 到扭转轴的弦向距离
%%%%%%%%%%%%%%%%%%%%%%%%% 翅膀自身惯性力—法向
% F_inert_y=-m_wing*(domega_z*x_com-domega_x*z_com+omega_y.*(omega_x*x_com+omega_z*z_com))*10^3; % mN
% M_inert_x=-z_com*F_inert_y*10^3;                    % 展向——翅膀自身惯性力矩;  初始逆时针(-)
%%%%%%%%%%%%%%%%%%%%%%%%%翅平面坐标系翅膀自身惯性力矩
% M_inert_x=-z_com*m_wing*(ddphi.*cos(psi)*x_com-ddpsi*z_com+dphi.^2.*(sin(2*psi)/2)*z_com); 
% M_inert_y=x_com*m_wing*(-ddphi.*sin(psi)*x_com-(dphi.*sin(psi)).^2*z_com-dpsi.^2*z_com); 
% M_inert_z=-x_com*m_wing*(ddphi.*cos(psi)*x_com-ddpsi*z_com+dphi.^2.*(sin(2*psi)/2)*z_com); 
M_inert_x=-k_inert*z_com*m_wing*(ddphi.*cos(x(1))*x_com+dphi.^2.*(sin(2*x(1))/2)*z_com);  % z_com^2*m_wing*ddpsi
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (7)  翅坐标系下: 被动扭转气动力矩{方向沿着扭转轴(x轴)}: ——针对右侧翅膀
% M_x=M_xtrans+M_xrd+M_xam+M_xrotcirc;  % ，没有弹性恢复力矩不可行哦
% M_x=M_xtrans+M_xrd+M_xam+M_hinge;
M_x=M_xtrans+M_xrd+M_xam+M_xrotcirc+M_inert_x+M_hinge;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (rad/s)^2*mg.mm^2=mg.mm/s^2.mm=10^-6mN.mm 
% fun_psi=[x(2);  M_x/I_xx+I_xy*ddphi.*cos(x(1))/I_xx+1/2*dphi.^2.*sin(2*x(1))];% 原始公式不可用于数值计算,需要转换下
% I_den=I_xx-k_am*I_xxam;     % 单位是: mg.mm^2;  注意这样处理，是因为虚质量项中domega_x=ddpsi;  
% fun_psi=[x(2);(M_x+I_xz*ddphi.*cos(x(1))+1/2*I_xx*dphi.^2.*sin(2*x(1)))/I_den];      % 4itemTorq1_1
% fun_psi=[x(2); (M_x+I_xz*(ddphi.*cos(x(1))-2*dphi.*x(2).*sin(x(1)))-1/2*I_xx*dphi.^2.*sin(2*x(1)))/I_den]; % 4itemTorq3
fun_psi=[x(2);(M_x-I_xz*ddphi.*cos(x(1))+1/2*I_xx*dphi.^2.*sin(2*x(1)))/I_den];            % 4itemTorq2
end

function stroke=stroke_exp(t)
% % 2014-Science-实测翅扭转角—拍打角函数 (注意：最终计算出来的扭转角，需要时间轴区间和实验时间轴配准)——需要求解和与核对
% % syms t                     % t1  自己设置时间轴
% % phi1=sym('3.9008+65.0445*cos(w*t)+4.2642*sin(w*t)+3.5806*cos(2*w*t)-2.9492*sin(2*w*t)+0.1319*cos(3*w*t)+0.3639*sin(3*w*t)+0.7844*cos(4*w*t)+0.2098*sin(4*w*t)');  %符号函数
% % dphi1=diff(phi1,t,1);
% % ddphi1=diff(phi1,t,2); 
% % dphi=inline(vectorize(dphi1),'w','t')            % 数值函数
% % ddphi=inline(vectorize(ddphi1),'w','t')  
w =1185.6;           %  角频率     %  f=188.7; % Hz——翅拍频率  % T=1/f;
phi=(3.9008+65.0445*cos(w*t)+4.2642*sin(w*t)+3.5806*cos(2*w*t)-2.9492*sin(2*w*t)...
         +0.1319*cos(3*w*t)+0.3639*sin(3*w*t)+0.7844*cos(4*w*t)+0.2098*sin(4*w*t))*pi/180;
dphi =(4.2642.*w.*cos(t.*w)-5.8984.*w.*cos(2*t.*w)+1.0917.*w.*cos(3*t.*w)+0.8392.*w.*cos(4*t.*w)...
           -65.0445.*w.*sin(t.*w)-7.1612.*w.*sin(2*t.*w)-0.3957.*w.*sin(3*t.*w)-3.1376.*w.*sin(4*t.*w))*pi/180;
ddphi =(11.7968.*w.^2.*sin(2*t.*w)-14.3224.*w.^2.*cos(2*t.*w)-1.1871.*w.^2.*cos(3*t.*w)-12.5504.*w.^2.*cos(4*t.*w)...
             -4.2642.*w.^2.*sin(t.*w)-65.0445.*w.^2.*cos(t.*w)-3.2751.*w.^2.*sin(3*t.*w)-3.3568.*w.^2.*sin(4*t.*w))*pi/180;         
stroke=[phi,dphi,ddphi];
end

function alpha_sim=alpha_cal(v_y_nonr,v_z_nonr)
% 由于alpha=atan2(cot(psi))=atan2(cot(pi/2-alpha))=atan2(tan(alpha)); %这里atan2给出象限正负值，尤其是alpha>pi/2时
% alpha=pi/2+x(1).*sign(dphi);      % (1*400)的行向量——几何攻角——弧度制   %输出——正确合理的几何攻角
% alpha_sim=atan2(omega_z,-omega_y)+pi/2;   % 这里输出: 正负几何攻角—不提倡使用XXX几何攻角—但是取绝对值是正确
alpha_sim=atan2(v_y_nonr,v_z_nonr);   % 这里输出: 正负几何攻角—不提倡使用XXX几何攻角—但是取绝对值是正确
end
end



