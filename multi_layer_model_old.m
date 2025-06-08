% Multilayer material model to calculate reverberation time and Qd
clear;close all;clc
% 设置频点和墙层数
N_layers = 3;
f=transpose(1e9:1e8:6e9); %——————设置频点/频带宽 为打印使为列矩阵
Length_frequency=length(f);
% 三层墙体结构
num_type_wall=0;
for type_wall = [1 2 3 4] % type_wall 依次计算承重墙，隔墙和地砖天花板，以及加了玻璃的窗户
    num_type_wall= num_type_wall+1;
%% 计算承重墙 厚度 220mm
% 承重墙结构为 石膏灰-混凝土-石膏灰
if type_wall ==1
        t=0;
for k=[15 19 15] % —— Gypsum plaster (5mm), Concrete with small gravel(220mm), Gypsum plaster (5mm)
ep=load('epsilon_paper_tabel.txt'); %20种材料的参数
e_INF=ep(k,1);e_s=ep(k,2);sigma_s=ep(k,3);r_t=ep(k,4);arf=ep(k,5); %对应材料的cole-cole模型的参数
w=2*pi*f;
miu0=(4*pi)*1e-7;
e0=1./(36*pi)*1e-9;
e_c2=e_INF+(e_s-e_INF)./(1+power(1i*w*r_t,1-arf))+sigma_s./(1i*w*e0); % Cole-Cole
RE=real(e_c2);
IMG=-imag(e_c2);
Tan_G=IMG./RE; 
miu_2r=f./f.*1;%相对磁导率/为写入文档保证矩阵长度一致为*（f./f的单位矩阵）
sigmma=f./f.*0;
t=t+1;
epsilon_r(:,t)=e_c2;
end
%% layers
atmosphere = ones(size(epsilon_r, 1),1);
epsilon_r=[atmosphere,epsilon_r,atmosphere];
for m=1:Length_frequency  
for n=1:N_layers+2
    omega(:,n)=2*pi*f; 
    epsilon(m,:)=e0.*epsilon_r(m,:);
    gama(m,n)=1i*omega(m,n).*sqrt(epsilon(m,n).*miu0);
    d(m,:)=[0,0.005,0.22,0.005,0];  % d 厚度，注意切换材料时这里需要更改 夹层是小砾石混凝土 No.19
end
end
for m=1:Length_frequency
for n=1:N_layers+2
    if n==1
        theta(m,n)=0/180*pi; 
    else
        theta(m,n)=asin(gama(m,n-1).*sin(theta(m,n-1))./gama(m,n));
    end
phi(m,n)=d(m,n).*gama(m,n).*cos(theta(m,n));
end
end

for m=1:Length_frequency   
for n=N_layers+2:-1:1
    if n==N_layers+2
    A(m,n)=1;
    C(m,n)=1;
    B(m,n)=0;
    D(m,n)=0;
    elseif n<N_layers+2 
Y(m,n+1)=cos(theta(m,n+1))./cos(theta(m,n)).*sqrt(epsilon(m,n+1)./epsilon(m,n));
Z(m,n+1)=cos(theta(m,n+1))./cos(theta(m,n)).*sqrt(epsilon(m,n)./epsilon(m,n+1));  
A(m,n)=exp(phi(m,n))./2.*(A(m,n+1).*(1+Y(m,n+1))+B(m,n+1).*(1-Y(m,n+1)));
B(m,n)=exp(-phi(m,n))./2.*(A(m,n+1).*(1-Y(m,n+1))+B(m,n+1).*(1+Y(m,n+1)));
C(m,n)=exp(phi(m,n))./2.*(C(m,n+1).*(1+Z(m,n+1))+D(m,n+1).*(1-Z(m,n+1)));
D(m,n)=exp(-phi(m,n))./2.*(C(m,n+1).*(1-Z(m,n+1))+D(m,n+1).*(1+Z(m,n+1)));
    end
end
end
 
RV=B(:,1)./A(:,1);
RH=D(:,1)./C(:,1);
TV=1./A(:,1);
TH=1./C(:,1);
 
Reflection_vertical(:,num_type_wall)=RV; 
Reflection_parallel(:,num_type_wall)=RH; 
Transmission_vertical(:,num_type_wall)=TV; 
Transmission_parallel(:,num_type_wall)=TH; 

for nn=1:1:Length_frequency 
    N_integral=1000;%积分N等分
    n_theta=(pi/2-0)/N_integral; %分成N等份 N+1次求和
    theta_i=0:n_theta:pi/2;%入射角用于积分求和
    SUM_Solid_Angle_mean_R(nn,num_type_wall)=sum((1-1/2*(abs(Reflection_vertical(nn,num_type_wall)).^2+abs(Reflection_parallel(nn,num_type_wall)).^2)).*cos(theta_i).*sin(theta_i)).*n_theta;
end

%% 计算隔墙
% 隔墙结构为 石膏灰-石膏板-石膏灰
elseif type_wall == 2
    t=0;
    epsilon_r = zeros(Length_frequency,3); % 清除上个材料的er值
    d=zeros(Length_frequency,5);
for k=[15 16 15] % —— Gypsum plaster (5mm), Plasterboard(18mm), Gypsum plaster (5mm)
ep=load('epsilon_paper_tabel.txt'); %20种材料的参数
e_INF=ep(k,1);e_s=ep(k,2);sigma_s=ep(k,3);r_t=ep(k,4);arf=ep(k,5); %对应材料的cole-cole模型的参数
w=2*pi*f;
miu0=(4*pi)*1e-7;
e0=1./(36*pi)*1e-9;
e_c2=e_INF+(e_s-e_INF)./(1+power(1i*w*r_t,1-arf))+sigma_s./(1i*w*e0); % Cole-Cole
RE=real(e_c2);
IMG=-imag(e_c2);
Tan_G=IMG./RE; 
miu_2r=f./f.*1;%相对磁导率/为写入文档保证矩阵长度一致为*（f./f的单位矩阵）
sigmma=f./f.*0;
t=t+1;
epsilon_r(:,t)=e_c2;
end
%% layers
atmosphere = ones(size(epsilon_r, 1),1);
epsilon_r=[atmosphere,epsilon_r,atmosphere];
for m=1:Length_frequency  
for n=1:N_layers+2
    omega(:,n)=2*pi*f; 
    epsilon(m,:)=e0.*epsilon_r(m,:);
    gama(m,n)=1i*omega(m,n).*sqrt(epsilon(m,n).*miu0);
    d(m,:)=[0,0.005,0.018,0.005,0];  % d 厚度，注意切换材料时这里需要更改 夹层是石膏板 No.16
end
end
for m=1:Length_frequency
for n=1:N_layers+2
    if n==1
        theta(m,n)=0/180*pi; 
    else
        theta(m,n)=asin(gama(m,n-1).*sin(theta(m,n-1))./gama(m,n));
    end
phi(m,n)=d(m,n).*gama(m,n).*cos(theta(m,n));
end
end

for m=1:Length_frequency   
for n=N_layers+2:-1:1
    if n==N_layers+2
    A(m,n)=1;
    C(m,n)=1;
    B(m,n)=0;
    D(m,n)=0;
    elseif n<N_layers+2 
Y(m,n+1)=cos(theta(m,n+1))./cos(theta(m,n)).*sqrt(epsilon(m,n+1)./epsilon(m,n));
Z(m,n+1)=cos(theta(m,n+1))./cos(theta(m,n)).*sqrt(epsilon(m,n)./epsilon(m,n+1));  
A(m,n)=exp(phi(m,n))./2.*(A(m,n+1).*(1+Y(m,n+1))+B(m,n+1).*(1-Y(m,n+1)));
B(m,n)=exp(-phi(m,n))./2.*(A(m,n+1).*(1-Y(m,n+1))+B(m,n+1).*(1+Y(m,n+1)));
C(m,n)=exp(phi(m,n))./2.*(C(m,n+1).*(1+Z(m,n+1))+D(m,n+1).*(1-Z(m,n+1)));
D(m,n)=exp(-phi(m,n))./2.*(C(m,n+1).*(1-Z(m,n+1))+D(m,n+1).*(1+Z(m,n+1)));
    end
end
end
 
RV=B(:,1)./A(:,1);
RH=D(:,1)./C(:,1);
TV=1./A(:,1);
TH=1./C(:,1);
 
Reflection_vertical(:,num_type_wall)=RV; 
Reflection_parallel(:,num_type_wall)=RH; 
Transmission_vertical(:,num_type_wall)=TV; 
Transmission_parallel(:,num_type_wall)=TH; 

for nn=1:1:Length_frequency 
    N_integral=1000;%积分N等分
    n_theta=(pi/2-0)/N_integral; %分成N等份 N+1次求和
    theta_i=0:n_theta:pi/2;%入射角用于积分求和
    SUM_Solid_Angle_mean_R(nn,num_type_wall)=sum((1-1/2*(abs(Reflection_vertical(nn,num_type_wall)).^2+abs(Reflection_parallel(nn,num_type_wall)).^2)).*cos(theta_i).*sin(theta_i)).*n_theta;
end

%% 计算天花板和地砖
% 天花板地砖结构为 石膏灰-混凝土-石膏灰
elseif type_wall == 3
    t=0;
    d=zeros(Length_frequency,5);
    epsilon_r = zeros(Length_frequency,3); % 清除上个材料的er值
for k=[15 19 15] % —— Gypsum plaster (5mm), Concrete with small gravel(100mm), Gypsum plaster (5mm)
ep=load('epsilon_paper_tabel.txt'); %20种材料的参数
e_INF=ep(k,1);e_s=ep(k,2);sigma_s=ep(k,3);r_t=ep(k,4);arf=ep(k,5); %对应材料的cole-cole模型的参数
w=2*pi*f;
miu0=(4*pi)*1e-7;
e0=1./(36*pi)*1e-9;
e_c2=e_INF+(e_s-e_INF)./(1+power(1i*w*r_t,1-arf))+sigma_s./(1i*w*e0); % Cole-Cole
RE=real(e_c2);
IMG=-imag(e_c2);
Tan_G=IMG./RE; 
miu_2r=f./f.*1;%相对磁导率/为写入文档保证矩阵长度一致为*（f./f的单位矩阵）
sigmma=f./f.*0;
%% 读取data矩阵数据计算等效的电长度
t=t+1;
epsilon_r(:,t)=e_c2;
end
%% layers
atmosphere = ones(size(epsilon_r, 1),1);
epsilon_r=[atmosphere,epsilon_r,atmosphere];
for m=1:Length_frequency  
for n=1:N_layers+2
    omega(:,n)=2*pi*f; 
    epsilon(m,:)=e0.*epsilon_r(m,:);
    gama(m,n)=1i*omega(m,n).*sqrt(epsilon(m,n).*miu0);
    d(m,:)=[0,0.005,0.100,0.005,0];  % d 厚度，注意切换材料时这里需要更改，小砾石混凝土 No.19
end
end
for m=1:Length_frequency
for n=1:N_layers+2
    if n==1
        theta(m,n)=0/180*pi; 
    else
        theta(m,n)=asin(gama(m,n-1).*sin(theta(m,n-1))./gama(m,n));
    end
phi(m,n)=d(m,n).*gama(m,n).*cos(theta(m,n));
end
end

for m=1:Length_frequency   
for n=N_layers+2:-1:1
    if n==N_layers+2
    A(m,n)=1;
    C(m,n)=1;
    B(m,n)=0;
    D(m,n)=0;
    elseif n<N_layers+2 
Y(m,n+1)=cos(theta(m,n+1))./cos(theta(m,n)).*sqrt(epsilon(m,n+1)./epsilon(m,n));
Z(m,n+1)=cos(theta(m,n+1))./cos(theta(m,n)).*sqrt(epsilon(m,n)./epsilon(m,n+1));  
A(m,n)=exp(phi(m,n))./2.*(A(m,n+1).*(1+Y(m,n+1))+B(m,n+1).*(1-Y(m,n+1)));
B(m,n)=exp(-phi(m,n))./2.*(A(m,n+1).*(1-Y(m,n+1))+B(m,n+1).*(1+Y(m,n+1)));
C(m,n)=exp(phi(m,n))./2.*(C(m,n+1).*(1+Z(m,n+1))+D(m,n+1).*(1-Z(m,n+1)));
D(m,n)=exp(-phi(m,n))./2.*(C(m,n+1).*(1-Z(m,n+1))+D(m,n+1).*(1+Z(m,n+1)));
    end
end
end
 
RV=B(:,1)./A(:,1);
RH=D(:,1)./C(:,1);
TV=1./A(:,1);
TH=1./C(:,1);
 
Reflection_vertical(:,num_type_wall)=RV; 
Reflection_parallel(:,num_type_wall)=RH; 
Transmission_vertical(:,num_type_wall)=TV; 
Transmission_parallel(:,num_type_wall)=TH; 


for nn=1:1:Length_frequency 
    N_integral=1000;%积分N等分
    n_theta=(pi/2-0)/N_integral; %分成N等份 N+1次求和
    theta_i=0:n_theta:pi/2;%入射角用于积分求和
    SUM_Solid_Angle_mean_R(nn,num_type_wall)=sum((1-1/2*(abs(Reflection_vertical(nn,num_type_wall)).^2+abs(Reflection_parallel(nn,num_type_wall)).^2)).*cos(theta_i).*sin(theta_i)).*n_theta;
end
%% 计算加装玻璃的窗户
% 玻璃的结构为了计算按照 glass-atmosphere-glass 计算
elseif type_wall == 4
    t=0;
    d=zeros(Length_frequency,5);
    epsilon_r = zeros(Length_frequency,3); % 清除上个材料的er值
for k=[13 21 13] % ——  Glass 6 - 12 单层玻璃，厚度取6mm, 在epsilon_paper_tabel.txt中增加了材料21作为空气夹层
ep=load('epsilon_paper_tabel.txt'); %20种材料的参数
e_INF=ep(k,1);e_s=ep(k,2);sigma_s=ep(k,3);r_t=ep(k,4);arf=ep(k,5); %对应材料的cole-cole模型的参数
w=2*pi*f;
miu0=(4*pi)*1e-7;
e0=1./(36*pi)*1e-9;
e_c2=e_INF+(e_s-e_INF)./(1+power(1i*w*r_t,1-arf))+sigma_s./(1i*w*e0); % Cole-Cole
RE=real(e_c2);
IMG=-imag(e_c2);
Tan_G=IMG./RE; 
miu_2r=f./f.*1;%相对磁导率/为写入文档保证矩阵长度一致为*（f./f的单位矩阵）
sigmma=f./f.*0;
%% 读取data矩阵数据计算等效的电长度
t=t+1;
epsilon_r(:,t)=e_c2;
end
%% layers
atmosphere = ones(size(epsilon_r, 1),1);
epsilon_r=[atmosphere,epsilon_r,atmosphere];
for m=1:Length_frequency  
for n=1:N_layers+2
    omega(:,n)=2*pi*f; 
    epsilon(m,:)=e0.*epsilon_r(m,:);
    gama(m,n)=1i*omega(m,n).*sqrt(epsilon(m,n).*miu0);
    d(m,:)=[0,0.004,0.012,0.004,0];  % d 厚度，注意切换材料时这里需要更改，玻璃 No.13
end
end
for m=1:Length_frequency
for n=1:N_layers+2
    if n==1
        theta(m,n)=0/180*pi; 
    else
        theta(m,n)=asin(gama(m,n-1).*sin(theta(m,n-1))./gama(m,n));
    end
phi(m,n)=d(m,n).*gama(m,n).*cos(theta(m,n));
end
end

for m=1:Length_frequency   
for n=N_layers+2:-1:1
    if n==N_layers+2
    A(m,n)=1;
    C(m,n)=1;
    B(m,n)=0;
    D(m,n)=0;
    elseif n<N_layers+2 
Y(m,n+1)=cos(theta(m,n+1))./cos(theta(m,n)).*sqrt(epsilon(m,n+1)./epsilon(m,n));
Z(m,n+1)=cos(theta(m,n+1))./cos(theta(m,n)).*sqrt(epsilon(m,n)./epsilon(m,n+1));  
A(m,n)=exp(phi(m,n))./2.*(A(m,n+1).*(1+Y(m,n+1))+B(m,n+1).*(1-Y(m,n+1)));
B(m,n)=exp(-phi(m,n))./2.*(A(m,n+1).*(1-Y(m,n+1))+B(m,n+1).*(1+Y(m,n+1)));
C(m,n)=exp(phi(m,n))./2.*(C(m,n+1).*(1+Z(m,n+1))+D(m,n+1).*(1-Z(m,n+1)));
D(m,n)=exp(-phi(m,n))./2.*(C(m,n+1).*(1-Z(m,n+1))+D(m,n+1).*(1+Z(m,n+1)));
    end
end
end
 
RV=B(:,1)./A(:,1);
RH=D(:,1)./C(:,1);
TV=1./A(:,1);
TH=1./C(:,1);
 
Reflection_vertical(:,num_type_wall)=RV; 
Reflection_parallel(:,num_type_wall)=RH; 
Transmission_vertical(:,num_type_wall)=TV; 
Transmission_parallel(:,num_type_wall)=TH; 


for nn=1:1:Length_frequency 
    N_integral=1000;%积分N等分
    n_theta=(pi/2-0)/N_integral; %分成N等份 N+1次求和
    theta_i=0:n_theta:pi/2;%入射角用于积分求和
    SUM_Solid_Angle_mean_R(nn,num_type_wall)=sum((1-1/2*(abs(Reflection_vertical(nn,num_type_wall)).^2+abs(Reflection_parallel(nn,num_type_wall)).^2)).*cos(theta_i).*sin(theta_i)).*n_theta;
end
%% 天花板和地砖按照单层计算，墙分承重墙和隔墙，分别按照[15 19 15] 和 [15 16 15] 计算
% RoomA 
d=3.28;
w=3.14;
h=2.782;

% % RoomB 
% d=6.8;
% w=3.5;
% h=2.78;

% % 2019房间长宽高(这一部分是文献验证，不与其余代码混合，统一注释)
% d= 9.1;
% w= 4.8;
% h= 4.1; 
% S_all=195;
% V = d*w*h; 
% S_floor_celling = 2*d*w;
% S_load_bearing_wall=0;
% S_partition_wall = S_all-S_floor_celling - S_load_bearing_wall;
% S_all_wall = S_partition_wall;
% CCS_opt_load_bearing_wall=(S_load_bearing_wall).*SUM_Solid_Angle_mean_R(:,1)./2;
% CCS_opt_S_partition_wall = S_partition_wall.*(SUM_Solid_Angle_mean_R(:,2))./2;
% CCS_opt_floor_celling = S_floor_celling.*(SUM_Solid_Angle_mean_R(:,3))./2;
% CCS_windows_door = 0; % 2019文献窗户和门的耦合损耗

%% Room A and B
V = d*w*h; 
S_floor_celling = 2*d*w;
S_all = 2*(d*h+d*w+h*w);   
%% 模型精细化                                     
% Room A 扣除飘窗屏蔽布的面积和门的屏蔽布面积
% S_cloth=2.5*1.75+0.87*2.12;                                                      % Room A cloth
S_cloth=0;                                                      % Room A cloth
% Room A 承重墙面积
S_load_bearing_wall = (1.23+1.8+1.27)*h;

% % Room B 分别扣除了阳台、过道门、小窗、厨房门和正门屏蔽布的面积
% S_cloth=2.35*2.3+1.08*2.46+0.85*1.36+1.57*2.1+1.54*2.37;  % Room B cloth
% % Room B 承重墙面积
% S_load_bearing_wall = (1.5+1.5+1.6+1)*h;

S_partition_wall = S_all-S_floor_celling-S_cloth-S_load_bearing_wall;

%% 计算隔墙的CCS 参考平面图，获取隔墙的面积 CCS=s*SUM*1/2，SUM 是立体角
CCS_opt_load_bearing_wall=(S_load_bearing_wall).*SUM_Solid_Angle_mean_R(:,1)./2;
CCS_opt_S_partition_wall = S_partition_wall.*(SUM_Solid_Angle_mean_R(:,2))./2;
CCS_opt_floor_celling = S_floor_celling.*(SUM_Solid_Angle_mean_R(:,3))./2;
CCS_windows = (2.5*1.75)*(SUM_Solid_Angle_mean_R(:,4))./2;
CCS_door = (2.12*0.87)/2; % Room A 窗户和门的耦合损耗
% 三部分CCS分开求和
% CCS_opt = CCS_opt_load_bearing_wall + CCS_opt_S_partition_wall + CCS_opt_floor_celling;
CCS_opt = CCS_opt_load_bearing_wall + CCS_opt_S_partition_wall + CCS_opt_floor_celling + CCS_windows + CCS_door;
%% Reverberation time and Qd
c=3e8; 
t_rev=1e9*V./CCS_opt./c;

% S_all_wall =S_all - S_cloth;
f_GHz=f./1e9;
% t_literature_2019=V*(-0.86.*f_GHz.^2+109.25.*f_GHz+29.49)./(2*pi.*f_GHz.*S_all_wall); % 单位 ns (1-40GHz)

% Qd
% Q_density_model = S_all_wall * (2*pi*f.*t_rev./1e9)./V;
Q_density_model = 195 * (2*pi*f.*t_rev./1e9)./160;
Q_density_literature_2019 = -0.86.*f_GHz.^2+109.25.*f_GHz+29.49;             % Ref[20](1-40GHz)
end
end