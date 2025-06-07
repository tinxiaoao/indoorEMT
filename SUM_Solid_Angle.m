function [SUM_Solid_Angle_mean_R,SUM_Solid_Angle_mean_T,SUM_Solid_Angle_mean_R_T,SUM_Solid_Angle_mean_r] = SUM_Solid_Angle(freq,tan_g,sigmma,miu_2r,epsilon_2r_Real,d_wall)
%UNTITLED3 此处显示有关此函数的摘要
%   此处显示详细说明
%% 利用反射系数求立体角得到品质因数 继而得到Pr/Pt
% freq;%频率Hz
% tan_g;%损耗角正切
% sigmma;%一般不同时赋值：tan_g
% miu_2r;%媒质2磁导率
% epsilon_2r_Real;%媒质2的e2实部 即相对介电常数
% d_wall;%墙体厚度
%%
miu_0=4*pi*1e-7;
epsilon_0=1/36/pi*1e-9;
lamuda=3e8./freq;%波长
omiga=2*pi*freq;
miu_1r=1;
epsilon_1r=1;
epsilon_2r_Imag=sigmma./omiga./epsilon_0+tan_g.*epsilon_2r_Real;%媒质2的e2虚部 即损耗角部分
epsilon_2r=epsilon_2r_Real-1i*epsilon_2r_Imag;
miu_1=miu_0.*miu_1r; 
miu_2=miu_0.*miu_2r;
epsilon_1=epsilon_0.*epsilon_1r;   

epsilon_2=epsilon_0.*epsilon_2r;
yita_1=sqrt(miu_1./epsilon_1);
yita_2=sqrt(miu_2./epsilon_2);
N=1000;%积分N等分
n_theta=(pi/2-0)/N; %分成N等份 N+1次求和
theta_i=0:n_theta:pi/2;%入射角用于积分求和
theta_t=asin(sin(theta_i).*sqrt(miu_1r.*epsilon_1r./miu_2r./epsilon_2r));%透射角
%%
k_1=omiga.*sqrt(miu_1.*epsilon_1);%介质1的波数 即真空中的波数
k_2=omiga.*sqrt(miu_2.*epsilon_2);%介质2的波数
%% a介质面的R||T.
Rp_vert_a=(yita_2.*cos(theta_i)-yita_1.*cos(theta_t))./(yita_2.*cos(theta_i)+yita_1.*cos(theta_t));%垂直极化波正方向入射a介质面上的反射系数
Rn_vert_a=(yita_1.*cos(theta_t)-yita_2.*cos(theta_i))./(yita_1.*cos(theta_t)+yita_2.*cos(theta_i));%垂直极化波负方向入射a介质面上的反射系数
Rp_para_a=(yita_1.*cos(theta_i)-yita_2.*cos(theta_t))./(yita_1.*cos(theta_i)+yita_2.*cos(theta_t));%平行极化波正方向入射a介质面上的反射系数
Rn_para_a=(yita_2.*cos(theta_t)-yita_1.*cos(theta_i))./(yita_2.*cos(theta_t)+yita_1.*cos(theta_i));%平行极化波负方向入射a介质面上的反射系数
Tp_vert_a=2*yita_2.*cos(theta_i)./(yita_2.*cos(theta_i)+yita_1.*cos(theta_t));%垂直极化波正方向入射a介质面上的透射系数
Tn_vert_a=2*yita_1.*cos(theta_t)./(yita_1.*cos(theta_t)+yita_2.*cos(theta_i));%垂直极化波负方向入射a介质面上的透射系数
Tp_para_a=2*yita_2.*cos(theta_i)./(yita_1.*cos(theta_i)+yita_2.*cos(theta_t));%平行极化波正方向入射a介质面上的透射系数
Tn_para_a=2*yita_1.*cos(theta_t)./(yita_2.*cos(theta_t)+yita_1.*cos(theta_i));%平行极化波负方向入射a介质面上的透射系数
%% b介质面的R||T.
Rp_vert_b=(yita_1.*cos(theta_t)-yita_2.*cos(theta_i))./(yita_1.*cos(theta_t)+yita_2.*cos(theta_i));%垂直极化波正方向入射b介质面上的反射系数
Rp_para_b=(yita_2.*cos(theta_t)-yita_1.*cos(theta_i))./(yita_2.*cos(theta_t)+yita_1.*cos(theta_i));%平行极化波正方向入射b介质面上的反射系数
Tp_vert_b=2*yita_1.*cos(theta_t)./(yita_1.*cos(theta_t)+yita_2.*cos(theta_i));%垂直极化波正方向入射b介质面上的透射系数
Tp_para_b=2*yita_1.*cos(theta_t)./(yita_2.*cos(theta_t)+yita_1.*cos(theta_i));%平行极化波正方向入射b介质面上的透射系数
%% 
theta_t=asin(k_1.*sin(theta_i)./real(k_2));
%% 总反射系数表达式
Reflection_vertical=Rp_vert_a+(Tp_vert_a.*Rp_vert_b.*Tn_vert_a.*exp(-2*1i*k_2.*d_wall./cos(theta_t)))./(1-Rn_vert_a.*Rp_vert_b.*exp(-2*1i*k_2.*d_wall./cos(theta_t)));
Reflection_parallel=Rp_para_a+(Tp_para_a.*Rp_para_b.*Tn_para_a.*exp(-2*1i*k_2.*d_wall./cos(theta_t)))./(1-Rn_para_a.*Rp_para_b.*exp(-2*1i*k_2.*d_wall./cos(theta_t)));
%% 总透射系数表达式
Transmission_vertical=Tp_vert_a.*Tp_vert_b.*exp(-1i.*k_2.*d_wall./cos(theta_t))./(1-Rp_vert_b.*Rn_vert_a.*exp(-2*1i*k_2.*d_wall./cos(theta_t)));
Transmission_parallel=Tp_para_a.*Tp_para_b.*exp(-1i.*k_2.*d_wall./cos(theta_t))./(1-Rp_para_b.*Rn_para_a.*exp(-2*1i*k_2.*d_wall./cos(theta_t)));
%% 
Reflection_parallel=-Reflection_parallel;%垂直极化波和平行极化波的参考方向不同
Transmission_parallel=-Transmission_parallel;%垂直极化波和平行极化波的参考方向不同
%% 立体角求和公式
SUM_Solid_Angle_mean_R=sum((1-1/2*(abs(Reflection_vertical).^2+abs(Reflection_parallel).^2)).*cos(theta_i).*sin(theta_i)).*n_theta;%只考虑反射
SUM_Solid_Angle_mean_T=sum((1/2*(abs(Transmission_vertical).^2+abs(Transmission_parallel).^2)).*cos(theta_i).*sin(theta_i)).*n_theta;%只考虑透射
SUM_Solid_Angle_mean_R_T=sum((1-1/2*(abs(Reflection_vertical).^2+abs(Reflection_parallel).^2+abs(Transmission_vertical).^2+abs(Transmission_parallel).^2)).*cos(theta_i).*sin(theta_i)).*n_theta;%考虑透射和反射
SUM_Solid_Angle_mean_r=sum((1-1/2*(abs(Rp_vert_a).^2+abs(Rp_para_a).^2)).*cos(theta_i).*sin(theta_i)).*n_theta;%只考虑反射Hill
end

 