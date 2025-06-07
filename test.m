
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