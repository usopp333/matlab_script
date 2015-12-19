%日期：2015-11-12
%进行EDG的设计
%以下各单位暂定为um
%这个中科院的设计是基于罗兰圆的设计，有的结构的设计是基于平场，但是那是光纤结构可以省去输入输出波导，硅光结构不用考虑；
%关于EDG光栅面的设计基本上文献中有两种方法：一点法和两点法；
%第一步：将中科院的EDG的模型进行复现（不能完全仿造因为有些结构只能靠自己判断）

n_si=3.472000;  %index of si
n_sio2=1.444000;  %index of sio2
n_s_eff=2.847330;%使用了cloudfdtd3D的精确计算
n_b_eff = 2.564554; %使用了cloudfdtd3D的精确计算
ng_b=3.753849301034707;  % 使用了cloudfdtd3D的精确计算并进行了相应的计算
pi = 3.1415926;
%% 设计参数的设置：
clambda = 1.550;  %器件的中心波长
dlambda = 0.02;   %信道宽度
N_in  = 1;  %输入通道
N_out = 4; %输出通道
theta_i = pi/4; %输入通道与大罗兰圆轴线的夹角（作为EDG的入射角）
d_io = 2;  %输入输出波导的芯芯间隔（关于输出端口的位置设置还不太清楚,暂且设置为输入和输出相邻间隔都是相同的）
wg_width = 0.500;  %输入输出波导的宽度
taper_l = 20;  %输入输出taper的长度

%% 相关结构参数的计算：
L_f = clambda * d_io / (2*dlambda*sin(theta_i)); %大罗兰圆的半径
theta_k(1) = theta_i - 2*asin(d_io/2/(L_f/2*cos(theta_i))); %第一个输出波导与大罗兰圆轴线的夹角，距离输入波导最近的输出波导
for k=2:N_out
      theta_k (k)= theta_k (k-1)- 2*asin(d_io/2/(L_f/2*cos(theta_k(k-1)))); 
end
theta_ck = (theta_k(4) - theta_k(1))/2 + theta_k(1);  %输出波导的中心位置与大罗兰圆轴线的夹角（可作为整个EDG的衍射角by：lzq）
fsr = N_out * dlambda;  %结合AWG的设计，FSR的至少应当是信道宽度和输出通道数目的乘积
m = ceil(clambda / fsr) ;  %根据公式FSR=lambda/m
d = m*clambda/((sin(theta_i)+sin(theta_ck))*n_s_eff);
Vc=2*pi*d_io/clambda*sqrt(n_si^2-n_sio2^2);  %归一化频率 中心频率
w0c=d_io*(0.321+2.1*Vc^(-1.5)+4*Vc^(-6));  %高斯波的W0，在输入taper的输入高斯的模场半径中心频率
theta_0=atan(clambda/(pi*w0c));  %利用的是高斯光的光束发散角公式
theta_total = 2*theta_0;  %入射宽口EDG展开角度（使用的方法是和AWG的方式一致）
%中科院EDG文档中说明的是根据罗兰圆的半径和光栅周期可以估算出中心光栅左右两侧的光栅数目；但是基本上的EDG的光栅周期都不是常数，况且说的太笼统，无法估算
