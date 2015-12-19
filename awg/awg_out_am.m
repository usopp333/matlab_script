%改程序主要是波分复用.ppt中的程序复现 
%程序中的基本的长度单位都为微米（um）
clc
clear all
close all
pi=3.141592654;
c=299792458;
n_si=3.47638;  %index of si
%n_si=3.47638;  %index of si
n_sio2=1.4482;  %index of sio2
%n_sio2=1.4482;  %index of sio2
dlambda=0.005;  %信道间隔  影响到阵列波导的长度差，使得水平阵列波导的尺寸变化，间隔约大，尺寸越小；
N_out=4;  %输出端口数目
%(da、d和d_out应该指的是相邻端口的总体间距，即可以理解为相邻的taper之间的中心间距)
%再根据输出的图像数据来判断taper的宽度和他们之间的相邻缝隙，减少相邻之间的串扰耦合为标准
da=2.4;  %罗兰圆输入端口宽度
d_in = da;
d=2;   %罗兰圆输入阵列波导端口宽度
d_out=2.4; %罗兰圆输出端口宽度 

dis_in=0.01;  %输入波导计算精度的划分
dis_out=0.01; %输出波导计算精度的划分
dis=0.01;     %阵列波导计算精度的划分
dis_t=0.01;   %阵列波导输出计算精度的划分
r=15;  %弯曲波导半径
lambda=1.550;  %中心波长
%n_b_eff=2.533617;  %直波导的有效折射率 cloudfdtd
n_b_eff=2.626644;   %cloudfdtd 脊波导
n_c_eff=2.423903;  %弯曲波导的有效折射率
n_s_eff=2.932861;  %平板自由传输区的有效折射率
ng_b=3.281933;  %直波导的群折射率
ng_s=3.701475;  %平板自由传输区的群折射率
ng_c=3.3;       %弯曲波导的群折射率（有待计算）

V=2*pi*da/lambda*sqrt(n_si^2-n_sio2^2); %归一化频率   输入波导
w0=da*(0.321+2.1*V^(-1.5)+4*V^(-6)); %输入波导高斯波的W0，在输入taper的输入高斯的模场半径
V1=2*pi*d/lambda*sqrt(n_si^2-n_sio2^2); %归一化频率  阵列波导
w01=d*(0.321+2.1*V1^(-1.5)+4*V1^(-6)); % 阵列输出波导高斯波的w0
lu = -3;
theta_max = sqrt(lu*lambda*sqrt(w01)/(8.7*n_si*pi));






fsr=1.5*N_out*dlambda;  %FSR的选取（取2.5*保证小的串扰）
%fsr=0.1;  %FSR的选取
m_r=lambda/fsr+1;  %修正的衍射级数
m=round(m_r*n_b_eff/ng_b);  %衍射级数
%L_f=d*d_out*n_s_eff/(m_r*dlambda);  %大罗兰圆的半径
L_f=da*d_in*n_s_eff/(m_r*dlambda);  %大罗兰圆的半径
%L_f = 60;
%求得相邻阵列波导的长度差
dL=m*lambda/n_b_eff;   %这里计算的dL，将弯曲波导和直波导的neff之间的差距忽略
%由于弯曲波导的长度差可以根据几何计算出来，真正和m有关的是直波导的长度（特别是l 2的长度）

k0=2*pi/lambda; %传播系数 K0

%对于高斯波在罗兰圆的发散角度的计算
%方法一：
% L_u=3; %非均匀性损耗（论文中是用来计算最外缘输出波导相对阵列波导中心的张角）
% we=sqrt(pi/2*w0);  %高斯远场分布的等效宽度
% theta0=lambda/n_s_eff/(we*sqrt(2*pi));
% theta_out=sqrt(L_u*theta0^2/8.7);  %最外缘输出波导相对阵列波导中心的张角
%theta_total=4*theta_out;

%方法二：利用高斯波的自由传输


wz=w0*sqrt(1+((lambda/n_si)*L_f/(pi*w0^2))^2);  %基膜高斯光束在L_f处的光斑半径；
R_z=L_f*(1+(pi*w0^2/((lambda/n_si)*L_f))^2);   %高斯光束的等相位曲率中心
%theta0=sqrt((wz/R_z)^2+((lambda/n_si)/(pi*wz))^2);  %输入的高斯波在罗兰圆半径处的发散角
theta0=atan(lambda/(pi*w0)); %利用的是高斯光的光束发散角公式
%theta0=2*asin(wz/L_f);
theta_d=2*asin(d/2/L_f);
theta_total=2*theta0;%+pi/6;  %留有余量（pi/4）
N_a=ceil(theta_total/theta_d);  %N_a的选取要根据在罗兰圆连接阵列波导端面处的高斯波的图像来判断，留有余量；

% invq0 = -j*lambda/n_si/(pi*w0^2);
% q0=1/invq0;  %求束腰处的复参数
% q=q0+L_f;
% invw2=-pi/(lambda/n_si)*imag(1./q);
% wz=sqrt(1./invw2);
% theta0=2*asin(wz/L_f);
% theta_d=2*asin(d/2/L_f);
% theta_total=theta0+pi/6;
% N_a=ceil(theta_total/theta_d);

%N_a=2*round(m_r)+1;   %阵列波导的个数（其中的一种算法）
%N_a=23;
%theta_total=90*pi/180;  %总的大罗兰圆的展开角度


    
%数值法模拟仿真AWG


i_in=floor(4*w0/dis_in);%模拟取样点的选取

for i=1:i_in;  %留一些余量，选择模拟精度为0.01um
    x1=-2*w0+(i-1)*dis_in;%输入taper与罗兰圆相接的端口各取样点的坐标
    f_in(i)=exp(-x1^2/w0^2); %输入高斯波
  
end
f_in=f_in/(sqrt(sum(abs(f_in).^2)));   %输入光场

figure;
title('输入光场');
plot(-2*w0:dis_in:-2*w0+(i_in-1)*dis_in,f_in);  %输入光场
%xlim([-40,40]);

%自由传输区域1
for i=1:theta_total/(dis/L_f);  %可以看做是阵列波导与罗兰圆相连接出的计算数目
    theta=-theta_total/2+(i-1)*(dis/L_f); %每个阵列波导与罗兰圆相连接处的位置角度（以水平坐标轴的正方向为0）
    %每个阵列波导相对应的自由传输区域的传播情况
    for n=1:i_in 
        x1=-2*w0+(n-1)*dis_in;
        %输入taper与罗兰圆相接的端口各取样点到达罗兰圆阵列波导输出端口取样点的距离
        r=sqrt((L_f*sin(theta)-x1)^2+(L_f*cos(theta))^2);
        f_a_in_2(n)=(1/sqrt(lambda)/j)*f_in(n)*(exp(j*k0*n_s_eff*r)/sqrt(r))*(1+L_f*cos(theta)/r)/2*dis_in;
        
    end
   
    f_a_in(i)=sum((f_a_in_2));
    I_a_in(i)=sum(abs(f_a_in(i))^2); %光强
    x1_fsrout(i)=[L_f*sin(theta)];
    
end
f_a_in_sum = sum(f_a_in);
figure;
plot(x1_fsrout,f_a_in);%阵列波导前端曲面的场强曲线
title('阵列波导的前端的场分布');
% 各个端口基膜
array_in=cell(N_a,1);  %在每一个阵列波导前端的模场分布
array_in_intensity=cell(N_a,1); %在每一个阵列波导前端的光强分布

figure; 
for i=1:N_a
    x_center=-(N_a-1)/2*d+(i-1)*d;
    for n=1:theta_total/(dis/L_f)
        x1=-(theta_total/(dis/L_f))/2*dis+(n-1)*dis;
        array_in{i}(n)=exp(-(x1-x_center)^2/w01^2);
        array_in_intensity{i}(n)=abs(array_in{i}(n))^2;
    end
    array_in{i}=array_in{i}/sqrt(sum(array_in_intensity{i})); %每个阵列波导的模场分布
    array_in_intensity{i}=array_in_intensity{i}/sum(array_in_intensity{i}); %每个阵列波导的光强分布
    
    title('阵列波导基膜');
    plot(-(theta_total/(dis/L_f))/2*dis:dis:-(theta_total/(dis/L_f))/2*dis+(n-1)*dis,array_in{i});hold on
end



%阵列波导光栅光场的分配
% C=1/sqrt(w01);  %C值大小和意义未知
C = 0;
for i=1:N_a;
    p_couple(i)=0;
    for n=1:theta_total/(dis/L_f)
       p_couple(i)=p_couple(i)+f_a_in(n)*array_in{i}(n)*(1-C)^2;
       %p_couple(i)=(f_a_in(n)/f_a_in_sum);
    end
    p_couplel(i)=(abs(p_couple(i)))^2;
end
In_gain=sum(p_couplel);
%阵列波导光栅的光场输入曲线的绘制
figure;
array_t=cell(N_a,1);
for i=1:N_a
    x_center=-(N_a-1)/2*d+(i-1)*d;
    for n=1:theta_total/(dis/L_f)
        x1=-(theta_total/(dis/L_f))/2*dis+(n-1)*dis;
    end
    array_t{i}=array_in{i}*(p_couplel(i)/In_gain); %每个阵列波导的模场分布
    %array_t{i}=array_in{i}*p_couplel(i); %每个阵列波导的模场分布
    title('阵列波导的光场输入');
    plot(-(theta_total/(dis/L_f))/2*dis:dis:-(theta_total/(dis/L_f))/2*dis+(n-1)*dis,array_t{i});hold on
 
% figure;
% plot(-(theta_total/(dis/L_f))/2*dis:dis:-(theta_total/(dis/L_f))/2*dis+(n-1)*dis,array_t{33});
end
% %阵列波导的输出端口的光场(对于非中心的阵列波导的输出端口进行坐标的旋转)
% %figure;
% z0=pi*w0^2*n_si/lambda;
% for n=1:theta_total/(dis/L_f)
%     x_center=-(N_a-1)/2*d+(i-1)*d;
%     theta_r=-theta_total/2+(i-1)*(d/L_f);
%     for i=1:N_a
%         x_center=-(N_a-1)/2*d+(i-1)*d;
%         theta_r=-theta_total/2+(i-1)*(d/L_f);
%         x=-(theta_total/(dis/L_f))/2*dis+(n-1)*dis;
%         z=sqrt((L_f/2)^2-x^2)-L_f/2;
%         array_r_out{n}(i)=p_couplel(i)/(w01*(x*sin(theta_r)+z*cos(theta_r)+L_f))*exp(-j*2*n_si*pi/lambda*(x*cos(theta_r)-z*sin(theta_r))^2)/2*(1/(L_f*(x*sin(theta_r)+z*cos(theta_r)+L_f))-j*(lambda/(pi*w01^2*(x*sin(theta_r)+z*cos(theta_r)+L_f))))*exp(-j*(2*n_si*pi/lambda*(x*sin(theta_r)+z*cos(theta_r)+L_f)-atan(x*sin(theta_r)+z*cos(theta_r)+L_f)/z0));
%     end
%     U_out(n)=sum(array_r_out{n});
%     I_out(n)=abs(U_out(n))^2;
% end
% %plot(-(theta_total/(dis/L_f))/2*dis:dis:-(theta_total/(dis/L_f))/2*dis+(n-1)*dis,abs(I_out));



%阵列波导的输出端口的光场（傅里叶变换）
% array_t_out=cell(N_a,1);
% array_t_out_intensity=cell(N_a,1);
% I_array_out=cell(N_a,1);

% 现在改为在阵列波导的每个输出都为阵列的基膜的输出
figure;
for n=1:theta_total/(dis_t/L_f)
    for i=1:N_a
        array_t_out{n}(i)=(array_in{i}(n)*exp(j*(2*n_si*pi*i*dL)/lambda));
        %array_t_out{n}(i)=(array_t{i}(n));
        %array_t_out_intensity{n}(i)=abs(array_t_out{n}(i))^2;
    end
    array_t_out{n}=sum(array_t_out{n});
    %I_array_out{i}=sum(abs(array_t_out{i})); 
    array_t_plot(n)=array_t_out{n};
end
title('阵列波导的输出光场');
plot(-(theta_total/(dis/L_f))/2*dis_t:dis_t:-(theta_total/(dis/L_f))/2*dis+(n-1)*dis_t,abs(array_t_plot));hold on
figure;
plot(-(theta_total/(dis/L_f))/2*dis_t:dis_t:-(theta_total/(dis/L_f))/2*dis+(n-1)*dis_t,phase(array_t_plot));hold on
title('阵列波导的输出光场相位');
%输出端口波导单模
array_out=cell(N_out,1); %在每个阵列波导前端的模场分布
array_out_intensity=cell(N_out,1); %光强
V2=2*pi*d_out/lambda*sqrt(n_si^2-n_sio2^2); %归一化频率
w02=d_out*(0.321+2.1*V2^(-1.5)+4*V2^(-6));
figure;
for i=1:N_out
    x_center=-(N_out-1)/2*d_out+(i-1)*d_out;
    for n=1:theta_total/(dis_out/L_f)
        x1=-(theta_total/(dis_out/L_f))/2*dis_out+(n-1)*dis_out;
        array_out{i}(n)=exp(-(x1-x_center)^2/w02^2);
        array_out_intensity{i}(n)=abs(array_in{i}(n))^2;
    end
    array_out{i}=array_out{i}/sqrt(sum(array_out{i}.^2));
    %array_out{i}=array_out{i}/abs(sum(array_out{i}));
    array_out_intensity{i}=array_out_intensity{i}/sum(array_out_intensity{i});
    title('输出端口基膜');
    plot(-(theta_total/(dis_out/L_f))/2*dis_out:dis_out:-(theta_total/(dis_out/L_f))/2*dis_out+(n-1)*dis_out,array_out{i});hold on
   %plot(f,array_out{i});hold on
end


%输出罗兰圆的汇聚光场，输出波导的输入光场
% for x=1:theta_total/(dis_t/L_f)
% 
%         U(x)=1/sqrt(lambda*L_f/n_si)*fft(array_t_plot(x));
%         I_U(x)=(U(x))^2;
% end
figure;
U=1/sqrt(lambda*L_f/n_si)*(fftshift(abs(fft(array_t_plot))));  %傅里叶变换(需要进行坐标的转化)
%U=1/sqrt(lambda*L_f/n_si)*((abs(fft(array_t_plot))));  
% ff = -(theta_total/(dis_t/L_f))/2*dis_t:dis_t:-(theta_total/(dis_t/L_f))/2*dis_t+(n-1)*dis_t;
ff = 1:theta_total/(dis_t/L_f);
f = (ff - (max(ff) + min(ff))/2 ) / (length(ff) *dis_t);
f=f*lambda*L_f/n_si;
df=f(2)-f(1);
plot(f,U);hold on
title('罗兰圆输出端的光场');



figure;
%输出端口的基膜和函数(这样就使得各个端口的相互的耦合抵消，会影响端口间隔的判断)
for  n=1:theta_total/(dis_out/L_f)
    for i=1:N_out
         array_out_s{n}(i) = array_out{i}(n);
    end
    array_out_s{n}=sum(array_out_s{n});
    array_out_plot(n)=array_out_s{n};
end
%array_out_plot=fftshift((abs(fft(array_out_plot))));
%plot(-(theta_total/(dis_out/L_f))/2*dis_out:dis_out:-(theta_total/(dis_out/L_f))/2*dis_out+(n-1)*dis_out,array_out_plot);hold on    
%plot(f,array_out_plot);
%plot(f,array_out_s);

%对于输出端口的基膜和函数采点选取，使得取点的间隔和进行了傅里叶变换后的U的取点间隔一致（采样有误差）
n1=0;
for n=1:round(df/dis_out):theta_total/(dis_out/L_f)
    n1=n1+1;
    array_out_p(n1)=array_out_plot(n);
end
x_f_n=n;
x_fsr=-(theta_total/(dis_out/L_f))/2*dis_out:round(df/dis_out)*dis_out:-(theta_total/(dis_out/L_f))/2*dis_out+(n-1)*dis_out;
plot(x_fsr,array_out_p);
title('输出端口基膜和函数-进行了采样取点后');

%对于输出端口的基膜分别进行采点选取，使得取点的间隔和进行了傅里叶变换后的U的取点间隔一致（采样有误差）
figure;
array_out1=cell(N_out,1);

for i=1:N_out
    n1=0;
    for n=1:round(df/dis_out):theta_total/(dis_out/L_f)
        n1=n1+1;
        array_out1{i}(n1)=array_out{i}(n);
    end
    plot(x_fsr,array_out1{i});hold on
end
title('输出端口的基膜-进行了采样取点后');

%对于罗兰圆的输出光场进行坐标的选取；使得和输出端的基膜的坐标相一致（有一定的误差）
%在这里我还进行了罗兰圆输出光场的归一化处理，这样的话，对于波导的输出端口来说，输出的总和为1
x_c=round(length(f)/2);
np=length(array_out_p);
n2=0;
for n=(x_c-floor(n1/2)):(x_c+floor(n1/2));
    n2=n2+1;
    if(n2>np)
        array_out_p(n2)=0;
        x_fsr=-(theta_total/(dis_out/L_f))/2*dis_out:round(df/dis_out)*dis_out:-(theta_total/(dis_out/L_f))/2*dis_out+(x_f_n-1)*dis_out+round(df/dis_out)*dis_out;
    end
    array_out_U(n2)=U(n);
    f_out(n2)=array_out_U(n2)*array_out_p(n2); %实际的输出端的输出光场，罗兰圆的输出光场和各个输出端口基膜和函数耦合
end
f_out=f_out/sum(f_out);  %罗兰圆输出光场的归一化处理

figure;
plot(x_fsr,array_out_U);
title('罗兰圆输出光场-进行了坐标选取后');
figure;
plot(x_fsr,f_out);
title('输出端输出光场-和基膜和函数耦合');



%输出端口的光场分配
for i=1:N_out;
    n2=0;
    p_couple_out(i)=0;
    for n=(x_c-floor(n1/2)):(x_c+floor(n1/2));
         n2=n2+1;
         if(n2>np)
           array_out1{i}(n2)=0;
        end
       p_couple_out(i)=p_couple_out(i)+U(n)*array_out1{i}(n2)*(1-C)^2;
        
    end
    p_couplel_out(i)=(abs(p_couple_out(i)))^2;
end
In_gain_out=sum(p_couplel_out);


%罗兰圆的输出光场和各个输出端口基膜的耦合（非和函数）
figure;
for i=1:N_out
    n2=0;
    for n=(x_c-floor(n1/2)):(x_c+floor(n1/2));
        n2=n2+1;
        if(n2>np)
           array_out1{i}(n2)=0;
           x_fsr=-(theta_total/(dis_out/L_f))/2*dis_out:round(df/dis_out)*dis_out:-(theta_total/(dis_out/L_f))/2*dis_out+(x_f_n-1)*dis_out+round(df/dis_out)*dis_out;
        end
        array_out_U(n2)=U(n);
        f_out1{i}(n2)=array_out_U(n2)*array_out1{i}(n2)*p_couplel_out(i)/In_gain_out; %实际的输出端的输出光场，罗兰圆的输出光场和各个输出端口基膜的耦合
    end
    
end
for i=1:N_out
    f_out1{i} = f_out1{i}/sum(f_out1{i});
    plot(x_fsr,f_out1{i});hold on
    title('输出端口输出光场-和各个基膜耦合');
end

% %模场卷积
% 
% p=1;
% for i=1:N_out
%    % F{i}(p)=0;
%     for n=1:theta_total/(dis_out/L_f)/2;
%         F{i}(n)=F{i}(n)+array_out{i}(n)^2*abs(U(n))^2;
%         %F{i}(p)=F{i}(p)+array_out{i}(n)^2*abs(f_a_in(n))^2;
%     end
%     
% end
    


