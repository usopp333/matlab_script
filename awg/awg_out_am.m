%�ĳ�����Ҫ�ǲ��ָ���.ppt�еĳ����� 
%�����еĻ����ĳ��ȵ�λ��Ϊ΢�ף�um��
clc
clear all
close all
pi=3.141592654;
c=299792458;
n_si=3.47638;  %index of si
%n_si=3.47638;  %index of si
n_sio2=1.4482;  %index of sio2
%n_sio2=1.4482;  %index of sio2
dlambda=0.005;  %�ŵ����  Ӱ�쵽���в����ĳ��Ȳʹ��ˮƽ���в����ĳߴ�仯�����Լ�󣬳ߴ�ԽС��
N_out=4;  %����˿���Ŀ
%(da��d��d_outӦ��ָ�������ڶ˿ڵ������࣬���������Ϊ���ڵ�taper֮������ļ��)
%�ٸ��������ͼ���������ж�taper�Ŀ�Ⱥ�����֮������ڷ�϶����������֮��Ĵ������Ϊ��׼
da=2.4;  %����Բ����˿ڿ��
d_in = da;
d=2;   %����Բ�������в����˿ڿ��
d_out=2.4; %����Բ����˿ڿ�� 

dis_in=0.01;  %���벨�����㾫�ȵĻ���
dis_out=0.01; %����������㾫�ȵĻ���
dis=0.01;     %���в������㾫�ȵĻ���
dis_t=0.01;   %���в���������㾫�ȵĻ���
r=15;  %���������뾶
lambda=1.550;  %���Ĳ���
%n_b_eff=2.533617;  %ֱ��������Ч������ cloudfdtd
n_b_eff=2.626644;   %cloudfdtd ������
n_c_eff=2.423903;  %������������Ч������
n_s_eff=2.932861;  %ƽ�����ɴ���������Ч������
ng_b=3.281933;  %ֱ������Ⱥ������
ng_s=3.701475;  %ƽ�����ɴ�������Ⱥ������
ng_c=3.3;       %����������Ⱥ�����ʣ��д����㣩

V=2*pi*da/lambda*sqrt(n_si^2-n_sio2^2); %��һ��Ƶ��   ���벨��
w0=da*(0.321+2.1*V^(-1.5)+4*V^(-6)); %���벨����˹����W0��������taper�������˹��ģ���뾶
V1=2*pi*d/lambda*sqrt(n_si^2-n_sio2^2); %��һ��Ƶ��  ���в���
w01=d*(0.321+2.1*V1^(-1.5)+4*V1^(-6)); % �������������˹����w0
lu = -3;
theta_max = sqrt(lu*lambda*sqrt(w01)/(8.7*n_si*pi));






fsr=1.5*N_out*dlambda;  %FSR��ѡȡ��ȡ2.5*��֤С�Ĵ��ţ�
%fsr=0.1;  %FSR��ѡȡ
m_r=lambda/fsr+1;  %���������伶��
m=round(m_r*n_b_eff/ng_b);  %���伶��
%L_f=d*d_out*n_s_eff/(m_r*dlambda);  %������Բ�İ뾶
L_f=da*d_in*n_s_eff/(m_r*dlambda);  %������Բ�İ뾶
%L_f = 60;
%����������в����ĳ��Ȳ�
dL=m*lambda/n_b_eff;   %��������dL��������������ֱ������neff֮��Ĳ�����
%�������������ĳ��Ȳ���Ը��ݼ��μ��������������m�йص���ֱ�����ĳ��ȣ��ر���l 2�ĳ��ȣ�

k0=2*pi/lambda; %����ϵ�� K0

%���ڸ�˹��������Բ�ķ�ɢ�Ƕȵļ���
%����һ��
% L_u=3; %�Ǿ�������ģ���������������������Ե�������������в������ĵ��Žǣ�
% we=sqrt(pi/2*w0);  %��˹Զ���ֲ��ĵ�Ч���
% theta0=lambda/n_s_eff/(we*sqrt(2*pi));
% theta_out=sqrt(L_u*theta0^2/8.7);  %����Ե�������������в������ĵ��Ž�
%theta_total=4*theta_out;

%�����������ø�˹�������ɴ���


wz=w0*sqrt(1+((lambda/n_si)*L_f/(pi*w0^2))^2);  %��Ĥ��˹������L_f���Ĺ�߰뾶��
R_z=L_f*(1+(pi*w0^2/((lambda/n_si)*L_f))^2);   %��˹�����ĵ���λ��������
%theta0=sqrt((wz/R_z)^2+((lambda/n_si)/(pi*wz))^2);  %����ĸ�˹��������Բ�뾶���ķ�ɢ��
theta0=atan(lambda/(pi*w0)); %���õ��Ǹ�˹��Ĺ�����ɢ�ǹ�ʽ
%theta0=2*asin(wz/L_f);
theta_d=2*asin(d/2/L_f);
theta_total=2*theta0;%+pi/6;  %����������pi/4��
N_a=ceil(theta_total/theta_d);  %N_a��ѡȡҪ����������Բ�������в������洦�ĸ�˹����ͼ�����жϣ�����������

% invq0 = -j*lambda/n_si/(pi*w0^2);
% q0=1/invq0;  %���������ĸ�����
% q=q0+L_f;
% invw2=-pi/(lambda/n_si)*imag(1./q);
% wz=sqrt(1./invw2);
% theta0=2*asin(wz/L_f);
% theta_d=2*asin(d/2/L_f);
% theta_total=theta0+pi/6;
% N_a=ceil(theta_total/theta_d);

%N_a=2*round(m_r)+1;   %���в����ĸ��������е�һ���㷨��
%N_a=23;
%theta_total=90*pi/180;  %�ܵĴ�����Բ��չ���Ƕ�


    
%��ֵ��ģ�����AWG


i_in=floor(4*w0/dis_in);%ģ��ȡ�����ѡȡ

for i=1:i_in;  %��һЩ������ѡ��ģ�⾫��Ϊ0.01um
    x1=-2*w0+(i-1)*dis_in;%����taper������Բ��ӵĶ˿ڸ�ȡ���������
    f_in(i)=exp(-x1^2/w0^2); %�����˹��
  
end
f_in=f_in/(sqrt(sum(abs(f_in).^2)));   %����ⳡ

figure;
title('����ⳡ');
plot(-2*w0:dis_in:-2*w0+(i_in-1)*dis_in,f_in);  %����ⳡ
%xlim([-40,40]);

%���ɴ�������1
for i=1:theta_total/(dis/L_f);  %���Կ��������в���������Բ�����ӳ��ļ�����Ŀ
    theta=-theta_total/2+(i-1)*(dis/L_f); %ÿ�����в���������Բ�����Ӵ���λ�ýǶȣ���ˮƽ�������������Ϊ0��
    %ÿ�����в������Ӧ�����ɴ�������Ĵ������
    for n=1:i_in 
        x1=-2*w0+(n-1)*dis_in;
        %����taper������Բ��ӵĶ˿ڸ�ȡ���㵽������Բ���в�������˿�ȡ����ľ���
        r=sqrt((L_f*sin(theta)-x1)^2+(L_f*cos(theta))^2);
        f_a_in_2(n)=(1/sqrt(lambda)/j)*f_in(n)*(exp(j*k0*n_s_eff*r)/sqrt(r))*(1+L_f*cos(theta)/r)/2*dis_in;
        
    end
   
    f_a_in(i)=sum((f_a_in_2));
    I_a_in(i)=sum(abs(f_a_in(i))^2); %��ǿ
    x1_fsrout(i)=[L_f*sin(theta)];
    
end
f_a_in_sum = sum(f_a_in);
figure;
plot(x1_fsrout,f_a_in);%���в���ǰ������ĳ�ǿ����
title('���в�����ǰ�˵ĳ��ֲ�');
% �����˿ڻ�Ĥ
array_in=cell(N_a,1);  %��ÿһ�����в���ǰ�˵�ģ���ֲ�
array_in_intensity=cell(N_a,1); %��ÿһ�����в���ǰ�˵Ĺ�ǿ�ֲ�

figure; 
for i=1:N_a
    x_center=-(N_a-1)/2*d+(i-1)*d;
    for n=1:theta_total/(dis/L_f)
        x1=-(theta_total/(dis/L_f))/2*dis+(n-1)*dis;
        array_in{i}(n)=exp(-(x1-x_center)^2/w01^2);
        array_in_intensity{i}(n)=abs(array_in{i}(n))^2;
    end
    array_in{i}=array_in{i}/sqrt(sum(array_in_intensity{i})); %ÿ�����в�����ģ���ֲ�
    array_in_intensity{i}=array_in_intensity{i}/sum(array_in_intensity{i}); %ÿ�����в����Ĺ�ǿ�ֲ�
    
    title('���в�����Ĥ');
    plot(-(theta_total/(dis/L_f))/2*dis:dis:-(theta_total/(dis/L_f))/2*dis+(n-1)*dis,array_in{i});hold on
end



%���в�����դ�ⳡ�ķ���
% C=1/sqrt(w01);  %Cֵ��С������δ֪
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
%���в�����դ�Ĺⳡ�������ߵĻ���
figure;
array_t=cell(N_a,1);
for i=1:N_a
    x_center=-(N_a-1)/2*d+(i-1)*d;
    for n=1:theta_total/(dis/L_f)
        x1=-(theta_total/(dis/L_f))/2*dis+(n-1)*dis;
    end
    array_t{i}=array_in{i}*(p_couplel(i)/In_gain); %ÿ�����в�����ģ���ֲ�
    %array_t{i}=array_in{i}*p_couplel(i); %ÿ�����в�����ģ���ֲ�
    title('���в����Ĺⳡ����');
    plot(-(theta_total/(dis/L_f))/2*dis:dis:-(theta_total/(dis/L_f))/2*dis+(n-1)*dis,array_t{i});hold on
 
% figure;
% plot(-(theta_total/(dis/L_f))/2*dis:dis:-(theta_total/(dis/L_f))/2*dis+(n-1)*dis,array_t{33});
end
% %���в���������˿ڵĹⳡ(���ڷ����ĵ����в���������˿ڽ����������ת)
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



%���в���������˿ڵĹⳡ������Ҷ�任��
% array_t_out=cell(N_a,1);
% array_t_out_intensity=cell(N_a,1);
% I_array_out=cell(N_a,1);

% ���ڸ�Ϊ�����в�����ÿ�������Ϊ���еĻ�Ĥ�����
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
title('���в���������ⳡ');
plot(-(theta_total/(dis/L_f))/2*dis_t:dis_t:-(theta_total/(dis/L_f))/2*dis+(n-1)*dis_t,abs(array_t_plot));hold on
figure;
plot(-(theta_total/(dis/L_f))/2*dis_t:dis_t:-(theta_total/(dis/L_f))/2*dis+(n-1)*dis_t,phase(array_t_plot));hold on
title('���в���������ⳡ��λ');
%����˿ڲ�����ģ
array_out=cell(N_out,1); %��ÿ�����в���ǰ�˵�ģ���ֲ�
array_out_intensity=cell(N_out,1); %��ǿ
V2=2*pi*d_out/lambda*sqrt(n_si^2-n_sio2^2); %��һ��Ƶ��
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
    title('����˿ڻ�Ĥ');
    plot(-(theta_total/(dis_out/L_f))/2*dis_out:dis_out:-(theta_total/(dis_out/L_f))/2*dis_out+(n-1)*dis_out,array_out{i});hold on
   %plot(f,array_out{i});hold on
end


%�������Բ�Ļ�۹ⳡ���������������ⳡ
% for x=1:theta_total/(dis_t/L_f)
% 
%         U(x)=1/sqrt(lambda*L_f/n_si)*fft(array_t_plot(x));
%         I_U(x)=(U(x))^2;
% end
figure;
U=1/sqrt(lambda*L_f/n_si)*(fftshift(abs(fft(array_t_plot))));  %����Ҷ�任(��Ҫ���������ת��)
%U=1/sqrt(lambda*L_f/n_si)*((abs(fft(array_t_plot))));  
% ff = -(theta_total/(dis_t/L_f))/2*dis_t:dis_t:-(theta_total/(dis_t/L_f))/2*dis_t+(n-1)*dis_t;
ff = 1:theta_total/(dis_t/L_f);
f = (ff - (max(ff) + min(ff))/2 ) / (length(ff) *dis_t);
f=f*lambda*L_f/n_si;
df=f(2)-f(1);
plot(f,U);hold on
title('����Բ����˵Ĺⳡ');



figure;
%����˿ڵĻ�Ĥ�ͺ���(������ʹ�ø����˿ڵ��໥����ϵ�������Ӱ��˿ڼ�����ж�)
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

%��������˿ڵĻ�Ĥ�ͺ����ɵ�ѡȡ��ʹ��ȡ��ļ���ͽ����˸���Ҷ�任���U��ȡ����һ�£���������
n1=0;
for n=1:round(df/dis_out):theta_total/(dis_out/L_f)
    n1=n1+1;
    array_out_p(n1)=array_out_plot(n);
end
x_f_n=n;
x_fsr=-(theta_total/(dis_out/L_f))/2*dis_out:round(df/dis_out)*dis_out:-(theta_total/(dis_out/L_f))/2*dis_out+(n-1)*dis_out;
plot(x_fsr,array_out_p);
title('����˿ڻ�Ĥ�ͺ���-�����˲���ȡ���');

%��������˿ڵĻ�Ĥ�ֱ���вɵ�ѡȡ��ʹ��ȡ��ļ���ͽ����˸���Ҷ�任���U��ȡ����һ�£���������
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
title('����˿ڵĻ�Ĥ-�����˲���ȡ���');

%��������Բ������ⳡ���������ѡȡ��ʹ�ú�����˵Ļ�Ĥ��������һ�£���һ������
%�������һ�����������Բ����ⳡ�Ĺ�һ�����������Ļ������ڲ���������˿���˵��������ܺ�Ϊ1
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
    f_out(n2)=array_out_U(n2)*array_out_p(n2); %ʵ�ʵ�����˵�����ⳡ������Բ������ⳡ�͸�������˿ڻ�Ĥ�ͺ������
end
f_out=f_out/sum(f_out);  %����Բ����ⳡ�Ĺ�һ������

figure;
plot(x_fsr,array_out_U);
title('����Բ����ⳡ-����������ѡȡ��');
figure;
plot(x_fsr,f_out);
title('���������ⳡ-�ͻ�Ĥ�ͺ������');



%����˿ڵĹⳡ����
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


%����Բ������ⳡ�͸�������˿ڻ�Ĥ����ϣ��Ǻͺ�����
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
        f_out1{i}(n2)=array_out_U(n2)*array_out1{i}(n2)*p_couplel_out(i)/In_gain_out; %ʵ�ʵ�����˵�����ⳡ������Բ������ⳡ�͸�������˿ڻ�Ĥ�����
    end
    
end
for i=1:N_out
    f_out1{i} = f_out1{i}/sum(f_out1{i});
    plot(x_fsr,f_out1{i});hold on
    title('����˿�����ⳡ-�͸�����Ĥ���');
end

% %ģ�����
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
    


