%�ĳ�����Ҫ�ǲ��ָ���.ppt�еĳ�����
%[f������Բ����ⳡ����Ҷ�任������꣬
% x_fsr�� ����˿ڵĺ������꣬
% U: ����Բ����ⳡ����Ҷ�任��
% f_out�� ����ⳡ�����������Ĥ��ˣ�
% array_out_U�� ����Բ����ⳡ����Ҷ�任�󣨽��������ѹ����
% f_out1�� �����ŵ��Ĺⳡ����ĳ˻�]
%  = AWG_lambda(
% lambda: �����仯��
% N_out_order: ����ŵ�����)
function [f,x_fsr_in,x_fsr,array_t_plot,U,f_out,array_out_U,f_out1,fsr,L_f,m,dL,N_a] = AWG_lambda(lambda,N_out_order,dlambda,da,d,d_out,fsrpar)
clambda=1.550;  %���Ĳ���
pi=3.141592654;
c=299792458;
n_si=3.472000;  %index of si
%n_si=3.47638;  %index of si
n_sio2=1.444000;  %index of sio2
%n_sio2=1.4482;  %index of sio2

%dlambda=0.008;  %�ŵ���� Ӱ�쵽���в����ĳ��Ȳʹ��ˮƽ���в����ĳߴ�仯�����Լ�󣬳ߴ�ԽС��

N_out=4;  %����˿���Ŀ
%(da��d��d_outӦ��ָ�������ڶ˿ڵ������࣬���������Ϊ���ڵ�taper֮������ļ��)
%�ٸ��������ͼ���������ж�taper�Ŀ�Ⱥ�����֮������ڷ�϶����������֮��Ĵ������Ϊ��׼
%da=2.4;  %����Բ����˿ڿ��
%d=2;  %����Բ�������в����˿ڿ��
%d_out=2.4; %����Բ����˿ڿ�� 

dis_in=0.01;  %���벨�����㾫�ȵĻ���
dis_out=0.01; %����������㾫�ȵĻ���
dis=0.01;     %���в������㾫�ȵĻ���
dis_t=0.01; % ���в���������˿ڼ��㾫�ȵĻ���
r=15;  %���������뾶
% %n_b_eff=2.533617;  %ֱ��������Ч������   ��ʹ����cloudfdtd����ļ�������
% n_b_eff=2.626644;   %cloudfdtd ������
% n_s_eff=2.932861;  %ƽ�����ɴ���������Ч������ ��ʹ����cloudfdtd����ļ�������
% ng_b=3.281933;  %ֱ������Ⱥ������
n_s_eff=2.847330;%ʹ����cloudfdtd3D�ľ�ȷ����
n_b_eff = 2.564554; %ʹ����cloudfdtd3D�ľ�ȷ����
ng_b=3.753849301034707;  % ʹ����cloudfdtd3D�ľ�ȷ���㲢��������Ӧ�ļ���


%���ڻ����Ĳ����ļ���
fsr=fsrpar*N_out*dlambda;  %FSR��ѡȡ
m_r=clambda/fsr+1;  %���������伶��
m=round(m_r*n_b_eff/ng_b);  %���伶��
d_in = da;    %����Բ������˿ڵļ��
%L_f=d*d_in*n_s_eff/(m_r*dlambda);  %������Բ�İ뾶%����������da����d����̫ȷ�����Ƚ�׼ȷ����d*d_out��������L_f�����������ȵ�Ӱ�첻���ر��Ƕ���������Ŀ�Ⱥ�λ�ü������䣻
L_f=d_in*d_out*n_s_eff/(m_r*dlambda);  %������Բ�İ뾶
%L_f=d*d*n_s_eff/(m_r*dlambda);  %������Բ�İ뾶
%L_f=2*2.4*n_s_eff/(m_r*dlambda);  %������Բ�İ뾶
%L_f = 60;
%����������в����ĳ��Ȳ�
dL=m*clambda/n_b_eff;   %��������dL��������������ֱ������neff֮��Ĳ�����
%�������������ĳ��Ȳ���Ը��ݼ��μ��������������m�йص���ֱ�����ĳ��ȣ��ر���l 2�ĳ��ȣ�
k0=2*pi/lambda; %����ϵ�� K0
    
%��ֵ��ģ�����AWG
V=2*pi*da/lambda*sqrt(n_si^2-n_sio2^2); %��һ��Ƶ��
w0=da*(0.321+2.1*V^(-1.5)+4*V^(-6)); %��˹����W0��������taper�������˹��ģ���뾶

Vc=2*pi*da/clambda*sqrt(n_si^2-n_sio2^2); %��һ��Ƶ�� ����Ƶ��
w0c=da*(0.321+2.1*Vc^(-1.5)+4*Vc^(-6)); %��˹����W0��������taper�������˹��ģ���뾶   ����Ƶ��
i_in=floor(4*w0/dis_in);%ģ��ȡ�����ѡȡ

for i=1:i_in;  %��һЩ������ѡ��ģ�⾫��Ϊ0.01um
    x1=-2*w0+(i-1)*dis_in;%����taper������Բ��ӵĶ˿ڸ�ȡ���������
    f_in(i)=exp(-x1^2/w0^2); %�����˹��
  
end
f_in=f_in/(sqrt(sum(abs(f_in).^2)));   %����ⳡ


%�����������ø�˹�������ɴ���
wz=w0c*sqrt(1+((clambda/n_si)*L_f/(pi*w0^2))^2);  %��Ĥ��˹������L_f���Ĺ�߰뾶��
%R_z=L_f*(1+(pi*w0^2/((clambda/n_si)*L_f))^2);   %��˹�����ĵ���λ��������
%theta0=sqrt((wz/R_z)^2+((lambda/n_si)/(pi*wz))^2);  %����ĸ�˹��������Բ�뾶���ķ�ɢ��
theta0=atan(clambda/(pi*w0c)); %���õ��Ǹ�˹��Ĺ�����ɢ�ǹ�ʽ
%theta0=2*asin(wz/L_f);
theta_d=2*asin(d/2/L_f);
theta_total=2*theta0;%+pi/6;  %��������  Ӱ�����в�������Ŀ
N_a=ceil(theta_total/theta_d);  %N_a��ѡȡҪ����������Բ�������в������洦�ĸ�˹����ͼ�����жϣ�����������



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

array_in=cell(N_a,1);  %��ÿһ�����в���ǰ�˵�ģ���ֲ�
array_in_intensity=cell(N_a,1); %��ÿһ�����в���ǰ�˵Ĺ�ǿ�ֲ�
V1=2*pi*d/lambda*sqrt(n_si^2-n_sio2^2); %��һ��Ƶ��
w01=d*(0.321+2.1*V1^(-1.5)+4*V1^(-6)); 
for i=1:N_a
    x_center=-(N_a-1)/2*d+(i-1)*d;
    for n=1:theta_total/(dis/L_f)
        x1=-(theta_total/(dis/L_f))/2*dis+(n-1)*dis;
        array_in{i}(n)=exp(-(x1-x_center)^2/w01^2);
        array_in_intensity{i}(n)=abs(array_in{i}(n))^2;
    end
    array_in{i}=array_in{i}/sqrt(sum(array_in_intensity{i})); %ÿ�����в�����ģ���ֲ�
    array_in_intensity{i}=array_in_intensity{i}/sum(array_in_intensity{i}); %ÿ�����в����Ĺ�ǿ�ֲ�
    
end

%���в�����դ�������
C=1/sqrt(w01);  %Cֵ��С������δ֪
for i=1:N_a;
    p_couple(i)=0;
    for n=1:theta_total/(dis/L_f)
       p_couple(i)=p_couple(i)+f_a_in(n)*array_in{i}(n)*(1-C)^2;
        
    end
    p_couplel(i)=(abs(p_couple(i)))^2;
end
In_gain=sum(p_couplel);
%���в�����դ�Ĺⳡ�������ߵĻ���

array_t=cell(N_a,1);
for i=1:N_a
    x_center=-(N_a-1)/2*d+(i-1)*d;
    for n=1:theta_total/(dis/L_f)
        x1=-(theta_total/(dis/L_f))/2*dis+(n-1)*dis;
    end
    array_t{i}=array_in{i}*p_couplel(i)/In_gain; %ÿ�����в�����ģ���ֲ�
end
% figure;
% plot(-(theta_total/(dis/L_f))/2*dis:dis:-(theta_total/(dis/L_f))/2*dis+(n-1)*dis,array_t{33});

%���в���������˿ڵĹⳡ(���ڷ����ĵ����в���������˿ڽ����������ת)
z0=pi*w0^2*n_si/lambda;
for n=1:theta_total/(dis/L_f)
    x_center=-(N_a-1)/2*d+(i-1)*d;
    theta_r=-theta_total/2+(i-1)*(d/L_f);
    for i=1:N_a
        x_center=-(N_a-1)/2*d+(i-1)*d;
        theta_r=-theta_total/2+(i-1)*(d/L_f);
        x=-(theta_total/(dis/L_f))/2*dis+(n-1)*dis;
        z=sqrt((L_f/2)^2-x^2)-L_f/2;
        array_r_out{n}(i)=p_couplel(i)/(w01*(x*sin(theta_r)+z*cos(theta_r)+L_f))*exp(-j*2*n_si*pi/lambda*(x*cos(theta_r)-z*sin(theta_r))^2)/2*(1/(L_f*(x*sin(theta_r)+z*cos(theta_r)+L_f))-j*(lambda/(pi*w01^2*(x*sin(theta_r)+z*cos(theta_r)+L_f))))*exp(-j*(2*n_si*pi/lambda*(x*sin(theta_r)+z*cos(theta_r)+L_f)-atan(x*sin(theta_r)+z*cos(theta_r)+L_f)/z0));
    end
    U_out(n)=sum(array_r_out{n});
    I_out(n)=abs(U_out(n))^2;
end

%���в���������˿ڵĹⳡ

for n=1:theta_total/(dis_t/L_f)
    for i=1:N_a
        array_t_out{n}(i)=(array_t{i}(n)*exp(j*(2*n_si*pi*i*dL)/lambda));
        
        %array_t_out{n}(i)=(array_in{i}(n)*exp(j*(2*n_si*pi*i*dL)/lambda));
        %array_t_out{n}(i)=(array_t{i}(n));
        %array_t_out_intensity{n}(i)=abs(array_t_out{n}(i))^2;
    end
    array_t_out{n}=sum(array_t_out{n});
    %I_array_out{i}=sum(abs(array_t_out{i})); 
    array_t_plot(n)=array_t_out{n};
end
x_fsr_in = -(theta_total/(dis/L_f))/2*dis_t:dis_t:-(theta_total/(dis/L_f))/2*dis+(n-1)*dis_t;
%����˿ڲ�����ģ
array_out=cell(N_out,1); %��ÿ�����в���ǰ�˵�ģ���ֲ�
array_out_intensity=cell(N_out,1); %��ǿ
V2=2*pi*d_out/lambda*sqrt(n_si^2-n_sio2^2); %��һ��Ƶ��
w02=d_out*(0.321+2.1*V2^(-1.5)+4*V2^(-6));
for i=1:N_out
    x_center=-(N_out-1)/2*d_out+(i-1)*d_out;
    for n=1:theta_total/(dis_out/L_f)
        x1=-(theta_total/(dis_out/L_f))/2*dis_out+(n-1)*dis_out;
        array_out{i}(n)=exp(-(x1-x_center)^2/w02^2);
        array_out_intensity{i}(n)=abs(array_in{i}(n))^2;
        
    end
    array_out{i}=array_out{i}/sqrt(sum(array_out{i}.^2));
    array_out_intensity{i}=array_out_intensity{i}/sum(array_out_intensity{i});
  
end


%�������Բ�Ļ�۹ⳡ���������������ⳡ  (�뾭������Ҷ�任)
U=1/sqrt(lambda*L_f/n_si)*(fftshift(abs(fft(array_t_plot))));  %����Ҷ�任(��Ҫ���������ת��)
%U=1/sqrt(lambda*L_f/n_si)*((abs(fft(array_t_plot))));  
% ff = -(theta_total/(dis_t/L_f))/2*dis_t:dis_t:-(theta_total/(dis_t/L_f))/2*dis_t+(n-1)*dis_t;
ff = 1:theta_total/(dis_t/L_f);
f = (ff - (max(ff) + min(ff))/2 ) / (length(ff) *dis_t);
f=f*lambda*L_f/n_si;
df=f(2)-f(1);

%����˿ڵĻ�Ĥ�ͺ���(������ʹ�ø����˿ڵ��໥����ϵ�������Ӱ��˿ڼ�����ж�)
for  n=1:theta_total/(dis_out/L_f)
    for i=1:N_out
         array_out_s{n}(i) = array_out{i}(n);
    end
    array_out_s{n}=sum(array_out_s{n});
    array_out_plot(n)=array_out_s{n};
end
n1=0;
%��������˿ڵĻ�Ĥ�ͺ����ɵ�ѡȡ��ʹ��ȡ��ļ���ͽ����˸���Ҷ�任���U��ȡ����һ�£���������
for n=1:round(df/dis_out):theta_total/(dis_out/L_f)
    n1=n1+1;
    array_out_p(n1)=array_out_plot(n);
end
x_f_n=n;
x_fsr=-(theta_total/(dis_out/L_f))/2*dis_out:round(df/dis_out)*dis_out:-(theta_total/(dis_out/L_f))/2*dis_out+(n-1)*dis_out;
%��������˿ڵĻ�Ĥ�ֱ���вɵ�ѡȡ��ʹ��ȡ��ļ���ͽ����˸���Ҷ�任���U��ȡ����һ�£���������

array_out1=cell(N_out,1);
for i=1:N_out
    n1=0;
    for n=1:round(df/dis_out):theta_total/(dis_out/L_f)
        n1=n1+1;
        array_out1{i}(n1)=array_out{i}(n);
    end

end

%��������Բ������ⳡ���������ѡȡ��ʹ�ú�����˵Ļ�Ĥ��������һ�£���һ������
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


%����Բ������ⳡ�͸�������˿ڻ�Ĥ����ϣ��Ǻͺ�����

i=N_out_order;
    n2=0;
    for n=(x_c-floor(n1/2)):(x_c+floor(n1/2));
        n2=n2+1;
        if(n2>np)
           array_out1{i}(n2)=0;
           x_fsr=-(theta_total/(dis_out/L_f))/2*dis_out:round(df/dis_out)*dis_out:-(theta_total/(dis_out/L_f))/2*dis_out+(x_f_n-1)*dis_out+round(df/dis_out)*dis_out;
        end
        array_out_U(n2)=U(n);
        f_out1(n2)=array_out_U(n2)*array_out1{i}(n2); %ʵ�ʵ�����˵�����ⳡ������Բ������ⳡ�͸�������˿ڻ�Ĥ�����
    end
    
%�����������������ⳡ���й�һ������
% f_out1 = f_out1/sum(f_out1);




    


