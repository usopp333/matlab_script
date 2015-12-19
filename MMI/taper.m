%% 设置taper函数，现在使用的是a*x^3+b*x^2+c*x+b的模型
close all
clear all
clc
taper_l = 2;
taper_w = 2;
wg_w = 0.5;

np = 1000;
% a = linspace(0,4,np);
% b = linspace(0,4,np);
x = linspace(0,taper_w/2,np);
z_taper = linspace(0,taper_w/2,np);
figure;
for a = 0.1:0.1:0.5
    for b = 0.1:0.1:2
        for c = 0.1:0.1:1
            d = 2*(taper_l - (a*(taper_w^4-wg_w^4)/16+b*(taper_w^3-wg_w^3)/8+c*(taper_w^2-wg_w^2)/4))/(taper_w-wg_w);
            e = -(a*((wg_w/2)^4)+b*((wg_w/2)^3)+c*(wg_w/2)^2+d*(wg_w/2));
            for n = 1:np
            z_taper(n) = a * x(n)^4 + b * x(n)^3 + c * x(n)^2 + d*x(n)+e;
            end
        end
        plot(z_taper,x,'color',[rand rand rand]);hold on
    end
end
title('taper函数');


%% 设置一个Y-junction使用函数曲线进行优化处理
clc

w = 0.5;
gap = 0.2;
mmi_l = 2;

%设置函数的类型为使用了论文里面的曲线的多项式拟合的函数

nz = 100;
z = linspace(0,mmi_l,nz);
x_z = linspace(0,w+gap/2,nz);
figure;
for a = -0.4577:0.0001:-0.4572 %-0.4574%
    for b = 2.8550:0.0001:2.8555 %2.8552%
        for c = -6.3251:0.0001:-6.3246 %-6.3248%
            for d = 5.5617:0.0001:5.5621 %5.5619%
                for e = -1.3213:0.0001:-1.3208 %-1.3210%
                    g = w/2;
                    f = (w+gap/2-a*(mmi_l)^6-b*(mmi_l)^5-c*(mmi_l)^4-d*(mmi_l)^3-e*(mmi_l)^2-g)/mmi_l; 
                    for n = 1:nz
                        x_z(n) = a*z(n)^6+b*z(n)^5+c*z(n)^4+d*z(n)^3+e*z(n)^2+f*z(n)+g;
                    end
                    plot(z,x_z,'color',[rand,rand,rand]);hold on
               end
            end
        end
    end 
end
plot([0,mmi_l],[0.25,0.25],'black--')
plot([0,mmi_l],[w+gap/2,w+gap/2],'black--')
title('使用了论文里面的曲线的多项式拟合函数的近似参数');

% 使用的函数参数是任意的参数
clc
clear all
w = 0.5;
gap = 0.2;
mmi_l = 2;
nz = 100;
z = linspace(0,mmi_l,nz);
x_z = linspace(w/2,w+gap/2,nz);
z_t = linspace(0,mmi_l,nz);
x_z_t = linspace(w/2,w+gap/2,nz);
line = linspace(w/2,w+gap/2,nz);
figure;
for n = 1:nz
    line(n) = (w/2+gap/2)/mmi_l*z(n)+w/2;
end
k = (w/2+gap/2)/mmi_l;
b = w/2;

for a = -0.4574 %-0.4574%
    for b = 2.8552 %2.8552%
        for c = -6.3248 %-6.3248%
            for d = 5.617:0.001:5.620 %5.5619%
                for e = -1.3212:0.01:-1.3209 %-1.3210%
                    num = 0;
                    g = w/2;
                    f = (w+gap/2-a*(mmi_l)^6-b*(mmi_l)^5-c*(mmi_l)^4-d*(mmi_l)^3-e*(mmi_l)^2-g)/mmi_l;                   
                    for n = 1:nz 
                        x_z(n) = a*z(n)^6+b*z(n)^5+c*z(n)^4+d*z(n)^3+e*z(n)^2+f*z(n)+g;
                       if((x_z(n)-line(n)) <= 0)
                           num = num+1;
                           if(num >= 0.9*nz)
                               for m = 1:nz
                                   z_t(m) = ((k-1/k)*z(m)+2*b-2*x_z(m))/(-1/k-k);
                                   x_z_t(m) = -1/k*(z_t(m)-z(m))+x_z(m);
                                   %x_z_t(m) = a*z(m)^6+b*z(m)^5+c*z(m)^4+d*z(m)^3+e*z(m)^2+f*z(m)+g;
                               end
                           end
                       end                         
                    end
                    plot(z,x_z,'color',[rand,rand,rand]);hold on
                   % plot(z_t,x_z_t,'color',[rand,rand,rand]);hold on
               end
            end
        end
    end 
end
plot(z,line,'r','LineWidth',2)
plot([0,mmi_l],[0.25,0.25],'black--')
plot([0,mmi_l],[w+gap/2,w+gap/2],'black--')
title('使用了任意变化的参数');

%% 测试函数 a/((z-z1)^2+b)....
% 使用的函数参数是任意的参数
clc
clear all
w = 0.5;
gap = 0.2;
mmi_l = 2;
nz = 100;
z = linspace(0,mmi_l,nz);
x_z = linspace(w/2,w+gap/2,nz);
z_t = linspace(0,mmi_l,nz);
x_z_t = linspace(w/2,w+gap/2,nz);
line = linspace(w/2,w+gap/2,nz);
figure;
for n = 1:nz
    line(n) = (w/2+gap/2)/mmi_l*z(n)+w/2;
end
% k = (w/2+gap/2)/mmi_l;
% b = w/2;
%首先使用三个均匀间隔极值点的曲线20151202
%三个极值点的位置
%现将所有b = 1;
z1 = mmi_l/6;
z2 = mmi_l/2;
z3 = 5*mmi_l/6;
%for z2 = mmi_l/6
for a1 = -0.5:0.1:0 
    a3 = (w+gap/2-w/2*(z2^2+1)/((mmi_l-z2)^2+1)-a1*(1/((mmi_l-z1)^2+1)-(z2^2+1)/(z1^2+1)/((mmi_l-z2)^2+1)))/(1/((mmi_l-z3)^2+1)-(z2^2+1)/(z3^2+1)/((mmi_l-z2)^2+1));
    a2 = (w/2-a3/(z3^2+1)-a1/(z1^2+1))*(z2^2+1);
           for n = 1:nz 
               x_z(n) = a1/((z(n)-z1)^2+1)+a2/((z(n)-z2)^2+1)+a3/((z(n)-z3)^2+1);
           end
            plot(z,x_z,'color',[rand,rand,rand]);hold on
%         end
%     end 
end
%end
plot(z,line,'r','LineWidth',2)
plot([0,mmi_l],[0.25,0.25],'black--')
plot([0,mmi_l],[w+gap/2,w+gap/2],'black--')
title('a/((z-z1)^2+b)....');

%% 使用高斯函数
% 使用的函数参数是任意的参数
clc
clear all
w = 0.5;
gap = 0.2;
mmi_l = 2;
z0 = 0;
x0 = w/2;
zt = mmi_l;
xt = w+gap/2;
nz = 100;
z = linspace(0,mmi_l,nz);
x_z = linspace(w/2,w+gap/2,nz);
z_t = linspace(0,mmi_l,nz);
x_z_t = linspace(w/2,w+gap/2,nz);
line = linspace(w/2,w+gap/2,nz);
figure;
for n = 1:nz
    line(n) = (w/2+gap/2)/mmi_l*z(n)+w/2;
end
% x_z(n) = a1*exp(-(z(n)-z1)^2/b1^2)+a2*exp(-(z(n)-z2)^2/b2^2)+a3*exp(-(z(n)-z3)^2/b3^2)+f4+f5....;
for b3 = 1
    for b2 = 1
        for b1 = 1
            for z3 = 0:1:2
                for z2 = 1
                    for z1 = 1
                        for a3 = -0.5:0.1:0.5
                            f11 = @(z)exp(-(z-z1)^2/b1^2);
                            f22 = @(z)exp(-(z-z2)^2/b2^2);
                            f3 = @(z)a3*exp(-(z-z3)^2/b3^2);
                            a2 = (xt - (f3(zt)  )- f11(zt)/f11(z0)*(x0-(f3(z0)  )))/(f22(zt)-f22(z0)*f11(zt)/f11(z0));
                            a1 = (x0 - a2*f22(z0)-(f3(z0)  ))/f11(z0);
                            if (a2 == Inf)
                                continue;
                            end
                            for n = 1:nz
                                x_z(n) = a1*exp(-(z(n)-z1)^2/b1^2)+a2*exp(-(z(n)-z2)^2/b2^2)+a3*exp(-(z(n)-z3)^2/b3^2);
                            end
                            plot(z,x_z,'color',[rand,rand,rand]);hold on
                        end
                    end
                end
            end
        end
    end
end          
plot(z,line,'r','LineWidth',2)
plot([0,mmi_l],[0.25,0.25],'black--')
plot([0,mmi_l],[w+gap/2,w+gap/2],'black--')
title('gaussian');


%% 将论文中的曲线进行多项式的拟合寻找合适的优化函数

clc
%close all
nx = 100;
x = linspace(0,2,13);
y = [0.5,0.5,0.6,0.7,0.9,1.26,1.4,1.4,1.4,1.4,1.31,1.2,1.2];
% y1 = (1-0.25)*rand(1,(nx-2));
% y = [0.25,y1,(0.5+0.1)];
 
figure;
plot(x,y/2,'r*');hold on;
for m = 6
    P=polyfit(x,y/2,m);
    xi=0:2/nx:2;
    yi=polyval(P,xi);
    plot(xi,yi,'color',[rand,rand,rand]);hold on
    title(num2str(m));
end



%% 使用自定义函数进行拟合
clear all
X=[65.31 44.88 39.92 34.73 29.51 24.17 18.20 12.33 6.02 0.64];Y=[0.98 2.02 2.89 2.87 2.81 5.02 5.11 4.75 5.584 5.17];
%fun=inline('a(1)*X+a(2)./X+a(3)','a','X');
fun=inline('a(1)+a(2)./(X^2)+a(3)./(X^4)','a','X');
b=[10 10 1];
a = nlinfit(X,Y,fun,b);
vpa(a,10)
A=a(1),B=a(2),C=a(3)
figure;
plot(X,Y)


nx = 100;
clear all
%x = linspace(0,2,13);
x = [0.00001 0.166666666666667 0.333333333333333 0.500000000000000 0.666666666666667 0.833333333333333 1 1.16666666666667 1.33333333333333 1.50000000000000 1.66666666666667 1.83333333333333 2];
y = [0.5 0.5 0.6 0.7 0.9 1.26 1.4 1.4 1.4 1.4 1.31 1.2 1.2];
fun = inline('A(1)+A(2)/x.^1+A(3)/x.^2','A','x');
b = [2 1 1];
A = nlinfit(x,y,fun,b);
vpa(A,13)
a = A(1),b = A(2),c = A(3)
figure;
plot(x,y/2,'r*');hold on;

%% 测试函数 a/((z-z1)^2+b)....
x = linspace(0,10,100);
y = linspace(0,10,100);
for n = 1:100
   %y(n) = 1/((x(n)-1)^2+1)+1/((x(n)-3)^2+1)+1/((x(n)-5)^2+1)+1/((x(n)-7)^2+1+1)/((x(n)-9)^2+1);
    y(n) = 1/((x(n)-1)^2+1);
end

for n = 1:100
   %y(n) = 1/((x(n)-1)^2+1)+1/((x(n)-3)^2+1)+1/((x(n)-5)^2+1)+1/((x(n)-7)^2+1+1)/((x(n)-9)^2+1);
    y1(n) = 2/((x(n)-1)^2+1);
end

for n = 1:100
   %y(n) = 1/((x(n)-1)^2+1)+1/((x(n)-3)^2+1)+1/((x(n)-5)^2+1)+1/((x(n)-7)^2+1+1)/((x(n)-9)^2+1);
    y2(n) = 1/((x(n)-1)^2+2);
end

for n = 1:100
   %y(n) = 1/((x(n)-1)^2+1)+1/((x(n)-3)^2+1)+1/((x(n)-5)^2+1)+1/((x(n)-7)^2+1+1)/((x(n)-9)^2+1);
    y3(n) = 0.5/((x(n)-1)^2+1);
end

figure;
plot(x,y, x, y1, x, y2, x, y3);

x = linspace(0,10,100);
y = linspace(0,10,100);
A = 1;
b = 1;
figure;
for b = 1:1:5
    for A = 1
        for n = 1:100
   %y(n) = 1/((x(n)-1)^2+1)+1/((x(n)-3)^2+1)+1/((x(n)-5)^2+1)+1/((x(n)-7)^2+1+1)/((x(n)-9)^2+1);
    y(n) = A*b^2/((x(n)-1)^2+b^2);
        end
    end
    plot(x,y,'color',[rand,rand,rand]);hold on
end

        













