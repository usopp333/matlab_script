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
 b4 = 0.5;
 z4 = 2;
for b3 = 0.5
    for b2 = 0.5
        for b1 = 0.5
            for z3 = 1.5
                for z2 = 1
                    for z1 = 0
                        for a4 = 1
                        for a3 = 1
                            f11 = @(z)exp(-(z-z1)^2/b1^2);
                            f22 = @(z)exp(-(z-z2)^2/b2^2);
                            f3 = @(z)a3*exp(-(z-z3)^2/b3^2);
                            f4 = @(z)a4*exp(-(z-z4)^2/b4^2);
                            a2 = (xt - (f3(zt)+f4(zt)  )- f11(zt)/f11(z0)*(x0-(f3(z0)+f4(z0)  )))/(f22(zt)-f22(z0)*f11(zt)/f11(z0));
%                             a2 = 1;
                            a1 = (x0 - a2*f22(z0)-(f3(z0)+f4(z0)  ))/f11(z0);
%                             a1 = 1;
                            if (a2 == Inf||a2==-Inf)
                                continue
                            end
                            for n = 1:nz
                                x_z(n) = a1*exp(-(z(n)-z1)^2/b1^2)+a2*exp(-(z(n)-z2)^2/b2^2)+a3*exp(-(z(n)-z3)^2/b3^2)+a4*exp(-(z(n)-z4)^2/b4^2);
                            end
                            plot(z,x_z,'color',[rand,rand,rand]);hold on
                        end
                        end                   
                    end
                end
            end           
        end
    end
end  

plot(z,line,'r','LineWidth',2)
plot([0,mmi_l],[0.25,0.25],'black--','LineWidth',2)
plot([0,mmi_l],[w+gap/2,w+gap/2],'black--','LineWidth',2)
title('gaussian');

gap = 0.2;
w = 0.5;
mmi_l = 2;
z0 = 0;
x0 = w/2;
zt = mmi_l;
xt = gap/2+w


