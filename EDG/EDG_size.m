%日期：2015-11-12
%进行EDG的设计
%以下各单位暂定为um
%这个中科院的设计是基于罗兰圆的设计，有的结构的设计是基于平场，但是那是光纤结构可以省去输入输出波导，硅光结构不用考虑；
%关于EDG光栅面的设计基本上文献中有两种方法：一点法和两点法；
%第一步：将中科院的EDG的模型进行复现（不能完全仿造因为有些结构只能靠自己判断）
%以上工作由于当时任务侧重点不在EDG暂时搁置：

%日期：2016-01-28
%由于海信公司的需要：（1）对中科院的EDG结构进行相应的优化
%                     (2)书写计算方法可以自定的生成EDG的各个的点坐标

clc
clear all
n_si=3.472000;  %index of si
n_sio2=1.444000;  %index of sio2
neff_factor = 0; 
%n_s_eff=2.959275 + neff_factor;%使用了cloudfdtd的精确计算（0.001um）
n_s_eff = 2.932935; % 使用0.05um的网格计算
n_b_eff = 2.632102; %使用了cloudfdtd3D的精确计算
%ng_b=3.753849301034707;  % 使用了cloudfdtd3D的精确计算并进行了相应的计算
pi = 3.1415926;
%% 设计参数的设置：
clambda = 1.310;  %器件的中心波长
dlambda = 0.02;   %信道宽度
N_in  = 1;  %输入通道
N_out = 4; %输出通道
theta_i = pi/4; %输入通道与大罗兰圆轴线的夹角（作为EDG的入射角）
d_io = 2;  %输入输出波导的芯芯间隔（关于输出端口的位置设置还不太清楚,暂且设置为输入和输出相邻间隔都是相同的）
wg_io = 1.85;
wg_width = 0.500;  %输入输出波导的宽度
taper_l = 12;  %输入输出taper的长度
wg_line = 5;
%% 计算群折射率
% lambda (1-20)
x_lambda = [1.200000e+000, 1.210526e+000, 1.221053e+000, 1.231579e+000, 1.242105e+000, 1.252632e+000, 1.263158e+000, 1.273684e+000, 1.284211e+000, 1.294737e+000, 1.305263e+000, 1.315789e+000, 1.326316e+000, 1.336842e+000, 1.347368e+000, 1.357895e+000, 1.368421e+000, 1.378947e+000, 1.389474e+000, 1.400000e+000];

% neff - neff(Abs) ((1e0)) (1-20)
neff_lambda = [2.746509e+000, 2.735879e+000, 2.725181e+000, 2.714418e+000, 2.703588e+000, 2.692690e+000, 2.681725e+000, 2.670693e+000, 2.659591e+000, 2.648421e+000, 2.637183e+000, 2.625874e+000, 2.614495e+000, 2.603047e+000, 2.591528e+000, 2.579937e+000, 2.568275e+000, 2.556541e+000, 2.544734e+000, 2.532854e+000];
p=polyfit(x_lambda,neff_lambda,1);
x1=linspace(min(x_lambda),max(x_lambda));
y1=polyval(p,x1);
% figure;
% plot(x_lambda,neff_lambda,'*',x1,y1);
% title('index 1310 width=0.5');
% ng_b = n_si - clambda*(abs(p(1)));


%% 相关结构参数的计算：
L_f = clambda * d_io/ (2*dlambda*sin(theta_i)); %大罗兰圆的半径
theta_k = linspace(0,90,N_out); %反射角（输出端口与轴线的夹角）
%theta_k(1) = ( theta_i - 2*asin(d_io/2/(L_f/2*cos(theta_i))) ); %第一个输出波导与大罗兰圆轴线的夹角，距离输入波导最近的输出波导
theta_k(1) = ( theta_i - 2*asin(d_io/2/(L_f*cos(theta_i))) ); %第一个输出波导与大罗兰圆轴线的夹角，距离输入波导最近的输出波导
for k=2:N_out
      %theta_k (k)= ( theta_k (k-1)- 2*asin(d_io/2/(L_f/2*cos(theta_k(k-1)))) ); 
      theta_k (k)= ( theta_k (k-1)- 2*asin(d_io/2/(L_f*cos(theta_k(k-1)))) ); 
end
%theta_ck = (theta_k(4) - theta_k(1))/2 + theta_k(1);  %输出波导的中心位置与大罗兰圆轴线的夹角（可作为整个EDG的衍射角by：lzq）
theta_ck = theta_k(2);  %把整个edg的衍射角的中心放在第二个中心波导（可作为整个EDG的衍射角by：lzq）
out_k = zeros(4,2); %输出波导中心在小罗兰圆上的坐标
for k = 1:N_out
    out_k(k,:) = [0-L_f*cos(theta_k(k))*sin(theta_k(k)),L_f-L_f*cos(theta_k(k))*cos(theta_k(k))];
end
fprintf('\n#入射角theta_i\n');
fprintf('theta_i = ');
fprintf(1,'%f',theta_i);
fprintf(';');
fprintf('\n#反射角theta_k\n'); 
fprintf('theta_k = [');
fprintf(1,'%f,',theta_k);
fprintf('];');
in_i = [-L_f/2,L_f/2];
close all
figure;
plot(out_k(:,1),out_k(:,2),'green+','LineWidth',4)
hold on
%L_f = clambda * d_io/ (dlambda*(sin(theta_i)+sin(theta_ck))); %大罗兰圆的半径
%fsr = N_out * dlambda;  %结合AWG的设计，FSR的至少应当是信道宽度和输出通道数目的乘积
fsr = 0.2;
%m = ceil(clambda / fsr) ;  %根据公式FSR=lambda/m
m = ceil(clambda / fsr)+1 ;  %根据公式FSR=lambda/m
d = m*clambda/((sin(theta_i)+sin(theta_ck))*n_s_eff); % 计算的光栅周期
k_wg = 2*pi/clambda;
Vc = k_wg*d_io*sqrt(n_si^2-n_sio2^2);  %归一化频率 中心频率
w0c = d_io*(0.321+2.1*Vc^(-1.5)+4*Vc^(-6));  %高斯波的W0，在输入taper的输入高斯的模场半径中心频率
%theta_0=atan(clambda/(pi*w0c));  %利用的是高斯光的光束发散角公式
theta_0 = 2*(clambda/(pi*w0c));  %利用的是高斯光的光束发散角公式
theta_total = theta_0;  %入射宽口EDG展开角度（使用的方法是和AWG的方式一致）
%theta_total = theta_0;  %入射宽口EDG展开角度（使用的方法是和AWG的方式一致）
%中科院EDG文档中说明的是根据罗兰圆的半径和光栅周期可以估算出中心光栅左右两侧的光栅数目；但是基本上的EDG的光栅周期都不是常数，况且说的太笼统，无法估算
%% 绘制相应的图形来确定罗兰圆上的光栅中心点坐标
%关于小罗兰圆的画法使用角度的绘制方法，可以使得小罗兰圆上的点坐标分布均匀，使得扫描的结果较好
nz = 10000; %计算步数 5000
alpha=0:2*pi/nz:2*pi;%角度[0,2*pi] 
%R=2;%半径 
small_circle = zeros(length(alpha),2);
x_half = (L_f/2)*cos(alpha); 
y_half = (L_f/2)*sin(alpha)+L_f/2;
small_circle(:,1) = x_half;
small_circle(:,2) = y_half;
plot(small_circle(:,1),small_circle(:,2),'r.')
%plot(x_half,y_half+L_f/2,'k') % 使用弧度函数绘制小罗兰圆
%syms z x_h;
z = linspace(-L_f,L_f,nz);
%x = linspace(-L_f,L_f,1000);
z_in_life = linspace(-L_f,L_f,nz);
z_in_right = linspace(-L_f,L_f,nz);
%x_h = solve('z^2+(x_h-L_f/2)^2=(L_f/2)^2','x_h');
%big_circle = linspace(0,L_f,nz);
big_circle = @(z) sqrt(L_f^2-z^2);
% 使用点坐标的函数来绘制小罗兰圆的时候会导致在小罗兰圆和水平轴线的连接处圆上的点坐标较少，使得扫描输入输出taper点位置时，误差较大
small_circle_up = @(z)(real(sqrt((L_f./2)^2-z^2)+L_f./2));
small_circle_down = @(z)(real(-sqrt((L_f./2)^2-z^2)+L_f./2));
%big_circle = @(z)(sqrt(L_f^2-z^2));
half_circle = linspace(0,L_f,nz);
% in_life =  linspace(0,L_f,nz);
% in_right =  linspace(0,L_f,nz);
big_circle_p = zeros(nz,2) ;
small_circle_up_p = zeros(nz,2);
small_circle_down_p = zeros(nz,2);
for i = 1:nz
    big_circle_p(i,:) = [z(i),big_circle(z(i))];
    small_circle_up_p(i,:) = [z(i),small_circle_up(z(i))];
    small_circle_down_p(i,:) = [z(i),small_circle_down(z(i))];
    %big_circle(i) = big_circle(z(i));
    %in_life(i) = (z(i)+L_f/2)*tan(theta_total/2+theta_i)+L_f/2;
end
in_life = (z_in_life+L_f/2)*tan(theta_total/2+pi/2-theta_i)+L_f/2;
%in_right = (z_in_right+L_f/2)*tan(theta_total/2-theta_i)+L_f/2;
if (theta_i>theta_total/2)
    in_right = (z_in_right+L_f/2)*tan(pi/2-theta_i-theta_total/2)+L_f/2;
end
if(theta_i<theta_total/2)
    in_right = (z_in_right+L_f/2)*tan(pi-(theta_total/2-(pi/2-theta_i)))+L_f/2;
end
out_line = cell(N_out,1); % 输出波导与光栅中心的连线
in_line = zeros(nz,2); % 输入波导与光栅中心的连线
for k = 1:N_out
    for j = 1:nz
        in_line(j,:) = [z(j),(z(j)-(-L_f/2))*((L_f-L_f/2)/(0-(-L_f/2)))+L_f/2];
        out_line{k}(j,1) = z(j);
        out_line{k}(j,2) = (z(j)-out_k(k,1))*((L_f-out_k(k,2))/(0-out_k(k,1)))+out_k(k,2);
    end
    plot(out_line{k}(:,1),out_line{k}(:,2)) %输出波导与光栅中心的连线
end
plot(in_line(:,1),in_line(:,2),'r') %输入波导与光栅中心的连线
% 求输入光夹角和大罗兰圆的交点坐标
for i = 1:nz
    dif_l =abs( big_circle(z(i))-in_life(i));
    dif_r = abs( big_circle(z(i))-in_right(i));
    if(dif_l<0.1);
        %fprintf(1, 'z_f=%f , x_f=%f\n',z(i),big_circle(z(i)));
        z_life = z(i);
        x_life = big_circle(z(i));
        life_edg = [z(i),big_circle(z(i))] ;
    end
    if(dif_r<0.1 && z(i)>0);
        %fprintf(1, 'z_r=%f , x_r=%f\n',z(i),big_circle(z(i)));
        z_right = z(i);
        x_right = big_circle(z(i));
        right_edg = [z(i),big_circle(z(i))] ;
    end
end
hold on
%plot(z_life,x_life,'red+','LineWidth',2)
%plot(z_right,x_right,'red+','LineWidth',2)
plot(life_edg(1,1),life_edg(1,2),'black*')
plot(right_edg(1,1),right_edg(1,2),'black*')
%theta_grating_center = (theta_i - theta_ck)/2+theta_ck;
%theta_d = 2*asin(d/2/L_f);  % 中心光栅的光栅面的长度
theta_d = d/L_f;   % 中心光栅的光栅面的长度
%theta_d = d*cos(theta_grating_center);  % 中心光栅的光栅面的长度
%光栅中心位置 计算简单起见，以EDG的中心轴线的光栅中点为起点，左右进行计算；
in_center = [-L_f/2,L_f/2]; %输入波导中心坐标
out_center = [-L_f*sin(theta_ck)*cos(theta_ck),L_f*sin(theta_ck)*sin(theta_ck)]; %输出波导中心坐标
%out_center = out_k(2,:); %输出波导中心 设置到第二个输出波导的中心
g_center = [0,L_f]; %中心轴线光栅的中心位置
%in_center
theta_life = asin(abs(z_life)/L_f);  % 左边展开角
theta_right = asin(abs(z_right)/L_f); %右边展开角
g_life_num =ceil( theta_life/theta_d); %左边光栅个数
g_right_num = ceil(theta_right/theta_d); %右边光栅个数

%% 1.左边光栅位置的计算 利用光栅公式 Ln - Ln-1 = m*lambda/neff;
l_l = linspace(0,1,g_life_num);
l = linspace(-L_f,L_f,nz/2);
%g_life_num = g_life_num-1; % 此处的-1是根据下面计算的结果，上面定义的右端的光栅的点数有些点是不满足光栅方程的，所以进行了相应的点缩放，要根据实际情况而定
l_f = @(z)(sqrt((in_center(1,1)-z)^2+(in_center(1,2)-big_circle(z))^2)+sqrt((out_center(1,1)-z)^2+(out_center(1,2)-big_circle(z))^2));
l_l(1) = l_f(g_center(1,1));
g_center_l = zeros(g_life_num,2) ;
%扫描大罗兰圆上的点坐标找到符合光栅方程的点
for j = 1:g_life_num
    for i = 1:(nz/2)
        l(i) = l_f(z(i));
        d_l = l_l(j)-l(i);
        %if (abs(d_l-m*clambda/n_s_eff)<0.04) % 2016-2-9
        if (abs(d_l-m*clambda/n_s_eff)<0.05)
            l_l(j+1) = l(i);
            %if(abs(l_l(j+1)-l_l(j))>1)
                if(big_circle(z(i))>0)
            %fprintf(1, 'l_l=%f',l_l);
            %fprintf(1, 'l_l=%f ,z = %f ,x = %f\n',l_l(j),z(i),big_circle(z(i)));
            g_center_l(g_life_num-j+1,:) = [z(i),big_circle(z(i))];
                end
           % end
        end
    end
end
plot(g_center_l(:,1),g_center_l(:,2),'red+','LineWidth',2)
% fprintf('\n左边光栅中心点(z)g_center_l(:,1)\n');
% fprintf(1, '%f,',g_center_l(:,1));
% fprintf('\n左边光栅中心点(x)g_center_l(:,2)\n');
% fprintf(1, '%f,',g_center_l(:,2));
%计算左边光栅中心点的光栅面和水平轴线的夹角
%根据中科院的PPT，要使得光栅面和输入光线垂直
%这里我的计算是 要使得光栅面和输入波导中心与输出波导中心的角平分线垂直
theta_g_i = linspace(0,2*pi,g_life_num);
theta_g_k = linspace(0,2*pi,g_life_num);
theta_p_l = linspace(0,2*pi,g_life_num+1);
g_center_l(g_life_num+1,:) =  g_center;
for i = 1:g_life_num
    theta_g_i(i) = atan(abs(g_center_l(i,1)-in_center(1,1))/abs(g_center_l(i,2)-in_center(1,2)));
    theta_g_k(i) = atan(abs(g_center_l(i,1)-out_center(1,1))/abs(g_center_l(i,2)-out_center(1,2)));
    theta_p_l(i) = (theta_g_i(i)-theta_g_k(i))/2+theta_g_k(i);
    %fprintf(1, 'theta_p_l=%f\n',theta_p_l);
end
% fprintf('\n左边光栅倾斜角theta_p_l\n');
% fprintf(1, '%f,',theta_p_l);
%计算左边光栅面的长度；
%根据每个光栅的光栅面要和相邻的光栅连接线垂直的几何关系可以一一推出每个光栅面的长度；
%中科院是根据相邻光栅面与大罗兰圆之间的几何关系进行的计算；
%由于上一步没有按照中科院的方法进行计算，所以这里算法也要做相应的修改；
%中心轴线的光栅面的长度根据一开始计算的光栅周期进行计算；
theta_p_center = (theta_i-theta_ck)/2+theta_ck;
d_g = linspace(0,1,g_life_num+1);
%d_g_center = d*cos(theta_p_center);
d_g_center = (d/2)/cos(theta_p_center);
%d_g(g_life_num+1) = d_g_center/2; %有一定的误差
d_g(g_life_num+1) = d_g_center; %有一定的误差
theta_p_l(g_life_num+1) = theta_p_center;
theta_g_g(g_life_num+1) = theta_p_center;
g_line = cell(g_life_num+1,1); % 左边光栅面的直线
for i = 1:g_life_num+1
    for j = 1:nz
        %g_line{j} = [z(j),(z(j)-g_center_l(i,1))*tan(pi-theta_p_l(i))+g_center_l(i,2)];
        g_line{i}(j,1) = z(j);
        g_line{i}(j,2) = (z(j)-g_center_l(i,1))*tan(pi-theta_p_l(i))+g_center_l(i,2);
    end
  %plot(g_line{i}(:,1),g_line{i}(:,2))
   
end
%利用扫描的方法找到求得光栅面垂直的连接线，既可以求出光栅面的长度
%注：并不是每个光栅面都和和它连接的线垂直，一般包保证连接线不能遮挡住入射光，所以一般是认为光栅面和左边的连接线是垂直的；
d_point_f_u = zeros(g_life_num+1,2);
d_point_f_d = zeros(g_life_num+1,2);
d_point_f_d(g_life_num+1,:) = [g_center_l(g_life_num+1,1)+d_g(g_life_num+1)/2*cos(theta_p_l(g_life_num+1)),g_center_l(g_life_num+1,2)-d_g(g_life_num+1)/2*sin(theta_p_l(g_life_num+1))];
g_link = cell(g_life_num+1,1);
for i = g_life_num+1:-1:1
    % 使用距离进行计算
    %d_point_f_u(i,:) = [g_center_l(i,1)-d_g(i)/2*cos(theta_p_l(i)),g_center_l(i,2)+d_g(i)/2*sin(theta_p_l(i))];
    d_point_f_u(i,:) = [2*g_center_l(i,1)-d_point_f_d(i,1),2*g_center_l(i,2)-d_point_f_d(i,2)];
    % 使用坐标进行计算
    %d_point_f_u(i,:) = [2*g_center_l(i,1)-d_point_f_d(i,1),2*g_center_l(i,2)-d_point_f_d(i,2)];
    for j = 1:nz
        g_link{i}(j,1) = z(j);
        %g_link{i}(j,2) = (z(j)-d_point_f_u(i,1))*tan(pi/2-theta_p_l(i))+d_point_f_u(i,2);
        % 要使得光栅面的连接面不遮挡光线
        g_link{i}(j,2) = (z(j)-d_point_f_u(i,1))*((d_point_f_u(i,2)-in_i(1,2))/(d_point_f_u(i,1)-in_i(1,1)))+d_point_f_u(i,2);
        if(i>1 && big_circle(z(j))>0)
        %if(abs(g_line{i-1}(j,2)-g_link{i}(j,2))<0.08 && abs(g_line{i-1}(j,1)-g_link{i}(j,1))<0.08)
        %if(abs(g_line{i-1}(j,2)-g_link{i}(j,2))<0.2 && abs(g_line{i-1}(j,1)-g_link{i}(j,1))<0.2)
        if(abs(g_line{i-1}(j,2)-g_link{i}(j,2))<0.05 && abs(g_line{i-1}(j,1)-g_link{i}(j,1))<0.05)
            d_g(i-1) = 2*sqrt((g_center_l(i-1,1)-g_line{i-1}(j,1))^2+(g_center_l(i-1,2)-g_line{i-1}(j,2))^2);
            % 这里使用g_link和g_line绘制出的点坐标不太相同（会导致整个光栅面的点坐标有所不同）
            % g_link保证不挡光；g_line保证准确的光栅面（扫描存在一定的误差，没有办法将两者合理的统一）
            %d_point_f_d(i-1,:) = [g_link{i}(j,1),g_link{i}(j,2)];
            d_point_f_d(i-1,:) = [g_line{i-1}(j,1),g_line{i-1}(j,2)];
        end
        end
    end
    plot(g_link{i}(:,1),g_link{i}(:,2),'r')
    k = 2*i;
    d_point_f (k-1,:)= d_point_f_u(i,:);
    d_point_f (k,:)= d_point_f_d(i,:);
    % 由于左边的光栅面和连接线使用现在的扫描方式误差较大，所以对之前定义的光栅面中心坐标和光栅面的倾斜角进行修正（右边的扫描很准确）
    g_center_l_revise(i,:) = [(d_point_f_u(i,1)+d_point_f_d(i,1))/2,(d_point_f_u(i,2)+d_point_f_d(i,2))/2]; 
    theta_p_l_revise(i) = atan(abs(d_point_f_u(i,2)-d_point_f_d(i,2))/abs(d_point_f_u(i,1)-d_point_f_d(i,1)));
    d_g_revise(i) = sqrt((d_point_f_u(i,1)-d_point_f_d(i,1))^2+(d_point_f_u(i,2)-d_point_f_d(i,2))^2);
end
%plot(d_point_f(:,1),d_point_f(:,2),'green+')



%% 1.右边光栅位置的计算 利用光栅公式 Ln - Ln-1 = m*lambda/neff;
l_l_r = linspace(0,1,g_right_num);
l = linspace(-L_f,L_f,nz/2);
l_f_r = @(z)(sqrt((in_center(1,1)-z)^2+(in_center(1,2)-big_circle(z))^2)+sqrt((out_center(1,1)-z)^2+(out_center(1,2)-big_circle(z))^2));
l_l_r(1) = l_f_r(g_center(1,1));
%g_right_num = g_right_num-10; % 此处的-10 是根据下面计算的结果，上面定义的右端的光栅的点数有些点是不满足光栅方程的，所以进行了相应的点缩放，要根据实际情况而定
g_center_r = zeros(g_right_num,2) ;
%扫描大罗兰圆上的点坐标找到符合光栅方程的点
zero_num = 0;
for j = 1:g_right_num
    %for i = nz/2-1:nz
    for i = nz/2+1:nz
        l(i) = l_f_r(z(i));
        d_l = -l_l_r(j)+l(i);
        %if (abs(d_l-m*clambda/n_s_eff)<0.04 )
        if (abs(d_l-m*clambda/n_s_eff)<0.05 )    
            if(z(i)<z_right)
            l_l_r(j+1) = l(i);
            if(abs(l_l_r(j+1)-l_l_r(j))>1)
            %fprintf(1, 'l_l_r=%f',l_l_r);
            %fprintf(1, 'l_l_r=%f ,z = %f ,x = %f\n',l_l_r(j),z(i),big_circle(z(i)));
            g_center_r(j,:) = [z(i),big_circle(z(i))];
            end
            end
        end
    end
    if(g_center_r(j,2) == 0)
        zero_num = zero_num+1;
    end
end
% 定义的g_center_r比实际上符合光栅方程的点要多，会出现多余的（0,0）点，使用下面的方法删除多余的点
for i = 1:zero_num
    g_center_r(g_right_num-zero_num+1,:)=[];
end
g_right_num = g_right_num-zero_num;

plot(g_center_r(:,1),g_center_r(:,2),'black+','LineWidth',2)
% fprintf('\n右边光栅中心点(z)g_center_r(:,1)\n');
% fprintf(1, '%f,',g_center_r(:,1));
% fprintf('\n右边光栅中心点(x)g_center_r(:,2)\n');
% fprintf(1, '%f,',g_center_r(:,2));
%计算右边光栅中心点的光栅面和水平轴线的夹角
%根据中科院的PPT，要使得光栅面和输入光线垂直
%这里我的计算是 要使得光栅面和输入波导中心与输出波导中心的角平分线垂直
theta_g_i = linspace(0,2*pi,g_right_num);
theta_g_k = linspace(0,2*pi,g_right_num);
theta_p_r = linspace(0,2*pi,g_right_num);
%g_center_r(1,:) =  g_center;
for i = 1:g_right_num
    theta_g_i(i) = atan(abs(g_center_r(i,1)-in_center(1,1))/abs(g_center_r(i,2)-in_center(1,2)));
    theta_g_k(i) = atan(abs(g_center_r(i,1)-out_center(1,1))/abs(g_center_r(i,2)-out_center(1,2)));
    theta_p_r(i) = (theta_g_i(i)-theta_g_k(i))/2+theta_g_k(i);  
end
% fprintf('\n右边光栅倾斜角theta_p_r\n');
% fprintf(1, '%f,',theta_p_r);
%计算右边光栅面的长度；
%根据每个光栅的光栅面要和相邻的光栅连接线垂直的几何关系可以一一推出每个光栅面的长度；
%中科院是根据相邻光栅面与大罗兰圆之间的几何关系进行的计算；
%由于上一步没有按照中科院的方法进行计算，所以这里算法也要做相应的修改；
%中心轴线的光栅面的长度根据一开始计算的光栅周期进行计算；
theta_p_center = (theta_i-theta_ck)/2+theta_ck;
d_g_r_center = d/cos(theta_p_center);
d_g_r = linspace(0,1,g_right_num);
%d_g_r(1) = d_g_r_center/2; %此处/2有待研究（/2效果好但是好像和我自己的计算不一致）
g_line_r = cell(g_right_num,1); % 右边光栅面的直线
for i = 1:g_right_num
    for j = 1:nz
        %g_line_r{j} = [z(j),(z(j)-g_center_r(i,1))*tan(pi-theta_p_r(i))+g_center_r(i,2)];
        g_line_r{i}(j,1) = z(j);
        g_line_r{i}(j,2) = (z(j)-g_center_r(i,1))*tan(pi-theta_p_r(i))+g_center_r(i,2);
    end
   % plot(g_line_r{i}(:,1),g_line_r{i}(:,2))
   
end
%利用扫描的方法找到求得光栅面垂直的连接线，既可以求出光栅面的长度
%注：并不是每个光栅面都和和它连接的线垂直，一般包保证连接线不能遮挡住入射光，所以一般是认为光栅面和右边的连接线是垂直的；
d_point_r_u = zeros(g_right_num,2);
d_point_r_d = zeros(g_right_num+1,2);
%d_point_r_u(1,:) = [g_center_l(g_life_num+1,1)-d_g(g_life_num+1)/2*cos(theta_p_l(g_life_num+1)),g_center_l(g_life_num+1,2)+d_g(g_life_num+1)/2*sin(theta_p_l(g_life_num+1))];
%d_point_r_d(1,:) = [g_center_r(1,1)+d_g_r(1)/2*cos(theta_p_r(1)),g_center_r(1,2)-d_g_r(1)/2*sin(theta_p_r(1))];
d_point_r_d(1,:) = d_point_f_d(g_life_num+1,:);
g_link_r = cell(g_right_num,1);
theta_g_r_link = linspace(0,1,g_right_num);
% for i = 1:g_right_num %  %利用距离的最小值来判断光栅面的上端点
%     for j = 1:nz
%         d_g_link_r(j) = sqrt((d_point_r_d(i,1)-g_line_r{i}(j,1))^2+(d_point_r_d(i,2)-g_line_r{i}(j,2))^2);
%     end
%     d_g_link_r_min(i) = min(d_g_link_r); %利用距离的最小值来判断光栅面的上端点，使得光栅面和右边的连接线是垂直的
%     %a = min(d_g_link_r);
%     min_index = find(d_g_link_r==d_g_link_r_min(i));
%     d_point_r_u(i,:) = [g_line_r{i}(min_index,1),g_line_r{i}(min_index,2)];
%     d_point_r_d(i+1,:) = [(g_center_r(i,1)*2-d_point_r_u(i,1)),(g_center_r(i,2)*2-d_point_r_u(i,2))];
%     d_g_r(i) = sqrt((d_point_r_u(i,1)- d_point_r_d(i+1,1) )^2+(d_point_r_u(i,2)- d_point_r_d(i+1,2))^2);
%     theta_g_r_link(i) = (atan(( (g_line_r{i}(min_index,2)-d_point_r_d(i,2))/(g_line_r{i}(min_index,1)-d_point_r_d(i,1)))));
%     g_link_r{i}(j,1) = z(j);
%     g_link_r{i}(j,2) = (z(j)-d_point_r_u(i,1))*tan(theta_g_r_link(i))+d_point_r_u(i,2);
%     % plot(g_link_r{i}(:,1),g_link_r{i}(:,2),'r')
%     k = 2*i;
%     d_point_r (k-1,:)= d_point_r_d(i,:);
%     d_point_r (k,:)= d_point_r_u(i,:);
%     % plot(g_link_r{i}(:,1),g_link_r{i}(:,2))
% end
for i = 1:g_right_num % 要使得光栅面的连接面不遮挡光线
    for j = 1:nz
        g_link_r{i}(j,1) = z(j);
        g_link_r{i}(j,2) = (z(j)-d_point_r_d(i,1))*((d_point_r_d(i,2)-in_i(1,2))/(d_point_r_d(i,1)-in_i(1,1)))+d_point_r_d(i,2); 
        %if(abs(g_line_r{i}(j,2)-g_link_r{i}(j,2))<0.15 && abs(g_line_r{i}(j,1)-g_link_r{i}(j,1))<0.15)
        if(abs(g_line_r{i}(j,2)-g_link_r{i}(j,2))<0.05 && abs(g_line_r{i}(j,1)-g_link_r{i}(j,1))<0.05)
            % 这里使用g_link和g_line绘制出的点坐标不太相同（会导致整个光栅面的点坐标有所不同）
            % g_link保证不挡光；g_line保证准确的光栅面（扫描存在一定的误差，没有办法将两者合理的统一）
            %d_point_r_u(i,:) = [g_link_r{i}(j,1),g_link_r{i}(j,2)];
            d_point_r_u(i,:) = [g_line_r{i}(j,1),g_line_r{i}(j,2)];
            d_point_r_d(i+1,:) = [(g_center_r(i,1)*2-d_point_r_u(i,1)),(g_center_r(i,2)*2-d_point_r_u(i,2))];
            d_g_r(i) = sqrt((d_point_r_u(i,1)-d_point_r_d(i+1,1))^2+(d_point_r_u(i,2)-d_point_r_d(i+1,2))^2);
        end
    end
    plot(g_link_r{i}(:,1),g_link_r{i}(:,2),'r')
    k = 2*i;
    d_point_r (k-1,:)= d_point_r_d(i,:);
    d_point_r (k,:)= d_point_r_u(i,:);
    g_center_r_revise(i,:) = [(d_point_r_u(i,1)+d_point_r_d(i+1,1))/2,(d_point_r_u(i,2)+d_point_r_d(i+1,2))/2]; 
    theta_p_r_revise(i) = atan(abs(d_point_r_u(i,2)-d_point_r_d(i+1,2))/abs(d_point_r_u(i,1)-d_point_r_d(i+1,1)));
end
        
        
d_point_r(k+1,:) = d_point_r_d(i+1,:);
%plot(d_point_r(:,1),d_point_r(:,2))

%d_g_all = [d_g;d_g_r];
d_point_all = [d_point_f;d_point_r]; %EDG的阶梯光栅的点数据
plot(d_point_all(:,1),d_point_all(:,2))
% fprintf('\nd_point_all(:,1)\n');
% fprintf(1,'%f,',d_point_all(:,1));
% fprintf('\nd_point_all(:,2)\n');
% fprintf(1,'%f,',d_point_all(:,2));

g_center_all_point = [g_center_l;g_center_r];
g_center_all_revise = [g_center_l_revise;g_center_r_revise];
theta_p_all = [theta_p_l,theta_p_r];
theta_p_all_revise = [theta_p_l_revise,theta_p_r_revise];
d_g_all = [d_g d_g_r]; % 光栅面的长度
d_g_all_revise = [d_g_revise d_g_r];
% fprintf('\n光栅中心点(z)g_center_all_point(:,1)\n');
% fprintf(1, '%f,',g_center_all_point(:,1));
% fprintf('\n光栅中心点(x)g_center_all_point(:,2)\n');
% fprintf(1, '%f,',g_center_all_point(:,2));
% fprintf('\n光栅倾斜角theta_p_all\n');
% fprintf(1, '%f,',theta_p_all);
% fprintf('\n光栅面的长度d_g_all\n');
% fprintf(1, '%f,',d_g_all);
fprintf('\n#修正的光栅中心点(z)g_center_all_revise(:,1)\n');
fprintf('g_center_all_z = [');
fprintf(1, '%f,',g_center_all_revise(:,1));
fprintf('];');
fprintf('\n#修正的光栅中心点(x)g_center_all_revise(:,2)\n');
fprintf('g_center_all_x = [');
fprintf(1, '%f,',g_center_all_revise(:,2));
fprintf('];');
fprintf('\n#修正的光栅倾斜角theta_p_all_revise\n');
fprintf('theta_p_all = [');
fprintf(1, '%f,',theta_p_all_revise);
fprintf('];');
fprintf('\n#修正的光栅面的长度d_g_all_revise\n');
fprintf('d_g_all = [');
fprintf(1, '%f,',d_g_all_revise);
fprintf('];');


%% 求解输入输出波导的参数和整体EDG结构的点坐标
%绘制每个波导的taper的直线函数，扫描求出taper直线和小罗兰圆函数的交点确定taper的点坐标
taper_iup_line = zeros(nz,2); % 输入波导taper上边直线
taper_idown_line = zeros(nz,2); % 输入波导taper下边直线
taper_kup_line = cell(N_out,1); % 输出波导的taper的直线
taper_kdown_line = cell(N_out,1);
taper_kup_small_point = cell(N_out,1);
taper_kdown_small_point = cell(N_out,1);
taper_kup_big_point = cell(N_out,1);
taper_kdown_big_point = cell(N_out,1);
wg_out_center = cell(N_out,1);
% 输入波导taper和波导连接处坐标
taper_iup_small_point = [-(L_f*cos(theta_i)+taper_l)*sin(theta_i)-wg_width/2*cos(theta_i),L_f-((L_f*cos(theta_i)+taper_l)*cos(theta_i)-wg_width/2*sin(theta_i))];
taper_idown_small_point = [-(L_f*cos(theta_i)+taper_l)*sin(theta_i)+wg_width/2*cos(theta_i),L_f-((L_f*cos(theta_i)+taper_l)*cos(theta_i)+wg_width/2*sin(theta_i))];
wg_in_center = [-(L_f*cos(theta_i)+taper_l+wg_line)*sin(theta_i),L_f-((L_f*cos(theta_i)+taper_l+wg_line)*cos(theta_i))];
% 输入波导taper端口的宽度坐标（并非和圆相连接）
taper_iup_big_point = [-(L_f*cos(theta_i))*sin(theta_i)-wg_io/2*cos(theta_i),L_f-((L_f*cos(theta_i))*sin(theta_i)-wg_io/2*sin(theta_i))];
taper_idown_big_point = [-(L_f*cos(theta_i))*sin(theta_i)+wg_io/2*cos(theta_i),L_f-((L_f*cos(theta_i))*sin(theta_i)+wg_io/2*sin(theta_i))];
plot(taper_iup_small_point(1,1),taper_iup_small_point(1,2),'red+')
plot(taper_idown_small_point(1,1),taper_idown_small_point(1,2),'red+')
plot(taper_iup_big_point(1,1),taper_iup_big_point(1,2),'red+')
plot(taper_idown_big_point(1,1),taper_idown_big_point(1,2),'red+')
% 输入taper直线
for i = 1:nz
    taper_iup_line(i,1) = z(i);
    taper_iup_line(i,2) = (z(i)-taper_iup_small_point(1,1))*((taper_iup_big_point(1,2)-taper_iup_small_point(1,2))/(taper_iup_big_point(1,1)-taper_iup_small_point(1,1)))+taper_iup_small_point(1,2);
    taper_idown_line(i,1) = z(i);
    taper_idown_line(i,2) = (z(i)-taper_idown_small_point(1,1))*((taper_idown_big_point(1,2)-taper_idown_small_point(1,2))/(taper_idown_big_point(1,1)-taper_idown_small_point(1,1)))+taper_idown_small_point(1,2);
end
plot(taper_iup_line(:,1),taper_iup_line(:,2),'green')
plot(taper_idown_line(:,1),taper_idown_line(:,2),'green')
%输出波导taper直线
for k = 1:N_out
    taper_kup_small_point{k}(1,:) = [-(L_f*cos(theta_k(k))+taper_l)*sin(theta_k(k))-wg_width/2*cos(theta_k(k)),L_f-((L_f*cos(theta_k(k))+taper_l)*cos(theta_k(k))-wg_width/2*sin(theta_k(k)))];
    taper_kdown_small_point{k}(1,:) = [-(L_f*cos(theta_k(k))+taper_l)*sin(theta_k(k))+wg_width/2*cos(theta_k(k)),L_f-((L_f*cos(theta_k(k))+taper_l)*cos(theta_k(k))+wg_width/2*sin(theta_k(k)))];
    taper_kup_big_point{k}(1,:) = [-(L_f*cos(theta_k(k)))*sin(theta_k(k))-wg_io/2*cos(theta_k(k)),L_f-((L_f*cos(theta_k(k)))*cos(theta_k(k))-wg_io/2*sin(theta_k(k)))];
    taper_kdown_big_point{k}(1,:) = [-(L_f*cos(theta_k(k)))*sin(theta_k(k))+wg_io/2*cos(theta_k(k)),L_f-((L_f*cos(theta_k(k)))*cos(theta_k(k))+wg_io/2*sin(theta_k(k)))];
    wg_out_center{k}(1,:) = [-(L_f*cos(theta_k(k))+taper_l+wg_line)*sin(theta_k(k)),L_f-((L_f*cos(theta_k(k))+taper_l+wg_line)*cos(theta_k(k)))];
    for i = 1:nz
        taper_kup_line{k}(i,1) = z(i);
        taper_kup_line{k}(i,2) = (z(i)-taper_kup_small_point{k}(1,1))*((taper_kup_big_point{k}(1,2)-taper_kup_small_point{k}(1,2))/(taper_kup_big_point{k}(1,1)-taper_kup_small_point{k}(1,1)))+taper_kup_small_point{k}(1,2);
        taper_kdown_line{k}(i,1) = z(i);
        taper_kdown_line{k}(i,2) = (z(i)-taper_kdown_small_point{k}(1,1))*((taper_kdown_big_point{k}(1,2)-taper_kdown_small_point{k}(1,2))/(taper_kdown_big_point{k}(1,1)-taper_kdown_small_point{k}(1,1)))+taper_kdown_small_point{k}(1,2);
    end
    plot(taper_kup_small_point{k}(1,1),taper_kup_small_point{k}(1,2),'red+')
    plot(taper_kdown_small_point{k}(1,1),taper_kdown_small_point{k}(1,2),'red+')
    plot(taper_kup_big_point{k}(1,1),taper_kup_big_point{k}(1,2),'red+')
    plot(taper_kdown_big_point{k}(1,1),taper_kdown_big_point{k}(1,2),'red+')
    plot(taper_kup_line{k}(:,1),taper_kup_line{k}(:,2),'green')
    plot(taper_kdown_line{k}(:,1),taper_kdown_line{k}(:,2),'green')
end
% 扫描小罗兰圆上的点坐标，求出taper与罗兰圆的交点坐标
small_circle_up = zeros(1,3); % EDG 小罗兰圆的
% 输入taper
for i = 1:nz
for j = 1:nz+1
%     dl_small_i_up = small_circle_up_p(i,2)-taper_iup_line(i,2);
%     dl_small_i_down = small_circle_down_p(i,2)-taper_idown_line(i,2);
    dl_small_i_up_z = small_circle(j,1)-taper_iup_line(i,1);
    dl_small_i_down_z = small_circle(j,1)-taper_idown_line(i,1);
    dl_small_i_up_x = small_circle(j,2)-taper_iup_line(i,2);
    dl_small_i_down_x = small_circle(j,2)-taper_idown_line(i,2);
    
    if(z(i)>-L_f/2 && z(i)<-L_f/4)
        if(abs(dl_small_i_up_z)<0.05 && abs(dl_small_i_up_x)<0.05)
            taper_iup_circle_point = [z(i+2),taper_iup_line(i+2,2)];
            small_circle_up(1,:) = [small_circle(j-40,1),small_circle(j-40,2),j-40];
        end
        if(abs(dl_small_i_down_z)<0.05 && abs(dl_small_i_down_x)<0.05)
            taper_idwon_circle_point = [z(i+2),taper_idown_line(i+2,2)];
        end
    end
end
end
taper_i_point = zeros(4,2);
taper_i_point(1,:) = taper_iup_small_point;
taper_i_point(2,:) = taper_idown_small_point;
taper_i_point(4,:) = taper_iup_circle_point;
taper_i_point(3,:) = taper_idwon_circle_point;
%taper_i_point(5,:) = taper_iup_small_point;
fprintf('\n#输入taper_i_point(:,1)\n');
fprintf('taper_in_z = [');
fprintf(1,'%f,',taper_i_point(:,1));
fprintf('];');
fprintf('\n#输入taper_i_point(:,2)\n');
fprintf('taper_in_x = [');
fprintf(1,'%f,',taper_i_point(:,2));
fprintf('];');
% fprintf('\n#输入wg_in_center(:,1)\n');
% fprintf(1,'%f,',wg_in_center(:,1));
% fprintf('\n#输入wg_in_center(:,2)\n');
% fprintf(1,'%f,',wg_in_center(:,2));
plot(taper_i_point(:,1),taper_i_point(:,2),'blue','LineWidth',2)
% 输出taper
taper_k_point = cell(N_out,1);
small_circle_down = zeros(1,3); % EDG 小罗兰圆的
for k = 1:N_out
for i = 1:nz
for j = 1:nz+1
%     dl_small_i_up = small_circle_up_p(i,2)-taper_iup_line(i,2);
%     dl_small_i_down = small_circle_down_p(i,2)-taper_idown_line(i,2);
    dl_small_k_up_z = small_circle(j,1)-taper_kup_line{k}(i,1);
    dl_small_k_down_z = small_circle(j,1)-taper_kdown_line{k}(i,1);
    dl_small_k_up_x = small_circle(j,2)-taper_kup_line{k}(i,2);
    dl_small_k_down_x = small_circle(j,2)-taper_kdown_line{k}(i,2);
    
    if(z(i)>-L_f/2 && z(i)<-L_f/4)
        if(abs(dl_small_k_up_z)<0.05 && abs(dl_small_k_up_x)<0.05)
            taper_kup_circle_point = [z(i+2),taper_kup_line{k}(i+2,2)];
        end
        if(abs(dl_small_k_down_z)<0.05 && abs(dl_small_k_down_x)<0.05)
            taper_kdwon_circle_point = [z(i+2),taper_kdown_line{k}(i+2,2)];
            small_circle_down(1,:) = [small_circle(j+40,1),small_circle(j+40,2),j+40];
        end
    end
end
end
taper_k_point{k}(1,1) = taper_kup_small_point{k}(1,1);
taper_k_point{k}(1,2) = taper_kup_small_point{k}(1,2);
taper_k_point{k}(2,1) = taper_kdown_small_point{k}(1,1);
taper_k_point{k}(2,2) = taper_kdown_small_point{k}(1,2);
taper_k_point{k}(4,1) = taper_kup_circle_point(1,1);
taper_k_point{k}(4,2) = taper_kup_circle_point(1,2);
taper_k_point{k}(3,1) = taper_kdwon_circle_point(1,1);
taper_k_point{k}(3,2) = taper_kdwon_circle_point(1,2);
%taper_k_point{k}(5,1) = taper_kup_small_point{k}(1,1);
%taper_k_point{k}(5,2) = taper_kup_small_point{k}(1,2);
fprintf('\n#输出taper_k_point_%f(:,1)\n',k);
fprintf(1,'taper_out_%d_z = [',k);
fprintf(1,'%f,',taper_k_point{k}(:,1));
fprintf('];');
fprintf('\n#输出taper_k_point_%f(:,2)\n',k);
fprintf(1,'taper_out_%d_x = [',k);
fprintf(1,'%f,',taper_k_point{k}(:,2));
fprintf('];');
% fprintf('\n#输出wg_out_center_%f(:,1)\n',k);
% fprintf(1,'%f,',wg_out_center{k}(:,1));
% fprintf('\n#输出wg_out_center_%f(:,2)\n',k);
% fprintf(1,'%f,',wg_out_center{k}(:,2));
plot(taper_k_point{k}(:,1),taper_k_point{k}(:,2),'blue','LineWidth',2)
end
%小罗兰圆需要绘制的部分
small_circle_edg = zeros(small_circle_down(1,3)-small_circle_up(1,3),2);
for c = small_circle_down(1,3):-1:small_circle_up(1,3)
    small_circle_edg((small_circle_down(1,3)-c+1),:)= small_circle(c,:) ;
end
% edg整体的点坐标
edg_all_point = [small_circle_edg;d_point_all];
fprintf('\n#整体(z)edg_all_point(:,1)\n');
fprintf('edg_z = [');
fprintf(1,'%f,',edg_all_point(:,1));
fprintf('];');
fprintf('\n#整体(x)edg_all_point(:,2)\n');
fprintf('edg_x = [');
fprintf(1,'%f,',edg_all_point(:,2));
fprintf('];');
% figure;
plot(small_circle_edg(:,1),small_circle_edg(:,2))
plot(edg_all_point(:,1),edg_all_point(:,2))
%% 绘图
hold on
plot(g_center(1,1),g_center(1,2),'blue*','LineWidth',2)
plot(g_center_all_revise(:,1),g_center_all_revise(:,2),'blue^','LineWidth',2)
plot(in_center(1,1),in_center(1,2),'blue*','LineWidth',2)
plot(out_center(1,1),out_center(1,2),'blue*','LineWidth',2)
plot(big_circle_p(:,1),big_circle_p(:,2),'black.')
plot(z_in_life,in_life,z_in_right,in_right)
%plot(small_circle_up_p(:,1),small_circle_up_p(:,2),'b.')
%plot(small_circle_down_p(:,1),small_circle_down_p(:,2),'b.')
plot([-L_f/2,L_f/2],[L_f/2,L_f/2],'black--')
plot([-L_f,L_f],[0,0],'black--')
plot([0,0],[0,L_f],'black--')
xlim([-L_f/2-20,L_f+20]);
ylim([-5,L_f+5]);
title('EDG-size');
fprintf('\nFinish\n');
