clc
close all
clear all
%close all
% CoreNum = 1; %调用的处理器个数  并行运算打开 （可选择核心数）
% if matlabpool('size')<=0  %之前没有打开
%     matlabpool('open','local',CoreNum);
% else  %之前已经打开
%     disp('matlab pool already started');
% end
N_out=1:4;
N_out_max=length(N_out);
in_lambda=1.545:0.005:1.555;
nmax=length(in_lambda);
%dlambda = 0.008;  %信道间隔
%da = 2.4;  %罗兰圆输入端口宽度
%d = 2;  %罗兰圆输入阵列波导端口宽度
%d_out = 2.4; %罗兰圆输出端口宽度 
%fsrpar = 1.8; %fsr=fsrpar*N_out*dlambda;  FSR的选取


%% 测量单一波长的谱线和结构参数的关系，调整罗兰圆最后最后汇聚场大小和位置等性质使得输出谱线带宽、强度和串扰更好
% 只研究中心频率
% 信道间隔

% figure;
% n = 1;
% lambda = 1.550;
% for dlambda = 0.001:0.002:0.01
%     names(n) = dlambda;
%     n = n + 1;
%     name =strcat( 'dlambda = ',num2str(dlambdal));
%     da = 2.4;  %罗兰圆输入端口宽度
%     d = 2;  %罗兰圆输入阵列波导端口宽度
%     d_out = 2.4; %罗兰圆输出端口宽度 
%     fsrpar = 1.6; %fsr=fsrpar*N_out*dlambda;  FSR的选取
%     N_out_order = 1;  %需要设置的函数输入值，就横向光场的输出曲线没有作用，但必须有的函数输入值；
%     [f,x_fsr_in,x_fsr,array_t_plot,U,f_out,array_out_U,f_out1,fsr,L_f,m,dL,N_a] = AWG_lambda(lambda,N_out_order,dlambda,da,d,d_out,fsrpar);
%     subplot(2,3,n-1);
%     plot(f,U,'color',[rand rand rand],'LineWidth',2);hold on ; %罗兰圆输出光场傅里叶变换后（不进行坐标的压缩）
%     xlabel('输出横向位置');
%     ylabel('光场分布');
%     xlim([-40,40]);
%     title(strcat('罗兰圆输出光场傅里叶变换后（不进行坐标的压缩）',name));
% end

% FSR的选取   fsr=fsrpar*N_out*dlambda;  这个参数有一定的要求，影响串扰和中心强度
% 从输出的场可以看出：1.中心强度随着fsrpar的增大而减少，反比；
%                     2.输出场的间隔随着fsrpar的增大而增大，正比；
%                     3.中心场宽与fsrpar影响不大；
%                     4.影响输出场的位置；

% n = 1;
% lambda = 1.550;
% %figure;
% for fsrpar = 1:0.2:2
%     names(n) = fsrpar;
%     n = n+1;
%     name = strcat('fsrpar:',num2str(names));
%     dlambda = 0.005; %信道间隔
%     da = 2;  %罗兰圆输入端口宽度
%     d = 2.4;  %罗兰圆输入阵列波导端口宽度
%     d_out = 2.4; %罗兰圆输出端口宽度 
%     N_out_order = 1;  %需要设置的函数输入值，就横向光场的输出曲线没有作用，但必须有的函数输入值；
%     [f,x_fsr_in,x_fsr,array_t_plot,U,f_out,array_out_U,f_out1,fsr,L_f,m,dL,N_a] = AWG_lambda(lambda,N_out_order,dlambda,da,d,d_out,fsrpar);
%     %subplot(2,2,1);
%     figure(1);
%     plot(fsrpar,L_f,'*');hold on;
%     title('罗兰圆的半径');
%     %subplot(2,2,2);
%     figure(2);
%     plot(fsrpar,N_a,'+');hold on;
%     title('阵列波导的个数');
%     %subplot(2,2,3);
%     figure(3);
%     plot(x_fsr_in,abs(array_t_plot)+0.001*(n-2),'color',[rand rand rand],'LineWidth',2);hold on ; %罗兰圆输出光场傅里叶变换后（不进行坐标的压缩）
%     xlabel('输出横向位置');
%     ylabel('光场分布');
%     title(strcat('罗兰圆输入光场',name));
%     %subplot(2,2,4);
%     figure(4);
%     plot(f,U,'color',[rand rand rand],'LineWidth',2);hold on ; %罗兰圆输出光场傅里叶变换后（不进行坐标的压缩）
%     xlabel('输出横向位置');
%     ylabel('光场分布');
%     %xlim([-40,40]);
%     title(strcat('罗兰圆输出光场傅里叶变换后',name));
% end
%     
%  阵列波导芯芯间隔的选取
% 这个问题因为一般的罗兰圆的半径是和阵列波导的芯芯间隔和输出输入芯芯间隔有关系总结出了几种形式：
%（1）L_f=d*d_out*n_s_eff/(m_r*dlambda);阵列波导和输出输入共同影响L_f（大部分论文，较为准确）；
%     曲线结论：只有阵列波导变化，阵列波导的间隔越窄，傅里叶变换以前的函数就越宽，傅里叶变换以后的波形就越窄，但是场宽和中心场的位置影响很小
%（2）L_f=da*d_out*n_s_eff/(m_r*dlambda);输入输出波导影响L_f；曲线结论：阵列波导的芯芯间隔对于场宽影响不大，间隔越大中心场的位置越偏离中心坐标
%（3）L_f=2*2.4*n_s_eff/(m_r*dlambda); 固定罗兰圆的半径；曲线结论：
%
% 以上的结论都是根据不同的罗兰圆的计算半径得出的；因为如果输出输入和阵列波导的芯芯间隔都有一样的尺寸，在罗兰圆的光场输出结果上会发现中心场的场宽
% 明显比输出端口的间隔要大得多，这样就会导致损耗变大；
% 利用这个L_f=d_in*d_out*n_s_eff/(m_r*dlambda);公式计算的罗兰圆的半径和阵列波导的芯芯间隔没有关系，可以发现缩小阵列波导的信心间隔，
% 输出光场的中心场宽会明显的变小，这个公式修改参数来调整整个AWG的输出更加容易；

% n=1;
% lambda = 1.550;
% figure;
% for d = 1.5:0.1:2                          
%     names(n) = d;
%     n = n+1;
%     name = strcat('d: ',num2str(names));
%     fsrpar = 1.6;
%     dlambda = 0.005;
%     da = 2.4;
%     d_out = 2.4;
%     N_out_order = 1;
%     [f,x_fsr_in,x_fsr,array_t_plot,U,f_out,array_out_U,f_out1,fsr,L_f,m,dL,N_a] = AWG_lambda(lambda,N_out_order,dlambda,da,d,d_out,fsrpar);
%     subplot(2,2,1);
%     plot(d,L_f,'*');hold on;
%     title('罗兰圆的半径');
%     subplot(2,2,2);
%     plot(d,N_a,'+');hold on;
%     title('阵列波导的数目');
%     subplot(2,2,3);
%     plot(x_fsr_in,abs(array_t_plot)+0.01*(n-2),'color',[rand rand rand],'LineWidth',2);hold on ; %罗兰圆输出光场傅里叶变换后（不进行坐标的压缩）
%     xlabel('输出横向位置');
%     ylabel('光场分布');
%  
%     title(strcat('罗兰圆输入光场',name));
%     subplot(2,2,4);
%     plot(f,U+0.2*(n-2),'color',[rand rand rand],'LineWidth',2);hold on ; %罗兰圆输出光场傅里叶变换后（不进行坐标的压缩）
%     xlabel('输出横向位置');
%     ylabel('光场分布');
%     xlim([-100,100]);
%     title(strcat('罗兰圆输出光场傅里叶变换后',name));
%     
%     
%  end

%输入输出端口在罗兰圆上的间隔
% n=1;
% lambda = 1.550;
% figure;
% for da = 2:0.1:2.4                         
%     names(n) = da;
%     n = n+1;
%     name = strcat('da: ',num2str(names));
%     fsrpar = 1.4;
%     dlambda = 0.005;
%     d = 2;
%     d_out = da;
%     N_out_order = 1;
%     [f,x_fsr_in,x_fsr,array_t_plot,U,f_out,array_out_U,f_out1,fsr,L_f,m,dL,N_a] = AWG_lambda(lambda,N_out_order,dlambda,da,d,d_out,fsrpar);
%     subplot(2,2,1);
%     plot(da,L_f,'*');hold on;
%     title('罗兰圆的半径');
%     subplot(2,2,2);
%     plot(da,N_a,'+');hold on;
%     title('阵列波导的数目');
%     subplot(2,2,3);
%     plot(x_fsr_in,abs(array_t_plot)+0.01*(n-2),'color',[rand rand rand],'LineWidth',2);hold on ; %罗兰圆输出光场傅里叶变换后（不进行坐标的压缩）
%     xlabel('输出横向位置');
%     ylabel('光场分布');
%  
%     title(strcat('罗兰圆输入光场',name));
%     subplot(2,2,4);
%     plot(f,U+0.2*(n-2),'color',[rand rand rand],'LineWidth',2);hold on ; %罗兰圆输出光场傅里叶变换后（不进行坐标的压缩）
%     xlabel('输出横向位置');
%     ylabel('光场分布');
%    % xlim([-100,100]);
%     title(strcat('罗兰圆输出光场傅里叶变换后',name));
%     
%     
%  end







% %% 横向光场

dlambda = 0.005;  %信道间隔
%da = 2.4;  %罗兰圆输入端口宽度
d = 2;  %罗兰圆输入阵列波导端口宽度
d_out = 2.4; %罗兰圆输出端口宽度 
fsrpar = 1.3; %fsr=fsrpar*N_out*dlambda;  FSR的选取
name = strcat('lambda:',num2str(in_lambda));
fn=0;
figure;
for da = 1:0.1:1.5
    na = da ;
    
    fn=fn+1;
    subplot(2,3,fn);
    namea = strcat('da:',num2str(na));
    for n=1:nmax  %并行运算
        nameall = strcat(namea,name);
        N_out_order = 1;  %需要设置的函数输入值，就横向光场的输出曲线没有作用，但必须有的函数输入值；
        lambda=in_lambda(n);
        %lambda = 1.550;
        [f,x_fsr_in,x_fsr,array_t_plot,U,f_out,array_out_U,f_out1,fsr,L_f,m,dL,N_a] = AWG_lambda(lambda,N_out_order,dlambda,da,d,d_out,fsrpar);
        lf(fn) = L_f;
        naa(fn) = N_a;
        %subplot(2,2,1); 
        %subplot(2,2,2);  
        plot(f,U,'color',[rand rand rand],'LineWidth',2) ; hold on; %罗兰圆输出光场傅里叶变换后（不进行坐标的压缩）
        xlabel('输出横向位置');
        ylabel('光场分布');
        xlim([-20,20]);      
    end
     title(strcat('总的输出光场',nameall));
end




%% 各个信道谱线
% sp_out = cell(N_out_max,1);
% sp_out_dB1 = cell(N_out_max,1);
% sp_out_dB2 = cell(N_out_max,1);
% dlambda = 0.008;  %信道间隔
% da = 2;  %罗兰圆输入端口宽度
% d = 1.8;  %罗兰圆输入阵列波导端口宽度
% d_out = 2.4; %罗兰圆输出端口宽度 
% fsrpar = 1.8; %fsr=fsrpar*N_out*dlambda;  FSR的选取
% for i = 1:N_out_max
%     N_out_order=N_out(i);
%     for  n = 1:nmax  %parfor 并行
%         lambda=in_lambda(n);
%         [f,x_fsr_in,x_fsr,U,f_out,array_out_U,f_out1,fsr,L_f,m,dL,N_a] = AWG_lambda(lambda,N_out_order,dlambda,da,d,d_out,fsrpar);
%         sp_in(n) = sum(U); 
%         sp(n) = sum(f_out1);
%     end
%     sp_out{i} = sp;
%     sp_out_dB1{i} = 10*log(sp/max(sp));
%     sp_out_dB2{i} = 10*log(sp/sum(sp));
%      
% end
% for i = 1:N_out_max
%     figure(1);
%     plot(in_lambda,sp_out{i},'color',[rand rand rand],'LineWidth',2);hold on
%     title('光谱');
%     xlabel('lambda');
%     ylabel('光场分布');
%     xlim([min(in_lambda)-0.001,max(in_lambda)+0.001]);
%     figure(2);
%     plot(in_lambda,sp_out_dB1{i},'color',[rand rand rand],'LineWidth',2);hold on
%     title('光谱');
%     xlabel('lambda');
%     ylabel('P (dB)');
%     xlim([min(in_lambda)-0.001,max(in_lambda)+0.001]);
%     figure(3);
%     plot(in_lambda,sp_out_dB2{i},'color',[rand rand rand],'LineWidth',2);hold on
%     title('光谱');
%     xlabel('lambda');
%     ylabel('P (dB)');
%     xlim([min(in_lambda)-0.001,max(in_lambda)+0.001]);
% end
% 
% 
% 
% % matlabpool close %并行运算关闭









