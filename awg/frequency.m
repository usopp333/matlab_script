clc
close all
clear all
%close all
% CoreNum = 1; %���õĴ���������  ��������� ����ѡ���������
% if matlabpool('size')<=0  %֮ǰû�д�
%     matlabpool('open','local',CoreNum);
% else  %֮ǰ�Ѿ���
%     disp('matlab pool already started');
% end
N_out=1:4;
N_out_max=length(N_out);
in_lambda=1.545:0.005:1.555;
nmax=length(in_lambda);
%dlambda = 0.008;  %�ŵ����
%da = 2.4;  %����Բ����˿ڿ��
%d = 2;  %����Բ�������в����˿ڿ��
%d_out = 2.4; %����Բ����˿ڿ�� 
%fsrpar = 1.8; %fsr=fsrpar*N_out*dlambda;  FSR��ѡȡ


%% ������һ���������ߺͽṹ�����Ĺ�ϵ����������Բ�������۳���С��λ�õ�����ʹ��������ߴ���ǿ�Ⱥʹ��Ÿ���
% ֻ�о�����Ƶ��
% �ŵ����

% figure;
% n = 1;
% lambda = 1.550;
% for dlambda = 0.001:0.002:0.01
%     names(n) = dlambda;
%     n = n + 1;
%     name =strcat( 'dlambda = ',num2str(dlambdal));
%     da = 2.4;  %����Բ����˿ڿ��
%     d = 2;  %����Բ�������в����˿ڿ��
%     d_out = 2.4; %����Բ����˿ڿ�� 
%     fsrpar = 1.6; %fsr=fsrpar*N_out*dlambda;  FSR��ѡȡ
%     N_out_order = 1;  %��Ҫ���õĺ�������ֵ���ͺ���ⳡ���������û�����ã��������еĺ�������ֵ��
%     [f,x_fsr_in,x_fsr,array_t_plot,U,f_out,array_out_U,f_out1,fsr,L_f,m,dL,N_a] = AWG_lambda(lambda,N_out_order,dlambda,da,d,d_out,fsrpar);
%     subplot(2,3,n-1);
%     plot(f,U,'color',[rand rand rand],'LineWidth',2);hold on ; %����Բ����ⳡ����Ҷ�任�󣨲����������ѹ����
%     xlabel('�������λ��');
%     ylabel('�ⳡ�ֲ�');
%     xlim([-40,40]);
%     title(strcat('����Բ����ⳡ����Ҷ�任�󣨲����������ѹ����',name));
% end

% FSR��ѡȡ   fsr=fsrpar*N_out*dlambda;  ���������һ����Ҫ��Ӱ�촮�ź�����ǿ��
% ������ĳ����Կ�����1.����ǿ������fsrpar����������٣����ȣ�
%                     2.������ļ������fsrpar��������������ȣ�
%                     3.���ĳ�����fsrparӰ�첻��
%                     4.Ӱ���������λ�ã�

% n = 1;
% lambda = 1.550;
% %figure;
% for fsrpar = 1:0.2:2
%     names(n) = fsrpar;
%     n = n+1;
%     name = strcat('fsrpar:',num2str(names));
%     dlambda = 0.005; %�ŵ����
%     da = 2;  %����Բ����˿ڿ��
%     d = 2.4;  %����Բ�������в����˿ڿ��
%     d_out = 2.4; %����Բ����˿ڿ�� 
%     N_out_order = 1;  %��Ҫ���õĺ�������ֵ���ͺ���ⳡ���������û�����ã��������еĺ�������ֵ��
%     [f,x_fsr_in,x_fsr,array_t_plot,U,f_out,array_out_U,f_out1,fsr,L_f,m,dL,N_a] = AWG_lambda(lambda,N_out_order,dlambda,da,d,d_out,fsrpar);
%     %subplot(2,2,1);
%     figure(1);
%     plot(fsrpar,L_f,'*');hold on;
%     title('����Բ�İ뾶');
%     %subplot(2,2,2);
%     figure(2);
%     plot(fsrpar,N_a,'+');hold on;
%     title('���в����ĸ���');
%     %subplot(2,2,3);
%     figure(3);
%     plot(x_fsr_in,abs(array_t_plot)+0.001*(n-2),'color',[rand rand rand],'LineWidth',2);hold on ; %����Բ����ⳡ����Ҷ�任�󣨲����������ѹ����
%     xlabel('�������λ��');
%     ylabel('�ⳡ�ֲ�');
%     title(strcat('����Բ����ⳡ',name));
%     %subplot(2,2,4);
%     figure(4);
%     plot(f,U,'color',[rand rand rand],'LineWidth',2);hold on ; %����Բ����ⳡ����Ҷ�任�󣨲����������ѹ����
%     xlabel('�������λ��');
%     ylabel('�ⳡ�ֲ�');
%     %xlim([-40,40]);
%     title(strcat('����Բ����ⳡ����Ҷ�任��',name));
% end
%     
%  ���в���оо�����ѡȡ
% ���������Ϊһ�������Բ�İ뾶�Ǻ����в�����оо������������оо����й�ϵ�ܽ���˼�����ʽ��
%��1��L_f=d*d_out*n_s_eff/(m_r*dlambda);���в�����������빲ͬӰ��L_f���󲿷����ģ���Ϊ׼ȷ����
%     ���߽��ۣ�ֻ�����в����仯�����в����ļ��Խխ������Ҷ�任��ǰ�ĺ�����Խ������Ҷ�任�Ժ�Ĳ��ξ�Խխ�����ǳ�������ĳ���λ��Ӱ���С
%��2��L_f=da*d_out*n_s_eff/(m_r*dlambda);�����������Ӱ��L_f�����߽��ۣ����в�����оо������ڳ���Ӱ�첻�󣬼��Խ�����ĳ���λ��Խƫ����������
%��3��L_f=2*2.4*n_s_eff/(m_r*dlambda); �̶�����Բ�İ뾶�����߽��ۣ�
%
% ���ϵĽ��۶��Ǹ��ݲ�ͬ������Բ�ļ���뾶�ó��ģ���Ϊ��������������в�����оо�������һ���ĳߴ磬������Բ�Ĺⳡ�������ϻᷢ�����ĳ��ĳ���
% ���Ա�����˿ڵļ��Ҫ��ö࣬�����ͻᵼ����ı��
% �������L_f=d_in*d_out*n_s_eff/(m_r*dlambda);��ʽ���������Բ�İ뾶�����в�����оо���û�й�ϵ�����Է�����С���в��������ļ����
% ����ⳡ�����ĳ�������Եı�С�������ʽ�޸Ĳ�������������AWG������������ף�

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
%     title('����Բ�İ뾶');
%     subplot(2,2,2);
%     plot(d,N_a,'+');hold on;
%     title('���в�������Ŀ');
%     subplot(2,2,3);
%     plot(x_fsr_in,abs(array_t_plot)+0.01*(n-2),'color',[rand rand rand],'LineWidth',2);hold on ; %����Բ����ⳡ����Ҷ�任�󣨲����������ѹ����
%     xlabel('�������λ��');
%     ylabel('�ⳡ�ֲ�');
%  
%     title(strcat('����Բ����ⳡ',name));
%     subplot(2,2,4);
%     plot(f,U+0.2*(n-2),'color',[rand rand rand],'LineWidth',2);hold on ; %����Բ����ⳡ����Ҷ�任�󣨲����������ѹ����
%     xlabel('�������λ��');
%     ylabel('�ⳡ�ֲ�');
%     xlim([-100,100]);
%     title(strcat('����Բ����ⳡ����Ҷ�任��',name));
%     
%     
%  end

%��������˿�������Բ�ϵļ��
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
%     title('����Բ�İ뾶');
%     subplot(2,2,2);
%     plot(da,N_a,'+');hold on;
%     title('���в�������Ŀ');
%     subplot(2,2,3);
%     plot(x_fsr_in,abs(array_t_plot)+0.01*(n-2),'color',[rand rand rand],'LineWidth',2);hold on ; %����Բ����ⳡ����Ҷ�任�󣨲����������ѹ����
%     xlabel('�������λ��');
%     ylabel('�ⳡ�ֲ�');
%  
%     title(strcat('����Բ����ⳡ',name));
%     subplot(2,2,4);
%     plot(f,U+0.2*(n-2),'color',[rand rand rand],'LineWidth',2);hold on ; %����Բ����ⳡ����Ҷ�任�󣨲����������ѹ����
%     xlabel('�������λ��');
%     ylabel('�ⳡ�ֲ�');
%    % xlim([-100,100]);
%     title(strcat('����Բ����ⳡ����Ҷ�任��',name));
%     
%     
%  end







% %% ����ⳡ

dlambda = 0.005;  %�ŵ����
%da = 2.4;  %����Բ����˿ڿ��
d = 2;  %����Բ�������в����˿ڿ��
d_out = 2.4; %����Բ����˿ڿ�� 
fsrpar = 1.3; %fsr=fsrpar*N_out*dlambda;  FSR��ѡȡ
name = strcat('lambda:',num2str(in_lambda));
fn=0;
figure;
for da = 1:0.1:1.5
    na = da ;
    
    fn=fn+1;
    subplot(2,3,fn);
    namea = strcat('da:',num2str(na));
    for n=1:nmax  %��������
        nameall = strcat(namea,name);
        N_out_order = 1;  %��Ҫ���õĺ�������ֵ���ͺ���ⳡ���������û�����ã��������еĺ�������ֵ��
        lambda=in_lambda(n);
        %lambda = 1.550;
        [f,x_fsr_in,x_fsr,array_t_plot,U,f_out,array_out_U,f_out1,fsr,L_f,m,dL,N_a] = AWG_lambda(lambda,N_out_order,dlambda,da,d,d_out,fsrpar);
        lf(fn) = L_f;
        naa(fn) = N_a;
        %subplot(2,2,1); 
        %subplot(2,2,2);  
        plot(f,U,'color',[rand rand rand],'LineWidth',2) ; hold on; %����Բ����ⳡ����Ҷ�任�󣨲����������ѹ����
        xlabel('�������λ��');
        ylabel('�ⳡ�ֲ�');
        xlim([-20,20]);      
    end
     title(strcat('�ܵ�����ⳡ',nameall));
end




%% �����ŵ�����
% sp_out = cell(N_out_max,1);
% sp_out_dB1 = cell(N_out_max,1);
% sp_out_dB2 = cell(N_out_max,1);
% dlambda = 0.008;  %�ŵ����
% da = 2;  %����Բ����˿ڿ��
% d = 1.8;  %����Բ�������в����˿ڿ��
% d_out = 2.4; %����Բ����˿ڿ�� 
% fsrpar = 1.8; %fsr=fsrpar*N_out*dlambda;  FSR��ѡȡ
% for i = 1:N_out_max
%     N_out_order=N_out(i);
%     for  n = 1:nmax  %parfor ����
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
%     title('����');
%     xlabel('lambda');
%     ylabel('�ⳡ�ֲ�');
%     xlim([min(in_lambda)-0.001,max(in_lambda)+0.001]);
%     figure(2);
%     plot(in_lambda,sp_out_dB1{i},'color',[rand rand rand],'LineWidth',2);hold on
%     title('����');
%     xlabel('lambda');
%     ylabel('P (dB)');
%     xlim([min(in_lambda)-0.001,max(in_lambda)+0.001]);
%     figure(3);
%     plot(in_lambda,sp_out_dB2{i},'color',[rand rand rand],'LineWidth',2);hold on
%     title('����');
%     xlabel('lambda');
%     ylabel('P (dB)');
%     xlim([min(in_lambda)-0.001,max(in_lambda)+0.001]);
% end
% 
% 
% 
% % matlabpool close %��������ر�









