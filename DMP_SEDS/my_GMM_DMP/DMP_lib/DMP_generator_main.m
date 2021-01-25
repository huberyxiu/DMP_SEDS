% 利用DMP生成新的相似轨迹，输入原始轨迹，输出原始轨迹和类似轨迹
clear all
close all

% load ('.\my_2Drecorded_motions\Trapezoid.mat','demos'); % 导入原始示教轨迹
% load ('.\my_3Drecorded_motions\my_demons_R.mat','demos'); % 导入原始示教轨迹
% demo = demos{1}; % 取示教轨迹的第一条作为初始轨迹 demo没有时间信息，只有x坐标和y坐标位置。
% load ('.\my_3Drecorded_motions\my_demons.mat','Pos2'); % 导入原始示教轨迹
% demo = Pos2; % 取示教轨迹的第一条作为初始轨迹 demo没有时间信息，只有x坐标和y坐标位置。
load ('.\my_3Drecorded_motions\my_demons.mat','Pos2'); % 导入原始示教轨迹
demo = Pos2; % 取示教轨迹的第一条作为初始轨迹 demo没有时间信息，只有x坐标和y坐标位置。
% 数据量有点大，需要较多时间处理，每三个取一个，最后插值处理
demo_num = size(demo,2);
% demo = demo(:,3:3:demo_num);
%% 用户设置基本参数
% 单维DMP系统
% ddy = -alpha*dy+alpha*beta(g-y)+f 
% dx = -gamma*x
para.dt = 0.1;  % 相邻轨迹点的坐标时间差
para.deltat = 0.1;  % 利用Euler单步法进行复现时采用的时间步长
para.BF = 500;  % 基本函数的个数
para.alpha = 8;  % 动态系统参数
para.beta = 2;   % 动态系统参数
para.gamma = 0.2;  % 标准系统参数
para.x0 = 1;     % 标准系统初值为1，终值趋向于0
para.c = 1:-1/(para.BF):1/(para.BF);   % 高斯函数的中心点
para.h = 1.5*para.BF./para.c;      % 基本函数的Gaussian Width

%% 
dt = para.dt;
Data = preprocess_DMP(demo,dt);
[Weight,ftarget] = DMP_Train(Data,para);
para.weight = Weight;
Data.ftarget = ftarget;
g_design = Data.g;
%%
% 利用新生成的轨迹生成其他考虑斥力的轨迹
N_gr = 4; % 轨迹的条数
G_des = zeros(Data.dimension,N_gr);
% Demo_gnr = zeros(N_gr,Data.dimension,Data.num);
% S曲线 ： G_des = [0,0];
G_des(:,1) = [0;0;0];
G_des(:,2) = 1.4*G_des(:,1);
G_des(:,3) = 1.5*G_des(:,2);
G_des(:,4) = 1.6*G_des(:,3);
Demo_gnr{1} = fliplr(Data.pos);
% 人工势场法的参数(Left_demons)
Data.fr = 20;
delta_fr = 20;
Data.rou0 = 10;
delta_rou = 5;
Data.rapid = 4;
for i=1:N_gr
    g_design = Data.g + G_des(:,i);
    Data_temp = DMP_Generator(Data,para,g_design);
    Demo_gnr{i+1} = Data_temp.pos(:,1:para.dt/para.deltat:end);    % 用于对轨迹进行筛选
    Demo_gnr{i+1} = fliplr(Demo_gnr{i+1});                          % 用于将轨迹变成和demo一样
    Data.fr = Data.fr + delta_fr;
    Data.rou0 = Data.rou0 + delta_rou;
end
% 对轨迹进行滤波，滤掉波动的轨迹
para.filter = 1;
if para.filter == 1
   for i=2:length(Demo_gnr)
       for j=1:Data.dimension
          Demo_gnr{i}(j,:) = smooth(Demo_gnr{i}(j,:),25); 
       end 
   end
end

% 对轨迹进行插补

% for i = 1:N_gr+1
%     x = 3:3:demo_num;
%     y1 = Demo_gnr{i}(1,:);
%     y2 = Demo_gnr{i}(2,:);
%     y3 = Demo_gnr{i}(3,:);
%     xi = 3:demo_num;
%     yi1 = interp1(x,y1,xi,'spline');
%     yi2 = interp1(x,y2,xi,'spline');
%     yi3 = interp1(x,y3,xi,'spline');
%     Demo_gnr_interp{i}(1,:) = yi1;
%     Demo_gnr_interp{i}(2,:) = yi2;
%     Demo_gnr_interp{i}(3,:) = yi3;
% end


%%
para.options = 1;
para.dimension = Data.dimension;
para.num = Data.num;
para.Ngr_demo = N_gr;
DMP_plot(para,Demo_gnr);
%% myplot
fig_option = [];
% 1表示绘制传统DMP的生成轨迹图(令frep=0)
% 2表示绘制改进DMP的生成轨迹图(frep不为0,Demo_gnr)
if ismember(1,fig_option)
    figure('name','DMPwithoutAPF');
    for i=1:length(Demo_gnr)
        if i==1
            plot(Demo_gnr{i}(1,:),Demo_gnr{i}(2,:),'r.');
        else
            plot(Demo_gnr{i}(1,:),Demo_gnr{i}(2,:),'b');
        end
        hold on
    end
    plot(0,0,'k*','markersize',12,'linewidth',2)
    box on
    xlabel('$\xi_1 /mm$','interpreter','latex','fontsize',9);
    ylabel('$\xi_2 /mm$','interpreter','latex','fontsize',9);
    set(get(gca,'xLabel'),'Fontname','Times New Roman','FontSize',9);
    set(get(gca,'yLabel'),'Fontname', 'Times New Roman','FontSize',9);
    set(gcf,'unit','centimeters','position',[10 5 4 4]) 
    set(gca,'Position',[0.3 0.27 0.5 0.65]); 
    ax = gca;
    ax.XLim = [min(demo(1,:))-(max(demo(1,:)-min(demo(1,:))))/10 max(demo(1,:))+(max(demo(1,:)-min(demo(1,:))))/10];
    ax.YLim = [min(demo(2,:))-(max(demo(2,:)-min(demo(2,:))))/10 max(demo(2,:))+(max(demo(2,:)-min(demo(2,:))))/10];
end
if ismember(2,fig_option)
    figure('name','DMPwithAPF');
    for i=1:length(Demo_gnr)
        if i==1
            plot(Demo_gnr{i}(1,:),Demo_gnr{i}(2,:),'r.');
        else
            plot(Demo_gnr{i}(1,:),Demo_gnr{i}(2,:),'b');
        end
        hold on
    end
    plot(0,0,'k*','markersize',12,'linewidth',2)
    set(gca,'Fontname', 'Times New Roman','FontSize',9);
    xlabel('$\xi_1 /mm$','interpreter','latex','fontsize',9);
    ylabel('$\xi_2 /mm$','interpreter','latex','fontsize',9);
    set(get(gca,'xLabel'),'Fontname','Times New Roman','FontSize',9);
    set(get(gca,'yLabel'),'Fontname', 'Times New Roman','FontSize',9);
    set(gcf,'unit','centimeters','position',[10 5 4 4]) 
    set(gca,'Position',[0.3 0.27 0.6 0.65]); 
    ax = gca;
    ax.XLim = [min(Data(1,:))-(max(Data(1,:)-min(Data(1,:))))/10 max(Data(1,:))+(max(Data(1,:)-min(Data(1,:))))/10];
    ax.YLim = [min(Data(2,:))-(max(Data(2,:)-min(Data(2,:))))/10 max(Data(2,:))+(max(Data(2,:)-min(Data(2,:))))/10];
end
%%
% 对算法加权重
% Demo_gnr{6} = Demo_gnr{1};
% Demo_gnr{7} = Demo_gnr{1};
% Demo_gnr{8} = Demo_gnr{1};
%%
demos = Demo_gnr;
% 保存所有的轨迹
if isempty(regexp(path,['my_generated_motions' pathsep], 'once'))
    addpath([pwd, '/my_generated_motions']);    % add my_generated_motions dir to path
end
save('my_generated_motions/my_pickcap_demo2.mat','demos');