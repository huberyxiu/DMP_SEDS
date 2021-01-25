
% This is a matlab script illustrating how to use SEDS_lib to learn
% an arbitrary model from a set of demonstrations.
%%
% To run this demo you need to provide the variable demos composed of all
% demosntration trajectories. To get more detailed information about the
% structure of the variable 'demo', type 'doc preprocess_demos' in the
% MATLAB command window
%% User Parameters and Setting
close all
clear all

% load('DMP_lib/recorded_motions/Trapezoid','demos')   % 从原始示教轨迹中得到轨迹
 load('.\DMP_lib\my_generated_motions\my_Trapezoid','demos');     % 从DMP中生成轨迹中得到轨迹
% demos(4) = []; demos(2) = [];
% the variable 'demos' composed of 3 demosntrations. Each demonstrations is
% recorded from Tablet-PC at 50Hz. Datas are in millimeters.

% % 测试一下2D_Multi_Model：
% % 结合Khamesh,Angle,Line
% % Khamesh = load('DMP_lib/recorded_motions/Khamesh','demos');
% % Angle = load('DMP_lib/recorded_motions/Angle','demos');
% % Line = load('DMP_lib/recorded_motions/Line','demos');
% % 生成的multi_models2D
% Khamesh = load('.\DMP_lib\my_generated_motions\my_Khamesh','demos');
% Angle = load('.\DMP_lib\my_generated_motions\my_Angle','demos');
% Line = load('.\DMP_lib\my_generated_motions\my_Line','demos');
% Khamesh.demos(4) = []; Khamesh.demos(2) = [];
% Angle.demos(4) = []; Angle.demos(2) = [];
% Line.demos(4) = []; Line.demos(2) = [];
% Multi_Model_2D = {...
%     Khamesh.demos{1},Khamesh.demos{2},Khamesh.demos{3},...
%     Angle.demos{1},Angle.demos{2},Angle.demos{3},...
%     Line.demos{1},Line.demos{2},Line.demos{3}
%     };
% demos = Multi_Model_2D;

% 测试一下示教的三维轨迹
% load('.\DMP_lib\my_3Drecorded_motions\my_demons_L','demos'); % 导入示教的三维轨迹
% load('.\DMP_lib\my_generated_motions\my_3DLeft_demon2','demos'); 

% 测试一下三维轨迹生成（Multi_Model)
% Left = load('.\DMP_lib\my_3Drecorded_motions\my_demons_L','demos'); 
% Right = load('.\DMP_lib\my_3Drecorded_motions\my_demons_R','demos'); 
% Left = load('.\DMP_lib\my_generated_motions\my_3DLeft_demon2','demos'); 
% Right = load('.\DMP_lib\my_generated_motions\my_3DRight_demon1','demos');   % 网球实验
% Left = load('.\DMP_lib\my_generated_motions\my_pickcap_demo2','demos'); 
% Right = load('.\DMP_lib\my_generated_motions\my_pickcap_demo1','demos');   % 帽子实验
% Right.demos(5) = []; Right.demos(4) = [];
% 
% Multi_Model_3D = {...
%     Left.demos{1},Left.demos{2},Left.demos{3},Left.demos{4},Left.demos{5},...
%     Right.demos{1},Right.demos{2},Right.demos{3},...
%     };
% demos = Multi_Model_3D;

% Pre-processing
dt = 0.1; %The time step of the demonstrations
tol_cutting = 1; % A threshold on velocity that will be used for trimming demos

% Training parameters
K = 7; %Number of Gaussian funcitons

% A set of options that will be passed to the solver. Please type 
% 'doc preprocess_demos' in the MATLAB command window to get detailed
% information about other possible options.
options.tol_mat_bias = 10^-6; % A very small positive scalar to avoid
                              % instabilities in Gaussian kernel [default: 10^-15]
                              
options.display = 1;          % An option to control whether the algorithm
                              % displays the output of each iterations [default: true]
                              
options.tol_stopping=10^-10;  % A small positive scalar defining the stoppping
                              % tolerance for the optimization solver [default: 10^-10]

options.max_iter = 500;       % Maximum number of iteration for the solver [default: i_max=1000]

options.objective = 'mse';    % 'likelihood': use likelihood as criterion to
                              % optimize parameters of GMM
                              % 'mse': use mean square error as criterion to
                              % optimize parameters of GMM
                              % 'direction': minimize the angle between the
                              % estimations and demonstrations (the velocity part)
                              % to optimize parameters of GMM                              
                              % [default: 'mse']

%% Putting GMR and SEDS library in the MATLAB Path
if isempty(regexp(path,['SEDS_lib' pathsep], 'once'))
    addpath([pwd, '/SEDS_lib']);    % add SEDS dir to path
end
if isempty(regexp(path,['GMR_lib' pathsep], 'once'))
    addpath([pwd, '/GMR_lib']);    % add GMR dir to path
end

%% SEDS learning algorithm
tic
% [x0 , xT, Data, index] = my_preprocess_demos(demos,dt,tol_cutting); %preprocessing datas
[x0 , xT, Data, index] = preprocess_demos(demos,dt,tol_cutting); %preprocessing datas
[Priors_0, Mu_0, Sigma_0] = initialize_SEDS(Data,K); %finding an initial guess for GMM's parameter
[Priors Mu Sigma]=SEDS_Solver(Priors_0,Mu_0,Sigma_0,Data,options); %running SEDS optimization solver

%% Simulation

% A set of options that will be passed to the Simulator. Please type 
% 'doc preprocess_demos' in the MATLAB command window to get detailed
% information about each option.
opt_sim.dt = 0.1;   
opt_sim.i_max = 3000;
opt_sim.tol = 0.1;
opt_sim.plot = 0;
d = size(Data,1)/2; %dimension of data


x0_all = [];

% 复现的起点
for i=1:length(demos)
    x0_all = [x0_all demos{i}(:,1)];
end

% 设计人工指定的起点和终点
% my_start = []; % 终点位置为0，这个是相对于终点的位置
% x0_all = [x0_all [-137.112871264805;-216.723229987012;-186.977171038150]];

% x0_all = Data(1:d,index(1:end-1)); %finding initial points of all demonstrations
fn_handle = @(x) GMR(Priors,Mu,Sigma,x,1:d,d+1:2*d);
[x xd]=Simulation(x0_all,[],fn_handle,opt_sim); %running the simulator
% [x xd]=Simulation(x0_gnr,[],fn_handle,opt_sim); %running the simulator
TrainingTime = toc;
%% 计算误差 误差公式见：《Learning Stable Nonlinear Dynamical Systems with GMMs》Page 8
% 计算原始数据的pos和vel的对应关系
Datap = Data(1:d,:);
Datad = Data(d+1:2*d,:);
% 计算生成数据的pos和vel的对应关系
Gnr_d = zeros(d,size(Datap,2));
for i=1:size(Datap,2)
    Gnr_d(:,i) = fn_handle(Datap(:,i));
end
% 计算所有轨迹的误差
error_record = [];
r = 0.4;
q = 0.6;
eps = 10^-10;
% 所有轨迹
error_count = 0;
for i=1:size(Datap,2)
    if norm(Datad(:,i)) > eps && norm(Gnr_d(:,i)) > eps
       error_record = [error_record , (r*((1-dot(Datad(:,i),Gnr_d(:,i)))/(norm(Datad(:,i))*norm(Gnr_d(:,i))+eps))^2 + ...
        q*(Datad(:,i)-Gnr_d(:,i))'*(Datad(:,i)-Gnr_d(:,i))/(norm(Datad(:,i))*norm(Datad(:,i))+eps))^0.5];
    end
    error_count = error_count+1;
end
% error = sum(error_record)/error_count;
error = mean(error_record);
errormax = max(error_record);
errorvar = var(error_record);
% % 单条轨迹
% error_single = 0;
% single_count = 0;
% for i=1:Sig_Num
%     if norm(Datad(:,i)) > eps && norm(Xd(:,i)) > eps
%        error_single = error_single + (r*((1-dot(Datad(:,i),Xd(:,i)))/(norm(Datad(:,i))*norm(Xd(:,i))+eps))^2 + ...
%         q*(Datad(:,i)-Xd(:,i))'*(Datad(:,i)-Xd(:,i))/(norm(Datad(:,i))*norm(Datad(:,i))+eps))^0.5;
%     end
%     single_count = single_count+1;
% end
% error2 = error_single/single_count;

% 只比较第一条轨迹的加速度
% xd_1 = xd(:,:,1);
% Xa_1 = diff(xd_1,1,2)/dt;
% Dataa = diff(Datad,1,2)/dt;
% Dataa_1 = Dataa(:,1:Sig_Num-1);
% normDa = [];
% for i=1:size(Dataa_1,2)
%     normDa = [normDa, norm(Dataa_1(:,i))];
% end
% Damax = max(normDa);
% normXa = [];
% for i=1:size(Xa_1,2)
%     normXa = [normXa, norm(Xa_1(:,i))];
% end
% Xamax = max(normXa);
%%
% plotting the result
figure('name','Results from Simulation','position',[265   200   520   720])
sp(1)=subplot(3,1,1);
hold on; box on
% plotGMM(Mu(1:2,:), Sigma(1:2,1:2,:), [0.6 1.0 0.6], 1,[0.6 1.0 0.6]);  
plot(Data(1,:),Data(2,:),'r.')
% plot(Data(1,1:Sig_Num),Data(2,1:Sig_Num),'r.')
xlabel('$\xi_1 (mm)$','interpreter','latex','fontsize',15);
ylabel('$\xi_2 (mm)$','interpreter','latex','fontsize',15);
title('Simulation Results')

sp(2)=subplot(3,1,2);
hold on; box on
% plotGMM(Mu([1 3],:), Sigma([1 3],[1 3],:), [0.6 1.0 0.6], 1,[0.6 1.0 0.6]);
plot(Data(1,:),Data(3,:),'r.')
xlabel('$\xi_1 (mm)$','interpreter','latex','fontsize',15);
ylabel('$\dot{\xi}_1 (mm/s)$','interpreter','latex','fontsize',15);

sp(3)=subplot(3,1,3);
hold on; box on
% plotGMM(Mu([2 4],:), Sigma([2 4],[2 4],:), [0.6 1.0 0.6], 1,[0.6 1.0 0.6]);
plot(Data(2,:),Data(4,:),'r.')
xlabel('$\xi_2 (mm)$','interpreter','latex','fontsize',15);
ylabel('$\dot{\xi}_2 (mm/s)$','interpreter','latex','fontsize',15);

for i=1:size(x,3)
    plot(sp(1),x(1,:,i),x(2,:,i),'linewidth',2)
    plot(sp(2),x(1,:,i),xd(1,:,i),'linewidth',2)
    plot(sp(3),x(2,:,i),xd(2,:,i),'linewidth',2)
    plot(sp(1),x(1,1,i),x(2,1,i),'ok','markersize',5,'linewidth',5)
    plot(sp(2),x(1,1,i),xd(1,1,i),'ok','markersize',5,'linewidth',5)
    plot(sp(3),x(2,1,i),xd(2,1,i),'ok','markersize',5,'linewidth',5)
end

for i=1:3
    axis(sp(i),'tight')
    ax=get(sp(i));
    axis(sp(i),...
        [ax.XLim(1)-(ax.XLim(2)-ax.XLim(1))/10 ax.XLim(2)+(ax.XLim(2)-ax.XLim(1))/10 ...
        ax.YLim(1)-(ax.YLim(2)-ax.YLim(1))/10 ax.YLim(2)+(ax.YLim(2)-ax.YLim(1))/10]);
    plot(sp(i),0,0,'k*','markersize',15,'linewidth',3)
    if i==1
        D = axis(sp(i));
    end
end


% plotting streamlines
% figure('name','Streamlines','position',[800   90   560   320])
% plotStreamLines(Priors,Mu,Sigma,D)
% hold on
% plot(Data(1,:),Data(2,:),'r.')
% plot(0,0,'k*','markersize',15,'linewidth',3)
% xlabel('$\xi_1 (mm)$','interpreter','latex','fontsize',15);
% ylabel('$\xi_2 (mm)$','interpreter','latex','fontsize',15);
% title('Streamlines of the model')
% set(gca,'position',[0.1300    0.1444    0.7750    0.7619])
%% 保存轨迹
save('my_designed_start.mat','x');
%% myPlot
fig_option = [2];
% 1表示绘制不含示教轨迹的流线图
% 2表示绘制handwriting的小图
% 3表示绘制三维图（只包含位置）
% 4表示绘制multi-model二维图（含示教轨迹）
% 5表示绘制
% 6表示只绘制示教轨迹图
% 7表示绘制示教和复现的大图
if ismember(1,fig_option) % 1表示绘制不含示教轨迹的流线图
    figure('name','Streamlines');
    plotStreamLines(Priors,Mu,Sigma,D)
    hold on
    plot(0,0,'k*','markersize',12,'linewidth',3)
    box on
    set(gca,'Fontname', 'Times New Roman','FontSize',9);
    xlabel('$\xi_1 /mm$','interpreter','latex','fontsize',9);
    ylabel('$\xi_2 /mm$','interpreter','latex','fontsize',9);
    set(get(gca,'xLabel'),'Fontname','Times New Roman','FontSize',9);
    set(get(gca,'yLabel'),'Fontname', 'Times New Roman','FontSize',9);
    set(gcf,'unit','centimeters','position',[10 5 8 4]) 
    set(gca,'Position',[0.15 0.24 0.75 0.65]);  
end
if ismember(2,fig_option) % 2表示绘制handwriting的小图
    figure('name','Handwriting');

    plot(Data(1,1:6:end),Data(2,1:6:end),'r.','markersize',5)
    hold on
    for i=1:length(demos)
        plot(x(1,1,i),x(2,1,i),'k.','markersize',15,'linewidth',2);
        hold on
        plot(x(1,:,i),x(2,:,i),'b','linewidth',1)
        hold on
    end


    plot(0,0,'k*','markersize',8,'linewidth',2)
%     set(gca,'xtick',[],'ytick',[])
%     box on
    axis off
    ax = gca;
    ax.XLim = [min(Data(1,:))-(max(Data(1,:)-min(Data(1,:))))/10 max(Data(1,:))+(max(Data(1,:)-min(Data(1,:))))/10];
    ax.YLim = [min(Data(2,:))-(max(Data(2,:)-min(Data(2,:))))/10 max(Data(2,:))+(max(Data(2,:)-min(Data(2,:))))/10];
    set(gcf,'unit','centimeters','position',[8 5 2.5 2]) 
end
if ismember(3,fig_option) % 2表示绘制正常三维生成图
    figure('name','3D_demons_reproduce');
    plot3(Data(1,1:15:end),Data(2,1:15:end),Data(3,1:15:end),'r.','markersize',5,'linewidth',1)
    hold on
    for i=1:length(demos)
%     for i=1:1
        plot3(x(1,1,i),x(2,1,i),x(3,1,i),'k.','markersize',15,'linewidth',1.5);
        hold on
        plot3(x(1,:,i),x(2,:,i),x(3,:,i),'b','linewidth',1.5)
        hold on
    end
    plot3(0,0,0,'k*','markersize',12,'linewidth',2)
    box on
    grid on
    set(gca,'Fontname', 'Times New Roman','FontSize',9);
    xlabel('$\xi_1 /mm$','interpreter','latex','fontsize',9);
    ylabel('$\xi_2 /mm$','interpreter','latex','fontsize',9);
    zlabel('$\xi_3 /mm$','interpreter','latex','fontsize',9);
    set(get(gca,'xLabel'),'Fontname','Times New Roman','FontSize',9);
    set(get(gca,'yLabel'),'Fontname', 'Times New Roman','FontSize',9);
    set(get(gca,'zLabel'),'Fontname', 'Times New Roman','FontSize',9);
    set(gcf,'unit','centimeters','position',[10 5 8 4]) 
    set(gca,'Position',[0.18 0.24 0.75 0.65]); 
end
if ismember(4,fig_option) % 4表示绘制multi-model二维图（含示教轨迹）
    figure('name','Multi-model 2D');
    plotStreamLines(Priors,Mu,Sigma,D)
    hold on
    plot(Data(1,1:6:end),Data(2,1:6:end),'r.')
    plot(0,0,'k*','markersize',10,'linewidth',2)
    box on
    set(gca,'Fontname', 'Times New Roman','FontSize',9);
    xlabel('$\xi_1 /mm$','interpreter','latex','fontsize',9);
    ylabel('$\xi_2 /mm$','interpreter','latex','fontsize',9);
    set(get(gca,'xLabel'),'Fontname','Times New Roman','FontSize',9);
    set(get(gca,'yLabel'),'Fontname', 'Times New Roman','FontSize',9);
    set(gcf,'unit','centimeters','position',[10 5 8 4]) 
    set(gca,'Position',[0.15 0.24 0.75 0.65]); 
end
if ismember(5,fig_option) % 5表示绘制x0_gnr三维生成图
    figure('name','3D_demons_generate');
    for i=1:size(x0_gnr,2)
        plot3(x(1,1,i),x(2,1,i),x(3,1,i),'k.','markersize',20,'linewidth',2);
        hold on
        plot3(x(1,:,i),x(2,:,i),x(3,:,i),'b','linewidth',2)
        hold on
    end
    plot3(0,0,0,'k*','markersize',12,'linewidth',2)
    box on
    set(gca,'Fontname', 'Times New Roman','FontSize',9);
    xlabel('$\xi_1 /mm$','interpreter','latex','fontsize',9);
    ylabel('$\xi_2 /mm$','interpreter','latex','fontsize',9);
    zlabel('$\xi_3 /mm$','interpreter','latex','fontsize',9);
    set(get(gca,'xLabel'),'Fontname','Times New Roman','FontSize',9);
    set(get(gca,'yLabel'),'Fontname', 'Times New Roman','FontSize',9);
    set(get(gca,'zLabel'),'Fontname', 'Times New Roman','FontSize',9);
    set(gcf,'unit','centimeters','position',[10 5 8 4]) 
    set(gca,'Position',[0.18 0.24 0.75 0.65]); 
end
 if ismember(6,fig_option) % 6表示只绘制示教轨迹图
    figure('name','Demonstrations');
    plot(Data(1,:),Data(2,:),'r.')
    hold on
    plot(0,0,'k*','markersize',12,'linewidth',2)
    box on
    xlabel('$\xi_1 /mm$','interpreter','latex','fontsize',9);
    ylabel('$\xi_2 /mm$','interpreter','latex','fontsize',9);
    set(get(gca,'xLabel'),'Fontname','Times New Roman','FontSize',9);
    set(get(gca,'yLabel'),'Fontname', 'Times New Roman','FontSize',9);
    set(gcf,'unit','centimeters','position',[10 5 4 4]) 
    set(gca,'Position',[0.3 0.27 0.5 0.65]); 
    ax = gca;
    ax.XLim = [min(Data(1,:))-(max(Data(1,:)-min(Data(1,:))))/10 max(Data(1,:))+(max(Data(1,:)-min(Data(1,:))))/10];
    ax.YLim = [min(Data(2,:))-(max(Data(2,:)-min(Data(2,:))))/10 max(Data(2,:))+(max(Data(2,:)-min(Data(2,:))))/10];
 end
if ismember(7,fig_option) % 7表示绘制示教和复现的半图
    figure('name','Reproduction');
    plot(Data(1,:),Data(2,:),'r.','markersize',5)
    hold on
    for i=1:length(demos)
        plot(x(1,1,i),x(2,1,i),'k.','markersize',20,'linewidth',2);
        hold on
        plot(x(1,:,i),x(2,:,i),'b','linewidth',1)
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