function [x0 , xT, Data, index] = my_preprocess_demos(demos,time,tol_cutting)
% 修改： Fred
% 对单次示教点采用多维高斯分布生成其他可能的示教点
% 输出的x0为单次示教的初始点，xT为单次示教的终点，Data为单次示教轨迹加上生成点
% index用于区分实际示教轨迹和生成点，表征一个划分
% This function preprocess raw data and put them in a format suitable for
% SEDS. The function computes the first time derivative of demonstrations,
% shift the final point of all demonstrations to the origin (this only
% simplify computation load in SEDS), and trim datas. The function can be
% called using: 
%
%          [x0 , xT, Data, index] = preprocess_demos(demos,time,tol_cutting)
%
% Inputs -----------------------------------------------------------------
%
%
%   o demos:   A variable containing all demonstrations (only
%              trajectories). The variable 'demos' should follow the
%              following format:
%              - demos{n}: d x T^n matrix representing the d dimensional
%                          trajectories. T^n is the number of datapoint in
%                          this demonstration (1 < n < N)
%
%   o time:    This variable can be provided in two ways. If time steps
%              between all demonstrations are the same, 'time' could be
%              given as a positive scalar (i.e. time = dt). If not, 'time'
%              should follow the following format:
%              - time{n}: 1 x T^n vector representing the time array of length
%                         T^n corresponding to the first demo  (1 < n < N)
%
%   o tol_cutting:  A small positive scalar that is used to trim data. It
%                   removes the redundant datapoint from the begining and
%                   the end of each demonstration that their first time
%                   derivative is less than 'tol_cutting'. Though this is
%                   not necessary for SEDS; however from practical point of
%                   view, it is very useful. There are always lots of noisy
%                   data at the begining (before the user starts the
%                   demosntration) and the end (after the user finished the
%                   demonstration) of each demosntration that are not
%                   useful.
%
% Outputs ----------------------------------------------------------------
%
%   o x0:      d x 1 array representing the mean of all demonstration's
%              initial points.
%
%   o xT:      d x 1 array representing the mean of all demonstration's
%              final point (target point).
%
%   o Data:    A 2d x N_Total matrix containing all demonstration data points.
%              Rows 1:d corresponds to trajectories and the rows d+1:2d
%              are their first time derivatives. Each column of Data stands
%              for a datapoint. All demonstrations are put next to each other 
%              along the second dimension. For example, if we have 3 demos
%              D1, D2, and D3, then the matrix Data is:
%                               Data = [[D1] [D2] [D3]]
%
%   o index:   A vector of N+1 components defining the initial index of each
%              demonstration. For example, index = [1 T1 T2 T3] indicates
%              that columns 1:T1-1 belongs to the first demonstration,
%              T1:T2-1 -> 2nd demonstration, and T2:T3-1 -> 3rd
%              demonstration.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%    Copyright (c) 2010 S. Mohammad Khansari-Zadeh, LASA Lab, EPFL,   %%%
%%%          CH-1015 Lausanne, Switzerland, http://lasa.epfl.ch         %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% The program is free for non-commercial academic use. Please contact the
% author if you are interested in using the software for commercial purposes.
% The software must not be modified or distributed without prior permission
% of the authors. Please acknowledge the authors in any academic publications
% that have made use of this code or part of it. Please use this BibTex
% reference:
% 
% S. M. Khansari-Zadeh and A. Billard, "Learning Stable Non-Linear Dynamical 
% Systems with Gaussian Mixture Models", IEEE Transaction on Robotics, 2011.
%
% To get latest upadate of the software please visit
%                          http://lasa.epfl.ch/khansari
%
% Please send your feedbacks or questions to:
%                           mohammad.khansari_at_epfl.ch

%%
%checking if a fixed time step is provided or not.
if length(time)==1
    dt = time;
end

d = size(demos{1},1); %dimensionality of demosntrations

i = 1; % 选择示教轨迹
clear tmp tmp_d

% de-noising data (not necessary)
for j=1:d
    tmp(j,:) = smooth(demos{i}(j,:),25); 
end

% computing the first time derivative
if length(time)==1
    tmp_d = diff(tmp,1,2)/dt;
else
    tmp_d = diff(tmp,1,2)./repmat(diff(time{i}),d,1);
end

% trimming demonstrations
ind = find(sqrt(sum(tmp_d.*tmp_d,1))>tol_cutting);
tmp = tmp(:,min(ind):max(ind)+1);
tmp_d = tmp_d(:,min(ind):max(ind));

% saving the initial point of each demo
% x0(:,i) = tmp(:,1);
x0 = tmp(:,1);
%saving the final point (target) of each demo
% xT(:,i) = demos{i}(:,end); 
xT = demos{i}(:,end); 

% shifting demos to the origin
tmp = tmp - repmat(xT,1,size(tmp,2));

% saving demos next to each other
Data_demo = [tmp;tmp_d zeros(d,1)];
index = size(Data_demo,2);
Data = Data_demo;

Gaussian_prepross = 1;
if Gaussian_prepross == 1
%% 高斯生成模块，基于Data_demo
% 此模块思路： 终止点的sigma为0，初始点的sigma较大，为一变量。
% 位置采用多维高斯分布，速度采用单维高斯分布，通过改变x,y方向速度分量改变合速度的方向和大小
% 多维高斯分布，速度方向的sigma由采样时间间隔和速度决定，与速度垂直方向的高斯分布各维可以独立也可以不独立
% 由同一个原示教轨迹的点生成的点采用和数据点接近的速度，单维高斯分布
    Data_gnr_smooth = [];
    for gnr_num=1:3 % 生成的轨迹条数
        Gauss_num = 1; % 只生成单点
        Sigma_max = 25;
        Sigma_delta = 0.1;
        Data_gnr = [];
        for i=1:index-1
            temp_p = Data_demo(1:2,i);
            temp_v = Data_demo(3:4,i);
            temp_oriv = Data_demo(3:4,i)/norm(Data_demo(3:4,i));
            temp_orin = [0,1;-1,0]'*temp_oriv;
            R = [temp_oriv';temp_orin'];
            sigmax = dt*norm(temp_v)/3;
            sigmay = min(Sigma_delta*(index-i),Sigma_max);
            Sigma_id = [sigmax.^2,0;0,sigmay.^2];
            temp_pos_delta = (mvnrnd(zeros(2,1),Sigma_id,Gauss_num))';
            temp_pos = temp_p+R*temp_pos_delta;
            Data_gnr = [Data_gnr temp_pos];
        end
        Data_gnr = [Data_gnr zeros(d,1)];
        % 对生成的轨迹进行相关处理
        for j=1:d
            gnr_tmp(j,:) = smooth(Data_gnr(j,:),25); 
        end
        gnr_tmp_d = diff(gnr_tmp,1,2)/dt;
        Data_gnr_smooth = [Data_gnr_smooth, [gnr_tmp;gnr_tmp_d zeros(d,1)]];
    end
    weight = 1;
    Data = repmat(Data_demo,1,weight);
    % Data = [repmat(Data_demo,1,weight) Data_gnr_smooth];
end
%% 画图调试模块
figure(1)
plot(Data(1,:),Data(2,:),'k.');
hold on
plot(Data_demo(1,:),Data_demo(2,:),'r.')

% figure(2)
% plot(Data(3,:),Data(4,:),'k.');
% hold on
% plot(Data_demo(3,:),Data_demo(4,:),'r.')
end



% xT = mean(xT,2); %we consider the mean value of all demonstraions' final point as the target
% x0 = mean(x0,2); %the mean value of all demonstrations' initial point