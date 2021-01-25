function Data_gnr=DMP_Generator(Data,para,g_design)
% 改变目标位置，利用已有的参数生成新的轨迹
% 输入原始轨迹参数，动态系统参数和目标位置
% 输出生成的新的轨迹点，新的轨迹点的个数与原始轨迹不一致

%  动力学方程为： dy = z;
%  动力学方程为： dz = alpha(belta(g-y)-z)+f
%  f是t的函数： f = sum(psiF)weight/sum(psiF)*x
%  Runge-Kutta法好像很复杂，采用Euler方法和小步长
deltat = para.deltat;
Durtime = Data.time(end);
Dim = Data.dimension;
Num = round(Durtime/deltat)+1;

Data_pos = zeros(Dim,Num);
Data_pos(:,1) = Data.y0;
Data_vel = zeros(Dim,Num);
Data_vel(:,1) = Data.vel(:,1);
Data_acc = zeros(Dim,Num);
ft = zeros(Dim,Num);
t = 0;
% 加入带斥力的人工势场法思想
% 由于采取了逆轨迹思想，原来的终点现在为起点，其附近的速度很小，斥力会很大，故需要乘以一个二次项将其压下来
% 二次项不能一直作用，加入一个消减因子，原来的vanish项
% 故斥力的表达式为： Frep =
% fr*vanish*(dist(x,start)^2)*(1/dist(x,nearest)-1/rou0)
% 实际上，发现斥力项必须要增大，因此fr需要取较大的值，但初始时fr不应该取很大，vanish项之前是为了保持稳定，但生成轨迹不需要稳定，故取消
% 而dist(x,start)^2项会急速增加，需要限制，先取消，在fr中体现
% 因此：fr*(1/dist(x,nearest)-1/rou0)
% 斥力方向为所有点的矢量和
fr = Data.fr; %控制斥力大小的常数
rou0 = Data.rou0; %斥力距离范围
rapid = Data.rapid;
% Frep = zeros(Dim,1);
% Frep = 20*[-1;0]; % 控制曲线的初始生成方向
% Frep = 20*[0;-1]; % 控制曲线的初始生成方向
Frep = 20*[0;0;1];
for j=2:Num
    for i=1:Dim
        ft(i,j-1) = f(t,para.gamma,para.BF,para.weight(i,:),para.h,para.c)*(g_design(i)-Data.y0(i)) + Frep(i);
        Data_pos(i,j) = Data_pos(i,j-1)+deltat*Data_vel(i,j-1);
        Data_acc(i,j-1) = para.alpha*(para.beta*(g_design(i)-Data_pos(i,j-1))-Data_vel(i,j-1))+ft(i,j-1);
        Data_vel(i,j) = Data_vel(i,j-1)+deltat*Data_acc(i,j-1);
    end
    vanish = x_t(t,para.gamma);
    for k=1:Data.num
        if norm((Data_pos(:,j)-Data.pos(:,k))) < 10^-4
            ori = zeros(3,1);
        else
            ori = (Data_pos(:,j)-Data.pos(:,k))/norm((Data_pos(:,j)-Data.pos(:,k)));  % 斥力的方向
        end
        
        distk = norm((Data_pos(:,j)-Data.pos(:,k)));   % 与原轨迹上的点的距离
        if distk>rou0||distk==0   % 计算单纯由距离从产生的斥力
            Ftemp = 0;
        else
            Ftemp = (1/distk-1/rou0);
        end
        Frep = Frep + ori*(fr*Ftemp*(1-vanish)^rapid);
    end
    t = t+deltat;
end
Data_gnr.pos = Data_pos;
Data_gnr.vel = Data_vel;
Data_gnr.acc = Data_acc;
% 验证ft和ftarget的量级
% plot(deltat:deltat:Num*deltat,ft(2,:),'r');
% hold on
% plot(Data.time,Data.ftarget(2,:),'b');
% 分别验证x和y的跟随情况
% plot(deltat:deltat:Num*deltat,Data_pos(2,:),'r');
% hold on
% plot(Data.time,Data.pos(2,:),'b');
end

