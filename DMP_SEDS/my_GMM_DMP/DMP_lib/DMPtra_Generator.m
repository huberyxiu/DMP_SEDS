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
Num = Durtime/deltat+1;

Data_pos = zeros(Dim,Num);
Data_pos(:,1) = Data.y0;
Data_vel = zeros(Dim,Num);
Data_vel(:,1) = Data.vel(:,1);
Data_acc = zeros(Dim,Num);
ft = zeros(Dim,Num);
t = 0;
for i=1:Dim
    for j=2:Num
        ft(i,j-1) = f(t,para.gamma,para.BF,para.weight(i,:),para.h,para.c)*(g_design(i)-Data.y0(i));
        Data_pos(i,j) = Data_pos(i,j-1)+deltat*Data_vel(i,j-1);
        Data_acc(i,j-1) = para.alpha*(para.beta*(g_design(i)-Data_pos(i,j-1))-Data_vel(i,j-1))+ft(i,j-1);
        Data_vel(i,j) = Data_vel(i,j-1)+deltat*Data_acc(i,j-1);
        t = t+deltat;
    end
    t = 0;
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

