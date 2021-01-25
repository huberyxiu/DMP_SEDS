function Data = preprocess_DMP(demo,dt)
% 输入待处理的数据和相邻轨迹间隔时间以及截断标志
% 输出处理之后的结果，Data为结构体

d = size(demo,1); % 示教轨迹的维度
num = size(demo,2);
% 对数据进行滤波
for j=1:d
    pos(j,:) = smooth(demo(j,:),25); 
end
pos = fliplr(pos);
% 利用差分计算速度和加速度
vel = diff(pos,1,2)/dt;
vel = [vel, zeros(d,1),];
acc = diff(vel,1,2)/dt;
acc = [zeros(d,1), acc];
% 轨迹的初始点
y0 = pos(:,1);
% 轨迹的目标点
g = pos(:,end); 
% 将轨迹移至原点
% pos = pos - repmat(g,1,size(pos,2));
% 生成时间序列
time = 0:dt:(size(pos,2)-1)*dt;

% 将所有的结果保存到Data结构体中
Data.dimension = d;
Data.num = num;
Data.time = time;
Data.pos = pos;
Data.vel = vel;
Data.acc = acc;
Data.y0 = y0;
Data.g = g;


