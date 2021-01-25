function Data_gnr=DMP_Generator(Data,para,g_design)
% �ı�Ŀ��λ�ã��������еĲ��������µĹ켣
% ����ԭʼ�켣��������̬ϵͳ������Ŀ��λ��
% ������ɵ��µĹ켣�㣬�µĹ켣��ĸ�����ԭʼ�켣��һ��

%  ����ѧ����Ϊ�� dy = z;
%  ����ѧ����Ϊ�� dz = alpha(belta(g-y)-z)+f
%  f��t�ĺ����� f = sum(psiF)weight/sum(psiF)*x
%  Runge-Kutta������ܸ��ӣ�����Euler������С����
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
% ������������˹��Ƴ���˼��
% ���ڲ�ȡ����켣˼�룬ԭ�����յ�����Ϊ��㣬�丽�����ٶȺ�С��������ܴ󣬹���Ҫ����һ���������ѹ����
% �������һֱ���ã�����һ���������ӣ�ԭ����vanish��
% �ʳ����ı��ʽΪ�� Frep =
% fr*vanish*(dist(x,start)^2)*(1/dist(x,nearest)-1/rou0)
% ʵ���ϣ����ֳ��������Ҫ�������fr��Ҫȡ�ϴ��ֵ������ʼʱfr��Ӧ��ȡ�ܴ�vanish��֮ǰ��Ϊ�˱����ȶ��������ɹ켣����Ҫ�ȶ�����ȡ��
% ��dist(x,start)^2��ἱ�����ӣ���Ҫ���ƣ���ȡ������fr������
% ��ˣ�fr*(1/dist(x,nearest)-1/rou0)
% ��������Ϊ���е��ʸ����
fr = Data.fr; %���Ƴ�����С�ĳ���
rou0 = Data.rou0; %�������뷶Χ
rapid = Data.rapid;
% Frep = zeros(Dim,1);
% Frep = 20*[-1;0]; % �������ߵĳ�ʼ���ɷ���
% Frep = 20*[0;-1]; % �������ߵĳ�ʼ���ɷ���
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
            ori = (Data_pos(:,j)-Data.pos(:,k))/norm((Data_pos(:,j)-Data.pos(:,k)));  % �����ķ���
        end
        
        distk = norm((Data_pos(:,j)-Data.pos(:,k)));   % ��ԭ�켣�ϵĵ�ľ���
        if distk>rou0||distk==0   % ���㵥���ɾ���Ӳ����ĳ���
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
% ��֤ft��ftarget������
% plot(deltat:deltat:Num*deltat,ft(2,:),'r');
% hold on
% plot(Data.time,Data.ftarget(2,:),'b');
% �ֱ���֤x��y�ĸ������
% plot(deltat:deltat:Num*deltat,Data_pos(2,:),'r');
% hold on
% plot(Data.time,Data.pos(2,:),'b');
end

