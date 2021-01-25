function Data = preprocess_DMP(demo,dt)
% �������������ݺ����ڹ켣���ʱ���Լ��ضϱ�־
% �������֮��Ľ����DataΪ�ṹ��

d = size(demo,1); % ʾ�̹켣��ά��
num = size(demo,2);
% �����ݽ����˲�
for j=1:d
    pos(j,:) = smooth(demo(j,:),25); 
end
pos = fliplr(pos);
% ���ò�ּ����ٶȺͼ��ٶ�
vel = diff(pos,1,2)/dt;
vel = [vel, zeros(d,1),];
acc = diff(vel,1,2)/dt;
acc = [zeros(d,1), acc];
% �켣�ĳ�ʼ��
y0 = pos(:,1);
% �켣��Ŀ���
g = pos(:,end); 
% ���켣����ԭ��
% pos = pos - repmat(g,1,size(pos,2));
% ����ʱ������
time = 0:dt:(size(pos,2)-1)*dt;

% �����еĽ�����浽Data�ṹ����
Data.dimension = d;
Data.num = num;
Data.time = time;
Data.pos = pos;
Data.vel = vel;
Data.acc = acc;
Data.y0 = y0;
Data.g = g;


