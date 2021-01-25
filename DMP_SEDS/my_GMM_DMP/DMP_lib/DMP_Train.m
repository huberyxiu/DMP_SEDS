function [Weight,ftarget]=DMP_Train(Data, para)
% ����DMP���ģ�ͺ�����
% �������LWR������õ�Ȩ��

% ��ÿһά����DMP��������ͬ�Ķ���ѧ�����ͱ�׼ϵͳ�����Լ���׼������ʽ�͸���

Num = length(Data.time);
Dim = Data.dimension;

% ����ftarget
ftarget = zeros(Dim,Num);
for i=1:Dim
    for j=1:Num
        ftarget(i,j) = Data.acc(i,j)-para.alpha*(para.beta*(Data.g(i)-Data.pos(i,j))-Data.vel(i,j));
    end
end

% ����s
s = zeros(Dim,Num);
for i=1:Dim
    for j=1:Num
        s(i,j) = x_t(Data.time(j),para.gamma)*(Data.g(i)-Data.y0(i));
    end
end
% plot(s(2,:)/(Data.g(2)-Data.y0(2)));
% ����Tau
Tau = zeros(Dim,para.BF,Num,Num);
for i=1:Dim
    for k=1:para.BF
        for j=1:Num
            Tau(i,k,j,j) = psiF(para.h,para.c,x_t(Data.time(j),para.gamma),k);
        end
    end
end
% ��֤PsiF������ʽ
% for k=1:para.BF
%     fplot(@(x) exp(-para.h(k)*(x-para.c(k))^2),[0,1]);
%     hold on
% end
% ����Ȩ��
Weight = zeros(Dim,para.BF);
% Tau_2d = zeros(Num,Num);
for i=1:Dim
    for j=1:para.BF
        Tau_2d(:,:) = Tau(i,j,:,:);
        Weight(i,j) = (s(i,:)*Tau_2d*ftarget(i,:)')/(s(i,:)*Tau_2d*s(i,:)');
    end
end
% ��֤Ȩ�صļ����Ƿ���ȷ
Psi = zeros(para.BF,Num);
for k=1:para.BF
    for j=1:Num
        Psi(k,j) = psiF(para.h,para.c,x_t(Data.time(j),para.gamma),k);
    end
end
J = zeros(Dim,para.BF);
for i=1:Dim
    for k=1:para.BF
        for j=1:Num
            J(i,k) = J(i,k) + Psi(k,j)*(ftarget(i,j)-Weight(i,k)*s(i,j))^2;
        end
    end
end
end