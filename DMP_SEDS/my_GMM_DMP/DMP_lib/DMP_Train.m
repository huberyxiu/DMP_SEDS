function [Weight,ftarget]=DMP_Train(Data, para)
% 输入DMP相关模型和数据
% 输出利用LWR方法求得的权重

% 对每一维进行DMP，采用相同的动力学参数和标准系统参数以及标准函数形式和个数

Num = length(Data.time);
Dim = Data.dimension;

% 计算ftarget
ftarget = zeros(Dim,Num);
for i=1:Dim
    for j=1:Num
        ftarget(i,j) = Data.acc(i,j)-para.alpha*(para.beta*(Data.g(i)-Data.pos(i,j))-Data.vel(i,j));
    end
end

% 计算s
s = zeros(Dim,Num);
for i=1:Dim
    for j=1:Num
        s(i,j) = x_t(Data.time(j),para.gamma)*(Data.g(i)-Data.y0(i));
    end
end
% plot(s(2,:)/(Data.g(2)-Data.y0(2)));
% 计算Tau
Tau = zeros(Dim,para.BF,Num,Num);
for i=1:Dim
    for k=1:para.BF
        for j=1:Num
            Tau(i,k,j,j) = psiF(para.h,para.c,x_t(Data.time(j),para.gamma),k);
        end
    end
end
% 验证PsiF函数形式
% for k=1:para.BF
%     fplot(@(x) exp(-para.h(k)*(x-para.c(k))^2),[0,1]);
%     hold on
% end
% 计算权重
Weight = zeros(Dim,para.BF);
% Tau_2d = zeros(Num,Num);
for i=1:Dim
    for j=1:para.BF
        Tau_2d(:,:) = Tau(i,j,:,:);
        Weight(i,j) = (s(i,:)*Tau_2d*ftarget(i,:)')/(s(i,:)*Tau_2d*s(i,:)');
    end
end
% 验证权重的计算是否正确
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