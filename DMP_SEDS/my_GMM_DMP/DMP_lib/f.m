function r = f(t,gamma,BF,weight,h,c)
% f������ǿ����������
% f = sum(psiF)weight/sum(psiF)*x
x = x_t(t,gamma);
% �������num
num = 0;
for i=1:BF
    num = num + psiF(h,c,x,i)*weight(i);
end
% �����ĸden
den = 0;
for i=1:BF
    den = den + psiF(h,c,x,i);
end
r = x*num/den;
end