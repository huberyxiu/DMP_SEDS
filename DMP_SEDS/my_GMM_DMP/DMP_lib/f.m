function r = f(t,gamma,BF,weight,h,c)
% f函数（强迫力函数）
% f = sum(psiF)weight/sum(psiF)*x
x = x_t(t,gamma);
% 计算分子num
num = 0;
for i=1:BF
    num = num + psiF(h,c,x,i)*weight(i);
end
% 计算分母den
den = 0;
for i=1:BF
    den = den + psiF(h,c,x,i);
end
r = x*num/den;
end