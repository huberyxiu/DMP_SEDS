% 定义基函数形式为高斯函数 PsiF = exp(-hi(x-ci)^2)
function r=psiF(h, c, x, i)
    r=exp(-h(i)*(x-c(i))^2); 
end