function LMDemo1
clc;
clear ;
close all;
 
Get_fx=@(x,a,b,y) ( a*exp(-b*x) -y);%误差函数
Get_Jx=@(x,a,b) ([ exp(-b*x),-a*exp(-b*x).*x ]);%雅可比矩阵
Get_A=@(x,a,b) ( transpose(Get_Jx(x,a,b))*Get_Jx(x,a,b) );
Get_g=@(x,a,b,y) ( transpose(Get_Jx(x,a,b))*Get_fx(x,a,b,y) );
 
Get_Fx=@(x,a,b,y) ( 0.5*transpose(Get_fx(x,a,b,y))*Get_fx(x,a,b,y) );%目标函数 0.5是为了方便处理系数
Get_L=@(x,a,b,y,h) ( Get_Fx(x,a,b,y)+h'*transpose(Get_Jx(x,a,b))*Get_fx(x,a,b,y)+...
    0.5*h'*Get_A(x,a,b)*h );%近似模型函数
 
% 拟合用数据。参见《数学试验》，p190，例2
X=[0.25 0.5 1 1.5 2 3 4 6 8]';
Y=[19.21 18.15 15.36 14.10 12.89 9.32 7.45 5.24 3.01]';
 
v=2;
tao=1e-10;
e1=1e-12;
e2=1e-12;
k=0;
k_max=100; %迭代次数
 
a0=10;b0=0.5; %设定初始值
t=[a0;b0];  % 变量赋值
Ndim=2;   %参数维度
A=Get_A(X,a0,b0);
g=Get_g(X,a0,b0,Y);
 
found=(norm(g)<=e1);
mou=tao*max(diag(A));   %μk
while (~found &&(k<k_max))%
    k=k+1;
    h_lm=-inv(A+mou*eye(Ndim))*g;
    if (norm(h_lm)<=e2*(norm(t)+e2))
        found=true;
    else
        t_new=t+h_lm;  %加上步长
        Fx=Get_Fx(X,t(1),t(2),Y); %新的
        Fx_h=Get_Fx(X,t_new(1),t_new(2),Y) ;
        L_0=Get_L(X,t(1),t(2),Y,zeros(Ndim,1));
        L_h=Get_L(X,t(1),t(2), Y,h_lm);
        rho=(Fx-Fx_h)./(L_0-L_h);
        if rho>0
            t=t_new;
            A=Get_A(X,t(1),t(2));
            g=Get_g(X,t(1),t(2),Y);
            found=(norm(g)<=e1);
            mou=mou*max([0.3333,1-(2*rho-1).^3]);
            v=2;
        else
            mou=mou*v;
            v=2*v;
        end
    end
end
k
PlotFitLine(X,Y,t);
end
 
 
function PlotFitLine(X,Y,t)
FitY=@(x,a,b) ( a*exp(-b*x) );
Y_hat=FitY(X,t(1),t(2));
figure;
plot(X,Y);
hold on;
plot(X,Y_hat,'g -*');
legend('原始数据','拟合数据');
title('LM最小二乘优化拟合');
end
