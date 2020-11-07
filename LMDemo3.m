function LMDemo3
 %估计表面摩擦系数 p=[us uc v0 aita] 动摩擦系数，静摩擦系数， 斯特里贝克曲线系数，表面粘性系数
 %输入数据：摩擦力f和法向力fn 速度曲线 V
% 计算函数f的雅克比矩阵，是解析式
syms a b c d v fn f y x real;
E=b+(a-b)*exp(-(v/c)^2)+d*v/fn-f/fn;%a*exp(-b*x);  %误差方程
Jsym=jacobian(E,[a b c d])   %求雅克比矩阵，偏导数矩阵
% [ exp(-v^2/c^2), 1 - exp(-v^2/c^2), (2*v^2*exp(-v^2/c^2)*(a - b))/c^3, v/fn]
% 拟合用数据。参见《数学试验》，p190，例2
global  V fn f n
n=9;
V=[0.25 0.25 0.25 0.25 0.25 0.25 0.25 0.25 0.25];
fn=[0.25 0.25 0.25 0.25 0.25 0.25 0.25 0.25 0.25];
f=[0.25 0.25 0.25 0.25 0.25 0.25 0.25 0.25 0.25];
% 2. LM算法
% 初始猜测s
% 数据个数
Ndata=9;
% 参数维数
Nparams=4;
% 迭代最大次数
n_iters=7;
% LM算法的阻尼系数初值
lamda=0.01;
% step1: 变量赋值
updateJ=1;
a_est=1;
b_est=5;
c_est=2;
d_est=1;
% step2: 迭代
for it=1:n_iters
    if updateJ==1
        % 根据当前估计值，计算雅克比矩阵
        J=zeros(Ndata,Nparams);
        for i=1:Ndata
            J(i,1)=exp(-(V(i)/c_est)^2);%[exp(-b_est*data_1(i)) -a_est*data_1(i)*exp(-b_est*data_1(i))]
            J(i,2)=1-exp(-(V(i)/c_est)^2);
            J(i,3)=exp(-(V(i)^2/c_est^2))*(a_est-b_est)*2*(V(i)^2)/(c_est^3);
            J(i,4)=V(i)/fn(i);
        end
        % 根据当前参数，得到函数值
        y_est= Error([a_est b_est c_est d_est]); %a_est*exp(-b_est*data_1);
        % 计算误差
        dd=y_est;
        % 计算（拟）海塞矩阵
        H=J'*J;
        % 若是第一次迭代，计算误差
        if it==1
            e=dot(dd,dd);
        end
    end
    % 根据阻尼系数lamda混合得到H矩阵
    H_lm=H+(lamda*eye(Nparams,Nparams)*diag(H));
    % 计算步长dp，并根据步长计算新的可能的\参数估计值
    dp=inv(H_lm)*(J'*dd(:));
    g = J'*dd(:);
    a_lm=a_est+dp(1);
    b_lm=b_est+dp(2);
    c_lm=c_est+dp(3);
    d_lm=d_est+dp(4);
    % 计算新的可能估计值对应的y和计算残差e
    y_est_lm =Error([a_lm b_lm  c_lm  d_lm ]) ;
    dd_lm=y_est_lm;
    e_lm=dot(dd_lm,dd_lm);
    % 根据误差，决定如何更新参数和阻尼系数
    if e_lm        
        lamda=lamda*10;
        a_est=a_lm;
        b_est=b_lm;
        c_est=c_lm;
        d_est=d_lm;
        e=e_lm;
        disp(e);
        updateJ=1;
    else
        updateJ=0;
        lamda=lamda*10;
    end
end
%显示优化的结果
a_est
b_est
c_est
d_est
end

 
function E=Error(t)
global   V fn f n
a_est=t(1);
b_est=t(2);
c_est=t(3);
d_est=t(4);
for i=1:n
    E=b_est+(a_est-b_est)*exp(-(V(i)/c_est)^2)+d_est*V(i)/fn(i)-f(i)/fn(i);%a*exp(-b*x);  %误差方程
end
end
