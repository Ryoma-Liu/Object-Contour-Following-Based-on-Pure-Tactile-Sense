function LMDemo6
 %估计表面摩擦系数
% 计算函数f的雅克比矩阵，是解析式
syms a b c d v fn f y x real;
E=b+(a-b)*exp(-(v/c)^2)+d*v/fn-f/fn;%a*exp(-b*x);  %误差方程
Jsym=jacobian(E,[a b c d])   %求雅克比矩阵，偏导数矩阵
% [ exp(-v^2/c^2), 1 - exp(-v^2/c^2), (2*v^2*exp(-v^2/c^2)*(a - b))/c^3, v/fn]
% 拟合用数据。参见《数学试验》，p190，例2
global  data_1 V fn f n
n=9;
data_1=[0.25 0.5 1 1.5 2 3 4 6 8];
V=[0.25 0.5 1 1.5 2 3 4 6 8];
fn=[0.25 0.5 1 1.5 2 3 4 6 8];
f=[0.25 0.5 1 1.5 2 3 4 6 8];
obs_1=[0 0 0 0 0 0 0 0 0];
% 2. LM算法
% 初始猜测s
a0=5; b0=0.1;  %设定初始值

% 数据个数
Ndata=length(obs_1);
% 参数维数
Nparams=4;
% 迭代最大次数
n_iters=10;
% LM算法的阻尼系数初值
lamda=0.01;
% step1: 变量赋值
updateJ=1;
a_est=a0;
b_est=b0;
c_est=a0;
d_est=b0;
y_init =  Error([a_est b_est c_est d_est]) ;   %计算初始输出
% step2: 迭代
for it=1:n_iters
    if updateJ==1
        % 根据当前估计值，计算雅克比矩阵
        J=zeros(Ndata,Nparams);
        for i=1:length(data_1)
            J(i,1)=exp(-(V(i)/c_est)^2);%[exp(-b_est*data_1(i)) -a_est*data_1(i)*exp(-b_est*data_1(i))]
            J(i,2)=1-exp(-(V(i)/c_est)^2);
            J(i,3)=(a_est-b_est)*2*(V(i)^2)/(c_est^3)*exp(-(V(i)/c_est)^2);
            J(i,4)=V(i)/fn(i);
        end
        % 根据当前参数，得到函数值
        y_est= Error([a_est b_est c_est d_est]); %a_est*exp(-b_est*data_1);
        % 计算误差
        d=obs_1-y_est;
        % 计算（拟）海塞矩阵
        H=J'*J;
        % 若是第一次迭代，计算误差
        if it==1
            e=dot(d,d);
        end
    end
    % 根据阻尼系数lamda混合得到H矩阵
    H_lm=H+(lamda*eye(Nparams,Nparams)*diag(H));
    % 计算步长dp，并根据步长计算新的可能的\参数估计值
    dp=inv(H_lm)*(J'*d(:));
    g = J'*d(:);
    a_lm=a_est+dp(1);
    b_lm=b_est+dp(2);
    c_lm=c_est+dp(3);
    d_lm=d_est+dp(4);
    % 计算新的可能估计值对应的y和计算残差e
    y_est_lm =Error([a_lm b_lm  c_lm  d_lm ]) ;
    d_lm=obs_1-y_est_lm;
    e_lm=dot(d_lm,d_lm);
    % 根据误差，决定如何更新参数和阻尼系数
    if e_lm        
        lamda=lamda/10;
        a_est=a_lm;
        b_est=b_lm;
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
global  data_1 V fn f n
a_est=t(1);
b_est=t(2);
c_est=t(3);
d_est=t(4);
E=zeros(1,4);
for i=1:n 
    E(i)=b_est+(a_est-b_est)*exp(-(V(i)/c_est)^2)+d_est*V(i)/fn(i)-f(i)/fn(i);%a*exp(-b*x);  %误差方程
end
end
