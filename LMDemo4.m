function LMDemo4
 %估计法向力，切向力
% 计算函数f的雅克比矩阵，是解析式
%syms x y z k fx fy fz mx my mz dz r ;
%E=(fx*z-fz*y-2*k*x-mx)^2+(fz*x-fx*z-2*k*y-my)^2+(fx*y-fy*x-2*k*(z-dz)-mz)^2+(x^2+y^2+(z-dz)^2-r^2);  %误差方程
%Jsym=jacobian(E,[x y z k])   %求雅克比矩阵，偏导数矩阵
% 
% Jsym =
%  
% [ 2*x + 2*fy*(mz + fy*x - fx*y + 2*k*(z-dz)) - 2*fz*(my - fz*x + fx*z + 2*k*y) + 4*k*(mx + fz*y - fx*z + 2*k*x),
%     2*y - 2*fx*(mz + fy*x - fx*y + 2*k*(z-dz)) + 2*fz*(mx + fz*y - fx*z + 2*k*x) + 4*k*(my - fz*x + fx*z + 2*k*y), 
%     2*z - 2*dz - 2*fx*(mx + fz*y - fx*z + 2*k*x) + 2*fx*(my - fz*x + fx*z + 2*k*y) + 4*k*(mz + fy*x - fx*y + 2*k*(z-dz)), 
%     4*x*(mx + fz*y - fx*z + 2*k*x) - 2*(2*dz - 2*z)*(mz + fy*x - fx*y - 2*k*(dz - z)) + 4*y*(my - fz*x + fx*z + 2*k*y)]
%  
%  雅克比矩阵
% [ 2*a_est + 2*fy(i)*(mz(i) + fy(i)*a_est - fx(i)*b_est + 2*d_est*(c_est-dz)) - 2*fz(i)*(my(i) - fz(i)*a_est + fx(i)*c_est + 2*d_est*b_est) + 4*d_est*(mx(i) + fz(i)*b_est - fx(i)*c_est + 2*d_est*a_est),
%     2*b_est - 2*fx(i)*(mz(i) + fy(i)*a_est - fx(i)*b_est + 2*d_est*(c_est-dz)) + 2*fz(i)*(mx(i) + fz(i)*b_est - fx(i)*c_est + 2*d_est*a_est) + 4*d_est*(my(i) - fz(i)*a_est + fx(i)*c_est + 2*d_est*b_est), 
%     2*c_est - 2*dz - 2*fx(i)*(mx(i) + fz(i)*b_est - fx(i)*c_est + 2*d_est*a_est) + 2*fx(i)*(my(i) - fz(i)*a_est + fx(i)*c_est + 2*d_est*b_est) + 4*d_est*(mz(i) + fy(i)*a_est - fx(i)*b_est + 2*d_est*(c_est-dz)), 
%     4*a_est*(mx(i) + fz(i)*b_est - fx(i)*c_est + 2*d_est*a_est)- 2*(2*dz - 2*c_est)*(mz(i) + fy(i)*a_est - fx(i)*b_est - 2*d_est*(dz - c_est)) + 4*b_est*(my(i) - fz(i)*a_est + fx(i)*c_est + 2*d_est*b_est)]
%  


% 拟合用数据。参见《数学试验》，p190，例2
global   fx fy fz mx my mz n  dz r
n=20;
fx=[0.004386168
0.004031294
0.003825956
0.003980974
0.004244776
0.004069917
0.003863509
0.003873473
0.004064788
0.003724719
0.003774387
0.003539226
0.004136775
0.004221406
0.003833911
0.003907866
0.004038007
0.004474937
0.003999941
0.003836547
];
fy= [0.007805193
0.007719117
0.007445964
0.007783121
0.008069525
0.007922116
0.007498033
0.007287029
0.007845691
0.007929229
0.007941548
0.007702964
0.008134819
0.00751478
0.007610855
0.007828148
0.007488279
0.007641008
0.0077917
0.007788647
];
fz= [0.012442049
0.013668551
0.013572009
0.014559381
0.013913077
0.010910856
0.012118932
0.014965001
0.012291541
0.014162754
0.014970316
0.012039643
0.014271995
0.013984738
0.014009763
0.012227078
0.012378971
0.012466386
0.010215771
0.014723545
];
mx= [0.061742466
0.058404379
0.060196154
0.059872173
0.061872028
0.060558859
0.058611743
0.057665724
0.060346752
0.056724686
0.059926193
0.05966923
0.062415965
0.058496088
0.059975754
0.060128696
0.058195189
0.060274411
0.063939914
0.059783969
];
my= [-0.069711529
-0.066272825
-0.06295722
-0.065182634
-0.069106884
-0.066382021
-0.062717423
-0.063188039
-0.066120192
-0.065547653
-0.062561147
-0.060187511
-0.065296039
-0.066600449
-0.065497071
-0.067707054
-0.066448651
-0.066487893
-0.065652072
-0.065122113
];
mz=[0.03050155
0.034349222
0.03162694
0.034221075
0.033109069
0.02167321
0.024908021
0.037691224
0.024675459
0.031150745
0.03658475
0.025872216
0.036092196
0.033306953
0.030374076
0.025616802
0.028371489
0.028421415
0.017899033
0.034938365
];
dz=12;r=18;
% 2. LM算法
% 初始猜测s
% 数据个数
Ndata=n;
% 参数维数
Nparams=4;
% 迭代最大次数
n_iters=30;
% LM算法的阻尼系数初值
lamda=0.01;
epulong=0.01;
% step1: 变量赋值
updateJ=1;
a_est=2;
b_est=2.5;
c_est=3;
d_est=0.32;
% step2: 迭代
for it=1:n_iters
    if updateJ==1
        % 根据当前估计值，计算雅克比矩阵
        J=zeros(Ndata,Nparams);
        for i=1:n
            J(i,1)= 2*a_est + 2*fy(i)*(mz(i) + fy(i)*a_est - fx(i)*b_est + 2*d_est*(c_est-dz)) - 2*fz(i)*(my(i) - fz(i)*a_est + fx(i)*c_est + 2*d_est*b_est) + 4*d_est*(mx(i) + fz(i)*b_est - fx(i)*c_est + 2*d_est*a_est);
            J(i,2)=  2*b_est - 2*fx(i)*(mz(i) + fy(i)*a_est - fx(i)*b_est + 2*d_est*(c_est-dz)) + 2*fz(i)*(mx(i) + fz(i)*b_est - fx(i)*c_est + 2*d_est*a_est) + 4*d_est*(my(i) - fz(i)*a_est + fx(i)*c_est + 2*d_est*b_est);
            J(i,3)=  2*c_est - 2*dz - 2*fx(i)*(mx(i) + fz(i)*b_est - fx(i)*c_est + 2*d_est*a_est) + 2*fx(i)*(my(i) - fz(i)*a_est + fx(i)*c_est + 2*d_est*b_est) + 4*d_est*(mz(i) + fy(i)*a_est - fx(i)*b_est + 2*d_est*(c_est-dz));
            J(i,4)=   4*a_est*(mx(i) + fz(i)*b_est - fx(i)*c_est + 2*d_est*a_est)- 2*(2*dz - 2*c_est)*(mz(i) + fy(i)*a_est - fx(i)*b_est - 2*d_est*(dz - c_est)) + 4*b_est*(my(i) - fz(i)*a_est + fx(i)*c_est + 2*d_est*b_est);
        end
        % 根据当前参数，得到函数值
        y_est= Error([a_est b_est c_est d_est]); %a_est*exp(-b_est*data_1);
        % 计算误差
        d=y_est;
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
    y_est_lm =Error([a_lm b_lm c_lm d_lm ]) ;
    f_lm=y_est_lm;
    e_lm=dot(f_lm,f_lm);
    % 根据误差，决定如何更新参数和阻尼系数
    if sum(chi2cdf([a_est b_est c_est d_est],4)-chi2cdf([a_lm b_lm c_lm d_lm],4))>epulong*dp'*(lamda*dp-J'*y_est_lm')
%e_lm        
        lamda=lamda/10;
        a_est=a_lm;
        b_est=b_lm;
        c_est=c_lm;
        d_est=d_lm;
        e=e_lm;
        disp(e);
        updateJ=1;
    else
        updateJ=0;
        disp(e)
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
global  fx fy fz mx my mz n dz r
x=t(1);
y=t(2);
z=t(3);
k=t(4);
E=zeros(1,n);
for i=1:n 
    E(i)=(fx(i)*z-fz(i)*y-2*k*x-mx(i))^2+(fz(i)*x-fx(i)*z-2*k*y-my(i))^2+(fx(i)*y-fy(i)*x-2*k*(z-dz)-mz(i))^2+(x^2+y^2+(z-dz)^2-r^2);  %误差方程
end
end
