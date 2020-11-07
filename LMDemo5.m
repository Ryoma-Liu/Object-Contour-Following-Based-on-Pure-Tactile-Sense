function LMDemo5
 %���Ʒ�������������
% ���㺯��f���ſ˱Ⱦ����ǽ���ʽ
%syms x y z k fx fy fz mx my mz dz r ;
%E=(fx*z-fz*y-2*k*x-mx)^2+(fz*x-fx*z-2*k*y-my)^2+(fx*y-fy*x-2*k*z-mz)^2+(x^2+y^2+(z-dz)^2-r^2);  %����
%Jsym=jacobian(E,[x y z k])   %���ſ˱Ⱦ���ƫ��������
% 
% Jsym =
%  
% [ 2*x + 2*fy*(mz + fy*x - fx*y + 2*k*z) - 2*fz*(my - fz*x + fx*z + 2*k*y) + 4*k*(mx + fz*y - fx*z + 2*k*x),
%     2*y - 2*fx*(mz + fy*x - fx*y + 2*k*z) + 2*fz*(mx + fz*y - fx*z + 2*k*x) + 4*k*(my - fz*x + fx*z + 2*k*y), 
%     2*z - 2*dz - 2*fx*(mx + fz*y - fx*z + 2*k*x) + 2*fx*(my - fz*x + fx*z + 2*k*y) + 4*k*(mz + fy*x - fx*y + 2*k*z), 
%     4*x*(mx + fz*y - fx*z + 2*k*x) + 4*y*(my - fz*x + fx*z + 2*k*y) + 4*z*(mz + fy*x - fx*y + 2*k*z)]
%  
%  �ſ˱Ⱦ���
% [ 2*a_est + 2*fy(i)*(mz(i) + fy(i)*a_est - fx(i)*b_est + 2*d_est*c_est) - 2*fz(i)*(my(i) - fz(i)*a_est + fx(i)*c_est + 2*d_est*b_est) + 4*d_est*(mx(i) + fz(i)*b_est - fx(i)*c_est + 2*d_est*a_est),
%     2*b_est - 2*fx(i)*(mz(i) + fy(i)*a_est - fx(i)*b_est + 2*d_est*c_est) + 2*fz(i)*(mx(i) + fz(i)*b_est - fx(i)*c_est + 2*d_est*a_est) + 4*d_est*(my(i) - fz(i)*a_est + fx(i)*c_est + 2*d_est*b_est), 
%     2*c_est - 2*dz - 2*fx(i)*(mx(i) + fz(i)*b_est - fx(i)*c_est + 2*d_est*a_est) + 2*fx(i)*(my(i) - fz(i)*a_est + fx(i)*c_est + 2*d_est*b_est) + 4*d_est*(mz(i) + fy(i)*a_est - fx(i)*b_est + 2*d_est*c_est), 
%     4*a_est*(mx(i) + fz(i)*b_est - fx(i)*c_est + 2*d_est*a_est) + 4*b_est*(my(i) - fz(i)*a_est + fx(i)*c_est + 2*d_est*b_est) + 4*c_est*(mz(i) + fy(i)*a_est - fx(i)*b_est + 2*d_est*c_est)]
%  

M = csvread('point3.csv');
M=M(1:465,:);
% ��������ݡ��μ�����ѧ���顷��p190����2
global  data_1 fx fy fz mx my mz n  k dz r
n=10;
for j=1:46
data_1=[0.25 0.5 1 1.5 2 3 4 6 8];
fx=M((j-1)*10+1:j*10,1);
fy= M((j-1)*10+1:j*10,2);
fz= M((j-1)*10+1:j*10,3);
mx= M((j-1)*10+1:j*10,4);
my= M((j-1)*10+1:j*10,5);
mz=M((j-1)*10+1:j*10,6);
dz=12;r=18;
obs_1=[0 0 0 0 0 0 0 0 0 0];
% 2. LM�㷨
% ��ʼ�²�s
a0=5; b0=0.1;  %�趨��ʼֵ
% ���ݸ���
Ndata=length(obs_1);
% ����ά��
Nparams=4;
% ����������
n_iters=3;
% LM�㷨������ϵ����ֵ
lamda=0.01;
epulong=0.01;
% step1: ������ֵ
updateJ=1;
a_est=a0;
b_est=b0;
c_est=a0;
d_est=b0;
y_init = Error([a_est b_est c_est d_est]) ;%a0*exp(-b0*data_1);   %�����ʼ���
% step2: ����
for it=1:n_iters
    if updateJ==1
        % ���ݵ�ǰ����ֵ�������ſ˱Ⱦ���
        J=zeros(Ndata,Nparams);
        for i=1:10
            J(i,1)= 2*a_est + 2*fy(i)*(mz(i) + fy(i)*a_est - fx(i)*b_est + 2*d_est*c_est) - 2*fz(i)*(my(i) - fz(i)*a_est + fx(i)*c_est + 2*d_est*b_est) + 4*d_est*(mx(i) + fz(i)*b_est - fx(i)*c_est + 2*d_est*a_est);
            J(i,2)= 2*b_est - 2*fx(i)*(mz(i) + fy(i)*a_est - fx(i)*b_est + 2*d_est*c_est) + 2*fz(i)*(mx(i) + fz(i)*b_est - fx(i)*c_est + 2*d_est*a_est) + 4*d_est*(my(i) - fz(i)*a_est + fx(i)*c_est + 2*d_est*b_est);
            J(i,3)= 2*c_est - 2*dz - 2*fx(i)*(mx(i) + fz(i)*b_est - fx(i)*c_est + 2*d_est*a_est) + 2*fx(i)*(my(i) - fz(i)*a_est + fx(i)*c_est + 2*d_est*b_est) + 4*d_est*(mz(i) + fy(i)*a_est - fx(i)*b_est + 2*d_est*c_est);
            J(i,4)= 4*a_est*(mx(i) + fz(i)*b_est - fx(i)*c_est + 2*d_est*a_est) + 4*b_est*(my(i) - fz(i)*a_est + fx(i)*c_est + 2*d_est*b_est) + 4*c_est*(mz(i) + fy(i)*a_est - fx(i)*b_est + 2*d_est*c_est);
        end
        % ���ݵ�ǰ�������õ�����ֵ
        y_est= Error([a_est b_est c_est d_est]); %a_est*exp(-b_est*data_1);
        % �������
        d=obs_1-y_est;
        % ���㣨�⣩��������
        H=J'*J;
        % ���ǵ�һ�ε������������
        if it==1
            e=dot(d,d);
        end
    end
    % ��������ϵ��lamda��ϵõ�H����
    H_lm=H+(lamda*eye(Nparams,Nparams)*diag(H));
    % ���㲽��dp�������ݲ��������µĿ��ܵ�\��������ֵ
    dp=-inv(H_lm)*(J'*d(:));
    g = J'*d(:);
    a_lm=a_est+dp(1);
    b_lm=b_est+dp(2);
    c_lm=c_est+dp(3);
    d_lm=d_est+dp(4);
    % �����µĿ��ܹ���ֵ��Ӧ��y�ͼ���в�e
    y_est_lm =Error([a_lm b_lm c_lm d_lm ]) ;
    f_lm=obs_1-y_est_lm;
    e_lm=dot(f_lm,f_lm);
    % ������������θ��²���������ϵ��
    if sum(chi2cdf([a_est b_est c_est d_est],4)-chi2cdf([a_lm b_lm  c_lm d_lm],4))>epulong*dp'*(lamda*dp-J'*y_est_lm')
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
        lamda=lamda*10;
    end
end
%��ʾ�Ż��Ľ��
x(j)=a_est;
y(j)=b_est;
z(j)=c_est;
kk(j)=d_est;
end
plot3(x,y,z,'o');
end

 
function E=Error(t)
global  fx fy fz mx my mz n k dz r
x=t(1);
y=t(2);
z=t(3);
k=t(4);
E=zeros(1,n);
for i=1:n 
    E(i)=(fx(i)*z-fz(i)*y-2*k*x-mx(i))^2+(fz(i)*x-fx(i)*z-2*k*y-my(i))^2+(fx(i)*y-fy(i)*x-2*k*z-mz(i))^2+(x^2+y^2+(z-dz)^2-r^2);  %����
end
end
