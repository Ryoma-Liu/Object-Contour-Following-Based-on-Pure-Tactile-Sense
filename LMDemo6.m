function LMDemo6
 %���Ʊ���Ħ��ϵ��
% ���㺯��f���ſ˱Ⱦ����ǽ���ʽ
syms a b c d v fn f y x real;
E=b+(a-b)*exp(-(v/c)^2)+d*v/fn-f/fn;%a*exp(-b*x);  %����
Jsym=jacobian(E,[a b c d])   %���ſ˱Ⱦ���ƫ��������
% [ exp(-v^2/c^2), 1 - exp(-v^2/c^2), (2*v^2*exp(-v^2/c^2)*(a - b))/c^3, v/fn]
% ��������ݡ��μ�����ѧ���顷��p190����2
global  data_1 V fn f n
n=9;
data_1=[0.25 0.5 1 1.5 2 3 4 6 8];
V=[0.25 0.5 1 1.5 2 3 4 6 8];
fn=[0.25 0.5 1 1.5 2 3 4 6 8];
f=[0.25 0.5 1 1.5 2 3 4 6 8];
obs_1=[0 0 0 0 0 0 0 0 0];
% 2. LM�㷨
% ��ʼ�²�s
a0=5; b0=0.1;  %�趨��ʼֵ

% ���ݸ���
Ndata=length(obs_1);
% ����ά��
Nparams=4;
% ����������
n_iters=10;
% LM�㷨������ϵ����ֵ
lamda=0.01;
% step1: ������ֵ
updateJ=1;
a_est=a0;
b_est=b0;
c_est=a0;
d_est=b0;
y_init =  Error([a_est b_est c_est d_est]) ;   %�����ʼ���
% step2: ����
for it=1:n_iters
    if updateJ==1
        % ���ݵ�ǰ����ֵ�������ſ˱Ⱦ���
        J=zeros(Ndata,Nparams);
        for i=1:length(data_1)
            J(i,1)=exp(-(V(i)/c_est)^2);%[exp(-b_est*data_1(i)) -a_est*data_1(i)*exp(-b_est*data_1(i))]
            J(i,2)=1-exp(-(V(i)/c_est)^2);
            J(i,3)=(a_est-b_est)*2*(V(i)^2)/(c_est^3)*exp(-(V(i)/c_est)^2);
            J(i,4)=V(i)/fn(i);
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
    dp=inv(H_lm)*(J'*d(:));
    g = J'*d(:);
    a_lm=a_est+dp(1);
    b_lm=b_est+dp(2);
    c_lm=c_est+dp(3);
    d_lm=d_est+dp(4);
    % �����µĿ��ܹ���ֵ��Ӧ��y�ͼ���в�e
    y_est_lm =Error([a_lm b_lm  c_lm  d_lm ]) ;
    d_lm=obs_1-y_est_lm;
    e_lm=dot(d_lm,d_lm);
    % ������������θ��²���������ϵ��
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
%��ʾ�Ż��Ľ��
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
    E(i)=b_est+(a_est-b_est)*exp(-(V(i)/c_est)^2)+d_est*V(i)/fn(i)-f(i)/fn(i);%a*exp(-b*x);  %����
end
end
