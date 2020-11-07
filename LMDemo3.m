function LMDemo3
 %���Ʊ���Ħ��ϵ�� p=[us uc v0 aita] ��Ħ��ϵ������Ħ��ϵ���� ˹���ﱴ������ϵ��������ճ��ϵ��
 %�������ݣ�Ħ����f�ͷ�����fn �ٶ����� V
% ���㺯��f���ſ˱Ⱦ����ǽ���ʽ
syms a b c d v fn f y x real;
E=b+(a-b)*exp(-(v/c)^2)+d*v/fn-f/fn;%a*exp(-b*x);  %����
Jsym=jacobian(E,[a b c d])   %���ſ˱Ⱦ���ƫ��������
% [ exp(-v^2/c^2), 1 - exp(-v^2/c^2), (2*v^2*exp(-v^2/c^2)*(a - b))/c^3, v/fn]
% ��������ݡ��μ�����ѧ���顷��p190����2
global  V fn f n
n=9;
V=[0.25 0.25 0.25 0.25 0.25 0.25 0.25 0.25 0.25];
fn=[0.25 0.25 0.25 0.25 0.25 0.25 0.25 0.25 0.25];
f=[0.25 0.25 0.25 0.25 0.25 0.25 0.25 0.25 0.25];
% 2. LM�㷨
% ��ʼ�²�s
% ���ݸ���
Ndata=9;
% ����ά��
Nparams=4;
% ����������
n_iters=7;
% LM�㷨������ϵ����ֵ
lamda=0.01;
% step1: ������ֵ
updateJ=1;
a_est=1;
b_est=5;
c_est=2;
d_est=1;
% step2: ����
for it=1:n_iters
    if updateJ==1
        % ���ݵ�ǰ����ֵ�������ſ˱Ⱦ���
        J=zeros(Ndata,Nparams);
        for i=1:Ndata
            J(i,1)=exp(-(V(i)/c_est)^2);%[exp(-b_est*data_1(i)) -a_est*data_1(i)*exp(-b_est*data_1(i))]
            J(i,2)=1-exp(-(V(i)/c_est)^2);
            J(i,3)=exp(-(V(i)^2/c_est^2))*(a_est-b_est)*2*(V(i)^2)/(c_est^3);
            J(i,4)=V(i)/fn(i);
        end
        % ���ݵ�ǰ�������õ�����ֵ
        y_est= Error([a_est b_est c_est d_est]); %a_est*exp(-b_est*data_1);
        % �������
        dd=y_est;
        % ���㣨�⣩��������
        H=J'*J;
        % ���ǵ�һ�ε������������
        if it==1
            e=dot(dd,dd);
        end
    end
    % ��������ϵ��lamda��ϵõ�H����
    H_lm=H+(lamda*eye(Nparams,Nparams)*diag(H));
    % ���㲽��dp�������ݲ��������µĿ��ܵ�\��������ֵ
    dp=inv(H_lm)*(J'*dd(:));
    g = J'*dd(:);
    a_lm=a_est+dp(1);
    b_lm=b_est+dp(2);
    c_lm=c_est+dp(3);
    d_lm=d_est+dp(4);
    % �����µĿ��ܹ���ֵ��Ӧ��y�ͼ���в�e
    y_est_lm =Error([a_lm b_lm  c_lm  d_lm ]) ;
    dd_lm=y_est_lm;
    e_lm=dot(dd_lm,dd_lm);
    % ������������θ��²���������ϵ��
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
%��ʾ�Ż��Ľ��
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
    E=b_est+(a_est-b_est)*exp(-(V(i)/c_est)^2)+d_est*V(i)/fn(i)-f(i)/fn(i);%a*exp(-b*x);  %����
end
end
