
function [sys,x0,str,ts]=funtmpl(t,x,u,flag)
%   t�ǲ���ʱ��
%   x��״̬����
%   u�Ǹ�simulinkģ�������
%   flag�Ƿ�����̵�״̬��־�������жϵ�ǰ״̬�ǳ�ʼ����ִ�С������
%   sys��flag�ĺ�������ͬ�ı�־ִ�в�ͬ��sys
%   x0��״̬�����ĳ�ʼֵ
%   str�Ǳ���������һ���ڳ�ʼ���ÿ�
%   ts��һ��1*2������������ts(1)�ǲ������ڣ�ts(2)��ƫ����
%   ��s-fun��Ҫ���ǵ���initialize��output����

switch flag,
    
%Initialization��flag=0ʱ�������ʼ����������35��
  case 0         
    [sys,x0,str,ts] = mdlInitializeSizes;

%Calculate  outputs
  case 3
    sys=mdlOutputs(t,x,u);

%Unused flags
  case {1,4,2,9}
    sys = []; % do nothing

  otherwise
         error('Simulink:blocks:unhandledFlag', num2str(flag));
end
%
%=============================================================================
% mdlInitializeSizes
% Return the sizes, initial conditions, and sample times for the S-function.
%=============================================================================
%��ʼ�������Ǳز����ٵģ����еĽṹ����ʽ�ǹ̶������
function [sys,x0,str,ts] = mdlInitializeSizes
sizes = simsizes;
%��������ģ������Ľṹ����simsizes����
sizes.NumContStates  = 0;
%ģ������״̬�����ĸ���
sizes.NumDiscStates  = 0;
%ģ����ɢ״̬�����ĸ���
sizes.NumOutputs     = 3;
%ģ����������ĸ���
sizes.NumInputs      = 7;
%ģ����������ĸ���
sizes.DirFeedthrough = 1;
%ģ���Ƿ����ֱͨ����������Ϊ1��ֱͨ����˼��������ֱ�ӿ������
sizes.NumSampleTimes = 1;
%ģ��Ĳ���ʱ�������������һ��

sys = simsizes(sizes);  %������󸳸�sys���
str = [];   %   str�Ǳ���������һ���ڳ�ʼ���ÿ�
x0  = [];   %   x0��״̬�����ĳ�ʼֵ����Ȼ��ʼֵ�ÿ�
ts  = [-1 0];   % sample time: [period, offset]

%start mdlOutputs
%�����е�mdlOutputs���ͻ����20*u��
function [sys]= mdlOutputs(t,x,u)
X=[u(1) u(2) u(3)];
Xd=[u(4) u(5) u(6)];
Xe=X-Xd;
Ks=u(7);
F=Ks*Xe;
sys=[F(1) F(2) F(3)];
% end mdlOutputs
