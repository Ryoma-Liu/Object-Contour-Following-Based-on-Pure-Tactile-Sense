
function [sys,x0,str,ts]=funtmpl(t,x,u,flag)
%   t是采样时间
%   x是状态变量
%   u是该simulink模块的输入
%   flag是仿真过程的状态标志，用来判断当前状态是初始化、执行、输出等
%   sys是flag的函数，不同的标志执行不同的sys
%   x0是状态变量的初始值
%   str是保留参数，一般在初始化置空
%   ts是一个1*2的向量，其中ts(1)是采样周期，ts(2)是偏移量
%   本s-fun主要就是调用initialize和output函数

switch flag,
    
%Initialization，flag=0时，具体初始化函数参照35行
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
%初始化函数是必不可少的，其中的结构体形式是固定不变的
function [sys,x0,str,ts] = mdlInitializeSizes
sizes = simsizes;
%用于设置模块参数的结构体用simsizes生成
sizes.NumContStates  = 0;
%模块连续状态变量的个数
sizes.NumDiscStates  = 0;
%模块离散状态变量的个数
sizes.NumOutputs     = 6;
%模块输出变量的个数
sizes.NumInputs      = 1;
%模块输入变量的个数
sizes.DirFeedthrough = 1;
%模块是否存在直通反馈，存在为1。直通的意思是输入能直接控制输出
sizes.NumSampleTimes = 1;
%模块的采样时间个数，至少是一个

sys = simsizes(sizes);  %设置完后赋给sys输出
str = [];   %   str是保留参数，一般在初始化置空
x0  = [];   %   x0是状态变量的初始值，显然初始值置空
ts  = [-1 0];   % sample time: [period, offset]

%start mdlOutputs
%当运行到mdlOutputs，就会输出20*u。
function [sys]= mdlOutputs(t,x,u)
load data.mat
j=1;
normal=zeros(1,3);
tforce=zeros(1,3);
 [loc, normal(1:3),tforce(1:3),forcemag(1:3)]=process_values(M(j,:));
sys=[normal tforce];
% end mdlOutputs
