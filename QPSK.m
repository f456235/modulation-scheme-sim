clear all;                  % 清除所有变量
close all;                  % 关闭所有窗口
clc;                        % 清屏
%% 基本参数
M=12;                       % 产生码元数    
L=100;                      % 每码元复制L次,每个码元采样次数
Ts=0.001;                   % 每个码元的宽度,即码元的持续时间
Rb=1/Ts;                    % 码元速率1K
dt=Ts/L;                    % 采样间隔
TotalT=M*Ts;                % 总时间
t=0:dt:TotalT-dt;           % 时间
TotalT2=(M/2)*Ts;           % 总时间2
t2=0:dt:TotalT2-dt;         % 时间2
Fs=1/dt;                    % 采样间隔的倒数即采样频率

%% 产生双极性波形
wave=[1,0,1,0,1,1,0,1,0,0,1,1];      % 产生二进制随机码,M为码元个数
wave=2*wave-1;              % 单极性变双极性
fz=ones(1,L);               % 定义复制的次数L,L为每码元的采样点数
x1=wave(fz,:);              % 将原来wave的第一行复制L次，称为L*M的矩阵
jidai=reshape(x1,1,L*M);    % 产生双极性不归零矩形脉冲波形，将刚得到的L*M矩阵，按列重新排列形成1*(L*M)的矩阵

%% I、Q路码元
% I路码元是基带码元的奇数位置码元，Q路码元是基带码元的偶数位置码元
I=[]; Q=[];
for i=1:M
    if mod(i, 2)~=0
        I = [I, wave(i)];
    else
        Q = [Q, wave(i)];
    end
end
x2 = I(fz,:);               % 将原来I的第一行复制L次，称为L*(M/2)的矩阵
I_lu = reshape(x2,1,L*(M/2));% 将刚得到的L*(M/2)矩阵，按列重新排列形成1*(L*(M/2))的矩阵
x3 = Q(fz,:);               % 将原来Q的第一行复制L次，称为L*(M/2)的矩阵
Q_lu = reshape(x3,1,L*(M/2));% 将刚得到的L*(M/2)矩阵，按列重新排列形成1*(L*(M/2))的矩阵
figure(1);                  % 绘制第1幅图
subplot(311);               % 窗口分割成3*1的，当前是第1个子图 
plot(t,jidai,'LineWidth',2);% 绘制基带码元波形，线宽为2
title('基带信号波形');      % 标题
xlabel('时间/s');           % x轴标签
ylabel('幅度');             % y轴标签
axis([0,TotalT,-1.1,1.1])   % 坐标范围限制

subplot(312);               % 窗口分割成3*1的，当前是第2个子图 
plot(t2,I_lu,'LineWidth',2);% 绘制基带码元波形，线宽为2
title('I路信号波形');       % 标题
xlabel('时间/s');           % x轴标签
ylabel('幅度');             % y轴标签
axis([0,TotalT2,-1.1,1.1])  % 坐标范围限制

subplot(313);               % 窗口分割成3*1的，当前是第3个子图 
plot(t2,Q_lu,'LineWidth',2);% 绘制基带码元波形，线宽为2
title('Q路信号波形');       % 标题
xlabel('时间/s');           % x轴标签
ylabel('幅度');             % y轴标签
axis([0,TotalT2,-1.1,1.1])  % 坐标范围限制
%% QPSK调制
fc=2000;                    % 载波频率2kHz       
zb1=cos(2*pi*fc*t2);        % 载波1
psk1=I_lu.*zb1;             % PSK1的调制 
zb2=sin(2*pi*fc*t2);        % 载波2
psk2=Q_lu.*zb2;             % PSK2的调制 
qpsk=psk1+psk2;             % QPSK的实现

figure(2);                  % 绘制第2幅图
subplot(211);               % 窗口分割成2*1的，当前是第1个子图 
plot(t2,qpsk,'LineWidth',2);% 绘制基带码元波形，线宽为2
title('QPSK信号波形');      % 标题
xlabel('时间/s');           % x轴标签
ylabel('幅度');             % y轴标签
axis([0,TotalT2,-2.1,2.1])  % 坐标范围限制

%% 信号经过高斯白噪声信道
tz=awgn(qpsk,200);           % 信号qpsk中加入白噪声，信噪比为SNR=20dB
subplot(212);               % 窗口分割成2*1的，当前是第2个子图 
plot(t2,tz,'LineWidth',2);  % 绘制2ASK信号加入白噪声的波形
axis([0,TotalT2,-2.5,2.5]); % 坐标范围设置
title('通过高斯白噪声信道后的信号');% 标题
xlabel('时间/s');           % x轴标签
ylabel('幅度');             % y轴标签
%% 解调部分
figure(3);
tz1=tz.*zb1;                % 相干解调，乘以相干载波
subplot(211)                % 窗口分割成2*1的，当前是第1个子图 
plot(t2,tz1,'LineWidth',2)  % 绘制I路乘以相干载波后的信号
axis([0,TotalT2,-2.5,2.5]); % 设置坐标范围
title("I路乘以相干载波后的信号")% 标题
xlabel('时间/s');           % x轴标签
ylabel('幅度');             % y轴标签

tz2=tz.*zb2;                % 相干解调，乘以相干载波
subplot(212)                % 窗口分割成2*1的，当前是第2个子图 
plot(t2,tz2,'LineWidth',2)  % 绘制Q路乘以相干载波后的信号
axis([0,TotalT2,-2.5,2.5]); % 设置坐标范围
title("Q路乘以相干载波后的信号")% 标题
xlabel('时间/s');           % x轴标签
ylabel('幅度');             % y轴标签
%% 加噪信号经过滤波器
% 低通滤波器设计
fp=2*Rb;                    % 低通滤波器截止频率，乘以2是因为下面要将模拟频率转换成数字频率wp=Rb/(Fs/2)
b=fir1(30, fp/Fs, boxcar(31));% 生成fir滤波器系统函数中分子多项式的系数
% fir1函数三个参数分别是阶数，数字截止频率，滤波器类型
% 这里是生成了30阶(31个抽头系数)的矩形窗滤波器
[h,w]=freqz(b, 1,512);      % 生成fir滤波器的频率响应
% freqz函数的三个参数分别是滤波器系统函数的分子多项式的系数，分母多项式的系数(fir滤波器分母系数为1)和采样点数(默认)512
lvbo1=fftfilt(b,tz1);       % 对信号进行滤波，tz1是等待滤波的信号，b是fir滤波器的系统函数的分子多项式系数
lvbo2=fftfilt(b,tz2);       % 对信号进行滤波，tz2是等待滤波的信号，b是fir滤波器的系统函数的分子多项式系数
figure(4);                  % 绘制第4幅图  
subplot(311);               % 窗口分割成3*1的，当前是第1个子图 
plot(w/pi*Fs/2,20*log(abs(h)),'LineWidth',2); % 绘制滤波器的幅频响应
title('低通滤波器的频谱');  % 标题
xlabel('频率/Hz');          % x轴标签
ylabel('幅度/dB');          % y轴标签

subplot(312)                % 窗口分割成3*1的，当前是第2个子图 
plot(t2,lvbo1,'LineWidth',2); % 绘制I路经过低通滤波器后的信号
axis([0,TotalT2,-2.1,2.1]);  % 设置坐标范围
title("I路经过低通滤波器后的信号");% 标题
xlabel('时间/s');           % x轴标签
ylabel('幅度');             % y轴标签

subplot(313)                % 窗口分割成3*1的，当前是第3个子图 
plot(t2,lvbo2,'LineWidth',2); % 绘制Q路经过低通滤波器后的信号
axis([0,TotalT2,-2.1,2.1]);  % 设置坐标范围
title("Q路经过低通滤波器后的信号");% 标题
xlabel('时间/s');           % x轴标签
ylabel('幅度');             % y轴标签

%% 抽样判决
figure(5);
k=0;                        % 设置抽样限值
pdst1=1*(lvbo1>0);          % 滤波后的向量的每个元素和0进行比较，大于0为1，否则为0
pdst2=1*(lvbo2>0);          % 滤波后的向量的每个元素和0进行比较，大于0为1，否则为0
subplot(311)                % 窗口分割成3*1的，当前是第1个子图 
plot(t2,pdst1,'LineWidth',2)% 画出经过抽样判决后的信号
axis([0,TotalT2,-0.1,1.1]); % 设置坐标范用
title("I路经过抽样判决后的信号")% 标题
xlabel('时间/s');           % x轴标签
ylabel('幅度');             % y轴标签

subplot(312)                % 窗口分割成3*1的，当前是第2个子图 
plot(t2,pdst2,'LineWidth',2)% 画出经过抽样判决后的信号
axis([0,TotalT2,-0.1,1.1]); % 设置坐标范用
title("Q路经过抽样判决后的信号")% 标题
xlabel('时间/s');           % x轴标签
ylabel('幅度');             % y轴标签
%% I、Q路合并
I_zong = [];
Q_zong = [];
% 取码元的中间位置上的值进行判决
for j=L/2:L:(L*M/2)
    if pdst1(j)>0
        I_zong=[I_zong,1];
    else
        I_zong=[I_zong,-1];
    end
end
% 取码元的中间位置上的值进行判决
for k=L/2:L:(L*M/2)
    if pdst2(k)>0
        Q_zong=[Q_zong,1];
    else
        Q_zong=[Q_zong,-1];
    end
end
code = [];
% 将I路码元为最终输出的奇数位置码元，将Q路码元为最终输出的偶数位置码元
for n=1:M
    if mod(n, 2)~=0
        code = [code, I_zong((n+1)/2)];
    else
        code = [code, Q_zong(n/2)];
    end
end

x4=code(fz,:);             % 将原来code的第一行复制L次，称为L*M的矩阵
dout=reshape(x4,1,L*M);    % 产生单极性不归零矩形脉冲波形，将刚得到的L*M矩阵，按列重新排列形成1*(L*M)的矩阵

subplot(313);              % 窗口分割成3*1的，当前是第3个子图 
plot(t,dout,'LineWidth',2);% 绘制基带码元波形，线宽为2
title('I、Q路合并信号波形'); % 标题
xlabel('时间/s');          % x轴标签
ylabel('幅度');            % y轴标签
axis([0,TotalT,-1.1,1.1])  % 坐标范围限制
%% 绘制频谱
%% 信源频谱
figure(6)                   % 绘制第6幅图
T=t(end);                   % 时间
df=1/T;                     % 频谱分辨率
N=length(jidai);            % 采样长度
f=(-N/2:N/2-1)*df;          % 频率范围
mf=fftshift(abs(fft(jidai)));%对信源信号采用快速傅里叶变换并移到矩阵中心
subplot(211);               % 窗口分割成2*1的，当前是第1个子图 
plot(f,mf,'LineWidth',2);   % 绘制信源频谱波形
title("基带信号频谱");      % 标题
xlabel('频率/Hz');          % x轴标签
ylabel('幅度');             % y轴标签
axis([-20000,20000,-inf,inf])% 坐标范围限制

%% QPSK信号频谱
T2=t2(end);                 % 时间2
df2=1/T2;                   % 频谱分辨率2
N2=length(qpsk);            % 采样长度2
f2=(-N2/2:N2/2-1)*df2;      % 频率范围
sf=fftshift(abs(fft(qpsk)));% 对QPSK信号采用快速傅里叶变换并将0-fs频谱移动到-fs/2-fs/2
subplot(212)                % 窗口分割成2*1的，当前是第2个子图 
plot(f2,sf,'LineWidth',2)   % 绘制QPSK调制信号频谱
title("QPSK信号频谱")       % 标题
xlabel('频率/Hz');          % x轴标签
ylabel('幅度');             % y轴标签
axis([-20000,20000,-inf,inf])   % 坐标范围限制

%% I路乘以相干载波后的频谱
mmf=fftshift(abs(fft(tz1)));% 对相干载波信号采用快速傅里叶变换并移到矩阵中心
figure(7)                   % 绘制第7幅图
subplot(211);               % 窗口分割成2*1的，当前是第1个子图 
plot(f2,mmf,'LineWidth',2)  % 画出乘以相干载波后的频谱
title("I路乘以相干载波后的频谱")
xlabel('频率/Hz');          % x轴标签
ylabel('幅度');             % y轴标签
axis([-20000,20000,-inf,inf])% 坐标范围限制

%% Q路乘以相干载波后的频谱
mmf2=fftshift(abs(fft(tz2)));% 对相干载波信号采用快速傅里叶变换并移到矩阵中心
subplot(212);               % 窗口分割成2*1的，当前是第2个子图 
plot(f2,mmf2,'LineWidth',2) % 画出乘以相干载波后的频谱
title("Q路乘以相干载波后的频谱")
xlabel('频率/Hz');          % x轴标签
ylabel('幅度');             % y轴标签
axis([-20000,20000,-inf,inf])% 坐标范围限制

%% 经过低通滤波后的频谱
figure(8);
dmf=fftshift(abs(fft(lvbo1)));%对低通滤波信号采用快速傅里叶变换并移到矩阵中心
subplot(211);               % 窗口分割成2*1的，当前是第1个子图 
plot(f2,dmf,'LineWidth',2)  % 画出经过低通滤波后的频谱
title("I路经过低通滤波后的频谱");
xlabel('频率/Hz');          % x轴标签
ylabel('幅度');             % y轴标签

dmf2=fftshift(abs(fft(lvbo2)));%对低通滤波信号采用快速傅里叶变换并移到矩阵中心
subplot(212);               % 窗口分割成2*1的，当前是第2个子图 
plot(f2,dmf2,'LineWidth',2) % 画出经过低通滤波后的频谱
title("Q路经过低通滤波后的频谱");
xlabel('频率/Hz');          % x轴标签
ylabel('幅度');             % y轴标签



