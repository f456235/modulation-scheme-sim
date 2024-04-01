clear all;                  % 清除所有变量
close all;                  % 关闭所有窗口
clc;                        % 清屏
%% 基本参数
M=12;                       % 产生码元数    
L=100;                      % 每码元复制L次,每个码元采样次数
Ts=0.001;                   % 每个码元的宽度,即码元的持续时间
Rb=1/Ts;                    % 码元速率1K
dt=Ts/L;                    % 采样间隔
TotalT=M*Ts;                % 总时间1
t=0:dt:TotalT-dt;           % 时间1
TotalT2=(M/2)*Ts;           % 总时间2
t2=0:dt:TotalT2-dt;         % 时间2
Fs=1/dt;                    % 采样间隔的倒数即采样频率

%% 产生单极性波形
wave=[1,0,1,0,1,1,0,1,0,0,1,1];      % 产生二进制随机码,M为码元个数
fz=ones(1,L);               % 定义复制的次数L,L为每码元的采样点数
fz2=ones(1,2*L);            % 定义复制的次数2L
x1=wave(fz,:);              % 将原来wave的第一行复制L次，称为L*M的矩阵
jidai=reshape(x1,1,L*M);    % 产生单极性不归零矩形脉冲波形，将刚得到的L*M矩阵，按列重新排列形成1*(L*M)的矩阵

%% 产生I、Q路码元
a=wave(1:2:end);            % I路码元是原始码元奇数位置的码元
b=wave(2:2:end);            % Q路码元是原始码元偶数位置的码元
% 码元电平转换，要将两两码元转换成四种电平，00转换成-1,01转换成-3,10转换成1,11转换成3
for i=1:length(a)/2
    if a(2*(i-1)+1:2*i)==[0 1]
        I(i)=-3;
    elseif a(2*(i-1)+1:2*i)==[0 0]
        I(i)=-1;
    elseif a(2*(i-1)+1:2*i)==[1 0]
        I(i)=1;
    elseif a(2*(i-1)+1:2*i)==[1 1]
        I(i)=3;
    end
        if b(2*(i-1)+1:2*i)==[0 1]
        Q(i)=-3;
    elseif b(2*(i-1)+1:2*i)==[0 0]
        Q(i)=-1;
    elseif b(2*(i-1)+1:2*i)==[1 0]
        Q(i)=1;
    elseif b(2*(i-1)+1:2*i)==[1 1]
        Q(i)=3;
    end
end
x2=a(fz,:);                    % 将原来a的第一行复制L次，称为L*(M/2)的矩阵
a_mayuan=reshape(x2,1,L*(M/2));% 产生单极性不归零矩形脉冲波形，将刚得到的L*(M/2)矩阵，按列重新排列形成1*(L*(M/2))的矩阵
x3=b(fz,:);                    % 将原来b的第一行复制L次，称为L*(M/2)的矩阵
b_mayuan=reshape(x3,1,L*(M/2));% 产生单极性不归零矩形脉冲波形，将刚得到的L*(M/2)矩阵，按列重新排列形成1*(L*(M/2))的矩阵
x4=I(fz2,:);                   % 将原来I的第一行复制2L次，称为2L*(M/4)的矩阵
I_mayuan=reshape(x4,1,(2*L)*(M/4));% 产生单极性不归零矩形脉冲波形，将刚得到的(2*L)*(M/4)矩阵，按列重新排列形成1*((2*L)*(M/4))的矩阵
x5=Q(fz2,:);                   % 将原来Q的第一行复制2L次，称为2L*(M/4)的矩阵
Q_mayuan=reshape(x5,1,(2*L)*(M/4));% 产生单极性不归零矩形脉冲波形，将刚得到的(2*L)*(M/4)矩阵，按列重新排列形成1*((2*L)*(M/4))的矩阵

%% 码元波形绘制
figure(1);                  % 绘制第1幅图
subplot(511);               % 窗口分割成5*1的，当前是第1个子图 
plot(t,jidai,'LineWidth',2);% 绘制基带码元波形，线宽为2
title('基带信号波形');      % 标题
xlabel('时间/s');           % x轴标签
ylabel('幅度');             % y轴标签
axis([0,TotalT,-0.1,1.1])   % 坐标范围限制

subplot(512);               % 窗口分割成5*1的，当前是第2个子图 
plot(t2,a_mayuan,'LineWidth',2);% I路码元波形，线宽为2
title('I路码元波形');       % 标题
xlabel('时间/s');           % x轴标签
ylabel('幅度');             % y轴标签
axis([0,TotalT2,-0.1,1.1])  % 坐标范围限制

subplot(513);               % 窗口分割成5*1的，当前是第3个子图 
plot(t2,b_mayuan,'LineWidth',2);% 绘制Q路码元波形，线宽为2
title('Q路码元波形');       % 标题
xlabel('时间/s');           % x轴标签
ylabel('幅度');             % y轴标签
axis([0,TotalT2,-0.1,1.1])  % 坐标范围限制

subplot(514);               % 窗口分割成5*1的，当前是第4个子图 
plot(t2,I_mayuan,'LineWidth',2);% 绘制转换后I路电平波形，线宽为2
title('电平转换后I路码元波形');% 标题
xlabel('时间/s');           % x轴标签
ylabel('幅度');             % y轴标签
axis([0,TotalT2,-4,4])      % 坐标范围限制

subplot(515);               % 窗口分割成5*1的，当前是第5个子图 
plot(t2,Q_mayuan,'LineWidth',2);% 绘制转换后Q路电平波形，线宽为2
title('电平转换后Q路码元波形');% 标题
xlabel('时间/s');           % x轴标签
ylabel('幅度');             % y轴标签
axis([0,TotalT2,-4,4])      % 坐标范围限制

%% QAM调制
fc=10000;                   % 载波频率2kHz       
zb1=cos(2*pi*fc*t2);        % 载波1
zb2=-sin(2*pi*fc*t2);       % 载波2
I_lu=I_mayuan.*zb1;         % I路波形 
Q_lu=Q_mayuan.*zb2;         % Q路波形 
qam=I_lu+Q_lu;              % QAM调制
figure(2);                  % 绘制第2幅图
subplot(411)                % 窗口分割成4*1的，当前是第1个子图 
plot(t2,I_lu,'LineWidth',2);% 绘制I路的波形 
title('I路信号波形')        % 标题
axis([0,TotalT2,-4,4]);     % 坐标范围限制
xlabel('时间/s');           % x轴标签
ylabel('幅度');             % y轴标签

subplot(412)                % 窗口分割成4*1的，当前是第2个子图 
plot(t2,Q_lu,'LineWidth',2);% 绘制Q路的波形 
title('Q路信号波形')        % 标题
axis([0,TotalT2,-4,4]);     % 坐标范围限制
xlabel('时间/s');           % x轴标签
ylabel('幅度');             % y轴标签

subplot(413)                % 窗口分割成4*1的，当前是第3个子图 
plot(t2,qam,'LineWidth',2); % 绘制QAM的波形 
title('QAM信号波形')        % 标题
axis([0,TotalT2,-7,7]);     % 坐标范围限制
xlabel('时间/s');           % x轴标签
ylabel('幅度');             % y轴标签
%% 信号经过高斯白噪声信道
tz=awgn(qam,20);            % 信号qam中加入白噪声，信噪比为SNR=20dB
subplot(414);               % 窗口分割成4*1的，当前是第4个子图 
plot(t2,tz,'LineWidth',2);  % 绘制QAM信号加入白噪声的波形
axis([0,TotalT2,-7,7]);     % 坐标范围设置
title('通过高斯白噪声信道后的信号');% 标题
xlabel('时间/s');           % x轴标签
ylabel('幅度');             % y轴标签

%% 解调部分
jt1=tz.*zb1;                % 相干解调，I路乘以相干载波
jt2=tz.*zb2;                % 相干解调，Q路乘以相干载波
figure(3);
subplot(511)                % 窗口分割成5*1的，当前是第1个子图 
plot(t2,jt1,'LineWidth',2)  % 绘制I路乘以相干载波后的信号
axis([0,TotalT2,-5,5]);     % 设置坐标范围
title("I路乘以相干载波后的信号")% 标题
xlabel('时间/s');           % x轴标签
ylabel('幅度');             % y轴标签

subplot(512)                % 窗口分割成5*1的，当前是第2个子图 
plot(t2,jt2,'LineWidth',2)  % 绘制Q路乘以相干载波后的信号
axis([0,TotalT2,-5,5]);     % 设置坐标范围
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
lvbo1=2*fftfilt(b,jt1);     % 对信号进行滤波，jt1是等待滤波的信号，b是fir滤波器的系统函数的分子多项式系数
lvbo2=2*fftfilt(b,jt2);     % 对信号进行滤波，jt2是等待滤波的信号，b是fir滤波器的系统函数的分子多项式系数

subplot(513);               % 窗口分割成5*1的，当前是第3个子图 
plot(w/pi*Fs/2,20*log(abs(h)),'LineWidth',2); % 绘制滤波器的幅频响应
title('低通滤波器的频谱');  % 标题
xlabel('频率/Hz');          % x轴标签
ylabel('幅度/dB');          % y轴标签

subplot(514)                % 窗口分割成5*1的，当前是第4个子图 
plot(t2,lvbo1,'LineWidth',2);% 绘制I路经过低通滤波器后的信号
axis([0,TotalT2,-5,5]);     % 设置坐标范围
title("I路经过低通滤波器后的信号");% 标题
xlabel('时间/s');           % x轴标签
ylabel('幅度');             % y轴标签

subplot(515)                % 窗口分割成5*1的，当前是第5个子图 
plot(t2,lvbo2,'LineWidth',2);% 绘制Q路经过低通滤波器后的信号
axis([0,TotalT2,-5,5]);     % 设置坐标范围
title("I路经过低通滤波器后的信号");% 标题
xlabel('时间/s');           % x轴标签
ylabel('幅度');             % y轴标签

%% 抽样判决
% 滤波后的波形进行判决，大于2的判决为3,0-2之间的判决为1，-2-0之间判决为-1，小于-2的判决为-3
for j=1:length(lvbo1)
   if lvbo1(j)>=2
       I_panjue(j)=3;
   elseif (lvbo1(j)>0 && lvbo1(j)<2)
       I_panjue(j)=1;
   elseif (lvbo1(j)>=-2 && lvbo1(j)<0)
       I_panjue(j)=-1;
   elseif lvbo1(j)<-2 
       I_panjue(j)=-3;
   end
end

for k=1:length(lvbo2)
   if lvbo2(k)>=2
       Q_panjue(k)=3;
   elseif (lvbo2(k)>0 && lvbo2(k)<2)
       Q_panjue(k)=1;
   elseif (lvbo2(k)>=-2 && lvbo2(k)<0)
       Q_panjue(k)=-1;
   elseif lvbo2(k)<-2 
       Q_panjue(k)=-3;
   end
end
figure(4);
subplot(611)                % 窗口分割成6*1的，当前是第1个子图 
plot(t2,I_panjue,'LineWidth',2);% 绘制经过判决的I路信号
axis([0,TotalT2,-5,5]);     % 设置坐标范围
title("I路经过判决后的信号");% 标题
xlabel('时间/s');           % x轴标签
ylabel('幅度');             % y轴标签

subplot(612)                % 窗口分割成6*1的，当前是第2个子图 
plot(t2,Q_panjue,'LineWidth',2);% 绘制经过判决的Q路信号
axis([0,TotalT2,-5,5]);     % 设置坐标范围
title("Q路经过判决后的信号");% 标题
xlabel('时间/s');           % x轴标签
ylabel('幅度');             % y轴标签

%% 判决之后再处理，取码元中间时刻判决
I_zong=[];
Q_zong=[];
for m=L:2*L:((2*L)*(M/4))
   if I_panjue(m)>=2
       I_zong=[I_zong,3];
   elseif (I_panjue(m)>0 && I_panjue(m)<2)
       I_zong=[I_zong,1];
   elseif (I_panjue(m)>=-2 && I_panjue(m)<0)
       I_zong=[I_zong,-1];
   elseif I_panjue(m)<-2 
       I_zong=[I_zong,-3];
   end
end

for n=L:2*L:((2*L)*(M/4))
   if Q_panjue(n)>=2
       Q_zong=[Q_zong,3];
   elseif (Q_panjue(n)>0 && Q_panjue(n)<2)
       Q_zong=[Q_zong,1];
   elseif (Q_panjue(n)>=-2 && Q_panjue(n)<0)
       Q_zong=[Q_zong,-1];
   elseif Q_panjue(n)<-2 
       Q_zong=[Q_zong,-3];
   end
end

x6=I_zong(fz2,:);                   % 将原来I_zong的第一行复制2L次，称为2L*(M/4)的矩阵
I_wave=reshape(x6,1,(2*L)*(M/4));% 产生单极性不归零矩形脉冲波形，将刚得到的(2*L)*(M/4)矩阵，按列重新排列形成1*((2*L)*(M/4))的矩阵
x7=Q_zong(fz2,:);                   % 将原来Q_zong的第一行复制2L次，称为2L*(M/4)的矩阵
Q_wave=reshape(x7,1,(2*L)*(M/4));% 产生单极性不归零矩形脉冲波形，将刚得到的(2*L)*(M/4)矩阵，按列重新排列形成1*((2*L)*(M/4))的矩阵

subplot(613)                % 窗口分割成6*1的，当前是第3个子图 
plot(t2,I_wave,'LineWidth',2);% 绘制经过中间时刻判决的I路波形
axis([0,TotalT2,-5,5]);     % 设置坐标范围
title("I路经过中间时刻判决后的信号");% 标题
xlabel('时间/s');           % x轴标签
ylabel('幅度');             % y轴标签

subplot(614)                % 窗口分割成6*1的，当前是第4个子图 
plot(t2,Q_wave,'LineWidth',2);% 绘制经过中间时刻判决的Q路波形
axis([0,TotalT2,-5,5]);     % 设置坐标范围
title("Q路经过中间时刻判决后的信号");% 标题
xlabel('时间/s');           % x轴标签
ylabel('幅度');             % y轴标签

%% I、Q两路恢复
% 电平转换成码元，-3转换成01，-1转换成00,1转换成10,3转换成11
I_yuan=[];
Q_yuan=[];
for i=1:length(I_zong)
    if I_zong(i)==3
        I_yuan=[I_yuan,1 1];
    elseif I_zong(i)==1
        I_yuan=[I_yuan,1 0];
    elseif I_zong(i)==-1
        I_yuan=[I_yuan,0 0];
    elseif I_zong(i)==-3
        I_yuan=[I_yuan,0 1];
    end
end

for i=1:length(Q_zong)
    if Q_zong(i)==3
        Q_yuan=[Q_yuan,1 1];
    elseif Q_zong(i)==1
        Q_yuan=[Q_yuan,1 0];
    elseif Q_zong(i)==-1
        Q_yuan=[Q_yuan,0 0];
    elseif Q_zong(i)==-3
        Q_yuan=[Q_yuan,0 1];
    end
end

code = [];
% 将I路码元为最终输出的奇数位置码元，将Q路码元为最终输出的偶数位置码元
for n=1:M
    if mod(n, 2)~=0
        code = [code, I_yuan((n+1)/2)];
    else
        code = [code, Q_yuan(n/2)];
    end
end

x8=wave(fz,:);              % 将原来wave的第一行复制L次，称为L*M的矩阵
dout=reshape(x8,1,L*M);     % 产生单极性不归零矩形脉冲波形，将刚得到的L*M矩阵，按列重新排列形成1*(L*M)的矩阵

subplot(615)                % 窗口分割成6*1的，当前是第5个子图 
plot(t,dout,'LineWidth',2); % 绘制恢复后的基带码元
axis([0,TotalT,-0.1,1.1]);  % 设置坐标范围
title("恢复后的基带码元");  % 标题
xlabel('时间/s');           % x轴标签
ylabel('幅度');             % y轴标签

subplot(616)                % 窗口分割成6*1的，当前是第6个子图 
plot(t,jidai,'LineWidth',2);% 绘制原始基带码元
axis([0,TotalT,-0.1,1.1]);  % 设置坐标范围
title("原始基带码元");      % 标题
xlabel('时间/s');           % x轴标签
ylabel('幅度');             % y轴标签

%% 绘制频谱
%% 信源频谱
T=t(end);                   % 时间
df=1/T;                     % 频谱分辨率
N=length(jidai);            % 采样长度
f=(-N/2:N/2-1)*df;          % 频率范围
mf=fftshift(abs(fft(jidai)));%对信源信号采用快速傅里叶变换并移到矩阵中心
figure(5);
subplot(211);               % 窗口分割成2*1的，当前是第1个子图 
plot(f,mf,'LineWidth',2);   % 绘制信源频谱波形
title("基带信号频谱");      % 标题
xlabel('频率/Hz');          % x轴标签
ylabel('幅度');             % y轴标签
axis([-20000,20000,-inf,inf]);% 坐标范围限制

%% QAM信号频谱
T2=t2(end);                 % 时间2
df2=1/T2;                   % 频谱分辨率2
N2=length(qam);             % 采样长度2
f2=(-N2/2:N2/2-1)*df2;      % 频率范围2
sf=fftshift(abs(fft(qam))); % 对QAM信号采用快速傅里叶变换并将0-fs频谱移动到-fs/2-fs/2

subplot(212)                % 窗口分割成2*1的，当前是第2个子图 
plot(f2,sf,'LineWidth',2)   % 绘制QAM调制信号频谱
title("QAM信号频谱")        % 标题
xlabel('频率/Hz');          % x轴标签
ylabel('幅度');             % y轴标签
axis([-20000,20000,-inf,inf]);% 坐标范围限制


%% 乘以相干载波后的频谱
mmf=fftshift(abs(fft(jt1))); % 对相干载波信号采用快速傅里叶变换并移到矩阵中心
figure(6)                   % 绘制第5幅图
subplot(211);               % 窗口分割成2*1的，当前是第1个子图 
plot(f2,mmf,'LineWidth',2)   % 画出I路乘以相干载波后的频谱
title("I路乘以相干载波后的频谱")
xlabel('频率/Hz');          % x轴标签
ylabel('幅度');             % y轴标签

mmf2=fftshift(abs(fft(jt2))); % 对相干载波信号采用快速傅里叶变换并移到矩阵中心
subplot(212);               % 窗口分割成2*1的，当前是第2个子图 
plot(f2,mmf2,'LineWidth',2)   % 画出Q路乘以相干载波后的频谱
title("Q路乘以相干载波后的频谱")
xlabel('频率/Hz');          % x轴标签
ylabel('幅度');             % y轴标签

%% 经过低通滤波后的频谱
dmf=fftshift(abs(fft(lvbo1)));%对低通滤波信号采用快速傅里叶变换并移到矩阵中心
figure(7);
subplot(211);               % 窗口分割成2*1的，当前是第1个子图 
plot(f2,dmf,'LineWidth',2)   % 画出I路经过低通滤波后的频谱
title("I路经过低通滤波后的频谱");
xlabel('频率/Hz');          % x轴标签
ylabel('幅度');             % y轴标签

dmf2=fftshift(abs(fft(lvbo2)));%对低通滤波信号采用快速傅里叶变换并移到矩阵中心
subplot(212);               % 窗口分割成2*1的，当前是第2个子图 
plot(f2,dmf2,'LineWidth',2)   % 画出Q路经过低通滤波后的频谱
title("Q路经过低通滤波后的频谱");
xlabel('频率/Hz');          % x轴标签
ylabel('幅度');             % y轴标签
