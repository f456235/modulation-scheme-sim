clear all;                  % 清除所有变量
close all;                  % 关闭所有窗口
clc;                        % 清屏

%% 基本参数
M = 12;                     % 产生码元数    
L = 100;                    % 每码元复制L次,每个码元采样次数
Ts = 0.001;                 % 每个码元的宽度,即码元的持续时间
Rb = 1 / Ts;                % 码元速率1K
dt = Ts / L;                % 采样间隔
TotalT = M * Ts;            % 总时间
t = 0:dt:TotalT-dt;         % 时间
Fs = 1 / dt;                % 采样间隔的倒数即采样频率

%% 产生单极性波形
wave = [1, 0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 1]; % 产生二进制随机码,M为码元个数
fz = ones(1, L);            % 定义复制的次数L,L为每码元的采样点数
x1 = wave(fz, :);           % 将原来wave的第一行复制L次，称为L*M的矩阵
jidai = reshape(x1, 1, L * M); % 产生单极性不归零矩形脉冲波形，将刚得到的L*M矩阵，按列重新排列形成1*(L*M)的矩阵

%% 2FSK调制
fc1 = 2000;                 % 载波1频率2kHz       
zb1 = cos(2 * pi * fc1 * t);% 载波1信号
fsk1 = jidai .* zb1;        % 基带信号1和载波1相乘 

fc2 = 10000;                % 载波2频率10kHz       
zb2 = cos(2 * pi * fc2 * t);% 载波2信号
fsk2 = (1 - jidai) .* zb2;  % 基带信号2和载波2相乘  

fsk = fsk1 + fsk2;          % 2FSK的调制

%% 解调部分
%% fsk信号经过第一个带通滤波后和载波1相乘
jt1 = fsk1 .* (-zb1);       % 相干解调，乘以相干载波

%% fsk信号经过第二个带通滤波后和载波2相乘
jt2 = fsk2 .* (-zb2);       % 相干解调，乘以相干载波

%% 抽样判决
pdst = 1 * (jt1 < jt2);     % 上面支路信号大于下面支路信号判决为1，否则为0 

%% 绘制输入信号、调制信号、载波和解调信号
figure;
subplot(511);               % 窗口分割成5*1的，当前是第1个子图 
plot(t, jidai, 'LineWidth', 2) % 绘制基带信号波形
title('Input Signal')       % 标题
xlabel('Time/s');           % x轴标签
ylabel('Amplitude');        % y轴标签
axis([0, TotalT, -0.1, 1.1])% 坐标范围限制

subplot(512);               % 窗口分割成5*1的，当前是第2个子图 
plot(t, fsk, 'LineWidth', 2);% 绘制2FSK的波形 
title('Modulated Signal')   % 标题
xlabel('Time/s');           % x轴标签
ylabel('Amplitude');        % y轴标签
axis([0, TotalT, -1.1, 1.1]);% 坐标范围限制

subplot(513);               % 窗口分割成5*1的，当前是第3个子图 
plot(t, zb1, 'r', 'LineWidth', 2);% 绘制载波1的波形
hold on;
plot(t, zb2, 'b', 'LineWidth', 2);% 绘制载波2的波形
title('Carrier Signals')   % 标题
legend('Carrier 1', 'Carrier 2'); % 图例
xlabel('Time/s');           % x轴标签
ylabel('Amplitude');        % y轴标签

subplot(514);               % 窗口分割成5*1的，当前是第4个子图 
plot(t, fsk1, 'r', 'LineWidth', 2);% 绘制解调信号1的波形
hold on;
plot(t, fsk2, 'b', 'LineWidth', 2);% 绘制解调信号2的波形
title('Demodulated Signals')% 标题
legend('Demodulated Signal 1', 'Demodulated Signal 2'); % 图例
xlabel('Time/s');           % x轴标签
ylabel('Amplitude');        % y轴标签

subplot(515);               % 窗口分割成5*1的，当前是第5个子图 
plot(t, pdst, 'b', 'LineWidth', 2); % Plotting the demodulated signal
hold on;
plot(t, jidai, 'r--', 'LineWidth', 2); % Plotting the original input signal
hold off;
title('Sampled and Demodulated Signal vs Original Signal');
xlabel('Time (s)');
ylabel('Amplitude');
legend('Demodulated Signal', 'Original Signal');