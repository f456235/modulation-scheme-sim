function[t,st]=lpf(f,sf,B)          % 频率LPF
df=f(2)-f(1);                       % 频率间隔
fN=length(f);                       % 采样点数
ym=zeros(1,fN);                     % 生成1行fN列的0向量
xm=floor(B/df);                     % 低频带宽频率除以间隔后的点数向下取整
xm_shift=[-xm:xm-1]+floor(fN/2);    % 因为前面做FFT将0频率搬移到中心处，
                                    % 因此，低通低频频率相应地搬移fN/2，才是对应的频率点
ym(xm_shift)=1;                     % 低通通过频率处幅度为1，其余为0，相当于理想低通
yf=ym.* sf;                         % FFT信号的频谱和对应低频带宽处频率值为1的行向量相乘
[t,st]=F2T(f,yf);                   % IFFT
