clear all;                  % Clear all variables
close all;                  % Close all windows
clc;                        % Clear command window

%% Basic parameters
M=12;                       % Number of symbols generated
L=100;                      % Each symbol repeats L times, sampling points per symbol
Ts=0.001;                   % Width of each symbol, symbol duration
Rb=1/Ts;                    % Symbol rate 1K
dt=Ts/L;                    % Sampling interval
TotalT=M*Ts;                % Total time
t=0:dt:TotalT-dt;           % Time vector
Fs=1/dt;                    % Sampling frequency

%% Generate unipolar waveform
wave=[1,0,1,0,1,1,0,1,0,0,1,1];      % Generate binary random code, M is the number of symbols
fz=ones(1,L);               % Define the number of copies L, L is the number of sampling points per symbol
x1=wave(fz,:);              % Repeat the first row of wave L times to form an L*M matrix
jidai=reshape(x1,1,L*M);    % Generate unipolar non-return-to-zero rectangular pulse waveform, rearrange the L*M matrix into a 1*(L*M) matrix

%% Unipolar to Bipolar
% Convert the baseband signal to bipolar, where jidai is 1 when jidai is 1; jidai is -1 when jidai is 0.
for n=1:length(jidai)
    if jidai(n)==1
        jidai(n)=1;
    else
        jidai(n)=-1;
    end
end

%% 2PSK modulation
fc=2000;                    % Carrier frequency 2kHz       
zb=sin(2*pi*fc*t);          % Carrier
psk=jidai.*zb;              % Analog modulation of 2PSK
figure(1);                  % Plot figure 1
subplot(211);               % Divide the window into 2*1, this is the first subplot
plot(t,jidai,'LineWidth',2);% Plot baseband symbol waveform, line width 2
title('Baseband signal waveform');      % Title
xlabel('Time/s');           % x-axis label
ylabel('Amplitude');        % y-axis label
axis([0,TotalT,-1.1,1.1])   % Coordinate range limit

subplot(212)                % Divide the window into 2*1, this is the second subplot
plot(t,psk,'LineWidth',2);  % Plot PASK waveform 
title('2PSK signal waveform')   % Title
axis([0,TotalT,-1.1,1.1]);  % Coordinate range limit
xlabel('Time/s');           % x-axis label
ylabel('Amplitude');        % y-axis label

%% Signal passing through Gaussian white noise channel
tz=awgn(psk,100);            % Add white noise to signal psk, SNR=15dB
figure(2);                  % Plot figure 2
subplot(211);               % Divide the window into 2*1, this is the first subplot 
plot(t,tz,'LineWidth',2);   % Plot waveform of 2PSK signal after adding white noise
axis([0,TotalT,-1.5,1.5]);  % Coordinate range setting
title('Signal after passing through Gaussian white noise channel');% Title
xlabel('Time/s');           % x-axis label
ylabel('Amplitude');        % y-axis label

%% Demodulation
tz=tz.*zb;                  % Coherent demodulation, multiply by coherent carrier
subplot(212)                % Divide the window into 2*1, this is the second subplot 
plot(t,tz,'LineWidth',1)    % Plot signal after multiplying coherent carrier
axis([0,TotalT,-1.5,1.5]);  % Set coordinate range
title("Signal after multiplying coherent carrier")% Title
xlabel('Time/s');           % x-axis label
ylabel('Amplitude');        % y-axis label

%% Filter for noise-added signal
% Low-pass filter design
fp=2*Rb;                    % Low-pass filter cut-off frequency, multiplied by 2 because the analog frequency is converted to digital frequency wp=Rb/(Fs/2)
b=fir1(30, fp/Fs, boxcar(31));% Generate coefficients of the numerator polynomial in the fir filter system function
% The three parameters of the fir1 function are order, digital cut-off frequency, and filter type.
% Here, a 30th-order (31 tap coefficient) rectangular window filter is generated
[h,w]=freqz(b, 1,512);      % Generate frequency response of fir filter
% The three parameters of the freqz function are the coefficients of the numerator polynomial of the filter system function, the coefficients of the denominator polynomial (the numerator coefficients of the fir filter are 1), and the number of sampling points (default) 512
lvbo=fftfilt(b,tz);         % Filter the signal, tz is the signal to be filtered, b is the coefficient of the numerator polynomial of the fir filter system function
figure(3);                  % Plot figure 3  
subplot(311);               % Divide the window into 3*1, this is the first subplot 
plot(w/pi*Fs/2,20*log(abs(h)),'LineWidth',2); % Plot the amplitude-frequency response of the low-pass filter
title('Frequency spectrum of low-pass filter');  % Title
xlabel('Frequency/Hz');          % x-axis label
ylabel('Amplitude/dB');          % y-axis label

subplot(312)                % Divide the window into 3*1, this is the second subplot 
plot(t,lvbo,'LineWidth',2); % Plot the signal after passing through the low-pass filter
axis([0,TotalT,-1.1,1.1]);  % Set coordinate range
title("Signal after low-pass filtering");% Title
xlabel('Time/s');           % x-axis label
ylabel('Amplitude');        % y-axis label

%% Sampling decision
k=0;                        % Set the sampling threshold
pdst=1*(lvbo>0);            % Compare each element of the filtered vector with 0, 1 if greater than 0, otherwise 0
subplot(313)                % Divide the window into 2*1, this is the third subplot 
plot(t,pdst,'LineWidth',2)  % Plot the signal after sampling decision
axis([0,TotalT,-0.1,1.1]);  % Set coordinate range
title("Signal after sampling decision")% Title
xlabel('Time/s');           % x-axis label
ylabel('Amplitude');        % y-axis label

%% Plot frequency spectrum
%% Frequency spectrum of 2PSK signal
T=t(end);                   % Time
df=1/T;                     % Frequency resolution
N=length(psk);              % Sampling length
f=(-N/2:N/2-1)*df;          % Frequency range
sf=fftshift(abs(fft(psk))); % Apply fast Fourier transform to 2PSK signal and move 0-Fs spectrum to -Fs/2-Fs/2
figure(4)                   % Plot figure 4
subplot(211)                % Divide the window into 2*1, this is the first subplot 
plot(f,sf,'LineWidth',2)    % Plot modulation signal spectrum
title("2PSK signal spectrum")       % Title
xlabel('Frequency/Hz');          % x-axis label
ylabel('Amplitude');             % y-axis label

%% Source spectrum
mf=fftshift(abs(fft(jidai)));% Apply fast Fourier transform to source signal and move to matrix center
subplot(212);               % Divide the window into 2*1, this is the second subplot 
plot(f,mf,'LineWidth',2);   % Plot source spectrum waveform
title("Baseband signal spectrum");      % Title
xlabel('Frequency/Hz');          % x-axis label
ylabel('Amplitude');             % y-axis label

%% Spectrum after multiplying coherent carrier
mmf=fftshift(abs(fft(tz))); % Apply fast Fourier transform to coherent carrier signal and move to matrix center
figure(5)                   % Plot figure 5
subplot(211);               % Divide the window into 2*1, this is the first subplot 
plot(f,mmf,'LineWidth',2)   % Plot spectrum after multiplying coherent carrier
title("Spectrum after multiplying coherent carrier")
xlabel('Frequency/Hz');          % x-axis label
ylabel('Amplitude');             % y-axis label

%% Spectrum after low-pass filtering
dmf=fftshift(abs(fft(lvbo)));% Apply fast Fourier transform to low-pass filtered signal and move to matrix center
subplot(212);               % Divide the window into 2*1, this is the second subplot 
plot(f,dmf,'LineWidth',2)   % Plot spectrum after low-pass filtering
title("Spectrum after low-pass filtering");
xlabel('Frequency/Hz');          % x-axis label
ylabel('Amplitude');             % y-axis label
