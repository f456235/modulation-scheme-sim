%% FM modulation and demodulation process with Digital-to-Analog Conversion (DAC)

%% Basic Parameters
clear all;
close all;
fm = 100;          % Baseband signal frequency
T = 2;             % Signal duration
fs = 20000;        % Sampling frequency (Nyquist rate is fs/2)
dt = 1/fs;         % Time sampling interval
N = T/dt;          % Number of samples
t = [0:N-1]*dt;    % Time vector

%% Baseband Signal Generation (Digital Signal)
digital_signal = [1 0 1 0 1 1 0 1 0 0 1 1];  % Example digital signal
Am = 1;                                       % Baseband signal amplitude
num_bits = numel(digital_signal);
T = num_bits / fm;                            % Update signal duration based on the number of bits
t = linspace(0, T, num_bits * fs);            % Update time vector based on the signal duration and sampling frequency

mt = zeros(1, numel(t));                      % Initialize baseband signal
for i = 1:num_bits
    mt((i-1)*fs+1:i*fs) = Am * digital_signal(i);
end

%% Carrier Signal Generation
fc = 1000;                                   % Carrier frequency
t_carrier = linspace(0, T, numel(t));       % Time vector for carrier signal
zaibo = cos(2*pi*fc*t_carrier);              % Carrier signal

%% Interpolation of Baseband Signal
mt_resampled = resample(mt, numel(t_carrier), numel(t));

%% Frequency Modulation (FM)
kf = 4000;                                      % Frequency deviation factor
min_length = min(length(mt_resampled), length(zaibo));
mt_resampled = mt_resampled(1:min_length);


SFM = cos(2*pi*fc*t + kf*cumsum(mt_resampled)*dt);

%% Digital-to-Analog Conversion (DAC)
analog_waveform = SFM;

%% Add Gaussian Noise
SNR = 20;                                     % Signal-to-Noise Ratio (SNR) in dB
analog_waveform = awgn(analog_waveform, SNR, 'measured');

%% Differentiation
for i=1:239999
    diff_SFM(i)=(SFM(i+1)-SFM(i))/dt;
end
demodulated_signal = abs((diff_SFM));
demodulated_signal = (demodulated_signal)/kf;
%% Envelope Detection


%% Plotting
figure;

% Plot Baseband Signal
subplot(3, 2, 1);
plot(t, mt, 'Linewidth', 2);
title('Baseband Signal');
xlabel('Time (s)');
ylabel('Amplitude');

% Plot Carrier Signal
subplot(3, 2, 2);
plot(t(1:10000), zaibo(1:10000), 'r', 'Linewidth', 2);
title('Carrier Signal');
xlabel('Time (s)');
ylabel('Amplitude');

% Plot Amplitude of FM Signal
subplot(3, 2, 3);
plot(t, SFM, 'Linewidth', 2);
title('FM Signal');
xlabel('Time (s)');
ylabel('Amplitude');

% Plot Demodulated Signal (after differentiation & after shift phase)
subplot(3, 2, 4);
plot(t(1:length(demodulated_signal)), demodulated_signal, 'b', 'Linewidth', 2);
title('Demodulated Signal');
xlabel('Time (s)');
ylabel('Amplitude');

% Plot Original Digital Signal
subplot(3, 2, 6);
stem(0:numel(digital_signal)-1, digital_signal, 'r', 'Linewidth', 2);
title('Original Digital Signal');
xlabel('Sample Index');
ylabel('Amplitude');


min_length = min(length(digital_signal), length(demodulated_signal));
errors = sum(digital_signal(1:min_length) ~= (demodulated_signal(1:min_length) > 0.5));
error_rate = errors / min_length;

% Compute data rate
data_rate = (min_length - errors) / T; % Correct bits per second

% Print out error rate and data rate
fprintf('Error Rate: %.4f\n', error_rate);
fprintf('Data Rate: %.2f bits per second\n', data_rate);
