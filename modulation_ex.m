% Generate analog signal
Fs = 2000; % Sampling frequency
t = 0:1/Fs:1; % Time vector
analog_signal = sin(2*pi*50*t); % Example analog signal (50 Hz sinusoidal signal)

% Generate message bits
message_bits = [1 0 1 0 1 1 0 1 0 0 1 1];

% AM modulation of analog signal
am_signal = ammod(analog_signal, 2000, 1000, 0);

% FM modulation of analog signal
fm_signal = fmmod(analog_signal, 2000, 1000, 50); % Example modulation index of 50

% FSK modulation of analog signal
fsk_signal = fskmod(analog_signal, 2000, 1000, 50); % Example frequency deviation of 50 Hz

% PSK modulation of message bits
psk_signal = pskmod(message_bits, 2, pi);

% QPSK modulation of message bits
qpsk_signal = pskmod(message_bits, 4, pi/4);

% π/4QPSK modulation of message bits
pi_4qpsk_signal = pi_4qpskmod(message_bits);

% QAM modulation of message bits
qam_signal = qammod(message_bits, 16);

% 16QAM modulation of message bits
qam16_signal = qammod(message_bits, 16);

% Plotting
figure;
subplot(3,2,1);
plot(t, analog_signal);
title('Analog Signal');
xlabel('Time');
ylabel('Amplitude');

subplot(3,2,2);
plot(t, am_signal);
title('AM Modulated Signal');
xlabel('Time');
ylabel('Amplitude');

subplot(3,2,3);
plot(t, fm_signal);
title('FM Modulated Signal');
xlabel('Time');
ylabel('Amplitude');

subplot(3,2,4);
plot(t, fsk_signal);
title('FSK Modulated Signal');
xlabel('Time');
ylabel('Amplitude');

subplot(3,2,5);
stem(message_bits, real(psk_signal));
title('PSK Modulated Signal');
xlabel('Time');
ylabel('Amplitude');

subplot(3,2,6);
stem(message_bits, real(qpsk_signal));
title('QPSK Modulated Signal');
xlabel('Time');
ylabel('Amplitude');