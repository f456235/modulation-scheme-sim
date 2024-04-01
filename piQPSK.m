clear all;
close all;

% The number of bits to send - Frame Length
N = 12;

% Generate a random bit stream
bit_stream = [1,0,1,0,1,1,0,1,0,0,1,1];

% 4 PHASE SHIFTS
P1 = pi/4;      % 45degrees phase shift
P2 = 3*pi/4;    % 135 degrees phase shift
P3 = 5*pi/4;    % 225 degree phase shift
P4 = 7*pi/4;    % 315 degree phase shift

% Frequency of Modulating Signal
f = 1;          % f --> time period

% Sampling rate of sine wave - This will define the resolution
fs = 100;

% Time for one bit
t = 0: 1/fs : 1;

% This time variable is just for plot
time = [];
time1 = [];
QPSK_signal = [];
Digital_signal = [];
carrier_signal=[];

for i = 1: length(bit_stream)
    current_bit = bit_stream(i);
    
    if current_bit == 0
        % If the current bit is 0, append a sequence of zeros
        segment = zeros(1, length(t));
    else
        % If the current bit is 1, append a sequence of ones
        segment = ones(1, length(t));
    end

    % Append the segment based on the next bit value
    
    Digital_signal = [Digital_signal, segment];
    time1 = [time1 t];
    t =  t + 1;
end

t = 0: 1/fs : 1;
for ii = 1:length(bit_stream)
  jj = ii + 1;

if ii == length(bit_stream)
    jj = ii;
end

%Code for generation of carrier signal
 carrier_signal=[carrier_signal (sin(2*pi*f*t))];
%Code for genearting QPSK signal modulated signal
if bit_stream(ii)==0
    if bit_stream(jj)==0
        bit00 = (bit_stream(ii)==0)*sin(2*pi*f*t + P1);
        QPSK_signal = [QPSK_signal (bit00)];
    else
       bit0 = (bit_stream(ii)==0)*sin(2*pi*f*t + P2);
       bit1 = (bit_stream(jj)==0)*sin(2*pi*f*t + P2);
       QPSK_signal = [QPSK_signal (bit0+bit1) ];
    end
end
if bit_stream(ii)==1
     if bit_stream(jj)==0
         bit1 = (bit_stream(ii)==0)*sin(2*pi*f*t + P3);
         bit0 = (bit_stream(jj)==0)*sin(2*pi*f*t + P3);
        QPSK_signal = [QPSK_signal (bit1+bit0) ];
     else
         bit11 = (bit_stream(jj)==1)*sin(2*pi*f*t + P4);
         QPSK_signal = [QPSK_signal (bit11) ];
     end
end
    time = [time t];
    t =  t + 1;
end

% Define the received Pi/4 QPSK signal
received_signal = QPSK_signal; % Assuming QPSK_signal is your received signal

T = 12; % Symbol duration (in seconds)

% Generate the two local oscillator signals (in-phase and quadrature)
t = 0:1/fs:1211/fs;
local_oscillator_I = cos(2*pi*f*t);
local_oscillator_Q = sin(2*pi*f*t);

% Perform coherent detection
in_phase_component = received_signal .* local_oscillator_I;
quadrature_component = received_signal .* local_oscillator_Q;

% Low-pass filtering (you can adjust the filter order as needed)
filter_order = 10;
lpf = ones(1, filter_order) / filter_order;
filtered_I = conv(in_phase_component, lpf, 'same');
filtered_Q = conv(quadrature_component, lpf, 'same');

% Combine the filtered in-phase and quadrature components to obtain the demodulated signal
demodulated_signal = filtered_I + filtered_Q;



% Plot the Original Digital Signal
subplot(7,1,1);
plot(time1,Digital_signal,'r','LineWidth',2);
xlabel('Time (bit period)');
ylabel('Amplitude');
title('Original Digital Signal');
axis([0 8 -0.5 1.5]);
grid on;

% Plot the carrier Signal
subplot(7,1,2);
plot(time,carrier_signal,'g','LineWidth',2);
xlabel('Time (bit period)');
ylabel('Amplitude');
title('Carrier Signal');
axis([0 time(end) -1.5 1.5]);
grid  on;

% Plot the QPSK Signal
subplot(7,1,3);
plot(time, QPSK_signal,'LineWidth',2);
xlabel('Time (bit period)');
ylabel('Amplitude');
title('QPSK Signal with two Phase Shifts');
axis([0 8 -1.5 1.5]);
grid  on;

% Plot the Demodulated Signal
subplot(7,1,4);
plot(demodulated_signal);
xlabel('Symbol Index');
ylabel('Decoded Symbol');
title('Demodulated Signal');
grid on;


