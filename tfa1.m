function tfa1(window_duration, frame_shift)
%%

% x is the sample amplitude of the signal
% fs is the sampling frequency of the signal
% spectrogram is the magnitude of STFT of the signal
% for fft purpose, it better to have window length in powers of 2

[x, fs] = audioread('tfa_assg3.wav');
x = x(:, 1); % first channel
% sound(x,fs)
% figure(1);
% plot(x);

%fs = fs*0.01;
Ts = 1/fs;
len_x = length(x); % number of samples in x 
duration_x = len_x * (1/fs); % 1/fs is the time period

% window_duration = 0.1; % frame size (in sec)
% frame_shift = 0.01; % given (in sec)

frame_shift_samples = floor(frame_shift * fs);
frame_size_samples = floor(window_duration * fs);
%num_overlapping_samples = floor(frame_size_samples - frame_shift_samples);
%num_overlapping_samples = max(256, 2^nextpow2(frame_size_samples));

% using default hamming window, default nfft, 
spectrogram(x, frame_size_samples, frame_shift_samples, 'yaxis');
