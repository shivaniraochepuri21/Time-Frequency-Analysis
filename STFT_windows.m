function tfa2(window)
% x is the sample amplitude of the signal
% fs is the sampling frequency of the signal
% spectrogram is the magnitude of STFT of the signal
% for fft purpose, it better to have window length in powers of 2

[x, fs] = audioread('tfa_assg3.wav');
x = x(:, 1); % first channel

% sound(x,fs)
% figure(1);
% plot(x);

Ts = 1/fs;
len_x = length(x); % number of samples in x 
duration_x = len_x * (1/fs); % 1/fs is the time period

window_duration = 20e-3; % frame size (in sec) - fixed
frame_shift = 10e-3; % given (in sec) - fixed

frame_shift_samples = frame_shift * fs;
frame_size_samples = window_duration * fs;
num_overlapping_samples = floor(frame_size_samples - frame_shift_samples);

% Using default nfft

% window is 1 implies use hamming window 
if (window == 1)
spectrogram(x,hamming(frame_size_samples),frame_shift_samples,[],'yaxis');
end

% window is 2 implies use hanning window 
if (window == 2)
spectrogram(x,hann(frame_size_samples),frame_shift_samples,[],'yaxis');
end

% window is 3 implies use rectangular window 
if (window == 3)
spectrogram(x,rectwin(frame_size_samples),frame_shift_samples,[],'yaxis');
end

% window is 4 implies use hamming window 
if (window == 4)
spectrogram(x,triang(frame_size_samples),frame_shift_samples,[],'yaxis');
end

end
