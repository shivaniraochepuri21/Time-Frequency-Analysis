% window_duratoin is the frame size, all in seconds

clc;close all;clear all;
window_duration1 = 100e-3; frame_shift1 = 10e-3;
tfa1(window_duration1, frame_shift1);

%%
clc;close all;clear all;
window_duration2 = 50e-3; frame_shift2 = 10e-3;
tfa1(window_duration2, frame_shift2);
%%
clc;close all;clear all;
window_duration3 = 20e-3; frame_shift3 = 10e-3;
tfa1(window_duration3, frame_shift3);

%%
clc;close all;clear all;
window_duration4 = 10e-3; frame_shift4 = 5e-3;
tfa1(window_duration4, frame_shift4);

%%
clc;close all;clear all;
tfa2(1); % hamming window
%%
clc;close all;clear all;
tfa2(2); % hanning window
%%
clc;close all;clear all;
tfa2(3); % rectandular window
%%
clc;close all;clear all;
tfa2(4); % triangular window

