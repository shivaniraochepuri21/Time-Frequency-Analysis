clc;
clear all; 
close all;
%%
%speech signal
[x,Fs] = audioread("arctic_a0005.wav");
fs = 8000;
ya = resample(x(:,1),fs,Fs);
% normalising the speech signal
ya = ya/abs(max(ya));
T = length(ya);
t = 1:T;

% DEGG - Derivative of electroglottographic signal
yb = resample(x(:,2),fs,Fs);
yb = yb/abs(max(yb));

%%
alpha = 80; %the balancing parameter of the data-fidelity constraint
tau = 0; %time-step of the dual ascent ( pick 0 for noise-slack ) 
tol = 1e-20; %tolerance of convergence criterion; typically around 1e-6 %
K = 2; %the number of modes to be recovered %NumIMF
DC = 0; %true if the first mode is put and kept at DC (0-freq)
init = 0; % 'InitializeMethod', init

%[imf,residual,info] = vmd(y,'NumIMF',K, 'AbsoluteTolerance', tol, 'PenaltyFactor', alpha, 'LMUpdateRate' ,tau);

%%
i = 1;
comp = 1;
y_iter = ya;
CF_arr = 0;

difference = 1e-6;
y_F0 = zeros(length(ya),1);
y_CF = 0;

while(1)
    
    [v1,v2,v3] = vmd(y_iter,'NumIMF',K);
    
    %selecting which vmf to use for the next iteration
    if v3.CentralFrequencies(1) <= v3.CentralFrequencies(2) 
        comp = 1;
    else
        comp = 2;
    end
    % choosing the VMF
    if comp == 1
        y_iter = v1(:,1);
    else
        y_iter = v1(:,2);
    end
    % appending the central freq of the vmfs we chose
    % to plot the central freq convergence graph
    if i == 1
        CF_arr = v3.CentralFrequencies(comp);
    else
        CF_arr = [CF_arr; v3.CentralFrequencies(comp)];
    end
    % comparing previous and current iteraiton's central freq
    if (i >= 2)
        if abs(CF_arr(i) - CF_arr(i-1)) <= difference
            y_F0 = y_iter;
            y_CF = CF_arr(i);
            break;
        end
        
        if(i>500)
            y_F0 = y_iter;
            y_CF = CF_arr(i);
            break;
        end
        
    end
        
    i = i + 1;
end

e = envelope(y_F0);

%thresholding - empirical for V/UV detection - given in paper
th = 0.2*max(e);

voiced_flag = zeros(1,length(e));
voiced = zeros(1,length(e));

for k=1:length(e)
    if e(k)>=th
        voiced(k)=e(k);
        voiced_flag(k)=0.01;
    end
end
% hilbert transform and analytical signal
yh = hilbert(voiced);
z = voiced + 1j*yh;

% finding phase of the obtained analytical signal
phi=zeros(1,length(z));
for j=1:length(z)
    phi(j)= atan(imag(z(j))/real(z(j)));
end

% omega=diff(phi)/(1/fs);

% differentiation of phase to get IF
omega = zeros(1,length(phi));
for i=1:length(phi)-1
    omega(i) = (phi(i+1) - phi(i))/(t(i+1) - t(i));
end
omega(length(phi)) = omega(length(phi)-1);

%% IF extraction using EMD (HHT)
% comment this section while running VMD based method
% uncomment this section to see EMD based results
% 
% [e1,e2,e3] = emd(ya);
% emd(ya)
% y_F0 = e1(:,9);
% e = envelope(y_F0);
% yh = hilbert(y_F0);
% z = y_F0 + 1j*yh;
% 
% phi = zeros(1,length(z));
% for j=1:length(z)
%     phi(j)= atan(imag(z(j))/real(z(j)));
% end
% 
% % differentiation of phase to get IF
% omega = diff(phi);
% 
% % hht
% hht(y_F0,fs)
% [hs,f,t,imfinsf,imfinse] = hht(y_F0);
% %figure('Name','Hilbert Spectrum');
% %plot(hs);

%% matlab function
ifq = instfreq(y_F0,fs);
% pitch or F0 estimation
f0 = pitch(ya,fs);

%% plotting signals

figure(1);
plot(ya);
grid on
title("Speech Signal");

figure(2);
plot(yb);
grid on
title("DEGG signal");

% histogram of the speech signal
figure(3);
histogram(ya, 'FaceColor', 'blue')

% 2d histogram of the speech signal
figure(4);
histogram2(t', ya, 'FaceColor', 'blue')
%

figure(5);
subplot(211)
plot(y_F0)
grid on
title("Fundamental Frequency component")
subplot(212)
plot(e);
grid on
title("Envelope of F0")

figure(6);
subplot(2,1,1)
plot(voiced)
hold on
plot(voiced_flag)
title("thresholding done")
subplot(2,1,2)
plot(e);
title("before thresholding- ie envelop of F0")

figure(7);
subplot(211)
plot(yh);
title("Hilbert Transform of F0");
subplot(212)
plot(z)
title("Analytical signal of F0");

t_v = 1:length(omega);
figure(8);
subplot(211)
plot(phi)
title("phase signal of F0")
subplot(212)
%scatter(t_v,omega)
stairs(omega)
title("Instantaneous Frequency plot")

figure(9);
periodogram(phi,[],[],Fs)

figure(10);
plot(ifq);
title("Inst Freq");

figure(11);
plot(f0);
title("F0 or pitch");

% uncomment this following lines to run vmd
% to get the plot for central frequency convergence
% figure(13);
% plot(CF_arr)

