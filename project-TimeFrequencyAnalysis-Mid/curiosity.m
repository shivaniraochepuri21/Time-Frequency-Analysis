clc;
clear all; 
close all;
%%
[x,Fs] = audioread("arctic_a0001.wav");
y = downsample(x,20); 
T = length(y);
t = 1:T;

% sound(y,Fs)
 figure(23);
 plot(y);
% title("Speech Signal");

alpha = 80; %the balancing parameter of the data-fidelity constraint
tau = 0; %time-step of the dual ascent ( pick 0 for noise-slack ) 
tol = 1e-20; %tolerance of convergence criterion; typically around 1e-6 %
K = 2; %the number of modes to be recovered %NumIMF
DC = 0; %true if the first mode is put and kept at DC (0-freq)
init = 0; % 'InitializeMethod', init

%[imf,residual,info] = vmd(y,'NumIMF',K, 'AbsoluteTolerance', tol, 'PenaltyFactor', alpha, 'LMUpdateRate' ,tau);

i = 1;
comp = 1;
y_iter = y;
CF_arr = 0;

diff = 0.001;
y_F0 = zeros(length(y),1);
y_CF = 0;

while(1)
    
    %[imf,residual,info] = vmd(y_iter,'NumIMF',K);
    %[u,u_hat,omega] = VMD(y_iter, alpha, tau, K, DC, init, tol);
    %if info.CentralFrequencies(1) >= info.CentralFrequencies(2)
    [v1,v2,v3] = vmd(y_iter,'NumIMF',K)

    figure(i);
    subplot(211);
    plot(t,v1(:,1));
    subplot(212);
    plot(t,v1(:,2));

    if v3.CentralFrequencies(1) <= v3.CentralFrequencies(2) 
        comp = 1;
    else
        comp = 2;
    end
    
    if comp == 1
        y_iter = v1(:,1);
    else
        y_iter = v1(:,2);
    end
    
    if i == 1
        CF_arr = v3.CentralFrequencies(comp);
    else
        CF_arr = [CF_arr; v3.CentralFrequencies(comp)];
    end
    
    if (i >= 2)
        if abs(CF_arr(i) - CF_arr(i-1)) <= diff
            y_F0 = y_iter;
            y_CF = CF_arr(i);
            break;
        end
        
        if(i>100)
            y_F0 = y_iter;
            y_CF = CF_arr(i);
            break;
        end
        
    end
        
    i = i + 1;
end

t_i = 1:i;
% figure();
% plot(t,y_F0);
% figure();
% plot(t_i, CF_arr)


e = envelope(y_F0);
figure(123);
plot(t,e);

yh = hilbert(y_F0);
ya = y_F0 + 1j*yh;
figure(223);
plot(t,yh);

figure(323);
plot(t,ya);

figure(4);
plot(t,y);
title("Speech Signal");