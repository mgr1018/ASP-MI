%% 1_1_a
clear all;
close all;
clc;

Fs = 1000; 
T = 1/Fs; 
L = 3000; 
t = (0:L-1)*T; 


S1 = zeros(size(t));
S1(1500) = 2;
Rxx1 = xcorr(S1, 'biased');
f1 = -Fs/2:Fs/length(Rxx1):Fs/2 - Fs/length(Rxx1);
P_rxx1 = abs(fftshift(fft(Rxx1)));
periodogram = (abs(fftshift(fft(S1))).^2)/L;
f2 = -Fs/2:Fs/length(periodogram):Fs/2 - Fs/length(periodogram);

S2 = 0.7 * sin(2*pi*0.5*t);
Rxx2 = xcorr(S2, 'biased');
P_rxx2 = abs(fftshift(fft(Rxx2)));
periodogram2 = (abs(fftshift(fft(S2))).^2)/L;

figure;
subplot(3, 1, 1);
plot([-L+1:L-1], Rxx2, 'LineWidth', 2);
hold on;
plot([-L+1:L-1], Rxx1, 'LineWidth', 2);
legend('Sinusoids', 'Impulse');
title('ACF of impulse and sinusoid');
xlabel('Lag');
ylabel('ACF');

subplot(3, 1, 2);
plot(f1, P_rxx2, 'LineWidth', 2);
hold on;
plot(f2, periodogram2,'LineWidth', 2);
legend('Definition 1', 'Definition 2');
title('PSD sinusoid');
xlabel('Normalised frequency (\pi rad/sample)');
ylabel('PSD');

subplot(3, 1, 3);
plot(f1, P_rxx1, 'LineWidth', 2);
hold on;
plot(f2, periodogram, '--', 'LineWidth', 2);
legend('Definition 1', 'Definition 2');
title('PSD impulse');
xlabel('Normalised frequency (\pi rad/sample)');
ylabel('PSD');

%% Part 1_2_a
close all;
clear all;
clc;
load sunspot.dat 

data = sunspot(:,2);
[periodogram_sunspot,w] = periodogram(data);
w = w/pi;
data_processed = detrend(data - mean(data));
periodogram_sunspot_processed = periodogram(data_processed);
log_data = log10(data);
log_data(isinf(log_data)) = 0;
log_data_1_processed = log_data - mean(log_data);
log_periodogram_1 = periodogram(log_data_1_processed);

figure
plot(w,10*log10(periodogram_sunspot));hold on
plot(w,10*log10(periodogram_sunspot_processed));hold on
plot(w,10*log10(log_periodogram_1));
ylim([-100,100])
title('sunspot series periodogram');
xlabel('normalized frequency (rad)');
ylabel('PSD');
legend('raw','centered, detrended','logarithmic data centered, detrended')

%% Part 1_2_b
close all;
clear all;
clc;

load EEG_Data\EEG_Data\EEG_Data_Assignment1.mat

L = length(POz);
T=1/fs;

pxx_10 = pwelch(POz,10*fs,0,L,fs);
pxx_5 = pwelch(POz,5*fs,0,L,fs);
pxx_1 = pwelch(POz,fs,0,L,fs);
[pxx,f] = pwelch(POz,L,0,L,fs);


figure
subplot(2,2,1)
plot(f,10*log10(pxx))
title('PSD');
xlabel('frequency');
ylabel('power');
xlim([1, 50])

subplot(2,2,2)
plot(f,10*log10(pxx_5))
title('PSD 5s window')
xlabel('freq (Hz)')
ylabel('Power')
xlim([1, 50])

subplot(2,2,3)
plot(f,10*log10(pxx_1))
title('PSD 1s window')
xlabel('freq (Hz)')
ylabel('Power')
xlim([1, 50])

subplot(2,2,4)
plot(f,10*log10(pxx_10))
title('PSD 10s window')
xlabel('freq (Hz)')
ylabel('Power')
xlim([1, 50])


%% Part 1_3
close all;
clear all;
clc;

t = [0:0.01:99.99]; 
fs = 10; 

noisy_sinusoid = 3*sin(2*5*pi*t) + randn(size(t)); 
noise = randn(size(t));
AR = arima('Constant', 0,'AR',{2.21, -2.94, 2.17, -0.96},'Variance',1);
filtered_noise = simulate(AR, length(t));

[PSD1, f] = correlogram(noisy_sinusoid, 'biased', fs);
[PSD1_unbiased, f] = correlogram(noisy_sinusoid, 'unbiased', fs);

[PSD2, f] = correlogram(noise, 'biased', fs);
[PSD2_unbiased, f] = correlogram(noise, 'unbiased', fs);

[PSD3, f] = correlogram(filtered_noise, 'biased', fs);
[PSD3_unbiased, f] = correlogram(filtered_noise, 'unbiased', fs);

f = 2*f/fs; 

figure(1);
subplot(2, 3, 1); plot(t, ACF(noisy_sinusoid, 'biased')); 
xlabel('Time'); 
ylabel('ACF'); title('ACF of Noisy Sinusoid (Biased)');

subplot(2, 3, 2); plot(t, ACF(noise, 'biased')); 
xlabel('Time (s)');
ylabel('ACF'); title('ACF of WGN (Biased)');

subplot(2, 3, 3); plot(t, ACF(filtered_noise, 'biased')); 
xlabel('Time');
ylabel('ACF'); title('ACF of AR(4) Process (Biased)');

subplot(2, 3, 4); plot(t, ACF(noisy_sinusoid, 'unbiased')); 
xlabel('Time'); 
ylabel('ACF'); title('ACF of Noisy Sinusoid (Unbiased)');

subplot(2, 3, 5); plot(t, ACF(noise, 'unbiased')); 
xlabel('Time');
ylabel('ACF'); title('ACF of WGN (Unbiased)');

subplot(2, 3, 6); plot(t, ACF(filtered_noise, 'unbiased')); 
xlabel('Time (s)');
ylabel('ACF'); title('ACF of AR(4) Process (Unbiased)');


figure(2);
subplot(3, 2, 1); plot(f, 10*log10(PSD1)); 
set(gca,'fontsize', 14); xlabel('Normalized Frequency'); 
ylabel('PSD(dB)'); title('Biased Noisy Sinusoid PSD');xlim([0,0.8])

subplot(3, 2, 3); plot(f, 10*log10(PSD2)); 
set(gca,'fontsize', 14); xlabel('Normalized Frequency');
ylabel('PSD(dB)'); title('Biased WGN PSD');xlim([0,0.8])

subplot(3, 2, 5); plot(f, 10*log10(PSD3)); 
set(gca,'fontsize', 14); xlabel('Normalized Frequency');
ylabel('PSD(dB)'); title('Biased Filtered WGN PSD');xlim([0,0.8])


subplot(3, 2, 2); plot(f, 10*log10(PSD1_unbiased)); 
set(gca,'fontsize', 14); xlabel('Normalized Frequency'); 
ylabel('PSD(dB)'); title('Unbiased Noisy Sinusoid PSD');xlim([0,0.8])

subplot(3, 2, 4); plot(f, 10*log10(PSD2_unbiased)); 
set(gca,'fontsize', 14); xlabel('Normalized Frequency');
ylabel('PSD(dB)'); title('Unbiased WGN PSD');xlim([0,0.8])

subplot(3, 2, 6); plot(f, 10*log10(PSD3_unbiased)); 
set(gca,'fontsize', 14); xlabel('Normalized Frequency');
ylabel('PSD(dB)'); title('Unbiased Filtered WGN PSD');

xlim([0,0.8])

%% Part 1_3_b/c 
clc; clear all; close all;

t = [0:0.01:9.99]; fs = 10; PSD_vect = [];

noisy_sinusoid = repmat(sin(2*pi*30*t), 1000, 1)' + 5*randn(1000, 1000);
figure(1);
for i = 1:100
    [PSD, f] = correlogram(noisy_sinusoid(:, i), 'biased', fs);
    f = 2*f/fs;
    subplot(2, 2, 1); hold on;
    plot(f, PSD, 'b'); hold off;

    subplot(2, 2, 3); hold on;
    plot(f, 10*log10(PSD), 'b'); hold off;
    PSD_vect = [PSD_vect; PSD]; 
end

subplot(2, 2, 1); hold on; plot(f, mean(PSD_vect, 1), 'r'); set(gca,'fontsize', 12);
xlabel('Normalized Frequency'); ylabel('PSD');
title('100 PSD estimates of a noisy sinusoid with their mean'); hold off;

subplot(2, 2, 2); hold on; plot(f, std(PSD_vect, 1)); set(gca,'fontsize', 12);
xlabel('Normalized Frequency'); ylabel('PSD');
title('Standard deviation of PSD estimate'); hold off;

subplot(2, 2, 3); hold on; plot(f, 10*log10(mean(PSD_vect, 1)), 'r'); set(gca,'fontsize', 12);
xlabel('Normalized Frequency'); ylabel('PSD(dB)');

title('100 PSD estimates of a noisy sinusoid (dB) with their mean');hold off;
subplot(2, 2, 4); hold on; plot(f, 10*log10(std(PSD_vect, 1))); set(gca,'fontsize', 12);
xlabel('Normalized Frequency'); ylabel('PSD(dB)');
title('Standard deviation of PSD estimate (dB)'); hold off;

%% part 1_3_d
clc; clear all; close all;
t = [0:0.5:50];
pts = [10, 50, 100];
PSDd = zeros(3,length(t));
x = zeros(1,length(t));

for i=1:length(pts)
    
    n = 1:pts(i);
    noise = 0.2*(randn(size(n))+1j*randn(size(n)));
    x(1:pts(i)) = 0.3* exp(1j*2*pi*0.2*n)+ 0.5*exp(1j*2*0.25*pi*n)+ noise;
    PSDd(i,:) = fft(x);  
end

normf = linspace(0, 1, length(PSDd(1,:)));
figure
plot(normf,abs(PSDd(1,:)))
hold on
plot(normf,abs(PSDd(2,:)))
plot(normf,abs(PSDd(3,:)))
xlim([0 0.6])
xlabel('Normalised Frequency')
ylabel('Magnitude')
set(gca, 'Fontsize', 15)
legend('n = 10','n = 50','n = 100')
title('PSD of complex exponential', 'Fontsize', 14)

%% part 1_3_e
clc; clear all; close all;

n = 0:30;
noise = 0.1*(randn(size(n))+1j*randn(size(n)));
x = 0.5*exp(1j*2*pi*0.5*n)+ 0.4*exp(1j*2*0.52*pi*n)+ noise;

tic
[X,R] = corrmtx(x, 14, 'mod'); 
[S,F] = pmusic(R, 2, [ ], 1, 'corr');
toc

tic
[Pxx, w] = periodogram(x);
toc 

figure(1); subplot(1, 2, 1);
plot(2*F,S); set(gca,'fontsize', 14);xlim([0.7,1.3]);
xlabel('Frequency (Hz)'); 
ylabel('Pseudospectrum'); title('MUSIC Pseudospectrum')
subplot(1, 2, 2); plot(w/pi, Pxx); set(gca,'fontsize', 14);xlim([0.7,1.3])
ylabel('PSD');
xlabel('Frequency (Hz)'); title('Standard PSD Estimate')

%% functions


function sol = ACF(x, biased)
    sol = []; N = length(x);
    for k=0:N-1
        if strcmp(biased, 'biased')
            a = (1/N);
        else
            a = (1/(N-k)); 
        end
        sol(k+1) = a * dot(x(k+1:N), conj(x(1:N-k)));
    end
end

function [PSD, f] = correlogram(s, biased, fs) 
    PSD = abs(fft(ACF(s, biased))); L = length(s);
    PSD = PSD(1:floor(L/2)+1);
    f = (fs/2)*[0:length(PSD)- 1]/length(PSD);
end