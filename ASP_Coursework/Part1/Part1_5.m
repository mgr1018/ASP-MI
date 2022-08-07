%% 1_5_a

close all; clear all; clc;
load ECG_data\Data\RAW-ECG.mat;
load ECG_data\Data\RRI-DATA.mat;

xRRI_1 = detrend(normalize(xRRI1));
xRRI_2 = detrend(normalize(xRRI2));
xRRI_3 = detrend(normalize(xRRI3));

RRI_data = {xRRI_1; xRRI_2; xRRI_3};
leg = {'Standard', '(W = 100s)', '(W = 500s)'};

figure(1);
for j = 1:3 
    subplot(3, 1, j); hold on;
    L = length(RRI_data{j});
    wl = [length(RRI_data{j}), 100, 500];
    for i = 1:3 % Window sizes
        [P_mean, w] = pwelch(RRI_data{j}, wl(i), 0, length(RRI_data{j}), 4);
        plot(w, 10*log10(P_mean),"LineWidth",1);
    end
    set(gca,'fontsize', 14);
    xlabel('Frequency (Hz)');
    ylabel('PSD(dB)'); 
    title(sprintf('RRI PSD Estimate (Trial %d)', j));
    legend(leg); hold off;
end

%% Part c
clc; clear all; close all;

load ECG_data\Data\RAW-ECG.mat;
load ECG_data\Data\RRI-DATA.mat;

xRRI_1 = detrend(normalize(xRRI1));
xRRI_2 = detrend(normalize(xRRI2));
xRRI_3 = detrend(normalize(xRRI3));


% 'Determining model order' code
leg = {};
figure(1); hold on;
for i = [1:10]
    [pxx, w] = pyulear(xRRI_2, i, 2048); % power spectrum estimate given AR model
    plot(w/pi, 10*log10(pxx));
    leg{end+1} = sprintf('Model Order %d', i);
end
xlabel('Normalized Frequency (\times \pi rad/sample)');
ylabel('Power/frequency (dB/(rad/sample))'); 
legend(leg); hold off;

% Trial 1 p=2
% Trial 2 p=10
% Trial 3 p=3

%%
% Plotting results
order = [2, 10, 3]; RRI_data = {xRRI_1; xRRI_2; xRRI_3};
figure(1); hold on;
for i = [1, 2, 3]
    [pxx, w] = pyulear(RRI_data{i}, order(i), 2048, 4); % power spectrum estimate given AR model
    plot(w, 10*log10(pxx));
end
set(gca,'fontsize', 12);
xlabel('Frequency (Hz)');
ylabel('PSD(dB)'); 
title('AR spectrum Estimate of RRI Data');
legend('Trial 1', 'Trial 2', 'Trial 3'); hold off;

