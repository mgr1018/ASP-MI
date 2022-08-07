%% 3_3_a
clc; clear variables; close all;

n_samples = 1500;
f = [100*ones(500,1); ...
    100*ones(500,1) + ([501:1000]' - 500)/2; ...
    100*ones(500,1) + (([1001:1500]' -  1000)/25).^2];
phi = cumsum(f);
fs = 1000;
variance = 0.05;
eta = wgn(n_samples,1,pow2db(variance),'complex');
y = exp(1i*2*pi*phi/fs) + eta;

time = 1:n_samples;
f_ax = 0:fs-1;
x = (1/fs)*exp(1i*2*pi*f_ax'*time./fs);

leaks = [0,0.03, 0.2, 0.5];
figure(1);
for idx = 1:4
    [w,~] = dft_clms(y.', x, leaks(idx));

    subplot(2,2,idx); hold on; set(gca,'fontsize', 16);
    mesh(time, f_ax(1:floor(fs/2)), abs(w(1:floor(fs/2),:)).^2);
    view(2);
    xlabel('Time (n)');
    ylabel('Frequency (Hz)');
    ylim([0,500]);
    title(strcat('DFT-CLMS Frequency Estimate ($\gamma=',num2str(leaks(idx)),'$)'),'Interpreter','Latex');
    hold off;
end

%% 3_3_d
clc; clear variables; close all;

load('EEG_Data_Assignment2.mat');

a = 1;
t_range = a:a+1199;
POz = detrend(POz(t_range));
n_samples = length(t_range);
f_ax = 0:fs-1;
x = (1/fs)*exp(1i*2*pi*f_ax'*t_range./fs);
[w,~] = dft_clms(POz, x, 0);

figure(1); hold on; set(gca,'fontsize', 16);
mesh(t_range, f_ax(1:floor(fs/2)), abs(w(1:floor(fs/2),:)).^2);
view(2);ylim([0,80]);
xlabel('Time (n)');
ylabel('Frequency (Hz)');
title('DFT-CLMS spectral frequency Estimate of POz');

%% fcts

function [h, error] = dft_clms(y, x, leak)

    [nft, n_samples] = size(x);
    h = zeros(nft, n_samples,'like',1i); 
    error = zeros(n_samples,1,'like',1i);
    lr = 1;
    for n = 1:n_samples
        error(n) = y(n) - h(:,n)' * x(:,n);
        if n < length(x)
            h(:, n+1) = (1-lr*leak)*h(:, n) + lr*conj(error(n))*x(:,n);
        end
    end
end

