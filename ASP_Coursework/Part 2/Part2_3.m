clc; clear variables; close all;
N = 1000;
realisations = 100;

% filter
b = [1 0 0.5];
a = 1;
order = length(b);

% params
mu = 0.01;
delta = 1:25;
M = [5 10 15 20];


MSPE = zeros(length(M),length(delta));
for m=1:length(M)
    xhat = zeros(length(delta),N,realisations);
    error = zeros(length(delta),N,realisations);
    w_anc = zeros(length(delta),N,realisations);
    xhat_anc = zeros(length(delta),N,realisations);
    A = zeros(length(delta),M(m),N,realisations);
    s = zeros(realisations,N);
    for d=1:length(delta)

        x = sin(0.01*pi*(1:N));
        for i=1:realisations
            v = normrnd(0, sqrt(1), 1, N);
            eta = filter(b, a, v);
            s(i,:) = x + eta;

            [xhat(d,:,i), error(d, :, i), A(d, :, :, i)] = lmsALE(s(i,:), mu, M(m), delta(d));
            %[w_anc(d,:,i), xhat_anc(d,:,i)] = anc_lms(s(i,:), v, mu, M(m));
        end
        MSPE(m,d) = mean(mean((repmat(x', 1, realisations)-squeeze(xhat(d,:,:))).^2));
    end
end



figure(1);subplot(1,2,1);hold on;
plot(MSPE','linewidth',2)
set(gca, 'Fontsize', 16)
legend('M=5','M=10','M=15','M=20')
xlabel('Delay ($\Delta$)', 'Interpreter', 'Latex');xlim([3, 25]);
ylabel('MSPE');
title('MSPE vs. Delay, varying filter orders', 'Interpreter','Latex');

subplot(1,2,2);
plot(M,MSPE(:,3),'linewidth',2)
xlabel('Filter Order (M)', 'Interpreter', 'Latex');
ylabel('MSPE');
set(gca, 'Fontsize', 16)
title('MSPE vs. Filter Order ($\Delta=3$)', 'Interpreter', 'Latex');


%% 2_3_d
clc; clear variables; close all;
load EEG_Data_Assignment2;

data = detrend(POz);

t = [0:length(data)-1];
lr = [0.01, 0.001, 0.0001];
lr_opt = 0.001;
M = [1, 5, 10, 20];
stdev = 2;

    mains = sin((2*pi*50/fs)*t) + (10^(-stdev))*randn(1, length(data));
    mains = mains';
    counter = 1;
    figure(1);
    sprintf('Variance: 1e-%d',stdev)

        for k = 1:length(M)
            [w,xhat] = anc_lms(POz, mains, lr_opt, M(k));

            subplot(2, 2, counter)
            ylim([0, 55]); xlim([1, 4.8]); hold on; set(gca,'fontsize', 16);
            spectrogram(xhat, rectwin(3540), round(0.3*(3540)), 16382, fs, 'yaxis');
            title(sprintf('M=%d ', M(k)),'Interpreter','Latex');
            hold off;
            counter = counter + 1;
        end

%% functs

function [xhat, error, A] = lmsALE(x, mu, M, delta)
    N = length(x);
    A = zeros(M,N);
    xhat = zeros(size(x));
    error = zeros(1,N);
    for i=delta+M:N
        xpast = x(i-delta:-1:i-delta-M+1);
        xhat(i) = A(:,i)'*xpast';
        error(i) = x(i)-xhat(i);
        A(:,i+1) = A(:,i) + mu*error(i)*xpast';
    end
    A = A(:,2:end);
end



function [w, xhat] = anc_lms(signal, sec_noise, lr, M)

    w = zeros(M, length(signal)); 
    eta = zeros(size(signal));
    xhat = zeros(size(signal));
    u = delayseq(repmat(sec_noise, 1, M), [0:M-1])';
 
    for n = 1:length(signal)
        eta(n) = dot(w(:, n), u(:, n));
        xhat(n) = signal(n)  - eta(n);
        w(:, n+1) = w(:, n) + lr*xhat(n)*u(:, n);
    end

end