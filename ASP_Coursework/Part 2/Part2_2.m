%% 2_2_a

clc; clear all; close all;
wo = 0.9; %real weight
std_eta = (0.5).^0.5;
input = std_eta*randn(1500, 1000); 
output = filter([1, wo], [1], input);
output = output(501:end, :); 
input = input(501:end, :);


rho = 0.001; %try changing
alpha = 0.7; %try changing

lr0 = [0.3, 0.1, 0.001]; %different learning rates
gass = {'standard', 'standard', 'ben', 'af', 'mx'};

figure(1);
for k = 1:3 %method
    coefs = zeros(5, 1000, 100); %method, trials, timestep
    lrs = [0.01, 0.1, lr0(k), lr0(k), lr0(k)]; 
    for j = 1:5
        for i = 1:100
            [a, coefs(j, :, i), e] = ...
                lms_arma(output(:, i), input(:, i), 0, 1, lrs(j), gass{j}, rho, alpha, 0);
        end
        
        param_error = -(squeeze(mean(coefs(j, :, :)-wo, 3)));
        subplot(3, 1, k); hold on; set(gca,'fontsize', 15); 
        plot([1:length(param_error)], param_error); ylim([0, 1]);
        hold off
    end
    subplot(3, 1, k); title(sprintf('LMS Weight Error Curves ($\\mu_{0}=%.3f$)',lr0(k)), 'Interpreter', 'Latex');
    lgd=legend('$\mu=0.01$','$\mu=0.1$', 'Benveniste', 'Ang \& Farhang', 'Matthews \& Xie', 'Interpreter', 'Latex'); 
    lgd.FontSize = 10;
    ylabel('w_{0} - w(n)'); xlabel('Time');
    
end

%% 2_2_c
clc; clear all; close all;

wo = 0.9; std_eta = (0.5).^0.5;
input = std_eta*randn(750, 10000);
output = filter([1, wo], [1], input);
output = output(501:end, :); 
input = input(501:end, :);

rho_g = 0.005; rho_b = 0.002;
lr_g = 1; lr_b = 0.1;

figure(1);
params_tot_g = zeros(250, 10000); error_tot_g = zeros(250, 10000);
params_tot_b = zeros(250, 10000); error_tot_b = zeros(250, 10000);
for i = 1:10000
    [~, params_tot_g(:, i), error_tot_g(:, i)] = ...
        gngd(output(:, i), input(:, i), 0, 1, lr_g, rho_g, 0);
    [~, params_tot_b(:, i), error_tot_b(:, i)] = ...
        lms_arma(output(:, i), input(:, i), 0, 1, lr_b,'ben', rho_b, 0, 0);
end

weight_g = (squeeze(mean(params_tot_g, 2)));
weight_b = (squeeze(mean(params_tot_b, 2)));
error_tot_g = squeeze(mean(mag2db(error_tot_g.^2), 2));
error_tot_b = squeeze(mean(mag2db(error_tot_b.^2), 2));

subplot(1,2,1)
plot([1:length(weight_g)], weight_g);hold on
plot([1:length(weight_b)], weight_b);
xlim([0 80]); set(gca,'fontsize', 16); 
title('Weight Estimations');
legend('GNGD', 'Benveniste');
ylabel('Weight value'); xlabel('Time');

subplot(1,2,2);
plot([1:length(error_tot_g)], error_tot_g);hold on
plot([1:length(error_tot_b)], error_tot_b);
xlim([0 80]); set(gca,'fontsize', 16); 
title('Squared Prediction Errors');
legend('GNGD', 'Benveniste');
ylabel('Error (dB)'); xlabel('Time');
hold off;


lr = 0.01;
eta = filter([1, 0, 0.5], [1], randn(1500, 1));
eta = eta(501:end);
x = sin([1:1000]*0.01*pi)';
s = x + eta;
sec_noise = 0.5*eta - 0.1*delayseq(eta, 2);
L = length(s);

MSPE = zeros(2, 2);

[~,xhat] = ale_lms(s, lr, 5, 5);
[~,xhat_anc] = anc_lms(s, sec_noise, lr, 10);


figure(2); subplot(2, 1, 1);
hold on; set(gca,'fontsize', 18);
plot([1:length(xhat)], xhat, 'r');
plot([1:length(xhat_anc)], xhat_anc, 'b');
plot([1:length(x)], x, 'k','LineWidth', 2);
xlabel('Time'); ylabel('Signal'); 
title('Signal Reconstruction with ALE and ANC filters');
lgd = legend('ALE', 'ANC', 'original');
lgd.FontSize = 10;
hold off;

MSPE_ale = 10*log10(movmean((x-xhat).^2, 30));
MSPE_anc = 10*log10(movmean((x-xhat_anc).^2, 30));

subplot(2, 1, 2);
hold on; set(gca,'fontsize', 18);
plot([1:length(MSPE_ale)], MSPE_ale, 'r');
plot([1:length(MSPE_anc)], MSPE_anc, 'b');
xlabel('Time'); ylabel('MSPE (dB)'); 
title('MSPE with ALE and ANC filters');
lgd = legend('ALE', 'ANC');
lgd.FontSize=10;
hold off;

%% functs
function [ar_params, ma_params, error] = lms_arma(output, input, p, q, lr, gass, rho, alpha, leak)

    params = zeros(p+q, length(output)); % n_parameters = order
    phi = zeros(p+q, length(output));
    lrs = lr*ones(size(output)); %adaptive learning rates
    error = ones(size(output));
    
    for i = max([p,q])+1:length(output)-1

        aug_dat = [flip(output(i-p:i-1)), flip(input(i-q:i-1))]';
        error(i) = output(i) - dot(aug_dat, params(:, i)); %signal error
        params(:, i+1) = (1-leak*lr)*params(:, i) + lrs(i)*(error(i))*aug_dat;
       
        if strcmp(gass, 'af')
            if i > max([p,q])+1
                phi(:, i) = alpha*phi(:, i-1) + error(i-1)*prev_aug_dat;
            end
            
        elseif strcmp(gass, 'ben')
            if i > max([p,q])+1
                phi(:, i) = (eye(length(prev_aug_dat))-lrs(i-1)*(prev_aug_dat(:)*prev_aug_dat(:).'))...
                    *phi(:, i-1) + error(i-1)*prev_aug_dat;
            end
        elseif strcmp(gass, 'mx')
            if i > max([p,q])+1
                phi(:, i) = error(i-1)*prev_aug_dat;
            end
            lrs(i+1) = lrs(i) + rho*error(i)*dot(aug_dat, phi(:, i));
        elseif strcmp(gass, 'standard')
        else
            error('wrong input parameters');
        end
        prev_aug_dat = aug_dat;
    end
    ar_params = params(1:p, :); ma_params = params(p+1:q, :);
end
%
function [params, error] = lms_estimator(data, order, lr)

params = zeros(order, length(data));
error = zeros(size(data)) ;

for i = order+1:length(data)
    current_error = data(i) - dot(flip(data(i-order:i-1)), params(:, i-1));
    error(i) = current_error;
    params(:, i) = params(:, i-1) + lr*(current_error)*flip(data(i-order:i-1));
end
end


function [params, error] = leaky_lms(data, order, lr, gamma)
params = zeros(order, length(data));
error = zeros(size(data)) ;

for i = order+1:length(data)
    current_error = data(i) - dot(flip(data(i-order:i-1)), params(:, i-1));
    error(i) = current_error;
    params(:, i) = (1-lr*gamma)*params(:, i-1) + lr*(current_error)*flip(data(i-order:i-1));
end
end

function [ar_params, ma_params, error] = gngd(output, input, p, q, lr, rho, leak)
    params = zeros(p+q, length(output));
    error = ones(size(output));
    reg = ones(size(output))/lr;
    
    for i = max([p,q])+1:length(output)-1

        aug_dat = [flip(output(i-p:i-1)), flip(input(i-q:i-1))]';
        error(i) = output(i) - dot(aug_dat, params(:, i)); 

        lr_now = lr/(reg(i) + dot(aug_dat, aug_dat));
        params(:, i+1) = (1-leak*lr)*params(:, i) + lr_now*(error(i))*aug_dat;

        if i >  max([p,q])+1
            num = rho*lr*error(i)*error(i-1)*dot(old_dat, aug_dat);
            den = ( reg(i-1) + dot(old_dat, old_dat) ).^2;
            reg(i+1) = reg(i) - num/den;
        end
        old_dat = aug_dat;
        
    end
    ar_params = params(1:p, :); ma_params = params(p+1:q, :);
end

function [w, xhat, error] = ale_lms(signal, lr, delta, M)
    w = zeros(M, length(signal)+1); 
    error = zeros(size(signal));
    xhat = signal(size(signal));
    
    for n = delta+M:length(signal)

        u = flip(signal(n-delta-M+1:n-delta));
        xhat(n) = dot(w(:, n), u);
        error(n) = signal(n)  - xhat(n);
        w(:, n+1) = w(:, n) + lr*error(n)*u;

    end
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
