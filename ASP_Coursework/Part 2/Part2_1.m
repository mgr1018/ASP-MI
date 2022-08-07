
%% part 2_1
clc; clear all; close all;

std_eta = 0.5;
x = filter([1], [1, -0.1, -0.8], std_eta*randn(1500, 100));
x = x(501:end, :);

[X1, error1] = lms_estimator(x(:, 1), 2, 0.05);
[X2, error2] = lms_estimator(x(:, 1), 2, 0.01);
figure(1); subplot(1, 2, 1); hold on; set(gca,'fontsize', 14);
plot([1:length(error1)], 10*log10(error1.^2));
plot([1:length(error2)], 10*log10(error2.^2));
xlabel('Time Step'); ylabel('Squared Prediction Error (dB)');
title('Squared Prediction Error of signal x(n)');
legend('$\mu=0.05$','$\mu=0.01$', 'Interpreter', 'Latex');
hold off;


error = zeros(2, 100, 1000);
lrs = [0.05, 0.01];

for i = 1:100
    for j = 1:2
        [a, error(j, i, :)] = lms_estimator(x(:, i), 2, lrs(j)); 
    end
end

error = squeeze(mean(10*log10(error.^2), 2)); 

figure(1); 
subplot(1, 2, 2); hold on; 
set(gca,'fontsize', 14);
plot([1:length(error(1, :))], error(1, :));
plot([1:length(error(2, :))], error(2, :));
xlabel('Time Step'); ylabel('Squared Prediction Error (dB)');
title('Learning Curve (100 realisations)');
legend('$\mu=0.05$','$\mu=0.01$', 'Interpreter', 'Latex');
hold off;

%% part 2_1_c
clc; clear all; close all;

std_eta = 0.5;
x = filter([1], [1, -0.1, -0.8], std_eta*randn(200500, 100));
x = x(501:end, :);

error = zeros(2, 100, 200000); lrs = [0.05, 0.01];
for i = 1:100
    for j = 1:2
        [a, error(j, i, :)] = lms_estimator(x(:, i), 2, lrs(j)); 
    end
end
error = squeeze(mean(error(:, 100, 1000:end).^2, [3, 2]));
M1 = (error(1) - std_eta.^2)/std_eta.^2
M2 = (error(2) - std_eta.^2)/std_eta.^2

%% part 2_1_d
clc; clear all; close all;

a1 = 0.1; a2 = 0.8;
std_eta = 0.5;
x = filter([1], [1, -0.1, -0.8], std_eta*randn(2000, 100));
x = x(501:end, :);

coefs = zeros(2, 2, 100, 1500); lrs = [0.05, 0.01];
for i = 1:100
    for j = 1:2
        [coefs(j, :, i, :), e] = lms_estimator(x(:, i), 2, lrs(j)); 
    end
end
coefs = squeeze(mean(coefs, 3)); 
params_tot1 = squeeze(coefs(1,:,:,:));
a1_est1 = sum(params_tot1(1,:),"all")/1500
a2_est1 = sum(params_tot1(2,:),"all")/1500

params_tot2 = squeeze(coefs(2,:,:,:));
a1_est2 = sum(params_tot2(1,:),"all")/1500
a2_est2 = sum(params_tot2(2,:),"all")/1500

%%
figure(1); subplot(1, 2, 1); hold on; set(gca,'fontsize', 18); ylim([0, 1]);
h(1)=plot([1:length(params_tot1(1, :))], params_tot1(1, :), 'r'); 
h(2)=plot([1:length(params_tot1(2, :))], params_tot1(2, :), 'b');
h(3)=plot([1:length(params_tot1(1, :))], a1*ones(1,length(params_tot1(1, :))), '--k');
h(4)=plot([1:length(params_tot1(2, :))], a2*ones(1,length(params_tot1(2, :))), '--m');
xlabel('Time Step'); ylabel('Parameter Estimate');
title('LMS Coefficient Estimation ($\mu=0.05$)', 'Interpreter', 'Latex');
legend(h,'$a_{1}$', '$a_{2}$', '$\hat{a}_{1}$','$\hat{a}_{2}$', 'Interpreter', 'Latex');
hold off;

subplot(1, 2, 2); hold on; set(gca,'fontsize', 18); ylim([0, 1]);
h(1)=plot([1:length(params_tot2(1, :))], params_tot2(1, :), 'r'); 
h(2)=plot([1:length(params_tot2(2, :))], params_tot2(2, :), 'b');
h(3)=plot([1:length(params_tot2(1, :))], a1*ones(1,length(params_tot2(1, :))), '--k');
h(4)=plot([1:length(params_tot2(2, :))], a2*ones(1,length(params_tot2(2, :))), '--m');
xlabel('Time Step'); ylabel('Parameter Estimate');
title('LMS Coefficient Estimation ($\mu=0.01$)', 'Interpreter', 'Latex');
legend(h,'$a_{1}$', '$a_{2}$', '$\hat{a}_{1}$','$\hat{a}_{2}$', 'Interpreter', 'Latex');
hold off;

%% 2_1_f
clc; clear all; close all;

a1 = 0.1; a2 = 0.8; std_eta = 0.5;
x = filter([1], [1, -0.1, -0.8], std_eta*randn(1500, 100));
x = x(501:end, :); 

lrs = [0.05, 0.01]; 
gammas = [0.2, 0.4, 0.8];
coefs = zeros(2, 1000, 100); %lr, timestep, trials
for j = 1:2 %lrs
    figure;
    for k = 1:3 %gamma
        for i = 1:100 %trials
            [coefs(:, :, i), e] = leaky_lms(x(:, i), 2, lrs(j), gammas(k)); 
        end
        params = squeeze(mean(coefs, 3));
        
        if(j==1)
            subplot(1, 2, 1); hold on;
            set(gca,'fontsize', 15); ylim([0, 1]); 
            if (k == 1)
                h(1)=plot([1:length(params(1, :))], a1*ones(1,length(params(1, :))), '--');
            end
                h(2)=plot([1:length(params(1, :))], params(1,:));

            title(sprintf('Leaky LMS Coefficient a1 ($\\mu=%.2f$)',...
                lrs(j)), 'Interpreter', 'Latex');
                    hleg = legend(sprintf('$a_{1}$'),sprintf('$estimation, \\gamma=%.1f$',gammas(1)),sprintf('$estimation, \\gamma=%.1f$',gammas(2)),sprintf('$estimation, \\gamma=%.1f$',gammas(3)),'Interpreter', 'Latex');
                    hleg.FontSize = 10;
                    xlabel('Time Step'); ylabel('Coefficient values');
           
    
            subplot(1, 2, 2); hold on;
            set(gca,'fontsize', 15); ylim([0, 1]); 
            if(k==1)
                h(1)=plot([1:length(params(2, :))], a2*ones(1,length(params(2, :))), '--');
            end
                h(2)=plot([1:length(params(2, :))], params(2,:));
    
            title(sprintf('Leaky LMS Coefficient a2 ($\\mu=%.2f$)',...
                lrs(1)), 'Interpreter', 'Latex');
                    hleg = legend(sprintf('$a_{2}$'),sprintf('$estimation, \\gamma=%.1f$',gammas(1)),sprintf('$estimation, \\gamma=%.1f$',gammas(2)),sprintf('$estimation, \\gamma=%.1f$',gammas(3)),'Interpreter', 'Latex');
                    hleg.FontSize = 10;
                    xlabel('Time Step'); ylabel('Coefficient values');

        else
            subplot(1, 2, 1); hold on;
            set(gca,'fontsize', 15); ylim([0, 1]); 
            if(k==1)
                h(1)=plot([1:length(params(1, :))], a1*ones(1,length(params(1, :))), '--');
            end
                h(2)=plot([1:length(params(1, :))], params(1, :));
    
            title(sprintf('Leaky LMS Coefficient a1 ($\\mu=%.2f$)',...
                lrs(2)), 'Interpreter', 'Latex');    
                    hleg = legend(sprintf('$a_{1}$'),sprintf('$estimation, \\gamma=%.1f$',gammas(1)),sprintf('$estimation, \\gamma=%.1f$',gammas(2)),sprintf('$estimation, \\gamma=%.1f$',gammas(3)),'Interpreter', 'Latex');
                    hleg.FontSize = 10;
                    xlabel('Time Step'); ylabel('Coefficient values');

            subplot(1, 2, 2); hold on;
            set(gca,'fontsize', 15); ylim([0, 1]); 
            if(k==1)
                h(1)=plot([1:length(params(2, :))], a2*ones(1,length(params(2, :))), '--');
            end
                h(2)=plot([1:length(params(2, :))], params(2,:));
    
            title(sprintf('Leaky LMS Coefficient a2 ($\\mu=%.2f$)',...
                lrs(2)), 'Interpreter', 'Latex');
                    hleg = legend(sprintf('$a_{2}$'),sprintf('$estimation, \\gamma=%.1f$',gammas(1)),sprintf('$estimation, \\gamma=%.1f$',gammas(2)),sprintf('$estimation, \\gamma=%.1f$',gammas(3)),'Interpreter', 'Latex');
                    hleg.FontSize = 10;
                    xlabel('Time Step'); ylabel('Coefficient values');
        end  
        
    end

   

end
hold off


%% functs

function [ar_params, ma_params, error] = lms_arma(output, input, p, q, lr, gass, rho, alpha, leak)

    params = zeros(p+q, length(output)); % n_parameters = order
    phi = zeros(p+q, length(output));
    lrs = lr*ones(size(output)); %adaptive learning rates
    error = ones(size(output));
    
    for i = max([p,q])+1:length(output)-1
        % Create augmented data vector with both y and x values
        aug_dat = [flip(output(i-p:i-1)), flip(input(i-q:i-1))]';
        error(i) = output(i) - dot(aug_dat, params(:, i)); %signal error
        % SimultaneousAR and MA updates
        params(:, i+1) = (1-leak*lr)*params(:, i) + lrs(i)*(error(i))*aug_dat;
        % Step-size update
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

params = zeros(order, length(data)); % n_parameters = order
error = zeros(size(data)) ;

for i = order+1:length(data)
    current_error = data(i) - dot(flip(data(i-order:i-1)), params(:, i-1));
    error(i) = current_error;
    params(:, i) = params(:, i-1) + lr*(current_error)*flip(data(i-order:i-1));
end
end


function [params, error] = leaky_lms(data, order, lr, gamma)
params = zeros(order, length(data)); % n_parameters = order
error = zeros(size(data)) ;

for i = order+1:length(data)
    current_error = data(i) - dot(flip(data(i-order:i-1)), params(:, i-1));
    error(i) = current_error;
    params(:, i) = (1-lr*gamma)*params(:, i-1) + lr*(current_error)*flip(data(i-order:i-1));
end
end

function [ar_params, ma_params, error] = gngd(output, input, p, q, lr, rho, leak)
    params = zeros(p+q, length(output)); % n_parameters = order
    error = ones(size(output));
    reg = ones(size(output))/lr;
    
    for i = max([p,q])+1:length(output)-1
        % Create augmented data vector with both y and x values
        aug_dat = [flip(output(i-p:i-1)), flip(input(i-q:i-1))]';
        error(i) = output(i) - dot(aug_dat, params(:, i)); %signal error
        % Simultaneous AR and MA updates
        lr_now = lr/(reg(i) + dot(aug_dat, aug_dat));
        params(:, i+1) = (1-leak*lr)*params(:, i) + lr_now*(error(i))*aug_dat;
        % Update regularization parameter
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
        % Define current u-vector
        u = flip(signal(n-delta-M+1:n-delta));
        xhat(n) = dot(w(:, n), u);
        error(n) = signal(n)  - xhat(n);
        % Update weights
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
        % Update weights
        w(:, n+1) = w(:, n) + lr*xhat(n)*u(:, n);
    end
end
