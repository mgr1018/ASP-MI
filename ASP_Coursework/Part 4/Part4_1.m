%% 4_1
clc; clear variables; close all;
load('time-series.mat');

y = detrend(y);
order = 4;
x = [zeros(order,1); y];
[~,error,y_hat] = lms_learning(y, x, order, 1e-5, 0, 'standard',0,false,0,false);

figure(1); 
set(gca,'fontsize', 14);
plot(1:length(y), y); hold on;
plot(1:length(y_hat), y_hat);
xlabel('Time (n)');
ylabel('Y'); title('signal and LMS Prediction');
legend('Signal', 'LMS'); hold off;

MSE = (error'*error)/1000;
Rp = pow2db(var(y_hat)/var(error));
fprintf(sprintf('MSE: %.4f, Rp: %.4f',MSE,Rp));

%% 2/3 
clc; clear variables; close all;
load('time-series.mat');
order = 4;
%y = detrend(y); % remove mean for part 2/3, comment out line for part 4
x = [zeros(order,1); y]; % input to AR(4) model

%part 2
%[params,error,y_hat] = lms_learning(y, x, order, 2e-7, 0, 'nonlinear',0,false,1,false);

%part 3
%[params,error,y_hat] = lms_learning(y, x, order, 2e-7, 0, 'nonlinear',0,false,50,false);

%part 4
%[params,error,y_hat] = lms_learning(y, x, order, 2e-7, 0, 'nonlinear',0,false,50,true);

%part 5
[params5,error,y_hat5] = lms_learning(y, x, order, 2e-7, 0, 'nonlinear',0,true,50,true);

figure(1); 
set(gca,'fontsize', 14);
plot(1:length(y), y);hold on;
%plot(1:length(y_hat), y_hat);hold on
plot(1:length(y_hat5), y_hat5);hold on %only for part 4-5
xlabel('Time (n)');xlim([0,150])
ylabel('Y'); title('Signal and its Dynamical-Perceptron Prediction (tanh non-linearity)');
legend('Signal', 'tanh activated LMS)','tanh activated LMS, pre-trained weights)'); hold off;


MSE = (error'*error)/1000;
Rp = pow2db(var(y_hat5)/var(error));
fprintf(sprintf('MSE: %.4f, Rp: %.4f',MSE,Rp));

%% fcts

function [params, err, y_hat] = lms_learning(y, x, order, lr, gamma, opt, lr_a, pt,a,bias)

    N = length(y);
    err = zeros(N,1);
    y_hat = zeros(N,1);

    if pt
        [params, a] = pretraining(y,x,100,20,order,lr,gamma,lr_a,opt);
    else
        if bias
         params= zeros(order+1,N);
        else
         params = zeros(order,N);
        end
    end

    for n = 1:N
        if strcmp(opt,'standard')
            if bias
                y_hat(n) = params(:,n).'*[1;x(n+order-1:-1:n)];
            else
                 y_hat(n) = params(:,n).'*x(n+order-1:-1:n);
            end
            act_factor = 1;
        elseif strcmp(opt,'nonlinear')
            if bias
                nonlinear = tanh(params(:,n).'*[1;x(n+order-1:-1:n)]);
            else
                nonlinear = tanh(params(:,n).'*x(n+order-1:-1:n));
            end
            y_hat(n) = a*nonlinear;
            act_factor = a*(1-(nonlinear.^2));
        else
            error("err");
        end
        err(n) = y(n) - y_hat(n);
        if n < N
            if bias
            params(:, n+1) = (1-lr*gamma)*params(:, n)...
                + lr*err(n)*act_factor*[1;x(n+order-1:-1:n)];
            else
                params(:, n+1) = (1-lr*gamma)*params(:, n)...
                + lr*err(n)*act_factor*x(n+order-1:-1:n);
            end
            if strcmp(opt,'nonlinear')
                a = a + lr_a*err(n)*nonlinear;
            end
        end
    end
end


function [params, a] = pretraining(y,x,n_epochs,K,order,lr,gamma,lr_a,opt)

    N = length(y);
    params = zeros(order+1,N); 
    weights = zeros(order+1,1);
    err = zeros(N,1);
    y_hat = zeros(N,1);
    a = 40; 
    for e = 1:n_epochs
        for k = 1:K
            if strcmp(opt,'standard')
                y_hat(k) = weights.'*[1;x(k+order-1:-1:k)];
                act_factor = 1;
            elseif strcmp(opt,'nonlinear')
                nonlinear = tanh(weights.'*[1;x(k+order-1:-1:k)]);
                y_hat(k) = a*nonlinear;
                act_factor = a*(1-(nonlinear.^2));
            else
                error("err");
            end
            err(k) = y(k) - y_hat(k);
            if k < K
                weights = (1-lr*gamma)*weights...
                    + lr*err(k)*act_factor*[1;x(k+order-1:-1:k)];
                if strcmp(opt,'nonlinear')
                    a = a + lr_a*err(k)*nonlinear;
                end
            end
        end
    end
    params(:,1:K) = repmat(weights, 1, K);
end
