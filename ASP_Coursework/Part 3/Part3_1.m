%% 3_1_a 
clc; clear variables; close all;
b1 = 1.5 + 1i; b2 = 2.5 - 0.5i;
realizations = 100;
input = wgn(1000,realizations,0,'complex');
data = zeros(1000,realizations,'like',1i);
error_store = zeros(2, 1000,realizations);

for idx = 1:realizations
    for jdx = 1:1000
        if jdx == 1
            data(jdx,idx) = 0i;
        else
            data(jdx,idx) = input(jdx,idx) + b1*input(jdx-1,idx) + b2*conj(input(jdx-1,idx));
        end
    end
    [~,error_store(1, :, idx),~] = aclms(data(:,idx),input(:,idx),1, 0.1, 0);
    [~,error_store(2, :, idx),~] = clms(data(:,idx),input(:,idx),1, 0.1, 0);
    
end

mean_error = mean(abs(error_store).^2,3);

figure(1); subplot(2,1,1); hold on; set(gca,'fontsize', 16);
plot(real(data(:)),imag(data(:)),'b.'); 
plot(real(input(:)),imag(input(:)), 'r.'); 
legend('WLMA','WGN')
title('Distribution of circular WGN and WLMA')
xlabel('Re(z)'); ylabel('Im(z)'); hold off;

subplot(2,1,2); hold on; set(gca,'fontsize', 18);
plot([1:1000], pow2db(mean_error(1,:)));
plot([1:1000], pow2db(mean_error(2,:)));
ylabel('$10 log(|e(n)|^2)$','Interpreter','latex'); xlabel('Samples');
title('Learning Curve'); legend('ACMLS', 'CMLS');
hold off;



%% 3_1_b
clc; clear variables; close all;
complex_v = zeros(3,5000,'like',1i);
names = {'high-wind','medium-wind','low-wind'};

figure(1);
for i = 1:3
    data = load(names{i});
    data_array = complex([data.v_east(:),data.v_north(:)]);
    complex_v(i,:) = complex(data.v_east,data.v_north);
    [eta,rho] = circularity(complex_v(i,:))
    str = split(names{i},'-');
    tit = strcat(upper(str{1}(1)), str{1}(2:end), sprintf(' Wind, circularity = %.3f',eta));

    subplot(1,3,i); hold on; set(gca,'fontsize', 16);
    plot(complex_v(i,:), '.');    
    xlabel('Re(z)'); ylabel('Im(z)');
    title(tit );
    hold off;

end


n_orders = 20;
MSPE = zeros(2,3,n_orders); 
lrs = [0.001, 0.01,0.1]; 
count = 0;
figure(2);
for j = 1:3
    count = count + 1 ;
    subplot(1, 3,count); hold on; set(gca,'fontsize', 18);
    for i = 1:2
        input = delayseq(complex_v(j,:).',1);
        for l = 1:n_orders
            if i == 1
                [~,err,~] = ...
                    clms(complex_v(j,:).',input,l-1, lrs(j), 0);
            else
                [~,err,~] = ...
                    aclms(complex_v(j,:).',input,l-1, lrs(j), 0);
            end
            sq_err = abs(err).^2;
            MSPE(i,j,l) = mean(sq_err);
        end
        plot(1:n_orders, 10*log10(squeeze(MSPE(i,j,:))));

    end
    str = split(names{j},'-');
    tit = strcat(upper(str{1}(1)), str{1}(2:end), sprintf(' Wind lr = %.3f',lrs(j)));
    legend('CLMS','ACLMS');
    xlabel('Filter Order'); ylabel('MSPE (dB)');
    title(tit);
    hold off;
end
%% 3_1_c
clc; clear variables; close all;
f0 = 50; 
fs = 5000;
V = 20; 
phi = 3*pi/7;
n = [1:10000]; 
v = sqrt(3/2)*V*exp(i*2*pi*(f0/fs)*n + phi);
[eta1,~] = circularity(v);
Va = 10;
eta = zeros(length(Va));
eta2 = zeros(length(Va));
delta_b = pi/3;


    A = (sqrt(6)/6)*(Va + 2*V);
    B = (sqrt(6)/6)*(Va + V*exp(-i*2*pi/3) + V*exp(i*2*pi/3));
    v_unbalanced1 = A*exp(i*2*pi*(f0/fs)*n + phi) + ...
    B*exp(-i*2*pi*(f0/fs)*n + phi);
    [eta,~] = circularity(v_unbalanced1);

    A = (sqrt(6)/6)*V*(2 + exp(i*delta_b));
    B = (sqrt(6)/6)*V*(1 + exp(-i*(delta_b+ 2*pi/3)) + exp(i*2*pi/3));
    v_unbalanced2 = A*exp(i*2*pi*(f0/fs)*n + phi) + ...
    B*exp(-i*2*pi*(f0/fs)*n + phi);
    [eta2,~] = circularity(v_unbalanced2);


figure(1); 
hold on; set(gca,'fontsize', 14);
plot(v,'b.'); plot(v_unbalanced1,'r.'); plot(v_unbalanced2,'c.');
xlabel('Re(z)'); ylabel('Im(z)');ylim([-100,150])
title('Circularity diagrams of complex Voltages');

lgd = legend(sprintf('Balanced system, circularity = %.3f',eta1), sprintf('Unbalanced magnitudes, circularity = %.3f',eta), sprintf('Unbalanced phasess, circularity = %.3f',eta2));
set(lgd, 'Interpreter','latex');

%% 3_1_e
clc; clear variables; close all;
f0 = 50; 
fs = 500; 
V = 10;
phi = 4*pi/7;
n = [1:100]; 
v = sqrt(3/2)*V*exp(i*2*pi*(f0/fs)*n + phi);
Va = 20; delta_b = pi/3;

A = (sqrt(6)/6)*(Va + 2*V);
B = (sqrt(6)/6)*(Va + V*exp(-i*2*pi/3) + V*exp(i*2*pi/3));
v_unbalanced1 = A*exp(i*2*pi*(f0/fs)*n + phi) + ...
B*exp(-i*2*pi*(f0/fs)*n + phi);
[eta(1),~] = circularity(v_unbalanced1);

A = (sqrt(6)/6)*V*(2 + exp(i*delta_b));
B = (sqrt(6)/6)*V*(1 + exp(-i*(delta_b + 2*pi/3)) + exp(i*2*pi/3));
v_unbalanced2 = A*exp(i*2*pi*(f0/fs)*n + phi) + ...
B*exp(-i*2*pi*(f0/fs)*n + phi);
[eta(2),~] = circularity(v_unbalanced2);

balanced_input = delayseq(v',1); 
[h_clms_balanced, err_clms_balanced, ~] = ...
    clms(v',balanced_input,0,0.0001,0);
[q_aclms_balanced, err_aclms_balanced, ~] = ...
    aclms(v',balanced_input,0,0.0001,0);

unbalanced_input1 = delayseq(v_unbalanced1', 1); 
[h_clms_unbalanced1, err_clms_unbalanced1, ~] = ...
    clms(v_unbalanced1',unbalanced_input1,0,0.0001,0);
[q_aclms_unbalanced1, err_aclms_unbalanced1, ~] = ...
    aclms(v_unbalanced1',unbalanced_input1,0,0.00001,0);

unbalanced_input2 = delayseq(v_unbalanced2', 1);
[h_clms_unbalanced2, err_clms_unbalanced2, ~] = ...
    clms(v_unbalanced2',unbalanced_input2,0,0.0001,0);
[q_aclms_unbalanced2, err_aclms_unbalanced2, ~] = ...
    aclms(v_unbalanced2',unbalanced_input2,0,0.0001,0);

f_clms_balanced = SLAR(h_clms_balanced,fs);
f_aclms_balanced = FWL(q_aclms_balanced,fs);
f_clms_unbalanced1 = SLAR(h_clms_unbalanced1,fs);
f_aclms_unbalanced1 = FWL(q_aclms_unbalanced1,fs);
f_clms_unbalanced2 = SLAR(h_clms_unbalanced2,fs);
f_aclms_unbalanced2 = FWL(q_aclms_unbalanced2,fs);

figure(1); subplot(1,3,1); hold on; set(gca,'fontsize', 16);
plot([1:length(f_clms_balanced)-2], f_clms_balanced(3:end),'LineWidth',1);
plot([1:length(f_aclms_balanced)-2],f_aclms_balanced(3:end),'LineWidth',1);
yticks([0:25:100]);
xlabel('Time (n)'); ylabel('Frequency (Hz)');
title('Balanced System','Interpreter','Latex'); 
legend('CLMS','ACLMS'); ylim([0,100]);

subplot(1,3,2); hold on; set(gca,'fontsize', 16);
plot([1:length(f_clms_unbalanced1)-2], f_clms_unbalanced1(3:end),'LineWidth',1);
plot([1:length(f_aclms_unbalanced1)-2],f_aclms_unbalanced1(3:end),'LineWidth',1);
yticks([0,25,50,75,100]);
xlabel('Time (n)'); ylabel('Frequency (Hz)');
title('Unbalanced Magnitudes','Interpreter','Latex'); 
legend('CLMS','ACLMS'); ylim([0,100]);

subplot(1,3,3); hold on; set(gca,'fontsize', 16);
plot([1:length(f_clms_unbalanced2)-2], f_clms_unbalanced2(3:end),'LineWidth',1);
plot([1:length(f_aclms_unbalanced2)-2],f_aclms_unbalanced2(3:end),'LineWidth',1);
yticks([0,25,50,75,100]);
xlabel('Time (n)'); ylabel('Frequency (Hz)');
title('Unbalanced Phases','Interpreter','Latex'); 
legend('CLMS','ACLMS'); ylim([0,100]);

%% functs

function f = SLAR(h,fs)
    f = (fs/(2*pi))*atan2(imag(h),real(h));
end


function f = FWL(q,fs)
    h = q(1,:); 
    g = q(2,:);
    f = (fs/(2*pi))*...
        atan2(real((sqrt(imag(h).^2 - abs(g).^2))),real(h));
end

function [q, error, y_hat] = aclms(y, x, M, lr, leak)
    q = zeros(2*(M+1), length(x),'like',1i);  
    error = zeros(size(x),'like',1i);
    y_hat = zeros(size(x),'like',1i);
    x_pad = [zeros(M,1,'like',1i); x]; % zero-padding
    for n = 1:length(x)
        x_a = [x_pad(n+M:-1:n); conj(x_pad(n+M:-1:n))];
        y_hat(n) = q(:,n)'*x_a;
        error(n) = y(n) - y_hat(n);
        if n < length(x)
            q(:, n+1) = (1-lr*leak)*q(:, n) + lr*conj(error(n))*x_a;
        end
    end
end

function [h, error, y_hat] = clms(y, x, M, lr, leak)

    h = zeros(M+1, length(x),'like',1i); 
    error = zeros(size(x),'like',1i);
    y_hat = zeros(size(x),'like',1i);
    x_pad = [zeros(M,1); x]; 
    for n = 1:length(x)
        y_hat(n) = h(:,n)'*x_pad(n+M:-1:n); 
        error(n) = y(n) - y_hat(n);
        if n < length(x)
            h(:, n+1) = (1-lr*leak)*h(:, n) + lr*conj(error(n))*x_pad(n+M:-1:n);
        end
    end
end

function [eta,rho] = circularity(data)
    p = mean(data.*data); c = mean(abs(data).^2);
    rho = p/c;
    eta = abs(rho);
end
