%% 1_6_a
close all; clear all; clc;
load 'PCR\PCAPCR.mat'

svdX = svd(X);
svdXN = svd(Xnoise);
figure
stem(svdX, 'x', 'linewidth', 2, 'markersize', 10);hold on
stem(svdXN, 'x', 'linewidth', 2, 'markersize', 10);
xlabel('Subspace dimension')
ylabel('Singular Value Magnitude')
set(gca, 'Fontsize', 12)
title('Singular Values for signal X and noisy signal Xnoise', 'Fontsize', 14)
legend('X','Xnoise')

figure
stem((svdX-svdXN).^2, 'x', 'linewidth', 2, 'markersize', 10);
xlabel('Subspace dimension')
ylabel('Squared error')
set(gca, 'Fontsize', 12)
title('Squared Error between SVs of X and Xnoise', 'Fontsize', 14)


%% 1_6_b

close all; clear all; clc;
load 'PCR\PCAPCR.mat'


[uClean, sClean, vClean] = svd(X);

rankClean = rank(X);

[uTrain, sTrain, vTrain] = svd(Xnoise);

rankTrain = rank(Xnoise);

xTrainDenoised = uTrain(:, 1: rankClean) * sTrain(1: rankClean, 1: rankClean) * vTrain(:, 1: rankClean)';

errorNoise = abs(vecnorm(X - Xnoise)) .^ 2;
errorDenoise = abs(vecnorm(X - xTrainDenoised)) .^ 2;

figure;
stem(errorNoise, 'x','linewidth', 2, 'markersize', 10);
hold on;
stem(errorDenoise, 'x','linewidth', 2, 'markersize', 10);
legend('Xnoise', 'X');
title('Squared error between noiseless input and reconstructed filtered noisy input', 'Fontsize', 14);
set(gca, 'Fontsize', 12)
xlabel('Subspace dimension');
ylabel('Squared error');

%% 1_6_c
close all; clear all; clc;
load 'PCR\PCAPCR.mat'

% Ordinary least squares (OLS)
% regression matrix by training data
coefOls = (Xnoise' * Xnoise) \ Xnoise' * Y;
% training and testing output by OLS regression model
yTrainOls = Xnoise * coefOls;
yTestOls = Xtest * coefOls;
% errors
errorTrainOls = abs(vecnorm(Y - yTrainOls)) .^ 2;
errorTestOls = abs(vecnorm(Ytest - yTestOls)) .^ 2;

% Principle component regression (PCR)
% rank of clean signal
rankClean = rank(X);
% singular value of noisy signals
[uTrain, sTrain, vTrain] = svd(Xnoise);
[uTest, sTest, vTest] = svd(Xtest);
% reconstruct denoised signals
xTrainDenoised = uTrain(:, 1: rankClean) * sTrain(1: rankClean, 1: rankClean) * vTrain(:, 1: rankClean)';
xTestDenoised = uTest(:, 1: rankClean) * sTest(1: rankClean, 1: rankClean) * vTest(:, 1: rankClean)';
% regression matrix by training data
coefPcr = vTrain(:, 1: rankClean) / sTrain(1: rankClean, 1: rankClean) * uTrain(:, 1: rankClean)' * Y;
% training and testing output by PCR regression model
yTrainPcrDenoised = xTrainDenoised * coefPcr;
yTestPcrDenoised = xTestDenoised * coefPcr;
% errors
errorTrainPcr = abs(vecnorm(Y - yTrainPcrDenoised)) .^ 2;
errorTestPcr = abs(vecnorm(Ytest - yTestPcrDenoised)) .^ 2;


figure;
subplot(2, 1, 1);
stem(errorTrainOls, 'o','linewidth', 2, 'markersize', 10);
hold on;
stem(errorTrainPcr, 'x','linewidth', 2, 'markersize', 7);
legend('OLS', 'PCR');
title('MSE between training data and estimations');
xlabel('Subspace dimension');
ylabel('MSE');
% testing set
subplot(2, 1, 2);
stem(errorTestOls, 'o','linewidth', 2, 'markersize', 10);
hold on;
stem(errorTestPcr, 'x','linewidth', 2, 'markersize', 7);
legend('OLS', 'PCR');
title('MSE between testing data and estimations');
xlabel('Subspace dimension');
ylabel('MSE');

%% 1_6_d


close all; clear all; clc;
load 'PCR\PCAPCR.mat'

nRps = 1e2;

errorOls = cell(nRps, 1);
errorPcr = cell(nRps, 1);

% regression matrix by OLS
coefOls = (Xnoise' * Xnoise) \ Xnoise' * Y;
% rank of clean signal
rankClean = rank(X);
% singular value of noisy signals
[uTrain, sTrain, vTrain] = svd(Xnoise);
% regression matrix by PCR
coefPcr = vTrain(:, 1: rankClean) / sTrain(1: rankClean, 1: rankClean) * uTrain(:, 1: rankClean)' * Y;

% generate test data and estimate based on regression matrix
for iRp = 1: nRps
    % OLS
    [yTest, yTestOls] = regval(coefOls);
    errorOls{iRp} = abs(vecnorm(yTest - yTestOls)) .^ 2;
    % PCR
    [yTest, yTestPcr] = regval(coefPcr);
    errorPcr{iRp} = abs(vecnorm(yTest - yTestPcr)) .^ 2;
end
errorOlsAvg = mean(cell2mat(errorOls));
errorPcrAvg = mean(cell2mat(errorPcr));

figure;
stem(errorOlsAvg, 'o','linewidth', 2, 'markersize', 10);
hold on;
stem(errorPcrAvg, 'x','linewidth', 2, 'markersize', 7);
legend('OLS', 'PCR');
title('MSE between data and estimations');
set(gca, 'Fontsize', 12)
xlabel('Subspace dimension');
ylabel('MSE');
