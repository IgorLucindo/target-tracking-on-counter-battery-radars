clear; close all; clc;

addpath('classes');
addpath('utils');
addpath('utils\prediction');
addpath('utils\simulation');

% variaveis para simulacao
numOfSim = 50;
% fixar P ou Ts
fixedVariable = 'P';


% variaveis de mruv
g = 9.81;
vel0 = [180 180 250]; a = [0 0 -g];
% chao
pos_floor = 0;
% intervalo de tempo de previsao
predTime = [1 5];

% entrada do sistema
u = a';
% matriz covariancia do ruido de processo
Q = cov(u');
% matriz covariancia do ruido de medicao
sigma2_v = [1e2; 1e2; 1e2];
R = [sigma2_v(1) 0 0;     % roubo
     0 sigma2_v(2) 0;
     0 0 sigma2_v(3)];
R = 0.1*R;                % R chutado para 10 vezes menor
% quero testar R obtido por autocovariance least-square

Ts_1ms = 0.001;
Ts_10ms = 0.01;
% carregar erros anteriores salvos em weights
switch fixedVariable
    case 'Ts'
        arrayLength = ceil(predTime(2)/Ts_10ms);
        errArraySize = [3 5 arrayLength];
        weights_path = 'weights\weights_comparingP.mat';
    case 'P'
        arrayLength = ceil(predTime(2)/Ts_1ms);
        errArraySize = [3 7 arrayLength];
        weights_path = 'weights\weights_comparingTs.mat';
end
[impErrAvgArray, shoErrAvgArray, numOfSimTotal] = loadWeights(weights_path, errArraySize);



% definir saida real
[y_true, impPt, shoPt] = setYTrue(Ts_1ms, vel0, g, u, pos_floor);
y_true_10ms = y_true(:, 10:10:end);



% simulacao
for i = 1:numOfSim
    tic;
    fprintf("\nrodando simulacao: " + i + "/" + numOfSim + "\n");
    
    % ruido de medicao
    v = generateNoise(sigma2_v, predTime(2), Ts_1ms);
    v_10ms = v(:, 10:10:end);

    switch fixedVariable
        case 'Ts'
            [impErrArray, shoErrArray] = runSimFixedTs(Ts_10ms, v_10ms, y_true_10ms, g, u, impPt, shoPt, pos_floor, predTime, Q, R);
        case 'P'
            [impErrArray, shoErrArray] = runSimFixedP(v, y_true, g, u, impPt, shoPt, pos_floor, predTime, Q, R);
    end

    % atualiza media dos erros e salva em weights
    impErrAvgArray = (impErrAvgArray*(numOfSimTotal-1) + impErrArray)/numOfSimTotal;
    shoErrAvgArray = (shoErrAvgArray*(numOfSimTotal-1) + shoErrArray)/numOfSimTotal;
    numOfSimTotal = numOfSimTotal + 1;
    save(weights_path, 'impErrAvgArray', 'shoErrAvgArray', 'numOfSimTotal');
    toc
end