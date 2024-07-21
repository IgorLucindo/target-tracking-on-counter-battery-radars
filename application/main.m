clear; close all; clc;

addpath('classes');
addpath(genpath('utils'));

% comparar P ou Ts
% COMPARE = 'P';
% COMPARE = 'Ts';
% COMPARE = 'Tsi';
COMPARE = 'Rad';
% numero de simulacoes
numOfSim = 13;
% ao mudar esses valores, delete o arquivo weights.mat
numOfP = 4;
numOfTs = 4;
Ts_ref = 0.005;


% variaveis de mruv
g = 9.81;
v0 = [180 180 250]; a = [0 0 -g];
% chao
p_floor = 0;
% intervalo de tempo de previsao
predTime = [1 5];

% entrada do sistema
u = a';
% matriz covariancia do ruido de processo
Q = cov(u');
% matriz covariancia do ruido de medicao
sigma2_n = [1e2; 1e2; 1e2];
R = [sigma2_n(1) 0 0;     % roubo
     0 sigma2_n(2) 0;
     0 0 sigma2_n(3)];
R = 0.1*R;                % R chutado para 10 vezes menor
% quero testar R obtido por autocovariance least-square

% periodo de amostragem
Ts_array = Ts_ref * (1:numOfTs);
Ts_10ms = 0.01;

% definir saida real
[y_true, impPt, shoPt] = setYTrue(Ts_ref, v0, g, u, p_floor);

% array de diferentes P
P_array = zeros(6, 6, numOfP);
for i = 1:numOfP
    P_array(:, :, i) = 10^(2*i-1) * eye(6);
end

% criar weights caso nao exista e definir labels
if exist('weights\weights.mat', 'file') == 0
    createWeights(Ts_array, Ts_10ms, P_array, predTime);
end

% carregar weights
data = load('weights\weights.mat');

% numero de simulacoes anteriores
numOfPrevSim = data.(['numOfSim' COMPARE]);

% simulacoes
for i = 1:numOfSim
    tic;
    fprintf("\nrodando simulacao: " + (numOfPrevSim+i) + "/" + (numOfPrevSim+numOfSim) + "\n");
    
    % radar
    switch COMPARE
        case 'Rad'
            [y, impPt, shoPt, Ts_rad] = getRadarTrajectoryDynamically(i);
        otherwise
            y = createTrajectory(y_true, sigma2_n, predTime, Ts_ref);

            rescale_10ms = Ts_10ms/Ts_ref;
            y_10ms = y(:, rescale_10ms:rescale_10ms:end);
    end

    % rodar simulacao
    switch COMPARE
        case 'P'
            [data.impErrP_array{end + 1}, data.shoErrP_array{end + 1}] = runSimP(Ts_10ms, P_array, y_10ms, g, u, impPt, shoPt, p_floor, predTime, Q, R);
        case 'Ts'
            [data.impErrTs_array{end + 1}, data.shoErrTs_array{end + 1}] = runSimTs(Ts_array, y, g, u, impPt, shoPt, p_floor, predTime, Q, R);
        case 'Tsi'
            [data.impErrTsi_array{end + 1}, data.shoErrTsi_array{end + 1}] = runSimTsi(Ts_array, y, g, u, impPt, shoPt, p_floor, predTime, Q, R);
        case 'Rad'
            [data.impErrRad_array{end + 1}, data.shoErrRad_array{end + 1}] = runSimTs([Ts_rad], y, g, u, impPt, shoPt, p_floor, predTime, Q, R);
    end
    data.(['numOfSim' COMPARE]) = numOfPrevSim + i;
    
    % salva erros em weights
    save('weights\weights.mat', '-struct', 'data');
    toc;
end