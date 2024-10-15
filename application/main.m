clear; close all; clc;

addpath('classes');
addpath(genpath('utils'));


% comparacao
comparisons = ["models", "P_kf", "P_ekf", "T_kf", "T_ekf", "interp_kf", "interp_ekf"];
% 1 - models      - kf e ekf em trajetoria simulada
% 2 - P_kf        - matriz inicial P diferentes para kf em trajetoria simulada
% 3 - P_ekf       - matriz inicial P diferentes para ekf em trajetoria simulada
% 4 - T_kf        - periodos de amostragem T diferentes para kf em trajetoria simulada
% 5 - T_ekf       - periodos de amostragem T diferentes para ekf em trajetoria simulada
% 6 - interp_kf   - kf com e sem interpolacao em trajetoria simulada
% 7 - interp_ekf  - ekf com e sem interpolacao em trajetoria simulada
% 8 - rad_kf      - kf em trajetoria real a partir de dados de radares (nao implementado)
% 9 - rad_ekf     - ekf em trajetoria real a partir de dados de radares (nao implementado)


% escolha a comparacao
compare = comparisons(3);

% variaveis de simulacao
% numero de simulacoes
numOfSim = 100;
% intervalo de tempo de previsao
predTime = [1 5];


% rodar simulacoes baseado na comparacao
runSimWithComparison(compare, numOfSim, predTime);