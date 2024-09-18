clear; close all; clc;

addpath('classes');
addpath(genpath('utils'));


% comparacao
% comparar kf e ekf em trajetoria simulada
compare = "models";


% variaveis de simulacao
% numero de simulacoes
numOfSim = 10;
% intervalo de tempo de previsao
predTime = [1 5];


% rodar simulacoes baseado na comparacao
runSimWithComparison(compare, numOfSim, predTime);