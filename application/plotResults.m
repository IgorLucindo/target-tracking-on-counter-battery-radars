clear; close all; clc;

addpath(genpath('utils'));

% metodo 1 - roda kalman filter para cada iteracao
% metodo 2 - rodar kalman filter, ida e volta, para cada iteracao
% metodo 3 - rodar kalman filter, ida e volta, multiplas vezes para cada iteracao

% plotar resultados
plotModel();         % plotar graficos para comparacao de filtragens kf e ekf
plotPKf();           % plotar graficos de kf para comparacao de diferentes valores de P
plotPEkf();          % plotar graficos de ekf para comparacao de diferentes valores de P
plotTKf();           % plotar graficos de kf para comparacao de diferentes valores de T
plotTEkf();          % plotar graficos de ekf para comparacao de diferentes valores de T
plotInterpKf();      % plotar graficos de kf para comparacao de interpolacao
plotInterpEkf();     % plotar graficos de ekf para comparacao de interpolacao