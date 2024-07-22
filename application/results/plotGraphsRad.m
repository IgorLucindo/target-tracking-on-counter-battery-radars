clear; close all; clc;

addpath(genpath('..\utils'));

% carregar weights
data = load('..\weights\weightsRad.mat');

% medias de erros
errRadSize = size(data.impErrRad_array{1});
[impErrTsAvg, shoErrTsAvg] = calculateErrorAverage(data.impErrTs_array, data.shoErrTs_array, errRadSize);
[impErrRadAvg, shoErrRadAvg] = calculateErrorAverage(data.impErrRad_array, data.shoErrRad_array, errRadSize);
impErrRadAvg = [impErrRadAvg, impErrTsAvg];
shoErrRadAvg = [shoErrRadAvg, shoErrTsAvg];
errRadSize = size(impErrRadAvg);

% tempo
predTime = data.predTime;
Ts_rad = 0.05;
t_rad = Ts_rad * (1:ceil((predTime(2) - predTime(1))/Ts_rad)) + 1;

tTs_cell = cell(errRadSize)';
for i = 1:errRadSize(2)
    tTs_cell{i, 1} = t_rad;
end
for i = 2:errRadSize(1)
    tTs_cell(:, i) = tTs_cell(:, 1);
end

% labels
labelsTs = data.labelsTs;
 

% comparando dados de radares
for i = 1:errRadSize(1)
    plotTitleImp(i) = "Erro de Impacto - " + (Ts_rad*1000) + "ms - " + labelsTs.method(i);
    plotTitleSho(i) = "Erro de Disparo - " + (Ts_rad*1000) + "ms - " + labelsTs.method(i);
end
plotLegend(1) = "Dados do Radar";
plotLegend(2) = "Dados Gerados (100 simulações)";
errorPlot(tTs_cell, impErrRadAvg, plotTitleImp, plotLegend, predTime, 1500);
errorPlot(tTs_cell, shoErrRadAvg, plotTitleSho, plotLegend, predTime, 120);