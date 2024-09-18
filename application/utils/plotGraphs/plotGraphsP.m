clear; close all; clc;

addpath(genpath('..\utils'));

% carregar resultados
data = load('..\weights\weights.mat');

% medias de erros
errPSize = size(data.impErrP_array{1});
[impErrPAvg, shoErrPAvg] = calculateErrorAverage(data.impErrP_array, data.shoErrP_array, errPSize);

% tempo
predTime = data.predTime;
T_ref = data.T_ref;

t_10ms = 0.01 * (1:ceil((predTime(2) - predTime(1))/0.01)) + 1;
tP_cell = cell(errPSize)';
for i = 1:errPSize(2)
    tP_cell{i, 1} = t_10ms;
end
for i = 2:errPSize(1)
    tP_cell(:, i) = tP_cell(:, 1);
end

% labels
labelsP = data.labelsP;


% comparando P
for i = 1:errPSize(1)
    plotTitleImp(i) = "Erro de Impacto - " + labelsP.method(i) + " e " + labelsP.T;
    plotTitleSho(i) = "Erro de Disparo - " + labelsP.method(i) + " e " + labelsP.T;
end
for i = 1:errPSize(2)
    plotLegend(i) = labelsP.P(i);
end
errorPlot(tP_cell, impErrPAvg, plotTitleImp, plotLegend, predTime, 120);
errorPlot(tP_cell, shoErrPAvg, plotTitleSho, plotLegend, predTime, 120);