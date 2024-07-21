clear; close all; clc;

addpath(genpath('utils'));

% carregar weights
data = load('weights\weights.mat');

% medias de erros
errTsiSize = size(data.impErrTsi_array{1});
[impErrTsiAvg, shoErrTsiAvg] = calculateErrorAverage(data.impErrTsi_array, data.shoErrTsi_array, errTsiSize);

% tempo
predTime = data.predTime;
Ts_ref = data.Ts_ref;

t_ref = Ts_ref * (1:ceil((predTime(2) - predTime(1))/Ts_ref)) + 1;
tTsi_cell = cell(errTsiSize(2), 1);
for i = 1:errTsiSize(2)
    tTsi_cell{i} = t_ref;
end

% labels
labelsTsi = data.labelsTs;
 

% comparando Ts
for i = 1:errTsiSize(1)
    plotTitleImp(i) = "Erro de Impacto Interpolado 5ms - " + labelsTsi.method(i);
    plotTitleSho(i) = "Erro de Disparo Interpolado 5ms - " + labelsTsi.method(i);
end
plotLegend(1) = labelsTsi.Ts(1) + " Não Interpolado";
for i = 2:errTsiSize(2)
    plotLegend(i) = labelsTsi.Ts(i) + " Interpolado para " + labelsTsi.Ts(1);
end
errorPlot(tTsi_cell, impErrTsiAvg, plotTitleImp, plotLegend, predTime, 120);
errorPlot(tTsi_cell, shoErrTsiAvg, plotTitleSho, plotLegend, predTime, 80);