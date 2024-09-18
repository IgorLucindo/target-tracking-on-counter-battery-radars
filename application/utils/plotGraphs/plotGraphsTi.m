clear; close all; clc;

addpath(genpath('..\utils'));

% carregar resultados
data = load('..\weights\weightsTiKF.mat');

% medias de erros
errTiSize = size(data.impErrTi_array{1});
[impErrTiAvg, shoErrTiAvg] = calculateErrorAverage(data.impErrTi_array, data.shoErrTi_array, errTiSize);
[impErrTAvg, shoErrTAvg] = calculateErrorAverage(data.impErrT_array, data.shoErrT_array, errTiSize);
impErrTAvg = [impErrTAvg(1, 2:end); impErrTiAvg(1, 2:end)]';
shoErrTAvg = [shoErrTAvg(1, 2:end); shoErrTiAvg(1, 2:end)]';
errSize = size(impErrTAvg);

% tempo
predTime = data.predTime;
T_ref = data.T_ref;

t_cell = cell(errSize)';
for i = 1:errSize(1)
    t_cell{1, i} = T_ref * (i+1) * (1:ceil((predTime(2) - predTime(1))/T_ref/(i+1))) + 1;
end
for i = 1:errSize(1)
    t_cell{2, i} = T_ref * (1:ceil((predTime(2) - predTime(1))/T_ref)) + 1;
end


% labels
labels = data.labelsT;


% comparando mesmo T
for i = 1:errSize(1)
    plotTitleImp(i) = "Erro de Impacto - " + labels.T(i+1);
    plotTitleSho(i) = "Erro de Disparo - " + labels.T(i+1);
end
plotLegend(1) = "Não Interpolado";
plotLegend(2) = "Interpolado para " + labels.T(1);
lines = ["b" "r--"];
errorPlot(t_cell, impErrTAvg, plotTitleImp, plotLegend, predTime, 120, lines);
errorPlot(t_cell, shoErrTAvg, plotTitleSho, plotLegend, predTime, 120, lines);