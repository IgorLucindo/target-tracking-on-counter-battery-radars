clear; close all; clc;

addpath(genpath('..\utils'));

% carregar weights
data = load('..\weights\weights10s.mat');

% medias de erros
errTsiSize = size(data.impErrTsi_array{1});
[impErrTsiAvg, shoErrTsiAvg] = calculateErrorAverage(data.impErrTsi_array, data.shoErrTsi_array, errTsiSize);
[impErrTsAvg, shoErrTsAvg] = calculateErrorAverage(data.impErrTs_array, data.shoErrTs_array, errTsiSize);
impErrTsAvg = [impErrTsAvg(1, :); impErrTsiAvg(1, :)]';
shoErrTsAvg = [shoErrTsAvg(1, :); shoErrTsiAvg(1, :)]';
errTsSize = size(impErrTsAvg);

% tempo
predTime = data.predTime;
Ts_ref = data.Ts_ref;

t_ref = Ts_ref * (1:ceil((predTime(2) - predTime(1))/Ts_ref)) + 1;
tTsi_cell = cell(errTsiSize)';
for i = 1:errTsiSize(2)
    tTsi_cell{i, 1} = t_ref;
end
for i = 2:errTsiSize(1)
    tTsi_cell(:, i) = tTsi_cell(:, 1);
end

tTs_cell = cell(errTsSize)';
for i = 1:errTsSize(1)
    tTs_cell{1, i} = Ts_ref * i * (1:ceil((predTime(2) - predTime(1))/Ts_ref/i)) + 1;
end
for i = 1:errTsSize(1)
    tTs_cell{2, i} = Ts_ref * (1:ceil((predTime(2) - predTime(1))/Ts_ref)) + 1;
end


% labels
labelsTsi = data.labelsTs;
 

% comparando Ts diferentes
for i = 1:errTsiSize(1)
    plotTitleImp(i) = "Erro de Impacto Interpolado para 5ms - " + labelsTsi.method(i);
    plotTitleSho(i) = "Erro de Disparo Interpolado para 5ms - " + labelsTsi.method(i);
end
plotLegend(1) = labelsTsi.Ts(1) + " Não Interpolado";
for i = 2:errTsiSize(2)
    plotLegend(i) = labelsTsi.Ts(i) + " Interpolado para " + labelsTsi.Ts(1);
end
errorPlot(tTsi_cell, impErrTsiAvg, plotTitleImp, plotLegend, predTime, 120);
errorPlot(tTsi_cell, shoErrTsiAvg, plotTitleSho, plotLegend, predTime, 80);


% comparando mesmo Ts
for i = 1:errTsSize(1)
    plotTitleImp(i) = "Erro de Impacto - " + labelsTsi.Ts(i) + " - metodo 1";
    plotTitleSho(i) = "Erro de Disparo - " + labelsTsi.Ts(i) + " - metodo 1";
end
plotLegend(1) = "Não Interpolado";
plotLegend(2) = "Interpolado para " + labelsTsi.Ts(1);
errorPlot(tTs_cell, impErrTsAvg, plotTitleImp, plotLegend, predTime, 120);
errorPlot(tTs_cell, shoErrTsAvg, plotTitleSho, plotLegend, predTime, 80);