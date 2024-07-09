clear; close all; clc;

addpath(genpath('utils'));

% carregar weights
data = load('weights\weights.mat');

% medias de erros
errPSize = size(data.impErrP_array{1});
errTsSize = size(data.impErrTs_array{1});
[impErrPAvg, shoErrPAvg] = errAverage(data.impErrP_array, data.shoErrP_array, errPSize);
[impErrTsAvg, shoErrTsAvg] = errAverage(data.impErrTs_array, data.shoErrTs_array, errTsSize);

% [impErrAvgArrayP, shoErrAvgArrayP, ~] = loadWeights('weights\weights_comparingP.mat');
% [impErrAvgArrayTs, shoErrAvgArrayTs, ~] = loadWeights('weights\weights_comparingTs.mat');
% impErrAvgArrayMethod = permute(impErrAvgArrayTs, [2, 1, 3]);
% shoErrAvgArrayMethod = permute(shoErrAvgArrayTs, [2, 1, 3]);

% tempo
predTime = data.predTime;
Ts_ref = data.Ts_ref;

t_10ms = 0.01 * (1:ceil((predTime(2) - predTime(1))/0.01)) + 1;
tP_cell = cell(errPSize(2), 1);
for i = 1:errPSize(2)
    tP_cell{i} = t_10ms;
end

tTs_cell = cell(errTsSize(2), 1);
for i = 1:errTsSize(2)
    tTs_cell{i} = Ts_ref * i * (1:ceil((predTime(2) - predTime(1))/(Ts_ref*i))) + 1;
end

% labels
labelsP = data.labelsP;
labelsTs = data.labelsTs;


% fullPlot(y, y_true, y_est);
% dynamicPlotRecord(y_true, y_est, impactPointPredArray, shootingPointPredArray);


% comparando P
for i = 1:errPSize(1)
    plotTitleImp(i) = "Erro de Impacto - " + labelsP.method(i) + " e " + labelsP.Ts;
    plotTitleSho(i) = "Erro de Disparo - " + labelsP.method(i) + " e " + labelsP.Ts;
end
for i = 1:errPSize(2)
    plotLegend(i) = labelsP.P(i);
end
errorPlot(tP_cell, impErrPAvg, plotTitleImp, plotLegend, predTime, 120);
errorPlot(tP_cell, shoErrPAvg, plotTitleSho, plotLegend, predTime, 80);
 

% comparando Ts
for i = 1:errTsSize(1)
    plotTitleImp(i) = "Erro de Impacto - " + labelsTs.method(i);
    plotTitleSho(i) = "Erro de Disparo - " + labelsTs.method(i);
end
for i = 1:errTsSize(2)
    plotLegend(i) = labelsTs.Ts(i);
end
errorPlot(tTs_cell, impErrTsAvg, plotTitleImp, plotLegend, predTime, 120);
errorPlot(tTs_cell, shoErrTsAvg, plotTitleSho, plotLegend, predTime, 80);



% comparando cada metodo
% plotTitle = ["Erro de Impacto - Ts = 1ms";
%              "Erro de Impacto - Ts = 5ms";
%              "Erro de Impacto - Ts = 10ms"];
% plotLegend = ["Metodo 1" "Metodo 2" "Metodo 3"];
% errorPlot(t_1ms, impErrAvgArrayMethod(1:3, :, :), plotTitle, plotLegend, 120);
% plotTitle = ["Erro de Impacto - Ts = 15ms";
%              "Erro de Impacto - Ts = 20ms";
%              "Erro de Impacto - Ts = 25ms";
%              "Erro de Impacto - Ts = 30ms"];
% plotLegend = ["Metodo 1" "Metodo 2" "Metodo 3"];
% errorPlot(t_1ms, impErrAvgArrayMethod(4:7, :, :), plotTitle, plotLegend, 120);
% plotTitle = ["Erro de Disparo - Ts = 1ms";
%              "Erro de Disparo - Ts = 5ms";
%              "Erro de Disparo - Ts = 10ms"];
% plotLegend = ["Metodo 1" "Metodo 2" "Metodo 3"];
% errorPlot(t_1ms, shoErrAvgArrayMethod(1:3, :, :), plotTitle, plotLegend, 80);
% plotTitle = ["Erro de Disparo - Ts = 15ms";
%              "Erro de Disparo - Ts = 20ms";
%              "Erro de Disparo - Ts = 25ms";
%              "Erro de Disparo - Ts = 30ms"];
% plotLegend = ["Metodo 1" "Metodo 2" "Metodo 3"];
% errorPlot(t_1ms, shoErrAvgArrayMethod(4:7, :, :), plotTitle, plotLegend, 80);