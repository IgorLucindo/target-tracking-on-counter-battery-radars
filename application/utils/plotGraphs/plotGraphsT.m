clear; close all; clc;

addpath(genpath('..\utils'));

% carregar resultados
data = load('..\weights\weights.mat');

% medias de erros
errSize = size(data.impErrT_array{1});
[impErrTAvg, shoErrTAvg] = calculateErrorAverage(data.impErrT_array, data.shoErrT_array, errSize);

% [impErrAvgArrayP, shoErrAvgArrayP, ~] = loadWeights('weights\weights_comparingP.mat');
% [impErrAvgArrayT, shoErrAvgArrayT, ~] = loadWeights('weights\weights_comparingT.mat');
% impErrAvgArrayMethod = permute(impErrAvgArrayTs, [2, 1, 3]);
% shoErrAvgArrayMethod = permute(shoErrAvgArrayTs, [2, 1, 3]);

% tempo
predTime = data.predTime;
T_ref = data.T_ref;

t_cell = cell(errSize)';
for i = 1:errSize(2)
    t_cell{i, 1} = T_ref * i * (1:ceil((predTime(2) - predTime(1))/(T_ref*i))) + 1;
end
for i = 2:errSize(1)
    t_cell(:, i) = t_cell(:, 1);
end

% labels
labelsT = data.labelsT;

% comparando T
for i = 1:errSize(1)
    plotTitleImp(i) = "Erro de Impacto - " + labelsT.method(i);
    plotTitleSho(i) = "Erro de Disparo - " + labelsT.method(i);
end
for i = 1:errSize(2)
    plotLegend(i) = labelsT.T(i);
end
lines = strings(1, errSize(2));
errorPlot(t_cell, impErrTAvg, plotTitleImp, plotLegend, predTime, 120, lines);
errorPlot(t_cell, shoErrTAvg, plotTitleSho, plotLegend, predTime, 120, lines);



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