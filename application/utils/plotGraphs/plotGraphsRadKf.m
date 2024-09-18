clear; close all; clc;

addpath(genpath('..\utils'));

% carregar resultados
data = load('..\weights\weightsRadKf.mat');

% medias de erros
errRadSize = size(data.impErrRad_array{1});
[impErrTAvg, shoErrTAvg] = calculateErrorAverage(data.impErrT_array, data.shoErrT_array, errRadSize);
[impErrRadAvg, shoErrRadAvg] = calculateErrorAverage(data.impErrRad_array, data.shoErrRad_array, errRadSize);
impErrAvg = [impErrRadAvg, impErrTAvg];
shoErrAvg = [shoErrRadAvg, shoErrTAvg];
errSize = size(impErrAvg);

% desvio padrao de erros
[impErrRadStd, shoErrRadStd] = calculateErrorStd(data.impErrRad_array, data.shoErrRad_array, impErrRadAvg, shoErrRadAvg, errRadSize);

% tempo
predTime = data.predTime;
T_rad = 0.05;
t_rad = T_rad * (1:ceil((predTime(2) - predTime(1))/T_rad)) + 1;

t_cell = cell(errSize)';
for i = 1:errSize(2)
    t_cell{i, 1} = t_rad;
end
for i = 2:errSize(1)
    t_cell(:, i) = t_cell(:, 1);
end

% labels
labelsT = data.labelsT;
 

% comparando dados de radares
for i = 1:errSize(1)
    plotTitleImp(i) = "Dados do Radar - Erro de Impacto - " + (T_rad*1000) + "ms - " + labelsT.method(i);
    plotTitleSho(i) = "Dados do Radar - Erro de Disparo - " + (T_rad*1000) + "ms - " + labelsT.method(i);
end
plotLegend(1) = "Dados Reais";
plotLegend(2) = "Dados Simulados";
lines = ["b" ""];
errorPlot(t_cell, impErrAvg, plotTitleImp, plotLegend, predTime, 1000, lines);
errorPlot(t_cell, shoErrAvg, plotTitleSho, plotLegend, predTime, 200, lines);


% desvio padrao de dados dos radares
lines = ["b"];
errorPlotStd(t_cell, impErrRadAvg, impErrRadStd, plotTitleImp, predTime, 1000, lines);
errorPlotStd(t_cell, shoErrRadAvg, shoErrRadStd, plotTitleSho, predTime, 200, lines);