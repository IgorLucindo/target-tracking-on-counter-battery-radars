clear; close all; clc;

addpath(genpath('..\utils'));

% carregar resultados
data = load('..\weights\weights.mat');

% medias de erros
errRadSize = size(data.impErrRad_array{1});
% [impErrTAvg, shoErrTAvg] = calculateErrorAverage(data.impErrT_array, data.shoErrT_array, errRadSize);
[impErrRadAvg, shoErrRadAvg] = calculateErrorAverage(data.impErrRad_array, data.shoErrRad_array, errRadSize);
[impErrRadEKFAvg, shoErrRadEKFAvg] = calculateErrorAverage(data.impErrRadEKF_array, data.shoErrRadEKF_array, errRadSize);
[impErrRadTiEKFAvg, shoErrRadTiEKFAvg] = calculateErrorAverage(data.impErrRadTiEKF_array, data.shoErrRadTiEKF_array, errRadSize);
ImpErrDesired = {[198 170 130]; [198 170 130]; [198 170 130]};
shoErrDesired = {[2.1 1.82 1.73]; [2.1 1.82 1.73]; [2.1 1.82 1.73]};
impErrRadAvg = [impErrRadAvg, impErrRadEKFAvg, impErrRadTiEKFAvg, ImpErrDesired];
shoErrRadAvg = [shoErrRadAvg, shoErrRadEKFAvg, shoErrRadTiEKFAvg, shoErrDesired];
% impErrRadAvg = [impErrRadAvg, impErrTAvg, ImpErrDesired];
% shoErrRadAvg = [shoErrRadAvg, shoErrTAvg shoErrDesired];
errRadSize = size(impErrRadAvg);

% tempo
predTime = data.predTime;
T_rad = 0.05;
t_rad = T_rad * (1:ceil((predTime(2) - predTime(1))/T_rad)) + 1;
t_Tirad = 0.005 * (1:ceil((predTime(2) - predTime(1))/0.005)) + 1;

tT_cell = cell(errRadSize)';
for i = 1:errRadSize(2)-2
    tT_cell{i, 1} = t_rad;
end
tT_cell{end - 1, 1} = t_Tirad;
tT_cell{end, 1} = [2 4 6];
for i = 2:errRadSize(1)
    tT_cell(:, i) = tT_cell(:, 1);
end

% labels
labelsT = data.labelsT;
 

% comparando dados de radares
for i = 1:errRadSize(1)
    plotTitleImp(i) = "Dados do Radar - Erro de Impacto - " + (T_rad*1000) + "ms - " + labelsT.method(i);
    plotTitleSho(i) = "Dados do Radar - Erro de Disparo - " + (T_rad*1000) + "ms - " + labelsT.method(i);
end
plotLegend(1) = "Kalman Filter";
plotLegend(2) = "Kalman Filter Extendido";
plotLegend(3) = "Kalman Filter Extendido Interpolado para 5ms";
plotLegend(4) = "Resultado Desejado (EKF)";
lines = ["" "" "" "--"];
errorPlot(tT_cell, impErrRadAvg, plotTitleImp, plotLegend, predTime, 3000, lines);
plotLegend(4) = "Resultado Desejado (Máxima Verossimilhança)";
errorPlot(tT_cell, shoErrRadAvg, plotTitleSho, plotLegend, predTime, 1100, lines);