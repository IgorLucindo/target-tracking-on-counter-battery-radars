clear; close all; clc;

addpath('utils');
addpath('utils\plot');

% intervalo de tempo de previsao
predTime = [1 5];

Ts_1ms = 0.001;
Ts_10ms = 0.01;
% carregar erros anteriores em weights
arrayLength = ceil(predTime(2)/Ts_10ms);
errArraySizeP = [3 5 arrayLength];
weights_path = 'weights\weights_comparingP.mat';
[impErrAvgArrayP, shoErrAvgArrayP, ~] = loadWeights(weights_path,errArraySizeP);

arrayLength = ceil(predTime(2)/Ts_1ms);
errArraySizeTs = [3 7 arrayLength];
weights_path = 'weights\weights_comparingTs.mat';
[impErrAvgArrayTs, shoErrAvgArrayTs, ~] = loadWeights(weights_path, errArraySizeTs);

impErrAvgArrayMethod = permute(impErrAvgArrayTs, [2, 1, 3]);
shoErrAvgArrayMethod = permute(shoErrAvgArrayTs, [2, 1, 3]);


% tempo
t = (1:arrayLength)*Ts_1ms;
t_1ms = t;
t_10ms = t(10:10:end);


% fullPlot(y, y_true, y_est);
% dynamicPlotRecord(y_true, y_est, impactPointPredArray, shootingPointPredArray);


% comparando P
plotTitle = ["Erro de Impacto - Metodo 1 e Ts = 10ms";
             "Erro de Impacto - Metodo 2 e Ts = 10ms";
             "Erro de Impacto - Metodo 3 e Ts = 10ms"];
plotLegend = ["P = 1e1" "P = 1e2" "P = 1e4" "P = 1e6" "P = 1e8"];
errorPlot(t_10ms, impErrAvgArrayP, plotTitle, plotLegend, 120);
plotTitle = ["Erro de Disparo - Metodo 1 e Ts = 10ms";
             "Erro de Disparo - Metodo 2 e Ts = 10ms";
             "Erro de Disparo - Metodo 3 e Ts = 10ms"];
plotLegend = ["P = 1e1" "P = 1e2" "P = 1e4" "P = 1e6" "P = 1e8"];
errorPlot(t_10ms, shoErrAvgArrayP, plotTitle, plotLegend, 80);
 

% comparando Ts
plotTitle = ["Erro de Impacto - Metodo 1";
             "Erro de Impacto - Metodo 2";
             "Erro de Impacto - Metodo 3"];
plotLegend = ["Ts = 1ms" "Ts = 5ms" "Ts = 10ms" "Ts = 15ms" "Ts = 20ms" "Ts = 25ms" "Ts = 30ms"];
errorPlot(t_1ms, impErrAvgArrayTs, plotTitle, plotLegend, 120);
plotTitle = ["Erro de Disparo - Metodo 1";
             "Erro de Disparo - Metodo 2";
             "Erro de Disparo - Metodo 3"];
plotLegend = ["Ts = 1ms" "Ts = 5ms" "Ts = 10ms" "Ts = 15ms" "Ts = 20ms" "Ts = 25ms" "Ts = 30ms"];
errorPlot(t_1ms, shoErrAvgArrayTs, plotTitle, plotLegend, 80);



% comparando cada metodo
plotTitle = ["Erro de Impacto - Ts = 1ms";
             "Erro de Impacto - Ts = 5ms";
             "Erro de Impacto - Ts = 10ms"];
plotLegend = ["Metodo 1" "Metodo 2" "Metodo 3"];
errorPlot(t_1ms, impErrAvgArrayMethod(1:3, :, :), plotTitle, plotLegend, 120);
plotTitle = ["Erro de Impacto - Ts = 15ms";
             "Erro de Impacto - Ts = 20ms";
             "Erro de Impacto - Ts = 25ms";
             "Erro de Impacto - Ts = 30ms"];
plotLegend = ["Metodo 1" "Metodo 2" "Metodo 3"];
errorPlot(t_1ms, impErrAvgArrayMethod(4:7, :, :), plotTitle, plotLegend, 120);
plotTitle = ["Erro de Disparo - Ts = 1ms";
             "Erro de Disparo - Ts = 5ms";
             "Erro de Disparo - Ts = 10ms"];
plotLegend = ["Metodo 1" "Metodo 2" "Metodo 3"];
errorPlot(t_1ms, shoErrAvgArrayMethod(1:3, :, :), plotTitle, plotLegend, 80);
plotTitle = ["Erro de Disparo - Ts = 15ms";
             "Erro de Disparo - Ts = 20ms";
             "Erro de Disparo - Ts = 25ms";
             "Erro de Disparo - Ts = 30ms"];
plotLegend = ["Metodo 1" "Metodo 2" "Metodo 3"];
errorPlot(t_1ms, shoErrAvgArrayMethod(4:7, :, :), plotTitle, plotLegend, 80);