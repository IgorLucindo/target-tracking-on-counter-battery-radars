% plotar graficos para comparacao de modelos
function plotModel()
    % retornar se arquivo de resultados nao existir
    resultsFilePath = "results\weights_model.mat";
    if exist(resultsFilePath, 'file') == 0
        return
    end

    % carregar resultados
    data = load(resultsFilePath);

    % medias de erros
    errTSize = size(data.impErrT_array{1});
    [impErrTAvg, shoErrTAvg] = calculateErrorAverage(data.impErrT_array, data.shoErrT_array, errTSize);
    [impErrTEKFAvg, shoErrTEKFAvg] = calculateErrorAverage(data.impErrTEKF_array, data.shoErrTEKF_array, errTSize);
    impErrAvg = [impErrTAvg(:, 1) impErrTEKFAvg(:, 1)];
    shoErrAvg = [shoErrTAvg(:, 1) shoErrTEKFAvg(:, 1)];
    errSize = size(impErrAvg);


    % tempo
    predTime = data.predTime;
    T_ref = data.T_ref;

    t_cell = cell(errSize)';
    for i = 1:errSize(2)
        t_cell{i, 1} = T_ref * (1:ceil((predTime(2) - predTime(1))/(T_ref))) + 1;
    end
    for i = 2:errSize(1)
        t_cell(:, i) = t_cell(:, 1);
    end

    % labels
    labelsT = data.labelsT;

    % comparando modelos
    for i = 1:errSize(1)
        plotTitleImp(i) = "Erro de Impacto - " + labelsT.T(1) + " - " + labelsT.method(i);
        plotTitleSho(i) = "Erro de Disparo - " + labelsT.T(1) + " - " + labelsT.method(i);
    end
    plotLegend(1) = "Kalman Filter";
    plotLegend(2) = "Kalman Filter Extendido";
    lines = strings(1, errSize(2));
    errorPlot(t_cell, impErrAvg, plotTitleImp, plotLegend, predTime, 500, lines);
    errorPlot(t_cell, shoErrAvg, plotTitleSho, plotLegend, predTime, 500, lines);
end