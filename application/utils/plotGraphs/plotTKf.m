% plotar graficos de kf para comparacao de diferentes valores de T
function plotTKf()
    % retornar se arquivo de resultados nao existir
    resultsFilePath = "results\weights_T_kf.mat";
    if exist(resultsFilePath, 'file') == 0
        return
    end

    % carregar resultados
    data = load(resultsFilePath);

    % medias de erros
    errSize = size(data.impErrT_array{1});
    numOfSim = size(data.impErrT_array, 2);
    [impErrAvg, shoErrAvg] = calculateErrorAverage(data.impErrT_array, data.shoErrT_array, errSize);

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
        plotTitle(i) = labelsT.method(i);
    end
    for i = 1:errSize(2)
        plotLegend(i) = "T = " + labelsT.T(i);
    end
    lines = strings(1, errSize(2));
    lim = [0 200];
    errorPlot(t_cell, impErrAvg, plotTitle, plotLegend, predTime, lim, lines);
    sgtitle("Erro de Impacto - Trajetória Simulada sem Arrasto - " + numOfSim + " Simulações");
    errorPlot(t_cell, shoErrAvg, plotTitle, plotLegend, predTime, lim, lines);
    sgtitle("Erro de Disparo - Trajetória Simulada sem Arrasto - " + numOfSim + " Simulações");
end