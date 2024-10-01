% plotar graficos de kf para comparacao de diferentes valores de P
function plotPKf()
    % retornar se arquivo de resultados nao existir
    resultsFilePath = "results\weights_P_kf.mat";
    if exist(resultsFilePath, 'file') == 0
        return
    end

    % carregar resultados
    data = load(resultsFilePath);

    % medias de erros
    errSize = size(data.impErrP_array{1});
    numOfSim = size(data.impErrP_array, 2);
    [impErrAvg, shoErrAvg] = calculateErrorAverage(data.impErrP_array, data.shoErrP_array, errSize);

    % tempo
    predTime = data.predTime;
    T_ref = data.T_ref;

    t = T_ref * (1:ceil((predTime(2) - predTime(1))/T_ref)) + 1;
    t_cell = cell(errSize)';
    for i = 1:errSize(2)
        t_cell{i, 1} = t;
    end
    for i = 2:errSize(1)
        t_cell(:, i) = t_cell(:, 1);
    end

    % labels
    labelsPKF = data.labelsPKF;


    % comparando P
    for i = 1:errSize(1)
        plotTitle(i) = labelsPKF.method(i);
    end
    for i = 1:errSize(2)
        plotLegend(i) = labelsPKF.P(i);
    end
    lines = strings(1, errSize(2));
    lim = [0 200];
    errorPlot(t_cell, impErrAvg, plotTitle, plotLegend, predTime, lim, lines);
    sgtitle("Erro de Impacto - Trajetória Simulada sem Arrasto - " + numOfSim + " Simulações - T = " + labelsPKF.T);
    errorPlot(t_cell, shoErrAvg, plotTitle, plotLegend, predTime, lim, lines);
    sgtitle("Erro de Disparo - Trajetória Simulada sem Arrasto - " + numOfSim + " Simulações - T = " + labelsPKF.T);
end