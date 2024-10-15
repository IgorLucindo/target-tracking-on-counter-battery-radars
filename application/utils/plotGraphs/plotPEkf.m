% plotar graficos de ekf para comparacao de diferentes valores de P
function plotPEkf()
    % retornar se arquivo de resultados nao existir
    resultsFilePath = "results\weights_P_ekf.mat";
    if exist(resultsFilePath, 'file') == 0
        return
    end

    % carregar resultados
    data = load(resultsFilePath);

    % medias de erros
    errSize = size(data.impErrPEKF_array{1});
    numOfSim = size(data.impErrPEKF_array, 2);
    [impErrAvg, shoErrAvg] = calculateErrorAverage(data.impErrPEKF_array, data.shoErrPEKF_array, errSize);

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
    labelsPEKF = data.labelsPEKF;


    % comparando P
    for i = 1:errSize(1)
        plotTitle(i) = labelsPKF.method(i);
    end
    for i = 1:errSize(2)
        plotLegend(i) = labelsPEKF.P(i);
    end
    lines = strings(1, errSize(2));
    errorPlot(t_cell, impErrAvg, plotTitle, plotLegend, predTime, [150 1200], lines);
    sgtitle("Erro do Ponto de Impacto - Trajetória Simulada com Arrasto - " + numOfSim + " Simulações - T = " + labelsPKF.T + " - γ_{inicial} = 200% γ_{real}");
    errorPlot(t_cell, shoErrAvg, plotTitle, plotLegend, predTime, [0 200], lines);
    sgtitle("Erro do Ponto de Disparo - Trajetória Simulada com Arrasto - " + numOfSim + " Simulações - T = " + labelsPKF.T + " - γ_{inicial} = 200% γ_{real}");
end