% plotar graficos de kf para comparacao de interpolacao
function plotInterpKf()
    % retornar se arquivo de resultados nao existir
    resultsFilePath = "results\weights_interp_kf.mat";
    if exist(resultsFilePath, 'file') == 0
        return
    end

    % carregar resultados
    data = load(resultsFilePath);

    % medias de erros
    errTiSize = size(data.impErrTi_array{1});
    numOfSim = size(data.impErrTi_array, 2);
    [impErrTiAvg, shoErrTiAvg] = calculateErrorAverage(data.impErrTi_array, data.shoErrTi_array, errTiSize);
    [impErrTAvg, shoErrTAvg] = calculateErrorAverage(data.impErrT_array, data.shoErrT_array, errTiSize);
    impErrTAvg = [impErrTAvg(1, 2:end); impErrTiAvg(1, 2:end)]';
    shoErrTAvg = [shoErrTAvg(1, 2:end); shoErrTiAvg(1, 2:end)]';
    errSize = size(impErrTAvg);

    % tempo
    predTime = data.predTime;
    T_ref = data.T_ref;

    t_cell = cell(errSize)';
    for i = 1:errSize(1)
        t_cell{1, i} = T_ref * (i+1) * (1:ceil((predTime(2) - predTime(1))/T_ref/(i+1))) + 1;
    end
    for i = 1:errSize(1)
        t_cell{2, i} = T_ref * (1:ceil((predTime(2) - predTime(1))/T_ref)) + 1;
    end


    % labels
    labels = data.labelsT;


    % comparando mesmo T
    for i = 1:errSize(1)
        plotTitle(i) = "T = " + labels.T(i+1);
    end
    plotLegend(1) = "Não Interpolado";
    plotLegend(2) = "Interpolado para " + labels.T(1);
    lim = [0 200];
    lines = ["b" "r--"];
    errorPlot(t_cell, impErrTAvg, plotTitle, plotLegend, predTime, lim, lines);
    sgtitle("Erro do Ponto de Impacto - Trajetória Simulada sem Arrasto - " + numOfSim + " Simulações");
    errorPlot(t_cell, shoErrTAvg, plotTitle, plotLegend, predTime, lim, lines);
    sgtitle("Erro do Ponto de Disparo - Trajetória Simulada sem Arrasto - " + numOfSim + " Simulações");
end