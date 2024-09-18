% Calcular desvio padrao de celulas
function [impErrStd, shoErrStd] = calculateErrorStd(impErr_array, shoErr_array, impErrAvg, shoErrAvg, errSize)
    % Obtém o número de simulações
    numOfSim = size(impErr_array, 2);

    % Inicializa as células para armazenar as médias
    impErrStd = cell(errSize);
    shoErrStd = cell(errSize);

    % Calcula o desvio padrão
    for i = 1:errSize(1)
        for j = 1:errSize(2)
            impErrSumSqDiff = 0;
            shoErrSumSqDiff = 0;
            for k = 1:numOfSim
                % Soma dos quadrados das diferenças
                impErrSumSqDiff = impErrSumSqDiff + (impErr_array{k}{i, j} - impErrAvg{i, j}).^2;
                shoErrSumSqDiff = shoErrSumSqDiff + (shoErr_array{k}{i, j} - shoErrAvg{i, j}).^2;
            end
            % Calcula o desvio padrão
            impErrStd{i, j} = sqrt(impErrSumSqDiff / numOfSim);
            shoErrStd{i, j} = sqrt(shoErrSumSqDiff / numOfSim);
        end
    end
end