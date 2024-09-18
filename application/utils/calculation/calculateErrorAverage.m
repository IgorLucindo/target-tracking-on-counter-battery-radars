% Calcular média de celulas
function [impErrAvg, shoErrAvg] = calculateErrorAverage(impErr_array, shoErr_array, errSize)
    % Obtém o número de simulações
    numOfSim = size(impErr_array, 2);

    % Inicializa as células para armazenar as médias
    impErrAvg = cell(errSize);
    shoErrAvg = cell(errSize);

    % Soma os valores de cada simulação e calcula a média
    for i = 1:errSize(1)
        for j = 1:errSize(2)
            impErrSum = 0;
            shoErrSum = 0;
            for k = 1:numOfSim
                impErrSum = impErrSum + impErr_array{k}{i, j};
                shoErrSum = shoErrSum + shoErr_array{k}{i, j};
            end
            impErrAvg{i, j} = impErrSum / numOfSim;
            shoErrAvg{i, j} = shoErrSum / numOfSim;
        end
    end
end