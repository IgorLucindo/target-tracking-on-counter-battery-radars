% rodar simulacoes com P fixo
function [impErrT, shoErrT] = runMultiSim(T_array, T_new, y, impPt, shoPt, predTime, R, filter, isInterp)
    numOfMethods = 3;
    numOfT = length(T_array);

    % metodo
    method = [1 2 3];

    % saida medida
    y_cell = cell(numOfT);
    for i = 1:numOfT
        y_cell{i} = y(:, i:i:end);
    end

    % iniciar errCell
    impErrT = cell(numOfMethods, numOfT);
    shoErrT = cell(numOfMethods, numOfT);

    for i = 1:numOfMethods
        for j = 1:numOfT
            [~, impErrT{i, j}, shoErrT{i, j}] = runSim(T_array(j), T_new, method(i), y_cell{j}, impPt, shoPt, predTime, R, filter, isInterp);
        end
    end
end