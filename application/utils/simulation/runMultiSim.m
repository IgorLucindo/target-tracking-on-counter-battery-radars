% rodar simulacoes com P fixo
function [impErr, shoErr] = runMultiSim(P_array, T_array, T_new, y, impPt, shoPt, predTime, R, model, filter, isInterp)
    % variaveis para multiplas simulacoes
    numOfT = length(T_array);
    numOfP = size(P_array, 3);
    method = [1 2 3];
    numOfMethods = 3;

    % saida medida
    y_cell = cell(numOfT);
    for i = 1:numOfT
        y_cell{i} = y(:, i:i:end);
    end

    % escolhe comparacao de T ou P
    errLength = max(numOfT, numOfP);

    % iniciar errCell
    impErr = cell(numOfMethods, errLength);
    shoErr = cell(numOfMethods, errLength);

    % simulacoes
    for i = 1:numOfMethods
        for j = 1:errLength
            [~, impErr{i, j}, shoErr{i, j}] = ...
                runSim(P_array(:, :, min(numOfP, j)), T_array(min(numOfT, j)), T_new, method(i), y_cell{min(numOfT, j)}, impPt, shoPt, predTime, R, model, filter, isInterp);
        end
    end
end