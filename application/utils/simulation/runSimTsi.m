% rodar simulacoes com P fixo
function [impErrTs, shoErrTs] = runSimTsi(Ts_array, y, g, u, impPt, shoPt, p_floor, predTime, Q, R)
    numOfMethods = 3;
    numOfTs = length(Ts_array);

    % metodo
    method = [1 2 3];

    % erro de covariancia
    P = 1e8*eye(6);

    % saida medida
    y_cell = cell(numOfTs);
    for i = 1:numOfTs
        y_cell{i} = y(:, i:i:end);
    end

    % iniciar errCell
    impErrTs = cell(numOfMethods, numOfTs);
    shoErrTs = cell(numOfMethods, numOfTs);

    for i = 1:numOfMethods
        for j = 1:numOfTs
            [~, impErrTs{i, j}, shoErrTs{i, j}] = runSimInterp(Ts_array(j), Ts_array(1), method(i), P, y_cell{j}, g, u, impPt, shoPt, p_floor, predTime, Q, R);
        end
    end
end