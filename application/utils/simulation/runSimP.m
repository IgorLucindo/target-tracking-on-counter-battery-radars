% rodar simulacoes comparando P
function [impErrP, shoErrP] = runSimP(T, P_array, y, g, u, impPt, shoPt, p_floor, predTime, R, plant, plantRev)
    numOfMethods = 3;
    numOfP = size(P_array, 3);

    % metodo
    method = [1 2 3];
    
    % iniciar errCell
    impErrP = cell(numOfMethods, numOfP);
    shoErrP = cell(numOfMethods, numOfP);

    for i = 1:numOfMethods
        for j = 1:numOfP
            [~, impErrP{i, j}, shoErrP{i, j}] = runSim(T, method(i), P_array(:, :, j), y, g, u, impPt, shoPt, p_floor, predTime, R);
        end
    end
end