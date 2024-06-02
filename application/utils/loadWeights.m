% carregar erros anteriores salvos em weights
function [impErrAvgArray, shoErrAvgArray, numOfSimTotal] = loadWeights(weights_path, errArraySize)
    % valores padroes se arquivo nao existir
    if exist(weights_path, 'file') == 0
        impErrAvgArray = zeros(errArraySize);
        shoErrAvgArray = zeros(errArraySize);
        numOfSimTotal = 1;
        return
    end
    
    % carregar
    data = load(weights_path);
    impErrAvgArray = data.impErrAvgArray;
    shoErrAvgArray = data.shoErrAvgArray;
    numOfSimTotal = data.numOfSimTotal;
end