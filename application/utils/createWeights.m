% criar weights caso nao exista e definir labels
function createWeights(resultsFilePath, T_array, T_P, P_array, predTime)
    numOfMethods = 3;
    method = [1 2 3];

    numOfP = size(P_array, 3);
    numOfT = length(T_array);

    data.labelsP.T = T_P/0.001 + "ms";
    data.labelsP.method = strings(numOfMethods);
    data.labelsP.P = strings(numOfP);

    data.labelsT.T = strings(numOfT);
    data.labelsT.method = strings(numOfMethods);

    data.impErrP_array = {}; data.shoErrP_array = {};
    data.impErrT_array = {}; data.shoErrT_array = {};
    data.impErrTEKF_array= {}; data.shoErrTEKF_array = {};
    data.impErrTi_array = {}; data.shoErrTi_array = {};
    data.impErrTiEKF_array = {}; data.shoErrTiEKF_array = {};
    data.impErrRad_array = {}; data.shoErrRad_array = {};
    data.impErrRadEKF_array = {}; data.shoErrRadEKF_array = {};
    data.impErrRadTiEKF_array = {}; data.shoErrRadTiEKF_array = {};

    data.numOfSimP = 0;
    data.numOfSimT = 0;
    data.numOfSimTi = 0;
    data.numOfSimRad = 0;

    for i = 1:numOfMethods
        data.labelsP.method(i) = "metodo " + method(i);
        data.labelsT.method(i) = "metodo " + method(i);
    end
    for i = 1:numOfP
        data.labelsP.P(i) = "P = 1e" + log10(P_array(1, 1, i)) + " I_{6}";
    end
    for i = 1:numOfT
        data.labelsT.T(i) = T_array(i)/0.001 + "ms";
    end

    data.T_ref = T_array(1);
    data.predTime = predTime;

    save(resultsFilePath, '-struct', 'data');
end