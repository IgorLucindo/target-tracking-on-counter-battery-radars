% criar weights caso nao exista e definir labels
function createWeights(resultsFilePath, T_array, P_array_kf, P_array_ekf, predTime)
    numOfMethods = 3;
    method = [1 2 3];

    numOfP = size(P_array_kf, 3);
    numOfT = length(T_array);

    T_ref = T_array(1);

    data.labelsPKF.T = T_ref/0.001 + "ms";
    data.labelsPKF.method = strings(numOfMethods);
    data.labelsPKF.P = strings(numOfP);

    data.labelsT.T = strings(numOfT);
    data.labelsT.method = strings(numOfMethods);

    data.impErrP_array = {}; data.shoErrP_array = {};
    data.impErrPEKF_array = {}; data.shoErrPEKF_array = {};
    data.impErrT_array = {}; data.shoErrT_array = {};
    data.impErrTEKF_array= {}; data.shoErrTEKF_array = {};
    data.impErrTi_array = {}; data.shoErrTi_array = {};
    data.impErrTiEKF_array = {}; data.shoErrTiEKF_array = {};
    data.impErrRad_array = {}; data.shoErrRad_array = {};
    data.impErrRadEKF_array = {}; data.shoErrRadEKF_array = {};
    data.impErrRadTiEKF_array = {}; data.shoErrRadTiEKF_array = {};

    data.numOfSimP_kf = 0; data.numOfSimP_ekf = 0;
    data.numOfSimT_kf = 0; data.numOfSimT_ekf = 0;
    data.numOfSimTi_kf = 0; data.numOfSimTi_ekf = 0;
    data.numOfSimRad_kf = 0; data.numOfSimRad_ekf = 0;

    for i = 1:numOfMethods
        data.labelsPKF.method(i) = "Método " + method(i);
        data.labelsT.method(i) = "Método " + method(i);
    end
    for i = 1:numOfP
        data.labelsPKF.P(i) = "P = 1e" + log10(P_array_kf(1, 1, i)) + " I_{6}";
        data.labelsPEKF.P(i) = "P(7, 7) = 1e" + log10(P_array_ekf(7, 7, i));
    end
    for i = 1:numOfT
        data.labelsT.T(i) = T_array(i)/0.001 + "ms";
    end

    data.T_ref = T_ref;
    data.predTime = predTime;

    save(resultsFilePath, '-struct', 'data');
end