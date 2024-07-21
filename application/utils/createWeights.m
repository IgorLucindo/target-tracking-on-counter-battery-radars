% criar weights caso nao exista e definir labels
function createWeights(Ts_array, Ts_P, P_array, predTime)
    numOfMethods = 3;
    method = [1 2 3];

    numOfP = size(P_array, 3);
    numOfTs = length(Ts_array);

    data.labelsP.Ts = Ts_P/0.001 + "ms";
    data.labelsP.method = strings(numOfMethods);
    data.labelsP.P = strings(numOfP);

    data.labelsTs.Ts = strings(numOfTs);
    data.labelsTs.method = strings(numOfMethods);

    data.impErrP_array = {}; data.shoErrP_array = {};
    data.impErrTs_array = {}; data.shoErrTs_array = {};
    data.impErrTsi_array = {}; data.shoErrTsi_array = {};
    data.impErrRad_array = {}; data.shoErrRad_array = {};

    data.numOfSimP = 0;
    data.numOfSimTs = 0;
    data.numOfSimTsi = 0;
    data.numOfSimRad = 0;

    for i = 1:numOfMethods
        data.labelsP.method(i) = "metodo " + method(i);
        data.labelsTs.method(i) = "metodo " + method(i);
    end
    for i = 1:numOfP
        data.labelsP.P(i) = "P = 1e" + log10(P_array(1, 1, i)) + " I_{6}";
    end
    for i = 1:numOfTs
        data.labelsTs.Ts(i) = Ts_array(i)/0.001 + "ms";
    end

    data.Ts_ref = Ts_array(1);
    data.predTime = predTime;

    save('weights\weights.mat', '-struct', 'data');
end