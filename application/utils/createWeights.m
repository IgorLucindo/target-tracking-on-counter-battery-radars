% criar weights caso nao exista e definir labels
function createWeights(Ts_array, Ts_P, P_array, predTime)
    numOfMethods = 3;
    method = [1 2 3];

    numOfP = size(P_array, 3);
    numOfTs = length(Ts_array);

    labelsP.Ts = Ts_P/0.001 + "ms";
    labelsP.method = strings(numOfMethods);
    labelsP.P = strings(numOfP);

    labelsTs.Ts = strings(numOfTs);
    labelsTs.method = strings(numOfMethods);

    impErrP_array = {}; shoErrP_array = {};
    impErrTs_array = {}; shoErrTs_array = {};

    numOfSimP = 0; numOfSimTs = 0;

    for i = 1:numOfMethods
        labelsP.method(i) = "metodo " + method(i);
        labelsTs.method(i) = "metodo " + method(i);
    end
    for i = 1:numOfP
        labelsP.P(i) = "P = 1e" + log10(P_array(1, 1, i)) + " I_{6}";
    end
    for i = 1:numOfTs
        labelsTs.Ts(i) = Ts_array(i)/0.001 + "ms";
    end

    Ts_ref = Ts_array(1);

    save('weights\weights.mat', 'impErrP_array', 'impErrTs_array', 'shoErrP_array', 'shoErrTs_array', 'labelsP', 'labelsTs', 'numOfSimP', 'numOfSimTs', 'Ts_ref', 'predTime');
end