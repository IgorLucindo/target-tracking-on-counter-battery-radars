% calcula media dos errCell
function [impErrAvg, shoErrAvg] = calculateErrorAverage(impErr_array, shoErr_array, errSize)
    numOfSim = size(impErr_array, 2);

    impErrAvg = cell(errSize);
    shoErrAvg = cell(errSize);

    for i = 1:errSize(1)
        for j = 1:errSize(2)
            impErrAvg{i, j} = impErr_array{1}{i, j};
            shoErrAvg{i, j} = shoErr_array{1}{i, j};
        end
    end
    for i = 1:errSize(1)
        for j = 1:errSize(2)
            for k = 2:numOfSim
                impErrAvg{i, j} = impErrAvg{i, j} + impErr_array{k}{i, j};
                shoErrAvg{i, j} = shoErrAvg{i, j} + shoErr_array{k}{i, j};
            end
            impErrAvg{i, j} = impErrAvg{i, j} / numOfSim;
            shoErrAvg{i, j} = shoErrAvg{i, j} / numOfSim;
        end
    end
end