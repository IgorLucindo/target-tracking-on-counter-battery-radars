% prever ponto de impacto e disparo de modelo nao linear
function [impPtPred, shoPtPred] = predImpShoPt(x, u, p_floor, plant, plantRev, modelType)
    if modelType == "kf"
        x = [x; 1e-4];
    end

    % ponto de impacto e disparo
    impPtPred = zeros(3, 1);
    shoPtPred = zeros(3, 1);

    % definir estado inicial
    plant = plant.setInitialState(x);
    plantRev = plantRev.setInitialState(x);

    % previsao
    while 1
        [plant, impPtPred] = plant.run(u);
        if impPtPred(3) < p_floor
            break
        end
    end
    while 1
        [plantRev, shoPtPred] = plantRev.run(u);
        if shoPtPred(3) < p_floor
            return
        end
    end
end