% prever ponto de impacto e disparo de modelo nao linear
function [impPtPred, shoPtPred] = predImpShoPt(x, u, p_floor, plant, plantRev, filter)
    if filter == "kf"
        gama = 1e-4;
        % gama = 0; % gama real inicial de 0
        % gama = 8.2873e-05; % gama real inicial de 1e-4
        % gama = 1.6797e-04; % gama real inicial de 2e-4
        % gama = 3.4314e-04; % gama real inicial de 4e-4
        x = [x; gama];
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