% definir saida real
function [y_true, impPt, shoPt] = setTrueTrajectory(x_true, T, u, p_floor)
    % periodo de amostragem pequeno para melhor resolucao de trajetoria
    T_5ms = 0.005;

    % planta
    [f, h, ~, ~, f_rev, ~] = getParamsEkf(T_5ms);
    plant = Plant(f, h, x_true);
    plantRev = Plant(f_rev, h, x_true);

    % saida real
    y_true = [];

    % ponto de impacto e disparo
    impPt = zeros(3, 1);
    shoPt = zeros(3, 1);

    while 1
        % planta
        [plant, y_true_i] = plant.run(u);
        y_true = [y_true y_true_i];
        
        % parar no impacto
        if y_true_i(3) < p_floor
            rescale = T/T_5ms;
            y_true = y_true(:, rescale:rescale:end);
            impPt = y_true_i;
            break;
        end
    end

    % prever ponto de disparo
    while 1
        [plantRev, shoPt] = plantRev.run(u);
        if shoPt(3) < p_floor
            return
        end
    end
end