% funcao nao usada mais
% previsao da trajetoria de impacto
function [y_pred] = setImpactTrajectory(plant_pred, y_pred, u, pos_floor)
    j = 0;
    while 1
        j = j + 1;
        [plant_pred, y_pred(:, j)] = plant_pred.run(u);
        if y_pred(3, j) < pos_floor
            y_pred = y_pred(:, 1:j);
            return
        end
    end
end