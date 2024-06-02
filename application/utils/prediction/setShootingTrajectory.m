% funcao nao usada mais
% previsao da trajetoria de disparo
function [y_pred_rev] = setShootingTrajectory(plant_pred_rev, y_pred_rev, u, pos_floor)
    j = 0;
    while 1
        j = j + 1;
        [plant_pred_rev, y_pred_rev(:, j)] = plant_pred_rev.run(u);
        if y_pred_rev(3, j) < pos_floor
            y_pred_rev = y_pred_rev(:, 1:j);
            return
        end
    end
end