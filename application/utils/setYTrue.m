% definir saida real
function [y_true, impPt, shoPt] = setYTrue(Ts, vel0, g, u, pos_floor)
    % matrizes de espaco de estados
    A = eye(6);
    A(1:3, 4:end) = Ts*eye(3);
    B = [Ts^2/2*eye(3); Ts*eye(3)];
    C = [eye(3) zeros(3)];

    % planta
    x_ss0 = [zeros(1, 3) vel0]';
    plant = Plant(A, B, C, x_ss0);

    % saida real
    y_true = [];

    % PI e PD
    impPt = zeros(3, 1);
    shoPt = zeros(3, 1);
    while 1
        % planta
        [plant, y_true_i] = plant.run(u);
        y_true = [y_true y_true_i];
        % parar no impacto
        if y_true_i(3) < pos_floor
            [impPt, shoPt] = setImpactShootingPoint(x_ss0, g, pos_floor);
            break;
        end
    end
end