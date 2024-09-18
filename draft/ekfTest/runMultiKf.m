% rodar kalman filter do inicio ate o fim e do fim ate o inicio,
% n vezes, para cada iteracao
function [ekf, ekf_rev, y_est, x_est] = runMultiKf(numOfFiltering, y, u, P, ekf, ekf_rev, i)
    % resetando valor de P (nao adequado)
    % ekf.P = P;
    
    y_est = zeros(3, i);
    for n = 1:numOfFiltering
        % filtragem direta
        for j = 2:i
            [ekf, y_est(:, j), x_est] = ekf.run(y(:, j), u);
        end
        ekf_rev.x = ekf.x;
        ekf_rev.P = ekf.P;
        if i < 10
            ekf.P
        end

        % filtragem reversa
        for j = 2:i
            l = i - j + 1;
            [ekf_rev, y_est(:, l), x_est] = ekf_rev.run(y(:, l), u);
        end
        ekf.x = ekf_rev.x;
        ekf.P = ekf_rev.P;
    end
end