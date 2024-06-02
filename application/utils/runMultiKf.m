% rodar kalman filter do inicio ate o fim e do fim ate o inicio,
% n vezes, para cada iteracao
function [kf, kf_rev, y_est, x_est] = runMultiKf(numOfFiltering, y_est, x_est, y, u, kf, kf_rev, i)
    for n = 1:numOfFiltering
        kf.x = x_est(:, 1);
        for j = 2:i
            [kf, y_est(:, j), x_est(:, j)] = kf.run(y(:, j), u);
        end
        kf_rev.x = x_est(:, i);
        for j = 2:i
            l = i - j + 1;
            [kf_rev, y_est(:, l), x_est(:, l)] = kf_rev.run(y(:, l), u);
        end
    end
end