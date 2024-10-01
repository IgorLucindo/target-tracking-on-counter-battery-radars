% rodar kalman filter do inicio ate o fim e do fim ate o inicio,
% k vezes, para cada iteracao
function [kf, kf_rev, y_est, x_est] = runMultiKf(numOfFiltering, i_change_method, filter, y, u, P, kf, kf_rev, i)
    y_est = zeros(3, i);
    P_77 = 1e-11;
    
    for k = 1:numOfFiltering
        % filtragem direta
        if i == i_change_method
            kf.P = P;
        else
            switch filter
                case "kf"
                    kf.P = P;
                case "ekf"
                    kf.P = kf_rev.P;
            end
        end
        kf.P(4:6, 1:3) = -kf.P(4:6, 1:3);
        kf.P(1:3, 4:6) = -kf.P(1:3, 4:6);
        if filter == "ekf"
            kf.P(4:6, 7) = -kf.P(4:6, 7);
            kf.P(7, 4:6) = -kf.P(7, 4:6);
            kf.P(7, 7) = P_77/1.2^(2*k-2);
        end
        % ida
        for j = 2:i
            [kf, y_est(:, j), ~] = kf.run(y(:, j), u);
        end
        [kf_rev, ~, ~] = kf_rev.setInitialState(kf.x);

        % filtragem reversa
        if i == i_change_method
            kf_rev.P = P;
        else
            switch filter
                case "kf"
                    kf_rev.P = P;
                case "ekf"
                    kf_rev.P = kf.P;
            end
        end
        kf_rev.P(4:6, 1:3) = -kf_rev.P(4:6, 1:3);
        kf_rev.P(1:3, 4:6) = -kf_rev.P(1:3, 4:6);
        if filter == "ekf"
            kf_rev.P(4:6, 7) = -kf_rev.P(4:6, 7);
            kf_rev.P(7, 4:6) = -kf_rev.P(7, 4:6);
            kf_rev.P(7, 7) = P_77/1.2^(2*k-1);
        end
        % volta
        for j = 2:i
            l = i - j + 1;
            [kf_rev, y_est(:, l), ~] = kf_rev.run(y(:, l), u);
        end
        [kf, ~, x_est] = kf.setInitialState(kf_rev.x);
    end
end