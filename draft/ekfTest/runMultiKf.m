% rodar kalman filter do inicio ate o fim e do fim ate o inicio,
% k vezes, para cada iteracao
function [kf, kf_rev, y_est, x_est] = runMultiKf(numOfFiltering, y, u, kf, kf_rev, P, i)
    y_est = zeros(3, i);
    P_77 = 1e-11/1.00^(i-20);
    for k = 1:numOfFiltering
        % filtragem direta
        if i == 60
            kf.P = P;
        else
            kf.P = kf_rev.P;
        end
        kf.P(4:6, 1:3) = -kf.P(4:6, 1:3);
        kf.P(1:3, 4:6) = -kf.P(1:3, 4:6);
        kf.P(4:6, 7) = -kf.P(4:6, 7);
        kf.P(7, 4:6) = -kf.P(7, 4:6);
        kf.P(7, 7) = P_77/1.2^(2*k-2);
        for j = 2:i
            [kf, y_est(:, j), ~] = kf.run(y(:, j), u);
        end
        [kf_rev, ~, ~] = kf_rev.setInitialState(kf.x);

        % filtragem reversa
        if i == 60
            kf_rev.P = P;
        else
            kf_rev.P = kf.P;
        end
        kf_rev.P(4:6, 1:3) = -kf_rev.P(4:6, 1:3);
        kf_rev.P(1:3, 4:6) = -kf_rev.P(1:3, 4:6);
        kf_rev.P(4:6, 7) = -kf_rev.P(4:6, 7);
        kf_rev.P(7, 4:6) = -kf_rev.P(7, 4:6);
        kf_rev.P(7, 7) = P_77/1.2^(2*k-1);
        for j = 2:i
            l = i - j + 1;
            [kf_rev, y_est(:, l), ~] = kf_rev.run(y(:, l), u);
        end
        [kf, ~, x_est] = kf.setInitialState(kf_rev.x);
    end
end