% definir P_array para comparacoes
function [P_array_kf, P_array_ekf] = setPArray(numOfP)
    P_array_kf = zeros(6, 6, numOfP);
    P_array_ekf = zeros(7, 7, numOfP);

    % se o tamanho for 1
    if numOfP == 1
        P_array_kf(:, :, 1) = 1e8*eye(6);
        P_array_ekf(:, :, 1) = 1e8*eye(7);
        P_array_ekf(7, 7, 1) = 1e-10;

        return;
    end

    % caso kf
    for i = 1:numOfP
        P_array_kf(:, :, i) = 10^i*eye(6);
    end

    % caso ekf
    for i = 1:numOfP
        P_array_ekf(:, :, i) = 1e8*eye(7);
        P_array_ekf(7, 7, i) = 10^-(7 + i);
    end
end