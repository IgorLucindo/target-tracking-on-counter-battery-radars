% rodar simulacao
function [y_est, impErr, shoErr] = runSim(Ts, method, P, y, g, u, impPt, shoPt, p_floor, predTime, Q, R)
    % matrizes de espaco de estados
    A = eye(6);
    A(1:3, 4:end) = Ts*eye(3);
    B = [Ts^2/2*eye(3); Ts*eye(3)];
    C = [eye(3) zeros(3)];
    A_rev = inv(A);
    B_rev = -A_rev*B;
    
    % Kalman filter
    kf = KalmanFilter(A, B, C, Q, R, P);
    kf_rev = KalmanFilter(A_rev, B_rev, C, Q, R, P);
    
    arrayLength = ceil(predTime(2)/Ts);
    % estado estimado
    x_est = zeros(length(A), arrayLength);
    
    ySize = [size(C, 1), arrayLength];
    % saida estimada
    y_est = zeros(ySize);

    arrayLength = ceil((predTime(2) - predTime(1))/Ts);
    % erros
    impErr = zeros(1, arrayLength);
    shoErr = zeros(1, arrayLength);
    
    impPtPredArray = zeros(ySize);
    shoPtPredArray = zeros(ySize);

    i = 0;
    j = 0;
    % loop
    while 1
        i = i + 1;
        currentTime = i*Ts;

        % parar quando tempo de execucao estiver fora do intervalo de predTime
        if currentTime > predTime(2)
            y_est = y_est(:, 1:i - 1);
            return
        end

        % kalman filter
        % primeira iteracao defini estado inicial de kalman filter
        if i == 1
            [kf, y_est(:, 1), x_est(:, 1)] = kf.setX0(y(:, 1));
            continue;
        end
        % metodo 1 - roda kalman filter para cada iteracao
        % metodo 2 - rodar kalman filter, ida e volta, para cada iteracao
        % metodo 3 - rodar kalman filter, ida e volta, multiplas vezes para cada iteracao
        switch method
            case 1
                [kf, y_est(:, i), x_est(:, i)] = kf.run(y(:, i), u);
            case 2
                numOfFiltering = 1;
                [kf, kf_rev, y_est, x_est] = runMultiKf(numOfFiltering, y_est, x_est, y, u, kf, kf_rev, i);
            case 3
                numOfFiltering = 5;
                [kf, kf_rev, y_est, x_est] = runMultiKf(numOfFiltering, y_est, x_est, y, u, kf, kf_rev, i);
        end

        % previsao da trasdajetoria de impacto e disparo no intervalo de predTime
        if currentTime > predTime(1) && currentTime <= predTime(2)
            j = j + 1;
            [impPtPred, shoPtPred] = setImpactShootingPoint(x_est(:, i), g, p_floor);
            impPtPredArray(:, j) = impPtPred;
            shoPtPredArray(:, j) = shoPtPred;
            % calcular erro dos pontos de impacto e de disparo
            [impErr(j), shoErr(j)] = getImpactShootingError(impPtPred, shoPtPred, impPt, shoPt);
        end
    end
end