% rodar simulacao
function [y_est, impErr, shoErr] = runSimInterp(Ts, Ts_new, method, P, y, g, u, impPt, shoPt, p_floor, predTime, Q, R)
    % matrizes de espaco de estados
    A = eye(6);
    A(1:3, 4:end) = Ts_new*eye(3);
    B = [Ts_new^2/2*eye(3); Ts_new*eye(3)];
    C = [eye(3) zeros(3)];
    A_rev = inv(A);
    B_rev = -A_rev*B;
    A2 = eye(6);
    A2(1:3, 4:end) = Ts*eye(3);
    B2 = [Ts^2/2*eye(3); Ts*eye(3)];
    
    % Kalman filter
    kf = KalmanFilter(A, B, C, Q, R, P);
    kf_rev = KalmanFilter(A_rev, B_rev, C, Q, R, P);
    kf2 = KalmanFilter(A2, B2, C, Q, R, P);
    
    arrayLength = ceil(predTime(2)/Ts_new);
    % estado estimado
    x_est = zeros(length(A), arrayLength);
    
    ySize = [size(C, 1), arrayLength];
    % saida estimada
    y_est = zeros(ySize);
    % saida interpolada
    y_new = zeros(ySize);
    if Ts == Ts_new
        y_new = y;
    end

    arrayLength = ceil((predTime(2) - predTime(1))/Ts_new);
    % erros
    impErr = zeros(1, arrayLength);
    shoErr = zeros(1, arrayLength);
    
    impPtPredArray = zeros(ySize);
    shoPtPredArray = zeros(ySize);

    i = 0; j = 0;
    % loop
    while 1
        i = i + 1;
        currentTime = (i - 1)*Ts_new + Ts;

        % obtem saida interpolada
        [y_new, kf2] = dynamicInterpolation(y, y_new, Ts, Ts_new, kf2, u, i);

        % parar quando tempo de execucao estiver fora do intervalo de predTime
        if currentTime > predTime(2)
            y_est = y_est(:, 1:i - 1);
            return
        end

        % kalman filter
        % primeira iteracao defini estado inicial de kalman filter
        if currentTime < 2*Ts
            continue
        elseif currentTime == 2*Ts
            [kf, y_est(:, 1), x_est(:, 1)] = kf.setX0(y_new(:, 1));
            for l = 2:i
                switch method
                    case 1
                        [kf, y_est(:, l), x_est(:, l)] = kf.run(y_new(:, l), u);
                    case 2
                        numOfFiltering = 1;
                        [kf, kf_rev, y_est, x_est] = runMultiKf(numOfFiltering, y_est, x_est, y_new, u, kf, kf_rev, l);
                    case 3
                        numOfFiltering = 5;
                        [kf, kf_rev, y_est, x_est] = runMultiKf(numOfFiltering, y_est, x_est, y_new, u, kf, kf_rev, l);
                end
            end
            continue;
        end
        % metodo 1 - roda kalman filter para cada iteracao
        % metodo 2 - rodar kalman filter, ida e volta, para cada iteracao
        % metodo 3 - rodar kalman filter, ida e volta, multiplas vezes para cada iteracao
        switch method
            case 1
                [kf, y_est(:, i), x_est(:, i)] = kf.run(y_new(:, i), u);
            case 2
                numOfFiltering = 1;
                [kf, kf_rev, y_est, x_est] = runMultiKf(numOfFiltering, y_est, x_est, y_new, u, kf, kf_rev, i);
            case 3
                numOfFiltering = 5;
                [kf, kf_rev, y_est, x_est] = runMultiKf(numOfFiltering, y_est, x_est, y_new, u, kf, kf_rev, i);
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