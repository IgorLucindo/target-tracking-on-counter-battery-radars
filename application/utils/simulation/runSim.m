% rodar simulacao
function [y_est, impErr, shoErr] = runSim(T, T_new, method, y, impPt, shoPt, predTime, R, filter, isInterp)
    % checar interpolacao
    if ~isInterp
        T_new = T;
    end
    % rescala de periodo de amostragem interpolado
    rescale = T/T_new;

    % gama inicial para filtragem
    gama = 8.2873e-05;

    % checar tipo de modelo
    switch filter
        case "kf"
            % matrizes de espaco de estados de KF
            [A, B, C, A_rev, B_rev] = getParamsKf(T_new);

            % matriz covariancia do ruido de processo
            Q = zeros(6);

            % valor inicial de matriz covariancia de estado
            P = 1e8*eye(6);

            % Kalman filter
            kf = KalmanFilter(A, B, C, Q, R, P);
            kf_rev = KalmanFilter(A_rev, B_rev, C, Q, R, P);

            % estado inicial para filtragem
            x0 = [zeros(5, 1); 200];

        case "ekf"
            % pegar funcoes de EKF
            [f, h, F, H, f_rev, F_rev] = getParamsEkf(T_new);

            % matriz covariancia do ruido de processo
            Q = zeros(7);

            % valor inicial de matriz covariancia de estado
            P = [1e3*eye(3) zeros(3) zeros(3, 1);
                 zeros(3) 1e2*eye(3) zeros(3, 1)
                 zeros(1, 3) zeros(1, 3) 0];

            % Kalman filter
            kf = ExtendedKalmanFilter(f, h, F, H, Q, R, P);
            kf_rev = ExtendedKalmanFilter(f_rev, h, F_rev, H, Q, R, P);

            % estado inicial para filtragem
            x0 = [zeros(5, 1); 200; gama];
    end

    % params das plantas
    T_100ms = 0.05;
    [f, h, ~, ~, f_rev, ~] = getParamsEkf(T_100ms);
    % planta da modelagem para obter ponto de impacto e disparo
    plant = Plant(f, h, zeros(7, 1));
    plantRev = Plant(f_rev, h, zeros(7, 1));

    % chao
    p_floor = 0;
    % gravidade
    g = 9.81;
    % entrada do sistema
    u = [0 0 -g]';
    
    arrayLength = ceil(predTime(2)/T_new);
    ySize = [3, arrayLength];
    % saida estimada
    y_est = zeros(ySize);
    % saida interpolada
    y_new = zeros(ySize);
    if T == T_new
        y_new = y;
    end

    arrayLength = ceil((predTime(2) - predTime(1))/T_new);
    % erros
    impErr = zeros(1, arrayLength);
    shoErr = zeros(1, arrayLength);
    
    impPtPredArray = zeros(ySize);
    shoPtPredArray = zeros(ySize);

    i = 0; j = 0;
    % loop
    while 1
        i = i + 1;
        currentTime = i*T;

        % parar quando tempo de execucao estiver fora do intervalo de predTime
        if currentTime > predTime(2)
            y_est = y_est(:, 1:end - 1);
            
            return
        end

        % pular caso amostras insuficientes
        if i < 3
            continue
        end

        % obtem saida interpolada
        y_new = dynamicInterp(y, y_new, T, T_new, i);

        % kalman filter
        % metodo 1 - roda kalman filter para cada iteracao
        % metodo 2 - rodar kalman filter, ida e volta, para cada iteracao
        % metodo 3 - rodar kalman filter, ida e volta, multiplas vezes para cada iteracao
        % filtragem das primeiras amostras
        if i == 3
            x0(1:3) = y_new(:, 1);
            % definir estado inicial
            [kf, y_est(:, 1), x_est] = kf.setInitialState(x0);
            for k = 2:rescale + 1
                switch method
                    case 1
                        [kf, y_est(:, k), x_est] = kf.run(y_new(:, k), u);
                    case 2
                        numOfFiltering = 1;
                        [kf, kf_rev, y_est(:, 1:k), x_est] = runMultiKf(numOfFiltering, y_new, u, P, kf, kf_rev, k);
                    case 3
                        numOfFiltering = 5;
                        [kf, kf_rev, y_est(:, 1:k), x_est] = runMultiKf(numOfFiltering, y_new, u, P, kf, kf_rev, k);
                end
            end
        end
        % filtragem das ultimas amostras
        for k = (i-2)*rescale + 2:(i-1)*rescale + 1
            switch method
                case 1
                    [kf, y_est(:, k), x_est] = kf.run(y_new(:, k), u);
                case 2
                    numOfFiltering = 1;
                    [kf, kf_rev, y_est(:, 1:k), x_est] = runMultiKf(numOfFiltering, y_new, u, P, kf, kf_rev, k);
                case 3
                    numOfFiltering = 5;
                    [kf, kf_rev, y_est(:, 1:k), x_est] = runMultiKf(numOfFiltering, y_new, u, P, kf, kf_rev, k);
            end
        end

        % previsao da trasdajetoria de impacto e disparo no intervalo de predTime
        if currentTime > predTime(1) && currentTime <= predTime(2)
            for k = (i-2)*rescale + 2:(i-1)*rescale + 1
                if j == arrayLength
                    continue
                end
                j = j + 1;

                % prever ponto de impacto e disparo
                % [impPtPred, shoPtPred] = predImpShoPtLinear(x_est, g, p_floor);
                [impPtPred, shoPtPred] = predImpShoPt(x_est, u, p_floor, plant, plantRev, filter);

                impPtPredArray(:, j) = impPtPred;
                shoPtPredArray(:, j) = shoPtPred;

                % calcular erro dos pontos de impacto e de disparo
                [impErr(j), shoErr(j)] = getImpactShootingError(impPtPred, shoPtPred, impPt, shoPt);
            end
        end
    end
end