% rodar simulacao
function [y_est, impErr, shoErr] = runSim(P, T, T_new, method, y, impPt, shoPt, predTime, R, model, filter, isInterp)
    % checar interpolacao
    if ~isInterp
        T_new = T;
    end
    % rescala de periodo de amostragem interpolado
    rescale = T/T_new;

    % gama inicial para filtragem
    gama = 1.242e-4; % +50% erro, gama real inicial de 1e-4
    % gama = 4.1435e-5; % -50% erro, gama real inicial de 1e-4
    % gama = 1.657e-4; % +100% erro, gama real inicial de 1e-4
    % gama = 2.486e-4; % +200% erro, gama real inicial de 1e-4
    % roubando
    % gama = 0; % gama real inicial de 0
    % gama = 8.2873e-05; % gama real inicial de 1e-4
    % gama = 1.6797e-04; % gama real inicial de 2e-4
    % gama = 3.4314e-04; % gama real inicial de 4e-4

    i_change_method = 3;
    % metodo escolhido
    method_ = method;


    % matriz covariancia do ruido de processo
    Q = zeros(length(P));

    % checar tipo de modelo
    switch filter
        case "kf"
            % matrizes de espaco de estados de KF
            [A, B, C, A_rev, B_rev] = getParamsKf(T_new);

            % Kalman filter
            kf = KalmanFilter(A, B, C, Q, R, P);
            kf_rev = KalmanFilter(A_rev, B_rev, C, Q, R, P);

            % estado inicial para filtragem
            x0 = [zeros(5, 1); 200];

        case "ekf"
            % pegar funcoes de EKF
            [f, h, F, H, f_rev, F_rev] = getParamsEkf(T_new);

            % Kalman filter
            kf = ExtendedKalmanFilter(f, h, F, H, Q, R, P);
            kf_rev = ExtendedKalmanFilter(f_rev, h, F_rev, H, Q, R, P);

            % estado inicial para filtragem
            x0 = [zeros(5, 1); 200; gama];

            % metodo inicial
            method_ = 1;

            % amostra para troca de metodo
            % evita P explodir no caso dos metodos 2 e 3
            changeMethodTime = 3;
            i_change_method = floor(changeMethodTime/T);
    end

    % params das plantas
    T_100ms = 0.1;
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

        % trocar de metodo para P nao explodir no caso do ekf
        if i == i_change_method && method ~= 1
            method_ = method;
            kf_rev.P = P;
            if filter == "ekf"
                kf_rev.P(7, 7) = 1e-11;
            end
            [kf_rev, ~, ~] = kf_rev.setInitialState(kf.x);
            for k = 2:i
                l = i - k + 1;
                [kf_rev, y_est(:, l), ~] = kf_rev.run(y(:, l), u);
            end
            [kf, ~, x_est] = kf.setInitialState(kf_rev.x);
        end

        % pular caso amostras insuficientes para interpolacao
        % isso nao altera calculos sem interpolacao
        if i < 3
            continue
        end

        % obtem saida interpolada
        y_new = dynamicInterp(y, y_new, T, T_new, i);

        % kalman filter
        % filtragem das primeiras amostras
        if i == 3
            x0(1:3) = y_new(:, 1);
            % definir estado inicial
            [kf, y_est(:, 1), x_est] = kf.setInitialState(x0);
            for k = 2:rescale + 1
                switch method_
                    case 1
                        [kf, y_est(:, k), x_est] = kf.run(y_new(:, k), u);
                    case 2
                        numOfFiltering = 1;
                        [kf, kf_rev, y_est(:, 1:k), x_est] = runMultiKf(numOfFiltering, i_change_method, filter, y_new, u, P, kf, kf_rev, k);
                    case 3
                        numOfFiltering = 5;
                        [kf, kf_rev, y_est(:, 1:k), x_est] = runMultiKf(numOfFiltering, i_change_method, filter, y_new, u, P, kf, kf_rev, k);
                end
            end
        end
        % filtragem das ultimas amostras
        for k = (i-2)*rescale + 2:(i-1)*rescale + 1
            switch method_
                case 1
                    [kf, y_est(:, k), x_est] = kf.run(y_new(:, k), u);
                case 2
                    numOfFiltering = 1;
                    [kf, kf_rev, y_est(:, 1:k), x_est] = runMultiKf(numOfFiltering, i_change_method, filter, y_new, u, P, kf, kf_rev, k);
                case 3
                    numOfFiltering = 5;
                    [kf, kf_rev, y_est(:, 1:k), x_est] = runMultiKf(numOfFiltering, i_change_method, filter, y_new, u, P, kf, kf_rev, k);
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
                switch model
                    case "linear"
                        [impPtPred, shoPtPred] = predImpShoPtLinear(x_est, g, p_floor);

                    case "nonLinear"
                        [impPtPred, shoPtPred] = predImpShoPt(x_est, u, p_floor, plant, plantRev, filter);
                end

                impPtPredArray(:, j) = impPtPred;
                shoPtPredArray(:, j) = shoPtPred;

                % calcular erro dos pontos de impacto e de disparo
                [impErr(j), shoErr(j)] = getImpactShootingError(impPtPred, shoPtPred, impPt, shoPt);
            end
        end
    end
end