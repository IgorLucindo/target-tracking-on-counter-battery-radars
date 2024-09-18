% rodar simulacoes dependendo da comparacao escolhida
function runSimWithComparison(compare, numOfSim, predTime)
    % definir flags de comparacao e escolher path do arquivo de resultado
    switch compare
        case "models"
            simFlags = dictionary(["P" "T" "Ti" "Rad"], [0 1 0 0]);
            resultsFilePath = "results\weights_model.mat";
        otherwise
            fprintf("comparison unrecognized\n");
            return;
    end


    % parametros para comparacao
    numOfP = 4;
    numOfT = 1;
    T_ref = 0.05;
    model = "nonLinear";
    detecThresh = ceil(10/T_ref);

    % variaveis do modelo
    g = 9.81;
    v0 = [180 180 250];
    a = [0 0 -g];
    % entrada do sistema
    u = a';
    switch model
        case "linear"
            gama = 0;
        case "nonLinear"
            gama = 1e-4;
    end

    % chao
    p_floor = 0;

    % matriz covariancia do ruido de medicao
    sigma2_n = 1e2;
    R = sigma2_n*eye(3);   % roubo
    R = 0.1*R;             % R chutado para 10 vezes menor
    % quero testar R obtido por autocovariance least-square

    % periodo de amostragem
    T_array = T_ref * (1:numOfT);
    T_10ms = 0.01;


    % definir trajetoria, ponto de impacto e ponto de disparo reais
    x_true = [zeros(1, 3) v0 gama]';
    [y_true, impPt, shoPt] = setTrueTrajectory(x_true, T_ref, u, p_floor);

    % array de diferentes P
    P_array = zeros(6, 6, numOfP);
    for i = 1:numOfP
        P_array(:, :, i) = 10^(10 - 2*i) * eye(6);
    end


    % criar weights caso nao exista e definir labels
    if exist(resultsFilePath, 'file') == 0
        createWeights(resultsFilePath, T_array, T_10ms, P_array, predTime);
    end

    % carregar weights
    data = load(resultsFilePath);

    % numero de simulacoes anteriores
    keys = simFlags.keys;
    numOfPrevSim = zeros(1, length(keys));
    for j = 1:length(keys)
        key = keys(j);
        numOfPrevSim(j) = data.('numOfSim' + keys(j));
    end


    % simulacoes
    for i = 1:numOfSim
        tic;
        % printar mensagem de simulacao
        for j = 1:length(keys)
            key = keys(j);
            if ~simFlags(key)
                continue
            end
            fprintf("simulacao - " + key + ": " + (numOfPrevSim(j)+i) + "/" + (numOfPrevSim(j)+numOfSim) + "    ");
        end

        % pegar trajetoria
        for j = 1:length(keys)
            key = keys(j);
            if ~simFlags(key)
                continue
            end
            switch key
                case 'P'
                    y_P = getSimulatedTrajectory(y_true, sigma2_n, predTime, T_10ms, detecThresh);
                case 'Rad'
                    [y_Rad, impPt, shoPt, T_rad] = getRadarTrajectoryDynamically(i);
                    T_array = [T_rad];
                otherwise
                    y = getSimulatedTrajectory(y_true, sigma2_n, predTime, T_ref, detecThresh);
            end
        end

        % simulacao dependendo da comparacao
        for j = 1:length(keys)
            key = keys(j);
            if ~simFlags(key)
                continue
            end

            % rodar simulacao
            isInterp = false;
            switch key
                case 'P'
                    [data.impErrP_array{end + 1}, data.shoErrP_array{end + 1}] = ...
                        runSimP(T_10ms, P_array, y_P, g, u, impPt, shoPt, p_floor, predTime, R);

                case 'T'
                    [data.impErrT_array{end + 1}, data.shoErrT_array{end + 1}] = ...
                        runMultiSim(T_array, T_ref, y, impPt, shoPt, predTime, R, "kf", isInterp);
                    [data.impErrTEKF_array{end + 1}, data.shoErrTEKF_array{end + 1}] = ...
                        runMultiSim(T_array, T_ref, y, impPt, shoPt, predTime, R, "ekf", isInterp);

                case 'Ti'
                    isInterp = true;
                    [data.impErrTi_array{end + 1}, data.shoErrTi_array{end + 1}] = ...
                        runMultiSim(T_array, T_ref, y, impPt, shoPt, predTime, R, "kf", isInterp);
                    [data.impErrTiEKF_array{end + 1}, data.shoErrTiEKF_array{end + 1}] = ...
                        runMultiSim(T_array, T_ref, y, impPt, shoPt, predTime, R, "ekf", isInterp);

                case 'Rad'
                    [data.impErrRad_array{end + 1}, data.shoErrRad_array{end + 1}] = ...
                        runMultiSim(T_array, T_ref, y_Rad, impPt, shoPt, predTime, R, "kf", isInterp);
                    [data.impErrRadEKF_array{end + 1}, data.shoErrRadEKF_array{end + 1}] = ...
                        runMultiSim(T_array, T_ref, y_Rad, impPt, shoPt, predTime, R, "ekf", isInterp);
                    isInterp = true;
                    [data.impErrRadTiEKF_array{end + 1}, data.shoErrRadTiEKF_array{end + 1}] = ...
                        runMultiSim(T_array, 0.005, y_Rad, impPt, shoPt, predTime, R, "ekf", isInterp);
            end

            data.('numOfSim' + keys(j)) = numOfPrevSim(j) + i;
        end

        % salva erros em weights
        save(resultsFilePath, '-struct', 'data');
        toc;
        fprintf("\n");
    end
end
