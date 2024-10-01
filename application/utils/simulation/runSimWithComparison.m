% rodar simulacoes dependendo da comparacao escolhida
function runSimWithComparison(compare, numOfSim, predTime)
    % definir parametros de simulacao a partir da comparacao escolhida
    [simFlags, resultsFilePath, model, T_ref, numOfT, numOfP] = setSimParams(compare);

    % parametros para comparacao
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

    % definir T_array para comparacoes
    T_array = T_ref * (1:numOfT);

    % definir P_array para comparacoes
    [P_array_kf, P_array_ekf] = setPArray(numOfP);

    % definir trajetoria, ponto de impacto e ponto de disparo reais
    x_true = [zeros(1, 3) v0 gama]';
    [y_true, impPt, shoPt] = setTrueTrajectory(x_true, T_ref, u, p_floor);


    % criar weights caso nao exista e definir labels
    if exist(resultsFilePath, 'file') == 0
        createWeights(resultsFilePath, T_array, P_array_kf, P_array_ekf, predTime);
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
                case {"Rad_kf", "Rad_ekf"}
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
            switch key
                case "P_kf"
                    [data.impErrP_array{end + 1}, data.shoErrP_array{end + 1}] = ...
                        runMultiSim(P_array_kf, T_array, T_ref, y, impPt, shoPt, predTime, R, model, "kf", false);

                case "P_ekf"
                    [data.impErrPEKF_array{end + 1}, data.shoErrPEKF_array{end + 1}] = ...
                        runMultiSim(P_array_ekf, T_array, T_ref, y, impPt, shoPt, predTime, R, model, "ekf", false);

                case "T_kf"
                    [data.impErrT_array{end + 1}, data.shoErrT_array{end + 1}] = ...
                        runMultiSim(P_array_kf, T_array, T_ref, y, impPt, shoPt, predTime, R, model, "kf", false);
                    
                case "T_ekf"
                    [data.impErrTEKF_array{end + 1}, data.shoErrTEKF_array{end + 1}] = ...
                        runMultiSim(P_array_ekf, T_array, T_ref, y, impPt, shoPt, predTime, R, model, "ekf", false);

                case "Ti_kf"
                    [data.impErrTi_array{end + 1}, data.shoErrTi_array{end + 1}] = ...
                        runMultiSim(P_array_kf, T_array, T_ref, y, impPt, shoPt, predTime, R, model, "kf", true);

                case "Ti_ekf"
                    [data.impErrTiEKF_array{end + 1}, data.shoErrTiEKF_array{end + 1}] = ...
                        runMultiSim(P_array_ekf, T_array, T_ref, y, impPt, shoPt, predTime, R, model, "ekf", true);

                case "Rad_kf"
                    [data.impErrRad_array{end + 1}, data.shoErrRad_array{end + 1}] = ...
                        runMultiSim(P_array_kf, T_array, T_ref, y_Rad, impPt, shoPt, predTime, R, model, "kf", false);

                case "Rad_ekf"
                    [data.impErrRadEKF_array{end + 1}, data.shoErrRadEKF_array{end + 1}] = ...
                        runMultiSim(P_array_ekf, T_array, T_ref, y_Rad, impPt, shoPt, predTime, R, model, "ekf", false);
            end

            data.('numOfSim' + keys(j)) = numOfPrevSim(j) + i;
        end

        % salva erros em weights
        save(resultsFilePath, '-struct', 'data');
        toc;
        fprintf("\n");
    end
end
