function [impErr, shoErr, state_error] = loop(kf, kf_rev, x0, T, predTime, f_model, f_model_rev, h, method, y, impPt, shoPt, u, p_floor, P, x_true_start, x_true_last)
    % ponto de impacto e disparo
    impPtPred = zeros(3, 1);
    shoPtPred = zeros(3, 1);

    arrayLength = ceil((predTime(2))/T);
    % trajetoria estimada
    y_est = zeros(3, arrayLength);

    arrayLength = ceil((predTime(2) - predTime(1))/T);
    % erro de impacto e disparo
    impErr = zeros(1, arrayLength);
    shoErr = zeros(1, arrayLength);

    method_ = 1;
    [kf, y_est(:, 1), x_est] = kf.setInitialState(x0);
    i = 1; j = 0;
    while 1
        i = i + 1;
        currentTime = i*T;

        % troca de metodo para n explodir
        time = 3;
        if i == floor(time/T) && method ~= 1
            method_ = method;
            kf_rev.P = P;
            kf_rev.P(7, 7) = 1e-11;
            [kf_rev, ~, x_est] = kf_rev.setInitialState(x_est);
            for k = 2:i
                l = i - k + 1;
                [kf_rev, y_est(:, l), x_est] = kf_rev.run(y(:, l), u);
            end
            [kf, ~, x_est] = kf.setInitialState(x_est);
        end

        % parar quando tempo de execucao estiver fora do intervalo de predTime
        if currentTime > predTime(2)
            switch method
                case 1
                    x_true = x_true_last;
                case {2, 3}
                    x_true = x_true_start;
            end
           state_error = abs((x_est - x_true)./x_true*100)';
            return
        end

        % extended kalman filter
        switch method_
            case 1
                % metodo 1
                [kf, y_est(:, i), x_est] = kf.run(y(:, i), u);
            case 2
                % metodo 2
                numOfFiltering = 1;
                [kf, kf_rev, y_est, x_est] = runMultiKf(numOfFiltering, y, u, kf, kf_rev, P, i);
            case 3
                % metodo 3
                numOfFiltering = 5;
                [kf, kf_rev, y_est, x_est] = runMultiKf(numOfFiltering, y, u, kf, kf_rev, P, i);
        end

        % previsao da trasdajetoria de impacto e disparo no intervalo de predTime
        if currentTime > predTime(1) && currentTime <= predTime(2)
            if j == arrayLength
                continue
            end
            j = j + 1;

            % prever ponto de impacto
            x_aux = x_est;
            while 1
                x_aux = f_model(x_aux, u);
                impPtPred = h(x_aux, u);
                if impPtPred(3) < p_floor
                    break
                end
            end
            
            % prever ponto de disparo
            x_aux = x_est;
            while 1
                x_aux = f_model_rev(x_aux, u);
                shoPtPred = h(x_aux, u);
                if shoPtPred(3) < p_floor
                    break
                end
            end

            % calcular erro dos pontos de impacto e de disparo
            impErr(j) = norm(impPtPred - impPt);
            shoErr(j) = norm(shoPtPred - shoPt);
        end
    end
end