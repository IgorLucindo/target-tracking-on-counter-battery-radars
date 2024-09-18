clear; clc; close all;

% escolher metodo
% metodo 1 - filtragem direta
% metodo 2 - filtragem direta e reversa 1 vez
method = 1;
% method = 2;
numOfsim = 1;

% variaveis de modelo
detecThresh = 0;
T = 0.05;
predTime = [1 5];
g = 9.81;
v0 = [180 180 250];
gama = 1e-4;
a = [0 0 -g];
p_floor = 0;
u = a';


% matriz covariancia do ruido de processo
Q = zeros(7);

% valor inicial de matriz covariancia de estado
% valores nao otimos
P = [1e8*eye(3) zeros(3) zeros(3, 1);
     zeros(3) 1e8*eye(3) zeros(3, 1);
     zeros(1, 3) zeros(1, 3) 0];

% matriz covariancia do ruido de medicao
sigma2_n = 1e0;
R = sigma2_n*eye(3);   % roubo
R = 0.1*R;             % R chutado para 10 vezes menor

% model params
[f_model, ~, ~, ~, f_model_rev, ~] = getParamsEkf(0.05);
% EKF params
[f, h, F, H, f_rev, F_rev] = getParamsEkf(T);

% Extended Kalman filter
ekf = ExtendedKalmanFilter(f, h, F, H, Q, R, P);
ekf_rev = ExtendedKalmanFilter(f_rev, h, F_rev, H, Q, R, P);


% trajetoria real
x_true = [1e3*ones(1, 3) v0 gama]';
x_aux = x_true;
i = 0;
while 1
    i = i + 1;
    x_aux = f(x_aux, u);
    if i == detecThresh + 1
        cheating_gama = x_aux(7);
    end
    y_true(:, i) = h(x_aux, u);
    if y_true(3, i) < p_floor
        break
    end
end
x_aux = x_true;
while 1
    x_aux = f_rev(x_aux, u);
    if x_aux(3) < p_floor
        break
    end
end


% ponto de impacto e disparo
impPt = y_true(:, end);
shoPt = x_aux(1:3);
impPtPred = zeros(3, 1);
shoPtPred = zeros(3, 1);


arrayLength = length(y_true);
% ruido de medicao
n = randn(3, arrayLength);
n = n - mean(n, 2)*ones(1, arrayLength);
n = n*sqrt(sigma2_n)./std(n, 0, 2);

% trajetoria medida
arrayLength = ceil(predTime(2)/T);
y = y_true(:, detecThresh + 1:detecThresh+arrayLength) + n(:, 1:arrayLength);

arrayLength = ceil((predTime(2))/T);
% trajetoria estimada
y_est = zeros(3, arrayLength);


arrayLength = ceil((predTime(2) - predTime(1))/T);
% erro de impacto e disparo
impErr = zeros(1, arrayLength);
shoErr = zeros(1, arrayLength);
impErrAvg = zeros(1, arrayLength);
shoErrAvg = zeros(1, arrayLength);


% simulacoes
for k = 1:numOfsim
    % loop
    x0 = [y(:, 1); 0; 0; 200; cheating_gama];
    [ekf, y_est(:, 1), x_est] = ekf.setInitialState(x0);
    % ekf.x = x0;
    i = 1; j = 0;
    while 1
        i = i + 1;
        currentTime = i*T;

        % parar quando tempo de execucao estiver fora do intervalo de predTime
        if currentTime > predTime(2)
            break
        end

        % extended kalman filter
        switch method
            case 1
                % metodo 1
                [ekf, y_est(:, i), x_est] = ekf.run(y(:, i), u);
            case 2
                % metodo 2
                numberOfFiltering = 1;
                [ekf, ekf_rev, y_est, x_est] = runMultiKf(numberOfFiltering, y, u, P, ekf, ekf_rev, i);

                % parar caso matriz P exploda
                if isnan(ekf.P(1))
                    "iteracao: " + i
                    ekf.P
                    ekf_rev.P
                    break
                end
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

    impErrAvg = impErrAvg + impErr;
    shoErrAvg = shoErrAvg + shoErr;
end


impErrAvg = impErrAvg/numOfsim;
shoErrAvg = shoErrAvg/numOfsim;

% plot das trajetorias
figure
plot(impErrAvg, 'LineWidth', 2)
title('impacto')
grid on
xlabel('amostras'), ylabel('erro (m)')
ylim([0 1e3])

figure
plot(shoErrAvg, 'LineWidth', 2)
title('disparo')
grid on
xlabel('amostras'), ylabel('erro (m)')
ylim([0 1e3])