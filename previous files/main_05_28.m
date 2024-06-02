clear; close all; clc;



% variaveis de mruv
g = 9.81;
vel0 = [180 180 250]; a = [0 0 -g];
% chao
pos_floor = 0;
% intervalo de tempo de previsao
predTime = [1 5];

% entrada do sistema
u = a';

% ruido de medicao
sigma2_v = [1e2; 1e2; 1e2];
Ts_5ms = 0.005;
v = randn(3, 60/Ts_5ms);
v = v - mean(v, 2)*ones(1, 60/Ts_5ms);
v = v.*sqrt(sigma2_v)./std(v, 0, 2);

% matriz covariancia do ruido de processo
Q = cov(u');

% matriz covariancia do ruido de medicao
R = [sigma2_v(1) 0 0;     % roubo
     0 sigma2_v(2) 0;
     0 0 sigma2_v(3)];
R = 0.1*R;                % R chutado para 10 vezes menos
% R = 10*R;               % e 10 vezes mais
% quero testar R obtido por autocovariance least-square



% variaveis para simulacao
% periodo de amostragem
Ts_5ms = 0.005;
Ts_10ms = 0.01;
Ts_15ms = 0.015;
Ts_20ms = 0.02;
Ts_25ms = 0.025;
Ts_30ms = 0.03;
Ts_35ms = 0.035;
% metodo
method1 = 1;
method2 = 2;
method3 = 3;
% erro de covariancia
P1 = 1e1*eye(6);
P2 = 1e2*eye(6);
P3 = 1e3*eye(6);
P4 = 1e4*eye(6);
P5 = 1e5*eye(6);
P6 = 1e6*eye(6);
P7 = 1e7*eye(6);
% ruido de medicao
v_5ms = v;
v_10ms = v(:, 2:2:end);
v_15ms = v(:, 3:3:end);
v_20ms = v(:, 4:4:end);
v_25ms = v(:, 5:5:end);
v_30ms = v(:, 6:6:end);
v_35ms = v(:, 7:7:end);


% simulacao
[~, ~, ~, ~, ~, shootingError_P2_m1_10ms] = runSim(Ts_10ms, method1, P2, v_10ms, g, vel0, u, pos_floor, predTime, Q, R);
[~, ~, ~, ~, ~, shootingError_P3_m1_10ms] = runSim(Ts_10ms, method1, P3, v_10ms, g, vel0, u, pos_floor, predTime, Q, R);
[~, ~, ~, ~, ~, shootingError_P4_m1_10ms] = runSim(Ts_10ms, method1, P4, v_10ms, g, vel0, u, pos_floor, predTime, Q, R);
[~, ~, ~, ~, ~, shootingError_P5_m1_10ms] = runSim(Ts_10ms, method1, P5, v_10ms, g, vel0, u, pos_floor, predTime, Q, R);
[~, ~, ~, ~, ~, shootingError_P6_m1_10ms] = runSim(Ts_10ms, method1, P6, v_10ms, g, vel0, u, pos_floor, predTime, Q, R);
[~, ~, ~, ~, ~, shootingError_P7_m1_10ms] = runSim(Ts_10ms, method1, P7, v_10ms, g, vel0, u, pos_floor, predTime, Q, R);

[~, ~, ~, ~, ~, shootingError_P2_m2_10ms] = runSim(Ts_10ms, method2, P2, v_10ms, g, vel0, u, pos_floor, predTime, Q, R);
[~, ~, ~, ~, ~, shootingError_P3_m2_10ms] = runSim(Ts_10ms, method2, P3, v_10ms, g, vel0, u, pos_floor, predTime, Q, R);
[~, ~, ~, ~, ~, shootingError_P4_m2_10ms] = runSim(Ts_10ms, method2, P4, v_10ms, g, vel0, u, pos_floor, predTime, Q, R);
[~, ~, ~, ~, ~, shootingError_P5_m2_10ms] = runSim(Ts_10ms, method2, P5, v_10ms, g, vel0, u, pos_floor, predTime, Q, R);
[~, ~, ~, ~, ~, shootingError_P6_m2_10ms] = runSim(Ts_10ms, method2, P6, v_10ms, g, vel0, u, pos_floor, predTime, Q, R);
[~, ~, ~, ~, ~, shootingError_P7_m2_10ms] = runSim(Ts_10ms, method2, P7, v_10ms, g, vel0, u, pos_floor, predTime, Q, R);

[~, ~, ~, ~, ~, shootingError_P2_m3_10ms] = runSim(Ts_10ms, method3, P2, v_10ms, g, vel0, u, pos_floor, predTime, Q, R);
[~, ~, ~, ~, ~, shootingError_P3_m3_10ms] = runSim(Ts_10ms, method3, P3, v_10ms, g, vel0, u, pos_floor, predTime, Q, R);
[~, ~, ~, ~, ~, shootingError_P4_m3_10ms] = runSim(Ts_10ms, method3, P4, v_10ms, g, vel0, u, pos_floor, predTime, Q, R);
[~, ~, ~, ~, ~, shootingError_P5_m3_10ms] = runSim(Ts_10ms, method3, P5, v_10ms, g, vel0, u, pos_floor, predTime, Q, R);
[~, ~, ~, ~, ~, shootingError_P6_m3_10ms] = runSim(Ts_10ms, method3, P6, v_10ms, g, vel0, u, pos_floor, predTime, Q, R);
[~, ~, ~, ~, ~, shootingError_P7_m3_10ms] = runSim(Ts_10ms, method3, P7, v_10ms, g, vel0, u, pos_floor, predTime, Q, R);

[~, ~, ~, t_5ms, ~, shootingError_P1_m1_5ms] = runSim(Ts_5ms, method1, P7, v_5ms, g, vel0, u, pos_floor, predTime, Q, R);
[~, ~, ~, ~, ~, shootingError_P1_m2_5ms] = runSim(Ts_5ms, method2, P7, v_5ms, g, vel0, u, pos_floor, predTime, Q, R);
[~, ~, ~, ~, ~, shootingError_P1_m3_5ms] = runSim(Ts_5ms, method3, P7, v_5ms, g, vel0, u, pos_floor, predTime, Q, R);
[~, ~, ~, t_10ms, ~, shootingError_P1_m1_10ms] = runSim(Ts_10ms, method1, P7, v_10ms, g, vel0, u, pos_floor, predTime, Q, R);
[~, ~, ~, ~, ~, shootingError_P1_m2_10ms] = runSim(Ts_10ms, method2, P7, v_10ms, g, vel0, u, pos_floor, predTime, Q, R);
[~, ~, ~, ~, ~, shootingError_P1_m3_10ms] = runSim(Ts_10ms, method3, P7, v_10ms, g, vel0, u, pos_floor, predTime, Q, R);
[~, ~, ~, t_15ms, ~, shootingError_P1_m1_15ms] = runSim(Ts_15ms, method1, P7, v_15ms, g, vel0, u, pos_floor, predTime, Q, R);
[~, ~, ~, ~, ~, shootingError_P1_m2_15ms] = runSim(Ts_15ms, method2, P7, v_15ms, g, vel0, u, pos_floor, predTime, Q, R);
[~, ~, ~, ~, ~, shootingError_P1_m3_15ms] = runSim(Ts_15ms, method3, P7, v_15ms, g, vel0, u, pos_floor, predTime, Q, R);
[~, ~, ~, t_20ms, ~, shootingError_P1_m1_20ms] = runSim(Ts_20ms, method1, P7, v_20ms, g, vel0, u, pos_floor, predTime, Q, R);
[~, ~, ~, ~, ~, shootingError_P1_m2_20ms] = runSim(Ts_20ms, method2, P7, v_20ms, g, vel0, u, pos_floor, predTime, Q, R);
[~, ~, ~, ~, ~, shootingError_P1_m3_20ms] = runSim(Ts_20ms, method3, P7, v_20ms, g, vel0, u, pos_floor, predTime, Q, R);
[~, ~, ~, t_25ms, ~, shootingError_P1_m1_25ms] = runSim(Ts_25ms, method1, P7, v_25ms, g, vel0, u, pos_floor, predTime, Q, R);
[~, ~, ~, ~, ~, shootingError_P1_m2_25ms] = runSim(Ts_25ms, method2, P7, v_25ms, g, vel0, u, pos_floor, predTime, Q, R);
[~, ~, ~, ~, ~, shootingError_P1_m3_25ms] = runSim(Ts_25ms, method3, P7, v_25ms, g, vel0, u, pos_floor, predTime, Q, R);
[~, ~, ~, t_30ms, ~, shootingError_P1_m1_30ms] = runSim(Ts_30ms, method1, P7, v_30ms, g, vel0, u, pos_floor, predTime, Q, R);
[~, ~, ~, ~, ~, shootingError_P1_m2_30ms] = runSim(Ts_30ms, method2, P7, v_30ms, g, vel0, u, pos_floor, predTime, Q, R);
[~, ~, ~, ~, ~, shootingError_P1_m3_30ms] = runSim(Ts_30ms, method3, P7, v_30ms, g, vel0, u, pos_floor, predTime, Q, R);
[~, ~, ~, t_35ms, ~, shootingError_P1_m1_35ms] = runSim(Ts_35ms, method1, P7, v_35ms, g, vel0, u, pos_floor, predTime, Q, R);
[~, ~, ~, ~, ~, shootingError_P1_m2_35ms] = runSim(Ts_35ms, method2, P7, v_35ms, g, vel0, u, pos_floor, predTime, Q, R);
[~, ~, ~, ~, ~, shootingError_P1_m3_35ms] = runSim(Ts_35ms, method3, P7, v_35ms, g, vel0, u, pos_floor, predTime, Q, R);


% plots
% fullPlot(y, y_true, y_est);
% dynamicPlotRecord(y_true, y_est, impactPointPredArray, shootingPointPredArray);

% mudando escala das diversas shootingError arrays
arrayLength_5ms = length(shootingError_P1_m1_5ms);
newLength = length(shootingError_P1_m1_10ms);
new_shootingError_P1_m1_10ms = interp1(linspace(1, newLength, newLength), shootingError_P1_m1_10ms, ...
                                       linspace(1, newLength, arrayLength_5ms), 'linear');
new_shootingError_P1_m2_10ms = interp1(linspace(1, newLength, newLength), shootingError_P1_m2_10ms, ...
                                       linspace(1, newLength, arrayLength_5ms), 'linear');
new_shootingError_P1_m3_10ms = interp1(linspace(1, newLength, newLength), shootingError_P1_m3_10ms, ...
                                       linspace(1, newLength, arrayLength_5ms), 'linear');
newLength = length(shootingError_P1_m1_15ms);
new_shootingError_P1_m1_15ms = interp1(linspace(1, newLength, newLength), shootingError_P1_m1_15ms, ...
                                       linspace(1, newLength, arrayLength_5ms), 'linear');
new_shootingError_P1_m2_15ms = interp1(linspace(1, newLength, newLength), shootingError_P1_m2_15ms, ...
                                       linspace(1, newLength, arrayLength_5ms), 'linear');
new_shootingError_P1_m3_15ms = interp1(linspace(1, newLength, newLength), shootingError_P1_m3_15ms, ...
                                       linspace(1, newLength, arrayLength_5ms), 'linear');
newLength = length(shootingError_P1_m1_20ms);
new_shootingError_P1_m1_20ms = interp1(linspace(1, newLength, newLength), shootingError_P1_m1_20ms, ...
                                       linspace(1, newLength, arrayLength_5ms), 'linear');
new_shootingError_P1_m2_20ms = interp1(linspace(1, newLength, newLength), shootingError_P1_m2_20ms, ...
                                       linspace(1, newLength, arrayLength_5ms), 'linear');
new_shootingError_P1_m3_20ms = interp1(linspace(1, newLength, newLength), shootingError_P1_m3_20ms, ...
                                       linspace(1, newLength, arrayLength_5ms), 'linear');
newLength = length(shootingError_P1_m1_25ms);
new_shootingError_P1_m1_25ms = interp1(linspace(1, newLength, newLength), shootingError_P1_m1_25ms, ...
                                       linspace(1, newLength, arrayLength_5ms), 'linear');
new_shootingError_P1_m2_25ms = interp1(linspace(1, newLength, newLength), shootingError_P1_m2_25ms, ...
                                       linspace(1, newLength, arrayLength_5ms), 'linear');
new_shootingError_P1_m3_25ms = interp1(linspace(1, newLength, newLength), shootingError_P1_m3_25ms, ...
                                       linspace(1, newLength, arrayLength_5ms), 'linear');
newLength = length(shootingError_P1_m1_30ms);
new_shootingError_P1_m1_30ms = interp1(linspace(1, newLength, newLength), shootingError_P1_m1_30ms, ...
                                       linspace(1, newLength, arrayLength_5ms), 'linear');
new_shootingError_P1_m2_30ms = interp1(linspace(1, newLength, newLength), shootingError_P1_m2_30ms, ...
                                       linspace(1, newLength, arrayLength_5ms), 'linear');
new_shootingError_P1_m3_30ms = interp1(linspace(1, newLength, newLength), shootingError_P1_m3_30ms, ...
                                       linspace(1, newLength, arrayLength_5ms), 'linear');
newLength = length(shootingError_P1_m1_35ms);
new_shootingError_P1_m1_35ms = interp1(linspace(1, newLength, newLength), shootingError_P1_m1_35ms, ...
                                       linspace(1, newLength, arrayLength_5ms), 'linear');
new_shootingError_P1_m2_35ms = interp1(linspace(1, newLength, newLength), shootingError_P1_m2_35ms, ...
                                       linspace(1, newLength, arrayLength_5ms), 'linear');
new_shootingError_P1_m3_35ms = interp1(linspace(1, newLength, newLength), shootingError_P1_m3_35ms, ...
                                       linspace(1, newLength, arrayLength_5ms), 'linear');

% comparando P
errorArray = [];
errorArray(1, :, :) = [shootingError_P1_m1_10ms; shootingError_P2_m1_10ms; shootingError_P3_m1_10ms; shootingError_P4_m1_10ms; ...
                       shootingError_P5_m1_10ms; shootingError_P6_m1_10ms; shootingError_P7_m1_10ms];
errorArray(2, :, :) = [shootingError_P1_m2_10ms; shootingError_P2_m2_10ms; shootingError_P3_m2_10ms; shootingError_P4_m2_10ms; ...
                       shootingError_P5_m2_10ms; shootingError_P6_m2_10ms; shootingError_P7_m2_10ms];
errorArray(3, :, :) = [shootingError_P1_m3_10ms; shootingError_P2_m3_10ms; shootingError_P3_m3_10ms; shootingError_P4_m3_10ms; ...
                       shootingError_P5_m3_10ms; shootingError_P6_m3_10ms; shootingError_P7_m3_10ms];
plotTitle = ["Erro de Disparo - Metodo 1 e Ts = 10ms";
             "Erro de Disparo - Metodo 2 e Ts = 10ms";
             "Erro de Disparo - Metodo 3 e Ts = 10ms"];
plotLegend = ["P = 1e1" "P = 1e2" "P = 1e3" "P = 1e4" "P = 1e5" "P = 1e6" "P = 1e7"];
errorPlot(t_10ms, errorArray, plotTitle, plotLegend);

% comparando Ts
errorArray = [];
errorArray(1, :, :) = [shootingError_P1_m1_5ms; new_shootingError_P1_m1_10ms; new_shootingError_P1_m1_15ms; ...
                       new_shootingError_P1_m1_20ms; new_shootingError_P1_m1_25ms; new_shootingError_P1_m1_30ms; new_shootingError_P1_m1_35ms];
errorArray(2, :, :) = [shootingError_P1_m2_5ms; new_shootingError_P1_m2_10ms; new_shootingError_P1_m2_15ms; ...
                       new_shootingError_P1_m2_20ms; new_shootingError_P1_m2_25ms; new_shootingError_P1_m2_30ms; new_shootingError_P1_m2_35ms];
errorArray(3, :, :) = [shootingError_P1_m3_5ms; new_shootingError_P1_m3_10ms; new_shootingError_P1_m3_15ms; ...
                       new_shootingError_P1_m3_20ms; new_shootingError_P1_m3_25ms; new_shootingError_P1_m3_30ms; new_shootingError_P1_m3_35ms];
plotTitle = ["Erro de Disparo - Metodo 1";
             "Erro de Disparo - Metodo 2";
             "Erro de Disparo - Metodo 3"];
plotLegend = ["Ts = 5ms" "Ts = 10ms" "Ts = 15ms" "Ts = 20ms" "Ts = 25ms" "Ts = 30ms" "Ts = 35ms"];
errorPlot(t_5ms, errorArray, plotTitle, plotLegend);

% comparando cada metodo
errorArray = [];
errorArray(1, :, :) = [shootingError_P1_m1_5ms; shootingError_P1_m2_5ms; shootingError_P1_m3_5ms];
errorArray(2, :, :) = [new_shootingError_P1_m1_10ms; new_shootingError_P1_m2_10ms; new_shootingError_P1_m3_10ms];
errorArray(3, :, :) = [new_shootingError_P1_m1_15ms; new_shootingError_P1_m2_15ms; new_shootingError_P1_m3_15ms];
plotTitle = ["Erro de Disparo - Ts = 5ms";
             "Erro de Disparo - Ts = 10ms";
             "Erro de Disparo - Ts = 15ms"];
plotLegend = ["Metodo 1" "Metodo 2" "Metodo 3"];
errorPlot(t_5ms, errorArray, plotTitle, plotLegend);
errorArray = [];
errorArray(1, :, :) = [new_shootingError_P1_m1_20ms; new_shootingError_P1_m2_20ms; new_shootingError_P1_m3_20ms];
errorArray(2, :, :) = [new_shootingError_P1_m1_25ms; new_shootingError_P1_m2_25ms; new_shootingError_P1_m3_25ms];
errorArray(3, :, :) = [new_shootingError_P1_m1_30ms; new_shootingError_P1_m2_30ms; new_shootingError_P1_m3_30ms];
errorArray(4, :, :) = [new_shootingError_P1_m1_35ms; new_shootingError_P1_m2_35ms; new_shootingError_P1_m3_35ms];
plotTitle = ["Erro de Disparo - Ts = 20ms";
             "Erro de Disparo - Ts = 25ms";
             "Erro de Disparo - Ts = 30ms";
             "Erro de Disparo - Ts = 35ms"];
plotLegend = ["Metodo 1" "Metodo 2" "Metodo 3"];
errorPlot(t_5ms, errorArray, plotTitle, plotLegend);



% simulacao
function [y, y_est, y_true, t, impactError, shootingError] = runSim(Ts, method, P, v, g, vel0, u, pos_floor, predTime, Q, R)
    arrayLength = ceil(predTime(2)/Ts);
    t = (1:arrayLength)*Ts;
    
    % amostra em que comeca a medida
    % referente ao tempo de 10 segundos
    detectionThreshold = ceil(10/Ts);
    
    % matrizes de espaco de estados
    A = eye(6);
    A(1:3, 4:end) = Ts*eye(3);
    B = [Ts^2/2*eye(3); Ts*eye(3)];
    C = [eye(3) zeros(3)];
    A_rev = inv(A);
    B_rev = -A_rev*B;
    
    % planta
    x_ss0 = [zeros(1,3) vel0]';
    plant = Plant(A, B, C, x_ss0);
    
    % Kalman filter
    kf = KalmanFilter(A, B, C, Q, R, P);
    kf_rev = KalmanFilter(A_rev, B_rev, C, Q, R, P);
    
    % estado estimado
    x_est = zeros(length(A), arrayLength);
    
    y_size = [size(C, 1), arrayLength];
    % saida real
    y_true = [];
    % saida medida
    y = zeros(y_size);
    % saida estimada
    y_est = zeros(y_size);
    
    impactError = zeros(1, arrayLength);
    shootingError = zeros(1, arrayLength);
    
    impactPoint = zeros(3, 1);
    shootingPoint = zeros(3, 1);

    impactPointPredArray = zeros(y_size);
    shootingPointPredArray = zeros(y_size);
    
    
    % gerar saida real
    while 1
        % planta
        [plant, y_true_i] = plant.run(u);
        y_true = [y_true y_true_i];
        % parar no impacto
        if y_true_i(3) < pos_floor
            [impactPoint, shootingPoint] = setImpactShootingPoint(x_ss0, g, pos_floor);
            break;
        end
    end
    
    
    i = 0;
    % loop
    while 1
        i = i + 1;
        currentTime = i*Ts;
    
        % radar
        y(:, i) = y_true(:, i+detectionThreshold) + v(:, i+detectionThreshold);
    
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
              numberOfFiltering = 1;
              [kf, kf_rev, y_est, x_est] = runMultipleKalmanFiltering(numberOfFiltering, y_est, x_est, y, u, kf, kf_rev, i);
           case 3
              numberOfFiltering = 5;
              [kf, kf_rev, y_est, x_est] = runMultipleKalmanFiltering(numberOfFiltering, y_est, x_est, y, u, kf, kf_rev, i);
        end
        % previsao da trajetoria de impacto e disparo no intervalo de predTime
        if currentTime > predTime(1) && currentTime <= predTime(2)
            [impactPointPred, shootingPointPred] = setImpactShootingPoint(x_est(:, i), g, pos_floor);
            impactPointPredArray(:, i) = impactPointPred;
            shootingPointPredArray(:, i) = shootingPointPred;
            % calcular erro dos pontos de impacto e de disparo
            [impactError(i), shootingError(i)] = getImpactShootingError(impactPointPred, shootingPointPred, impactPoint, shootingPoint);

        % parar quando tempo de execucao estiver fora do intervalo de predTime
        elseif currentTime > predTime(2)
            y = y(:, 1:i);
            y_est = y_est(:, 1:i);
            break
        end


        % dynamic plot
        % dynamicPlot(y_true, y_est, impactPointPred, shootingPointPred, i);
    end
end



% rodar kalman filter do inicio ate o fim e do fim ate o inicio,
% n vezes, para cada iteracao
function [kf, kf_rev, y_est, x_est] = runMultipleKalmanFiltering(numberOfFiltering, y_est, x_est, y, u, kf, kf_rev, i)
    for n = 1:numberOfFiltering
        kf.x = x_est(:, 1);
        for j = 2:i
            [kf, y_est(:, j), x_est(:, j)] = kf.run(y(:, j), u);
        end
        kf_rev.x = x_est(:, i);
        for j = 2:i
            l = i - j + 1;
            [kf_rev, y_est(:, l), x_est(:, l)] = kf_rev.run(y(:, l), u);
        end
    end
end



% funcao nao usada mais
% previsao da trajetoria de impacto
function [y_pred] = setImpactTrajectory(plant_pred, y_pred, u, pos_floor)
    j = 0;
    while 1
        j = j + 1;
        [plant_pred, y_pred(:, j)] = plant_pred.run(u);
        if y_pred(3, j) < pos_floor
            y_pred = y_pred(:, 1:j);
            return
        end
    end
end



% funcao nao usada mais
% previsao da trajetoria de disparo
function [y_pred_rev] = setShootingTrajectory(plant_pred_rev, y_pred_rev, u, pos_floor)
    j = 0;
    while 1
        j = j + 1;
        [plant_pred_rev, y_pred_rev(:, j)] = plant_pred_rev.run(u);
        if y_pred_rev(3, j) < pos_floor
            y_pred_rev = y_pred_rev(:, 1:j);
            return
        end
    end
end



% previsao do ponto de impacto e disparo
function [impactPointPred, shootingPointPred] = setImpactShootingPoint(x_est, g, pos_floor)
    pos0x = x_est(1); pos0y = x_est(2); pos0z = x_est(3);
    vel0x = x_est(4); vel0y = x_est(5); vel0z = x_est(6);

    ti = (vel0z + sqrt(vel0z^2 + 2*g*(pos0z - pos_floor)))/g;
    td = (vel0z - sqrt(vel0z^2 + 2*g*(pos0z - pos_floor)))/g;

    pos1x = pos0x + vel0x*ti; pos2x = pos0x + vel0x*td;
    pos1y = pos0y + vel0y*ti; pos2y = pos0y + vel0y*td;

    impactPointPred = [pos1x; pos1y; pos_floor];
    shootingPointPred = [pos2x; pos2y; pos_floor];
end



% calcular erro dos pontos de impacto e de disparo
function [impactError, shootingError] = getImpactShootingError(impactPointPred, shootingPointPred, impactPoint, shootingPoint)
    impactError = norm(impactPointPred - impactPoint);
    shootingError = norm(shootingPointPred - shootingPoint);
end



function dynamicPlot(y_true, y_est, impactPointPred, shootingPointPred, i)
    if i == 1
        figure
    end
    % remove the zeros points
    zero_columns = all(y_est == 0, 1);
    y_est = y_est(:, ~zero_columns);

    plot3(y_true(1, :), y_true(2, :), y_true(3, :), 'b', ...
         y_est(1, :), y_est(2, :), y_est(3, :), 'rx', ...
         impactPointPred(1, :), impactPointPred(2, :), impactPointPred(3, :), 'r*',...
         shootingPointPred(1, :), shootingPointPred(2, :), shootingPointPred(3, :), 'r*', 'LineWidth', 2)
    axis([-1e3, 12e3, -1e3, 12e3, 0, 4e3])
    xlabel('x (m)'), ylabel('y (m)'), zlabel('z (m)')
    title('Posição do Projétil')
    legend('Medida', 'Filtrada', 'Ponto de Impacto', 'Ponto de Disparo')
    grid on
    set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1], 'Color', [1 1 1]);

    pause(0);
end



function dynamicPlotRecord(y_true, y_est, impactPointPredArray, shootingPointPredArray)
    video = VideoWriter('plot.mp4', 'MPEG-4');
    video.FrameRate = 60;
    open(video)

    zero_columns = all(impactPointPredArray == 0, 1);
    num_zero_columns = sum(zero_columns);

    figure

    plot(y_true(2, :), y_true(3, :), 'b', 'LineWidth', 2)
    hold on

    h = animatedline('LineStyle', '-', 'Color', 'r', 'LineWidth', 4);
    p1 = animatedline('Marker', 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'r', 'MarkerSize', 5);
    p2 = animatedline('Marker', 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'r', 'MarkerSize', 5);

    axis([-1e3, 12e3, 0, 4e3])
    %view(45, 15);
    xlabel('x (m)'), ylabel('z (m)'), zlabel('z (m)')
    title('Posição do Projétil')
    legend('Medida', 'Filtrada', 'Ponto de Impacto', 'Ponto de Disparo')
    grid on
    set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1], 'Color', [1 1 1]);

    for i = 1:num_zero_columns
        addpoints(h, y_est(2, i), y_est(3, i));
    end
    for i = num_zero_columns + 1:length(impactPointPredArray)
        clearpoints(p1);
        clearpoints(p2);
        addpoints(h, y_est(2, i), y_est(3, i));
        addpoints(p1, impactPointPredArray(2, i), impactPointPredArray(3, i));
        addpoints(p2, shootingPointPredArray(2, i), shootingPointPredArray(3, i));
        drawnow
        frame = getframe(gcf);
        writeVideo(video, frame)
    end

    close(video)
end



function fullPlot(y, y_true, y_est)
    figure

    plot3(y_true(1, :), y_true(2, :), y_true(3, :), 'b', ...
          y(1, :), y(2, :), y(3, :), 'k', ...
          y_est(1, :), y_est(2, :), y_est(3, :), 'rx', 'LineWidth', 2)
    xlabel('x (m)'), ylabel('y (m)'), zlabel('z (m)')
    title('Posição do Projétil')
    legend('Real', 'Medida', 'Filtrada', 'ItemHitFcn', @cb_legend)
    axis tight
    grid on
    set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1], 'Color', [1 1 1]);
end



function errorPlot(t, errorArray, plotTitle, plotLegend)
    figure
    for i = 1:size(errorArray, 1)
        subplot(1, size(errorArray, 1), i)
        hold on
        for j = 1:size(errorArray, 2)
            error = reshape(errorArray(i, j, :), 1, []);
            plot(t, error, 'LineWidth', 2)
        end
        title(plotTitle(i))
        ylim([0 80])
        xlabel('tempo (s)'), ylabel('erro (m)')
        legend(plotLegend, 'ItemHitFcn', @cb_legend)
        grid on
    end
    set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1], 'Color', [1 1 1]);
end



function cb_legend(~, evt)
    if strcmp(evt.Peer.Visible, 'on')
        evt.Peer.Visible = 'off';
    else 
        evt.Peer.Visible = 'on';
    end
end