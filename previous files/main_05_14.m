clear; close all; clc;



% variaveis de mruv
g = 9.81;
vel0 = [180 180 250]; a = [0 0 -g];
% chao
floor = 0;
% intervalo de tempo de previsao
predTime = [1 5];

% entrada do sistema
u = a';

% ruido de medicao
sigma2_v = [1e2; 1e2; 1e2];
Ts_10ms = 0.01;
v = randn(3, 60/Ts_10ms);
v = v - mean(v, 2)*ones(1, 60/Ts_10ms);
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
Ts_10ms = 0.01;
Ts_20ms = 0.02;
% metodo
method1 = 1;
method2 = 2;
method3 = 3;
% erro de covariancia
P1 = 1e3*eye(6);
P2 = 1e2*eye(6);
P3 = 1e1*eye(6);
% ruido de medicao
v_10ms = v;
v_20ms = v(2:2:end);

% simulacao
[y, y_est, y_true, t, impactError_P1_m1_20ms, shootingError_P1_m1_20ms] = runSim(Ts_20ms, method1, P1, v_20ms, g, vel0, u, floor, predTime, Q, R);
[y, y_est, y_true, t, impactError_P2_m1_20ms, shootingError_P2_m1_20ms] = runSim(Ts_20ms, method1, P2, v_20ms, g, vel0, u, floor, predTime, Q, R);
[y, y_est, y_true, t, impactError_P3_m1_20ms, shootingError_P3_m1_20ms] = runSim(Ts_20ms, method1, P3, v_20ms, g, vel0, u, floor, predTime, Q, R);
[y, y_est, y_true, t, impactError_P1_m2_20ms, shootingError_P1_m2_20ms] = runSim(Ts_20ms, method2, P1, v_20ms, g, vel0, u, floor, predTime, Q, R);
[y, y_est, y_true, t, impactError_P2_m2_20ms, shootingError_P2_m2_20ms] = runSim(Ts_20ms, method2, P2, v_20ms, g, vel0, u, floor, predTime, Q, R);
[y, y_est, y_true, t, impactError_P3_m2_20ms, shootingError_P3_m2_20ms] = runSim(Ts_20ms, method2, P3, v_20ms, g, vel0, u, floor, predTime, Q, R);
[y, y_est, y_true, t, impactError_P1_m3_20ms, shootingError_P1_m3_20ms] = runSim(Ts_20ms, method3, P1, v_20ms, g, vel0, u, floor, predTime, Q, R);
[y, y_est, y_true, t, impactError_P2_m3_20ms, shootingError_P2_m3_20ms] = runSim(Ts_20ms, method3, P2, v_20ms, g, vel0, u, floor, predTime, Q, R);
[y, y_est, y_true, t, impactError_P3_m3_20ms, shootingError_P3_m3_20ms] = runSim(Ts_20ms, method3, P3, v_20ms, g, vel0, u, floor, predTime, Q, R);
t_20ms = t;

[y, y_est, y_true, t, impactError_P1_m1_10ms, shootingError_P1_m1_10ms] = runSim(Ts_10ms, method1, P1, v_10ms, g, vel0, u, floor, predTime, Q, R);
[y, y_est, y_true, t, impactError_P2_m1_10ms, shootingError_P2_m1_10ms] = runSim(Ts_10ms, method1, P2, v_10ms, g, vel0, u, floor, predTime, Q, R);
[y, y_est, y_true, t, impactError_P3_m1_10ms, shootingError_P3_m1_10ms] = runSim(Ts_10ms, method1, P3, v_10ms, g, vel0, u, floor, predTime, Q, R);
[y, y_est, y_true, t, impactError_P1_m2_10ms, shootingError_P1_m2_10ms] = runSim(Ts_10ms, method2, P1, v_10ms, g, vel0, u, floor, predTime, Q, R);
[y, y_est, y_true, t, impactError_P2_m2_10ms, shootingError_P2_m2_10ms] = runSim(Ts_10ms, method2, P2, v_10ms, g, vel0, u, floor, predTime, Q, R);
[y, y_est, y_true, t, impactError_P3_m2_10ms, shootingError_P3_m2_10ms] = runSim(Ts_10ms, method2, P3, v_10ms, g, vel0, u, floor, predTime, Q, R);
[y, y_est, y_true, t, impactError_P1_m3_10ms, shootingError_P1_m3_10ms] = runSim(Ts_10ms, method3, P1, v_10ms, g, vel0, u, floor, predTime, Q, R);
[y, y_est, y_true, t, impactError_P2_m3_10ms, shootingError_P2_m3_10ms] = runSim(Ts_10ms, method3, P2, v_10ms, g, vel0, u, floor, predTime, Q, R);
[y, y_est, y_true, t, impactError_P3_m3_10ms, shootingError_P3_m3_10ms] = runSim(Ts_10ms, method3, P3, v_10ms, g, vel0, u, floor, predTime, Q, R);
t_10ms = t;

% plots
fullPlot(y, y_true, y_est);

% comparando P
errorArray = [];
errorArray(1, :, :) = [shootingError_P1_m1_20ms; shootingError_P2_m1_20ms; shootingError_P3_m1_20ms];
errorArray(2, :, :) = [shootingError_P1_m2_20ms; shootingError_P2_m2_20ms; shootingError_P3_m2_20ms];
errorArray(3, :, :) = [shootingError_P1_m3_20ms; shootingError_P2_m3_20ms; shootingError_P3_m3_20ms];
plotTitle = ["Erro de Disparo - Metodo 1 e Ts = 20ms";
             "Erro de Disparo - Metodo 2 e Ts = 20ms";
             "Erro de Disparo - Metodo 3 e Ts = 20ms"];
plotLegend = ["P = 1000" "P = 100" "P = 10"];
errorPlot(t_20ms, errorArray, plotTitle, plotLegend);

% comparando Ts
errorArray = [];
errorArray(1, :, :) = [repelem(shootingError_P1_m1_20ms, Ts_20ms/Ts_10ms); shootingError_P1_m1_10ms];
errorArray(2, :, :) = [repelem(shootingError_P3_m2_20ms, Ts_20ms/Ts_10ms); shootingError_P3_m2_10ms];
errorArray(3, :, :) = [repelem(shootingError_P3_m3_20ms, Ts_20ms/Ts_10ms); shootingError_P3_m3_10ms];
plotTitle = ["Erro de Disparo - Metodo 1 e Melhor P Associado";
             "Erro de Disparo - Metodo 2 e Melhor P Associado";
             "Erro de Disparo - Metodo 3 e Melhor P Associado"];
plotLegend = ["Ts = 20ms" "Ts = 10ms"];
errorPlot(t_10ms, errorArray, plotTitle, plotLegend);

% comparando cada metodo
errorArray = [];
errorArray(1, :, :) = [repelem(shootingError_P1_m1_20ms, Ts_20ms/Ts_10ms); repelem(shootingError_P3_m2_20ms, Ts_20ms/Ts_10ms); repelem(shootingError_P3_m3_20ms, Ts_20ms/Ts_10ms)];
errorArray(2, :, :) = [shootingError_P1_m1_10ms; shootingError_P3_m2_10ms; shootingError_P3_m3_10ms];
plotTitle = ["Erro de Disparo - Ts = 20ms";
             "Erro de Disparo - Ts = 10ms"];
plotLegend = ["Metodo 1" "Metodo 2" "Metodo 3"];
errorPlot(t, errorArray, plotTitle, plotLegend);



% simulacao
function [y, y_est, y_true, t, impactError, shootingError] = runSim(Ts, method, P, v, g, vel0, u, floor, predTime, Q, R)
    arrayLength = predTime(2)/Ts;
    t = (1:arrayLength)*Ts;
    
    % amostra em que comeca a medida
    % referente ao tempo de 10 segundos
    detectionThreshold = 10/Ts;
    
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
    
    impactPoint = zeros(2, 1);
    shootingPoint = zeros(2, 1);
    
    
    % gerar saida real
    while 1
        % planta
        [plant, y_true_i] = plant.run(u);
        y_true = [y_true y_true_i];
        % parar no impacto
        if y_true_i(3) < floor
            [impactPoint, shootingPoint] = setImpactShootingPoint(x_ss0, g, floor);
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
              [y_est, x_est] = runMultipleKalmanFiltering(numberOfFiltering, y_est, x_est, y, u, kf, kf_rev, i);
           case 3
              numberOfFiltering = 5;
              [y_est, x_est] = runMultipleKalmanFiltering(numberOfFiltering, y_est, x_est, y, u, kf, kf_rev, i);
        end
        % previsao da trajetoria de impacto e disparo no intervalo de predTime
        if currentTime > predTime(1) && currentTime <= predTime(2)
            [impactPointPred, shootingPointPred] = setImpactShootingPoint(x_est(:, i), g, floor);
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
function [y_est, x_est] = runMultipleKalmanFiltering(numberOfFiltering, y_est, x_est, y, u, kf, kf_rev, i)
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
function [y_pred] = setImpactTrajectory(plant_pred, y_pred, u, floor)
    j = 0;
    while 1
        j = j + 1;
        [plant_pred, y_pred(:, j)] = plant_pred.run(u);
        if y_pred(3, j) < floor
            y_pred = y_pred(:, 1:j);
            return
        end
    end
end



% funcao nao usada mais
% previsao da trajetoria de disparo
function [y_pred_rev] = setShootingTrajectory(plant_pred_rev, y_pred_rev, u, floor)
    j = 0;
    while 1
        j = j + 1;
        [plant_pred_rev, y_pred_rev(:, j)] = plant_pred_rev.run(u);
        if y_pred_rev(3, j) < floor
            y_pred_rev = y_pred_rev(:, 1:j);
            return
        end
    end
end



% previsao do ponto de impacto e disparo
function [impactPointPred, shootingPointPred] = setImpactShootingPoint(x_est, g, floor)
    pos0x = x_est(1); pos0y = x_est(2); pos0z = x_est(3);
    vel0x = x_est(4); vel0y = x_est(5); vel0z = x_est(6);

    ti = (vel0z + sqrt(vel0z^2 + 2*g*(pos0z - floor)))/g;
    td = (vel0z - sqrt(vel0z^2 + 2*g*(pos0z - floor)))/g;

    pos1x = pos0x + vel0x*ti; pos2x = pos0x + vel0x*td;
    pos1y = pos0y + vel0y*ti; pos2y = pos0y + vel0y*td;

    impactPointPred = [pos1x; pos1y; floor];
    shootingPointPred = [pos2x; pos2y; floor];
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

    plot3(y_true(1, :), y_true(2, :), y_true(3, :), 'b', ...
         y_est(1, :), y_est(2, :), y_est(3, :), 'rx', ...
         impactPointPred(1, :), impactPointPred(2, :), impactPointPred(3, :), 'r*',...
         shootingPointPred(1, :), shootingPointPred(2, :), shootingPointPred(3, :), 'g*', 'LineWidth', 1)
    axis([-1e3, 12e3, -1e3, 12e3, 0, 4e3])
    legend('Medida', 'Filtrada', 'Ponto de Impacto', 'Ponto de Disparo')
    grid on

    pause(0);
    fig=gcf;
    fig.Units='normalized';
    fig.OuterPosition=[0 0 1 1];
end



function fullPlot(y, y_true, y_est)
    figure

    plot3(y_true(1, :), y_true(2, :), y_true(3, :), 'b', ...
          y(1, :), y(2, :), y(3, :), 'k', ...
          y_est(1, :), y_est(2, :), y_est(3, :), 'r--', 'LineWidth', 2)
    xlabel('x (m)'), ylabel('y (m)'), zlabel('z (m)')
    title('Posicao do projetil')
    legend('Real', 'Medida', 'Filtrada')
    axis tight
    grid on

    fig=gcf;
    fig.Units='normalized';
    fig.OuterPosition=[0 0 1 1];
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
        ylim([0 40])
        xlabel('tempo (s)'), ylabel('erro (m)')
        legend(plotLegend)
        grid on
    end
    fig=gcf;
    fig.Units='normalized';
    fig.OuterPosition=[0 0 1 1];
end


