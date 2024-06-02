clear; close all; clc;



% ====================== periodo de amostragem 100ms ========================



% variaveis de mruv
vel0 = [250 250]; a = [0 -9.81];

% período de amostragem
Ts = 0.1;
k_100ms = (0:1:60/Ts)';
k = k_100ms;

% matrizes de espaco de estados
A = [1 0 Ts 0;
     0 1 0  Ts;
     0 0 1  0;
     0 0 0  1];

B = [Ts^2/2 0;
     0      Ts^2/2;
     Ts     0;
     0      Ts];

C = [eye(2) zeros(2)];

A_rev = inv(A);

B_rev = -A_rev*B;


% entrada do sistema
u = a'*ones(1, length(k));

% ruido de medicao
sigma2_v = [1e2; 2e4];
v = randn(size(C, 1), 60/0.001);
v = v - mean(v, 2)*ones(1, 60/0.001);
v = v.*sqrt(sigma2_v)./std(v, 0, 2);
v_100ms = v(100:100:end);


% matriz covariancia do ruido de processo
Q = cov(u');
% matriz covariancia do ruido de medicao
R = [sigma2_v(1) 0;       % roubo
     0 sigma2_v(2)];
R = 0.1*R;                % R chutado para 10 vezes menos
% R = 10*R;               % e 10 vezes mais
% quero testar R obtido por autocovariance least-square

% erro de covariancia
% valor de P interfere muito para cada metodo
P1 = 1e3*eye(length(A));
P2 = 1e1*eye(length(A));
P3 = 1e1*eye(length(A));

% planta
x_ss0 = [0; 0; vel0(1); vel0(2)];
plant = Plant(A, B, C, x_ss0);
plant_pred = Plant(A, B, C, zeros(length(A), 1));
plant_pred_rev = Plant(A_rev, B_rev, C, zeros(length(A), 1));

% Kalman filter
kf1 = KalmanFilter(A, B, C, Q, R, P1);
kf2 = KalmanFilter(A, B, C, Q, R, P2);
kf3 = KalmanFilter(A, B, C, Q, R, P3);
kf_rev2 = KalmanFilter(A_rev, B_rev, C, Q, R, P2);
kf_rev3 = KalmanFilter(A_rev, B_rev, C, Q, R, P3);

x_size = [length(A), length(k)];
% estado estimado
x_est1 = zeros(x_size);
x_est2 = zeros(x_size);
x_est3 = zeros(x_size);

y_size = [size(C, 1), length(k)];
% saida real
y_true = zeros(y_size);
% saida medida
y = zeros(y_size);
% saida estimada
y_est1 = zeros(y_size);
y_est2 = zeros(y_size);
y_est3 = zeros(y_size);
% saida prevista
y_pred1 = zeros(y_size);
y_pred2 = zeros(y_size);
y_pred3 = zeros(y_size);
y_pred_rev1 = zeros(y_size);
y_pred_rev2 = zeros(y_size);
y_pred_rev3 = zeros(y_size);

impactError1_100ms = zeros(1, length(k));
impactError2_100ms = zeros(1, length(k));
impactError3_100ms = zeros(1, length(k));
shootingError1_100ms = zeros(1, length(k));
shootingError2_100ms = zeros(1, length(k));
shootingError3_100ms = zeros(1, length(k));

detectionThreshold = 1e4/100;
impactPoint = zeros(2, 1);
shootingPoint = zeros(2, 1);



% simulacao
for i = 1:length(k)
    % planta
    [plant, y_true(:, i)] = plant.run(u(:, i));

    % parar no impacto
    if y_true(2, i) < 0
        y_true = y_true(:, 1:i);
        impactPoint = y_true(: , end);
        shootingPoint = y_true(: , 1);
        break;
    end
end



"rodando loop para peridodo de amostragem 100ms"
i = 0;
% loop
while 1
    i = i + 1;
    % radar
    y(:, i) = y_true(:, i+detectionThreshold) + v_100ms(:, i+detectionThreshold);


    % kalman filter
    % primeira iteracao defini estado inicial de kalman filter
    if i == 1
        [kf1, y_est1(:, i), x_est1(:, i)] = kf1.setX0(y(:, i));
        [kf2, y_est2(:, i), x_est2(:, i)] = kf2.setX0(y(:, i));
        [kf3, y_est3(:, i), x_est3(:, i)] = kf3.setX0(y(:, i));
        continue;
    end
    % metodo 1
    % roda kalman filter para cada iteracao
    [kf1, y_est1(:, i), x_est1(:, i)] = kf1.run(y(:, i), u(:, i));
    % metodo 2
    % rodar kalman filter, ida e volta, para cada iteracao
    numberOfFiltering = 1;
    [y_est2, x_est2] = runMultipleKalmanFiltering(numberOfFiltering, y_est2, x_est2, y, u, kf2, kf_rev2, i);
    % metodo 3
    % rodar kalman filter, ida e volta, multiplas vezes para cada iteracao
    numberOfFiltering = 5;
    [y_est3, x_est3] = runMultipleKalmanFiltering(numberOfFiltering, y_est3, x_est3, y, u, kf3, kf_rev3, i);

    % previsao da trajetoria de impacto e disparo
    % metodo 1
    plant_pred = plant_pred.setX0(x_est1(:, i));
    plant_pred_rev = plant_pred_rev.setX0(x_est1(:, i));
    [y_pred1, impactError1_100ms(i)] = setImpactTrajectory(plant_pred, y_pred1, u, impactPoint, k);
    [y_pred_rev1, shootingError1_100ms(i)] = setShootingTrajectory(plant_pred_rev, y_pred_rev1, u, shootingPoint, k);
    % metodo 2
    plant_pred = plant_pred.setX0(x_est2(:, 1));
    plant_pred_rev = plant_pred_rev.setX0(x_est2(:, 1));
    [y_pred2, impactError2_100ms(i)] = setImpactTrajectory(plant_pred, y_pred2, u, impactPoint, k);
    [y_pred_rev2, shootingError2_100ms(i)] = setShootingTrajectory(plant_pred_rev, y_pred_rev2, u, shootingPoint, k);
    % metodo 3
    plant_pred = plant_pred.setX0(x_est3(:, 1));
    plant_pred_rev = plant_pred_rev.setX0(x_est3(:, 1));
    [y_pred3, impactError3_100ms(i)] = setImpactTrajectory(plant_pred, y_pred3, u, impactPoint, k);
    [y_pred_rev3, shootingError3_100ms(i)] = setShootingTrajectory(plant_pred_rev, y_pred_rev3, u, shootingPoint, k);
    

    % parar no impacto
    if y_true(2, i+detectionThreshold) < 0
        y = y(:, 1:i);
        y_est1 = y_est1(:, 1:i);
        y_est2 = y_est2(:, 1:i);
        y_est3 = y_est3(:, 1:i);
        break
    end


    % dynamic plot
    % dynamicPlot(y_true, y_est, y_pred1, y_pred2, impactError, shootingError, i);
end
"feito"



% ====================== periodo de amostragem 20ms ========================



% período de amostragem
Ts = 0.02;
k_20ms = (0:1:60/Ts)';
k = k_20ms;

% matrizes de espaco de estados
A = [1 0 Ts 0;
     0 1 0  Ts;
     0 0 1  0;
     0 0 0  1];

B = [Ts^2/2 0;
     0      Ts^2/2;
     Ts     0;
     0      Ts];

A_rev = inv(A);

B_rev = -A_rev*B;


% entrada do sistema
u = a'*ones(1, length(k));

% ruido de medicao
v_20ms = v(20:20:end);


% planta
plant = Plant(A, B, C, x_ss0);
plant_pred = Plant(A, B, C, zeros(length(A), 1));
plant_pred_rev = Plant(A_rev, B_rev, C, zeros(length(A), 1));

% Kalman filter
kf1 = KalmanFilter(A, B, C, Q, R, P1);
kf2 = KalmanFilter(A, B, C, Q, R, P2);
kf3 = KalmanFilter(A, B, C, Q, R, P3);
kf_rev2 = KalmanFilter(A_rev, B_rev, C, Q, R, P2);
kf_rev3 = KalmanFilter(A_rev, B_rev, C, Q, R, P3);

x_size = [length(A), length(k)];
% estado estimado
x_est1 = zeros(x_size);
x_est2 = zeros(x_size);
x_est3 = zeros(x_size);

y_size = [size(C, 1), length(k)];
% saida real
y_true = zeros(y_size);
% saida medida
y = zeros(y_size);
% saida estimada
y_est1 = zeros(y_size);
y_est2 = zeros(y_size);
y_est3 = zeros(y_size);
% saida prevista
y_pred1 = zeros(y_size);
y_pred2 = zeros(y_size);
y_pred3 = zeros(y_size);
y_pred_rev1 = zeros(y_size);
y_pred_rev2 = zeros(y_size);
y_pred_rev3 = zeros(y_size);

impactError1_20ms = zeros(1, length(k));
impactError2_20ms = zeros(1, length(k));
impactError3_20ms = zeros(1, length(k));
shootingError1_20ms = zeros(1, length(k));
shootingError2_20ms = zeros(1, length(k));
shootingError3_20ms = zeros(1, length(k));

detectionThreshold = 1e4/20;
impactPoint = zeros(2, 1);
shootingPoint = zeros(2, 1);



% simulacao
for i = 1:length(k)
    % planta
    [plant, y_true(:, i)] = plant.run(u(:, i));

    % parar no impacto
    if y_true(2, i) < 0
        y_true = y_true(:, 1:i);
        impactPoint = y_true(: , end);
        shootingPoint = y_true(: , 1);
        break;
    end
end



"rodando loop para peridodo de amostragem 20ms"
i = 0;
% loop
while 1
    i = i + 1;
    % radar
    y(:, i) = y_true(:, i+detectionThreshold) + v_20ms(:, i+detectionThreshold);


    % kalman filter
    % primeira iteracao defini estado inicial de kalman filter
    if i == 1
        [kf1, y_est1(:, i), x_est1(:, i)] = kf1.setX0(y(:, i));
        [kf2, y_est2(:, i), x_est2(:, i)] = kf2.setX0(y(:, i));
        [kf3, y_est3(:, i), x_est3(:, i)] = kf3.setX0(y(:, i));
        continue;
    end
    % metodo 1
    % roda kalman filter para cada iteracao
    [kf1, y_est1(:, i), x_est1(:, i)] = kf1.run(y(:, i), u(:, i));
    % metodo 2
    % rodar kalman filter, ida e volta, para cada iteracao
    numberOfFiltering = 1;
    [y_est2, x_est2] = runMultipleKalmanFiltering(numberOfFiltering, y_est2, x_est2, y, u, kf2, kf_rev2, i);
    % metodo 3
    % rodar kalman filter, ida e volta, multiplas vezes para cada iteracao
    numberOfFiltering = 5;
    [y_est3, x_est3] = runMultipleKalmanFiltering(numberOfFiltering, y_est3, x_est3, y, u, kf3, kf_rev3, i);

    % previsao da trajetoria de impacto e disparo
    % metodo 1
    plant_pred = plant_pred.setX0(x_est1(:, i));
    plant_pred_rev = plant_pred_rev.setX0(x_est1(:, i));
    [y_pred1, impactError1_20ms(i)] = setImpactTrajectory(plant_pred, y_pred1, u, impactPoint, k);
    [y_pred_rev1, shootingError1_20ms(i)] = setShootingTrajectory(plant_pred_rev, y_pred_rev1, u, shootingPoint, k);
    % metodo 2
    plant_pred = plant_pred.setX0(x_est2(:, 1));
    plant_pred_rev = plant_pred_rev.setX0(x_est2(:, 1));
    [y_pred2, impactError2_20ms(i)] = setImpactTrajectory(plant_pred, y_pred2, u, impactPoint, k);
    [y_pred_rev2, shootingError2_20ms(i)] = setShootingTrajectory(plant_pred_rev, y_pred_rev2, u, shootingPoint, k);
    % metodo 3
    plant_pred = plant_pred.setX0(x_est3(:, 1));
    plant_pred_rev = plant_pred_rev.setX0(x_est3(:, 1));
    [y_pred3, impactError3_20ms(i)] = setImpactTrajectory(plant_pred, y_pred3, u, impactPoint, k);
    [y_pred_rev3, shootingError3_20ms(i)] = setShootingTrajectory(plant_pred_rev, y_pred_rev3, u, shootingPoint, k);


    % parar no impacto
    if y_true(2, i+detectionThreshold) < 0
        y = y(:, 1:i);
        y_est1 = y_est1(:, 1:i);
        y_est2 = y_est2(:, 1:i);
        y_est3 = y_est3(:, 1:i);
        break
    end


    % dynamic plot
    % dynamicPlot(y_true, y_est, y_pred, y_pred_rev, impactError, shootingError, i);
end
"feito"



% ====================== periodo de amostragem 1ms ========================
% 
% 
% 
% % período de amostragem
% Ts = 0.001;
% k_1ms = (0:1:60/Ts)';
% k = k_1ms;
% 
% % matrizes de espaco de estados
% A = [1 0 Ts 0;
%      0 1 0  Ts;
%      0 0 1  0;
%      0 0 0  1];
% 
% B = [Ts^2/2 0;
%      0      Ts^2/2;
%      Ts     0;
%      0      Ts];
% 
% A_rev = inv(A);
% 
% B_rev = -A_rev*B;
% 
% 
% % entrada do sistema
% u = a'*ones(1, length(k));
% 
% % ruido de medicao
% v_1ms = v;
% 
% 
% % planta
% plant = Plant(A, B, C, x_ss0);
% plant_pred = Plant(A, B, C, zeros(length(A), 1));
% plant_pred_rev = Plant(A_rev, B_rev, C, zeros(length(A), 1));
% 
% % Kalman filter
% kf1 = KalmanFilter(A, B, C, Q, R, P1);
% kf2 = KalmanFilter(A, B, C, Q, R, P2);
% kf3 = KalmanFilter(A, B, C, Q, R, P3);
% kf_rev2 = KalmanFilter(A_rev, B_rev, C, Q, R, P2);
% kf_rev3 = KalmanFilter(A_rev, B_rev, C, Q, R, P3);
% 
% x_size = [length(A), length(k)];
% % estado estimado
% x_est1 = zeros(x_size);
% x_est2 = zeros(x_size);
% x_est3 = zeros(x_size);
% 
% y_size = [size(C, 1), length(k)];
% % saida real
% y_true = zeros(y_size);
% % saida medida
% y = zeros(y_size);
% % saida estimada
% y_est1 = zeros(y_size);
% y_est2 = zeros(y_size);
% y_est3 = zeros(y_size);
% % saida prevista
% y_pred1 = zeros(y_size);
% y_pred2 = zeros(y_size);
% y_pred3 = zeros(y_size);
% y_pred_rev1 = zeros(y_size);
% y_pred_rev2 = zeros(y_size);
% y_pred_rev3 = zeros(y_size);
% 
% impactError1_1ms = zeros(1, length(k));
% impactError2_1ms = zeros(1, length(k));
% impactError3_1ms = zeros(1, length(k));
% shootingError1_1ms = zeros(1, length(k));
% shootingError2_1ms = zeros(1, length(k));
% shootingError3_1ms = zeros(1, length(k));
% 
% detectionThreshold = 1e4;
% impactPoint = zeros(2, 1);
% shootingPoint = zeros(2, 1);
% 
% 
% 
% % simulacao
% for i = 1:length(k)
%     % planta
%     [plant, y_true(:, i)] = plant.run(u(:, i));
% 
%     % parar no impacto
%     if y_true(2, i) < 0
%         y_true = y_true(:, 1:i);
%         impactPoint = y_true(: , end);
%         shootingPoint = y_true(: , 1);
%         break;
%     end
% end
% 
% 
%
% "rodando loop para peridodo de amostragem 1ms"
% i = 0;
% % loop
% while 1
%     i = i + 1;
%     % radar
%     y(:, i) = y_true(:, i+detectionThreshold) + v_1ms(:, i+detectionThreshold);
% 
% 
%     % kalman filter
%     % primeira iteracao defini estado inicial de kalman filter
%     if i == 1
%         [kf1, y_est1(:, i), x_est1(:, i)] = kf1.setX0(y(:, i));
%         [kf2, y_est2(:, i), x_est2(:, i)] = kf2.setX0(y(:, i));
%         [kf3, y_est3(:, i), x_est3(:, i)] = kf3.setX0(y(:, i));
%         continue;
%     end
%     % metodo 1
%     % roda kalman filter para cada iteracao
%     [kf1, y_est1(:, i), x_est1(:, i)] = kf1.run(y(:, i), u(:, i));
%     % metodo 2
%     % rodar kalman filter, ida e volta, para cada iteracao
%     numberOfFiltering = 1;
%     [y_est2, x_est2] = runMultipleKalmanFiltering(numberOfFiltering, y_est2, x_est2, y, u, kf2, kf_rev2, i);
%     % metodo 3
%     % rodar kalman filter, ida e volta, multiplas vezes para cada iteracao
%     numberOfFiltering = 5;
%     [y_est3, x_est3] = runMultipleKalmanFiltering(numberOfFiltering, y_est3, x_est3, y, u, kf3, kf_rev3, i);
% 
%     % previsao da trajetoria de impacto e disparo
%     % metodo 1
%     plant_pred = plant_pred.setX0(x_est1(:, i));
%     plant_pred_rev = plant_pred_rev.setX0(x_est1(:, i));
%     [y_pred1, impactError1_1ms(i)] = setImpactTrajectory(plant_pred, y_pred1, u, impactPoint, k);
%     [y_pred_rev1, shootingError1_1ms(i)] = setShootingTrajectory(plant_pred_rev, y_pred_rev1, u, shootingPoint, k);
%     % metodo 2
%     plant_pred = plant_pred.setX0(x_est2(:, 1));
%     plant_pred_rev = plant_pred_rev.setX0(x_est2(:, 1));
%     [y_pred2, impactError2_1ms(i)] = setImpactTrajectory(plant_pred, y_pred2, u, impactPoint, k);
%     [y_pred_rev2, shootingError2_1ms(i)] = setShootingTrajectory(plant_pred_rev, y_pred_rev2, u, shootingPoint, k);
%     % metodo 3
%     plant_pred = plant_pred.setX0(x_est3(:, 1));
%     plant_pred_rev = plant_pred_rev.setX0(x_est3(:, 1));
%     [y_pred3, impactError3_1ms(i)] = setImpactTrajectory(plant_pred, y_pred3, u, impactPoint, k);
%     [y_pred_rev3, shootingError3_1ms(i)] = setShootingTrajectory(plant_pred_rev, y_pred_rev3, u, shootingPoint, k);
% 
% 
%     % parar no impacto
%     if y_true(2, i+detectionThreshold) < 0
%         y = y(:, 1:i);
%         y_est1 = y_est1(:, 1:i);
%         y_est2 = y_est2(:, 1:i);
%         y_est3 = y_est3(:, 1:i);
%         break
%     end
% 
% 
%     % dynamic plot
%     % dynamicPlot(y_true, y_est, y_pred, y_pred_rev, impactError, shootingError, i);
% end
% "feito"



% full plot
% fullplot(y, y_true, y_est);

figure

subplot(2, 1, 1)
plot(k_100ms, impactError1_100ms, 'r', ...
     k_100ms, impactError2_100ms, 'g', ...
     k_100ms, impactError3_100ms, 'b', ...
     'LineWidth', 2)
xlabel('k-detectionThreshold'), ylabel('erro (m)')
title('Erro de Impacto em 100ms')
legend('Metodo 1', 'Metodo 2', 'Metodo 3')
grid on

subplot(2, 1, 2)
plot(k_20ms, impactError1_20ms, 'r', ...
     k_20ms, impactError2_20ms, 'g', ...
     k_20ms, impactError3_20ms, 'b', ...
     'LineWidth', 2)
xlabel('k-detectionThreshold'), ylabel('erro (m)')
title('Erro de Impacto em 20ms')
legend('Metodo 1', 'Metodo 2', 'Metodo 3')
grid on

fig=gcf;
fig.Units='normalized';
fig.OuterPosition=[0 0 1 1];



figure

subplot(2, 1, 1)
plot(k_100ms, shootingError1_100ms, 'r', ...
     k_100ms, shootingError2_100ms, 'g', ...
     k_100ms, shootingError3_100ms, 'b', ...
     'LineWidth', 2)
xlabel('k-detectionThreshold'), ylabel('erro (m)')
title('Erro de Disparo em 100ms')
legend('Metodo 1', 'Metodo 2', 'Metodo 3')
grid on

subplot(2, 1, 2)
plot(k_20ms, shootingError1_20ms, 'r', ...
     k_20ms, shootingError2_20ms, 'g', ...
     k_20ms, shootingError3_20ms, 'b', ...
     'LineWidth', 2)
xlabel('k-detectionThreshold'), ylabel('erro (m)')
title('Erro de Disparo em 20ms')
legend('Metodo 1', 'Metodo 2', 'Metodo 3')
grid on

fig=gcf;
fig.Units='normalized';
fig.OuterPosition=[0 0 1 1];



% rodar kalman filter do inicio ate o fim e do fim ate o inicio,
% n vezes, para cada iteracao
function [y_est, x_est] = runMultipleKalmanFiltering(numberOfFiltering, y_est, x_est, y, u, kf, kf_rev, i)
    for n = 1:numberOfFiltering
        kf.x = x_est(:, 1);
        for j = 2:i
            [kf, y_est(:, j), x_est(:, j)] = kf.run(y(:, j), u(:, j));
        end
        kf_rev.x = x_est(:, i);
        for j = 2:i
            l = i - j + 1;
            [kf_rev, y_est(:, l), x_est(:, l)] = kf_rev.run(y(:, l), u(:, l));
        end
    end
end



% previsao da trajetoria de impacto
function [y_pred1, impactError] = setImpactTrajectory(plant_pred, y_pred1, u, impactPoint, k)
    for j = 1:length(k)
        [plant_pred, y_pred1(:, j)] = plant_pred.run(u(:, j));
        if y_pred1(2, j) < 0
            y_pred1 = y_pred1(:, 1:j);
            impactError = abs(mean(y_pred1(:, j) - impactPoint));
            return
        end
    end
    impactError = 1e4;
end



% previsao da trajetoria de impacto
function [y_pred2, shootingError] = setShootingTrajectory(plant_pred_rev, y_pred2, u, shootingPoint, k)
    for j = 1:length(k)
        [plant_pred_rev, y_pred2(:, j)] = plant_pred_rev.run(u(:, j));
        if y_pred2(2, j) < 0
            y_pred2 = y_pred2(:, 1:j);
            shootingError = abs(mean(y_pred2(:, j) - shootingPoint));
            return
        end
    end
    shootingError = 1e4;
end



function dynamicPlot(y_true, y_est, y_pred1, y_pred2, impactError, shootingError, i)
    if i == 1
        figure
    end

    subplot(2, 2, [1 2])
    plot(y_true(1, :), y_true(2, :), 'b*', ...
         y_est(1, :), y_est(2, :), 'rx', ...
         y_pred1(1, :), y_pred1(2, :), 'r--',...
         y_pred2(1, :), y_pred2(2, :), 'g--', 'LineWidth', 1)
    axis([-1e3, 14e3, -5e2, 4e3])
    legend('Medida', 'Filtrada', 'Trajetoria de Impacto', 'Trajetoria de Disparo')
    grid on

    subplot(2, 2, 3)
    plot(impactError, 'LineWidth', 2)
    xlabel('k'), ylabel('erro (m)')
    title('Erro de Impacto')
    grid on

    subplot(2, 2, 4)
    plot(shootingError, 'LineWidth', 2)
    xlabel('k'), ylabel('erro (m)')
    title('Erro de Disparo')
    grid on

    pause(0);
    fig=gcf;
    fig.Units='normalized';
    fig.OuterPosition=[0 0 1 1];
end



function fullplot(y, y_true, y_est)
    figure

    plot(y(1, :), y(2, :), 'k', ...
        y_true(1, :), y_true(2, :), 'b', ...
        y_est(1, :), y_est(2, :), 'r--', 'LineWidth', 2)
    xlabel('x (m)'), ylabel('y (m)')
    title('Resposta do Kalman filter')
    legend('Medida', 'Real', 'Filtrada')
    axis tight
    grid on

    fig=gcf;
    fig.Units='normalized';
    fig.OuterPosition=[0 0 1 1];
end


