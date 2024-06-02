clear; close all; clc;


% variaveis de mruv
vel0 = [250 250]; a = [0 -9.81];

% período de amostragem
Ts = 0.1;
k = (0:1:60/Ts)';

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
v = randn(size(C, 1), length(k));
v = v - mean(v, 2)*ones(1, length(k));
v = v.*sqrt(sigma2_v)./std(v, 0, 2);


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
P = 1e1*eye(length(A));

% planta
x_ss0 = [0; 0; vel0(1); vel0(2)];
plant = Plant(A, B, C, x_ss0);
plant_pred = Plant(A, B, C, zeros(length(A), 1));
plant_pred_rev = Plant(A_rev, B_rev, C, zeros(length(A), 1));

% Kalman filter
kf = KalmanFilter(A, B, C, Q, R, P);
kf_rev = KalmanFilter(A_rev, B_rev, C, Q, R, P);

% estado estimado
x_est = zeros(length(A), length(k));

y_size = [size(C, 1), length(k)];
% saida real
y_true = zeros(y_size);
% saida medida
y = zeros(y_size);
% saida estimada
y_est = zeros(y_size);
% saida prevista
y_pred1 = zeros(y_size);
y_pred2 = zeros(y_size);

impactError = zeros(1, length(k));
shootingError = zeros(1, length(k));

detectionThreshold = 100;



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



i = 0;
% loop
while 1
    i = i + 1;
    % radar
    y(:, i) = y_true(:, i+detectionThreshold) + v(:, i);


    % kalman filter
    % primeira iteracao defini estado inicial de kalman filter
    if i == 1
        [kf, y_est(:, i), x_est(:, i)] = kf.setX0(y(:, i));
        continue;
    end
    % metodo 1
    % roda kalman filter para cada iteracao
    % [kf, y_est(:, i), x_est(:, i)] = kf.run(y(:, i), u(:, i));
    % metodo 2 e metodo 3
    % rodar kalman filter multiplas vezes para cada iteracao
    numberOfFiltering = 5;
    [y_est, x_est] = runMultipleKalmanFiltering(numberOfFiltering, y_est, x_est, y, u, kf, kf_rev, i);


    % previsao da trajetoria de impacto e de disparo
    % metodo 1
    % plant_pred = plant_pred.setX0(x_est(:, i));
    % plant_pred_rev = plant_pred_rev.setX0(x_est(:, i));
    % metodo 2 e metodo 3
    plant_pred = plant_pred.setX0(x_est(:, 1));
    plant_pred_rev = plant_pred_rev.setX0(x_est(:, 1));
    
    [y_pred1, impactError(i)] = setImpactTrajectory(plant_pred, y_pred1, u, impactPoint, k);
    [y_pred2, shootingError(i)] = setShootingTrajectory(plant_pred_rev, y_pred2, u, shootingPoint, k);
    

    % parar no impacto
    if y_est(2, i) < 0 || i == length(y_true)-detectionThreshold
        y = y(:, 1:i);
        y_est = y_est(:, 1:i);
        break
    end


    % dynamic plot
    % dynamicPlot(y_true, y_est, y_pred1, y_pred2, impactError, shootingError, i);
end

% full plot
% fullplot(y, y_true, y_est);

impacterror_2s_3s_4s = [impactError(2/Ts), impactError(3/Ts), impactError(4/Ts)]
shootingerror_2s_3s_4s = [shootingError(2/Ts), shootingError(3/Ts), shootingError(4/Ts)]


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


