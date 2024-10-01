clear; clc; close all;

% variaveis de modelo
detecThresh = 200;
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
     zeros(1, 3) zeros(1, 3) 1e-10];

% matriz covariancia do ruido de medicao
sigma2_n = 1e2;
R = sigma2_n*eye(3);   % roubo
R = 0.1*R;             % R chutado para 10 vezes menor

% model params
[f_model, ~, ~, ~, f_model_rev, ~] = getParamsEkf(0.1);
% EKF params
[f, h, F, H, f_rev, F_rev] = getParamsEkf(T);

% Extended Kalman filter
kf = ExtendedKalmanFilter(f, h, F, H, Q, R, P);
kf_rev = ExtendedKalmanFilter(f_rev, h, F_rev, H, Q, R, P);


% trajetoria real
x_true = [0*ones(1, 3) v0 gama]';
x_aux = x_true;
i = 0;
while 1
    i = i + 1;
    x_aux = f(x_aux, u);
    if i == detecThresh + 1
        x_true_start = x_aux;
    elseif i == detecThresh + 100
        x_true_last = x_aux;
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

arrayLength = length(y_true);

numOfSim = 10;
for i = 1:numOfSim
    % ruido de medicao
    n = randn(3, arrayLength);
    n = n - mean(n, 2)*ones(1, arrayLength);
    n = n*sqrt(sigma2_n)./std(n, 0, 2);

    % trajetoria medida
    arrayLength = ceil(predTime(2)/T);
    y = y_true(:, detecThresh + 1:detecThresh+arrayLength) + n(:, 1:arrayLength);

    % estado estimado inicial
    x0 = [y(:, 1); 0; 0; 200; 2e-4];

    % metodo 1
    method = 1;
    [impErr1, shoErr1, state_error1] = loop(kf, kf_rev, x0, T, predTime, f_model, f_model_rev, h, method, y, impPt, shoPt, u, p_floor, P, x_true_start, x_true_last);

    % metodo 2
    method = 2;
    [impErr2, shoErr2, state_error2] = loop(kf, kf_rev, x0, T, predTime, f_model, f_model_rev, h, method, y, impPt, shoPt, u, p_floor, P, x_true_start, x_true_last);

    % metodo 3
    method = 3;
    [impErr3, shoErr3, state_error3] = loop(kf, kf_rev, x0, T, predTime, f_model, f_model_rev, h, method, y, impPt, shoPt, u, p_floor, P, x_true_start, x_true_last);

    % calcular media de erro
    if i == 1
        impErrAvg1 = impErr1; shoErrAvg1 = shoErr1;
        state_error1_avg = state_error1;

        impErrAvg2 = impErr2; shoErrAvg2 = shoErr2;
        state_error2_avg = state_error2;

        impErrAvg3 = impErr3; shoErrAvg3 = shoErr3;
        state_error3_avg = state_error3;
    else
        impErrAvg1 = impErrAvg1 + impErr1; shoErrAvg1 = shoErrAvg1 + shoErr1;
        state_error1_avg = state_error1_avg + state_error1;

        impErrAvg2 = impErrAvg2 + impErr2; shoErrAvg2 = shoErrAvg2 + shoErr2;
        state_error2_avg = state_error2_avg + state_error2;

        impErrAvg3 = impErrAvg3 + impErr3; shoErrAvg3 = shoErrAvg3 + shoErr3;
        state_error3_avg = state_error3_avg + state_error3;
    end
end
impErrAvg1 = impErrAvg1/numOfSim; shoErrAvg1 = shoErrAvg1/numOfSim;
state_error1_avg = state_error1_avg/numOfSim

impErrAvg2 = impErrAvg2/numOfSim; shoErrAvg2 = shoErrAvg2/numOfSim;
state_error2_avg = state_error2_avg/numOfSim

impErrAvg3 = impErrAvg3/numOfSim; shoErrAvg3 = shoErrAvg3/numOfSim;
state_error3_avg = state_error3_avg/numOfSim


% plot dos erros
figure
plot(impErrAvg1, 'LineWidth', 2)
hold on
plot(impErrAvg2, 'LineWidth', 2)
hold on
plot(impErrAvg3, 'LineWidth', 2)
title('impacto')
legend('metodo 1', 'metodo 2', 'metodo 3')
grid on
xlabel('amostras'), ylabel('erro (m)')
ylim([0 1e3])
% set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1], 'Color', [1 1 1])

figure
plot(shoErrAvg1, 'LineWidth', 2)
hold on
plot(shoErrAvg2, 'LineWidth', 2)
hold on
plot(shoErrAvg3, 'LineWidth', 2)
title('disparo')
legend('metodo 1', 'metodo 2', 'metodo 3')
grid on
xlabel('amostras'), ylabel('erro (m)')
ylim([0 1e3])
% set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1], 'Color', [1 1 1])