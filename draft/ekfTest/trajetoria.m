clear; clc; close all;

% variaveis de modelo
T = 0.05;
detecThresh = ceil(40/T);
predTime = [1 5];
g = 9.81;
v0 = [180 180 250];
gama = 0;
a = [0 0 -g];
p_floor = 0;
u = a';
testState = 7;
numOfFiltering = 5;


% matriz covariancia do ruido de processo
Q = zeros(7);

% valor inicial de matriz covariancia de estado
% valores nao otimos
P = [1e8*eye(3) zeros(3) zeros(3, 1);
     zeros(3) 1e8*eye(3) zeros(3, 1);
     zeros(1, 3) zeros(1, 3) 0];

% matriz covariancia do ruido de medicao
sigma2_n = 1e2;
R = sigma2_n*eye(3);   % roubo
R = 0.1*R;             % R chutado para 10 vezes menor

% EKF params
[f, h, F, H, f_rev, F_rev] = getParamsEkf(T);

% Extended Kalman filter
kf = ExtendedKalmanFilter(f, h, F, H, Q, R, P);

% trajetoria real
x_true = [0*ones(1, 3) v0 gama]';
x_aux = x_true;
i = 0;
while 1
    i = i + 1;
    x_aux = f(x_aux, u);
    y_true(:, i) = h(x_aux, u);
    if y_true(3, i) < p_floor
        break
    end
end

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

[kf, y_est(:, 1), x_est] = kf.setInitialState([y(:, 1)' 0 0 200 0]);
for i = 2:arrayLength
    [kf, y_est(:, i), x_est] = kf.run(y(:, i), u);
end


figure
plot3(y_true(1, :), y_true(2, :), y_true(3, :), 'b', ...
      y_est(1, :), y_est(2, :), y_est(3, :), 'r*', 'LineWidth', 2)
lgd = legend("Real", "Filtrada", 'FontSize', 15);
title("Trajetórias Real e Filtrada com Medições Próximas ao Ponto de Impacto", 'FontSize', 16)
xlabel('x (m)', 'FontSize', 15), ylabel('y (m)', 'FontSize', 15), zlabel('z (m)', 'FontSize', 15)
grid on
axis tight
set(lgd, 'ItemTokenSize', [30, 50]);
set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1], 'Color', [1 1 1])