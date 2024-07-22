clear; clc; close all;

addpath(genpath('utils'));

axis = 1;

Ts_ref = 0.005;
g = 9.81;
v0 = [180 180 250];
u = [0 0 -g]';
p_floor = 0;
predTime = [1 5];
sigma2_n = [1e2; 1e2; 1e2];

% trajetoria medida 5ms
[y_true, ~, ~] = setYTrue(Ts_ref, v0, g, u, p_floor);
y_5ms = createTrajectory(y_true, sigma2_n, predTime, Ts_ref);
% trajetoria medida 20ms
y_20ms = y_5ms(:, 1:4:end);
% trajetoria de 20ms interpolada para 5ms
y_20ms_Interp_5ms = interp1(1:length(y_20ms(axis, :)), y_20ms(axis, :), 1:1/4:length(y_20ms(axis, :)), 'linear');

% espectro
Y_5ms = fft(y_5ms(axis, :));
Y_20ms = fft(y_20ms(axis, :));
Y_20ms_Interp_5ms = fft(y_20ms_Interp_5ms);
% magnitude
Y_5ms_mag = abs(Y_5ms);
Y_20ms_mag = abs(Y_20ms);
Y_20ms_Interp_5ms_mag = abs(Y_20ms_Interp_5ms);

% plot
figure
subplot(2, 1, 1)
plot(Y_5ms_mag, 'LineWidth', 2)
hold on
plot(Y_20ms_mag, 'LineWidth', 2)
title("Espectro do Sinal Medido de Diferentes Ts")
xlabel("frequência (Hz)"), ylabel("magnitude")
legend("5ms", "20ms", "ItemHitFcn", @cb_legend)
grid on

subplot(2, 1, 2)
plot(Y_5ms_mag, 'LineWidth', 2)
hold on
plot(Y_20ms_Interp_5ms, 'LineWidth', 2)
title("Espectro do Sinal Medido Com e Sem Interpolação")
xlabel("frequência (Hz)"), ylabel("magnitude")
legend("5ms Não Interpolado", "20ms Interpolado para 5ms", "ItemHitFcn", @cb_legend)
grid on

set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1], 'Color', [1 1 1])