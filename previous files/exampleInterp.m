clear; clc; close all;


t_max = 4*4;
t = 0.06*(0:t_max);
y_new = -t.^2 + t;

figure
plot(t(1:end - 1), y_new(1:end - 1), 'b')
hold on
plot(t(1:end - 1), y_new(1:end - 1), 'b|', ...
     t(1:4:end - 1), y_new(1:4:end - 1), 'g*', ...
     t(end), y_new(end), 'r*', 'LineWidth', 2)
title('Posicao Interpolada')
legend('Interpolada 5ms', 'Interpolada 5ms', 'Medida 20ms', 'estimada (kalman filter)', 'ItemHitFcn', @cb_legend)
grid on
set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1], 'Color', [1 1 1]);