clear; clc; close all;

% variaveis de modelo
T = 0.05;
g = 9.81;
v0 = [180 180 250];
gama1 = 0;
gama2 = 1e-4;
a = [0 0 -g];
p_floor = 0;
u = a';

% model params
[f, h, ~, ~, ~, ~] = getParamsEkf(T);

% estado inicial real
x_true1 = [0*ones(1, 3) v0 gama1]';
x_true2 = [0*ones(1, 3) v0 gama2]';
% trajetoria real
x_aux = x_true1;
i = 0;
while 1
    i = i + 1;
    x_aux = f(x_aux, u);
    y_true1(:, i) = h(x_aux, u);
    if y_true1(3, i) < p_floor
        break
    end
end
x_aux = x_true2;
i = 0;
while 1
    i = i + 1;
    x_aux = f(x_aux, u);
    y_true2(:, i) = h(x_aux, u);
    if y_true2(3, i) < p_floor
        break
    end
end

plot3(y_true1(1, :), y_true1(2, :), y_true1(3, :), 'b', 'LineWidth', 2)
hold on
plot3(y_true2(1, :), y_true2(2, :), y_true2(3, :), 'r', 'LineWidth', 2)
grid on
title('Trajetórias Simuladas')
legend('sem arrasto', 'com arrasto')
set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1], 'Color', [1 1 1])