clear; clc; close all;

% periodo de amostragem
T = 0.05;
% tempo de simulacao
t = 0:T:30;

% entrada de controle
u = [0; 0; -9.81];


% saida real
y_true = zeros(3, length(t));
y_true_rev = zeros(3, length(t));

% obter trajetoria real diretamente
x_true = [0; 0; 0; 100; 100; 180; 0.000245];
for k = 2:length(t)
    x_true = f(x_true, u, T);
    y_true(:, k) = h(x_true, u);
end

% obter trajetoria real reversamente
x_true_rev = x_true;
for k = 1:length(t)
    y_true_rev(:, k) = h(x_true_rev, u);
    x_true_rev = f_rev(x_true_rev, u, T);
end

% erro entre trajetorias
y_error = zeros(3, length(t));
for k = 1:length(t)
    l = length(t) - k + 1;
    y_error(:, k) = y_true(:, k) - y_true_rev(:, l);
end

% plotar trajetoria
figure
plot3(y_true(1, :), y_true(2, :), y_true(3, :), 'b', 'LineWidth', 2)
hold on
plot3(y_true_rev(1, :), y_true_rev(2, :), y_true_rev(3, :), 'r--', 'LineWidth', 2)
title("Modelagem Direta e Reversa Nao Linear")
legend("Modelagem Direta (+T)", "Modelagem Reversa (-T)")
xlabel("x (m)"), ylabel("y (m)"), zlabel("z (m)")
grid on
set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1], 'Color', [1 1 1])

% plotar erro entre modelagens
figure
subplot(3, 1, 1)
plot(t, y_error(1, :), 'k', 'LineWidth', 2)
title("Erro Entre Modelagem Direta e Reversa - Eixo x")
xlabel('tempo (s)'), ylabel('erro (m)')
grid on

subplot(3, 1, 2)
plot(t, y_error(2, :), 'k', 'LineWidth', 2)
title("Erro Entre Modelagem Direta e Reversa - Eixo y")
xlabel('tempo (s)'), ylabel('erro (m)')
grid on

subplot(3, 1, 3)
plot(t, y_error(3, :), 'k', 'LineWidth', 2)
title("Erro Entre Modelagem Direta e Reversa - Eixo z")
xlabel('tempo (s)'), ylabel('erro (m)')
grid on

set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1], 'Color', [1 1 1])


% funcao transicao de estado
function x_pred = f(x, u, T)
    px = x(1); py = x(2); pz = x(3); vx = x(4); vy = x(5); vz = x(6); gama = x(7);
    g = -u(3);

    v = sqrt(vx^2 + vy^2 + vz^2);
    h = 1.01*1e-4;

    x_pred = [px + vx*T - gama*v*vx*T^2/4;
         py + vy*T - gama*v*vy*T^2/4;
         pz + vz*T - gama*v*vz*T^2/4 - g*T^2/2;
         vx + (T^2*((gama^2*vx*vy^2)/4 + (gama*vx*((gama*v)/2 + (gama*vx^2)/(2*v))*v)/2 + (gama*vx*vz*(g + (gama*vz*v)/2))/(2*v) + (gama*h*vx*vz*v)/2))/2 - (T*gama*vx*v)/2;
         vy + (T^2*((gama^2*vx^2*vy)/4 + (gama*vy*((gama*v)/2 + (gama*vy^2)/(2*v))*v)/2 + (gama*vy*vz*(g + (gama*vz*v)/2))/(2*v) + (gama*h*vy*vz*v)/2))/2 - (T*gama*vy*v)/2;
         vz + (T^2*(((gama*v)/2 + (gama*vz^2)/(2*v))*(g + (gama*vz*v)/2) + (gama^2*vx^2*vz)/4 + (gama^2*vy^2*vz)/4 + (gama*h*vz^2*v)/2))/2 - T*(g + (gama*vz*v)/2);
         gama + (T^2*(gama*h^2*vz^2 + gama*h*(g + (gama*vz*v)/2)))/2 - T*gama*h*vz];
end

% funcao transicao de estado
function x_pred = f_rev(x, u, T)
    px = x(1); py = x(2); pz = x(3); vx = x(4); vy = x(5); vz = x(6); gama = x(7);
    g = -u(3);

    v = sqrt(vx^2 + vy^2 + vz^2);
    h = 1.01*1e-4;

    T = -T;

    x_pred = [px + vx*T - gama*v*vx*T^2/4;
         py + vy*T - gama*v*vy*T^2/4;
         pz + vz*T - gama*v*vz*T^2/4 - g*T^2/2;
         vx + (T^2*((gama^2*vx*vy^2)/4 + (gama*vx*((gama*v)/2 + (gama*vx^2)/(2*v))*v)/2 + (gama*vx*vz*(g + (gama*vz*v)/2))/(2*v) + (gama*h*vx*vz*v)/2))/2 - (T*gama*vx*v)/2;
         vy + (T^2*((gama^2*vx^2*vy)/4 + (gama*vy*((gama*v)/2 + (gama*vy^2)/(2*v))*v)/2 + (gama*vy*vz*(g + (gama*vz*v)/2))/(2*v) + (gama*h*vy*vz*v)/2))/2 - (T*gama*vy*v)/2;
         vz + (T^2*(((gama*v)/2 + (gama*vz^2)/(2*v))*(g + (gama*vz*v)/2) + (gama^2*vx^2*vz)/4 + (gama^2*vy^2*vz)/4 + (gama*h*vz^2*v)/2))/2 - T*(g + (gama*vz*v)/2);
         gama + (T^2*(gama*h^2*vz^2 + gama*h*(g + (gama*vz*v)/2)))/2 - T*gama*h*vz];
end

function y_est = h(x, u)
    y_est = x(1:3);
end