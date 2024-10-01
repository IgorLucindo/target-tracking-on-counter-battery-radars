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
     zeros(1, 3) zeros(1, 3) 1e-8];

% matriz covariancia do ruido de medicao
sigma2_n = 1e2;
R = sigma2_n*eye(3);   % roubo
R = 0.1*R;             % R chutado para 10 vezes menor

% model params
[f_model, ~, ~, ~, f_model_rev, ~] = getParamsEkf(0.1);
% EKF params
[f, h, F, H, f_rev, F_rev] = getParamsEkf(T);

% Extended Kalman filter
ekf = ExtendedKalmanFilter(f, h, F, H, Q, R, P);
ekf_rev = ExtendedKalmanFilter(f_rev, h, F_rev, H, Q, R, P);


% gama array
gama_array_true = [];

% trajetoria real
x_true = [0*ones(1, 3) v0 gama]';
x_aux = x_true;
i = 0;
while 1
    i = i + 1;
    x_aux = f(x_aux, u);
    if i == detecThresh + 1
        cheating_gama = x_aux(7);
    end
    gama_array_true = [gama_array_true x_aux(7)];
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
impPtPred = zeros(3, 1);
shoPtPred = zeros(3, 1);


gama_array_avg = [];

numOfFiltering = 5;
numOfSim = 100;
for l = 1:numOfSim
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

    gama_array = [];

    n = 0;
    x_est = [y(:, 1); 0; 0; 200; 0.1e-4];
    % loop
    for k = 1:numOfFiltering
        % ida
        if k == 1
            ekf.P = P;
        elseif k > 1
            ekf.P = ekf_rev.P;
            ekf.P(4:6, 1:3) = -ekf.P(4:6, 1:3);
            ekf.P(1:3, 4:6) = -ekf.P(1:3, 4:6);
            ekf.P(4:6, 7) = -ekf.P(4:6, 7);
            ekf.P(7, 4:6) = -ekf.P(7, 4:6);
        end
        ekf.P(7, 7) = P(7, 7)/2^(2*k-2);
    
        [ekf, y_est(:, 1), x_est] = ekf.setInitialState(x_est);
        aux_array = [];
        for i = 2:arrayLength
            aux_array = [aux_array x_est(7)];
            [ekf, y_est(:, i), x_est] = ekf.run(y(:, i), u);
        end
        gama_array = [gama_array; aux_array];
    
        % volta
        if k == 1
            ekf_rev.P = P;
        elseif k > 1
            ekf_rev.P = ekf.P;
            ekf_rev.P(4:6, 1:3) = -ekf_rev.P(4:6, 1:3);
            ekf_rev.P(1:3, 4:6) = -ekf_rev.P(1:3, 4:6);
            ekf_rev.P(4:6, 7) = -ekf_rev.P(4:6, 7);
            ekf_rev.P(7, 4:6) = -ekf_rev.P(7, 4:6);
        end
        ekf_rev.P(7, 7) = P(7, 7)/2^(2*k-1);
    
        [ekf_rev, ~, x_est] = ekf_rev.setInitialState(x_est);
        aux_array = [];
        for i = 2:arrayLength
            j = arrayLength - i + 1;
            aux_array = [aux_array x_est(7)];
            [ekf_rev, ~, x_est] = ekf_rev.run(y(:, j), u);
        end
        gama_array = [gama_array; flip(aux_array)];
    end

    % calcula media de gama
    if length(gama_array_avg) == 0
        gama_array_avg = gama_array;
    else
        gama_array_avg = gama_array_avg + gama_array;
    end
end
gama_array_avg = gama_array_avg/numOfSim;

% plot de gama
arrayLength = length(gama_array_avg(1, :));
t = T:T:T*arrayLength;
plotlegend = ["real"];

figure
plot(t, gama_array_true(detecThresh + 1:detecThresh + arrayLength), 'LineWidth', 2)
for k = 1:numOfFiltering
    hold on
    plot(t, gama_array_avg(2*k-1, :), 'LineWidth', 2)
    hold on
    plot(t, gama_array_avg(2*k, :), 'LineWidth', 2)
    plotlegend(2*k) = "ida " + k;
    plotlegend(2*k + 1) = "volta " + k;
end
legend(plotlegend, 'ItemHitFcn', @cb_legend)
title('estimacao de gama')
grid on
ylim([0.5e-4 2.2e-4])
set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1], 'Color', [1 1 1])

function cb_legend(~, evt)
    if strcmp(evt.Peer.Visible, 'on')
        evt.Peer.Visible = 'off';
    else 
        evt.Peer.Visible = 'on';
    end
end