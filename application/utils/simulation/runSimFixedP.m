% rodar simulacoes com P fixo
function [impErrArray, shoErrArray] = runSimFixedP(v, y_true, g, u, impPt, shoPt, pos_floor, predTime, Q, R)
    % metodo
    method1 = 1;
    method2 = 2;
    method3 = 3;
    % erro de covariancia
    P = 1e8*eye(6);
    % periodo de amostragem
    Ts_1ms = 0.001;
    Ts_5ms = 0.005;
    Ts_10ms = 0.01;
    Ts_15ms = 0.015;
    Ts_20ms = 0.02;
    Ts_25ms = 0.025;
    Ts_30ms = 0.03;
    Ts_base = Ts_1ms;

    [~, ~, impErrArray(1, 1, :), shoErrArray(1, 1, :)] = runSim(Ts_base, Ts_1ms, method1, P, v, y_true, g, u, impPt, shoPt, pos_floor, predTime, Q, R);
    [~, ~, impErrArray(1, 2, :), shoErrArray(1, 2, :)] = runSim(Ts_base, Ts_5ms, method1, P, v, y_true, g, u, impPt, shoPt, pos_floor, predTime, Q, R);
    [~, ~, impErrArray(1, 3, :), shoErrArray(1, 3, :)] = runSim(Ts_base, Ts_10ms, method1, P, v, y_true, g, u, impPt, shoPt, pos_floor, predTime, Q, R);
    [~, ~, impErrArray(1, 4, :), shoErrArray(1, 4, :)] = runSim(Ts_base, Ts_15ms, method1, P, v, y_true, g, u, impPt, shoPt, pos_floor, predTime, Q, R);
    [~, ~, impErrArray(1, 5, :), shoErrArray(1, 5, :)] = runSim(Ts_base, Ts_20ms, method1, P, v, y_true, g, u, impPt, shoPt, pos_floor, predTime, Q, R);
    [~, ~, impErrArray(1, 6, :), shoErrArray(1, 6, :)] = runSim(Ts_base, Ts_25ms, method1, P, v, y_true, g, u, impPt, shoPt, pos_floor, predTime, Q, R);
    [~, ~, impErrArray(1, 7, :), shoErrArray(1, 7, :)] = runSim(Ts_base, Ts_30ms, method1, P, v, y_true, g, u, impPt, shoPt, pos_floor, predTime, Q, R);

    [~, ~, impErrArray(2, 1, :), shoErrArray(2, 1, :)] = runSim(Ts_base, Ts_1ms, method2, P, v, y_true, g, u, impPt, shoPt, pos_floor, predTime, Q, R);
    [~, ~, impErrArray(2, 2, :), shoErrArray(2, 2, :)] = runSim(Ts_base, Ts_5ms, method2, P, v, y_true, g, u, impPt, shoPt, pos_floor, predTime, Q, R);
    [~, ~, impErrArray(2, 3, :), shoErrArray(2, 3, :)] = runSim(Ts_base, Ts_10ms, method2, P, v, y_true, g, u, impPt, shoPt, pos_floor, predTime, Q, R);
    [~, ~, impErrArray(2, 4, :), shoErrArray(2, 4, :)] = runSim(Ts_base, Ts_15ms, method2, P, v, y_true, g, u, impPt, shoPt, pos_floor, predTime, Q, R);
    [~, ~, impErrArray(2, 5, :), shoErrArray(2, 5, :)] = runSim(Ts_base, Ts_20ms, method2, P, v, y_true, g, u, impPt, shoPt, pos_floor, predTime, Q, R);
    [~, ~, impErrArray(2, 6, :), shoErrArray(2, 6, :)] = runSim(Ts_base, Ts_25ms, method2, P, v, y_true, g, u, impPt, shoPt, pos_floor, predTime, Q, R);
    [~, ~, impErrArray(2, 7, :), shoErrArray(2, 7, :)] = runSim(Ts_base, Ts_30ms, method2, P, v, y_true, g, u, impPt, shoPt, pos_floor, predTime, Q, R);

    [~, ~, impErrArray(3, 1, :), shoErrArray(3, 1, :)] = runSim(Ts_base, Ts_1ms, method3, P, v, y_true, g, u, impPt, shoPt, pos_floor, predTime, Q, R);
    [~, ~, impErrArray(3, 2, :), shoErrArray(3, 2, :)] = runSim(Ts_base, Ts_5ms, method3, P, v, y_true, g, u, impPt, shoPt, pos_floor, predTime, Q, R);
    [~, ~, impErrArray(3, 3, :), shoErrArray(3, 3, :)] = runSim(Ts_base, Ts_10ms, method3, P, v, y_true, g, u, impPt, shoPt, pos_floor, predTime, Q, R);
    [~, ~, impErrArray(3, 4, :), shoErrArray(3, 4, :)] = runSim(Ts_base, Ts_15ms, method3, P, v, y_true, g, u, impPt, shoPt, pos_floor, predTime, Q, R);
    [~, ~, impErrArray(3, 5, :), shoErrArray(3, 5, :)] = runSim(Ts_base, Ts_20ms, method3, P, v, y_true, g, u, impPt, shoPt, pos_floor, predTime, Q, R);
    [~, ~, impErrArray(3, 6, :), shoErrArray(3, 6, :)] = runSim(Ts_base, Ts_25ms, method3, P, v, y_true, g, u, impPt, shoPt, pos_floor, predTime, Q, R);
    [~, ~, impErrArray(3, 7, :), shoErrArray(3, 7, :)] = runSim(Ts_base, Ts_30ms, method3, P, v, y_true, g, u, impPt, shoPt, pos_floor, predTime, Q, R);
end