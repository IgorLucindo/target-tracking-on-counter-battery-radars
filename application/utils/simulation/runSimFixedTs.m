% rodar simulacoes com Ts fixo
function [impErrArray, shoErrArray] = runSimFixedTs(Ts, v, y_true, g, u, impPt, shoPt, pos_floor, predTime, Q, R)
    % metodo
    method1 = 1;
    method2 = 2;
    method3 = 3;
    % erro de covariancia
    P1 = 1e1*eye(6);
    P2 = 1e2*eye(6);
    P3 = 1e4*eye(6);
    P4 = 1e6*eye(6);
    P5 = 1e8*eye(6);

    [~, ~, impErrArray(1, 1, :), shoErrArray(1, 1, :)] = runSim(Ts, Ts, method1, P1, v, y_true, g, u, impPt, shoPt, pos_floor, predTime, Q, R);
    [~, ~, impErrArray(1, 2, :), shoErrArray(1, 2, :)] = runSim(Ts, Ts, method1, P2, v, y_true, g, u, impPt, shoPt, pos_floor, predTime, Q, R);
    [~, ~, impErrArray(1, 3, :), shoErrArray(1, 3, :)] = runSim(Ts, Ts, method1, P3, v, y_true, g, u, impPt, shoPt, pos_floor, predTime, Q, R);
    [~, ~, impErrArray(1, 4, :), shoErrArray(1, 4, :)] = runSim(Ts, Ts, method1, P4, v, y_true, g, u, impPt, shoPt, pos_floor, predTime, Q, R);
    [~, ~, impErrArray(1, 5, :), shoErrArray(1, 5, :)] = runSim(Ts, Ts, method1, P5, v, y_true, g, u, impPt, shoPt, pos_floor, predTime, Q, R);

    [~, ~, impErrArray(2, 1, :), shoErrArray(2, 1, :)] = runSim(Ts, Ts, method2, P1, v, y_true, g, u, impPt, shoPt, pos_floor, predTime, Q, R);
    [~, ~, impErrArray(2, 2, :), shoErrArray(2, 2, :)] = runSim(Ts, Ts, method2, P2, v, y_true, g, u, impPt, shoPt, pos_floor, predTime, Q, R);
    [~, ~, impErrArray(2, 3, :), shoErrArray(2, 3, :)] = runSim(Ts, Ts, method2, P3, v, y_true, g, u, impPt, shoPt, pos_floor, predTime, Q, R);
    [~, ~, impErrArray(2, 4, :), shoErrArray(2, 4, :)] = runSim(Ts, Ts, method2, P4, v, y_true, g, u, impPt, shoPt, pos_floor, predTime, Q, R);
    [~, ~, impErrArray(2, 5, :), shoErrArray(2, 5, :)] = runSim(Ts, Ts, method2, P5, v, y_true, g, u, impPt, shoPt, pos_floor, predTime, Q, R);

    [~, ~, impErrArray(3, 1, :), shoErrArray(3, 1, :)] = runSim(Ts, Ts, method3, P1, v, y_true, g, u, impPt, shoPt, pos_floor, predTime, Q, R);
    [~, ~, impErrArray(3, 2, :), shoErrArray(3, 2, :)] = runSim(Ts, Ts, method3, P2, v, y_true, g, u, impPt, shoPt, pos_floor, predTime, Q, R);
    [~, ~, impErrArray(3, 3, :), shoErrArray(3, 3, :)] = runSim(Ts, Ts, method3, P3, v, y_true, g, u, impPt, shoPt, pos_floor, predTime, Q, R);
    [~, ~, impErrArray(3, 4, :), shoErrArray(3, 4, :)] = runSim(Ts, Ts, method3, P4, v, y_true, g, u, impPt, shoPt, pos_floor, predTime, Q, R);
    [~, ~, impErrArray(3, 5, :), shoErrArray(3, 5, :)] = runSim(Ts, Ts, method3, P5, v, y_true, g, u, impPt, shoPt, pos_floor, predTime, Q, R);
end