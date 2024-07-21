% criar sinal de trajetoria de radar
function y = createTrajectory(y_true, sigma2_n, predTime, Ts_ref)
    % ruido de medicao
    n = createNoise(sigma2_n, predTime(2), Ts_ref);

    % amostra em que comeca a medida referente ao tempo de 10 segundos
    detecThresh = ceil(10/Ts_ref);
    arrayLength = ceil(predTime(2)/Ts_ref);

    % medicao
    y = y_true(:, detecThresh:detecThresh+arrayLength) + n(:, 1:1+arrayLength);
end