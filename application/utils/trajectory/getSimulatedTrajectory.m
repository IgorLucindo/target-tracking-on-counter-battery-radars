% criar trajetoria simulada
function y = getSimulatedTrajectory(y_true, sigma2_n, predTime, T_ref, detecThresh)
    % ruido de medicao
    n = createNoise(sigma2_n, predTime(2), T_ref);

    arrayLength = ceil(predTime(2)/T_ref);
    % medicao
    y = y_true(:, detecThresh + 1:detecThresh+arrayLength) + n(:, 1:arrayLength);
end