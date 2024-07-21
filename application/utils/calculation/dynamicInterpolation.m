% calcula interpolacao dynamicamente de um periodo de amostragem para outro a partir da estimacao de kalman filter
function [y_new, kf] = dynamicInterpolation(y, y_new, Ts, Ts_new, kf, u, i)
    rescale = Ts/Ts_new;
    if mod(i - 1, rescale) == 0 && Ts ~= Ts_new
        j = (i - 1)/rescale + 1;

        % kalman filter
        if j == 1
            [kf, ~, ~] = kf.setX0(y(:, 1));
            return
        end
        [kf, y_est, ~] = kf.run(y(:, j), u);

        % interpolacao
        y_new(1, (j-2)*rescale + 1:j*rescale) = interp1(0:2, [y(1, j-1:j) y_est(1)], 0:1/rescale:2 - 1/rescale, 'linear');
        y_new(2, (j-2)*rescale + 1:j*rescale) = interp1(0:2, [y(2, j-1:j) y_est(2)], 0:1/rescale:2 - 1/rescale, 'linear');
        y_new(3, (j-2)*rescale + 1:j*rescale) = interp1(0:2, [y(3, j-1:j) y_est(3)], 0:1/rescale:2 - 1/rescale, 'spline');
    end
end