% calcula interpolacao dynamicamente de um periodo de amostragem para outro a partir da estimacao de kalman filter
function y_new = dynamicInterp(y, y_new, T, T_new, i)
    if T == T_new
        return
    end
    % valor inicial
    if i == 3
        y_new(:, 1) = y(:, 1);
    end

    rescale = T/T_new;
    % interpolacao de 3 em 3 pontos
    y_new(1, (i-3)*rescale + 2:(i-1)*rescale + 1) = interp1(0:2, [y(1, i-2:i)], 1/rescale:1/rescale:2, 'linear');
    y_new(2, (i-3)*rescale + 2:(i-1)*rescale + 1) = interp1(0:2, [y(2, i-2:i)], 1/rescale:1/rescale:2, 'linear');
    y_new(3, (i-3)*rescale + 2:(i-1)*rescale + 1) = interp1(0:2, [y(3, i-2:i)], 1/rescale:1/rescale:2, 'spline');
end