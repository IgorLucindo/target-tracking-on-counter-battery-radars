% previsao do ponto de impacto e disparo
function [impPtPred, shoPtPred] = setImpactShootingPoint(x_est, g, p_floor)
    p0x = x_est(1); p0y = x_est(2); p0z = x_est(3);
    v0x = x_est(4); v0y = x_est(5); v0z = x_est(6);

    ti = (v0z + sqrt(v0z^2 + 2*g*(p0z - p_floor)))/g;
    td = (v0z - sqrt(v0z^2 + 2*g*(p0z - p_floor)))/g;

    p1x = p0x + v0x*ti; p2x = p0x + v0x*td;
    p1y = p0y + v0y*ti; p2y = p0y + v0y*td;

    impPtPred = [p1x; p1y; p_floor];
    shoPtPred = [p2x; p2y; p_floor];
end