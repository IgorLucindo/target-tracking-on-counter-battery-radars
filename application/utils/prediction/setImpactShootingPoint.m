% previsao do ponto de impacto e disparo
function [impPtPred, shoPtPred] = setImpactShootingPoint(x_est, g, pos_floor)
    pos0x = x_est(1); pos0y = x_est(2); pos0z = x_est(3);
    vel0x = x_est(4); vel0y = x_est(5); vel0z = x_est(6);

    ti = (vel0z + sqrt(vel0z^2 + 2*g*(pos0z - pos_floor)))/g;
    td = (vel0z - sqrt(vel0z^2 + 2*g*(pos0z - pos_floor)))/g;

    pos1x = pos0x + vel0x*ti; pos2x = pos0x + vel0x*td;
    pos1y = pos0y + vel0y*ti; pos2y = pos0y + vel0y*td;

    impPtPred = [pos1x; pos1y; pos_floor];
    shoPtPred = [pos2x; pos2y; pos_floor];
end