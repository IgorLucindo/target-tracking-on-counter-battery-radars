% calcular erro dos pontos de impacto e de disparo
function [impErr, shoErr] = getImpactShootingError(impPtPred, shoPtPred, impPt, shoPt)
    impErr = norm(impPtPred - impPt);
    shoErr = norm(shoPtPred - shoPt);
end