% pegar sinal de trajetoria do radar dinamicamente
function [y, impPt, shoPt, Ts] = getRadarTrajectoryDynamically(i)
    % arquivo
    radarFiles = dir('radarData');
    radarFile = radarFiles(i + 2);

    % pegar trajetoria
    [y, impPt, shoPt, Ts] = getRadarTrajectory("radarData/" + radarFile.name);
end