% pegar sinal de trajetoria do radar dinamicamente
function [y, impPt, shoPt, T] = getRadarTrajectoryDynamically(i)
    % arquivo
    radarFiles = dir('radarData');
    radarFile = radarFiles(i + 2);

    % pegar trajetoria
    [y, impPt, shoPt, T] = getRadarTrajectory("radarData/" + radarFile.name);
end