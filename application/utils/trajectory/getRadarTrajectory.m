% pegar sinal de trajetoria do radar
function [y, impPt, shoPt, T] = getRadarTrajectory(filename)
    % dados do radar
    data = readmatrix(filename);

    % pontos de impacto e disparo
    impPt = readmatrix(filename, 'Range', 'B12:D12')';
    shoPt = [0; 0; 2.5];

    % vetor tempo t das amostras do radar
    latency = data(:,4);
    t_top_str = readcell(filename, 'Range', 'A10:A10');
    t_0_str = readcell(filename, 'Range', 'B16:B16');
    t_top_str = t_top_str{1}(6:end);
    t_0_str = t_0_str{1};
    deltaT = (datenum(t_0_str, 'HH:MM:SS:FFF') - datenum(t_top_str, 'HH:MM:SS.FFF')).*(24*60*60);
    t = data(:, 3) - data(1, 3) - latency + deltaT;

    % sinal de trajetoria do radar
    y_esf_radar = data(:, 5:7);
    y_esf_radar(:, 1:2) = y_esf_radar(:, 1:2)*180/pi;
    % corrigir amostragem do sinal de trajetoria
    T = 0.05;
    y_esf = resample(y_esf_radar, t, 1/T)';
    % passar sinal de trajetoria para coordenadas cartesianas
    y = zeros(size(y_esf));
    [y(1, :), y(2, :), y(3, :)] = aer2enu(y_esf(2, :), y_esf(1, :), y_esf(3, :));
end