clear; clc; close all;

addpath('classes');
addpath(genpath('utils'));

[y, impPt, shoPt, Ts_rad] = getRadarTrajectory("radarData/81mm1.xlsx");

method = 1;
P = 1e8*eye(6);
g = 9.81;
u = [0; 0; -g];
p_floor = 0;
predTime = [1 5];
Q = cov(u');
sigma2_n = [1e2; 1e2; 1e2];
R = [sigma2_n(1) 0 0;
     0 sigma2_n(2) 0;
     0 0 sigma2_n(3)];
[~, impPterr, shoErr] = runSim(Ts_rad, method, P, y, g, u, impPt, shoPt, p_floor, predTime, Q, R);

plot3(y(1, :), y(2, :), y(3, :), 'k.')
hold on
plot3(shoPt(1), shoPt(2), shoPt(3), 'g*')
grid on