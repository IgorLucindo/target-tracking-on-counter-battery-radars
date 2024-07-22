clear; clc; close all;

addpath('classes');
addpath(genpath('utils'));

[y, impPt, shoPt, Ts_rad] = getRadarTrajectory("radarData/81mm1.xlsx");

figure
plot3(y(1, :), y(2, :), y(3, :), 'k.')
hold on
plot3(impPt(1), impPt(2), impPt(3), 'r*')
hold on
plot3(shoPt(1), shoPt(2), shoPt(3), 'g*')
grid on