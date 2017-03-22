% ME 112 - Eni Asebiomo 
% Final Project
clear all
close all

run BirdFootLinkageAnalysis.m

% Foot Velocity from Linkage Analysis
vx = x_vel_ms;
vy = y_vel_ms;

%Read in Force Plate Data from excel file
filename = 'CuratedForcePlateData.xlsx';

%Forces (Normal, Shear)
force_plate = xlsread(filename, 'A2:C601');
time = force_plate(:,1);
fnormal = force_plate(:,2);
fshear = force_plate(:,3);

fnavg = mean(fnormal(54:569)); % Trim data to only include actual step

m = 0.528; % Bird Mass in kg
g = 9.81;
fncalc = m*g;

alpha = fncalc/fnavg;
fy = fnormal .* alpha; % Adjusted normal and shear forces to 
fx = fshear .* alpha;  % account for variability in data collection

figure
plot(time, fnormal, time, fshear);
legend('Normal Force','Shear Force');
title('Force Plate Data');
xlabel('time (ms)');
ylabel('force (N)');
figure
plot(time, vx, time, vy);
legend( 'X Velocity', 'Y Velocity');

fmag = sqrt(fnormal.^2 + fshear.^2);
fdir = atan(fnormal./fshear);

fx = fx'; % Transposing matricies in order to calculate dot product
fy = fy';

power = abs(fx.*vx) + abs(fy.*vy);
avgpower = mean(power(54:569));

figure
plot(time, power);
title('Output power during stride');
xlabel('time (ms)');
ylabel('power (W)');

avgvx = mean(vx);
mu = 0.8;
Pcalc = mu*m*g*avgvx;

plotfixer

