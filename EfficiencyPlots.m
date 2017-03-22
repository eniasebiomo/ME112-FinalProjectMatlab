%% CRAWLER Motor Characterization
clear all
close all
%% a. Stall Torque - plot voltage vs stall current

voltage = 3; %volts
I_stall = 2.28; %amps

R = voltage/I_stall;

%% b. No Load - plot of net Voltage (V - i_nL*R) versus omega_nL
%Convert to radians/sec

% No-Load	
V = 3;
I_nl = 0.17;
RPM = 10380;

%Convert to radians/sec
Rads = RPM .* ((2*pi)/60);

Vnet = V - (I_nl * R);

% coeffs = polyfit(Rads, Vnet, 1);

k = Vnet/Rads;

%We chose to use the Tf with higher voltage because we expect to run our
%motor also at a higher voltage
Tf = k*I_nl;

%% Motor-Winch Test Losses and Efficiency 


m = 0.080; %mass of bolt [kg]
g = 9.81;

I_winch = [0.28];
RPM = [11995];

%Convert to radians/sec
Rads = RPM .* ((2*pi)/60);

height = [0.61]; %distance needed to be raised [m]
time = [22.92]; %time required to raise mass  [s]
velocity = height./time;

Tm = k.*I_winch;

P_in = (Tm - Tf).*Rads;

Fg = m*g;
P_out = Fg*velocity;

Loss = P_in - P_out;

%Transmission Efficiency
Efficiency = P_out./P_in;

%% Graphs of motor efficiency, output power, and torque vs. speed
i_nl = Tf/k;

V3 = 2.74;

istall3 = V3/R;

wnl3 = (V3 - (i_nl*R))/k;

i3 = linspace(i_nl, istall3, 20);

omega3 = (V3-i3*R)/k;

T_l3 = k*i3 - Tf;

P_out3 = T_l3.*omega3;

figure(1)
plot(omega3, T_l3, 'r')
ylabel('Motor Load Torque (Nm)')
xlabel('Angular Veloity (rad/s)')
title('Motor Load Torque V = 3V')

figure(2)
plot(omega3, P_out3, 'r')
ylabel('Motor Output Power (W)')
xlabel('Angular Veloity (rad/s)')
title('Motor Output Power at V = 3V')

P_in3 = V3*i3;

eta_3 = P_out3./P_in3;

figure(3)
plot(i3, eta_3, 'r')
ylabel('Motor Efficiency')
xlabel('Current (Amps)')
title('Motor Efficiency at V = 3V')

plotfixer
