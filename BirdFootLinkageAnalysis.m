close all;
clear all;
% Andrew Edoimioya, 2/27/2017
% Adapted from HalfAngleMethod by Mark Cutkosky

% This method uses tangent half-angle formulas to simplify catching
% the 2 inversions (see Wikipedia "tangent half angle")
% Assumes link 1 is aligned with X axis.
%See '4bar_assemblies_diagram.pdf' in this folder for explanation
%of notation

% Start with equations (9.29) - (9.31) per Stanisic.
% Vector loop is R2 + R3 = R1 + R4, where R1 is horizontal
% from left to right. All angles anticlockwise from horizontal.
% Theta2 is input crank angle.

%Define link lengths and inversion
%(Alternatively one could read in initial joint xy locations
%from a file and use those to get the link lengths)
R1 = 3;
R2 = 1.125;
R3 = 3.75;
R4 = 2.5;
inversion = 1;    % 1 if inverted solution, else 0

%Add provisions for a coupler location
%(Modify this to suit your design.)
gammac = -1*pi*(50/180); %units of radians
r5=2.92;


%Define range of angles to plot and step size
theta2start = 0;
theta2end = 2*pi;
numsteps = 10;
theta2step = (theta2end-theta2start)/numsteps;
theta2s = [theta2start:theta2step:theta2end];

% Empirically determined crank angular velocity
omega_L = 2.18; %[rad/s]

%Initialize plot
clf;
%Guess a reasonable plot area (may need to tweak these)
leftlim = -2;
rightlim = 6;
bottomlim = -4;
toplim = 4;
axis([leftlim rightlim bottomlim toplim],'equal');
grid('on')

%The fixed ground link
X1 = 0;
Y1 = 0;
X4 = R1;
Y4 = 0;

hold on;
%Ground link rotation angle and rotation matrix
theta_gnd = pi*(35/180);
rot_mat = [cos(theta_gnd) -sin(theta_gnd) 0; sin(theta_gnd) cos(theta_gnd) 0;...
    0 0 1];

homo_gnd = [X1 Y1 1; X4 Y4 1]';
tran_gnd = rot_mat*homo_gnd;

%Plot rotated ground link X1             Y4             Y1
line([tran_gnd(1,2);tran_gnd(1,1)],[tran_gnd(2,2);tran_gnd(2,1)],'Color','k','LineWidth',1);

hold on;

for i=1:length(theta2s)

cosq2 = cos(theta2s(i));
sinq2 = sin(theta2s(i));
  
% Define some substitutions (9.29-9.31)
A = 2*R4*(R1 - R2*cosq2);
B = -2*R2*R4*sinq2;
C = R3^2 - R2^2 - R4^2 - R1^2 + 2*R2*R1*cosq2;

%Solve for the roots of half angle equation (9.35)
%If general
u41 = (B + sqrt(A^2 + B^2 - C^2))/(A+C);
q4 = 2*atan(u41);

%If inversion
if(inversion)
u42 = (B - sqrt(A^2 + B^2 - C^2))/(A+C);
q4 = 2*atan(u42);
end        %end if

%Solve for theta3 (9.37-9.39)
%range = -PI to +PI
cosq3 = (-R2*cosq2 + R4*cos(q4) + R1)/R3;
sinq3 = (-R2*sinq2 + R4*sin(q4))/R3;
q3 = atan2(sinq3,cosq3);  %don't really need it for plotting

%Compute locations of joints 2,3
X2(i) = R2 * cosq2;
Y2(i) = R2 * sinq2;
X3(i) = X2(i)+R3*cosq3;
Y3(i) = Y2(i)+R3*sinq3;


%Coupler location
xca1(i) = R2*cos(theta2s(i))+r5*cos(q3+gammac);
yca1(i) = R2*sin(theta2s(i))+r5*sin(q3+gammac);
theta_c(i) = q3+gammac;

if(i>1)
   dx = xca1(i) - xca1(i-1);
   dy = yca1(i) - yca1(i-1);
   x_velocities(i) = (dx/theta2step)*omega_L;
   y_velocities(i) = (dy/theta2step)*omega_L;
end

% pause(1);  %Omit if you don't want it to wait for key press at each loop
end    %end for i

%Affine Transformations
homo_2 = [X2; Y2; ones(1, length(X2))];
tran_2 = rot_mat*homo_2;

homo_3 = [X3; Y3; ones(1, length(X3))];
tran_3 = rot_mat*homo_3;

homo_c = [xca1; yca1; ones(1, length(xca1))];
tran_c = rot_mat*homo_c;

for i=1:length(theta2s)
    %Now plot linkages 2 and 3 transformed
    line([tran_gnd(1,1);tran_2(1,i)],[tran_gnd(2,1);tran_2(2,i)],'Color','g','LineWidth',1);
    line([tran_2(1,i);tran_3(1,i)],[tran_2(2,i);tran_3(2,i)],'Color','r','LineWidth',1);
    line([tran_3(1,i);tran_gnd(1,2)],[tran_3(2,i);tran_gnd(2,2)],'Color','b','LineWidth',1);
    
    %Plot Coupler
    plot(tran_c(1,i),tran_c(2,i),'*');
    line([tran_2(1,i);tran_c(1,i)],[tran_2(2,i);tran_c(2,i)],'Color','y','LineWidth',1);
    pause(1);
end

%Now we want to plot the foot onto the coupler curve
%We first want to combine the data of the coupler point and the orientation of
%the coupler foot
couplerdata = [xca1; yca1; theta_c]'; 

polygon = load('PointsData.txt');

%Design 1 (Lex's Design) - Uncomment when using design 1
% RotPoly = Rotq(1.55*pi/2);
% RotGlobal = Rotq(3.25*pi/4);

%Design 2 (Andrew's Design) - Uncomment when using design 2
% RotPoly = Rotq(1.55*pi/2);
% RotGlobal = Rotq(2.75*pi/4);

%Design 3 - Uncomment when using design 3
% RotPoly = Rotq(0.48*pi/2);
% RotGlobal = Rotq(1.2*pi/4);

%Design 4 - Uncomment when using design 3
RotPoly = Rotq(0.45*pi/2);
RotGlobal = Rotq(0.9*pi/4);

homogenous = [polygon ones(length(polygon),1)]';
%Initialize plot
figure(2)
%Guess a reasonable plot area (may need to tweak these)
% Design 1
% leftlim = -6;
% rightlim = 0;
% bottomlim = -6;
% toplim = 0;

% Design 2
% leftlim = -4;
% rightlim = 0;
% bottomlim = -4;
% toplim = 0;

% Design 3
% leftlim = -2;
% rightlim = 4;
% bottomlim = -4;
% toplim = 0;

%Design 4
leftlim = -1;
rightlim = 4;
bottomlim = -4;
toplim = 0;
axis([leftlim rightlim bottomlim toplim],'equal');
grid('on')

for i=1:length(couplerdata)
    x = couplerdata(i, 1);
    y = couplerdata(i, 2);
    theta = couplerdata(i, 3);
    Rot = Rotq(theta);
    Tran = Tranxy(x,y);
    points = RotGlobal*Tran*RotPoly*Rot*homogenous;
    plotlines(points, 'b', 'r');
%     pause(1);
    hold on
end

%Plot velocity information in m/s
x_vel_ms = x_velocities.*0.0254;
y_vel_ms = y_velocities.*0.0254;

vel = sqrt(x_vel_ms.^2 + y_vel_ms.^2);

figure(3)
plot(theta2s, x_vel_ms, theta2s, y_vel_ms, theta2s, vel)
xlabel('Crank angle, \theta [radians]')
ylabel('Foot Velocity [m/s]')
legend('X-velocity', 'Y-velocity', 'Total Speed')
leftlim = 0;
rightlim = 6.28;
bottomlim = -0.15;
toplim = 0.15;
axis([leftlim rightlim bottomlim toplim]);
plotfixer

%Compute average, max and min velocities
avg_vel = mean(vel);
max_vel = max(vel);
min_vel = min(vel(vel>0));