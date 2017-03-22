%ME 112 - Andrew Edoimioya
%Final Project
%Deflection Calculations


%% Change these to fit your settings

E = 3.2e9; %Young's Modulus 

% Hardboard (duron): 3.1-5.52 GPa [Pa = N/m^2]
% Acrylic: 3.2 GPa
% PLA plastic: 3.5 GPa

b = 1.5; %[cm]
h = 0.59; %[cm]

I = (b*h^3)/12; %Moment of inertia [cm^4)

P = 2.5; %Force [N]
L = 10; %Length of link [m]

%% Maximun deflection

nu_max = (-P*L^3)/(3*E*I); %[cm]
