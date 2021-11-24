clc;
close all;
clear all;


% Inertia moments
I_x = 0.04; %[kg*m^2]
I_y = 0.06; %[kg*m^2]
I_z = 0.08; %[kg*m^2]


%definisco orbita in coordinate cartesiane
r=[-5799.7381; 4578.0767;  2953.8806]*10^3;
v=[-5.085; -5.301; 0.0391]*10^3;
muP = astroConstants(13)*10^9;
G = 6.67259e-11; %[Nm^2/kg^2]
m_t = 5.9726*10^24; %[kg]

%trasformo in coordinate perifocali
[a, e, i, OM, om, th] = car2kep(r, v, muP);
% a = 8.625e6;
% e=0.1195;
% i = 0.3815;

e = norm(e);

th_0 = 0;
h = cross(r,v);
h = norm(h);
n_S = sqrt(muP/a^3);

% Initial conditions of angular velocities
w_x_0 = 1e-6; %[rad/s]
w_y_0 = 1e-6; %[rad/s]
w_z_0 = n_S; %[rad/s]

%Orbital period
T = 3*[2*pi*sqrt(a^3/muP)];

%Defining angular velocity of earth and inclination

eps = deg2rad(11.5); %inclination
w_E = (2*pi)/(24*3600); %angular veloity of earth


%Defining terms for computing H0 of magnetic field (see slides of Lab 7).
%First number refers to the pedice and second to the apice
g_1_0 = -29682*10^-9;
g_1_1 = -1789*10^-9;
h_1_1 = 5318*10^-9;
H_0 = sqrt((g_1_0)^2 + (g_1_1)^2 + (h_1_1)^2);

% Defining the Areas of solar panels (SP) and the main body of the S/C (MB)
A_
