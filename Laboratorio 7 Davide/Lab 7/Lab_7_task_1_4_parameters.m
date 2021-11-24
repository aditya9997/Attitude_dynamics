clc;
close all;
clear all;

%NB: tutte le unità di misura sono quelle del sistema internazionale
% GENERAL ELLIPTIC ORBIT WITH ONLY GRAVITY GRADIENT


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

% 
% out = sim('Lab_7_task_1_4');
% out2 = sim('Lab_7_task_1_4');
% 
% n = out2.simout2;







