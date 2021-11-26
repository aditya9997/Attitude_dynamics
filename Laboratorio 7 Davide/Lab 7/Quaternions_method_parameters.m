clear all
close all
clc

%% definisco dati

%costanti
   %GG
G = 6.67259 * 1e-11; %[N*m^2*kg^-2]
Mt = 5.9726 * 1e24; %[kg]
muP = G*Mt;

%definisco orbita in coordinate cartesiane
r=[-5799.7381; 4578.0767;  2953.8806]*10^3;
v=[-5.085; -5.301; 0.0391]*10^3;

%trasformo in coordinate perifocali
[a, e, i, OM, om, th] = car2kep(r, v, muP);
e=norm(e);
ii = i;
theta00 = 0;

a = 8.6251e+06;
e = 0.1195;
ii = 0.3815;

n = sqrt(muP/a^3);  

%SRP
RS_MB = 0.5;
RS_SP = 0.1;
RD_MB = 0.1*ones(6,1);
RD_SP = 1;
A_MB = 0.5; %[m^2]
A_SP = 1; %[m^2]
c = 299792*10^3; %[m/s]
n_Sun = 2*pi/(365*24*60*60);
eps = deg2rad(23.45);
w_E = (2*pi)/(24*3600); %angular veloity of earth

%Defining terms for computing H0 of magnetic field (see slides of Lab 7).
%First number refers to the pedice and second to the apice
g_1_0 = -29682*10^-9;
g_1_1 = -1789*10^-9;
h_1_1 = 5318*10^-9;
H_0 = sqrt((g_1_0)^2 + (g_1_1)^2 + (h_1_1)^2);

%vettore inerzia [kg*m^2]
Ix = 0.04;
Iy = 0.06;
Iz = 0.08; 

%condizioni iniziali
om_x = 1e-6;
om_y = 1e-6;
om_z = n;
omega_0 = [om_x, om_y, om_z];

%tempo
Time = (2*pi/n);