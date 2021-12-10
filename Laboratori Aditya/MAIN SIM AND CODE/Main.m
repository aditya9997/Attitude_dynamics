clear all
close all
clc

%% definisco dati

%Orbit definition
G = 6.67259 * 1e-11; %[N*m^2*kg^-2]
Mt = 5.9726 * 1e24; %[kg]
muP = G*Mt;
a = 8.6251e+06;
e = 0.1195;
i = 0;
n = sqrt(muP/a^3);  

%SRP
RS_MB = 0.5;
RS_SP = 0.1;

RD = 0.1*ones(1,10);

A_MB = 0.5*ones(1,6); %[m^2]

RD_MB = 0.1*ones(6,1);
RD_SP = 0.1*ones(4,1);


A_SP = 1; %[m^2]
c = 299792*10^3; %[m/s]
n_Sun = 2*pi/(365*24*60*60);
eps = degtorad(23.45);
w_E = (2*pi)/(24*3600); %angular veloity of earth

%Defining terms for computing H0 of magnetic field (see slides of Lab 7).
%First number refers to the pedice and second to the apice
g_1_0 = -29682*10^-9;
g_1_1 = -1789*10^-9;
h_1_1 = 5318*10^-9;
H_0 = sqrt((g_1_0)^2 + (g_1_1)^2 + (h_1_1)^2);

%condizioni iniziali
om_x = 0;
om_y = 0;
om_z = n;
omega_0 = [om_x, om_y, om_z];
theta0 = 0;

%tempo
Time = 1*(2*pi/n);

%Mass and three dimensions of the S/C main body(MB) and solar panel(SP)
MB = [2000; 0.3; 1; 0.12].*1e-2;
SP = [2000; 0.3; 1; 0.12].*1e-2;

%magnetometer
p_acc = 0.5; %[degree], look at slide 12 of attitude sensors, bernelli's slide
f_mag = 18; %[Hz] by https://www.cubesatshop.com/product/nss-magnetometer/, in alternative we can put f=5 since bernelli's slides
