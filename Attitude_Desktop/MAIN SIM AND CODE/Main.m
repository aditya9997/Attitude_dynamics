clear all
close all
clc

%re = 6371*1e3;
%a = 7195 *1e3; %semi major axis in meters
%i = deg2rad(98.7);
%rp = 818.3*1e3 + re;
%ra = 830.8*1e3 + re;
%e = (ra-rp)/(ra+rp);
%% definisco dati

%Orbit definition
G = 6.67259 * 1e-11; %[N*m^2*kg^-2]
Mt = 5.9726 * 1e24; %[kg]
muP = G*Mt;
a = 8.6251e06;
e = 0.1195;
i = 0.3815;
n = sqrt(muP/a^3);  
T = 2*pi/n;


%initial condition A_BN
A_BN_0 =  [1/3, 2/3, 2/3; ...
           2/3, 1/3, -2/3; ...
          -2/3, 2/3,-1/3];

%SRP
rho_d_MB = 0.1;
rho_s_MB = 0.5;
rho_s_SP = 0.1;
rho_d_SP = 0.1;

c = 299792*10^3; %[m/s]
n_Sun = 2*pi/(365*24*60*60);
eps = deg2rad(23.45);
w_E = (2*pi)/(24*3600); %angular veloity of earth

%condizioni iniziali
om_x = 2;
om_y = 2;
om_z = 2;
omega_0 = [om_x, om_y, om_z];
theta0 = 0;

%tempo
Time = 1*(2*pi/n);

%Mass and three dimensions of the S/C main body(MB) and solar panel(SP)
MB = [10; 0.3; 0.3; 0.3];
SP = [30; 0.3; 1; 0.1];

%magnetometer
p_acc = 0.5; %[degree], look at slide 12 of attitude sensors, bernelli's slide
f_mag = 18; %[Hz] by https://www.cubesatshop.com/product/nss-magnetometer/, in alternative we can put f=5 since bernelli's slides

% Reaction Wheels (using 4 RW with the 3 axis + diagonal model)
A_rw = [1 0 0 1/sqrt(3);...
        0 1 0 1/sqrt(3);
        0 0 1 1/sqrt(3)];

A_rw_star =  [5/6 -1/6 -1/6;...
             -1/6  5/6 -1/6;
             -1/6 -1/6  5/6;
              1/(2*sqrt(3)) 1/(2*sqrt(3)) 1/(2*sqrt(3))];


% %Gyro
t_sample = 0.1;
sig_b = 0.3*pi/180/3600/sqrt(t_sample);
sig_n = 0.15*pi/180/3600/sqrt(t_sample);
sig_b = sig_b^2;

mr = 1.2*1e-4;
p_acc = deg2rad(5);
f_mag = 5;


rho_d_MB = 0.1;
rho_s_MB = 0.5;

[gn, gm, gvali, gsvi] = textread('igrfSg.txt','%f %f %f %f');
[hn, hm, hvali, hsvi] = textread('igrfSh.txt','%f %f %f %f');
G = [gn, gm, gvali, gsvi];
H = [hn, hm, hvali, hsvi];

%slew manoevre A_BN desired
A_slew = [0.4777 0.0636 0.8762;
    0.7262 0.5327 0.4345;
    0.4944 0.8439 0.2083];