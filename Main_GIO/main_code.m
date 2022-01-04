clear all
close all
clc

%% Orbit data

% Constants
G = 6.67259 * 1e-11;   %[m^3/kg/s^2]    Universal gravitational constant
Mt = 5.9726 * 1e24;    %[kg]            Earth's mass
mu_P = G * Mt;         %[m^3/s^2]       Earth gravitational constant

% Orbit initial data
a = 8.6251e06;         %[m]             Semi-major axis
e = 0.1195;            %[-]             Eccentricity
i = 0.3815;            %[rad]           Inclination

% Orbit parametres computation
n = sqrt(mu_P / a^3);  %[rad/s]         Mean angular velocity
T = 2*pi / n;          %[s]             Orbit period


%% Initial state

% True anomaly
theta_0 = 0;                    %[rad]          Initial true anomaly

% Angular velocity
w_0x = 0.12e-5;                       %[rad/s]        omega x
w_0y = 0.18e-5;                       %[rad/s]        omega y
w_0z = 0.15e-5;                       %[rad/s]        omega z
omega_0 = [w_0x; w_0y; w_0z];   %[rad/s] [3x1]  omega vector

% Attitude
A_BN_0 = eye(3);                %[-]            DCM

%% Spacecraft characteristics

MB = [1; 0.3; 0.3; 0.3];       %[kg; m; m; m] Mass, a, b, h
SP = [3; 0.3; 1; 0.1];         %[kg; m; m; m] Mass, a, b, h


%% Environment

% Gravity gradient
[gn, gm, gvali, gsvi] = textread('igrfSg.txt','%f %f %f %f');
[hn, hm, hvali, hsvi] = textread('igrfSh.txt','%f %f %f %f');
G = [gn, gm, gvali, gsvi];
H = [hn, hm, hvali, hsvi];

% Solar radiation
c = 2.99792*10^8;               %[m/s]      Light velocity
rho_d_MB = 0.1;
rho_s_MB = 0.5;
rho_s_SP = 0.1;
rho_d_SP = 0.1;

%% Sensors

%magnetometer
p_acc = 0.5; %[degree]
f_mag = 18; %[Hz] 

%Gyro
t_sample = 1;
sig_b = 0.3*pi/180/3600/sqrt(t_sample);
sig_n = 0.15*pi/180/3600/sqrt(t_sample);
sig_b = sig_b^2;
mr = 1.2*1e-4;

%% Reaction Wheels (using 4 RW with the 3 axis + diagonal model)

A_rw = [1 0 0 1/sqrt(3);...
        0 1 0 1/sqrt(3);
        0 0 1 1/sqrt(3)];

A_rw_star =  [5/6 -1/6 -1/6;...
             -1/6  5/6 -1/6;
             -1/6 -1/6  5/6;
              1/(2*sqrt(3)) 1/(2*sqrt(3)) 1/(2*sqrt(3))];

%% Simulation data

Time = 1*T;         %[s]

%% plot
out = sim('Main_simulation');
    %Slew
      
        %error
figure(1)
plot(out.error)
title("Attitude Error")
xlabel("Time [s]")
ylabel("Attitude Error [-]")
legend("error")
grid on

        %omega
figure(2)
plot(out.omega)
title("ω_BN")
xlabel("Time [s]")
ylabel("ω [rad]")
legend("ω_1", "ω_2", "ω_3")
grid on

    %A_BN
figure(3)
plot(out.A_BN)
title("A_BN")
xlabel("Time [s]")
ylabel("A_BN [rad]")
legend("A11", "A12", "A13", "A21", "A22", "A23", "A31", "A32", "A33")
grid on



