clear all
close all
clc

%% Orbit data

% Constants
G = 6.67259 * 1e-11;   %[m^3/kg/s^2]    Universal gravitational constant
Mt = 5.9726 * 1e24;    %[kg]            Earth's mass
mu_P = G * Mt;         %[m^3/s^2]       Earth gravitational constant
c = 2.99792*10^8;      %[m/s]           Light velocity

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
w_0x = 0.12e-5;                     %[rad/s]        omega x
w_0y = 0.18e-5;                     %[rad/s]        omega y
w_0z = 0.15e-5;                     %[rad/s]        omega z
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

t_sample = 1;

% Sun sensor
acc_S = deg2rad(0.25);

% Magnetometre
acc_m = deg2rad(5);
mr = 1.2*1e-4;

% Gyroscopes
sig_b = 0.3*pi/180/3600/sqrt(t_sample);
sig_n = 0.15*pi/180/3600/sqrt(t_sample);
sig_b = sig_b^2;

%% Reaction Wheels (using 4 RW with the 3 axis + diagonal model)

A_rw = [1 0 0 1/sqrt(3);...
        0 1 0 1/sqrt(3);
        0 0 1 1/sqrt(3)];

A_rw_star =  [5/6 -1/6 -1/6;...
             -1/6  5/6 -1/6;
             -1/6 -1/6  5/6;
              1/(2*sqrt(3)) 1/(2*sqrt(3)) 1/(2*sqrt(3))];

%% Simulation data

Time = 1*T;                     %[s]
control = 2;                    %[-] Set control to choose the maneouvre, 1 for detumbling, 2 for trajectory tracking, 3 for slew

%% Maneouvres

if control == 1
    
    %Detumbling

    % Angular velocity
    w_0x = 1;                       %[rad/s]        omega x
    w_0y = 2;                       %[rad/s]        omega y
    w_0z = 3;                       %[rad/s]        omega z 
    omega_0 = [w_0x; w_0y; w_0z];   %[rad/s] [3x1]  omega vector

    % Attitude
    A_BN_0 = eye(3);                %[-]            DCM
    
    out = sim('Main_simulation_davide');

    %omega
    figure(1)
    subplot(3,1,1)
    plot(out.w_B1)
    title("ω_Bx")
    xlabel("Time [s]")
    ylabel("ω [rad]")
    legend("ω_1")
    grid on
    subplot(3,1,2)
    plot(out.w_B2)
    title("ω_By")
    xlabel("Time [s]")
    ylabel("ω [rad]")
    legend("ω_2")
    grid on    
    subplot(3,1,3)
    plot(out.w_B3)
    title("ω_Bz")
    xlabel("Time [s]")
    ylabel("ω [rad]")
    legend("ω_3")
    grid on

elseif control == 3
    %Slew

    % Angular velocity
    w_0x = 0.7;                       %[rad/s]        omega x
    w_0y = 0.5;                       %[rad/s]        omega y
    w_0z = 0.3;                       %[rad/s]        omega z
    omega_0 = [w_0x; w_0y; w_0z];     %[rad/s] [3x1]  omega vector

    % Attitude
    A_BN_0 = eye(3);                %[-]            DCM
    %A_BN_0 =  [1/3, 2/3, 2/3; 2/3, 1/3, -2/3; -2/3, 2/3,-1/3];
   
    out = sim('Main_simulation_davide');
    

    %error
    figure(1)
    plot(out.contr_err)
    title("Attitude Error")
    xlabel("Time [s]")
    ylabel("Attitude Error [-]")
    legend("error")
    grid on
    
    %omega
    figure(2)
    plot(out.w_B)
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

elseif control == 2
    
    % Trajectory tracking

    % Angular velocity
    w_0x = 0.12e-5;                       %[rad/s]        omega x
    w_0y = 0.18e-5;                       %[rad/s]        omega y
    w_0z = 0.15e-5;                       %[rad/s]        omega z
    omega_0 = [w_0x; w_0y; w_0z];         %[rad/s] [3x1]  omega vector

    % Attitude
    A_BN_0 = eye(3);                %[-]            DCM
    %A_BN_0 =  [1/3, 2/3, 2/3; 2/3, 1/3, -2/3; -2/3, 2/3,-1/3];
   
    out = sim('Main_simulation_davide');
    out1 = sim('Main_simulation_davide.slx');

    %error7
    figure(1)
    plot(out.contr_err)
    title("Attitude Error")
    xlabel("Time [s]")
    ylabel("Attitude Error [-]")
    legend("error")
    grid on
    
    %omega
    figure(2)
    plot(out.w_B)
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

end
%% Plot
figure(1);
subplot(3, 1, 1)
plot(out1.w_B1, 'r', 'Linewidth', 2);
hold on
plot(out.w_B1, 'b', 'Linewidth', 2);
xlabel('Time [s]');
ylabel('w_x [rad/s]');
legend( 'Closed loop', 'Open loop');
grid on
hold off

subplot(3, 1, 2)
plot(out1.w_B2, 'r', 'Linewidth', 2);
hold on 
plot(out.w_B2, 'b', 'Linewidth', 2);
xlabel('Time [s]');
ylabel('w_y [rad/s]');
legend( 'Closed loop', 'Open loop');
grid on
hold off

subplot(3, 1, 3)
plot(out1.w_B3, 'r', 'Linewidth', 2);
hold on 
plot(out.w_B3, 'b', 'Linewidth', 2);
xlabel('Time [s]');
ylabel('w_z [rad/s]');
legend( 'Closed loop', 'Open loop');
grid on
hold off

figure(2);
subplot(3, 1, 1)
plot(out1.err1, 'r', 'Linewidth', 2);
hold on 
plot(out.err1, 'b', 'Linewidth', 2);
xlabel('Time [s]');
ylabel('Err [deg] ');
legend( 'Closed loop', 'Open loop');
grid on
hold off

subplot(3, 1, 2)
plot(out1.err2, 'r', 'Linewidth', 2);
hold on 
plot(out.err2, 'b', 'Linewidth', 2);
xlabel('Time [s]');
ylabel('Err2 [deg] ');
legend( 'Closed loop', 'Open loop');
grid on
hold off

subplot(3, 1, 3)
plot(out1.err3, 'r', 'Linewidth', 2);
hold on 
plot(out.err3, 'b', 'Linewidth', 2);
xlabel('Time [s]');
ylabel('Err3 [deg] ');
legend( 'Closed loop', 'Open loop');
grid on
hold off