close all
clear all
clc

%% Relative attitide error
I_x = 0.04 ;
I_y = 0.06 ;
I_z = 0.08 ;
w_0x = 0 ;
w_0y = 0 ;
w_0z = 2*pi ;

I = diag( [I_x I_y I_z ] ) ;
I_inv = inv( I ) ;
w_0 = [ w_0x w_0y w_0z ] ;

% Orbit characheristics (LVLH frame)
R = 100000 ;
mu= 398600 ;
n = sqrt(mu/R^3);
wl_0x = 1e-6 ;
wl_0y = 1e-6 ;
wl_0z = n ;
omega_L = [wl_0x; wl_0y; wl_0z] ;

M = [0 0 0]';

SimTime = 2000;
EOM = sim('EquationOfM.slx');
omega = EOM.w;
A0 = eye(3);
DCM_N_to_B = sim('DCM_N_to_B.slx');
DCM = DCM_N_to_B.DCM ;

LVLH = sim('LVLH_error.slx');

figure
plot(LVLH.attitude_err)

%% Gravity gradient
I_x = 0.04 ;
I_y = 0.06 ;
I_z = 0.08 ;
w_0x = 0 ;
w_0y = 0 ;
w_0z = 2*pi ;

I = diag( [I_x I_y I_z ] ) ;
I_inv = inv( I ) ;
w_0 = [ w_0x w_0y w_0z ] ;

R = 100000 ;
mu = 398600 ;
n = sqrt(mu/R^3);
wl_0x = 1e-6 ;
wl_0y = 1e-6 ;
wl_0z = n ;
omega_L = [wl_0x; wl_0y; wl_0z] ;

SimTime = 2000;

% LVLH = sim('LVLH_error.slx');
% A_BL = LVLH.A_BL;
A_BL = eye(3);

GG = sim('Gravity_gradient.slx');
M = Gravity_gradient.M;
EOM = sim('EquationOfM.slx');
omega = EOM.w;
