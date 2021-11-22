close all
clear all
clc

%% INPUT DATA
I_x = 0.07 ;
I_y = 0.0504 ;
I_z = 0.0109 ;
w_0x = 0.45 ;
w_0y = 0.52 ;
w_0z = 0.55 ;

I = diag( [I_x I_y I_z ] ) ;
I_inv = inv( I ) ;
w_0 = [ w_0x w_0y w_0z ] ;

EOM = sim('EquationOfM.slx');
omega = EOM.w;

%% Euler angles






%% Pointing error of a stabilized S/C

I_x = 0.07 ;
I_y = 0.0504 ;
I_z = 0.0109 ;

C = 0.2; % [0.2;2*pi]
w_0x = C ;
w_0y = 0.1 ;
w_0z = 0.1 ;

I = diag( [I_x I_y I_z ] ) ;
I_inv = inv( I ) ;
w_0 = [ w_0x w_0y w_0z ] ;







