close all
clear all
clc

%% INPUT DATA
I_x = 0.07 ;
I_y = 0.0504 ;
I_z = 0.0109 ;
w_0x = 0.5 ;
w_0y = 0 ;
w_0z = 0 ;

SimTime = 100;
I = diag( [I_x I_y I_z ] ) ;
I_inv = inv( I ) ;
w_0 = [ w_0x w_0y w_0z ] ;

EOM = sim('EquationOfM.slx');
omega = EOM.w;
figure
plot(omega)

%% DCM

SimTime = 100;
%A0 = eye(3);
A0 = [sqrt(2) sqrt(2) 0
      -sqrt(2) sqrt(2) 0
      0 0 1];
DCM = sim('DCM_N_to_B.slx');

%plot
figure
plot(DCM.DCM)
figure
plot(DCM.EA)
% figure
% sc = importGeometry('spacecraft.stl');
% pdegplot(sc); hold on
% quiver3(A0(1,:)); hold on
% quiver3(A0(2,:)); hold on
% quiver3(A0(3,:)); hold on
% plot3(DCM.DCM(1,:), DCM.DCM(2,:), DCM.DCM(2,:)); hold on
% plot3(DCM.DCM(2,1,:), DCM.DCM(2,2,:), DCM.DCM(2,2,:)); hold on
% plot3(DCM.DCM(3,1,:), DCM.DCM(3,2,:), DCM.DCM(3,2,:));

%% quaternions

q0 = [0 0 1 0];
SimTime = 100;
q = sim('Quaternions.slx');

%plot
figure
plot(q.q)
figure
plot(q.DCM)

