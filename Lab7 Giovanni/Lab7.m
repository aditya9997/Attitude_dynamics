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
Rs_MB = ones(1,6)*0.5;
Rs_SP = ones(1,4)*0.1;
Rd = ones(1,10)*0.1;
A_MB = ones(1,6)*0.5; %[m^2]
A_SP = ones(1,4)*1; %[m^2]
c = ones(1,10)*299792*10^3; %[m/s]

n_Sun = 2*pi/(365*24*60*60);
eps = degtorad(23.45);

%vettore inerzia [kg*m^2]
% Ixx = 3.3311;
% Ixy = -0.0119;
% Ixz = 0.0027;
% Iyx = Ixy;
% Iyy = 2.2024;
% Iyz = 0.0161;
% Izy = Iyz;
% Izx = Ixz;
% Izz = 1.8033;
% 
% Imat = [Ixx, Ixy, Ixz;
%         Iyx, Iyy, Iyz;
%         Izx, Izy, Izz];
%     
% I = eig(Imat);
% Iinv = inv(diag(I));
% 
% Ix = I(1);
% Iy = I(2);
% Iz = I(3); 

Ix = 0.04;
Iy = 0.06;
Iz = 0.08; 

%condizioni iniziali
om_x = 1e-6;
om_y = 1e-6;
om_z = n;
omega_0 = [om_x, om_y, om_z];

%tempo
 
Time = 5*(2*pi/n);

% out = sim('Quaternions_Method');
% om_BL = out.Omega_BL.data;
% s = size(om_BL, 1);
% t=linspace(0, Time, s);
% 
% 
% plot(t, om_BL(:,:))


