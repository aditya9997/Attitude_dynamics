clear all
close all
clc

%% EXERCISE 1
%inertia moments
I_x = 0.07 ;
I_y = 0.0504 ;
I_z = 0.0109 ;
I = diag( [I_x I_y I_z ] ) ;

%generic spin axis
w_0x = 0.45 ;
w_0y = 0.52 ;
w_0z = 0.55 ;
% w_0x = 0 ;
% w_0y = 0 ;
% w_0z = 0.55 ;

I_inv = inv( I ) ;
w_0 = [ w_0x w_0y w_0z ] ;

sim('EquationOfM.slx')

% plot
figure
subplot(1,2,1)
plot(ans.w)
ylabel('w [rad/s]')
legend( 'w_x', 'w_y', 'w_z' )
grid on

subplot(1,2,2)
plot(ans.wdot)
ylabel('wdot [rad/s^2]')
legend( 'wdot_x', 'wdot_y', 'wdot_z' )
grid on

%% EXERCISE 2
%simulation of symmetric s/c
I(1,1) = I_y ;
ans1 = sim('EquationOfM.slx') ;

%analytic solution
lambda = ( (I_z - I_y) / I_y ) * w_0z ;
wrealx = @(t) w_0x * cos( lambda*t ) - w_0y * sin( lambda*t ) ; 
wrealy = @(t) w_0x * sin( lambda*t ) + w_0y * cos( lambda*t ) ;
wrealz = @(t) w_0z ;

%plot
figure
subplot(1,2,1)
t = 0:0.1:10 ;
plot(wrealx(t)), hold on
plot(wrealy(t)), hold on
plot(wrealz(t))
ylabel('w real [rad/s]')
legend( 'w_x', 'w_y', 'w_z' )
grid on

subplot(1,2,2)
plot(ans1.w)
ylabel('w [rad/s]')
legend( 'w_x', 'w_y', 'w_z' )
grid on

%% EXERCISE 3
I(1,1) = 0.01 ;
I(2,2) = 0.05 ;
I(3,3) = 0.07 ;
w_0 = [0.1 0.1 0.1];

spinaxis = 3;

switch spinaxis
    case 1
        w_0(1) = 2*pi ;
    case 2
        w_0(2) = 2*pi ;
    case 3
        w_0(3) = 2*pi ;
end


ans2 = sim('EquationOfM.slx') ;
figure
plot(ans2.w)
ylabel('w [rad/s]')
legend( 'w_x', 'w_y', 'w_z' )
grid on

%% EXERCISE 4
I(1,1) = 0.07 ;
I(2,2) = 0.0504 ;
I(3,3) = 0.0109 ;
Iw = 0.005 ;

w_0 = [1e-6 1e-6 0.02] ;

spinaxis = 3;

% switch spinaxis
switch spinaxis
    case 1
        w_0(1) = 2*pi ;
    case 2
        w_0(2) = 2*pi ;
    case 3
        w_0(3) = 2*pi ;
end

%plot
SimTime = 50;
ans4 = sim('EulerWheel.slx') ;
figure
plot(ans4.Omega)
ylabel('w [rad/s]')
legend( 'w_x', 'w_y', 'w_z' )
grid on

%% Ellipsoid plot

%angular momentum's ellipsoid
h = I*w_0';
h_norm = norm(h);
[Xh, Yh, Zh] = ellipsoid(0, 0, 0, h_norm/I_x, h_norm/I_y, h_norm/I_z );

%kinetic energy and its ellipsoid
T = 1/2.*w_0*I*w_0';
[XT, YT, ZT] = ellipsoid(0, 0, 0, sqrt(2*T/I_x), sqrt(2*T/I_y), sqrt(2*T/I_z));

%plot
figure
surf(Xh, Yh, Zh); hold on;
alpha 0.05
surf(XT, YT, ZT);
alpha 0.05






