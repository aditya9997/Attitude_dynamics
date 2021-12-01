%% 6U Cubesat : detumbling,slew manoeuvre and Sun pointing
%
% AUTHOR:
%
%   NAME:                                   ID code :      Badge number:
%
%	Andrea De Vittori                       10467544       898352
% 
% DATE:
%
%   20/11/2018, MATLAB
%
% -------------------------------------------------------------------------
%% DETUMBLING:

% Data:
close all
clear 
clc
disp('---> de-tumbling...')
J=zeros(1,3); 

% VAR                 DESCRIPTION              UNIT OF MEASURE:

% Satellite Inertia:

J(1)=0.218;          %Ixx                     [kg*m^2]
J(2)=0.166;          %Iyy                     [kg*m^2]
J(3)=0.082;          %Izz                     [kg*m^2]

% CG thrusters parameters:

x=0.02518/2;         % in plane position      [m]
c=0.183;             % half lenght of the sc  [m]
F=10e-3;             % Thrust of a thruster   [N]
T_R=10e-3;           % rise time              [s]
T_F=50e-3;           % fall time              [s]
T_d=5e-3;            % delay time             [s]

%Parameters for gyro:

f=262;
Ts=1/f;

sigma_n=0.15*pi/180/60/sqrt(Ts);
sigma_b=0.3*pi/180/3600/sqrt(Ts);

% parameters for the differential Riccati equation:
R=eye(3)*1e-8;
Q=eye(3)*1e-8;
S_0=eye(3);

%Parameters of magnetometer:

fM=18; 
Tm=1/fM;% refresh period of the magnetometer
pointing_accuracy=0.5;% degrees

%Initial Conditions for detumbling:

om_0=[0.22 0.26 0.22];
A0=[0.5335 0.808 0.25;....
    -0.808 0.3995 0.433;...
    0.25 -0.433 0.866];
% Time of integration:

T_det=180;

%Parameters for the Kalman filter:

C=1*eye(3); %  matrix of the sensors

% matrix for bang bang control :

T_max_1=2*c*F*sind(15);
T_max_2=2*c*F*cosd(15);
T_max_3=2*(x*F*cosd(15)-x*F*sind(15));
T_max=diag([T_max_1;T_max_2;T_max_3]);

%Recalling the Simulink file:

sim('detumbling')

% To workspace :

om_x_f = om_det_filter.Data(:,1:3:end-2); 
om_y_f = om_det_filter.Data(:,2:3:end-1);
om_z_f = om_det_filter.Data(:,3:3:end);
om_x_m = om_det_m.Data(:,1:3:end-2); 
om_y_m = om_det_m.Data(:,2:3:end-1);
om_z_m = om_det_m.Data(:,3:3:end);
Mc_x = Control_action.Data(:,1:3:end-2);
Mc_y = Control_action.Data(:,2:3:end-1);
Mc_z = Control_action.Data(:,3:3:end);
CG=CG.Data(:,:);
A_end=[A.Data(1,end-2) A.Data(1,end-1) A.Data(1,end);...
    A.Data(2,end-2) A.Data(2,end-1) A.Data(2,end);...
    A.Data(3,end-2) A.Data(3,end-1) A.Data(3,end)];
Time=  om_det_filter.Time;
Timecg=linspace(0,T_det,length(CG(:,1)));
figure()

% Plot angular velocities filtered:

subplot(1,2,1)
Legend=['$\omega_x$';'$\omega_y$';'$\omega_z$'];
PLOT(Time,'Time [s]','$\omega$ [rad/s]','$\omega_{det_{filtered}}$ ',1.2,[om_x_f;om_y_f;om_z_f],Legend)

% Plot angular velocities not filtered:
subplot(1,2,2)
Legend=['$\omega_x$';'$\omega_y$';'$\omega_z$'];
PLOT(Time,'Time [s]','$\omega$ [rad/s]','$\omega_{det_{noisy}}$ ',1.2,[om_x_m;om_y_m;om_z_m],Legend)

% Plot momentum due to CG thrusters:
figure()
Legend=['$M_x$';'$M_y$';'$M_z$'];
PLOT(Timecg,'Time [s]','Torque[Nm]','',1.2,[CG(:,1)';CG(:,2)';CG(:,3)'],Legend)
disp('---> end of de-detumbling...')
%% SLEW MANOEUVRE:
% Initial conditions:
 disp('--->  slew...')
om_0=[om_x_f(1,end) om_y_f(1,end) om_z_f(1,end)];
A0=A_end;

% Time fo integration and Delta_t shift related to detambling:
T_slew=140;
Delta_t=T_det;

% Maximum torque of a single reaction wheel:
m_max=1e-3;
momentum_storage=1e-2;
 sim('slew')
 
% To workspace :
om_x_s = om_slew.Data(:,1); 
om_y_s = om_slew.Data(:,2);
om_z_s = om_slew.Data(:,3);
SRP_x = SRP_slew.Data(:,1);
SRP_y = SRP_slew.Data(:,2);
SRP_z = SRP_slew.Data(:,3);
GG_x = GG_slew.Data(:,1);
GG_y = GG_slew.Data(:,2);
GG_z = GG_slew.Data(:,3);
MDT_x = MDT_slew.Data(:,1);
MDT_y = MDT_slew.Data(:,2);
MDT_z = MDT_slew.Data(:,3);
A_s11=A_slew.Data(1,1:3:end-2);A_s12= A_slew.Data(1,2:3:end-1);A_s13= A_slew.Data(1,3:3:end);
A_s21=A_slew.Data(2,1:3:end-2);A_s22= A_slew.Data(2,2:3:end-1);A_s23= A_slew.Data(2,3:3:end);
A_s31=A_slew.Data(3,1:3:end-2);A_s32= A_slew.Data(3,2:3:end-1);A_s33= A_slew.Data(3,3:3:end);

A_ss11=A_slew1.Data(1,1:3:end-2);A_ss12= A_slew1.Data(1,2:3:end-1);A_ss13= A_slew1.Data(1,3:3:end);
A_ss21=A_slew1.Data(2,1:3:end-2);A_ss22= A_slew1.Data(2,2:3:end-1);A_ss23= A_slew1.Data(2,3:3:end);
A_ss31=A_slew1.Data(3,1:3:end-2);A_ss32= A_slew1.Data(3,2:3:end-1);A_ss33= A_slew1.Data(3,3:3:end);

A_d11=A_d_slew.Data(1,1:3:end-2);A_d12= A_d_slew.Data(1,2:3:end-1);A_d13= A_d_slew.Data(1,3:3:end);
A_d21=A_d_slew.Data(2,1:3:end-2);A_d22= A_d_slew.Data(2,2:3:end-1);A_d23= A_d_slew.Data(2,3:3:end);
A_d31=A_d_slew.Data(3,1:3:end-2);A_d32= A_d_slew.Data(3,2:3:end-1);A_d33= A_d_slew.Data(3,3:3:end);

A_n_s11=A_noise_slew.Data(1,1:3:end-2);A_n_s12= A_noise_slew.Data(1,2:3:end-1);A_n_s13= A_noise_slew.Data(1,3:3:end);
A_n_s21=A_noise_slew.Data(2,1:3:end-2);A_n_s22= A_noise_slew.Data(2,2:3:end-1);A_n_s23= A_noise_slew.Data(2,3:3:end);
A_n_s31=A_noise_slew.Data(3,1:3:end-2);A_n_s32= A_noise_slew.Data(3,2:3:end-1);A_n_s33= A_noise_slew.Data(3,3:3:end);

Tot_dist=[GG_x'+SRP_x'+MDT_x';GG_y'+SRP_y'+MDT_y';GG_z'+SRP_z'+MDT_z'];
Tot_dist_1=zeros(1,length(Tot_dist(1,:)));
for i=1:length(Tot_dist(1,:))
    Tot_dist_1(i)=norm( Tot_dist(:,i));
end

SUN1=SUN.Data;
sc1=sc.Data;

Time_slew=  om_slew.Time;

figure()

% Plot angular velocities 
Legend=['$\omega_x$';'$\omega_y$';'$\omega_z$'];
PLOT(Time_slew','Time [s]','$\omega$ [rad/s]','$\omega_{slew}$ ',1.2,[om_x_s';om_y_s';om_z_s'],Legend)

% plot SRP:
figure()
Legend=['$SRP_x$';'$SRP_y$';'$SRP_z$'];
PLOT(Time_slew,'Time [s]','Torque [Nm]','$SRP_{slew}$',1.2,[SRP_x';SRP_y';SRP_z'],Legend)

% magnitude of the disturbing torques:
figure()
Legend=[''];
PLOT(Time_slew,'Time [s]','Torque [Nm]','Magnitude of the disturbing torques',1.2,[Tot_dist_1],Legend)

% plot A_slew:
figure()
% subplot(1,2,1)
Legend=['$A_{11}$';'$A_{12}$';'$A_{13}$';'$A_{21}$';'$A_{22}$';'$A_{23}$';...
        '$A_{31}$';'$A_{32}$';'$A_{33}$'];
PLOT(Time_slew,'Time [s]','','$A_{slew}$',1.2,[A_s11;A_s12;A_s13;A_s21;A_s22;A_s23;A_s31;A_s32;A_s33],Legend)

% plot A_slew_desired:
figure()
% subplot(1,2,2)
PLOT(Time_slew,'Time [s]','','$A_{slew_{d}}$',1.2,[A_d11;A_d12;A_d13;A_d21;A_d22;A_d23;A_d31;A_d32;A_d33],Legend)
%A_noisy:
figure()
PLOT(Time_slew,'Time [s]','','$A_{slew_{noisy}}$',1.2,[A_n_s11;A_n_s12;A_n_s13;A_n_s21;A_n_s22;A_n_s23;A_n_s31;A_n_s32;A_n_s33],Legend)

As_end=[A_ss11(end) A_ss12(end) A_ss13(end); ...
        A_ss21(end) A_ss22(end) A_ss23(end);...
         A_ss31(end) A_ss32(end) A_ss33(end)];
     
 As1_end=[A_s11(end) A_s12(end) A_s13(end); ...
         A_s21(end) A_s22(end) A_s23(end);...
         A_s31(end) A_s32(end) A_s33(end)];
% plot of the error:
figure()
Legend='err';
PLOT(Time_slew,'Time [s]','','Attitude error',1.2,err.Data',Legend)
% Plot angular momentum: 
figure()
Legend=['$h_1$';'$h_2$';'$h_3$';'$h_4$'];
PLOT(Time_slew,'Time [s]','$[Nms] $','$Angular \ momentum_{slew}$ ',1.2,[h_slew.Data(:,1)';h_slew.Data(:,2)';h_slew.Data(:,3)';h_slew.Data(:,4)'],Legend)

% Plot of the sphere:
figure()
background()
T=1:floor(length(A_s11)/50):length(A_s11);
[X,Y,Z] = sphere(60);
SPHERE_PLOT=surf(X,Y,Z);
 set( SPHERE_PLOT,'FaceAlpha',0.3);
       shading interp;
 quiver3(0,0,0,1,0,0,'color','w','linewidth',3);
quiver3(0,0,0,0,1,0,'color','w','linewidth',3);
quiver3(0,0,0,0,0,1,'color','w','linewidth',3);
text(0,0,1,texlabel('NZ'),'Color','w','FontSize',10,'FontWeight','bold');
text(0,1,0,texlabel('NY'),'Color','w','FontSize',10,'FontWeight','bold');
text(1,0,0,texlabel('NX'),'Color','w','FontSize',10,'FontWeight','bold');
text(sind(45),cosd(45),0,texlabel('Equatorial plane'),'Color','w')
text(A_s11(end)+0.1,A_s12(end)+0.1,A_s13(end)-0.02,texlabel('Ecliptic plane'),'Color','y')
circlePlane3D([0,0,0],[0,0,1], 1, 0.2, [ 255 255 255]/255)
hold on
circlePlane3D([0,0,0],[0,-sind(23.45),cosd(23.45)], 1, 0.2, [ 255 255 0]/255)
text(A_s11(1),A_s12(1),A_s13(1),texlabel('X_{in}'),'Color','r','FontSize',12,'FontWeight','bold');
text(A_s21(1),A_s22(1),A_s23(1),texlabel('Y_{in}'),'Color','g','FontSize',12,'FontWeight','bold');
text(A_s31(1),A_s32(1),A_s33(1),texlabel('Z_{in}'),'Color','b','FontSize',12,'FontWeight','bold');
for i=T
    x1=(SUN1(i,:)+sc1(i,:))/norm(SUN1(i,:)+sc1(i,:));
    
if i<T(end)
quiver3(0,0,0,A_s11(i),A_s12(i),A_s13(i),'color','r');
hold on
quiver3(0,0,0,A_s21(i),A_s22(i),A_s23(i),'color','g');
quiver3(0,0,0,A_s31(i),A_s32(i),A_s33(i),'color','b');
else
 quiver3(0,0,0,A_s11(i),A_s12(i),A_s13(i),'color','r','linewidth',2);
hold on
quiver3(0,0,0,A_s21(i),A_s22(i),A_s23(i),'color','g','linewidth',2);
quiver3(0,0,0,A_s31(i),A_s32(i),A_s33(i),'color','b','linewidth',2); 
end
h=scatter3(x1(1),x1(2),x1(3),80,'y','filled');
axis equal
drawnow
if i<T(end)
delete(h)
end
end
text(A_s11(end),A_s12(end),A_s13(end),texlabel('X_{end}'),'Color','r','FontSize',12,'FontWeight','bold');
text(A_s21(end),A_s22(end),A_s23(end),texlabel('Y_{end}'),'Color','g','FontSize',12,'FontWeight','bold');
text(A_s31(end),A_s32(end),A_s33(end),texlabel('Z_{end}'),'Color','b','FontSize',12,'FontWeight','bold');

 disp('---> end of slew...')
 %% SUN POINTING:
  disp('---> tracking...')
 % Initial Conditions:
 
om_0_p=[om_x_s(end,1) om_y_s(end,1) om_z_s(end,1) ];
A0_p=As_end;% for the ideal attitude
A0_p1=As1_end;% for the reconstructed state
h_0_p=[h_slew.Data(end,1);h_slew.Data(end,2);h_slew.Data(end,3);...
    h_slew.Data(end,4)];% initial conditon for the reaction wheels dynamics = last value in slew

% Time of integration  integration and Delta_t shift related to detambling+slew:
T_pointing=86400;
Delta_t_p=T_det+T_slew;
m_max_p=1e-3;

sim('Sun_pointing1')

% To workspace

om_x_p = om_p.Data(:,1); 
om_y_p = om_p.Data(:,2);
om_z_p = om_p.Data(:,3);
SRP_x_p = SRP_p.Data(:,1);
SRP_y_p = SRP_p.Data(:,2);
SRP_z_p = SRP_p.Data(:,3);
GG_x_p = GG_p.Data(:,1);
GG_y_p = GG_p.Data(:,2);
GG_z_p = GG_p.Data(:,3);
MDT_x_p = MDT_p.Data(:,1);
MDT_y_p = MDT_p.Data(:,2);
MDT_z_p = MDT_p.Data(:,3);
A_p11=A_p.Data(1,1:3:end-2);A_p12= A_p.Data(1,2:3:end-1);A_p13= A_p.Data(1,3:3:end);
A_p21=A_p.Data(2,1:3:end-2);A_p22= A_p.Data(2,2:3:end-1);A_p23= A_p.Data(2,3:3:end);
A_p31=A_p.Data(3,1:3:end-2);A_p32= A_p.Data(3,2:3:end-1);A_p33= A_p.Data(3,3:3:end);

A_d11=A_d_p.Data(1,1:3:end-2);A_d12= A_d_p.Data(1,2:3:end-1);A_d13= A_d_p.Data(1,3:3:end);
A_d21=A_d_p.Data(2,1:3:end-2);A_d22= A_d_p.Data(2,2:3:end-1);A_d23= A_d_p.Data(2,3:3:end);
A_d31=A_d_p.Data(3,1:3:end-2);A_d32= A_d_p.Data(3,2:3:end-1);A_d33= A_d_p.Data(3,3:3:end);

A_n_p11=A_noise_p.Data(1,1:3:end-2);A_n_p12= A_noise_p.Data(1,2:3:end-1);A_n_p13= A_noise_p.Data(1,3:3:end);
A_n_p21=A_noise_p.Data(2,1:3:end-2);A_n_p22= A_noise_p.Data(2,2:3:end-1);A_n_p23= A_noise_p.Data(2,3:3:end);
A_n_p31=A_noise_p.Data(3,1:3:end-2);A_n_p32= A_noise_p.Data(3,2:3:end-1);A_n_p33= A_noise_p.Data(3,3:3:end);

Tot_dist_p=[GG_x_p'+SRP_x_p'+MDT_x_p';GG_y_p'+SRP_y_p'+MDT_y_p';GG_z_p'+SRP_z_p'+MDT_z_p'];
Tot_dist_1_p=zeros(1,length(Tot_dist_p(1,:)));
for i=1:length(Tot_dist_p(1,:))
    Tot_dist_1_p(i)=norm( Tot_dist_p(:,i));
end
SUN=SUN.Data;
sc=sc.Data;

Time_p= A_d_p.Time;

% Plot angular velocities 
figure()
Legend=['$\omega_x$';'$\omega_y$';'$\omega_z$'];
PLOT(Time_p','Time [s]','$\omega$ [rad/s]','$\omega_{pointing}$ ',1.2,[om_x_p';om_y_p';om_z_p'],Legend)

% plot SRP:
figure()
Legend=['$SRP_x$';'$SRP_y$';'$SRP_z$'];
PLOT(Time_p,'Time [s]','Torque [Nm]','$SRP_{pointing}$',1.2,[SRP_x_p';SRP_y_p';SRP_z_p'],Legend)

% plot GG torque:
figure()
Legend=['$GG_x$';'$GG_y$';'$GG_z$'];
PLOT(Time_p,'Time [s]','Torque [Nm]','$GG_{pointing}$',1.2,[GG_x_p';GG_y_p';GG_z_p'],Legend)

% plot magnetic disturbance:
figure()
Legend=['$MDT_x$';'$MDT_y$';'$MDT_z$'];
PLOT(Time_p,'Time [s]','Torque [Nm]','$Magnetic \ disturbing\ torque_{pointing}$',1.2,[MDT_x_p';MDT_y_p';MDT_z_p'],Legend)
% magnitude of the disturbing torques:
figure()
Legend=[''];
PLOT(Time_p,'Time [s]','Torque [Nm]','Magnitude of the disturbing torques',1.2,[Tot_dist_1_p],Legend)

% plot A_pointing:
figure()
 subplot(1,2,1)
Legend=['$A_{11}$';'$A_{12}$';'$A_{13}$';'$A_{21}$';'$A_{22}$';'$A_{23}$';...
        '$A_{31}$';'$A_{32}$';'$A_{33}$'];
PLOT(Time_p,'Time [s]','','$A_{pointing}$',0.6,[A_p11;A_p12;A_p13;A_p21;A_p22;A_p23;A_p31;A_p32;A_p33],Legend)

% plot A_pointing_desired:

 subplot(1,2,2)
PLOT(Time_p,'Time [s]','','$A_{pointing_{d}}$',0.6,[A_d11;A_d12;A_d13;A_d21;A_d22;A_d23;A_d31;A_d32;A_d33],Legend)

% plot A_noisy:
figure()
PLOT(Time_p,'Time [s]','','$A_{pointing_{noisy}}$',0.6,[A_n_p11;A_n_p12;A_n_p13;A_n_p21;A_n_p22;A_n_p23;A_n_p31;A_n_p32;A_n_p33],Legend)

% plot angular momentum for pointing:
figure()
Legend=['$h_1$';'$h_2$';'$h_3$';'$h_4$'];
PLOT(Time_p,'Time [s]','$[Nms] $','$Angular \ momentum_{pointing}$ ',1.2,[h_p.Data(:,1)';h_p.Data(:,2)';h_p.Data(:,3)';h_p.Data(:,4)'],Legend)
% PLOT of the sphere:
figure()
background()

T=1:floor(length(A_p11)/50):length(A_p11);
[X,Y,Z] = sphere(60);
SPHERE_PLOT=surf(X,Y,Z);
 set( SPHERE_PLOT,'FaceAlpha',0.3);
       shading interp;
 quiver3(0,0,0,1,0,0,'color','w','linewidth',2);
quiver3(0,0,0,0,1,0,'color','w','linewidth',2);
quiver3(0,0,0,0,0,1,'color','w','linewidth',2);
text(0,0,1,texlabel('NZ'),'Color','w','FontSize',12,'FontWeight','bold');
text(0,1,0,texlabel('NY'),'Color','w','FontSize',12,'FontWeight','bold');
text(1,0,0,texlabel('NX'),'Color','w','FontSize',12,'FontWeight','bold');
text(sind(45),cosd(45),0,texlabel('Equatorial plane'),'Color','w')
text(A_p11(1),A_p12(1),A_p13(1),texlabel('Ecliptic plane'),'Color','y')
circlePlane3D([0,0,0],[0,0,1], 1, 0.2, [ 255 255 255]/255)
hold on
circlePlane3D([0,0,0],[0,-sind(23.45),cosd(23.45)], 1, 0.2, [ 255 255 0]/255)

for i=T
    x1=(SUN(i,:)+sc(i,:))/norm(SUN(i,:)+sc(i,:));
    
if i<T(end)
 quiver3(0,0,0,A_p11(i),A_p12(i),A_p13(i),'color','r');
hold on
quiver3(0,0,0,A_p21(i),A_p22(i),A_p23(i),'color','g');
quiver3(0,0,0,A_p31(i),A_p32(i),A_p33(i),'color','b');
else
quiver3(0,0,0,A_p11(i),A_p12(i),A_p13(i),'color','r','linewidth',2);
hold on
quiver3(0,0,0,A_p21(i),A_p22(i),A_p23(i),'color','g','linewidth',2);
quiver3(0,0,0,A_p31(i),A_p32(i),A_p33(i),'color','b','linewidth',2);
end
h=scatter3(x1(1),x1(2),x1(3),80,'y','filled');
axis equal
drawnow
if i<T(end)
delete(h)
end
end
text(A_p11(end),A_p12(end),A_p13(end),texlabel('X'),'Color','r','FontSize',12,'FontWeight','bold');
text(A_p21(end),A_p22(end),A_p23(end),texlabel('Y'),'Color','g','FontSize',12,'FontWeight','bold');
text(A_p31(end),A_p32(end),A_p33(end),texlabel('Z'),'Color','b','FontSize',12,'FontWeight','bold');
 disp('---> end of tracking...')
%% PLOT FUNCTION:
 function PLOT(time,Xlabel,Ylabel,Title,width,xData,Legend)
 for i=1:length(xData(:,1))
 h=plot(time,xData(i,:));
 hold on
set(h,'linewidth',width); 
if i==1
h = title(Title,'interpreter','latex');
set(h,'FontSize',20);
h = ylabel(Ylabel,'interpreter','latex');
set(h,'FontSize',15)
h = xlabel(Xlabel,'interpreter','latex');
set(h,'FontSize',15)
end

 end
 h=legend(Legend,'interpreter','latex');
 set(h,'FontSize',13)
 end
 %% BACKGROUND:
 function[]= background
% background.m - function which allows the user to set the background
% appearence
% 
% DESCRIPTION:
% 
%       1- set the properties of the background of a 2D/3D plot
%
%   1- CHANGE THE PROPERTIES OF THE BACKGROUND:
%   First of all the function creates the 'background' axes.Then 
%   it moves them  to the bottom. Load in a background  the image that 
%    you saved previously in the current folder of you code and display it
%    using the correct colors.
%    Subsequently it Turns the handlevisibility off so that we don't 
%    inadvertently plot into the axes again.
%     Also,it makes the axes invisible,hiding the frame system.
%
% INPUT:
%	the function has no input         
%
% OUTPUT:
%   the function has no output
%
% AUTHOR:
%   Andrea De  vittori, 28/12/2017, MATLAB, lbackgroundambertMR.m
%

% -------------------------------------------------------------------------

%  1- CHANGE THE PROPERTIES OF THE BACKGROUND:
ha = axes('units','normalized', ...
            'position',[0 0 1 1]);

uistack(ha,'bottom');
%load the image:
I=imread('nero.jpg');
imagesc(I);
set(ha,'handlevisibility','off', ...
            'visible','off')
% hide the axis
        set(gcf, 'Color', 'None')
        axis('off')
        hold on
 end
 %% circle plane
 function H = circlePlane3D( center, normal, radious, theintv, color  )
%circlePlane3D.m -this function generates a circle plane in a 3D space.
%
% DESCRIPTION:
% findcolor is made by 2 parts
%       1-case of a plane with normal [0 0 +/-1]
%       2-generic case:
%
%   1- CASE NORMAL [0 0 +/-1]
%      the reason why there's this first part of the code is due to the
%      presence of a cross product between the z-axis and the normal to
%      the plane,that in case of parallel vectors gives as output an error.
%      the other part of this portion of the function works as the generic
%      case
%
%   2-GENERIC CASE
%     steps:
%     1)define an angle between0 and 2pi
%     2)generate circle polygon
%     3)define he angle between the norma and the z axis
%     4)apply the rotation matrix
%     5)translate the centre
% INPUT:
%   VARIABLE:      DIM:   DESCRIPTION:                                        UNIT OF MEASURE:
%
%	center         [3]    position of the centre of the cicle plane            [km]                         
%   normal         [3]    direction of the normal to the circle plane          
%   radious        [1]    radious of the circle plane                          [km]
%   theintv        [1]    interval theta which allow you to control 
%                         your polygon
%   color         [char]  color of the circle plane 3D                                
%,
% OUTPUT:
%   VARIABLE:      DIM:   DESCRIPTION:                                        UNIT OF MEASURE:
%
%	H            [figure]  plot of the circl plane 3D                                     
%
% REFERENCES:
%   - 
% AUTHOR:
%   Andrea de vittori, 13/01/2017, MATLAB, circlePlane3D.m
%--------------------------------------------------------------------------
%%   1- CASE NORMAL [0 0 +/-1]
%generate circle polygon
if isequal(normal,[0 0 1])==1||isequal(normal,[0 0 -1])==1
    t = 0:theintv:2*pi;
x = radious*cos(t);
y = radious*sin(t);
z = zeros(size(x));
%compute rotate theta and axis
axis = [0 0 1];
normal = normal/norm(normal);
ang = acos(dot(axis,normal));
% A skew symmetric representation of the normalized axis 
axis_skewed = [ 0 -axis(3) axis(2) ; axis(3) 0 -axis(1) ; -axis(2) axis(1) 0]; 
%  rotation matrix 
R = eye(3) + sin(ang)*axis_skewed + (1-cos(ang))*axis_skewed*axis_skewed;
fx = R(1,1)*x + R(1,2)*y + R(1,3)*z;
fy = R(2,1)*x + R(2,2)*y + R(2,3)*z;
fz = R(3,1)*x + R(3,2)*y + R(3,3)*z;
%translate center
fx = fx+center(1);
fy = fy+center(2);
fz = fz+center(3);
H = fill3(fx, fy, fz, color);
set(H,'FaceAlpha',0.25)
else
t = 0:theintv:2*pi;
x = radious*cos(t);
y = radious*sin(t);
z = zeros(size(x));
%compute rotate theta and axis
zaxis = [0 0 1];
normal = normal/norm(normal);
ang = acos(dot(zaxis,normal));
axis = cross(zaxis, normal)/norm(cross(zaxis, normal));
% A skew symmetric representation of the normalized axis 
axis_skewed = [ 0 -axis(3) axis(2) ; axis(3) 0 -axis(1) ; -axis(2) axis(1) 0]; 
%  rotation matrix 
R = eye(3) + sin(ang)*axis_skewed + (1-cos(ang))*axis_skewed*axis_skewed;
fx = R(1,1)*x + R(1,2)*y + R(1,3)*z;
fy = R(2,1)*x + R(2,2)*y + R(2,3)*z;
fz = R(3,1)*x + R(3,2)*y + R(3,3)*z;
%translate center
fx = fx+center(1);
fy = fy+center(2);
fz = fz+center(3);
H = fill3(fx, fy, fz, color);
set(H,'FaceAlpha',0.25)

end

end