function [a, e, i, OM, om, th] = car2kep(r, v, mu)
% car 2 kep m Conversion from Cartesian coordinates to Keplerian elements
%
% PROTOTYPE
%
%[e, i OM, om, th car 2 kep(r, v, mu)
%
% DESCRIPTION
% Conversion from Cartesian coordinates to Keplerian elements Angles in radians
%
% INPUT
% r [3 x 1] Position vector [km]
% v [3 x 1] Velocity vector [km/s]
% mu [1 x 1] Gravitational parameter [km^3/s^2]
%
% OUTPUT
% a [1 x 1] Semi major axis [km]
% e [1 x 1] Eccentricity [-]
% i [1 x 1] Inclination [rad]
% OM [1 x 1] RAAN [rad]
% om [1 x 1] Pericentre anomaly [rad]
% th [1 x 1] True anomaly [rad]

% calcolo la norma di r
r_norm=norm(r);
%calcolo la norma di v
v_norm=norm(v);
%calcolo il momento angolare specifico
%nota: cross è il prodotto vettore
h=cross(r, v);
%calcolo il modulo del momento angolare specifico
h_norm=norm(h);
%calcolo l'inclinazione
    %calcolo l'elemento in posizione z di h, cioe' hz
    h_z=h(3,1);
%calcolo ora l'inclinazione
i = acos(h_z/h_norm);
%calcolo il vettore eccentricita'
e= 1/mu * [((v_norm)^2 - (mu/r_norm)).*r - (dot(r, v).*v)];
%calcolo il modulo di e
e_norm= norm(e);
%calcolo energia meccanica
E=0.5*(v_norm)^2 - (mu/r_norm);
%calcolo il semiasse maggiore
a= -(mu/(2*E));
%Linea dei nodi
    %definisco versore K
    K=[0; 0; 1];
%calcolo vettore linea dei nodi
N=cross(K, h);
%calcolo norma di N
N_norm= norm(N);
%Calcolo l'ascensione retta del nodo ascendente
    %Calcolo le componenti Nx e Ny
    N_x=N(1, 1);
    N_y=N(2, 1);
%Calcolo ascensione retta del nodo ascendente (OM)
%Metto un if per evitare l'ambguita' dovuto all'
%arcocoseno
if(N_y >=0)
    OM=acos(N_x/N_norm);
else
    OM=(2*pi-acos(N_x/N_norm));
end
%Calcolo l'anomalia del pericentro
    %estraggo la componente e_z
    e_z= e(3, 1);
%calcolo l'anomalia del pericentro (om) 
%Metto un if per evitare l'ambiguita' dovuta all'
%arcocoseno
if(e_z >= 0)
    om=acos((dot(N, e)/(N_norm*e_norm)));
else
    om=2*pi - acos((dot(N, e)/(N_norm*e_norm)));
end
%calcolo la velocita' radiale
%nota: se Vr > 0, ci stiamo allontanando dal pericentro
%Se Vr < 0, ci stiamo avvicinando al pericentro
%perchè la lunghezza del raggio vettore diminuisce 
%nel tempo
v_r=(dot(r, v)/r_norm);
%calcolo l'anomalia vera (th), uso un if per evitare 
%l'ambiguita' dovuta all'arcoseno
if(v_r >=0)
    th=acos(dot(e, r)/(e_norm*r_norm));
else
    th=2*pi - acos(dot(e, r)/(e_norm*r_norm));
end
end

