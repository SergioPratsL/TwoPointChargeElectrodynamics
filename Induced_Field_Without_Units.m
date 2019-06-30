function [Eo, Bo] = Induced_Field_Without_Units(R, v)
% Induced field based on the Lienard Wiechert potentials
% http://en.wikipedia.org/wiki/Li%C3%A9nard%E2%80%93Wiechert_potential

% R is the retarded distance 
% v is the charge velocity with units such as the speed of light is 1

% Ligth speed is normalized
% c = 2.9979 * 10^8;
c = 1;   % :P

% Permitivity (Farads by meter)
% E_0 = 8.8542 * 10^-12;
% I will use        E_0 = 1;        % :P

v_norm = v / c;

R_norm = R / norm(R);

Sigma = 1 / sqrt( 1 - norm(v_norm)^2);

% Coef_Eo = 1/(4*pi*E_0) * q * (1/Sigma^2) * 1/( 1 - dot(R_norm, v_norm) )^3 * 1 / norm(R)^2;
% Constants out!!
Coef_Eo = (1/Sigma^2) * 1/( 1 - dot(R_norm, v_norm) )^3 * 1 / norm(R)^2;

Eo = Coef_Eo * (R_norm - v_norm);

Bo = cross( R_norm,  Eo);

