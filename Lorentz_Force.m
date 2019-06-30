function F = Lorentz_Force(v, E, B)
% Lorentz force with speed of light, charge, permeativity and other things
% normalized.

F = E + cross(v, B);

