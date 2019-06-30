function [Eo, Bo] = RadiatedFieldAcceleratedCharge(Position, Velocity, Acceleration)
% Calculates the radiated field caused by a charge of unit charge
% under an 'Acceleration' acceleration, moving at 'Velocity' velocity
% evaluated at a point at 'Position' measured from the charge at the
% retarded position (Position is a 3-dim space vector).

% This calculations are done ignoring permitivity and 4*pi factors 
% and normalizing the speed of light, c=1.

Distance = norm(Position);

VectorProduct = cross( Position, cross((Position - Distance*Velocity), Acceleration));

Eo = 1 / (Distance - dot(Position,Velocity))^3 * VectorProduct;

Bo = cross( Position, Eo ) / Distance;


end
