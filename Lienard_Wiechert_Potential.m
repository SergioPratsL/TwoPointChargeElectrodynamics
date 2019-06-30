function [ V, A ] = Lienard_Wiechert_Potential( R, v )
% R is the retarded distance 
% v is the charge velocity
% All what can be normalized (c, permitivity...) has been normalized

V = 1 /(norm(R) - dot(R,v));

A = v  / (norm(R) - dot(R,v));

end

