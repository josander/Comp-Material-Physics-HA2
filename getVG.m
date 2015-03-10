function VG = getVG(a,G)
%GETVG Function that gets the potental
%   a =  lattice parameter in a.u
%   G = G-point 

%Defining constants
d = [1 1 1; -1 -1 -1].*(a/8);
vbase = [-0.0560 0.0138 0.0181];

%Calculate S for given G
S = exp(-1i*G*d(1,:)') + exp(-1i*G*d(2,:)');

%Setting vG
G = G./(2*pi/a);
if round(norm(G)^2) == 3
    VG = vbase(1)*S;
elseif round(norm(G)^2) == 8
    VG = vbase(2)*S;
elseif round(norm(G)^2) == 11
    VG = vbase(3)*S;
else
    VG = 0;
end

