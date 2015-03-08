function [ Heye ] = getHeye( K, G )
%GETHEYE Summary of this function goes here
%   Detailed explanation goes here
hbar = 1;
me = 1;

factor = hbar^2/(2*me);

Gsize = size(G);

Heye = eye(Gsize(1));

for i = 1:Gsize(1) 
    Heye(i,i) = factor*norm(K + G(i,:))^2; 
end

end

