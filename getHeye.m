function [ Heye ] = getHeye( K, G )
%GETHEYE Summary of this function goes here
%   Detailed explanation goes here
hbar = 1;
m_e = 1;

factor = hbar^2/(2*m_e);

G_size = size(G);

Heye = eye(G_size(1));

for i = 1:G_size(1) 
    Heye(i,i) = factor*norm(K + G(i,:))^2; 
end

end

