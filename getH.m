function [ H ] = getH(a,  K, G )
%GETHEYE Function that builds a hamiltonian
%   a =  lattice parameter in a.u
%   K = k-vector
%   G = G-matrix 
hbar = 1;
me = 1;

factor = hbar^2/(2*me);
Gsize = size(G);
H = zeros(Gsize(1), Gsize(1));

V = zeros(Gsize(1),Gsize(1));

for i = 1:Gsize(1) 
    H(i,i) = factor*norm(K + G(i,:))^2; 
    for j = i:Gsize(1)
        Gij = G(i,:) - G(j,:);
        V(i,j) = getVG(a,Gij);
        V(j,i) = conj(V(i,j));
    end
end

H = H + V;


end

