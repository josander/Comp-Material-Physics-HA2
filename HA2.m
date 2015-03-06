%% Home assignment 2
clc
clf
clear all

% Lattice parameter [Å]
a = 5.43;   

% Basis vectors
d(1,:) = a/8*[1 1 1];         
d(2,:) = -a/8*[1 1 1]; 

% Define G-vectors for which the structure factor > zero
G(1,:) = 2*pi/a*[0 0 0];
G(2,:) = 2*pi/a*[1 1 1];
G(3,:) = 2*pi/a*[1 1 3];
G(4,:) = 2*pi/a*[1 3 3];

% Get sizes of the arrays
Gsize = size(G);
dSize = size(d);

% Get the structure factor in reciprocal space
for i = 1:Gsize(1)
    for j = 1:dSize(1)
        S_G(i) = cos(norm(G(i,:).*d(j,:)));
    end
end

% Define form factor for the different G-vectors
v_G = [0 -0.056 0.0138 0.0181];

% Get potential in reciprocal space
V_G = v_G.*S_G;

% Get vector for the 3D space
for i = 1:3
    r(i,:) = linspace(1,10,100);
end

% Calc potential in space. j: distance from origo in one direction. i:sum
% over different G-vectors
for j = 1:length(r)
    for i = 1:Gsize(1)
        V_r(i,j) =  V_G(1,i) * cos(norm(G(i,:).*r(:,j)'));
    end
end

%%
for i = 1:length(r)-1
    plot3(V_r(1,i), V_r(2,i+1), V_r(3,:))
    hold on
end

hold off

xlabel('X')
ylabel('Y')
zlabel('Z')



