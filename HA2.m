%% Home assignment 2
% Task 1: Plot the periodic potential in the plane (-1 1 0)

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

% Number of grid points
nPoints = 100;

% Get vector for the 3D space
for i = 1:3
    r(i,:) = linspace(1,10,nPoints);
end

% Initialise with zeros
V_r = zeros(nPoints, nPoints, nPoints);

% Calc potential in space
for k = 1:length(r)
    for l = 1:length(r)
        for m = 1:length(r)
            for i = 1:Gsize(1)
                V_r(k,l,m) =  V_r(k,l,m) + V_G(1,i) * cos(norm(G(i,:).*[r(1,k) r(2,l) r(3,m)]));
            end
        end
    end
end

% Get the potential in the plane
for i = 1:length(r)-1
    potInPlan(i,:) = V_r(i,i+1,:); 
end

% Plot the potential in the plane
surf(potInPlan)

xlabel('X')
ylabel('Y')
zlabel('Z')

shg

%% Task 2: Get band structure


k(1,:) = [0 0 0];
k(2,:) = pi/a*[-1 1 1];
k(3,:) = 2*pi/a*[0 0 1];