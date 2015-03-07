%% Home assignment 2
% Task 1: Plot the periodic potential in the plane (-1 1 0)

clc
clear all

% Lattice parameter [Å]
a = 5.43;   

% Basis vectors
d(1,:) = a/8*[1 1 1];         
d(2,:) = -a/8*[1 1 1]; 

% Define G-vectors for which the structure factor > zero
G = getG(a);

biggest = norm(G(end,:));

% Check if G vector is wrong
for u = 1:length(G)-1
    if( norm(G(u,:)) > biggest)
        output = ['Error in G, row ', num2str(u)];
    	disp(output);
    end
end

% Define form factor for the different G-vectors
v_G = getFormFact(G, a);


%%
% Get sizes of the arrays
Gsize = size(G);
dSize = size(d);

% Initialize with zeros
S_G = zeros(1,Gsize(1));

% Get the structure factor in reciprocal space
for i = 1:Gsize(1)
    for j = 1:dSize(1)
        %S_G(i) = S_G(i) + exp(-i*(G(i,:)*d(j,:)'));
        S_G(i) = S_G(i) + cos(G(i,:)*d(j,:)');
    end
end


% Get potential in reciprocal space
V_G = v_G.*S_G;

% Number of grid points
nPoints = 100;

% Get vector for the 3D space
for i = 1:3
    r(i,:) = linspace(-5,5,nPoints);
end

% Initialise with zeros
V_r = zeros(nPoints, nPoints, nPoints);

% Calc potential in space
for l = 1:length(r)
    for m = 1:length(r)
        for n = 1:length(r)
            for i = 1:Gsize(1)
                %V_r(l,m,n) =  V_r(l,m,n) + V_G(1,i) * real(exp(i*G(i,:)*[r(1,l) r(2,m) r(3,n)]'));
                V_r(l,m,n) =  V_r(l,m,n) + V_G(1,i) * cos(G(i,:)*[r(1,l) r(2,m) r(3,n)]');
            end
        end
    end
end

% Get the potential in the plane
for i = 1:length(r)
    potInPlan(i,:) = V_r(i,i,:); 
end

figure(1)

% Plot the potential in the plane
surf(r(1,:), r(1,:), potInPlan)
%contourf(V_r)

xlabel('X')
ylabel('Y')
zlabel('Z')

shg

%% Task 2: Get band structure

clc

% Generate k-vectors
k(1,:) = [0 0 0];
k(2,:) = pi/a*[-1 1 1];
k(3,:) = 2*pi/a*[0 0 1];

% Generate a lot of G-vectors
% FCC: all ever or all odd
G(1,:) = 2*pi/a*[1 1 1];
G(2,:) = 2*pi/a*[1 1 3];
G(3,:) = 2*pi/a*[1 3 1];
G(4,:) = 2*pi/a*[3 1 1];
G(5,:) = 2*pi/a*[2 2 2];
G(6,:) = 2*pi/a*[3 3 1];
G(7,:) = 2*pi/a*[1 3 3];
G(8,:) = 2*pi/a*[3 1 3];
G(9,:) = 2*pi/a*[2 2 4];
G(10,:) = 2*pi/a*[2 4 2];
G(11,:) = 2*pi/a*[4 2 2];
G(12,:) = 2*pi/a*[3 3 3];
G(13,:) = 2*pi/a*[1 5 1];
G(14,:) = 2*pi/a*[1 1 5];
G(15,:) = 2*pi/a*[5 1 1];
G(16,:) = 2*pi/a*[1 3 5];
G(17,:) = 2*pi/a*[5 3 1];
G(18,:) = 2*pi/a*[1 5 3];
G(19,:) = 2*pi/a*[5 1 3];
G(20,:) = 2*pi/a*[3 1 5];
G(21,:) = 2*pi/a*[3 5 1];
G(22,:) = 2*pi/a*[4 2 4];
G(23,:) = 2*pi/a*[2 4 4];
G(24,:) = 2*pi/a*[4 4 2];
G(25,:) = 2*pi/a*[5 3 3];
G(26,:) = 2*pi/a*[3 3 5];
G(27,:) = 2*pi/a*[3 5 3];

biggest = norm(G(length(G),:));

% Check if G vector is wrong
for u = 1:length(G)-1
    if( norm(G(u,:)) > biggest)
        output = ['Error in G, row ', num2str(u)];
    	disp(output);
    end
end

Gamma = 2*pi/a*[0 0 0];
X = 2*pi/a*[1 0 0];
W = 2*pi/a*[1 0.5 1];
L = 2*pi/a*[0.5 0.5 0.5];
K = 2*pi/a*[0.75 0.75 0];

%% Task 3: The valence electron charge density


