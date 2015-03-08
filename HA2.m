%% Home assignment 2
% Task 1: Plot the periodic potential in the plane (-1 1 0)

clc
clear all

% Lattice parameter [Å]
a = 5.43;   

% Basis vectors
d(1,:) = a/8*[1 1 1];         
d(2,:) = -a/8*[1 1 1]; 

% Maxvalue in the G-vector
maxValue = 3;

% Define G-vectors for which the structure factor > zero
[G v_G] = constructG(a, maxValue);

% Get sizes of the arrays
Gsize = size(G);
dSize = size(d);

% Initialize with zeros
S_G = zeros(1,Gsize(1));

% Get the structure factor in reciprocal space
for i = 1:Gsize(1)
        S_G(i) = 2*cos(G(i,:)*d(1,:)');
end

% Get potential in reciprocal space
V_G = v_G.*S_G;

% Number of grid points
nPoints = 50;

% Get vector for the 3D space
r = linspace(-8,8,nPoints);

% Initialise with zeros
V_r = zeros(nPoints, nPoints, nPoints);

tic

% Calc potential in space
for l = 1:length(r)
    for m = 1:length(r)
        for n = 1:length(r)
            for i = 1:Gsize(1)
                V_r(l,m,n) =  V_r(l,m,n) + V_G(1,i) * cos(G(i,:)*[r(l) r(m) r(n)]');
            end
        end
    end
end

toc

% Get the potential in the plane
for i = 1:length(r)
    potInPlan(i,:) = real(V_r(i,i,:)); 
end

% Plot the potential in the plane
contourf(r(:)*sqrt(2), r(:), potInPlan',20,'LineStyle','none')
colorbar

xlabel('Z, l = 0')
ylabel('X = Y, h = k')

%% Task 1: Nice plot

set(gcf,'renderer','painters','PaperPosition',[0 0 12 7]);
contourf(r(:)/sqrt(2), r(:), potInPlan',20,'LineStyle','none')
colorbar

plotTickLatex2D

title('Empirical pseudopotential in Si (-1 1 0)','Interpreter','latex', 'fontsize', 14);
X = xlabel('X = Y,  h = k [\AA]', 'Interpreter','latex', 'fontsize', 12);
Y = ylabel('Z, l = 0 [\AA]','Interpreter','latex', 'fontsize', 12);
set(Y, 'Units', 'Normalized', 'Position', [-0.09, 0.5, 0]);
set(X, 'Units', 'Normalized', 'Position', [0.5, -0.05, 0]);

print(gcf,'-depsc2','task1.eps')

%% Task 2: Get band structure

clc

% Lattice parameter [Å]
a = 5.43;   

% Constants
hbar = 1;                   % [au]
me = 1;                     % [au]

% Generate k-vectors
k(1,:) = [0 0 0];
k(2,:) = pi/a*[-1 1 1];
k(3,:) = 2*pi/a*[0 0 1];

% Basis vectors
d(1,:) = a/8*[1 1 1];         
d(2,:) = -a/8*[1 1 1]; 

% Maxvalue in the G-vector
maxValue = 2;

% Define E cut off
Ecut = 1000;              % [au]

% Define G-vectors for which the structure factor > zero
[G] = constructGbig(a, maxValue, k(3,:), Ecut);

Gsize = size(G)
(2*maxValue+1)^3+1

dG = [];

tic
for m = 1:Gsize(1)
    for n = 1:Gsize(1)
        dG(m,n) = norm(G(m,:) - G(n,:));
    end
end
toc

%%
tic
V_dG = getV_dG(V_r, dG, r);
toc


Ham = 1;


% The symmetry points in k-space
Gamma = 2*pi/a*[0 0 0];
X = 2*pi/a*[1 0 0];
W = 2*pi/a*[1 0.5 1];
L = 2*pi/a*[0.5 0.5 0.5];
K = 2*pi/a*[0.75 0.75 0];

%% Task 3: The valence electron charge density


