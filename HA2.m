%% Home assignment 2
% Task 1: Plot the periodic potential in the plane (-1 1 0)
set(0, 'defaultTextInterpreter', 'latex');

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
r = linspace(-6,6,nPoints);

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

set(gcf,'renderer','painters','PaperPosition',[0 0 12 8]);
contourf(r(:)/sqrt(2), r(:), potInPlan',35,'LineStyle','none')
colorbar

plotTickLatex2D

title('Empirical pseudopotential in Si (-1 1 0)','Interpreter','latex', 'fontsize', 14);
X = xlabel('X = Y [\AA],  h = k', 'Interpreter','latex', 'fontsize', 12);
Y = ylabel('Z [\AA], l = 0','Interpreter','latex', 'fontsize', 12);
set(Y, 'Units', 'Normalized', 'Position', [-0.09, 0.5, 0]);
set(X, 'Units', 'Normalized', 'Position', [0.5, -0.055, 0]);

print(gcf,'-depsc2','task1.eps')

%% Task 2: Get band structure
% Find convergence

clc
clear all

% Lattice parameter [au]
a = 5.43/0.529177;   

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
maxValue = 3;

% Define E cut off
EcutInitial = 4;
EcutFinal = 20;              % [au]
dE = 0.25;

for Ecut = EcutInitial:dE:EcutFinal
    
    % Define G-vectors for which the structure factor > zero
    [G1] = constructGbig(a, maxValue, k(1,:), Ecut);
    [G2] = constructGbig(a, maxValue, k(2,:), Ecut);
    [G3] = constructGbig(a, maxValue, k(3,:), Ecut);

    % Get dG and form factor for dG
    [dG1, v_dG1] = getDG(G1);
    [dG2, v_dG2] = getDG(G2);
    [dG3, v_dG3] = getDG(G3);

    % Get sizes of the dG-matrices
    dGsize1 = size(dG1);
    dGsize2 = size(dG2);
    dGsize3 = size(dG3);
    
    % Allocate memory for struc factor
    S_dG1 = zeros(dGsize1(1),dGsize1(1));
    S_dG2 = zeros(dGsize2(1),dGsize2(1));
    S_dG3 = zeros(dGsize3(1),dGsize3(1));

    % Get the structure factor in reciprocal space
    for i = 1:dGsize1(1)
        for j = 1:dGsize1(1)
            S_dG1(i,j) = 2*cos([dG1(i,j,1) dG1(i,j,2) dG1(i,j,3)]*d(1,:)');
        end
    end
    
    for i = 1:dGsize2(1)
        for j = 1:dGsize2(1)
            S_dG2(i,j) = 2*cos([dG2(i,j,1) dG2(i,j,2) dG2(i,j,3)]*d(1,:)');
        end
    end
    
    for i = 1:dGsize3(1)
        for j = 1:dGsize3(1)
            S_dG3(i,j) = 2*cos([dG3(i,j,1) dG3(i,j,2) dG3(i,j,3)]*d(1,:)');
        end
    end

    % Calculate the potential 
    V_dG1 = v_dG1 * S_dG1;
    V_dG2 = v_dG2 * S_dG2;
    V_dG3 = v_dG3 * S_dG3;

    % Get the Hamiltonian for the different k-vectors
    H1 = getHeye(k(1,:), G1) + V_dG1;
    H2 = getHeye(k(2,:), G2) + V_dG2;
    H3 = getHeye(k(3,:), G3) + V_dG3;

    % Find eigenvalue-matrices and eigenvectors
    [eigVecs1, e1] = eig(H1);
    [eigVecs2, e2] = eig(H2);
    [eigVecs3, e3] = eig(H3);

    % Take out the eigenvalues
    eigs1 = real(diag(e1));
    eigs2 = real(diag(e2));
    eigs3 = real(diag(e3));
    
    % Find index of the minimal eigenvalue
    index1 = find(eigs1 == min(eigs1));
    index2 = find(eigs2 == min(eigs2));
    index3 = find(eigs3 == min(eigs3));
    
    % Get the number of the iteration
    index = (Ecut-EcutInitial)/dE+1;

    % Get the minimal eigenvalue in Hartree energy
    minEig1(index) = eigs1(index1);
    minEig2(index) = eigs2(index2);
    minEig3(index) = eigs3(index3);
        
    % Save the energy cut off
    EnergyCut(index) = Ecut;
    
    % Print the loop-number in the command window
    out = ['Index ', num2str(index), ' out of ', num2str((EcutFinal-EcutInitial)/dE+1)];
    disp(out);
    
end

clf
plot(EnergyCut, minEig1, EnergyCut, minEig2, EnergyCut, minEig3);

%% Task 2: Nice plot of convergence

set(gcf,'renderer','painters','PaperPosition',[0 0 12 7.5]);
plot(EnergyCut, minEig1, EnergyCut, minEig2, EnergyCut, minEig3);
xlim([EcutInitial EcutFinal]);

plotTickLatex2D

title('Convergence in eigenvalues with respect to $E_{cut}$','Interpreter','latex', 'fontsize', 14);
X = xlabel('$E_{cut} [au]$', 'Interpreter','latex', 'fontsize', 12);
Y = ylabel('Energy [au]','Interpreter','latex', 'fontsize', 12);
set(Y, 'Units', 'Normalized', 'Position', [-0.09, 0.5, 0]);
set(X, 'Units', 'Normalized', 'Position', [0.5, -0.06, 0]);

print(gcf,'-depsc2','convergence.eps')

%% Task 2: Plot band structure

clc

% Lattice parameter [au]
a = 5.43/0.529177;   

% Constants
hbar = 1;                   % [au]
me = 1;                     % [au]

% Define how many of the lowest bands that should be plotted
plotNumBands = 8;

% Basis vectors
d(1,:) = a/8*[1 1 1];         
d(2,:) = -a/8*[1 1 1]; 

% Maxvalue in the G-vector
maxValue = 5;

% Define E cut off
Ecut = 15;              % [au]

% The symmetry points in k-space
Gamma = 2*pi/a*[0 0 0];
X = 2*pi/a*[1 0 0];
W = 2*pi/a*[1 0.5 1];
L = 2*pi/a*[0.5 0.5 0.5];
K = 2*pi/a*[0.75 0.75 0];


tickLable = {'$\Gamma$','X','W', 'L','$\Gamma$', 'K', '$\Gamma$'};
%Define the path in k-space
sPoints =[Gamma ;X ; W ;L;Gamma; K ;Gamma]; 

step = 0.001;
% Get the path in k-space
[kMat, tickPoint] = getkMat(sPoints, step);

% Allocate memory
minEig = zeros(length(kMat),1);

for kNum = 1:length(kMat)
   
    
    % Define G-vectors for which the structure factor > zero
    [G] = constructGbig(a, maxValue, kMat(kNum,:), Ecut);

    % Get dG ang form factor for dG
    [dG, v_dG] = getDG(G);

    % Get size of dG
    dGsize = size(dG);
    
    S_dG = zeros(dGsize(1),dGsize(1));

    % Get the structure factor in reciprocal space
    for i = 1:dGsize(1)
        for j = 1:dGsize(1)
            S_dG(i,j) = 2*cos([dG(i,j,1) dG(i,j,2) dG(i,j,3)]*d(1,:)');
        end
    end

    % Get potential in reciprocal space for G_m - G_m'
    V_dG = v_dG * S_dG;

    % Get Hamiltonian
    H = getHeye(kMat(kNum,:), G) + V_dG;

    % Get eigenvalues and eigenvectors
    [eigVecs, eigs] = eig(H);

    eigs = diag(eigs);
    
    for nBands = 1:plotNumBands
        
        % Find index of the minimal eigenvalue
        index = find(eigs == min(eigs));

        % Get the minimal eigenvalue in Hartree energy
        minEig(kNum, nBands) = eigs(index);
        
        % Set the lowest value to something big
        eigs(index) = 100;
    
    end
    
    disp(num2str(kNum));
    
end


set(gcf,'renderer','painters','PaperPosition',[0 0 12 7], 'PaperUnits', 'Centimeters');

plot(minEig);

xlabel('k-vectors');
ylabel('Energy [Hartree]');
set(0, 'defaultTextInterpreter', 'latex');

%% Task 2: Nice plot of band structure

plot(abs(minEig));
set(gca,'XTick',tickPoint)
set(gca,'XTickLabel',tickLable)
set(gca,'XGrid','on')
axis([1 length(kMat) min(min(minEig)) max(max(minEig))])
xlabel('\textbf{k}-vectors','fontsize', 14);
ylabel('Energy [Hartree]');

title('Band structure');

plotTickLatex2D

print(gcf,'-depsc2','bandstruc.eps')

%% Task 3: The valence electron charge density

clc
clear all

w = [1 4 3]/8;

