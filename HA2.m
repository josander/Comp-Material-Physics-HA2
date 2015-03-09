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
set(get(hcb,'Title'),'String','A Title')
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
a = 1*5.43/0.529177;   

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
maxValue = 5;

% Define E cut off
EcutInitial = 2;
EcutFinal = 15;              % [au]
dE = 0.1;

for Ecut = EcutInitial:dE:EcutFinal
    
    % Define G-vectors for which the structure factor > zero
    [G1] = constructGbig(a, maxValue, k(1,:), Ecut);
    [G2] = constructGbig(a, maxValue, k(2,:), Ecut);
    [G3] = constructGbig(a, maxValue, k(3,:), Ecut);

    % Get the Hamiltonian for the different k-vectors
    H1 = getH(a, k(1,:), G1);
    H2 = getH(a, k(1,:), G2);
    H3 = getH(a, k(1,:), G3);

    % Find eigenvalue-matrices and eigenvectors
    [eigVecs1, e1] = eig(H1);
    [eigVecs2, e2] = eig(H2);
    [eigVecs3, e3] = eig(H3);
    %[eigVecs1, e1] = eigs(H1, 5,'sr');
    %[eigVecs2, e2] = eigs(H2, 5,'sr');
    %[eigVecs3, e3] = eigs(H3, 5,'sr');

    
    % Take out the eigenvalues
    eigs1 = real(diag(e1));
    eigs2 = real(diag(e2));
    eigs3 = real(diag(e3));

    % Find index of the minimal eigenvalue
    index1 = find(eigs1 == min(eigs1));
    index2 = find(eigs2 == min(eigs2));
    index3 = find(eigs3 == min(eigs3));
    
    % Get the number of the iteration
    index = round((Ecut-EcutInitial)/dE+1);

    % Get the minimal eigenvalue in Hartree energy
    minEig1(index) = eigs1(index1(1));
    minEig2(index) = eigs2(index2(1));
    minEig3(index) = eigs3(index3(1));
    
    eigs1(index1) = 100;
    eigs2(index2) = 100;
    eigs3(index3) = 100;
    
    index1 = find(eigs1 == min(eigs1));
    index2 = find(eigs2 == min(eigs2));
    index3 = find(eigs3 == min(eigs3));
    
    
    min2Eig1(index) = eigs1(index1(1));
    min2Eig2(index) = eigs2(index2(1));
    min2Eig3(index) = eigs3(index3(1));
    
        
    
    % Save the energy cut off
    EnergyCut(index) = Ecut;
    
    % Print the loop-number in the command window
    out = ['Index ', num2str(index), ' out of ', num2str((EcutFinal-EcutInitial)/dE+1)];
    disp(out);
    
end

clf
plot(EnergyCut, minEig1, EnergyCut, minEig2, EnergyCut, minEig3);

%% Task 2: Nice plot of convergence


set(gcf,'renderer','painters','PaperPosition',[0 0 12 7], 'PaperUnits', 'Centimeters');
plot(EnergyCut, minEig1, EnergyCut, minEig2, EnergyCut, minEig3);
xlim([EcutInitial EcutFinal]);

plotTickLatex2D

title('Convergence in eigenvalues with respect to $E_{cut}$','Interpreter','latex', 'fontsize', 14);
X = xlabel('$E_{cut} [au]$', 'Interpreter','latex', 'fontsize', 12);
Y = ylabel('Energy [au]','Interpreter','latex', 'fontsize', 12);
set(Y, 'Units', 'Normalized', 'Position', [-0.12, 0.5, 0]);
set(X, 'Units', 'Normalized', 'Position', [0.5, -0.05, 0]);

print(gcf,'-depsc2','convergence.eps')

%% Task 2: Plot band structure

clc

% Lattice parameter [au]
a = 5.43/0.529177;   

% Constants
hbar = 1;                   % [au]
me = 1;                     % [au]

% Define how many of the lowest bands that should be plotted
plotNumBands = 6;

% Basis vectors
d(1,:) = a/8*[1 1 1];         
d(2,:) = -a/8*[1 1 1]; 

% Maxvalue in the G-vector
maxValue = 4;

% Define E cut off
Ecut = 10;              % [au]

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

  
    % Get Hamiltonian
    H = getH(a, kMat(kNum,:), G1);


    % Get eigenvalues and eigenvectors
    [eigVecs, eigs] = eig(H);

    eigs = diag(eigs);
    
    for nBands = 1:plotNumBands
        
        % Find index of the minimal eigenvalue
        index = find(eigs == min(eigs));

        % Get the minimal eigenvalue in Hartree energy
        minEig(kNum, nBands) = eigs(index(1));
        
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
X = xlabel('\textbf{k}-vectors','fontsize', 12);
Y = ylabel('Energy [Hartree]');

title('Band structure','fontsize', 14);

plotTickLatex2D
set(Y, 'Units', 'Normalized', 'Position', [-0.09, 0.5, 0]);
set(X, 'Units', 'Normalized', 'Position', [0.5, -0.05, 0]);

print(gcf,'-depsc2','bandstruc.eps')

save('band.mat','minEig', 'kMat', 'minEig3')
%% Task 3: The valence electron charge density

clc
clear all

w = [1 4 3]/8;

