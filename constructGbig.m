function [G] = constructGbig(a, maxValue, kVec, Ecut)

% Constants
hbar = 1;                   % [au]
me = 1;                     % [au]

% Basis vectors in the reciprocal space
b = [-1 1 1; 1 -1 1; 1 1 -1];

% First G-vector
G = [];

% Find different G-vectors where h,k,l all odd or all even
for i = -maxValue:maxValue

   for j = -maxValue:maxValue

        for k = -maxValue:maxValue

            temp = i*b(1,:) + j*b(2,:) + k*b(3,:);

            if mod(temp(1),2) == 0 && mod(temp(2),2) == 0 && mod(temp(3),2) == 0
                
                E = hbar^2/(2*me)*norm(kVec + 2*pi/a.*temp)^2;

                if E < Ecut
                
                    G(end+1,:) = temp;
                
                end

           elseif mod(temp(1),2) == 1 && mod(temp(2),2) == 1 && mod(temp(3),2) == 1

                E = hbar^2/(2*me)*norm(kVec + 2*pi/a.*temp)^2;

                if E < Ecut
                
                    G(end+1,:) = temp;
                
                end

           end

        end

    end

end

% Include the factor
G = 2*pi/a*G;

end