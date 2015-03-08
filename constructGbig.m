function [G] = constructGbig(a, maxValue)
    
b = [-1 1 1; 1 -1 1; 1 1 -1];

G(1,:) = [0 0 0];

for i = -maxValue:maxValue
    
    for j = -maxValue:maxValue
        
        for k = -maxValue:maxValue
            
            temp = i*b(1,:) + j*b(2,:) + k*b(3,:);
            
            if mod(temp(1),2) == 0 && mod(temp(2),2) == 0 && mod(temp(3),2) == 0
            
                G(end+1,:) = temp;
                
            elseif mod(temp(1),2) == 1 && mod(temp(2),2) == 1 && mod(temp(3),2) == 1
            
                G(end+1,:) = temp;

            end
            
        end
        
    end
    
end


G = 2*pi/a*G;

end