function [G, v_G] = constructG(a, maxValue)
    
b = [-1 1 1; 1 -1 1; 1 1 -1];

G(1,:) = [0 0 0];
v_G(1) = 0;

for i = -maxValue:maxValue  
    for j = -maxValue:maxValue        
        for k = -maxValue:maxValue         
            temp = i*b(1,:) + j*b(2,:) + k*b(3,:);            
            if round(norm(temp)^2) == 3
            
                G(end+1,:) = temp;
                v_G(end+1) = -0.056;
                
            elseif round(norm(temp)^2) == 8
            
                G(end+1,:) = temp;
                v_G(end+1) = 0.0138;
                
            elseif round(norm(temp)^2) == 11
            
                G(end+1,:) = temp;
                v_G(end+1) = 0.0181;
                
            end
            
        end
        
    end
    
end


G = 2*pi/a*G;

end