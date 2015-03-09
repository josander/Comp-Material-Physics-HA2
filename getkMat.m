function [ kMat, tickPoint ] = getkMat( sPoints, step )
%GETKMAT Summary of this function goes here
%   Detailed explanation goes here


kMat = [0 0 0];
n = 1;

tickPoint = [n];

for i = 2:length(sPoints)
    
    while norm(sPoints(i,:) - kMat(n,:)) > 5*step
        wVector = (sPoints(i,:) - kMat(n,:))/norm(sPoints(i,:) - kMat(n,:));    


        kMat = [kMat; kMat(n,:) + wVector*step];
        
        n = n + 1;  
   
    end
    tickPoint = [tickPoint (n-1)];
    
end


end

