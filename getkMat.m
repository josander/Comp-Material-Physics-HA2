function [ kMat, tickPoint ] = getkMat( sPoints, step )
%GETKMAT Function that construct a path though sPoints 
%   sPoints = n*3 matrix with path-points
%   step = step length of the path

kMat = [0 0 0];
n = 1;

tickPoint = [n];

for i = 2:length(sPoints)
    while norm(sPoints(i,:) - kMat(n,:)) > 1.5*step
        % Unit-vector pointing the direction of the next path point
        wVector = (sPoints(i,:) - kMat(n,:))/norm(sPoints(i,:) - kMat(n,:));    
        kMat = [kMat; kMat(n,:) + wVector*step];
        
        n = n + 1;     
    end
    kMat = [kMat; sPoints(i,:)];
    n = n + 1; 
    % Vector of the indices of the pathpoints in kMat
    tickPoint = [tickPoint (n-1)];
    
end
end

