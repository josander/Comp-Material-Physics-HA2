function [G] = getG(a)

    G(1,:) = [0 0 0];
    
    G(end+1,:) = [1 1 1];

    v = [-1 1 1];
    P = perms(v);
    sizeP = size(P);
    G(end+1:end+sizeP(1),:) = P;

    v = [-1 -1 1];
    P = perms(v);
    sizeP = size(P);
    G(end+1:end+sizeP(1),:) = P;

    v = [-1 -1 -1];
    P = perms(v);
    sizeP = size(P);
    G(end+1:end+sizeP(1),:) = P;

    v = [0 2 2];
    P = perms(v);
    sizeP = size(P);
    G(end+1:end+sizeP(1),:) = P;

    v = [0 2 -2];
    P = perms(v);
    sizeP = size(P);
    G(end+1:end+sizeP(1),:) = P;

    v = [0 -2 -2];
    P = perms(v);
    sizeP = size(P);
    G(end+1:end+sizeP(1),:) = P;

    v = [1 1 3];
    P = perms(v);
    sizeP = size(P);
    G(end+1:end+sizeP(1),:) = P;

    v = [-1 1 3];
    P = perms(v);
    sizeP = size(P);
    G(end+1:end+sizeP(1),:) = P;

    v = [-1 -1 3];
    P = perms(v);
    sizeP = size(P);
    G(end+1:end+sizeP(1),:) = P;

    v = [-1 1 -3];
    P = perms(v);
    sizeP = size(P);
    G(end+1:end+sizeP(1),:) = P;

    v = [-1 -1 -3];
    P = perms(v);
    sizeP = size(P);
    G(end+1:end+sizeP(1),:) = P;


    G = 2*pi/a*G;

end