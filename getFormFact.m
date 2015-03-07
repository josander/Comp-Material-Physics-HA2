function [v_G] = getFormFact(G, a)

    sizeG = size(G);

    G = G*a/(2*pi);

    G = abs(G)

    v_G = zeros(sizeG(1));

    for i = 1:sizeG(1)

       if norm(G(i,:)) == sqrt(3)
           v_G(i) = -0.056;
       elseif norm(G(i,:)) == sqrt(8)
           v_G(i) = 0.0138;
       elseif norm(G(i,:)) == sqrt(11)
           v_G(i) = 0.0181;
       end

    end

end