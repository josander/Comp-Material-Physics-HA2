function [v_G] = getFormFact(G, a)

    sizeG = size(G);

    G = G*a/(2*pi);

    G = abs(G)

    v_G = zeros(1,sizeG(1));

    for i = 1:sizeG(1)

       if round(norm(G(i,:))^2) == 3
           v_G(i) = -0.056;
       elseif round(norm(G(i,:))^2) == 8
           v_G(i) = 0.0138;
       elseif round(norm(G(i,:))^2) == 11
           v_G(i) = 0.0181;
       end

    end

end