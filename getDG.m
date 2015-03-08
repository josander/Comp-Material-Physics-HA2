function [dG, v_dG] = getDG(G)

Gsize = size(G);

dG = zeros(Gsize(1), Gsize(1), 3);
v_dG = zeros(Gsize(1), Gsize(1));

for m = 1:Gsize(1)
    for n = m:Gsize(1)
        temp = norm(G(m,:) - G(n,:))^2;

        if round(temp) == 3
            
            dG(m, n,:) = G(m,:) - G(n,:);
            dG(n, m,:) = -dG(m, n, :);
            v_dG(m,n) = -0.056;
            v_dG(n,m) = -0.056;
                
        elseif round(temp) == 8
            
            dG(m, n,:) = G(m,:) - G(n,:);
            dG(n, m,:) = -dG(m, n,:);
            v_dG(m,n) = 0.0138;
            v_dG(n,m) = 0.0138;
                
        elseif round(temp) == 11
            
            dG(m,n,:) = G(m,:) - G(n,:);
            dG(n, m, :) = -dG(m, n, :);
            v_dG(m,n) = 0.0181;
            v_dG(n,m) = 0.0181;
                
        end
    end
end


end