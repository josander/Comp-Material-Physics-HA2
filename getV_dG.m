function V_dG = getV_dG(V_r, dG, r)

V_dG = zeros(size(V_r));

for l = 1:length(r)
    for m = 1:length(r)
        for n = 1:length(r)
            V_dG = V_dG + V_r(l,m,n)*exp(-i*dG*[r(l) r(m) r(n)]');
        end
    end
end


end