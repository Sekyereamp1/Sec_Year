function [psi,distance_psi] = Inv_dist_interp(indxg,psi,Pr_z,P)
I = cell(P.Nz);
for z = 1:P.Nz
    I{z} = zeros(P.Na,P.Na);
    for a = 1:P.Na
        
        loc = floor((indxg(a,z)-1)/P.int)+1;
        if loc >= 1 && loc <= P.Na-1
            alpha = (indxg(a,z)-P.int*(floor((indxg(a,z)-1)/P.int))-1)/P.int;
            I{z}(loc,a) = I{z}(loc,a)+(1-alpha);
            I{z}(loc+1,a) = I{z}(loc+1,a)+alpha;
        else
            if loc == 0
                
                I{z}(1,a) = I{z}(1,a)+1;
            else
                I{z}(P.Na,a) = I{z}(P.Na,a)+1;
            end
        end
        
    end
end



distance_psi = 1;
tol_psi      = 10e-9;
ite_psi      = 0;
ite_psi_max  = 10000;



while distance_psi > tol_psi && ite_psi < ite_psi_max
    ite_psi = ite_psi + 1;
    
    tpsi = zeros(P.Na,P.Nz);

    for z = 1:P.Nz
        tpsi = tpsi + I{z}*psi(:,z)*Pr_z(z,:);
    end

    
    distance_psi = max(max(abs(psi-tpsi)));
    psi = tpsi;
    
end


end

