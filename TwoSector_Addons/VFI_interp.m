function [V,indxg,distance_VFI,c_sol,s_sol,it_sol,wt_sol] = ...
    VFI_interp(skip_VFI,inc_sol,r,V,indxg,P,Agrid,Agrid_int,Pr_z)


Util_int = ones(P.Na_int,P.Na,P.Nz);
TV = V;

for z = 1:P.Nz
    for a = 1:P.Na
        temp(:,1,1) = (Agrid(a)*(1+r)+inc_sol(a,z) ...
            - Inc_Tax_Func(Agrid(a)*r+inc_sol(a,z),P)...
            - Agrid_int - P.tau_w*Agrid(a))./(1+P.tau_c);
        Util_int(:,a,z) = ((temp.*(temp>0)+eps).^(1-P.sigma)-1)/(1-P.sigma);
    end
end
Util_int(P.Na_int+1,:,:) = -inf;

distance_VFI = 1;
tol_VFI      = 0.000001;
ite_VFI      = 0;
max_ite_VFI  = 500;
while distance_VFI > tol_VFI && ite_VFI < max_ite_VFI
    ite_VFI = ite_VFI + 1;
    
    
    % Interpolation for the value function

    A_int = linspace(1,P.Na,P.Na_int);
    V_int = ScaleTime(V,A_int);
    
    temp = zeros(P.Na_int,P.Nz);
    for z = 1:P.Nz
        temp(:,z) = P.beta*V_int*Pr_z(z,:)';
    end
    temp(P.Na_int+1,:) = -inf;
    
    for z = 1:P.Nz
        for a = 1:P.Na
            if a == 1
                ap = 1;
            else
                ap = indxg(a-1,z);
            end

            for i = ap:P.Na_int
                U2 = Util_int(i+1,a,z)+temp(i+1,z);
                U1 = Util_int(i,a,z)+temp(i,z);
                if U2 < U1
                    indxg(a,z) = i;
                    TV(a,z)    = U1;
                    break;
                end
            end
        end
    end
    
    distance_VFI = max(max(abs(TV-V)));
    V = TV;
    
    if skip_VFI == 1
        distance_VFI = 0;
        for a = 1:P.Na
            indxg(a,:) = (a-1)*P.int+1;
        end
    end
    
    
    
end

c_sol  = zeros(P.Na,P.Nz);
s_sol  = zeros(P.Na,P.Nz);
it_sol = zeros(P.Na,P.Nz);
wt_sol = zeros(P.Na,P.Nz);
for a = 1:P.Na
    for z = 1:P.Nz
        c_sol(a,z) = ((1+r)*Agrid(a)+inc_sol(a,z)...
            -Inc_Tax_Func(Agrid(a)*r+inc_sol(a,z),P)...
            - Agrid_int(indxg(a,z)) - Agrid(a)*P.tau_w)/(1+P.tau_c);
        s_sol(a,z) = Agrid_int(indxg(a,z))-Agrid(a);
        it_sol(a,z) = Inc_Tax_Func(Agrid(a)*r+inc_sol(a,z),P);
        wt_sol(a,z) = P.tau_w*Agrid(a);
    end
end



end

