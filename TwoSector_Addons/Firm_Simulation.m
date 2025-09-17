function [a_sim,z_sim,k_sim,y_sim,n_sim,d_sim,pi_sim,MPK_sim,dist_a] = ...
    Firm_Simulation(y_sol,k_sol,n_sol,pi_sol,k_sol_FB,Agrid,P,psi)


psi_sim      = zeros(P.Na*P.Nz,1);
psicdf_sim   = zeros(P.Na*P.Nz,1);

psi_a_sim    = zeros(P.Na*P.Nz,1);
psi_z_sim    = zeros(P.Na*P.Nz,1);
psi_k_sim    = zeros(P.Na*P.Nz,1);
psi_y_sim    = zeros(P.Na*P.Nz,1);
psi_n_sim    = zeros(P.Na*P.Nz,1);
psi_pi_sim   = zeros(P.Na*P.Nz,1);
psi_k_FB_sim = zeros(P.Na*P.Nz,1);

i = 0;

for z = 1:P.Nz
    for a = 1:P.Na
        i = i + 1;
        psi_sim(i) = psi(a,z);
        if i > 1
            psicdf_sim(i) = psicdf_sim(i-1) + psi_sim(i);
        else
            psicdf_sim(i) = psi_sim(i);
        end
        psi_a_sim(i) = a;
        psi_z_sim(i) = z;
        
        psi_k_sim(i) = k_sol(a,z);
        psi_y_sim(i) = y_sol(a,z);
        psi_n_sim(i) = n_sol(a,z);
        psi_pi_sim(i)   = pi_sol(a,z);
        psi_k_FB_sim(i) = k_sol_FB(a,z);
    end
end




dim_sim  = P.Na*P.Nz;


seed    = linspace(0,1,P.N_sim);
loc     = 1;


z_sim   = zeros(P.N_sim,1);
a_sim   = zeros(P.N_sim,1);
k_sim   = zeros(P.N_sim,1);
y_sim   = zeros(P.N_sim,1);
n_sim   = zeros(P.N_sim,1);
pi_sim  = zeros(P.N_sim,1);
d_sim   = zeros(P.N_sim,1);
k_FB_sim = zeros(P.N_sim,1);

for i = 1:P.N_sim
    
    flag_exit = 0;
    while flag_exit == 0
        if seed(i) <= psicdf_sim(loc) || loc == dim_sim

            a_sim(i)   = psi_a_sim(loc);
            z_sim(i)   = psi_z_sim(loc);

            pi_sim(i)  = psi_pi_sim(loc);
            k_sim(i)   = psi_k_sim(loc);
            y_sim(i)   = psi_y_sim(loc);
            n_sim(i)   = psi_n_sim(loc);
            d_sim(i)   = max(k_sim(i)-Agrid(psi_a_sim(loc)),0);
            k_FB_sim(i) = psi_k_FB_sim(loc);
            flag_exit  = 1;
        else
            loc = loc + 1;
        end
    end
end


MPK_sim = y_sim./k_sim*P.alpha*P.gamma;
dist_a  = k_FB_sim - k_sim;



end

