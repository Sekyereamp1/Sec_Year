function [k_sol,n_sol,y_sol,pi_sol,A_sol,k_sol_FB] = Firm_Problem(zgrid,col_cons,r,w,P)

% This code solves for the firm's profit maximization problem.

y_sol = zeros(P.Na,P.Nz);
k_sol = zeros(P.Na,P.Nz);
n_sol = zeros(P.Na,P.Nz);
A_sol = zeros(P.Na,P.Nz);
pi_sol = zeros(P.Na,P.Nz);
k_sol_FB = zeros(P.Na,P.Nz);


% Part 1: the case without financial constraint

A_FB = zgrid;

k_FB = A_FB.*P.gamma^(1/(1-P.gamma))*(P.alpha/(r+P.delta))^((1-P.gamma*(1-P.alpha))/(1-P.gamma))...
    *((1-P.alpha)/w)^((1-P.alpha)*P.gamma/(1-P.gamma));
n_FB = A_FB.*P.gamma^(1/(1-P.gamma))*(P.alpha/(r+P.delta))^((P.alpha*P.gamma)/(1-P.gamma))...
    *((1-P.alpha)/w)^((1-P.alpha*P.gamma)/(1-P.gamma));
y_FB = A_FB.*P.gamma^(P.gamma/(1-P.gamma))*(P.alpha/(r+P.delta))^(P.alpha*P.gamma/(1-P.gamma))...
    *((1-P.alpha)/w)^((1-P.alpha)*P.gamma/(1-P.gamma));

% y_FB_2 = A_FB.^(1-P.gamma).*k_FB.^(P.alpha*P.gamma).*n_FB.^((1-P.alpha)*P.gamma);


pi_FB = (1-P.gamma)*y_FB;

% (y_FB - w*n_FB - (r+P.delta)*k_FB)./A_FB



% Part 2: the case with financial constraint

for z = 1:P.Nz
    for a = 1:P.Na
        
        k_sol_FB(a,z) = k_FB(z);
        if k_FB(z) < col_cons(a,z)
            y_sol(a,z) = y_FB(z);
            k_sol(a,z) = k_FB(z);
            n_sol(a,z) = n_FB(z);
            pi_sol(a,z) = pi_FB(z);
            A_sol(a,z) = A_FB(z);
        else
            k_sol(a,z) = col_cons(a,z);
            A_sol(a,z) = zgrid(z);
            n_sol(a,z) = ((1-P.alpha)*P.gamma/w)^(1/(1-(1-P.alpha)*P.gamma))...
                *A_sol(a,z)^((1-P.gamma)/(1-(1-P.alpha)*P.gamma))...
                *k_sol(a,z)^(P.alpha*P.gamma/(1-(1-P.alpha)*P.gamma));
            y_sol(a,z) = A_sol(a,z)^(1-P.gamma)*k_sol(a,z)^(P.alpha*P.gamma)*n_sol(a,z)^((1-P.alpha)*P.gamma);
            pi_sol(a,z) = y_sol(a,z)-(r+P.delta)*k_sol(a,z)-w*n_sol(a,z);
        end
    end
end
            


end

