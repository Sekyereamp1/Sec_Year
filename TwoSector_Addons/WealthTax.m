%% Notes

% This code calculates the effect of taxing assets

clear
clc
tic



%% Step 1: Define the parameters


P.alpha   = 0.4714; % Capital share
P.delta   = 0.06;   % capital depreciation rate
P.gamma   = 0.7;    % span-of-control parameter
P.sigma   = 2;      % CRRA parameter
P.beta    = 0.920280; % Discount factor
P.Nw      = 5.433; % 15.5% of individuals are entrepreneurs


% Ability distribution
P.rho     = 0.85; % Hazard rate
P.sigma_z = 0.75; % std dev of log-normal productivity distribution
P.Nz      = 21;     % Grid points on z (entrepreneur's ability)
[zgrid,Pr_z] = tauchen(P.Nz,0,P.rho,P.sigma_z,3);
zgrid = exp(zgrid);


% Asset holdings
P.Na      = 201; % Grid points on a (wealth)
P.int     = 10;  % Interpolation grid numbers
P.Na_int  = (P.Na-1)*P.int + 1;
P.Aupper  = 2500;

Agrid     = linspace(log(1.1),log(P.Aupper+1),P.Na)';
Agrid     = exp(Agrid)-1;


% Taxation
P.tau_c   = 0.057;
P.tau_ic  = 0.151;
P.tau_il  = 0.841;
P.tau_w   = 0;

% P.tau_ic  = 0.000;
% P.tau_il  = 1;
% P.tau_w   = 0.062;
% 
% P.tau_ic  = 0.000;
% P.tau_il  = 0.80;
% P.tau_w   = 0.000;
% 
% P.tau_ic  = 0.000;
% P.tau_il  = 1.000;
% P.tau_w   = 0.036;
% 
% P.tau_c   = 0.057;
% P.tau_ic  = 0.000;
% P.tau_il  = 1.000;
% P.tau_w   = 0.000;


Agrid_int = ones(P.Na_int,1);
for i = 1:P.Na-1
    for j = 1:P.int
        if j == 1
            Agrid_int((i-1)*P.int+j) = Agrid(i);
        else
            Agrid_int((i-1)*P.int+j) = Agrid(i) + (Agrid(i+1)-Agrid(i))/P.int*(j-1);
        end
    end
end
Agrid_int((P.Na-1)*P.int+1) = Agrid(P.Na);


% Frictions

P.phi     = 2.009; % Financial friction
% P.phi     = 1000;
mf        = 0.00; % Micro-finance: to solve the problem of zero-asset entrepreneurs
col_cons  = Agrid*ones(1,P.Nz)*P.phi+mf;



% !!!!! Initial guess of r

% Initialize some functions


V     = ones(P.Na,P.Nz);
indxg = ones(P.Na,P.Nz);
psi   = ones(P.Na,P.Nz)/(P.Na*P.Nz);



r          = 0.05;
rh         = 0.075;
rl         = -0.01;
distance_r = 1;
tol_r      = 0.001;
ite_r      = 0;
max_ite_r  = 20;


% betah = 0.97;
% betal = 0.86;

while (abs(distance_r) > tol_r) && (ite_r <= max_ite_r)
    
    ite_r = ite_r + 1;
    
    % !!!!! Initial guess of w
    w          = 0.2;
    wh         = 2;
    wl         = 0.01;
    tol_w      = 0.0001;
    distance_w = 1;
    ite_w      = 0;
    max_ite_w  = 25;
    
    
    %% Step 2: Solve the firm's PMP and occupational choice problem
    



    while (abs(distance_w) > tol_w) && (ite_w <= max_ite_w)
        
        error_VFI = 0;
        error_psi = 0;

        ite_w = ite_w + 1;
        
        
        [k_sol,n_sol,y_sol,pi_sol,A_sol,k_sol_FB] = Firm_Problem(zgrid,col_cons,r,w,P);
        
        inc_sol = pi_sol;
        
        inc_sol_int = ones(P.Na_int,P.Nz);
        for z = 1:P.Nz
            inc_sol_int(:,z) = interp1(Agrid,inc_sol(:,z),Agrid_int);
        end
        
 
        %% Step 3: Solve the Dynamic Programming Problem for Entrepreneurs
        
        
        [V,indxg,distance_VFI,c_sol,s_sol,it_sol,wt_sol] = VFI_interp(0,inc_sol,r,V,indxg,P,Agrid,Agrid_int,Pr_z);
        if distance_VFI > 0.000001
            error_VFI = 1;
        end
        
        %% Step 4: Solve for the invariant distribution.

        [psi,distance_psi] = Inv_dist_interp(indxg,psi,Pr_z,P);
        if distance_psi > 10e-9
            error_psi = 1;
        end

        
        %% Step 5: Solve the market clearing condition and update prices
        
        
        AD_n = sum(sum(n_sol.*psi));
        AS_n = P.Nw; % No occupational choice
        distance_w = (AD_n-AS_n);
        
        [ w,wh,wl ] = Adjust_Price( AD_n,AS_n,w,wh,wl,tol_w );


    end
    
    AD_k = sum(sum(k_sol.*psi));
    AS_k = 0;
    for a = 1:P.Na
        for z = 1:P.Nz
            AS_k = AS_k + Agrid(a)*psi(a,z);
        end
    end
    
    distance_r = AD_k - AS_k;
    



    fprintf('======================================================== \n');
    fprintf('error_VFI, %1i, error_psi, %1i  \n',error_VFI,error_psi);
    fprintf('Distance N: %5.4f, w: %5.4f \n',distance_w,w);
    fprintf('Distance K: %5.4f, r: %5.4f, Iteration: %2i \n',distance_r,r,ite_r);
    

    % distance_r = 0;
    
    % [ P.beta,betah,betal ] = Adjust_Price( AD_k,AS_k,P.beta,betah,betal,tol_r );
    % fprintf('    New beta: %7.6f \n',P.beta);

    [ r,rh,rl ] = Adjust_Price( AD_k,AS_k,r,rh,rl,tol_r );
    fprintf('    New r: %5.4f \n',r);

    
    
end

%% Step 6: Calculate the statistics for calibration



% Capital-output ratio

AS_y = sum(sum(y_sol.*psi));


G    = sum(sum(it_sol.*psi)) + AD_k.*P.tau_w + sum(sum(c_sol.*P.tau_c.*psi));

GDP_i = sum(sum(pi_sol.*psi)) + w*sum(sum(psi.*n_sol)) + (r+P.delta)*AD_k;
GDP_o = sum(sum(y_sol.*psi));
GDP_e = sum(sum((c_sol+s_sol).*psi)) + w*sum(sum(psi.*n_sol))...
    + P.delta*sum(sum(k_sol.*psi)) + G;



K    = AD_k;


K_share  = AD_k*(r+P.delta)/GDP_i;
L_share  = w*P.Nw/GDP_i;
pi_share = sum(sum(pi_sol.*psi))/GDP_i;

% Taxation

it_share = sum(sum(it_sol.*psi))/GDP_i;
wt_share = sum(sum(wt_sol.*psi))/GDP_i;
ct_share = sum(sum(c_sol.*P.tau_c.*psi))/GDP_i;


% Simulate some firms to calculate model moments

P.N_sim  = 1000000;

[a_sim,z_sim,k_sim,y_sim,n_sim,d_sim,pi_sim,MPK_sim,dist_a] = ...
    Firm_Simulation(y_sol,k_sol,n_sol,pi_sol,k_sol_FB,Agrid,P,psi);


disp_MPK = std(log(MPK_sim));


z_sim_new = y_sim./(k_sim.^P.alpha.*n_sim.^(1-P.alpha)).^P.gamma;


dummy_fin = (dist_a>0);

dkratio = d_sim./k_sim;
X = [z_sim_new,y_sim,dkratio,dummy_fin,Agrid(a_sim),k_sim,n_sim];
X = sortrows(X,1);

quantile_size = round(size(X,1)*0.25);

for i = 1:3
    dkratio_group(i) = mean(X((i-1)*quantile_size+1:i*quantile_size,3));
    fin_group(i)     = mean(X((i-1)*quantile_size+1:i*quantile_size,4));
    asset_group(i)   = mean(X((i-1)*quantile_size+1:i*quantile_size,5));
    capital_group(i) = mean(X((i-1)*quantile_size+1:i*quantile_size,6));
end

dkratio_group(4) = mean(X(3*quantile_size+1:size(X,1),3));
fin_group(4)     = mean(X(3*quantile_size+1:size(X,1),4));
asset_group(4)   = mean(X(3*quantile_size+1:size(X,1),5));
capital_group(4) = mean(X(3*quantile_size+1:size(X,1),6));

Y_var = Occu_Separate((z_sim_new>prctile(z_sim_new,75)),z_sim_new,1);
X_var = Occu_Separate((z_sim_new>prctile(z_sim_new,75)),k_sim,1);
rankr = corr(Y_var,X_var,'type','spearman');



TFP = AS_y/(AD_k^0.33*AD_n^(1-0.33));
LP  = AS_y/AD_n;



fprintf('======================================================== \n');
fprintf('Equilibrium Model Moments: \n');
fprintf('Capital Income Share:        %+5.3f \n', K_share);
fprintf('Capital-Output Ratio:        %+5.3f \n', AD_k/AS_y);
fprintf('Std Log Output:              %+5.3f \n', std(log(y_sim)));
fprintf('Aggregate Debt-Output Ratio: %+5.3f \n', sum(d_sim)/sum(y_sim));
fprintf('GDP:                         %+5.3f \n', GDP_e);
fprintf('Gov size:                    %+5.3f \n', G);
fprintf('Gov size rel to GDP:         %+5.3f \n', G/GDP_e);
fprintf('Aggregate TFP:               %+5.3f \n',TFP);
fprintf('Aggregate Labor Prod:        %+5.3f \n',LP);
fprintf('======================================================== \n');
fprintf('Percentage of Financially Constrained: \n');
fprintf('    Q1:   %5.3f \n',fin_group(1));
fprintf('    Q2:   %5.3f \n',fin_group(2));
fprintf('    Q3:   %5.3f \n',fin_group(3));
fprintf('    Q4:   %5.3f \n',fin_group(4));
fprintf('======================================================== \n');
toc

