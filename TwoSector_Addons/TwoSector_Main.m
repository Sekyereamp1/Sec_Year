function results = TwoSector_Main()
% TW0SECTOR_MAIN
% -------------------------------------------------------------------------
% Two-sector extension based on the structure/logic of the uploaded code:
% - Price search via Adjust_Price
% - Firm problems solved on (a,z) grids
% - DP with interpolation and occupational choice (Worker, Farmer[A], Entrepreneur[N])
% - Invariant distribution via Inv_dist_interp
% - Stone–Geary (log) utility over two goods (A and N)
% - Agriculture uses (k, t) with land market; Non-ag uses (k, n) with labor
% - No land-backed collateral (constraints: k_s <= phi_s * a)
% -------------------------------------------------------------------------

clc

%% PARAMETERS (follow your previous style; only add what's needed)
P = struct();

% --- From WealthTax-style baseline (kept) ---
P.alpha   = 0.4714;   % baseline capital share (used as fallback)
P.delta   = 0.06;     % depreciation
P.gamma   = 0.7;      % span-of-control
P.sigma   = 2;        % CRRA in older code (not used here, kept for compat)
P.beta    = 0.920280; % discount

% Ability processes (Tauchen grids for zA and zN)
P.rhoA     = 0.85;    % persistence ag
P.sigma_zA = 0.75;    % shock std (log space)
P.NzA      = 11;      % grid size ag
P.rhoN     = 0.85;    % persistence non-ag
P.sigma_zN = 0.75;    % shock std (log space)
P.NzN      = 11;      % grid size non-ag

% Asset grid consistent with earlier files
P.Na      = 201;
P.int     = 10;
P.Na_int  = (P.Na-1)*P.int + 1;
P.Aupper  = 2500;

Agrid     = linspace(log(1.1),log(P.Aupper+1),P.Na)';
Agrid     = exp(Agrid)-1;

% Build fine asset grid for interpolation (as in earlier files)
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

% Financial frictions (sector-specific; fallback to original phi if absent)
P.phi  = 2.009;
P.phiA = P.phi;
P.phiN = P.phi;

% Land endowment
P.Tbar = 1.0;

% Stone–Geary preferences (new)
P.psiA  = 0.35;
P.psiN  = 1 - P.psiA;
P.cbarA = 0.05;
P.cbarN = 0.05;

% Sector-specific production parameters (fallbacks to baseline alpha,gamma)
P.alphaA = P.alpha;
P.gammaA = P.gamma;
P.alphaN = P.alpha;
P.gammaN = P.gamma;

% INITIAL PRICES (same bracket-search logic as your files)
p      = 1.0;  ph      = 4.0;  pl      = 0.05;   % ag good relative price
w      = 0.8;  wh      = 4.0;  wl      = 0.01;   % wage
rhoL   = 1.0;  rhoH    = 6.0;  rhol    = 0.01;   % land rent (varrho)
r      = 0.05; rh      = 0.10; rl      = -0.01;  % bond rate
tol_p  = 1e-3; tol_w   = 1e-4; tol_rho = 1e-3; tol_r = 1e-3;

% Storage
V     = zeros(P.Na, P.NzA*P.NzN);
indxg = ones(P.Na, P.NzA*P.NzN);
psi   = ones(P.Na, P.NzA*P.NzN) / (P.Na * P.NzA * P.NzN);

% Build (zA,zN) grid and transition
[zA_grid, PrA, zN_grid, PrN, z_pair, PrZ] = build_z_grids(P);

%% Outer loop: interest rate
distance_r = inf; ite_r = 0; max_ite_r = 30;
while (abs(distance_r) > tol_r) && (ite_r <= max_ite_r)
    ite_r = ite_r + 1;
    R = r + P.delta;

    %% Inner loop: wage, land rent, relative price
    distance_w = inf; distance_rho = inf; distance_p = inf;
    ite_w = 0; max_ite_w = 30;

    while ( (abs(distance_w) > tol_w) || (abs(distance_rho) > tol_rho) || (abs(distance_p) > tol_p) ) ...
            && (ite_w <= max_ite_w)

        ite_w = ite_w + 1;

        % --- 1) Firms by sector (arrays on (a,z)) ---
        col_cons_A = (Agrid * ones(1,P.NzA)) * P.phiA;
        col_cons_N = (Agrid * ones(1,P.NzN)) * P.phiN;

        [kA, tA, yA, piA] = Firm_Problem_A(zA_grid, col_cons_A, R, rhoL, p, P);
        [kN, nN, yN, piN] = Firm_Problem_N(zN_grid, col_cons_N, R, w, P);

        % --- 2) DP with occupation choice + Stone–Geary ---
        piA_full = kron(ones(1,P.NzN), piA);     % Na x (NzA*NzN)
        piN_full = kron(piN,           ones(1,P.NzA));

        Lambda = rhoL * P.Tbar;  % land rent rebate per capita

        [V, indxg, pol_s, m_disc] = VFI_OC_SG_interp( ...
            0, piA_full, piN_full, r, w, p, Lambda, V, indxg, P, Agrid, Agrid_int, PrZ);

        % --- 3) Invariant distribution on (a,zA,zN) ---
        [psi, ~] = Inv_dist_interp(indxg, psi, PrZ, struct('Nz', P.NzA*P.NzN, 'Na', P.Na, 'int', P.int));

        % --- 4) Aggregation and market clearing ---
        AGG = aggregates_2Sectors(psi, pol_s, Agrid, P, ...
                zA_grid, zN_grid, ...
                kA, tA, yA, kN, nN, yN, ...
                p, w, rhoL, R, m_disc);

        % Residuals
        exK = AGG.K_demand - AGG.K_supply;
        exN = AGG.N_demand - AGG.N_supply;
        exT = AGG.T_demand - P.Tbar;
        exG = (p*AGG.YA + AGG.YN) - (p*AGG.CA + AGG.CN + P.delta*AGG.K_demand);

        fprintf(['it_r %02d | it_in %02d | exK=%+.3e exN=%+.3e exT=%+.3e exG=%+.3e | ' ...
                 'p=%.3f w=%.3f rhoL=%.3f r=%.3f\n'], ...
                 ite_r, ite_w, exK, exN, exT, exG, p, w, rhoL, r);

        if max(abs([exK, exN, exT, exG])) < 5e-3
            break;
        end

        % --- 5) Adjust prices using your bracket logic ---
        [w,  wh,  wl]  = Adjust_Price( AGG.N_demand, AGG.N_supply,  w,  wh,  wl,  tol_w );
        [rhoL, rhoH, rhol] = Adjust_Price( AGG.T_demand, P.Tbar,       rhoL, rhoH, rhol, tol_rho );

        % Relative price p by value-clearing residual
        if exG > 0
            pl = p; p = 0.5*p + 0.5*ph;
        elseif exG < 0
            ph = p; p = 0.5*p + 0.5*pl;
        end

        distance_w   = AGG.N_demand - AGG.N_supply;
        distance_rho = AGG.T_demand - P.Tbar;
        distance_p   = exG;
    end

    % Capital market residual (outer loop)
    distance_r = exK;
    [r, rh, rl] = Adjust_Price( AGG.K_demand, AGG.K_supply, r, rh, rl, tol_r );
end

results = struct();
results.prices = [p, w, rhoL, r];
results.V      = V;
results.psi    = psi;
results.pol_s  = pol_s;
results.aggregates = AGG;
results.params = P;
results.msg    = 'Two-sector equilibrium (approx).';

fprintf('\nFinal prices: p=%.4f, w=%.4f, rhoL=%.4f, r=%.4f\n', p, w, rhoL, r);

end
