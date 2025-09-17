function AGG = aggregates_2Sectors(psi, pol_s, Agrid, P, ...
                                   zA_grid, zN_grid, ...
                                   kA, tA, yA, kN, nN, yN, ...
                                   p, w, varrho, R, m_disc)
% Aggregate factors and goods using stationary distribution and occupation map.

Na = numel(Agrid); NzA = numel(zA_grid); NzN = numel(zN_grid);
Nz = NzA * NzN;

kA_full = kron(ones(1,NzN), kA);     % Na x Nz
tA_full = kron(ones(1,NzN), tA);
yA_full = kron(ones(1,NzN), yA);

kN_full = kron(kN, ones(1,NzA));     % Na x Nz
nN_full = kron(nN, ones(1,NzA));
yN_full = kron(yN, ones(1,NzA));

I_W = (pol_s == 1);
I_A = (pol_s == 2);
I_N = (pol_s == 3);

A3 = repmat(Agrid, 1, Nz);
K_supply = sum(sum( psi .* A3 ));
K_demand = sum(sum( psi(I_A) .* kA_full(I_A) )) + sum(sum( psi(I_N) .* kN_full(I_N) ));

N_supply = sum( psi(I_W), 'all' );
N_demand = sum( psi(I_N) .* nN_full(I_N), 'all' );

T_demand = sum( psi(I_A) .* tA_full(I_A), 'all' );

YA = sum( psi(I_A) .* yA_full(I_A), 'all' );
YN = sum( psi(I_N) .* yN_full(I_N), 'all' );

m_pos = max(m_disc, 0);
CA = sum( psi(:) .* ( P.cbarA + (P.psiA/p) * m_pos(:) ) );
CN = sum( psi(:) .* ( P.cbarN + (P.psiN)   * m_pos(:) ) );

AGG = struct();
AGG.K_supply = K_supply;  AGG.K_demand = K_demand;
AGG.N_supply = N_supply;  AGG.N_demand = N_demand;
AGG.T_demand = T_demand;
AGG.YA = YA; AGG.YN = YN; AGG.CA = CA; AGG.CN = CN;
end
