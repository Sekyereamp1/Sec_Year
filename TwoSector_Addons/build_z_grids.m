function [zA_grid, PrA, zN_grid, PrN, z_pair, PrZ] = build_z_grids(P)
% Build Tauchen grids for z_A and z_N, combine via Kronecker.

[zA, PrA] = tauchen(P.NzA, 0, P.rhoA, P.sigma_zA, 3);
[zN, PrN] = tauchen(P.NzN, 0, P.rhoN, P.sigma_zN, 3);
zA_grid = exp(zA);
zN_grid = exp(zN);

z_pair = zeros(P.NzA*P.NzN, 2);
idx = 0;
for iA = 1:P.NzA
    for iN = 1:P.NzN
        idx = idx + 1;
        z_pair(idx, :) = [zA_grid(iA), zN_grid(iN)];
    end
end

PrZ = kron(PrN, PrA);
end
