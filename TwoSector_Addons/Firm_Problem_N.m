function [kN, nN, yN, piN] = Firm_Problem_N(zN_grid, col_cons_N, R, w, P)
% Non-ag: y = z^(1-gN) (k^aN n^(1-aN))^gN, profit = y - Rk - w n

Na = size(col_cons_N,1); NzN = numel(zN_grid);
kN = zeros(Na, NzN);
nN = zeros(Na, NzN);
yN = zeros(Na, NzN);
piN= zeros(Na, NzN);

aN = P.alphaN; gN = P.gammaN;

for j = 1:NzN
    z = zN_grid(j);
    Bz = z^(1-gN);

    theta = (1-aN)/aN * (R/w);  % n/k if mu=0
    c1 = aN * gN * Bz * theta^((1-aN)*gN);
    expo = 1/(gN - 1);
    k_uc = ( R / c1 )^expo;
    n_uc = theta * k_uc;

    for i = 1:Na
        kcap = col_cons_N(i,j);
        if k_uc <= kcap
            k = k_uc; n = n_uc;
        else
            k = kcap;
            e = (1-aN)*gN - 1;
            base = w / (gN*(1-aN)*Bz * k^(aN*gN));
            if abs(e) < 1e-12
                n = theta * k;
            else
                n = base^(1/e);
            end
        end

        y  = z^(1-gN) * (k^aN * n^(1-aN))^gN;
        pi = y - R*k - w*n;

        kN(i,j) = k; nN(i,j) = n; yN(i,j) = y; piN(i,j) = pi;
    end
end
end
