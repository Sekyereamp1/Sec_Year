function [kA, tA, yA, piA] = Firm_Problem_A(zA_grid, col_cons_A, R, varrho, p, P)
% Agriculture: y = z^(1-gA) (k^aA t^(1-aA))^gA, profit = p*y - Rk - varrho t

Na = size(col_cons_A,1); NzA = numel(zA_grid);
kA = zeros(Na, NzA);
tA = zeros(Na, NzA);
yA = zeros(Na, NzA);
piA= zeros(Na, NzA);

aA = P.alphaA; gA = P.gammaA;

for j = 1:NzA
    z = zA_grid(j);
    Az = p * z^(1-gA);

    theta = (1-aA)/aA * (R/varrho);  % t/k if mu=0
    c1 = Az * aA * gA * theta^((1-aA)*gA);
    expo = 1/(gA - 1);
    k_uc = ( R / c1 )^expo;
    t_uc = theta * k_uc;

    for i = 1:Na
        kcap = col_cons_A(i,j);
        if k_uc <= kcap
            k = k_uc; t = t_uc;
        else
            k = kcap;
            e = (1-aA)*gA - 1;
            base = varrho / (Az*(1-aA)*gA * k^(aA*gA));
            if abs(e) < 1e-12
                t = theta * k;
            else
                t = base^(1/e);
            end
        end

        y  = z^(1-gA) * (k^aA * t^(1-aA))^gA;
        pi = p*y - R*k - varrho*t;

        kA(i,j) = k; tA(i,j) = t; yA(i,j) = y; piA(i,j) = pi;
    end
end
end
