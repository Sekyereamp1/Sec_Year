function [V, indxg, pol_s, m_disc] = VFI_OC_SG_interp(~, piA, piN, r, w, p, Lambda, V, indxg, P, Agrid, Agrid_int, PrZ)
% DP with occupational choice (Worker=1, Farmer[A]=2, Entrepreneur[N]=3)
% Stone–Geary (log) over (cA,cN). Interpolation over a' grid via ScaleTime.

Na = numel(Agrid); Nz = size(PrZ,1);
V_new  = V;
pol_s  = ones(Na, Nz);
m_disc = zeros(Na, Nz);

tol_VFI = 1e-6; max_ite_VFI = 400; distance = inf; ite = 0;
while (distance > tol_VFI) && (ite < max_ite_VFI)
    ite = ite + 1;

    A_int_idx = linspace(1, Na, P.Na_int);
    V_int = ScaleTime(V, A_int_idx);

    EV = zeros(P.Na_int, Nz);
    for z = 1:Nz
        EV(:,z) = P.beta * V_int * PrZ(z,:)';
    end

    for z = 1:Nz
        for a = 1:Na
            eW = (1+r)*Agrid(a) + w + Lambda;
            eA = (1+r)*Agrid(a) + piA(a,z) + Lambda;
            eN = (1+r)*Agrid(a) + piN(a,z) + Lambda;

            mW = eW - Agrid_int - p*P.cbarA - P.cbarN;
            mA = eA - Agrid_int - p*P.cbarA - P.cbarN;
            mN = eN - Agrid_int - p*P.cbarA - P.cbarN;

            UW = -inf(P.Na_int,1); IA = (mW>0); UW(IA) = log(mW(IA));
            UA = -inf(P.Na_int,1); IB = (mA>0); UA(IB) = log(mA(IB));
            UN = -inf(P.Na_int,1); IC = (mN>0); UN(IC) = log(mN(IC));

            VW_all = UW + EV(:,z);
            VA_all = UA + EV(:,z);
            VN_all = UN + EV(:,z);

            [vW, iw]  = max(VW_all);
            [vA, ia2] = max(VA_all);
            [vN, in2] = max(VN_all);

            [Vbest, choice] = max([vW, vA, vN]);
            V_new(a,z) = Vbest;
            pol_s(a,z) = choice;
            switch choice
                case 1
                    idx = iw;  m_disc(a,z) = mW(idx);
                case 2
                    idx = ia2; m_disc(a,z) = mA(idx);
                case 3
                    idx = in2; m_disc(a,z) = mN(idx);
            end

            % store coarse-grid index for Inv_dist_interp
            if idx <= 1
                indxg(a,z) = 1;
            elseif idx >= P.Na_int
                indxg(a,z) = P.Na_int-1;
            else
                indxg(a,z) = idx;
            end
        end
    end

    distance = max(abs(V_new(:) - V(:)));
    V = V_new;
end

end
