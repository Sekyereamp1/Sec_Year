Two-Sector Add-On (MATLAB)
==========================

Drop these files into the same folder as your existing codebase that includes:
Adjust_Price.m, Inv_dist_interp.m, ScaleTime.m, tauchen.m, VFI_interp.m, Firm_Problem.m, etc.

New files:
 - TwoSector_Main.m        : main driver (price loops, firms, DP, distribution, aggregates)
 - build_z_grids.m         : Tauchen grids for z_A and z_N + Kronecker transition
 - Firm_Problem_A.m        : agriculture firm block (k,t)
 - Firm_Problem_N.m        : non-agriculture firm block (k,n)
 - VFI_OC_SG_interp.m      : DP with occupation choice + Stone–Geary (log)
 - aggregates_2Sectors.m   : aggregation and market clearing
 - Inc_Tax_Func.m          : stub returns zero (compat)

Run:
  addpath(pwd);
  results = TwoSector_Main();

Outputs:
  results.prices  -> [p, w, rhoL, r]
  results.aggregates (K, N, T, YA, YN, CA, CN)
  results.psi, results.pol_s, results.V

Notes:
 - No land-backed collateral; k_s <= phi_s * a only.
 - Stone–Geary feasibility requires m > 0. If you hit many states with m<=0,
   reduce cbarA/cbarN or adjust grids/prices.
