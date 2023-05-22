# First, read decomposition
irf_factor_1 = load(data_path_1*"irf_factor.jld","irf_factor")
irf_factor = load(data_path*"irf_factor.jld","irf_factor")
Policy_S_0  = load(data_path*"Policies_S_0.jld","other_variables")
State_S_0 = load(data_path*"States_S_0.jld","state_variables") 
S_0 = [0; 0; 0; State_S_0; Policy_S_0]
Policy_S_Star  = load(data_path*"Policies_S_Star.jld","other_variables")
State_S_Star = load(data_path*"States_S_Star.jld","state_variables") 
S_Star = [0; 0; 0; State_S_Star; Policy_S_Star];
Policy_S_1  = load(data_path*"Policies_S_1.jld","other_variables")
State_S_1 = load(data_path*"States_S_1.jld","state_variables") 
S_1 = [0; 0; 0; State_S_1; Policy_S_1];
CF_Policies_A = load(data_path*"Counterfactual_Policies_A.jld","cf_policy_mat")
CF_Policies_B = load(data_path*"Counterfactual_Policies_B.jld","cf_policy_mat")
# Objects needed: 
# [1] State S0, S1, S star [three objects]
# [2] counterfactual policies, CF_Policies_A, CF_Policies_B [two objects]
# Will write these in once we have data filtered

# The lines above [reading in objects] correspond to lines 1-77 of calc_decomposition.m

# The next starts from line 78 of the Matlab file


## 2. Implement Decomposition
# -------------------------------------------------- #
# [1] Derive the individual policies at S Star
# -------------------------------------------------- #
sav_a_0 = S_Star[jx_sav_a]; # [S_Star[jx_aggrw] .* S_Star[jx_tht_a]] + S_Star[jx_w] .* S_Star[jx_l] .* prm_vec[param_id].labor_alloc_a - prm_vec[param_id].lmbd_a * S_Star[jx_c_a];
k_a_0   = S_Star[jx_k_a]; #sav_a_0 .* [1-S_Star[jx_sh_a]] ./ S_Star[jx_q];
b_a_0   = S_Star[jx_sav_a] - S_Star[jx_k_a]*S_Star[jx_q]; #sav_a_0 .* [S_Star[jx_sh_a]];
c_a_0   = prm_vec[param_id].lmbd_a * S_Star[jx_c_a];

sav_b_0 = S_Star[jx_sav_b]; # [S_Star[jx_aggrw] .* S_Star[jx_tht_a]] + S_Star[jx_w] .* S_Star[jx_l] .* prm_vec[param_id].labor_alloc_a - prm_vec[param_id].lmbd_a * S_Star[jx_c_a];
k_b_0   = S_Star[jx_k_b]; #sav_a_0 .* [1-S_Star[jx_sh_a]] ./ S_Star[jx_q];
b_b_0   = S_Star[jx_sav_b] - S_Star[jx_k_b]*S_Star[jx_q]; #sav_a_0 .* [S_Star[jx_sh_a]];
c_b_0   = prm_vec[param_id].lmbd_b * S_Star[jx_c_b];

sav_c_0 = S_Star[jx_sav_c]; # [S_Star[jx_aggrw] .* S_Star[jx_tht_a]] + S_Star[jx_w] .* S_Star[jx_l] .* prm_vec[param_id].labor_alloc_a - prm_vec[param_id].lmbd_a * S_Star[jx_c_a];
k_c_0   = S_Star[jx_k_c]; #sav_a_0 .* [1-S_Star[jx_sh_a]] ./ S_Star[jx_q];
b_c_0   = S_Star[jx_sav_c] - S_Star[jx_k_c]*S_Star[jx_q]; #sav_a_0 .* [S_Star[jx_sh_a]];
c_c_0   = prm_vec[param_id].lmbd_c * S_Star[jx_c_c];



# -------------------------------------------------- #
# [2] Derive the individual policies at S1
# -------------------------------------------------- #
sav_a_1 = S_1[jx_sav_a]; # [S_Star[jx_aggrw] .* S_Star[jx_tht_a]] + S_Star[jx_w] .* S_Star[jx_l] .* prm_vec[param_id].labor_alloc_a - prm_vec[param_id].lmbd_a * S_Star[jx_c_a];
k_a_1   = S_1[jx_k_a]; #sav_a_0 .* [1-S_Star[jx_sh_a]] ./ S_Star[jx_q];
b_a_1   = S_1[jx_sav_a] - S_1[jx_k_a]*S_1[jx_q]; #sav_a_0 .* [S_Star[jx_sh_a]];
c_a_1   = prm_vec[param_id].lmbd_a * S_1[jx_c_a];

sav_b_1 = S_1[jx_sav_b]; # [S_Star[jx_aggrw] .* S_Star[jx_tht_a]] + S_Star[jx_w] .* S_Star[jx_l] .* prm_vec[param_id].labor_alloc_a - prm_vec[param_id].lmbd_a * S_Star[jx_c_a];
k_b_1   = S_1[jx_k_b]; #sav_a_0 .* [1-S_Star[jx_sh_a]] ./ S_Star[jx_q];
b_b_1   = S_1[jx_sav_b] - S_1[jx_k_b]*S_1[jx_q]; #sav_a_0 .* [S_Star[jx_sh_a]];
c_b_1   = prm_vec[param_id].lmbd_b * S_1[jx_c_b];

sav_c_1 = S_1[jx_sav_c]; # [S_Star[jx_aggrw] .* S_Star[jx_tht_a]] + S_Star[jx_w] .* S_Star[jx_l] .* prm_vec[param_id].labor_alloc_a - prm_vec[param_id].lmbd_a * S_Star[jx_c_a];
k_c_1   = S_1[jx_k_c]; #sav_a_0 .* [1-S_Star[jx_sh_a]] ./ S_Star[jx_q];
b_c_1   = S_1[jx_sav_c] - S_1[jx_k_c]*S_1[jx_q]; #sav_a_0 .* [S_Star[jx_sh_a]];
c_c_1   = prm_vec[param_id].lmbd_c * S_1[jx_c_c];

# Here is line 142 of the matlab file

# -------------------------------------------------- #
# [3] Scale Counterfactual Policies
# -------------------------------------------------- #
c_a_CFA = CF_Policies_A[1,1] * prm_vec[param_id].lmbd_a;
c_b_CFA = CF_Policies_A[2,1] * prm_vec[param_id].lmbd_b;
c_c_CFA = CF_Policies_A[3,1] * prm_vec[param_id].lmbd_c;

b_a_CFA = CF_Policies_A[1,3] * prm_vec[param_id].lmbd_a;
b_b_CFA = CF_Policies_A[2,3] * prm_vec[param_id].lmbd_b;
b_c_CFA = CF_Policies_A[3,3] * prm_vec[param_id].lmbd_c;

k_a_CFA = CF_Policies_A[1,4] * prm_vec[param_id].lmbd_a;
k_b_CFA = CF_Policies_A[2,4] * prm_vec[param_id].lmbd_b;
k_c_CFA = CF_Policies_A[3,4] * prm_vec[param_id].lmbd_c;

c_a_CFB = CF_Policies_B[1,1] * prm_vec[param_id].lmbd_a;
c_b_CFB = CF_Policies_B[2,1] * prm_vec[param_id].lmbd_b;
c_c_CFB = CF_Policies_B[3,1] * prm_vec[param_id].lmbd_c;

b_a_CFB = CF_Policies_B[1,3] * prm_vec[param_id].lmbd_a;
b_b_CFB = CF_Policies_B[2,3] * prm_vec[param_id].lmbd_b;
b_c_CFB = CF_Policies_B[3,3] * prm_vec[param_id].lmbd_c;

k_a_CFB = CF_Policies_B[1,4] * prm_vec[param_id].lmbd_a;
k_b_CFB = CF_Policies_B[2,4] * prm_vec[param_id].lmbd_b;
k_c_CFB = CF_Policies_B[3,4] * prm_vec[param_id].lmbd_c;

# here, line 170 of matlab

# -------------------------------------------------- #
# [4] Compute Direct Effects
# 3 by 3 Table, a, b, c
# -------------------------------------------------- #

Direct_Effects_k = [k_a_CFA-k_a_0, k_b_CFA-k_b_0, k_c_CFA-k_c_0];
Direct_Effects_b = [b_a_CFA-b_a_0, b_b_CFA-b_b_0, b_c_CFA-b_c_0];
Direct_Effects_c = [c_a_CFA-c_a_0, c_b_CFA-c_b_0, c_c_CFA-c_c_0];

Return_Effects_k = [k_a_CFA-k_a_CFB, k_b_CFA-k_b_CFB, k_c_CFA-k_c_CFB];
Return_Effects_b = [b_a_CFA-b_a_CFB, b_b_CFA-b_b_CFB, b_c_CFA-b_c_CFB];
Return_Effects_c = [c_a_CFA-c_a_CFB, c_b_CFA-c_b_CFB, c_c_CFA-c_c_CFB];


# -------------------------------------------------- #
# [5] Compute Total Effects
# -------------------------------------------------- #
Total_Effects_k = [k_a_1-k_a_0, k_b_1-k_b_0, k_c_1-k_c_0];
Total_Effects_b = [b_a_1-b_a_0, b_b_1-b_b_0, b_c_1-b_c_0];
Total_Effects_c = [c_a_1-c_a_0, c_b_1-c_b_0, c_c_1-c_c_0];

# -------------------------------------------------- #
# [6] Load MPK/MPB/MPC at S_Star
# -------------------------------------------------- #
MPK_a = S_Star[jx_MPK1] / S_Star[jx_q];
MPK_b = S_Star[jx_MPK2] / S_Star[jx_q];
MPK_c = S_Star[jx_MPK3] / S_Star[jx_q];
MPC_a = S_Star[jx_MPC1];
MPC_b = S_Star[jx_MPC2];
MPC_c = S_Star[jx_MPC3];
MPB_a = 1 - S_Star[jx_MPC1] - S_Star[jx_MPK1];
MPB_b = 1 - S_Star[jx_MPC2] - S_Star[jx_MPK2];
MPB_c = 1 - S_Star[jx_MPC3] - S_Star[jx_MPK3];

MP_k = [MPK_a; MPK_b; MPK_c];
MP_b = [MPB_a; MPB_b; MPB_c];
MP_c = [MPC_a; MPC_b; MPC_c];

# -------------------------------------------------- #
# [7] Load Net Worth Changes
# -------------------------------------------------- #
n_a_0 = S_Star[jx_n_a]; #S_Star[jx_aggrw] .* S_Star[jx_tht_a];
n_b_0 = S_Star[jx_n_b]; #S_Star[jx_aggrw] .* S_Star[jx_tht_b];
n_c_0 = S_Star[jx_n_c]; #S_Star[jx_aggrw] - n_a_0 - n_b_0;


aggr_wealth_1 = (S_0[jx_sav_a] + S_0[jx_sav_b] + S_0[jx_sav_c])*(S_1[jx_pi] + (1-prm_vec[param_id].delta)*S_1[jx_q]);

n_a_1 = ( prm_vec[param_id].s_bar_a*prm_vec[param_id].xi*aggr_wealth_1 + (1-prm_vec[param_id].xi)*S_0[jx_sav_a]*  
((1-S_0[jx_sh_a])*(S_1[jx_pi] + (1-prm_vec[param_id].delta)*S_1[jx_q]) + S_0[jx_sh_a]*S_0[jx_nom_i]/S_1[jx_infl]) ); # S_1[jx_n_a); #S_1[jx_aggrw) .* S_1[jx_tht_a);

n_b_1 = ( (1- prm_vec[param_id].s_bar_c - prm_vec[param_id].s_bar_a)*prm_vec[param_id].xi*aggr_wealth_1 + (1-prm_vec[param_id].xi)*S_0[jx_sav_b]*  
((1-S_0[jx_sh_b])*(S_1[jx_pi] + (1-prm_vec[param_id].delta)*S_1[jx_q]) + S_0[jx_sh_b]*S_0[jx_nom_i]/S_1[jx_infl]) );  # S_1[jx_n_c); #S_1[jx_aggrw) - n_a_1 - n_b_1;

n_c_1 = ( prm_vec[param_id].s_bar_c*prm_vec[param_id].xi*aggr_wealth_1 + (1-prm_vec[param_id].xi)*S_0[jx_sav_c]* 
((1-S_0[jx_sh_c])*(S_1[jx_pi] + (1-prm_vec[param_id].delta)*S_1[jx_q]) + S_0[jx_sh_c]*S_0[jx_nom_i]/S_1[jx_infl]) ); 
# -------------------------------------------------- #
# (8) Indirect Effects and Residuals
# -------------------------------------------------- #
Indirect_Effects_k = [MPK_a*(n_a_1-n_a_0); 
                      MPK_b*(n_b_1-n_b_0); 
                      MPK_c*(n_c_1-n_c_0)];
Indirect_Effects_b = [MPB_a*(n_a_1-n_a_0); 
                      MPB_b*(n_b_1-n_b_0); 
                      MPB_c*(n_c_1-n_c_0)];
Indirect_Effects_c = [MPC_a*(n_a_1-n_a_0); 
                      MPC_b*(n_b_1-n_b_0); 
                      MPC_c*(n_c_1-n_c_0)];

# Note: these are vectors

                  
# -------------------------------------------------- #
# (9) Net Worth Component: Bond
# Given bond positions b_0 in state S_0
# What would be the value in S_Star?
# What would be the value in S_1?
# The difference is the bond term
# -------------------------------------------------- #
# Bond and Capital Positions in S_0
sav_a_orig = S_0[jx_sav_a]; 
k_a_orig   = S_0[jx_k_a]; 
b_a_orig   = S_0[jx_sav_a] - S_0[jx_k_a]*S_0[jx_q]; 
sav_b_orig = S_0[jx_sav_b]; 
k_b_orig   = S_0[jx_k_b]; 
b_b_orig   = S_0[jx_sav_b] - S_0[jx_k_b]*S_0[jx_q]; 
sav_c_orig = S_0[jx_sav_c]; 
k_c_orig   = S_0[jx_k_c]; 
b_c_orig   = S_0[jx_sav_c] - S_0[jx_k_c]*S_0[jx_q]; 

# here corr to line 262 of Matlab

# Would-be bond value in S_Star
bv_a_0 = b_a_orig * S_0[jx_nom_i] / S_Star[jx_infl];
bv_b_0 = b_b_orig * S_0[jx_nom_i] / S_Star[jx_infl];
bv_c_0 = b_c_orig * S_0[jx_nom_i] / S_Star[jx_infl];

# Would-be bond value in S_1
bv_a_1 = b_a_orig * S_0[jx_nom_i] / S_1[jx_infl];
bv_b_1 = b_b_orig * S_0[jx_nom_i] / S_1[jx_infl];
bv_c_1 = b_c_orig * S_0[jx_nom_i] / S_1[jx_infl];

# Bond Components
Bond_Component = [bv_a_1 - bv_a_0, bv_b_1 - bv_b_0, bv_c_1 - bv_c_0];

# -------------------------------------------------- #
# (10) Net Worth Component: Capital, Profit
# -------------------------------------------------- #
Profit_Component = [k_a_orig * (S_1[jx_pi]-S_Star[jx_pi] ); 
                    k_b_orig * (S_1[jx_pi]-S_Star[jx_pi] ); 
                    k_c_orig * (S_1[jx_pi]-S_Star[jx_pi])];

# -------------------------------------------------- #
# (11) Capital Price Component
# -------------------------------------------------- #
Capital_Price_Component =(1-prm_vec[param_id].delta) .* 
    [k_a_orig * (S_1[jx_q]-S_Star[jx_q]);
     k_b_orig * (S_1[jx_q]-S_Star[jx_q]); 
     k_c_orig * (S_1[jx_q]-S_Star[jx_q])] 
# -------------------------------------------------- #
# (12) Transfer Component
# The rest must come from transfer
# No need to compute again
# -------------------------------------------------- #
NetWorth_Change    = [n_a_1 - n_a_0; n_b_1 - n_b_0; n_c_1 - n_c_0];
NetWorth_0         = [n_a_0; n_b_0; n_c_0];
NetWorth_1         = [n_a_1; n_b_1; n_c_1];
Transfer_Component = @. NetWorth_Change - Capital_Price_Component - Profit_Component - Bond_Component;

# -------------------------------------------------- #
# (13) Compute the separate d(wl)/d(veps) term
# -------------------------------------------------- #                
Labor_Income_0 = @. S_Star[jx_w] * S_Star[jx_l]*[prm_vec[param_id].labor_alloc_a*prm_vec[param_id].lmbd_a, 
     prm_vec[param_id].labor_alloc_b*prm_vec[param_id].lmbd_b, 
     prm_vec[param_id].labor_alloc_c*prm_vec[param_id].lmbd_c];
Labor_Income_1 = S_1[jx_w] .* S_1[jx_l] .*[prm_vec[param_id].labor_alloc_a*prm_vec[param_id].lmbd_a, 
     prm_vec[param_id].labor_alloc_b*prm_vec[param_id].lmbd_b, 
     prm_vec[param_id].labor_alloc_c*prm_vec[param_id].lmbd_c];
Labor_Income_Change = S_1[jx_w] * S_1[jx_l] - S_Star[jx_w] * S_Star[jx_l];                
Labor_Income_Component = ( Labor_Income_Change .* 
    [prm_vec[param_id].labor_alloc_a*prm_vec[param_id].lmbd_a, 
     prm_vec[param_id].labor_alloc_b*prm_vec[param_id].lmbd_b, 
     prm_vec[param_id].labor_alloc_c*prm_vec[param_id].lmbd_c] );

# -------------------------------------------------- #
# (14) Adjust Direct and Indirect Effect
# Add Labor_Income_Component * MPR to Indirect Effect
# Take them away from direct effect
# -------------------------------------------------- #  
Direct_Effects_k = @. Direct_Effects_k - Labor_Income_Component * MP_k;
Direct_Effects_b = @. Direct_Effects_b - Labor_Income_Component * MP_b;
Direct_Effects_c = @. Direct_Effects_c - Labor_Income_Component * MP_c;

Indirect_Effects_k = @. Indirect_Effects_k + Labor_Income_Component * MP_k;
Indirect_Effects_b = @. Indirect_Effects_b + Labor_Income_Component * MP_b;
Indirect_Effects_c = @. Indirect_Effects_c + Labor_Income_Component * MP_c;

NetWorth_Component = @. NetWorth_Change + Labor_Income_Component;

# -------------------------------------------------- #
# (15) Compute Residuals
# -------------------------------------------------- #  
Residual_k = @. Total_Effects_k - Direct_Effects_k - Indirect_Effects_k;
Residual_b = @. Total_Effects_b - Direct_Effects_b - Indirect_Effects_b;
Residual_c = @. Total_Effects_c - Direct_Effects_c - Indirect_Effects_c;


Other_Effects_c = @. Direct_Effects_c - Return_Effects_c;
Other_Effects_b = @. Direct_Effects_b - Return_Effects_b;
Other_Effects_k = @. Direct_Effects_k - Return_Effects_k;

## 6. Further Decompose Return Effect on Capital
# -------------------------------------------------- #
# E(rk) Change
# -------------------------------------------------- #  
Erk_change = (S_1[jx_Erk] - S_Star[jx_Erk]);
Erk_component_k = zeros(3)
Erk_component_k[1] = Erk_change * S_Star[jx_dk_dErk1];
Erk_component_k[2] = Erk_change * S_Star[jx_dk_dErk2];
Erk_component_k[3] = Erk_change * S_Star[jx_dk_dErk3];

# -------------------------------------------------- #
# E(rf) Change
# -------------------------------------------------- #  
Erf_change = (S_1[jx_Erf] - S_Star[jx_Erf]);
Erf_component_k = zeros(3)
Erf_component_k[1] = Erf_change * S_Star[jx_dk_dErf1];
Erf_component_k[2] = Erf_change * S_Star[jx_dk_dErf2];
Erf_component_k[3] = Erf_change * S_Star[jx_dk_dErf3];

decomp_mat[:,1,ccc]  = ones(1,3)*irf_factor_1*sum(Total_Effects_k)/sum([k_a_0 ;k_b_0; k_c_0]);
decomp_mat[:,2,ccc]  = irf_factor.*log.(([k_a_0; k_b_0; k_c_0] + Total_Effects_k)./[k_a_0; k_b_0; k_c_0]);
decomp_mat[:,3,ccc]  = @. NetWorth_0./(S_0[[jx_sav_a ;jx_sav_b ;jx_sav_c]].*(1 .-S_0[[jx_sh_a ;jx_sh_b ;jx_sh_c]]));
decomp_mat[:,4,ccc]  = MP_k;
decomp_mat[:,5,ccc]  = @. irf_factor*log.((NetWorth_0 + Labor_Income_0 + NetWorth_Component)/(NetWorth_0 + Labor_Income_0));
decomp_mat[:,6,ccc]  = @. irf_factor*(Direct_Effects_k)./[k_a_0; k_b_0; k_c_0];
decomp_mat[:,7,ccc]  = decomp_mat[:,2,ccc] - decomp_mat[:,3,ccc].*decomp_mat[:,4,ccc].*decomp_mat[:,5,ccc] - decomp_mat[:,6,ccc];

  
decomp_mat2[1,ccc]  = irf_factor_1*NetWorth_0[1]/sum(NetWorth_0)*(NetWorth_1[1]/NetWorth_0[1]-sum(NetWorth_1)/sum(NetWorth_0)); 
decomp_mat2[2,ccc]  = irf_factor_1*NetWorth_0[1]/sum(NetWorth_0)*(Labor_Income_Component[1]/NetWorth_0[1] - sum(Labor_Income_Component)/sum(NetWorth_0));
decomp_mat2[3,ccc] = irf_factor_1*NetWorth_0[1]/sum(NetWorth_0)*(Bond_Component[1]/NetWorth_0[1] - sum(Bond_Component)/sum(NetWorth_0));
decomp_mat2[4,ccc]  = irf_factor_1*NetWorth_0[1]/sum(NetWorth_0)*(Profit_Component[1]/NetWorth_0[1] - sum(Profit_Component)/sum(NetWorth_0));
decomp_mat2[5,ccc]  = irf_factor_1*NetWorth_0[1]/sum(NetWorth_0)*(Capital_Price_Component[1]/NetWorth_0[1] - sum(Capital_Price_Component)/sum(NetWorth_0));
decomp_mat2[6,ccc]  = irf_factor_1*NetWorth_0[1]/sum(NetWorth_0)*(Transfer_Component[1]/NetWorth_0[1] - sum(Transfer_Component)/sum(NetWorth_0));
