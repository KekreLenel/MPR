## Produce Table 5 and Table 6 of the Paper: Targeted and Untargeted Moments

jx_z_shock  = 1;
jx_m_shock  = 2;
jx_p_shock  = 3;

jx_k        = 4;
jx_s_a      = 5;
jx_s_c      = 6;
jx_m        = 7;
jx_w_s      = 8;
jx_p        = 9;

jx_y        = 10;
jx_l        = 11;
jx_c_a      = 12;
jx_c_b      = 13;
jx_c_c      = 14;
jx_inv      = 15;
jx_pi       = 16;
jx_Erf      = 17;
jx_Erk      = 18;
jx_q        = 19;
jx_sh_a     = 20;
jx_sh_b     = 21;
jx_sh_c     = 22;
jx_sav      = 23;
jx_aggrw    = 24;
jx_w        = 25;
jx_Eexc     = 26;
jx_infl     = 27;
jx_nom_i    = 28;
jx_ksav     = 29;
jx_kap      = 30;
jx_v1       = 31;
jx_v2       = 32;
jx_v3       = 33;
jx_ql1      = 34;
jx_ql2      = 35;
jx_ql3      = 36;
jx_MPC1     = 37;
jx_MPC2     = 38;
jx_MPC3     = 39;
jx_MPS1     = 40;
jx_MPS2     = 41;
jx_MPS3     = 42;
jx_MPK1     = 43;
jx_MPK2     = 44;
jx_MPK3     = 45;
jx_dc_dErk1 = 46;
jx_dc_dErk2 = 47;
jx_dc_dErk3 = 48;
jx_db_dErk1 = 49;
jx_db_dErk2 = 50;
jx_db_dErk3 = 51;
jx_dk_dErk1 = 52;
jx_dk_dErk2 = 53;
jx_dk_dErk3 = 54;
jx_dc_dErf1 = 55;
jx_dc_dErf2 = 56;
jx_dc_dErf3 = 57;
jx_db_dErf1 = 58;
jx_db_dErf2 = 59;
jx_db_dErf3 = 60;
jx_dk_dErf1 = 61;
jx_dk_dErf2 = 62;
jx_dk_dErf3 = 63;
jx_k_a      = 64;
jx_k_b      = 65;
jx_k_c      = 66;
jx_sav_a    = 67;
jx_sav_b    = 68;
jx_sav_c    = 69;
jx_n_a      = 70;
jx_n_b      = 71;
jx_n_c      = 72;
jx_qcon     = 73;
jx_cona     = 74;
jx_conb     = 75;
jx_conc     = 76;
jx_mca      = 77;
jx_mcb      = 78;
jx_mcc      = 79;
jx_q1       = 80;
jx_q2       = 81;
jx_q3       = 82;
jx_q4       = 83;



# Shocks
z_shocks = temp_series[:,jx_z_shock];
m_shocks = temp_series[:,jx_m_shock];
p_shocks = temp_series[:,jx_p_shock];


# adjustment for government debt in the data
gov_adj = prm_vec[param_id].gov_debt/(1.0 - prm_vec[param_id].gov_debt);

# Cumiulative Growth
z_series = exp.(cumsum(z_shocks));

# Only Disaster Series
dis_series        = zeros(size(z_shocks));
dis_series[(z_shocks .>= -prm_vec[param_id].varphi_p .- sqrt(eps())) .& (z_shocks .<= -prm_vec[param_id].varphi_p .+ sqrt(eps()))] = z_shocks[(z_shocks .>= -prm_vec[param_id].varphi_p .- sqrt(eps())) .& (z_shocks .<= -prm_vec[param_id].varphi_p .+ sqrt(eps()))];
z_shocks_no_dis   = z_shocks-dis_series;

# Capital Stock
k_series        = temp_series[:,jx_k] .* z_series;
k_series_unadj  = temp_series[:,jx_k];
k_chosen_series = k_series  ./ exp.(dis_series); # Capital Chosen at t
k_chosen_series = @. [k_chosen_series[2:end]; k_chosen_series[end]];

# Wealth Shares
s_a_series = temp_series[:,jx_s_a];
s_c_series = temp_series[:,jx_s_c];
s_a_series_adj = @. (s_a_series + prm_vec[param_id].s_trgt_a*gov_adj)/(1+gov_adj);
s_c_series_adj = @. (s_c_series + prm_vec[param_id].s_trgt_c*gov_adj)/(1+gov_adj);

# Monetary Policy
m_series = temp_series[:,jx_m];

# Real Wage Chosen Last Period
w_state_series = temp_series[:,jx_w_s];

# Other states
p_series     = exp.(temp_series[:,jx_p]);

# Labor, Output
l_series = temp_series[:,jx_l];
y_series = temp_series[:,jx_y] .* z_series;

# Expected Returns
E_rf_series      = temp_series[:,jx_Erf];
E_rk_series      = temp_series[:,jx_Erk];
E_exc_series     = temp_series[:,jx_Eexc];

# Price of Capital
q_series = temp_series[:,jx_q];
qcon_series = temp_series[:,jx_qcon];

# Inflation, Nominal Rate
infl_series     = temp_series[:,jx_infl];
nom_i_series    = temp_series[:,jx_nom_i];
pi_series       = temp_series[:, jx_pi];

# Wages
w_series = temp_series[:,jx_w] .* z_series;

# Investment
inv_series = temp_series[:,jx_inv] .* z_series;

# Consumnption
c_a_series = @. prm_vec[param_id].lmbd_a * temp_series[:,jx_c_a] .* z_series; 
c_b_series = @. prm_vec[param_id].lmbd_b * temp_series[:,jx_c_b] .* z_series; 
c_c_series = @. (1-prm_vec[param_id].lmbd_a-prm_vec[param_id].lmbd_b) * temp_series[:,jx_c_c].* z_series; 
c_series   = @. c_a_series + c_b_series + c_c_series;

# Portfolio Shares
sh_a_series = temp_series[:,jx_sh_a];
sh_b_series = temp_series[:,jx_sh_b];
sh_c_series = temp_series[:,jx_sh_c];

# adjust portfolio share 
sh_a_adj_series = @. (temp_series[:,jx_sh_a]+prm_vec[param_id].s_trgt_a ./s_a_series*gov_adj) ./(1+prm_vec[param_id].s_trgt_a ./s_a_series*gov_adj);
sh_b_adj_series = @. (temp_series[:,jx_sh_b]+(1 - prm_vec[param_id].s_trgt_c - prm_vec[param_id].s_trgt_a) ./(1 -s_a_series - s_c_series)*gov_adj) ./(1+(1 - prm_vec[param_id].s_trgt_c - prm_vec[param_id].s_trgt_a) ./(1 -s_a_series - s_c_series)*gov_adj);
sh_c_adj_series = @. (temp_series[:,jx_sh_c]+prm_vec[param_id].s_trgt_c ./s_c_series*gov_adj) ./(1+prm_vec[param_id].s_trgt_c ./s_c_series*gov_adj);

# Others
aggr_wealth_series = temp_series[:,jx_aggrw] .* z_series;

# Savings; Depolyed Capital
sav_series      = temp_series[:,jx_sav] .* z_series;
ksav_series     = temp_series[:,jx_ksav] .* z_series;
h_kap_series    = temp_series[:,jx_kap] .* z_series;

# Derived Return Series
rf_series           = @. log.(nom_i_series[1:end.-1]  ./ infl_series[2:end]);
rf_series           = @. [rf_series[1]; rf_series];
rk_series           = @. log.(exp.(dis_series[2:end]).*(pi_series[2:end] + (1.0-prm_vec[param_id].delta)*q_series[2:end]) ./q_series[1:end-1]);
rk_series           = @. [rk_series[1]; rk_series];
exc_ret_series      = @. rk_series - rf_series;
rA_series           = @. log((1 + prm_vec[param_id].debt_to_equity )*exp.(rk_series) - prm_vec[param_id].debt_to_equity*exp.(rf_series));
exc_retA_series     = @. (1 + prm_vec[param_id].debt_to_equity )*(exp.(rk_series) - exp.(rf_series));
E_rA_series         = @. (1 + prm_vec[param_id].debt_to_equity )*E_rk_series - prm_vec[param_id].debt_to_equity*E_rf_series;
E_excA_series       = @. (1 + prm_vec[param_id].debt_to_equity )*(E_rk_series - E_rf_series);


# unlevered claim on capital, dividend flow as given by investment series
divk_series       = k_series[1:end-1].*(pi_series[1:end-1] + (1 - prm_vec[param_id].delta)*q_series[1:end-1]) 
                        - q_series[1:end-1].*k_series[2:end];
divk_series       = [divk_series; divk_series[end]]; 
divk_price_series = divk_series./(q_series.*k_chosen_series); 

# dividend series : aggregate profits, financed with fixed leverage, aggregate investment
div_series    = k_series[2:end-1].*(pi_series[2:end-1] + (1 - prm_vec[param_id].delta)*q_series[2:end-1] - prm_vec[param_id].debt_to_equity/(1.0+prm_vec[param_id].debt_to_equity).*q_series[1:end-2].*exp.(rf_series[2:end-1]) ) - 1/(1 + prm_vec[param_id].debt_to_equity).*q_series[2:end-1].*k_series[3:end];
div_series    = [div_series[1]; div_series; div_series[end]];
qE_series     = 1/(1 + prm_vec[param_id].debt_to_equity).*q_series.*k_chosen_series;

div_price_series = div_series./qE_series; 
div_price_smoothed_series = movsum(div_series,3)./qE_series;
div_price_smoothed_series[1:3] .= div_price_smoothed_series[1]*4;


MPC1 = temp_series[:,jx_MPC1];
MPC2 = temp_series[:,jx_MPC2];
MPC3 = temp_series[:,jx_MPC3];
MPK1 = temp_series[:,jx_MPK1];
MPK2 = temp_series[:,jx_MPK2];
MPK3 = temp_series[:,jx_MPK3];
MPS1 = temp_series[:,jx_MPS1];
MPS2 = temp_series[:,jx_MPS2];
MPS3 = temp_series[:,jx_MPS3];

MPR1 = MPK1  ./ MPS1;
MPR2 = MPK2  ./ MPS2;
MPR3 = MPK3  ./ MPS3;

q1_series       = temp_series[:,jx_q1];
q2_series       = temp_series[:,jx_q2];
q3_series       = temp_series[:,jx_q3];
q4_series       = temp_series[:,jx_q4];
yield1_series = 1 ./q1_series.-1;
yield2_series = (1 ./q2_series).^(1/2).-1;
yield3_series = (1 ./q3_series).^(1/3).-1;
yield4_series = (1 ./q4_series).^(1/4).-1;
treas_1y_series = yield4_series*4;

logq_series = log.(q_series);
logw_series = log.(w_series);
logl_series = log.(l_series);
loginv_series = log.(inv_series);
logc_series = log.(c_series);
logy_series = log.(y_series);

sav_a_series = temp_series[:,jx_sav_a].*z_series; 
sav_b_series = temp_series[:,jx_sav_b].*z_series; 
sav_c_series = temp_series[:,jx_sav_c].*z_series;
sav_ret_a_series = @. sav_a_series[1:end-1].*(sh_a_series[1:end-1].*exp.(rf_series[2:end]) + (1 - sh_a_series[1:end-1]).*exp.(rk_series[2:end]));
sav_ret_b_series = @. sav_b_series[1:end-1].*(sh_b_series[1:end-1].*exp.(rf_series[2:end]) + (1 - sh_b_series[1:end-1]).*exp.(rk_series[2:end]));
sav_ret_c_series = @. sav_c_series[1:end-1].*(sh_c_series[1:end-1].*exp.(rf_series[2:end]) + (1 - sh_c_series[1:end-1]).*exp.(rk_series[2:end]));
sav_ret_a_series = [sav_ret_a_series[1]; sav_ret_a_series];
sav_ret_b_series = [sav_ret_b_series[1]; sav_ret_b_series];
sav_ret_c_series = [sav_ret_c_series[1]; sav_ret_c_series];
z_series_previous = [1; z_series[1:end-1]];
bg_prm = gov_adj*mean(q_series.*k_chosen_series./z_series);
transfer_a_series = @. prm_vec[param_id].xi*prm_vec[param_id].s_bar_a*aggr_wealth_series./y_series -  prm_vec[param_id].xi*sav_ret_a_series./y_series - prm_vec[param_id].s_trgt_a*bg_prm*(z_series_previous.*exp.(rf_series) - z_series)./y_series;
transfer_b_series = @. prm_vec[param_id].xi*(1-prm_vec[param_id].s_bar_a - prm_vec[param_id].s_bar_c)*aggr_wealth_series./y_series -  prm_vec[param_id].xi*sav_ret_b_series./y_series - (1- prm_vec[param_id].s_trgt_a - prm_vec[param_id].s_trgt_c)*bg_prm*(z_series_previous.*exp.(rf_series) - z_series)./y_series;
transfer_c_series = @. prm_vec[param_id].xi*prm_vec[param_id].s_bar_c*aggr_wealth_series./y_series -  prm_vec[param_id].xi*sav_ret_c_series./y_series - prm_vec[param_id].s_trgt_c*bg_prm*(z_series_previous.*exp.(rf_series) - z_series)./y_series;
bg_series =  bg_prm*z_series;




