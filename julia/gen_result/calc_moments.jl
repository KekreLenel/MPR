## HP filtered y series
HP_LAMBDA = 1600;
num_periods = 1000;
x = log.(y_series[1:num_periods]);
include("hp_filter.jl");
y_filtered = xhp;

corr_yhp_Erk   = cor(y_filtered[2:num_periods], E_rk_series[2:num_periods]);
corr_yhp_Erf   = cor(y_filtered[2:num_periods], E_rf_series[2:num_periods]);
corr_yhp_EecxA = cor(y_filtered[2:num_periods], E_excA_series[2:num_periods]);
corr_dyhp_dErk   = cor(diff(y_filtered[2:num_periods]), diff(E_rk_series[2:num_periods]));
corr_dyhp_dErf   = cor(diff(y_filtered[2:num_periods]), diff(E_rf_series[2:num_periods]));
corr_dyhp_dEecxA = cor(diff(y_filtered[2:num_periods]), diff(E_excA_series[2:num_periods]));
corr_dyhp_Erk   = cor(diff(y_filtered[2:num_periods]), (E_rk_series[2:num_periods-1]));
corr_dyhp_Erf   = cor(diff(y_filtered[2:num_periods]), (E_rf_series[2:num_periods-1]));
corr_dyhp_EecxA = cor(diff(y_filtered[2:num_periods]), (E_excA_series[2:num_periods-1]));


## HP filtered y serie
x = log.(y_series[1:num_periods]);
y_detrend = ols_detrend(x);

corr_ydt_Erk   = cor(y_detrend[2:num_periods], E_rk_series[2:num_periods]);
corr_ydt_Erf   = cor(y_detrend[2:num_periods], E_rf_series[2:num_periods]);
corr_ydt_EecxA = cor(y_detrend[2:num_periods], E_excA_series[2:num_periods]);
corr_dydt_dErk   = cor(diff(y_detrend[2:num_periods]), diff(E_rk_series[2:num_periods]));
corr_dydt_dErf   = cor(diff(y_detrend[2:num_periods]), diff(E_rf_series[2:num_periods]));
corr_dydt_dEecxA = cor(diff(y_detrend[2:num_periods]), diff(E_excA_series[2:num_periods]));


std_k       = std(log.(k_series[2:end]./k_series[1:end-1]));
std_w       = std(log.(w_series[2:end]./w_series[1:end-1]));
std_l       = std(log.(l_series[2:end]./l_series[1:end-1]));
std_y       = std(log.(y_series[2:end]./y_series[1:end-1]));
std_c       = std(log.(c_series[2:end]./c_series[1:end-1]));
std_ca      = std(log.(c_a_series[2:end]./c_a_series[1:end-1]));
std_cb      = std(log.(c_b_series[2:end]./c_b_series[1:end-1]));
std_cc      = std(log.(c_c_series[2:end]./c_c_series[1:end-1]));
std_inv     = std(log.(inv_series[2:end]./inv_series[1:end-1]));


corr_cw     = cor(log.(c_series[2:end]./c_series[1:end-1]),log.(w_series[2:end]./w_series[1:end-1]));
corr_cl     = cor(log.(c_series[2:end]./c_series[1:end-1]),log.(l_series[2:end]./l_series[1:end-1]));
corr_cy     = cor(log.(c_series[2:end]./c_series[1:end-1]),log.(y_series[2:end]./y_series[1:end-1]));
corr_cinv   = cor(log.(c_series[2:end]./c_series[1:end-1]),log.(inv_series[2:end]./inv_series[1:end-1]));

corr_yw     = cor(log.(y_series[2:end]./y_series[1:end-1]),log.(w_series[2:end]./w_series[1:end-1]));
corr_yl     = cor(log.(y_series[2:end]./y_series[1:end-1]),log.(l_series[2:end]./l_series[1:end-1]));
corr_yinv  = cor(log.(y_series[2:end]./y_series[1:end-1]),log.(inv_series[2:end]./inv_series[1:end-1]));

skew_k      = skewness(log.(k_series[2:end]./k_series[1:end-1]));
skew_w      = skewness(log.(w_series[2:end]./w_series[1:end-1]));
skew_l      = skewness(log.(l_series[2:end]./l_series[1:end-1]));
skew_y      = skewness(log.(y_series[2:end]./y_series[1:end-1]));
skew_c      = skewness(log.(c_series[2:end]./c_series[1:end-1]));
skew_inv    = skewness(log.(inv_series[2:end]./inv_series[1:end-1]));

std_k_lev   = std(log.(k_series[1:end]));
std_w_lev   = std(log.(w_series[1:end]));
std_l_lev   = std(log.(l_series[1:end]));
std_y_lev   = std(log.(y_series[1:end]));
std_c_lev   = std(log.(c_series[1:end]));
std_inv_lev = std(log.(inv_series[1:end]));

std_sa      = std(s_a_series[1:end]);
std_sc      = std(s_c_series[1:end]);
std_divprice = std(div_price_series);
std_divprice_smoothed = std(div_price_smoothed_series);

std_infl    = std((infl_series[1:end].-1));
std_Erf      = std((E_rf_series[1:end].-1));
std_rf      = std((rf_series[1:end].-1));
std_Erk      = std((E_rk_series[1:end].-1));
std_Eexc     = std((E_exc_series[1:end]));
std_EexcA    = std((E_excA_series[1:end]));
std_realized_excA = std(exp.(rA_series[1:end]).-exp.(rf_series[1:end]));
std_nom_i   = std(nom_i_series[1:end].-1);

ac_excA      = cor(E_excA_series[3:end-1], E_excA_series[2:end-2]);
ac_Erf      = cor(E_rf_series[3:end-1], E_rf_series[2:end-2]);
ac_rf      = cor(exp.(rf_series[3:end-1]).-1, exp.(rf_series[2:end-2]).-1);
ac_divprice  = cor( div_price_series[3:end-1],  div_price_series[2:end-2]);

yoy_ac_divprice  = cor( div_price_series[6:end-1],  div_price_series[2:end-1-4]);
ac_divprice_smoothed  = cor(div_price_smoothed_series[5:end], div_price_smoothed_series[4:end-1]);
ac_yoy_divprice_smoothed  = cor(div_price_smoothed_series[8:end], div_price_smoothed_series[4:end-4]);



global io *= @sprintf("\nNUMERICAL VERIFICATION \n");
global io *= @sprintf("-----------------------\n\n");
global io *= @sprintf("STATES MEAN AND STD \n");
global io *= @sprintf("-----------------------\n");
global io *= @sprintf("        avrg.     grid_mean  std.  grid_dev\n");
global io *= @sprintf("k       %6.2f    %6.2f   %6.2f    %6.2f\n",       mean(k_series_unadj), k_grid_mean,          std(k_series_unadj),     k_grid_dev);
global io *= @sprintf("s_a     %6.2f%%   %6.2f%%  %6.2f%%   %6.2f%%\n", 100*mean(s_a_series),   100*s_a_grid_mean,  100*std(s_a_series), 100*s_a_grid_dev);
global io *= @sprintf("s_c     %6.2f%%   %6.2f%%  %6.2f%%   %6.2f%%\n", 100*mean(s_c_series),   100*s_c_grid_mean,  100*std(s_c_series), 100*s_c_grid_dev);
global io *= @sprintf("p       %6.2f%%   %6.2f%%  %6.2f%%   %6.2f%%\n", 100*mean(p_series),       100*dis_grid_mean,      100*std(p_series),     100*dis_grid_dev);
global io *= @sprintf("w       %6.2f    %6.2f   %6.2f    %6.2f\n", mean(w_state_series), w_grid_mean,      std(w_state_series), w_grid_dev);
global io *= @sprintf("m       %6.2f    %6.2f   %6.2f    %6.2f\n",       mean(m_series),       m_grid_mean,          std(m_series),           m_grid_dev);



global io *= @sprintf("FIRST MOMENTS \n");
global io *= @sprintf("-----------------------\n\n");
global io *= @sprintf("REAL OUTCOMES \n");
global io *= @sprintf("-----------------------\n");
global io *= @sprintf("c        %10.6f\n", mean(c_series[2:end-1]./z_series[2:end-1]));
global io *= @sprintf("k        %10.6f\n", mean(k_series[2:end-1]./z_series[2:end-1]));
global io *= @sprintf("w        %10.6f\n", mean(w_series[2:end-1]./z_series[2:end-1]));
global io *= @sprintf("y        %10.6f\n", mean(y_series[2:end-1]./z_series[2:end-1]));
global io *= @sprintf("l        %10.6f\n", mean(l_series[2:end-1]));
global io *= @sprintf("inv      %10.6f\n", mean(inv_series[2:end-1]./z_series[2:end-1]));
global io *= @sprintf("\nPRICES (ANNUALIZED) \n");
global io *= @sprintf("-----------------------\n");
global io *= @sprintf("infl     %10.6f%%\n", 400*mean(log.(infl_series[2:end-1])));
global io *= @sprintf("E[rf]    %10.6f%%\n", 400*mean(exp.(rf_series[2:end-1]).-1));
global io *= @sprintf("E[rk]    %10.6f%%\n", 400*mean(exp.(rk_series[2:end-1]).-1));
global io *= @sprintf("E[rk-rf] %10.6f%%\n", 400*mean(exp.(rk_series[2:end-1]) .- exp.(rf_series[2:end-1])));
global io *= @sprintf("E[rA-rf] %10.6f%%\n", 400*mean(exp.(rA_series[2:end-1]) .- exp.(rf_series[2:end-1])));
global io *= @sprintf("\nAGENTS WEALTH AND PORTFOLIOS\n");
global io *= @sprintf("-----------------------\n");
global io *= @sprintf("s_a      %10.6f\n", mean(s_a_series_adj[2:end-1]));
global io *= @sprintf("s_c      %10.6f\n", mean(s_c_series_adj[2:end-1]));
global io *= @sprintf("share a  %10.6f\n", mean(sh_a_adj_series[2:end-1]));
global io *= @sprintf("share b  %10.6f\n", mean(sh_b_adj_series[2:end-1]));
global io *= @sprintf("share c  %10.6f\n", mean(sh_c_adj_series[2:end-1]));
global io *= @sprintf("\nMPCs and MPK\n");
global io *= @sprintf("-----------------------\n");
global io *= @sprintf("MPC A        %10.6f\n", mean(MPC1[2:end-1]));
global io *= @sprintf("MPC B        %10.6f\n", mean(MPC2[2:end-1]));
global io *= @sprintf("MPC C        %10.6f\n", mean(MPC3[2:end-1]));
global io *= @sprintf("MPS A        %10.6f\n", mean(MPS1[2:end-1]));
global io *= @sprintf("MPS B        %10.6f\n", mean(MPS2[2:end-1]));
global io *= @sprintf("MPS C        %10.6f\n", mean(MPS3[2:end-1]));
global io *= @sprintf("MPK A        %10.6f\n", mean(MPK1[2:end-1]));
global io *= @sprintf("MPK B        %10.6f\n", mean(MPK2[2:end-1]));
global io *= @sprintf("MPK C        %10.6f\n", mean(MPK3[2:end-1]));
global io *= @sprintf("\n-----------------------\n");
global io *= @sprintf("\nDIV/PRICE\n");
global io *= @sprintf("-----------------------\n");
global io *= @sprintf("D/P        %10.6f\n", mean(div_price_series[2:end-1]));
global io *= @sprintf("-----------------------\n\n");
global io *= @sprintf("bg/y         %10.6f%%\n", 100*mean(bg_series[2:end-1]./y_series[2:end-1]));
global io *= @sprintf("bg/(qk + bg) %10.6f%%\n", 100*mean(bg_series[2:end-1]./(aggr_wealth_series[2:end-1] .+ bg_series[2:end-1])));
global io *= @sprintf("bg/z         %10.6f%%\n", 100*mean(bg_series[2:end-1]./(z_series[2:end-1])));
global io *= @sprintf("-----------------------\n\n");



global io *= @sprintf("BUSINESS CYCLE MOMENTS \n");
global io *= @sprintf("-----------------------\n\n");
global io *= @sprintf("STD  log. growth \n");
global io *= @sprintf("-----------------------\n\n");
global io *= @sprintf("std(k)        %10.6f%%\n", 100*mean(std_k));
global io *= @sprintf("std(w)        %10.6f%%\n", 100*mean(std_w));
global io *= @sprintf("std(y)        %10.6f%%\n", 100*mean(std_y));
global io *= @sprintf("std(c)        %10.6f%%\n", 100*mean(std_c));
global io *= @sprintf("std(c_a)      %10.6f%%\n", 100*mean(std_ca));
global io *= @sprintf("std(c_b)      %10.6f%%\n", 100*mean(std_cb));
global io *= @sprintf("std(c_c)      %10.6f%%\n", 100*mean(std_cc));
global io *= @sprintf("std(l)        %10.6f%%\n", 100*mean(std_l));
global io *= @sprintf("std(inv)      %10.6f%%\n", 100*mean(std_inv));
global io *= @sprintf("\n-----------------------\n\n");

global io *= @sprintf("CORR log growth \n");
global io *= @sprintf("-----------------------\n\n");
global io *= @sprintf("cor(c,w)        %10.6f%%\n", 100*mean(corr_cw));
global io *= @sprintf("cor(c,y)        %10.6f%%\n", 100*mean(corr_cy));
global io *= @sprintf("cor(c,l)        %10.6f%%\n", 100*mean(corr_cl));
global io *= @sprintf("cor(c,inv)      %10.6f%%\n", 100*mean(corr_cinv));
global io *= @sprintf("-----------------------\n\n");
global io *= @sprintf("cor(y,w)        %10.6f%%\n", 100*mean(corr_yw));
global io *= @sprintf("cor(y,l)        %10.6f%%\n", 100*mean(corr_yl));
global io *= @sprintf("cor(y,inv)      %10.6f%%\n", 100*mean(corr_yinv));
global io *= @sprintf("\n-----------------------\n\n");


global io *= @sprintf("CORR hp filtered log y and returns \n");
global io *= @sprintf("-----------------------\n\n");
global io *= @sprintf("cor(y,E(rf))     %10.6f%%\n", 100*mean(corr_yhp_Erf));
global io *= @sprintf("cor(y,E(rk)      %10.6f%%\n", 100*mean(corr_yhp_Erk));
global io *= @sprintf("cor(y,E(exc)     %10.6f%%\n", 100*mean(corr_yhp_EecxA));
global io *= @sprintf("-----------------------\n");
global io *= @sprintf("cor(dy,E(rf))   %10.6f%%\n", 100*mean(corr_dyhp_Erf));
global io *= @sprintf("cor(dy,E(rk)    %10.6f%%\n", 100*mean(corr_dyhp_Erk));
global io *= @sprintf("cor(dy,E(exc)   %10.6f%%\n", 100*mean(corr_dyhp_EecxA));
global io *= @sprintf("-----------------------\n");
global io *= @sprintf("cor(dy,dE(rf))   %10.6f%%\n", 100*mean(corr_dyhp_dErf));
global io *= @sprintf("cor(dy,dE(rk)    %10.6f%%\n", 100*mean(corr_dyhp_dErk));
global io *= @sprintf("cor(dy,dE(exc)   %10.6f%%\n", 100*mean(corr_dyhp_dEecxA));
global io *= @sprintf("-----------------------\n\n");

global io *= @sprintf("CORR detrended log. y and returns \n");
global io *= @sprintf("-----------------------\n\n");
global io *= @sprintf("cor(y,E(rf))     %10.6f%%\n", 100*mean(corr_ydt_Erf));
global io *= @sprintf("cor(y,E(rk)      %10.6f%%\n", 100*mean(corr_ydt_Erk));
global io *= @sprintf("cor(y,E(exc)     %10.6f%%\n", 100*mean(corr_ydt_EecxA));
global io *= @sprintf("-----------------------\n");
global io *= @sprintf("cor(dy,dE(rf))   %10.6f%%\n", 100*mean(corr_dydt_dErf));
global io *= @sprintf("cor(dy,dE(rk)    %10.6f%%\n", 100*mean(corr_dydt_dErk));
global io *= @sprintf("cor(dy,dE(exc)   %10.6f%%\n", 100*mean(corr_dydt_dEecxA));
global io *= @sprintf("-----------------------\n\n");

global io *= @sprintf("SKEW log. growth\n");
global io *= @sprintf("-----------------------\n\n");
global io *= @sprintf("skew(k)        %10.6f%%\n", 100*mean(skew_k));
global io *= @sprintf("skew(w)        %10.6f%%\n", 100*mean(skew_w));
global io *= @sprintf("skew(y)        %10.6f%%\n", 100*mean(skew_y));
global io *= @sprintf("skew(c)        %10.6f%%\n", 100*mean(skew_c));
global io *= @sprintf("skew(l)        %10.6f%%\n", 100*mean(skew_l));
global io *= @sprintf("skew(inv)      %10.6f%%\n", 100*mean(skew_inv));
global io *= @sprintf("\n-----------------------\n\n");

global io *= @sprintf("PRICE MOMENTS \n");
global io *= @sprintf("-----------------------\n\n");
global io *= @sprintf("STD  returns (ANNUALIZED) \n");
global io *= @sprintf("-----------------------\n\n");
global io *= @sprintf("std(infl)     %10.6f%%\n", 400*mean(std_infl));
global io *= @sprintf("std(rf)       %10.6f%%\n", 400*mean(std_rf));
global io *= @sprintf("std(E[rf])    %10.6f%%\n", 400*mean(std_Erf));
global io *= @sprintf("std(E[rk])    %10.6f%%\n", 400*mean(std_Erk));
global io *= @sprintf("std(E[rk-rf]) %10.6f%%\n", 400*mean(std_Eexc));
global io *= @sprintf("std(E[rA-rf]) %10.6f%%\n", 400*mean(std_EexcA));
global io *= @sprintf("std(rA-rf])   %10.6f%%\n\n", 400*mean(std_realized_excA));


global io *= @sprintf("\nDIVIDEND PRICE STD\n");
global io *= @sprintf("std(d/p) %10.6f%%\n", 100*mean(std_divprice));
global io *= @sprintf("\nSMOOTHED DIVIDEND PRICE STD\n");
global io *= @sprintf("std(d/p) %10.6f%%\n", 100*mean(std_divprice_smoothed));

global io *= @sprintf("\nAUTOcor \n");
global io *= @sprintf("-----------------------\n\n");
global io *= @sprintf("ac(E[rA-rf])   %10.6f%%\n", 100*mean(ac_excA));
global io *= @sprintf("ac(E[rf])      %10.6f%%\n", 100*mean(ac_Erf));
global io *= @sprintf("ac(rf)         %10.6f%%\n", 100*mean(ac_rf));
global io *= @sprintf("ac(d/p)        %10.6f%%\n", 100*mean(ac_divprice));
global io *= @sprintf("\nAUTOcor YEAR OVER YEAR  \n");
global io *= @sprintf("-----------------------\n\n");
global io *= @sprintf("ac(d/p)        %10.6f%%\n", 100*mean(yoy_ac_divprice));
global io *= @sprintf("\nAUTOcor SMOOTH \n");
global io *= @sprintf("ac(d/p) smooth %10.6f%%\n", 100*mean(ac_divprice_smoothed));
global io *= @sprintf("\nAUTOcor SMOOTH YEAR OVER YEAR \n");
global io *= @sprintf("ac(d/p) smooth %10.6f%%\n", 100*ac_yoy_divprice_smoothed);

global io *= @sprintf("\nWEALTH SHARE STD\n");
global io *= @sprintf("std(sa)       %10.6f%%\n", 100*std_sa);
global io *= @sprintf("std(sc)       %10.6f%%\n", 100*std_sc);
global io *= @sprintf("\n-----------------------\n\n");

