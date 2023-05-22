% -------------------------------------------------------------------------
% calc_moments.m: calculate simulation moments 
% -------------------------------------------------------------------------
% authors:         Rohan Kekre and Moritz Lenel
% for updates see: https://github.com/KekreLenel/MPR
% -------------------------------------------------------------------------
% HP filtered y series
HP_LAMBDA = 1600;
x = log(y_series(1:10000));
hp_filter;
y_filtered = xhp;

corr_yhp_Erk     = corr(y_filtered(2:10000), E_rk_series(2:10000));
corr_yhp_Erf     = corr(y_filtered(2:10000), E_rf_series(2:10000));
corr_yhp_EecxA   = corr(y_filtered(2:10000), E_excA_series(2:10000));
corr_dyhp_dErk   = corr(diff(y_filtered(2:10000)), diff(E_rk_series(2:10000)));
corr_dyhp_dErf   = corr(diff(y_filtered(2:10000)), diff(E_rf_series(2:10000)));
corr_dyhp_dEecxA = corr(diff(y_filtered(2:10000)), diff(E_excA_series(2:10000)));
corr_dyhp_Erk    = corr(diff(y_filtered(2:10000)), (E_rk_series(2:10000-1)));
corr_dyhp_Erf    = corr(diff(y_filtered(2:10000)), (E_rf_series(2:10000-1)));
corr_dyhp_EecxA  = corr(diff(y_filtered(2:10000)), (E_excA_series(2:10000-1)));


% HP filtered y serie
x = log(y_series(1:10000));
y_detrend = detrend(x);

corr_ydt_Erk     = corr(y_detrend(2:10000), E_rk_series(2:10000));
corr_ydt_Erf     = corr(y_detrend(2:10000), E_rf_series(2:10000));
corr_ydt_EecxA   = corr(y_detrend(2:10000), E_excA_series(2:10000));
corr_dydt_dErk   = corr(diff(y_detrend(2:10000)), diff(E_rk_series(2:10000)));
corr_dydt_dErf   = corr(diff(y_detrend(2:10000)), diff(E_rf_series(2:10000)));
corr_dydt_dEecxA = corr(diff(y_detrend(2:10000)), diff(E_excA_series(2:10000)));


std_k       = std(log(k_series(2:end)./k_series(1:end-1)));
std_w       = std(log(w_series(2:end)./w_series(1:end-1)));
std_l       = std(log(l_series(2:end)./l_series(1:end-1)));
std_y       = std(log(y_series(2:end)./y_series(1:end-1)));
std_c       = std(log(c_series(2:end)./c_series(1:end-1)));
std_ca      = std(log(c_a_series(2:end)./c_a_series(1:end-1)));
std_cb      = std(log(c_b_series(2:end)./c_b_series(1:end-1)));
std_cc      = std(log(c_c_series(2:end)./c_c_series(1:end-1)));
std_inv     = std(log(inv_series(2:end)./inv_series(1:end-1)));


corr_cw     = corr(log(c_series(2:end)./c_series(1:end-1)),log(w_series(2:end)./w_series(1:end-1)));
corr_cl     = corr(log(c_series(2:end)./c_series(1:end-1)),log(l_series(2:end)./l_series(1:end-1)));
corr_cy     = corr(log(c_series(2:end)./c_series(1:end-1)),log(y_series(2:end)./y_series(1:end-1)));
corr_cinv   = corr(log(c_series(2:end)./c_series(1:end-1)),log(inv_series(2:end)./inv_series(1:end-1)));

corr_yw     = corr(log(y_series(2:end)./y_series(1:end-1)),log(w_series(2:end)./w_series(1:end-1)));
corr_yl     = corr(log(y_series(2:end)./y_series(1:end-1)),log(l_series(2:end)./l_series(1:end-1)));
corr_yinv   = corr(log(y_series(2:end)./y_series(1:end-1)),log(inv_series(2:end)./inv_series(1:end-1)));

skew_k      = skewness(log(k_series(2:end)./k_series(1:end-1)));
skew_w      = skewness(log(w_series(2:end)./w_series(1:end-1)));
skew_l      = skewness(log(l_series(2:end)./l_series(1:end-1)));
skew_y      = skewness(log(y_series(2:end)./y_series(1:end-1)));
skew_c      = skewness(log(c_series(2:end)./c_series(1:end-1)));
skew_inv    = skewness(log(inv_series(2:end)./inv_series(1:end-1)));

std_k_lev   = std(log(k_series(1:end)));
std_w_lev   = std(log(w_series(1:end)));
std_l_lev   = std(log(l_series(1:end)));
std_y_lev   = std(log(y_series(1:end)));
std_c_lev   = std(log(c_series(1:end)));
std_inv_lev = std(log(inv_series(1:end)));

std_sa       = std(s_a_series(1:end));
std_sc       = std(s_c_series(1:end));
std_divprice = std(div_price_series);
std_divprice_smoothed = std(div_price_smoothed_series);

std_infl     = std((infl_series(1:end)-1));
std_Erf      = std((E_rf_series(1:end)-1));
std_rf       = std((exp(rf_series(1:end))-1));
std_Erk      = std((E_rk_series(1:end)-1));
std_Eexc     = std((E_exc_series(1:end)));
std_EexcA    = std((E_excA_series(1:end)));
std_realized_excA = std((exp(rA_series(1:end))-exp(rf_series(1:end))));
std_nom_i    = std(nom_i_series(1:end)-1);

ac_excA      = corr(E_excA_series(3:end-1), E_excA_series(2:end-2));
ac_Erf       = corr(E_rf_series(3:end-1), E_rf_series(2:end-2));
ac_rf        = corr(exp(rf_series(3:end-1))-1, exp(rf_series(2:end-2))-1);
ac_divprice  = corr(div_price_series(3:end-1), div_price_series(2:end-2));

yoy_ac_divprice = corr(div_price_series(6:end-1), div_price_series(2:end-1-4));
ac_divprice_smoothed = corr(div_price_smoothed_series(5:end), div_price_smoothed_series(4:end-1));
ac_yoy_divprice_smoothed = corr(div_price_smoothed_series(8:end), div_price_smoothed_series(4:end-4));

fprintf(fileID,'\nNUMERICAL VERIFICATION \n');
fprintf(fileID,'-----------------------\n\n');
fprintf(fileID,'STATES MEAN AND STD \n');
fprintf(fileID,'-----------------------\n');
fprintf(fileID,'        avrg.     grid_mean  std.  grid_dev\n');
fprintf(fileID,'k       %6.2f    %6.2f   %6.2f    %6.2f\n',       mean(k_series_unadj), k_grid_mean,          std(k_series_unadj),     k_grid_dev);
fprintf(fileID,'s_a     %6.2f%%   %6.2f%%  %6.2f%%   %6.2f%%\n', 100*mean(s_a_series),   100*s_a_grid_mean,  100*std(s_a_series), 100*s_a_grid_dev);
fprintf(fileID,'s_c     %6.2f%%   %6.2f%%  %6.2f%%   %6.2f%%\n', 100*mean(s_c_series),   100*s_c_grid_mean,  100*std(s_c_series), 100*s_c_grid_dev);
fprintf(fileID,'p       %6.2f%%   %6.2f%%  %6.2f%%   %6.2f%%\n', 100*mean(p_series),       100*dis_grid_mean,      100*std(p_series),     100*dis_grid_dev);
fprintf(fileID,'w       %6.2f    %6.2f   %6.2f    %6.2f\n', mean(w_state_series), w_grid_mean,      std(w_state_series), w_grid_dev);
fprintf(fileID,'m       %6.2f    %6.2f   %6.2f    %6.2f\n',       mean(m_series),       m_grid_mean,          std(m_series),           m_grid_dev);


fprintf(fileID,'FIRST MOMENTS \n');
fprintf(fileID,'-----------------------\n\n');
fprintf(fileID,'REAL OUTCOMES \n');
fprintf(fileID,'-----------------------\n');
fprintf(fileID,'c        %10.6f\n', mean(c_series(2:end-1)./z_series(2:end-1)));
fprintf(fileID,'k        %10.6f\n', mean(k_series(2:end-1)./z_series(2:end-1)));
fprintf(fileID,'w        %10.6f\n', mean(w_series(2:end-1)./z_series(2:end-1)));
fprintf(fileID,'y        %10.6f\n', mean(y_series(2:end-1)./z_series(2:end-1)));
fprintf(fileID,'l        %10.6f\n', mean(l_series(2:end-1)));
fprintf(fileID,'inv      %10.6f\n', mean(inv_series(2:end-1)./z_series(2:end-1)));
fprintf(fileID,'\nPRICES (ANNUALIZED) \n');
fprintf(fileID,'-----------------------\n');
fprintf(fileID,'infl     %10.6f%%\n', 400*mean(infl_series(2:end-1) - 1.0));
fprintf(fileID,'rf       %10.6f%%\n', 400*mean(exp(rf_series(2:end-1))-1));
fprintf(fileID,'rk       %10.6f%%\n', 400*mean(exp(rk_series(2:end-1))-1));
fprintf(fileID,'rk-rf    %10.6f%%\n', 400*mean(exp(rk_series(2:end-1)) - exp(rf_series(2:end-1))));
fprintf(fileID,'rA-rf    %10.6f%%\n', 400*mean(exp(rA_series(2:end-1)) - exp(rf_series(2:end-1))));
fprintf(fileID,'\nAGENTS WEALTH AND PORTFOLIOS\n');
fprintf(fileID,'-----------------------\n');
fprintf(fileID,'s_a      %10.6f\n', mean(s_a_series_adj(2:end-1)));
fprintf(fileID,'s_c      %10.6f\n', mean(s_c_series_adj(2:end-1)));
fprintf(fileID,'share a  %10.6f\n', mean(sh_a_adj_series(2:end-1)));
fprintf(fileID,'share b  %10.6f\n', mean(sh_b_adj_series(2:end-1)));
fprintf(fileID,'share c  %10.6f\n', mean(sh_c_adj_series(2:end-1)));
fprintf(fileID,'\nMPCs and MPK\n');
fprintf(fileID,'-----------------------\n');
fprintf(fileID,'MPC A        %10.6f\n', mean(MPC1(2:end-1)));
fprintf(fileID,'MPC B        %10.6f\n', mean(MPC2(2:end-1)));
fprintf(fileID,'MPC C        %10.6f\n', mean(MPC3(2:end-1)));
fprintf(fileID,'MPS A        %10.6f\n', mean(MPS1(2:end-1)));
fprintf(fileID,'MPS B        %10.6f\n', mean(MPS2(2:end-1)));
fprintf(fileID,'MPS C        %10.6f\n', mean(MPS3(2:end-1)));
fprintf(fileID,'MPK A        %10.6f\n', mean(MPK1(2:end-1)));
fprintf(fileID,'MPK B        %10.6f\n', mean(MPK2(2:end-1)));
fprintf(fileID,'MPK C        %10.6f\n', mean(MPK3(2:end-1)));
fprintf(fileID,'\n-----------------------\n');
fprintf(fileID,'\nDIV/PRICE\n');
fprintf(fileID,'-----------------------\n');
fprintf(fileID,'D/P        %10.6f\n', mean(div_price_series(2:end-1)));
fprintf(fileID,'-----------------------\n\n');
fprintf(fileID,'bg/y         %10.6f%%\n', 100*mean(bg_series(2:end-1)./y_series(2:end-1)));
fprintf(fileID,'bg/(qk + bg) %10.6f%%\n', 100*mean(bg_series(2:end-1)./(aggr_wealth_series(2:end-1) + bg_series(2:end-1))));
fprintf(fileID,'bg/z         %10.6f%%\n', 100*mean(bg_series(2:end-1)./(z_series(2:end-1))));
fprintf(fileID,'-----------------------\n\n');



fprintf(fileID,'BUSINESS CYCLE MOMENTS \n');
fprintf(fileID,'-----------------------\n\n');
fprintf(fileID,'STD  log growth \n');
fprintf(fileID,'-----------------------\n\n');
fprintf(fileID,'std(k)        %10.6f%%\n', 100*std_k);
fprintf(fileID,'std(w)        %10.6f%%\n', 100*std_w);
fprintf(fileID,'std(y)        %10.6f%%\n', 100*std_y);
fprintf(fileID,'std(c)        %10.6f%%\n', 100*std_c);
fprintf(fileID,'std(c_a)      %10.6f%%\n', 100*std_ca);
fprintf(fileID,'std(c_b)      %10.6f%%\n', 100*std_cb);
fprintf(fileID,'std(c_c)      %10.6f%%\n', 100*std_cc);
fprintf(fileID,'std(l)        %10.6f%%\n', 100*std_l);
fprintf(fileID,'std(inv)      %10.6f%%\n', 100*std_inv);
fprintf(fileID,'\n-----------------------\n\n');

fprintf(fileID,'CORR log growth \n');
fprintf(fileID,'-----------------------\n\n');
fprintf(fileID,'corr(c,w)        %10.6f%%\n', 100*corr_cw);
fprintf(fileID,'corr(c,y)        %10.6f%%\n', 100*corr_cy);
fprintf(fileID,'corr(c,l)        %10.6f%%\n', 100*corr_cl);
fprintf(fileID,'corr(c,inv)      %10.6f%%\n', 100*corr_cinv);
fprintf(fileID,'-----------------------\n\n');
fprintf(fileID,'corr(y,w)        %10.6f%%\n', 100*corr_yw);
fprintf(fileID,'corr(y,l)        %10.6f%%\n', 100*corr_yl);
fprintf(fileID,'corr(y,inv)      %10.6f%%\n', 100*corr_yinv);
fprintf(fileID,'\n-----------------------\n\n');


fprintf(fileID,'CORR hp filtered log y and returns \n');
fprintf(fileID,'-----------------------\n\n');
fprintf(fileID,'corr(y,E(rf))     %10.6f%%\n', 100*corr_yhp_Erf);
fprintf(fileID,'corr(y,E(rk)      %10.6f%%\n', 100*corr_yhp_Erk);
fprintf(fileID,'corr(y,E(exc)     %10.6f%%\n', 100*corr_yhp_EecxA);
fprintf(fileID,'-----------------------\n');
fprintf(fileID,'corr(dy,E(rf))   %10.6f%%\n', 100*corr_dyhp_Erf);
fprintf(fileID,'corr(dy,E(rk)    %10.6f%%\n', 100*corr_dyhp_Erk);
fprintf(fileID,'corr(dy,E(exc)   %10.6f%%\n', 100*corr_dyhp_EecxA);
fprintf(fileID,'-----------------------\n');
fprintf(fileID,'corr(dy,dE(rf))   %10.6f%%\n', 100*corr_dyhp_dErf);
fprintf(fileID,'corr(dy,dE(rk)    %10.6f%%\n', 100*corr_dyhp_dErk);
fprintf(fileID,'corr(dy,dE(exc)   %10.6f%%\n', 100*corr_dyhp_dEecxA);
fprintf(fileID,'-----------------------\n\n');

fprintf(fileID,'CORR detrended log y and returns \n');
fprintf(fileID,'-----------------------\n\n');
fprintf(fileID,'corr(y,E(rf))     %10.6f%%\n', 100*corr_ydt_Erf);
fprintf(fileID,'corr(y,E(rk)      %10.6f%%\n', 100*corr_ydt_Erk);
fprintf(fileID,'corr(y,E(exc)     %10.6f%%\n', 100*corr_ydt_EecxA);
fprintf(fileID,'-----------------------\n');
fprintf(fileID,'corr(dy,dE(rf))   %10.6f%%\n', 100*corr_dydt_dErf);
fprintf(fileID,'corr(dy,dE(rk)    %10.6f%%\n', 100*corr_dydt_dErk);
fprintf(fileID,'corr(dy,dE(exc)   %10.6f%%\n', 100*corr_dydt_dEecxA);
fprintf(fileID,'-----------------------\n\n');

fprintf(fileID,'SKEW log growth\n');
fprintf(fileID,'-----------------------\n\n');
fprintf(fileID,'skew(k)        %10.6f%%\n', 100*skew_k);
fprintf(fileID,'skew(w)        %10.6f%%\n', 100*skew_w);
fprintf(fileID,'skew(y)        %10.6f%%\n', 100*skew_y);
fprintf(fileID,'skew(c)        %10.6f%%\n', 100*skew_c);
fprintf(fileID,'skew(l)        %10.6f%%\n', 100*skew_l);
fprintf(fileID,'skew(inv)      %10.6f%%\n', 100*skew_inv);
fprintf(fileID,'\n-----------------------\n\n');

fprintf(fileID,'PRICE MOMENTS \n');
fprintf(fileID,'-----------------------\n\n');
fprintf(fileID,'STD  returns (ANNUALIZED) \n');
fprintf(fileID,'-----------------------\n\n');
fprintf(fileID,'std(infl)     %10.6f%%\n', 400*std_infl);
fprintf(fileID,'std(rf)       %10.6f%%\n', 400*std_rf);
fprintf(fileID,'std(E[rf])    %10.6f%%\n', 400*std_Erf);
fprintf(fileID,'std(E[rk])    %10.6f%%\n', 400*std_Erk);
fprintf(fileID,'std(E[rk-rf]) %10.6f%%\n', 400*std_Eexc);
fprintf(fileID,'std(E[rA-rf]) %10.6f%%\n', 400*std_EexcA);
fprintf(fileID,'std(rA-rf])   %10.6f%%\n\n', 400*std_realized_excA);


fprintf(fileID,'\nDIVIDEND PRICE STD\n');
fprintf(fileID,'std(d/p) %10.6f%%\n', 100*std_divprice);
fprintf(fileID,'\nSMOOTHED DIVIDEND PRICE STD\n');
fprintf(fileID,'std(d/p) %10.6f%%\n', 100*std_divprice_smoothed);

fprintf(fileID,'\nAUTOCORR \n');
fprintf(fileID,'-----------------------\n\n');
fprintf(fileID,'ac(E[rA-rf])   %10.6f%%\n', 100*ac_excA);
fprintf(fileID,'ac(E[rf])      %10.6f%%\n', 100*ac_Erf);
fprintf(fileID,'ac(rf)         %10.6f%%\n', 100*ac_rf);
fprintf(fileID,'ac(d/p)        %10.6f%%\n', 100*ac_divprice);
fprintf(fileID,'\nAUTOCORR YEAR OVER YEAR  \n');
fprintf(fileID,'-----------------------\n\n');
fprintf(fileID,'ac(d/p)        %10.6f%%\n', 100*yoy_ac_divprice);
fprintf(fileID,'\nAUTOCORR SMOOTH \n');
fprintf(fileID,'ac(d/p) smooth %10.6f%%\n', 100*ac_divprice_smoothed);
fprintf(fileID,'\nAUTOCORR SMOOTH YEAR OVER YEAR \n');
fprintf(fileID,'ac(d/p) smooth %10.6f%%\n', 100*ac_yoy_divprice_smoothed);

fprintf(fileID,'\nWEALTH SHARE STD\n');
fprintf(fileID,'std(sa)       %10.6f%%\n', 100*std_sa);
fprintf(fileID,'std(sc)       %10.6f%%\n', 100*std_sc);
fprintf(fileID,'\n-----------------------\n\n');

