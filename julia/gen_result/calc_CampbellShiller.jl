# Extract Campbell-Shiller

#  Unlevered Campbell-Shiller (capital excess returns)
global io *= @sprintf("Monetary policy shock - effects on impact \n");
global io *= @sprintf("-----------------------\n");
global io *= @sprintf("log(inv)         %10.6fbp\n", 10000*(loginv_series[3]-loginv_series[2]));
global io *= @sprintf("log(c)           %10.6fbp\n", 10000*(logc_series[3]-logc_series[2]));
global io *= @sprintf("log(y)           %10.6fbp\n", 10000*(logy_series[3]-logy_series[2]));
global io *= @sprintf("log(excret)           %10.6fbp\n", 10000*(exc_retA_series[3])-exc_retA_series[2]);
global io *= @sprintf("-----------------------\n\n");

M_effects[:,ccc] = [ 10000*(loginv_series[3]-loginv_series[2]); 10000*(logc_series[3]-logc_series[2]) ;10000*(logy_series[3]-logy_series[2])];

#  Method 1 (Individual)
steady_state_DP_ratio = divk_price_series[2];
rho_dp = @. 1/(1+steady_state_DP_ratio);
helper = rho_dp*ones(n_irf_periods-3,1);
helper = cumprod(helper,dims =1);
realized_return= (rk_series[3]- rk_series[2]);
cash_flow_news= log.(divk_series[3]/divk_series[2]) + sum(helper[1:end-1].*log.(divk_series[4:end-1]./divk_series[3:end-2]) );
rf_news= sum(helper[1:end-1].*((rf_series[4:end-1]) .- (rf_series[2]) ));
exc_news= sum(helper[1:end-1].*(exc_ret_series[4:end-1] .- exc_ret_series[2] ));

#  Levered Campbell-Shiller (equity excess returns)

#  Method 1 (Individual)
steady_state_DP_ratioA = div_price_series[2];
rho_dpA = @. 1/(1+steady_state_DP_ratioA);
helperA = rho_dpA*ones(n_periods-3,1);
helperA = cumprod(helperA,dims =1);
realized_returnA = (rA_series[3]- rA_series[2]);
cash_flow_newsA = log.(div_series[3]/div_series[2]) + sum(helperA[1:end-1].*(log.(div_series[4:end-1]./div_series[3:end-2]))) ;
rf_newsA =  sum(helperA[1:end-1].*(rf_series[4:end-1] .- rf_series[2]));
exc_newsA = sum(helperA[1:end-1].*(exc_retA_series[4:end-1] .- exc_retA_series[2] ));
SUM = cash_flow_newsA .- rf_newsA .- exc_newsA;
CS_Decomposition[:,ccc] = [cash_flow_newsA/SUM, -rf_newsA/SUM, -exc_newsA/SUM];



global io *= @sprintf("Campbell Shiller Decomposition - not levered \n");
global io *= @sprintf("-----------------------\n\n");
global io *= @sprintf("ExcRet           %10.6f%%\n", 100*realized_return);
global io *= @sprintf("CF News       %10.6f%% (%10.6f%%)\n", 100*cash_flow_news,  100*cash_flow_news/realized_return);
global io *= @sprintf("rF News       %10.6f%% (%10.6f%%)\n", 100*rf_news, -100*rf_news/realized_return);
global io *= @sprintf("Ex News       %10.6f%% (%10.6f%%)\n", 100*exc_news, -100*exc_news/realized_return);
global io *= @sprintf("\n-----------------------\n");
global io *= @sprintf("SUM           %10.6f%% (%10.6f%%)\n", 100*cash_flow_news - 100*rf_news- 100*exc_news, 
100*cash_flow_news/realized_return  -100*rf_news/realized_return -100*exc_news/realized_return);
global io *= @sprintf("\n-----------------------\n\n");

global io *= @sprintf("Campbell Shiller Decomposition - levered \n");
global io *= @sprintf("-----------------------\n\n");
global io *= @sprintf("ExcRet           %10.6f%%\n", 100*realized_returnA);
global io *= @sprintf("CF News       %10.6f%% (%10.6f%%)\n", 100*cash_flow_newsA,  100*cash_flow_newsA/realized_returnA);
global io *= @sprintf("rF News       %10.6f%% (%10.6f%%)\n", 100*rf_newsA, -100*rf_newsA/realized_returnA);
global io *= @sprintf("Ex News       %10.6f%% (%10.6f%%)\n", 100*exc_newsA, -100*exc_newsA/realized_returnA);
global io *= @sprintf("\n-----------------------\n");
global io *= @sprintf("SUM           %10.6f%% (%10.6f%%)\n", 100*cash_flow_newsA - 100*rf_newsA- 100*exc_newsA, 
100*cash_flow_newsA/realized_returnA  -100*rf_newsA/realized_returnA -100*exc_newsA/realized_returnA);
global io *= @sprintf("\n-----------------------\n\n");
