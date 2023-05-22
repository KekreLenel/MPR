% -------------------------------------------------------------------------
% calc_CampbellShiller.m: calculate campbell shiller decomposition 
% for levered and unlevered equity claim 
% -------------------------------------------------------------------------
% authors:         Rohan Kekre and Moritz Lenel
% for updates see: https://github.com/KekreLenel/MPR
% -------------------------------------------------------------------------

fprintf(fileID,'Monetary policy shock - effects on impact \n');
fprintf(fileID,'-----------------------\n');
fprintf(fileID,'log(inv)         %10.6fbp\n', 10000*(loginv_series(3)-loginv_series(2)));
fprintf(fileID,'log(c)           %10.6fbp\n', 10000*(logc_series(3)-logc_series(2)));
fprintf(fileID,'log(y)           %10.6fbp\n', 10000*(logy_series(3)-logy_series(2)));
fprintf(fileID,'log(excret)      %10.6fbp\n', 10000*(exc_retA_series(3)-exc_retA_series(2)));
fprintf(fileID,'-----------------------\n\n');

M_effects(:,ccc) = [ 10000*(loginv_series(3)-loginv_series(2)); 10000*(logc_series(3)-logc_series(2)) ;10000*(logy_series(3)-logy_series(2))];

% Unlevered Campbell-Shiller (capital excess returns)
steady_state_DP_ratio = divk_price_series(2);
rho_dp = 1/(1+steady_state_DP_ratio);
helper = rho_dp*ones(n_irf_periods-3,1);
helper = cumprod(helper);
realized_return= (rk_series(3)- rk_series(2));
cash_flow_news= log(divk_series(3)/divk_series(2)) + sum(helper(1:end-1).*log(divk_series(4:end-1)./divk_series(3:end-2)) );
rf_news= sum(helper(1:end-1).*((rf_series(4:end-1)) - (rf_series(2)) ));
exc_news= sum(helper(1:end-1).*(exc_ret_series(4:end-1) - exc_ret_series(2) ));

% Levered Campbell-Shiller (equity excess returns)
steady_state_DP_ratioA = div_price_series(2);
rho_dpA = 1/(1+steady_state_DP_ratioA);
helperA = rho_dpA*ones(n_periods-3,1);
helperA = cumprod(helperA);
realized_returnA = (rA_series(3)- rA_series(2));
cash_flow_newsA = log(div_series(3)/div_series(2)) + sum(helperA(1:end-1).*(log(div_series(4:end-1)./div_series(3:end-2)))) ;
rf_newsA =  sum(helperA(1:end-1).*(rf_series(4:end-1) - rf_series(2) ));
exc_newsA = sum(helperA(1:end-1).*(exc_retA_series(4:end-1) - exc_retA_series(2) ));
SUM = cash_flow_newsA - rf_newsA- exc_newsA;

CS_Decomposition(:,ccc) = [cash_flow_newsA/SUM, -rf_newsA/SUM, -exc_newsA/SUM];

fprintf(fileID,'Campbell Shiller Decomposition - not levered \n');
fprintf(fileID,'-----------------------\n\n');
fprintf(fileID,'ExcRet        %10.6f%%\n', 100*realized_return);
fprintf(fileID,'CF News       %10.6f%% (%10.6f%%)\n', 100*cash_flow_news,  100*cash_flow_news/realized_return);
fprintf(fileID,'rF News       %10.6f%% (%10.6f%%)\n', 100*rf_news, -100*rf_news/realized_return);
fprintf(fileID,'Ex News       %10.6f%% (%10.6f%%)\n', 100*exc_news, -100*exc_news/realized_return);
fprintf(fileID,'\n-----------------------\n');
fprintf(fileID,'SUM           %10.6f%% (%10.6f%%)\n', 100*cash_flow_news - 100*rf_news- 100*exc_news, ...
100*cash_flow_news/realized_return  -100*rf_news/realized_return -100*exc_news/realized_return);
fprintf(fileID,'\n-----------------------\n\n');

fprintf(fileID,'Campbell Shiller Decomposition - levered \n');
fprintf(fileID,'-----------------------\n\n');
fprintf(fileID,'ExcRet        %10.6f%%\n', 100*realized_returnA);
fprintf(fileID,'CF News       %10.6f%% (%10.6f%%)\n', 100*cash_flow_newsA,  100*cash_flow_newsA/realized_returnA);
fprintf(fileID,'rF News       %10.6f%% (%10.6f%%)\n', 100*rf_newsA, -100*rf_newsA/realized_returnA);
fprintf(fileID,'Ex News       %10.6f%% (%10.6f%%)\n', 100*exc_newsA, -100*exc_newsA/realized_returnA);
fprintf(fileID,'\n-----------------------\n');
fprintf(fileID,'SUM           %10.6f%% (%10.6f%%)\n', 100*cash_flow_newsA - 100*rf_newsA- 100*exc_newsA, ...
100*cash_flow_newsA/realized_returnA  -100*rf_newsA/realized_returnA -100*exc_newsA/realized_returnA);
fprintf(fileID,'\n-----------------------\n\n');
