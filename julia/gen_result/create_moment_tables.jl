if (ccc == ix_bm || ccc == ix_idio_bm || ccc == ix_interm_bm)  

    # Moment Names
    targeted_names = ["\$\\sigma(\\Delta\\log c)\$" 
                      "\$\\sigma(\\Delta\\log x)\$" 
                      "\$4r_{+1}\$" 
                      "\$4[r^e_{+1}-r_{+1}]\$" 
                      "\$\\sigma(4\\mathbb{E}r_{+1})\$" 
                      "\$\\rho(\\mathbb{E}r_{+1})\$" 
                      "\$q k^a/a^a\$" 
                      "\$q k^c/a^c\$" 
                      "\$\\lambda^aa^a/\\sum_i\\lambda^ia^i\$" 
                      "\$\\lambda^ca^c/\\sum_i\\lambda^ia^i\$" 
                      "\$-\\sum_i\\lambda^ib^i/\\sum_i\\lambda^ia^i\$"]; 
    
    if ccc == ix_idio_bm 
    param_names = ["\$\\sigma^z\$" 
                   "\$\\chi^x\$" 
                   "\$\\beta\$" 
                   "\$\\gamma^a = \\gamma^b = \\gamma^c\$" 
                   "\$\\sigma^p\$"
                   "\$\\rho^p\$" 
                   "\$\\eta^b\$" 
                   "\$\\underline{k}\$" 
                   "\$\\xi\\bar{s}^a\$" 
                   "\$\\xi\\bar{s}^c\$" 
                   "\$b^g\$"]; 
    else 
    param_names = ["\$\\sigma^z\$" 
                   "\$\\chi^x\$" 
                   "\$\\beta\$" 
                   "\$\\gamma^b\$" 
                   "\$\\sigma^p\$" 
                   "\$\\rho^p\$" 
                   "\$\\gamma^a\$" 
                   "\$\\underline{k}\$" 
                   "\$\\xi\\bar{s}^a\$" 
                   "\$\\xi\\bar{s}^c\$" 
                   "\$b^g\$"]; 
    
    end
    
    if ccc == ix_idio_bm  
    param_desc = ["std. dev. prod." 
                  "capital adj. cost" 
                  "discount factor" 
                  "RRA" 
                  "std. dev. log dis. prob." 
                  "persist. log dis. prob." 
                  "idio. risk \$b\$" 
                  "lower bound \$k^i\$" 
                  "newborn endowment \$a\$" 
                  "newborn endowment \$c\$" 
                  "real value govt. bonds"];           
          
    param_value = [100*prm_vec[param_id].sig_z 
                    prm_vec[param_id].chiX 
                    prm_vec[param_id].bbeta_a 
                    prm_vec[param_id].gma_b 
                    prm_vec[param_id].disast_std 
                    prm_vec[param_id].rho_p 
                    prm_vec[param_id].idio_risk_c^2 
                    prm_vec[param_id].kbar 
                    100*prm_vec[param_id].xi*prm_vec[param_id].s_bar_a 
                    100*prm_vec[param_id].xi*prm_vec[param_id].s_bar_c 
                    -mean(bg_series./z_series)];
    else 
    param_desc = ["std. dev. prod." 
                  "capital adj. cost" 
                  "discount factor" 
                  "RRA \$b\$" 
                  "std. dev. log dis. prob." 
                  "persist. log dis. prob." 
                  "RRA \$a\$" 
                  "lower bound \$k^i\$" 
                  "newborn endowment \$a\$" 
                  "newborn endowment \$c\$" 
                  "real value govt. bonds"];   

    param_value = [100*prm_vec[param_id].sig_z 
                    prm_vec[param_id].chiX 
                    prm_vec[param_id].bbeta_a 
                    prm_vec[param_id].gma_b 
                    prm_vec[param_id].disast_std 
                    prm_vec[param_id].rho_p 
                    prm_vec[param_id].gma_a 
                    prm_vec[param_id].kbar 
                    100*prm_vec[param_id].xi*prm_vec[param_id].s_bar_a 
                    100*prm_vec[param_id].xi*prm_vec[param_id].s_bar_c 
                    -mean(bg_series./z_series)];
    end
              
    # Moment Targets
    if (ccc == ix_idio_bm)  
    targeted_values = [0.5, 2.1, 1.3, 7.3, 2.2, 0.79, 2.0, 1.1, 18, 23, -10];
    elseif (ccc == ix_interm_bm)  
    targeted_values = [0.5, 2.1, 1.3, 7.3, 2.2, 0.79, 4.4, 1.1, 2, 23, -10];
    else
    targeted_values = [0.5, 2.1, 1.3, 7.3, 2.2, 0.79, 2.0, 1.1, 18, 23, -10];
    end              
    
    # Moment Values
    targeted_moments    = zeros(12,1);
    targeted_moments[1] = 100*std(log.(c_series[3:end-1]./ c_series[2:end-2]));
    targeted_moments[2] = 100*std(log.(inv_series[3:end-1]./ inv_series[2:end-2]));
    targeted_moments[3] = 400*mean(exp.(rf_series[2:end-1]).-1);
    targeted_moments[4] = 400*mean(exp.(rA_series[2:end-1])) - 400*mean(exp.(rf_series[2:end-1]));
    targeted_moments[5] = 400*std(E_rf_series[2:end-1]);
    targeted_moments[6] = cor( E_rf_series[3:end-1], E_rf_series[2:end-2]);
    targeted_moments[7] = (1-mean(sh_a_adj_series[2:end-1]));
    targeted_moments[8] = (1-mean(sh_c_adj_series[2:end-1]));
    targeted_moments[9] = 100*mean(s_a_series_adj[2:end-1]);
    targeted_moments[10] = 100*mean(s_c_series_adj[2:end-1]);
    targeted_moments[11] = -100*prm_vec[param_id].gov_debt; 
    

    # Write Targeted Moments
    T = ""
    T *= @sprintf("& Description & Value & Moment & Target & Model \\\\ \n");
    T *= @sprintf("\\hline \n");
    T *= @sprintf("%s & %s & %6.1f\\%% & %s & %6.1f\\%% & %6.1f\\%% \\\\ \n", param_names[1], param_desc[1], param_value[1], targeted_names[1],  targeted_values[1],  targeted_moments[1]);
    T *= @sprintf("%s & %s & %6.1f & %s & %6.1f\\%% & %6.1f\\%% \\\\ \n", param_names[2], param_desc[2], param_value[2], targeted_names[2],  targeted_values[2],  targeted_moments[2]);
    T *= @sprintf("%s & %s & %6.3f & %s & %6.1f\\%% & %6.1f\\%% \\\\ \n", param_names[3], param_desc[3], param_value[3], targeted_names[3],  targeted_values[3],  targeted_moments[3]);
    T *= @sprintf("%s & %s & %6.0f & %s & %6.1f\\%% & %6.1f\\%% \\\\ \n", param_names[4], param_desc[4], param_value[4], targeted_names[4],  targeted_values[4],  targeted_moments[4]);
    T *= @sprintf("%s & %s & %6.2f\\%% & %s & %6.1f\\%% & %6.1f\\%% \\\\ \n", param_names[5], param_desc[5], param_value[5], targeted_names[5],  targeted_values[5],  targeted_moments[5]);
    T *= @sprintf("%s & %s & %6.2f & %s & %6.2f & %6.2f \\\\ \n", param_names[6], param_desc[6], param_value[6],   targeted_names[6],  targeted_values[6],  targeted_moments[6]);
    T *= @sprintf("%s & %s & %6.0f & %s & %6.1f & %6.1f \\\\ \n", param_names[7], param_desc[7], param_value[7],   targeted_names[7],  targeted_values[7],  targeted_moments[7]);
    T *= @sprintf("%s & %s & %6.1f\$\\cdot k^{ss}\$ & %s & %6.1f & %6.1f \\\\ \n", param_names[8], param_desc[8], param_value[8],   targeted_names[8],  targeted_values[8],  targeted_moments[8]);
    T *= @sprintf("%s & %s & %6.0f\\%% & %s & %6.1f\\%% & %6.1f\\%% \\\\ \n", param_names[9], param_desc[9], param_value[9], targeted_names[9],  targeted_values[9],  targeted_moments[9]);
    T *= @sprintf("%s & %s & %6.0f\\%% & %s & %6.1f\\%% & %6.1f\\%% \\\\ \n", param_names[10], param_desc[10], param_value[10], targeted_names[10], targeted_values[10], targeted_moments[10]);
    T *= @sprintf("%s & %s & See Below & %s & %6.1f\\%% & %6.1f\\%% \\\\ \n", param_names[11], param_desc[11],   targeted_names[11], targeted_values[11], targeted_moments[11]);
    # T *= @sprintf("%s & %s & %6.1f\\%% & %s & %6.2f\\%% & %6.2f\\%% \\\\ \n", param_names(12), param_desc(12), param_value(12), targeted_names(12), targeted_values(12), targeted_moments(12));
    T *= @sprintf("\\hline \n");
    T *= @sprintf("\\multicolumn{6}{@{}l}{\\scriptsize Disutility parameters \$\\\bar{\\nu}^i\$ are set to \$(%6.2f, %6.2f, %6.2f)\$ to jointly match the average labor and the steady state labor income shares.} \\\\ \n", mpr.thtbar_vec[1], mpr.thtbar_vec[2], mpr.thtbar_vec[3]);
    T *= @sprintf("\\hline \n");
    write(tab_path*"Targeted_Moments_"*string(ccc)*".tex", T);
end
    # 6. Write Untargeted Moments (Paper Table 6)
    # Moment Names

if ccc == ix_bm
    untargeted_names = ["\$\\sigma(\\Delta\\log y)\$"
                        "\$\\sigma(\\Delta\\log \\ell)\$"
                        "\$\\sigma(d/p)\$"
                        "\$\\sum_i\\lambda^impr^i\$"
                        "\$mpr^a\$"
                        "\$mpr^b\$"
                        "\$mpr^c\$"
                        "\$\\sum_i\\lambda^impc^i\$"
                        "\$mpc^a\$"
                        "\$mpc^b\$"
                        "\$mpc^c\$"];
                        # "\$\\sigma(\\mathbb{E}[r_{+1}])\$"

    # Data Values
    untargeted_data = [ 0.8, 0.8,  0.2, 0.2, 0.2];

    # Simulated Values
    untargeted_moments = zeros(11,1);
    untargeted_moments[1] = 100*std(log.(y_series[3:end-1]./ y_series[2:end-2]));
    untargeted_moments[2] = 100*std(log.(l_series[3:end-1]./ l_series[2:end-2]));
    untargeted_moments[3] = 100*std(div_price_smoothed_series[2:end-1])/4;
    untargeted_moments[4] = mean(prm_vec[param_id].lmbd_a*MPR1 .+ prm_vec[param_id].lmbd_b*MPR2 .+ prm_vec[param_id].lmbd_c*MPR3);
    untargeted_moments[5] = mean(MPR1);
    untargeted_moments[6] = mean(MPR2);
    untargeted_moments[7] = mean(MPR3);
    untargeted_moments[8] = mean(prm_vec[param_id].lmbd_a*MPC1 .+ prm_vec[param_id].lmbd_b*MPC2 .+ prm_vec[param_id].lmbd_c*MPC3);
    untargeted_moments[9] = mean(MPC1);
    untargeted_moments[10] = mean(MPC2);
    untargeted_moments[11] = mean(MPC3);

    # Write Unargeted Moments

    T = ""
    T *= @sprintf("Moment & Data & Model \\\\ \n");
    T *= @sprintf("\\hline \n");
    T *= @sprintf("%s & %6.1f\\%% & %6.1f\\%% \\\\ \n", untargeted_names[1],  untargeted_data[1],  untargeted_moments[1]);
    T *= @sprintf("%s & %6.1f\\%% & %6.1f\\%% \\\\ \n", untargeted_names[2],  untargeted_data[2],  untargeted_moments[2]);
    T *= @sprintf("%s & %6.1f\\%% & %6.1f\\%% \\\\ \n", untargeted_names[3],  untargeted_data[3],  untargeted_moments[3]);
    T *= @sprintf("\\hline \n");
    T *= @sprintf("%s & \\approx %6.1f\\%% & %6.1f\\%% \\\\ \n", untargeted_names[4],  untargeted_data[4],  untargeted_moments[4]);
    T *= @sprintf("%s & & %6.1f  \\\\ \n", untargeted_names[5], untargeted_moments[5]);
    T *= @sprintf("%s & & %6.1f \\\\ \n", untargeted_names[6], untargeted_moments[6]);
    T *= @sprintf("%s & & %6.1f \\\\ \n", untargeted_names[7], untargeted_moments[7]);
    T *= @sprintf("\\hline \n");
    T *= @sprintf("%s & \\approx %6.1f\\%% & %6.2f\\%% \\\\ \n", untargeted_names[8],  untargeted_data[5],  untargeted_moments[8]);
    T *= @sprintf("%s & & %6.2f \\\\ \n", untargeted_names[9], untargeted_moments[9]);
    T *= @sprintf("%s & & %6.2f \\\\ \n", untargeted_names[10], untargeted_moments[10]);
    T *= @sprintf("%s & & %6.2f \\\\ \n", untargeted_names[11], untargeted_moments[11]);
    T *= @sprintf("\\hline \n");
    write(tab_path*"Untargeted_Moments_"*string(ccc)*".tex", T);
end
