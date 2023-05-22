per_person_wealth = zeros(n_I)
c_ss_vec = zeros(n_I)
s_ss_vec = zeros(n_I)
# average preference parameters 
bbeta_avrg = sum(s_trgt_vec.*bbeta_vec) 
gma_avrg   = sum(s_trgt_vec.*gma_vec)
ies_avrg   = sum(s_trgt_vec.*ies_vec)


# corresponding returns
rk_ss = 1.0/bbeta_avrg  - 1.0 
rf_ss = 1.0/bbeta_avrg  - 1.0

# steady state labor
l_ss = l_target 

# Capital, wage, output in stationary economy
k_ss = l_ss * ((rk_ss + ddelta) / aalpha)^(1.0 / (aalpha - 1.0))
w_ss = (1.0 - aalpha) * (k_ss^aalpha) * (l_ss^(-aalpha))
y_ss = (k_ss^aalpha) * l_ss^(1-aalpha)

# Wealth shares, consumption
s_ss_vec = s_trgt_vec
for iii = 1:n_I
    per_person_wealth[iii] = s_ss_vec[iii]/lmbd_vec[iii]
    c_ss_vec[iii] = w_ss*l_target*labor_alloc_vec[iii] + per_person_wealth[iii] * rk_ss * k_ss
end

for iii = 1:n_I
    if (labor_alloc_vec[iii] > sqrt(eps())) 
        thtbar = 0.5
        diff = 1.0
        iter = 0
        while (diff > 1E-12)

            iter = iter + 1
            thtbar_new = ies_vec[iii]*w_ss/c_ss_vec[iii] * (l_ss*labor_alloc_vec[iii])^(-1.0/tht) * 
                    (1.0 + (1.0/ies_vec[iii]-1.0) * thtbar * ((l_ss*labor_alloc_vec[iii])^(1.0+1.0/tht)) / (1.0+1.0/tht))
            diff = abs(thtbar-thtbar_new)
            thtbar = thtbar + 0.5*(thtbar_new-thtbar)

            if (iter > 1000) 
                error( "ERROR: No convergence in finding thtbar.")
            end

        end
        thtbar_vec[iii] = thtbar 
    else
        thtbar_vec[iii] = 0.0
    end
end

# ----------------------------------------------------------------- #
# Value function; labor endowment prices
# ----------------------------------------------------------------- #
for iii = 1:n_I

    v_normalization_vec[iii] = (1-bbeta_vec[iii])
    
    # value function
    v_ss_vec[iii] = ( v_normalization_vec[iii]/(1.0-bbeta_vec[iii]) )^(1.0/(1.0 - 1.0/ies_vec[iii])) *  
                    c_ss_vec[iii]*(1+(1/ies_vec[iii]-1)*thtbar_vec[iii]*((l_ss*labor_alloc_vec[iii])^(1.0 + 1.0/tht))/(1.0 + 1.0/tht))^( (1.0/ies_vec[iii])/(1.0-1.0/ies_vec[iii]) )
    
    # labor endowment price
    q_l_ss_vec[iii] = w_ss * (l_ss*labor_alloc_vec[iii]) / (1.0-bbeta_vec[iii] )
end  

# total individual wealth including endowment claim
tot_wealth_ss_vec = q_l_ss_vec + per_person_wealth*k_ss*(1.0+rk_ss)

# ----------------------------------------------------------------- #
# Scaling parameters of wage and capital grid 
# ----------------------------------------------------------------- #
k_grid_mean   = k_grid_adj*k_ss
w_grid_mean   = w_grid_adj*w_ss

k_grid_dev   = k_dev_param*k_grid_mean
w_grid_dev   = w_dev_param*w_grid_mean

# ----------------------------------------------------------------- #
# Print SS
# ----------------------------------------------------------------- #
println("---------------------------------------------")
println( "STEADY STATE")
println("---------------------------------------------")
println(                     " k_ss   =            " , k_ss)
println(              " k_grid_mean   =            " , k_grid_mean)
println(                     " y_ss   =            " , y_ss)
println(                     " l_ss   =            " , l_ss)
println(                     " w_ss   =            " , w_ss)
println("")
println( "---------------------------------------------")
println( " q_l_ss =          " , q_l_ss_vec[1], " ", q_l_ss_vec[2], " ", q_l_ss_vec[3])
println( " thtbar =            " , thtbar_vec[1], " ", thtbar_vec[2], " ", thtbar_vec[3])
println( " c_ss =            " , c_ss_vec[1], " ", c_ss_vec[2], " ", c_ss_vec[3])
println( " v_ss =            " , v_ss_vec[1], " ", v_ss_vec[2], " ", v_ss_vec[3])
println("---------------------------------------------")
println("")

# Export calibrated thtbar for result files
save(results_folder*"extra_data.jld","thtbar_vec", thtbar_vec);
