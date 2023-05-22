function calc_unexpected_transition( nxt_mat, n_nxt, nxt_mat_2, sss,
                                       c_vec, l_aggr, q_current, q_l_current, infl,
                                       nom_i, share_vec, trans_shock_vec,
                                       s_nxt, param::mpr_economy)
    @unpack_mpr_economy param
    brent_delta = 5E-16
    min_cons_sav = 1E-08
    c_temp=zeros(n_I)
    v_temp = zeros(1, n_I)
    mc_temp = zeros(1, n_I)
    survivor_wealth_vec=zeros(1, n_I)
    savings_vec=zeros(n_I)
    q_l_nxt=zeros(1, n_I)

    # =============================================================================
    # CURRENT STATE ===============================================================
    # =============================================================================
    k_aggr = state_grid[idx_k, sss]

    # =============================================================================
    # CURRENT PERIOD GUESSES ======================================================
    # =============================================================================
    pi_current   = aalpha * (l_aggr/k_aggr).^(1.0-aalpha)    
    aggr_wealth  = k_aggr*((1-ddelta)*q_current + pi_current)   
    y_current    = k_aggr.^aalpha * l_aggr.^(1.0-aalpha)
    w_choice = (1.0-aalpha)*y_current/l_aggr

    # =============================================================================
    # NEXT PERIOD =================================================================
    # =============================================================================
    z_temp = trans_shock_vec[sidx_z]

    k_nxt      = nxt_mat_2

    # value functions
    for iii in 1 : n_I
        v_temp[:,iii] = nxt_mat[:,iii]
        mc_temp[:,iii] = nxt_mat[:,n_I+iii]
    end

    # tomorrow's price, labor choice, infl and price of time endowment
    q_nxt               = nxt_mat[:,2*n_I + 1]
    l_aggr_nxt          = nxt_mat[:,2*n_I + 2]
    infl_nxt            = nxt_mat[:,2*n_I + 3]
    q_l_nxt[:,1]        .= nxt_mat[:,2*n_I + 4] .* exp.(z_temp)
    q_l_nxt[:,2]        .= nxt_mat[:,2*n_I + 5] .* exp.(z_temp)
    q_l_nxt[:,3]        .= nxt_mat[:,2*n_I + 6] .* exp.(z_temp)


    # profit next period per unit of capital
    pi_nxt = @. aalpha*exp(z_temp*(1.0-aalpha))*(l_aggr_nxt^(1.0-aalpha))*k_nxt^(aalpha-1.0)

    # capital return, adjustment in disaster state
    rk_vec = ( (1.0 .- ddelta) .* q_nxt .+ pi_nxt) ./ q_current
    rf_vec = nom_i ./ infl_nxt

    #! =============================================================================
    #! STORE VALUE FUNCTION, NEXT CAPITAL, NEXT WEALTH SHARES ======================
    #! =============================================================================
    c_temp    = c_vec
    for iii = 1 : n_I
        # Read parameters for agent

        wealth = aggr_wealth .* wealth_share_grid[iii,sss]
        labor  = l_aggr .* labor_alloc_vec[iii]

        # Enforce minimal spending and mininal saving
        consumption =  minimum( [maximum([c_temp[iii], min_cons_sav]),
                        wealth + w_choice*labor-min_cons_sav] )

        # Get Next Period Wealth Share
        r_alpha_vec = share_vec[iii] .* rf_vec .+ (1.0 .- share_vec[iii]) .* rk_vec
        savings = wealth .+ w_choice .* labor .- consumption

        # update wealth transitions accordingly
        survivor_wealth_vec[:,iii] .= lmbd_vec[iii] .* savings .* r_alpha_vec
        savings_vec[iii] = savings
    end

    # (4.2) Next period wealth share
    temp_nxt_aggr_wealth = sum(sum(survivor_wealth_vec,dims=2))

    # The original code is sum(survivor_wealth_vec,2), I assume they add along the second row.
    # This part needs to be revisited

    for iii = 1: n_I
        s_nxt[iii] = (((1.0-xi)* survivor_wealth_vec[1,iii] + xi * s_bar_vec[iii] * temp_nxt_aggr_wealth) / (temp_nxt_aggr_wealth))
    end

    return s_nxt

end
