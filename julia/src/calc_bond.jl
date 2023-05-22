using Revise
using UnPack, LinearAlgebra, StructArrays
using Interpolations
using ForwardDiff

function calc_bond_prices(nxt_mat, n_nxt, nxt_mat_2, nxt_bond_prices, sss, c_vec, l_aggr,
    q_current,q_l_current,infl,nom_i, share_vec, constraint_binding,param)
    #used input: results_vec
    @unpack_mpr_economy param
    min_cons_sav = 1E-08
    n_wealth = 11
    nom_i_grid = zeros(n_wealth)
    saving_grid = zeros(n_wealth)
    capital_grid = zeros(n_wealth)
    bond_grid = zeros(n_wealth)
    q_grid = zeros(n_wealth)
    E_rk_grid=zeros(n_wealth)
    E_rf_grid=zeros(n_wealth)
    v_temp = zeros(n_quad,n_I)
    mc_temp = zeros(n_quad,n_I)
    n_wealth = 11; sqrt_eps = sqrt(eps())
    ratio_grid = collect(LinRange(-1.0, 1.0, n_wealth))
    q_inc = collect(LinRange(0.0025, -0.0025, n_wealth))
    nom_i_inc = collect(LinRange(-0.0025, 0.0025, n_wealth))
    cons_grid = similar(ratio_grid); share_grid = similar(ratio_grid)
    wealth_grid = similar(ratio_grid)
    dc_dErk_mat = zeros(n_I); 
    db_dErk_mat = zeros(n_I); 
    dk_dErk_mat= zeros(n_I);
    dc_dErf_mat = zeros(n_I); 
    db_dErf_mat = zeros(n_I); 
    dk_dErf_mat= zeros(n_I);
    savings_vec = zeros(n_I);
    v_store = zeros(n_I);
    mc_store = zeros(n_I);
    MPC_mat = zeros(n_I);
    MPS_mat = zeros(n_I);
    MPK_mat = zeros(n_I);

    v_idio =  zeros(n_quad*n_idio, n_I)
    mc_idio =  zeros(n_quad*n_idio, n_I)
    q_l_nxt_idio = zeros(n_quad*n_idio, n_I)
    dz_vec_idio = zeros(n_quad*n_idio)
    rf_vec_idio = zeros(n_quad*n_idio)
    big_weight_vec_idio = zeros(n_quad*n_idio)
    r_alpha_vec_idio = zeros(n_quad*n_idio)


    # same as in equilibrium calculation
    k_aggr = state_grid[idx_k, sss] # current capital stock, disaster prob & wage
    p_dis = exp(state_grid[idx_dis,sss])

    big_weight_vec = zeros(n_quad) # initialize
    big_weight_vec[1:n_quad-1] .= quad_weight_vec .* (1.0-p_dis)
    big_weight_vec[n_quad] = p_dis

    pi_current = aalpha * (l_aggr/k_aggr)^(1.0-aalpha)
    aggr_wealth = k_aggr * (pi_current + (1.0-ddelta)*q_current)
    y_current = k_aggr^aalpha * l_aggr^(1.0-aalpha)
    w_choice = (1.0-aalpha)*y_current/l_aggr 

    k_nxt = copy(nxt_mat_2)
    for iii in 1:n_I
        v_temp[:,iii] .= nxt_mat[:,iii]
        mc_temp[:,iii] .= nxt_mat[:,n_I+iii]
    end
    q_nxt = nxt_mat[:,2*n_I+1]
    l_aggr_nxt = nxt_mat[:,2*n_I+2]
    infl_nxt = nxt_mat[:,2*n_I+3]
    q_l_nxt = nxt_mat[:,(2*n_I+4):(2*n_I+6)] .* exp.(dz_vec) # more convenient compared to original fortran code

    pi_nxt = aalpha*exp.( (dz_vec).*(1.0-aalpha)).*(l_aggr_nxt.^(1.0-aalpha)).*k_nxt.^(aalpha-1.0)

    # gross capital return (adjust in disaster state)
    rk_vec =( (1.0 .- ddelta) .* q_nxt .+ pi_nxt ) ./ q_current
    rk_vec[n_quad] *= exp(dz_vec[n_quad])

    # Gross risk free return
    rf_vec = nom_i ./ infl_nxt

    # basic idiosyncratic risk expansion
    if (use_idio_risk == 1)
        for kkk in 1:n_idio
            for iii in 1:n_I
                v_idio[ ((kkk-1)*n_quad+1):(kkk*n_quad),iii].=v_temp[:,iii]
                mc_idio[ ((kkk-1)*n_quad+1):(kkk*n_quad),iii].=mc_temp[:,iii]
                q_l_nxt_idio[ ((kkk-1)*n_quad+1):(kkk*n_quad),iii] .= q_l_nxt[:,iii].*exp(idio_shock_grid[kkk,iii])
            end
            dz_vec_idio[((kkk-1)*n_quad+1):(kkk*n_quad)] .= dz_vec
            rf_vec_idio[((kkk-1)*n_quad+1):(kkk*n_quad)] .= rf_vec
            big_weight_vec_idio[((kkk-1)*n_quad+1):(kkk*n_quad)] .= big_weight_vec .* idio_weight_vec[kkk]
        end
    end

    # calculate pricing kernel
    M_vec_mat =zeros(n_quad,n_I)
    # marg_util_c_vec = zeros(n_I)

    for iii in 1:n_I
        # (7.1) Read parameters and quantities (same as in part (4) above)
        ies   = ies_vec[iii]
        bbeta = bbeta_vec[iii]
        gma   = gma_vec[iii]
        thtbar  = thtbar_vec[iii]
        # v_normalization = v_normalization_vec[iii]

        wealth = aggr_wealth .* wealth_share_grid[iii,sss]
        labor  = l_aggr .* labor_alloc_vec[iii]

        consumption =  min(max(c_vec[iii], min_cons_sav),
                          (wealth + w_choice *labor -min_cons_sav) )
        r_alpha_vec = share_vec[iii] .* rf_vec .+ (1.0 .- share_vec[iii]) .* rk_vec
        savings     = wealth .+ w_choice .* labor .- consumption
        savings_vec[iii] = savings

        # (7.2) Calculate (Realized) Pricing Kernel
        if (use_idio_risk == 1)
            # ------------------------------------------------------------ !
            # With Idiosyncratic Risk
            # ------------------------------------------------------------ !
            # Expand portfolio return
            for kkk in 1 : n_idio
                r_alpha_vec_idio[((kkk-1)*n_quad+1):(kkk*n_quad)] = share_vec[iii].*rf_vec .+ (1.0 - share_vec[iii]).*rk_vec.*exp.(idio_shock_grid[kkk,iii])
            end
           # Same question as before

            # Other Operations
            next_period_share_idio = (savings .* r_alpha_vec_idio .+ q_l_nxt_idio[:,iii]) ./ (tot_wealth_ss_vec[iii])
            v_vec_twisted_idio = (v_idio[:,iii] .* next_period_share_idio) .^ (1.0 .- gma)
            EV = sum(big_weight_vec_idio .* v_vec_twisted_idio)

            util, util_c_deriv, labor_part = util_fun(consumption, labor, iii, thtbar_vec,  ies_vec, tht)

            #Need to be revisited

            objf = ( (1.0-bbeta) .* util .+ bbeta .* (EV .^( (1.0 .- (1.0 ./ ies)) ./ (1.0 .- gma)))) .^ (1.0 ./ (1.0-1.0 ./ ies))
            
            mc_store[iii] = util_c_deriv
            v_store[iii] = objf*( (wealth + q_l_current[iii])/tot_wealth_ss_vec[iii])^(-1)

            # SDF: first, the full form
            M_idio = bbeta* ( next_period_share_idio .^ (-gma) ) .* (EV .^( (gma - (1.0/ies))/(1.0-gma) )) .* (v_idio[:,iii].^((1.0/ies) -gma )) .* mc_idio[:,iii] ./ util_c_deriv 

            # Then for each of the n_quad aggregate nodes, average across the n_idio Idiosyncratic nodes
            M_vec_mat[:,iii]  .= 0.0
            for kkk in 1 : n_idio
                M_vec_mat[:,iii] .= M_vec_mat[:,iii] .+ M_idio[((kkk-1) .* n_quad .+1):(kkk .* n_quad)] .* idio_weight_vec[kkk]
            end

        else
            # ------------------------------------------------------------ !
            # No Idiosyncratic Risk
            # ------------------------------------------------------------ !
            next_period_share = (savings .* r_alpha_vec .+ q_l_nxt[:,iii]) ./ (tot_wealth_ss_vec[iii])
            v_vec_twisted = (v_temp[:,iii] .* next_period_share) .^ (1.0 .- gma)
            EV = sum(big_weight_vec .* v_vec_twisted)

            util, util_c_deriv, labor_part = util_fun(consumption, labor, iii,thtbar_vec,  ies_vec, tht)

            objf = ( (1.0-bbeta) .* util .+  bbeta .*  (EV .^( (1.0 .- (1.0 ./ ies))/(1.0 .- gma)))) .^( 1.0 ./ (1.0 .- 1.0 ./ies))

            mc_store[iii] = util_c_deriv
            v_store[iii] = objf*( (wealth + q_l_current[iii])/tot_wealth_ss_vec[iii])^(-1)

            M_vec_mat[:,iii] .=  bbeta* ( next_period_share.^ (-gma) ) .* (EV .^( (gma - (1.0/ies))/(1.0-gma) )) .* (v_temp[:,iii].^((1.0/ies) -gma )) .* mc_temp[:,iii] ./ util_c_deriv 
        end


    end

    # This directly uses matrix multiplication to find the max and arg max
    # big_weight_vec is of size (n_quad, 1), M_vec_mat is size (n_quad,n_I)
    # nxt_bond_prices, infl_nxt are size (n_quad, 1)
    bond_price_prelim = vec(M_vec_mat' * (nxt_bond_prices./infl_nxt.*big_weight_vec)) # a three element vector
    bond_prices_out, who_prices_bond = findmax(bond_price_prelim)

    E_rf = dot(big_weight_vec, rf_vec); E_rk = dot(big_weight_vec, rk_vec)

    # Calculate MPC and MPR: no labor endowment trades
    for iii in 1:n_I
        ies   = ies_vec[iii]
        bbeta = bbeta_vec[iii]
        gma   = gma_vec[iii]
        thtbar  = thtbar_vec[iii]

        wealth = aggr_wealth .* wealth_share_grid[iii,sss]
        labor  = l_aggr .* labor_alloc_vec[iii]
        labor_part =  1.0 .+ (1.0 ./ ies .- 1.0) .* thtbar .* tht ./ (1.0 .+ tht) .* labor .^((1.0 .+ tht) ./ tht)
        # initialize
        wealth_grid .= ratio_grid.+ wealth 
        cons_grid .= c_vec[iii]
        share_grid .= share_vec[iii]
        # loop update

        for nnn in 1:n_wealth
            diff_inner = 1.0
            while (diff_inner > 1E-5)
                # update consump and portfolio choices
                wealth = wealth_grid[nnn]
                curr_share = share_grid[nnn]
                r_alpha_vec = curr_share .* rf_vec .+ (1. - curr_share).*rk_vec
                consumption = cons_grid[nnn]
                savings = wealth + w_choice * labor - consumption

                if (use_idio_risk==1)
                    for kkk in 1 : n_idio
                        r_alpha_vec_idio[((kkk-1)*n_quad+1):(kkk*n_quad)] .=  curr_share.*rf_vec .+ (1.0 - curr_share).*rk_vec.*exp.(idio_shock_grid[kkk,iii])
                    end
                   # Same question as before
     
                    # Other Operations
                    next_period_share_idio = (savings .* r_alpha_vec_idio .+ q_l_nxt_idio[:,iii]) ./ (tot_wealth_ss_vec[iii])
                    v_vec_twisted_idio = (v_idio[:,iii] .* next_period_share_idio) .^ (1.0 .- gma)
                    EV = sum(big_weight_vec_idio .* v_vec_twisted_idio)

                    deriv_helper_idio = bbeta .* (v_idio[:,iii].^((1.0/ies) - gma )) .* (next_period_share_idio.^(-gma)).* (EV.^( (gma - (1.0/ies))/(1.0-gma) )) .* mc_idio[:,iii]

                    if (constr_agents[iii] == 1 && constraint_binding[iii] == 1)
                        temp = sum(deriv_helper_idio .* rf_vec_idio .* big_weight_vec_idio)
                    else
                        temp = sum(deriv_helper_idio .* r_alpha_vec_idio .* big_weight_vec_idio)
                    end
                else # No Idiosyncratic risk
                    next_period_share = ( savings .* r_alpha_vec .+ q_l_nxt[:,iii]) ./ (tot_wealth_ss_vec[iii])
                    v_vec_twisted = (v_temp[:,iii] .* next_period_share) .^(1.0 .- gma)
                    EV = sum(big_weight_vec .* v_vec_twisted)


                    deriv_helper = bbeta .* (v_temp[:,iii].^((1.0/ies) - gma )) .* (next_period_share.^(-gma)).* (EV.^( (gma - (1.0/ies))/(1.0-gma) )) .* mc_temp[:,iii]

                    if (constr_agents[iii] == 1 && constraint_binding[iii] == 1)
                        temp = sum(deriv_helper .* rf_vec .* big_weight_vec)
                    else
                        temp = sum(deriv_helper .* r_alpha_vec .* big_weight_vec)
                    end
                end

                cons_update = labor_part .*  temp.^(-ies)
                cons_update = min(max(cons_update, min_cons_sav), (wealth .+ w_choice .*labor .-min_cons_sav) )

                # update portfolio choice

                if (constr_agents[iii] == 1 && constraint_binding[iii] == 1)
                    share_update = 1.0 - kbar_constr[iii]*q_current/savings
                else
                    # find natural leverage limit to prevent negative payoffs
                    temp_vec = @. -rk_vec/(rf_vec-rk_vec + sqrt_eps )
                    temp_vec[temp_vec.>=0.0] .= -1E4

                    share_low = max(low_guess_fixed, (maximum(temp_vec)+sqrt_eps))

                    # find natural leverage limit to prevent negative payoffs
                    temp_vec = @. -rk_vec/(rf_vec-rk_vec + sqrt_eps )
                    temp_vec[temp_vec.<=0.0] .= 1E4
                    share_high = min(high_guess_fixed, (minimum(temp_vec)-sqrt_eps))

                    share_update = curr_share
                    share_update = calc_portfolio_share(1.0, share_low, share_high, iii, big_weight_vec, 
                    savings, labor, cons_update, 0.0, q_l_nxt[:,iii], 
                    v_temp[:,iii],mc_temp[:,iii], rf_vec, rk_vec, q_current, share_update, param::mpr_economy)
                end

                # Update values
                diff_inner = @. abs(consumption - cons_update)
                cons_grid[nnn] = @. consumption + 0.5*(cons_update - consumption)
                share_grid[nnn] = share_update
            end
        end

        # Compute MPC, MPR at center n_GH_points
        saving_grid .= wealth_grid .+ w_choice .* labor .- cons_grid
        wealth = aggr_wealth * wealth_share_grid[iii,sss]

        interp = interpolate(wealth_grid, cons_grid, FritschCarlsonMonotonicInterpolation()) 
        MPC_mat[iii] = ForwardDiff.derivative(extrapolate(interp, Flat()), wealth)

        interp = interpolate(wealth_grid, saving_grid, FritschCarlsonMonotonicInterpolation())
        MPS_mat[iii] = ForwardDiff.derivative(interp, wealth)

        capital_grid .= @. (1. - share_grid) * saving_grid 
        interp = interpolate(wealth_grid, capital_grid, FritschCarlsonMonotonicInterpolation())
        MPK_mat[iii] = ForwardDiff.derivative(extrapolate(interp, Flat()), wealth)

        if (constr_agents[iii] == 1 && constraint_binding[iii] == 1)
            MPK_mat[iii] = 0.0
        end
        
    end
    
    # ============================================================ !
    # Derivatives to Price of Capital
    # ============================================================ !
    for iii=1 : n_I
        ies   = ies_vec[iii]
        bbeta = bbeta_vec[iii]
        gma   = gma_vec[iii]
        thtbar  = thtbar_vec[iii]

        wealth = aggr_wealth .* wealth_share_grid[iii,sss]
        labor  = (l_aggr .* labor_alloc_vec[iii])
        labor_part =  1.0 .+ (1.0/ies - 1.0) .* thtbar .* tht ./ (1.0 .+ tht) .* labor .^ ( (1.0 .+ tht) ./ tht)
        
        # Initialize
        q_grid .= q_current .+ q_inc
        cons_grid .= c_vec[iii]
        share_grid .= share_vec[iii]
        # Loop Update
        for nnn = 1 : n_wealth
            diff_inner = 1.0
             while(diff_inner > 1E-05)
                # ============================================================ !
                # Update Consumption Based on Portfolio Choice
                # ============================================================ !
                curr_q = q_grid[nnn]
                curr_share = share_grid[nnn]

                rk_vec_use = ( (1.0 .- ddelta) .* q_nxt .+ pi_nxt ) ./ curr_q
                rk_vec_use[n_quad] = exp.(dz_vec[n_quad]) .* rk_vec_use[n_quad]

                r_alpha_vec = curr_share .* rf_vec .+ (1.0 .- curr_share) .* rk_vec_use
                consumption = cons_grid[nnn]
                savings = wealth + w_choice*labor - consumption

                if (use_idio_risk==1)
                    for kkk in 1 : n_idio
                        r_alpha_vec_idio[((kkk-1)*n_quad+1):(kkk*n_quad)] .= curr_share.*rf_vec .+ (1.0 .- curr_share).*rk_vec_use.*exp.(idio_shock_grid[kkk,iii])
                    end
                   # Same question as before
     
                    # Other Operations
                    next_period_share_idio = (savings .* r_alpha_vec_idio .+ q_l_nxt_idio[:,iii]) ./ (tot_wealth_ss_vec[iii] )
                    v_vec_twisted_idio = (v_idio[:,iii] .* next_period_share_idio) .^ (1.0 .- gma)
                    EV = sum(big_weight_vec_idio .* v_vec_twisted_idio)

                    deriv_helper_idio =  bbeta .* (v_idio[:,iii].^((1.0/ies) - gma )) .* (next_period_share_idio.^(-gma)).* (EV.^( (gma - (1.0/ies))/(1.0-gma) )) .* mc_idio[:,iii]


                    if (constr_agents[iii] == 1 && constraint_binding[iii] == 1)
                        temp = sum(deriv_helper_idio .* rf_vec_idio .* big_weight_vec_idio)
                    else
                        temp = sum(deriv_helper_idio .* r_alpha_vec_idio .* big_weight_vec_idio)
                    end
                else # that is no idio risk
                    next_period_share = ( savings .* r_alpha_vec .+ q_l_nxt[:,iii]) ./ (tot_wealth_ss_vec[iii])
                    v_vec_twisted = (v_temp[:,iii] .* next_period_share) .^(1.0 .- gma)
                    EV = sum(big_weight_vec .* v_vec_twisted)


                    deriv_helper = bbeta .* (v_temp[:,iii].^((1.0/ies) - gma )) .* (next_period_share.^(-gma)).* (EV.^( (gma - (1.0/ies))/(1.0-gma) )) .* mc_temp[:,iii]

                    if (constr_agents[iii] == 1 && constraint_binding[iii] == 1)
                        temp = sum(deriv_helper .* rf_vec .* big_weight_vec)
                    else
                        temp = sum(deriv_helper .* r_alpha_vec .* big_weight_vec)
                    end
                end
                cons_update = labor_part .* (temp.^(-ies) ) 
                cons_update = min(max(cons_update, min_cons_sav), (wealth .+ w_choice .*labor .-min_cons_sav) )

                # update portfolio choice

                if (constr_agents[iii] == 1 && constraint_binding[iii] == 1)
                    share_update = 1.0 - kbar_constr[iii]*k_ss*q_current/savings
                else
                    # find natural leverage limit to prevent negative payoffs
                    temp_vec = @. -rk_vec_use/(rf_vec-rk_vec_use + sqrt_eps )
                    temp_vec[temp_vec.>=0.0] .= -1E4

                    share_low = max(low_guess_fixed, (maximum(temp_vec)+sqrt_eps))

                    # find natural leverage limit to prevent negative payoffs
                    temp_vec = @. -rk_vec_use/(rf_vec-rk_vec_use + sqrt_eps )
                    temp_vec[temp_vec.<=0.0] .= 1E4
                    share_high = min(high_guess_fixed, (minimum(temp_vec)-sqrt_eps))

                    share_update = curr_share
                    share_update = calc_portfolio_share(1.0, share_low, share_high, iii, big_weight_vec, 
                    savings, labor, cons_update, 0.0, q_l_nxt[:,iii], 
                    v_temp[:,iii],mc_temp[:,iii], rf_vec, rk_vec_use, q_current, share_update, param::mpr_economy)
                end

                # Update values
                diff_inner = @. abs(consumption - cons_update)
                cons_grid[nnn] = consumption .+ 0.5.*(cons_update .- consumption)
                share_grid[nnn] = share_update
            end
        end
        # Compute MPC, MPR at the center point


        saving_grid .= wealth + w_choice * labor .- cons_grid
        E_rk_grid .= E_rk .* q_current ./ q_grid

        interp = interpolate(E_rk_grid, cons_grid, FritschCarlsonMonotonicInterpolation()) 
        dc_dErk_mat[iii] = ForwardDiff.derivative(interp, E_rk)

        bond_grid .= share_grid .* saving_grid
        interp = interpolate(E_rk_grid, bond_grid, FritschCarlsonMonotonicInterpolation())
        db_dErk_mat[iii] = ForwardDiff.derivative(extrapolate(interp, Flat()), E_rk)

        capital_grid .= @. (1. - share_grid) * saving_grid / q_current
        interp = interpolate(E_rk_grid, capital_grid, FritschCarlsonMonotonicInterpolation())
        dk_dErk_mat[iii] = ForwardDiff.derivative(extrapolate(interp, Flat()), E_rk)

        if (constr_agents[iii] == 1 && constraint_binding[iii] == 1)
            dk_dErk_mat[iii] = 0.0
        end
    end

    # ============================================================ !
    # Derivatives to Nominal Rate
    # ============================================================ !
    for iii=1 : n_I
        ies   = ies_vec[iii]
        bbeta = bbeta_vec[iii]
        gma   = gma_vec[iii]
        thtbar  = thtbar_vec[iii]

        wealth = aggr_wealth*wealth_share_grid[iii,sss]
        labor  = (l_aggr*labor_alloc_vec[iii])
        labor_part =  1.0 + (1.0/ies - 1.0) * thtbar*tht/(1.0 + tht) * labor^( (1.0+tht)/tht)

        # Initialize
        nom_i_grid .= nom_i .+ nom_i_inc
        cons_grid .= c_vec[iii]
        share_grid .= share_vec[iii]

        for nnn = 1 : n_wealth
            diff_inner = 1.0
            while(diff_inner > 1E-05)
                # ============================================================ !
                # Update Consumption Based on Portfolio Choice
                # ============================================================ !
                curr_nom_i = nom_i_grid[nnn]
                curr_share = share_grid[nnn]

                rf_vec_use = curr_nom_i ./ infl_nxt

                r_alpha_vec = curr_share .* rf_vec_use .+ (1.0 .- curr_share) .* rk_vec
                consumption = cons_grid[nnn]
                savings = wealth + w_choice*labor - consumption

                if (use_idio_risk==1)
                    for kkk in 1 : n_idio
                        r_alpha_vec_idio[((kkk-1)*n_quad+1):(kkk*n_quad)] .=  curr_share.*rf_vec_use .+ (1.0 .- curr_share).*rk_vec.*exp.(idio_shock_grid[kkk,iii])
                        rf_vec_idio[((kkk-1)*n_quad+1):(kkk*n_quad)] = rf_vec_use
                    end
                   # Same question as before
     
                    # Other Operations
                    next_period_share_idio = (savings .* r_alpha_vec_idio .+ q_l_nxt_idio[:,iii]) ./ (tot_wealth_ss_vec[iii])
                    v_vec_twisted_idio = (v_idio[:,iii] .* next_period_share_idio) .^ (1.0 .- gma)
                    EV = sum(big_weight_vec_idio .* v_vec_twisted_idio)

                    deriv_helper_idio = bbeta .* (v_idio[:,iii].^((1.0/ies) - gma )) .* (next_period_share_idio.^(-gma)).* (EV.^( (gma - (1.0/ies))/(1.0-gma) )) .* mc_idio[:,iii]

                    if (constr_agents[iii] == 1 && constraint_binding[iii] == 1)
                        temp = sum(deriv_helper_idio .* rf_vec_idio .* big_weight_vec_idio)
                    else
                        temp = sum(deriv_helper_idio .* r_alpha_vec_idio .* big_weight_vec_idio)
                    end
                else # that is no idio risk
                    next_period_share = ( savings .* r_alpha_vec .+ q_l_nxt[:,iii]) ./ (tot_wealth_ss_vec[iii] .* exp.(dz_vec))
                    v_vec_twisted = (v_temp[:,iii] .* next_period_share) .^(1.0 .- gma)
                    EV = sum(big_weight_vec .* v_vec_twisted)


                    deriv_helper = bbeta .* (v_temp[:,iii].^((1.0/ies) - gma )) .* (next_period_share.^(-gma)).* (EV.^( (gma - (1.0/ies))/(1.0-gma) )) .* mc_temp[:,iii]

                    if (constr_agents[iii] == 1 && constraint_binding[iii] == 1)
                        temp = sum(deriv_helper .* rf_vec_use .* big_weight_vec)
                    else
                        temp = sum(deriv_helper .* r_alpha_vec .* big_weight_vec)
                    end
                end
                cons_update = labor_part .* ( temp.^(-ies))
                cons_update = min(max(cons_update, min_cons_sav), (wealth .+ w_choice .*labor .-min_cons_sav) )

                # update portfolio choice

                if (constr_agents[iii] == 1 && constraint_binding[iii] == 1)
                    share_update = 1.0 - kbar_constr[iii]*q_current/savings
                else
                    # find natural leverage limit to prevent negative payoffs
                    temp_vec = @. -rk_vec/(rf_vec_use-rk_vec + sqrt_eps )
                    temp_vec[temp_vec.>=0.0] .= -1E4

                    share_low = max(low_guess_fixed, (maximum(temp_vec)+sqrt_eps))

                    # find natural leverage limit to prevent negative payoffs
                    temp_vec = @. -rk_vec/(rf_vec_use-rk_vec + sqrt_eps )
                    temp_vec[temp_vec.<=0.0] .= 1E4
                    share_high = min(high_guess_fixed, (minimum(temp_vec)-sqrt_eps))

                    share_update = curr_share
                    share_update = calc_portfolio_share(1.0, share_low, share_high, iii, big_weight_vec, 
                    savings, labor, cons_update, 0.0, q_l_nxt[:,iii], 
                    v_temp[:,iii], mc_temp[:,iii], rf_vec_use, rk_vec, q_current, share_update, param::mpr_economy)
                end

                # Update values
                diff_inner = @. abs(consumption - cons_update)
                cons_grid[nnn] = @. consumption + 0.5*(cons_update - consumption)
                share_grid[nnn] = share_update
            end
        end
        # Compute MPC, MPR at the center point
        saving_grid .= @. wealth + w_choice*labor - cons_grid
        E_rf_grid .= @. E_rf * nom_i_grid / nom_i

        # MPC

        interp = interpolate(E_rf_grid, cons_grid, FritschCarlsonMonotonicInterpolation()) 
        dc_dErf_mat[iii] = ForwardDiff.derivative(extrapolate(interp, Flat()), E_rf)
        

        # MPB Actual (Bond)
        bond_grid .= @. share_grid * saving_grid
        interp = interpolate(E_rf_grid, bond_grid, FritschCarlsonMonotonicInterpolation()) 
        db_dErf_mat[iii] = ForwardDiff.derivative(extrapolate(interp, Flat()), E_rf)

       
        # MPK Actual (Capital)
        capital_grid .= @. (1-share_grid) * saving_grid / q_current

        interp = interpolate(E_rf_grid, capital_grid, FritschCarlsonMonotonicInterpolation()) 
        dk_dErf_mat[iii] = ForwardDiff.derivative(extrapolate(interp, Flat()), E_rf)

       
        if (constr_agents[iii] == 1 && constraint_binding[iii] == 1) 
            dk_dErf_mat[iii] = 0.0
        end
    end


    E_rf = sum(big_weight_vec .* log.(rf_vec) )   
    E_rk = sum(big_weight_vec .* log.(rk_vec) )   



    # store all variables

    # ============================================================ !
    # Store all relevant variables
    # ============================================================ !
    results_vec = zeros(70)

    results_vec[1] = y_current
    results_vec[2] = l_aggr
    results_vec[3] = c_vec[1]
    results_vec[4] = c_vec[2]
    results_vec[5] = c_vec[3]
    results_vec[6] = k_nxt[1] - (1.0-ddelta)*k_aggr
    results_vec[7] = pi_current
    results_vec[8] = E_rf
    results_vec[9] = E_rk
    results_vec[10] = q_current
    results_vec[11] = share_vec[1]
    results_vec[12] = share_vec[2]
    results_vec[13] = share_vec[3]
    results_vec[14] = sum(lmbd_vec .* savings_vec)
    results_vec[15] = aggr_wealth
    results_vec[16] = w_choice
    results_vec[17] = E_rk - E_rf
    results_vec[18] = infl
    results_vec[19] = nom_i
    results_vec[20] = sum(lmbd_vec .*savings_vec .* (1.0 .- share_vec))
    results_vec[21] = k_aggr
    results_vec[22] = v_store[1]
    results_vec[23] = v_store[2]
    results_vec[24] = v_store[3]
    results_vec[25] = q_l_current[1]
    results_vec[26] = q_l_current[2]
    results_vec[27] = q_l_current[3]
    results_vec[28] = MPC_mat[1]
    results_vec[29] = MPC_mat[2]
    results_vec[30] = MPC_mat[3]
    results_vec[31] = MPS_mat[1]
    results_vec[32] = MPS_mat[2]
    results_vec[33] = MPS_mat[3]
    results_vec[34] = MPK_mat[1]
    results_vec[35] = MPK_mat[2]
    results_vec[36] = MPK_mat[3]
    results_vec[37] = dc_dErk_mat[1] #
    results_vec[38] = dc_dErk_mat[2]
    results_vec[39] = dc_dErk_mat[3]
    results_vec[40] = db_dErk_mat[1] # 
    results_vec[41] = db_dErk_mat[2]
    results_vec[42] = db_dErk_mat[3]
    results_vec[43] = dk_dErk_mat[1] 
    results_vec[44] = dk_dErk_mat[2]
    results_vec[45] = dk_dErk_mat[3]
    results_vec[46] = dc_dErf_mat[1]
    results_vec[47] = dc_dErf_mat[2]
    results_vec[48] = dc_dErf_mat[3]
    results_vec[49] = db_dErf_mat[1]
    results_vec[50] = db_dErf_mat[2]
    results_vec[51] = db_dErf_mat[3]
    results_vec[52] = dk_dErf_mat[1]
    results_vec[53] = dk_dErf_mat[2]
    results_vec[54] = dk_dErf_mat[3]
    results_vec[55] = lmbd_vec[1]*savings_vec[1]*(1.0-share_vec[1])/q_current
    results_vec[56] = lmbd_vec[2]*savings_vec[2]*(1.0-share_vec[2])/q_current
    results_vec[57] = lmbd_vec[3]*savings_vec[3]*(1.0-share_vec[3])/q_current
    results_vec[58] = lmbd_vec[1]*savings_vec[1]
    results_vec[59] = lmbd_vec[2]*savings_vec[2]
    results_vec[60] = lmbd_vec[3]*savings_vec[3]
    results_vec[61] = lmbd_vec[1]*aggr_wealth*wealth_share_grid[1,sss]
    results_vec[62] = lmbd_vec[2]*aggr_wealth*wealth_share_grid[2,sss]
    results_vec[63] = lmbd_vec[3]*aggr_wealth*wealth_share_grid[3,sss]
    results_vec[64] = sum(big_weight_vec.*M_vec_mat[:,3].*rk_vec)
    results_vec[65] = constraint_binding[1] 
    results_vec[66] = constraint_binding[2]
    results_vec[67] = constraint_binding[3]
    results_vec[68] = mc_store[1] 
    results_vec[69] = mc_store[2]
    results_vec[70] = mc_store[3]

    return (who_prices_bond, bond_prices_out, results_vec)
end
