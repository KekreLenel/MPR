function calc_equilibrium_and_update(nxt_mat::Array{Float64, 2},  nxt_mat_2::Array{Float64, 1}, sss::Int64,
    c_vec::Array{Float64, 1}, l_aggr::Float64, q_current::Float64, q_l_current::Array{Float64, 1},
    infl::Float64, nom_i::Float64, share_vec::Array{Float64, 1}, monetary_shock::Float64, constraint_binding::Array{Int64,1}, param::mpr_economy;outer_iter::Int64=1)
    
    @unpack_mpr_economy param # unpack all the values
    constraint_binding_new = Array{Int64,1}(undef,n_I);
    constraint_binding_new .= 0;
    M_idio_mat = zeros(n_quad*n_idio, n_I)
    q_l_nxt_idio =  zeros(n_quad*n_idio, n_I)
    q_l_nxt = zeros(n_quad, n_I)
    v_idio =  zeros(n_quad*n_idio, n_I)
    mc_idio =  zeros(n_quad*n_idio, n_I)
    dz_vec_idio =  zeros(n_quad*n_idio)
    big_weight_vec_idio =  zeros(n_quad*n_idio)
    r_alpha_vec_idio =  zeros(n_quad*n_idio)
    next_period_share_idio =  zeros(n_quad*n_idio)
    v_vec_twisted_idio =  zeros(n_quad*n_idio)
    rf_vec_idio =  zeros(n_quad*n_idio)
    deriv_helper_idio =  zeros(n_quad*n_idio)
    M_idio =  zeros(n_quad*n_idio)
    c_vec_new = zeros(n_I)
    q_l_new = zeros(n_I)
    min_cons_sav = 1E-8 

    # declare some local variables
    v_temp = zeros(n_quad,n_I)
    mc_temp = zeros(n_quad, n_I)

    #         ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ !
    # Load Current State
    #         ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ !
    k_aggr = state_grid[idx_k, sss] # current capital stock, disaster prob & wage
    p_dis = exp(state_grid[idx_dis,sss])
    w_current = state_grid[idx_w, sss]

    # weighting over n_quad quadrature points, account for disaster
    big_weight_vec = zeros(n_quad) # initialize
    big_weight_vec[1:n_quad-1] = quad_weight_vec .* (1.0-p_dis)
    big_weight_vec[n_quad] = p_dis

    pi_current = aalpha * (l_aggr/k_aggr)^(1.0-aalpha)
    aggr_wealth = k_aggr * (pi_current + (1.0-ddelta)*q_current)
    y_current = k_aggr^aalpha * l_aggr^(1.0-aalpha)
    w_choice = (1.0-aalpha)*y_current/l_aggr 

    #         ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ !
    # Basic next-period updates
    #         ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ !

    # load next-period values: those all have dimensions n_quad
    k_nxt = copy(nxt_mat_2)
    for iii in 1:n_I
        v_temp[:,iii]  = nxt_mat[:,      iii] 
        mc_temp[:,iii] = nxt_mat[:,n_I + iii] 
    end
    q_nxt        = nxt_mat[:, 2*n_I+1]
    l_aggr_nxt   = nxt_mat[:, 2*n_I+2]
    infl_nxt     = nxt_mat[:, 2*n_I+3]
    q_l_nxt[:,1] = nxt_mat[:, 2*n_I+4] .* exp.(dz_vec)
    q_l_nxt[:,2] = nxt_mat[:, 2*n_I+5] .* exp.(dz_vec)
    q_l_nxt[:,3] = nxt_mat[:, 2*n_I+6] .* exp.(dz_vec)
    # Price of capital (exported)
    q_new = (k_nxt[1]/k_aggr)^(chiX)


    # Output
    y_next = @. exp((dz_vec)*(1.0-aalpha))*l_aggr_nxt*((k_nxt/l_aggr_nxt)^aalpha)

    # Wage
    w_next_choice = @. (1-aalpha)*y_next/l_aggr_nxt

    # Profit
    pi_nxt = @. aalpha*exp( (dz_vec)*(1.0-aalpha)) *(l_aggr_nxt^(1.0-aalpha))*k_nxt^(aalpha-1.0)


    # gross capital return (adjust in disaster state)
    rk_vec = @. ( (1.0 .- ddelta) .* q_nxt .+ pi_nxt ) ./ q_current
    rk_vec[n_quad] *= exp(dz_vec[n_quad])

    # Gross risk free return
    rf_vec = nom_i ./ infl_nxt

    #         ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ !
    # Solve for nominal rate
    #         ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ !
    # at this point inflation & consumption choices are held low_guess_fixed
    # excess demand solved for aggregate bond portfolios

    # current consump & rate
    c_temp, rf_vec_temp, share_temp = copy(c_vec), copy(rf_vec), copy(share_vec)

    # get excess bond

    # the following needs more work
   
    
    share_temp, k_temp_vec, b_temp_vec, excess_b, constraint_binding_new = calc_excess_bond_nom(aggr_wealth, v_temp, mc_temp, rf_vec_temp, rk_vec, q_current, w_choice, 
    l_aggr, k_nxt[1]*q_current, q_l_nxt, sss, big_weight_vec, c_temp, share_temp, constraint_binding, param::mpr_economy)

    # initial guess
    nom_i_guess = nom_i
    nom_i_A = nom_i_guess
    excess_b_A = excess_b

    # if no excess bond demand, just return the guess
    
    if (excess_b==0.0)
        nom_i_new = nom_i_guess
        excess_b_S = excess_b
        
    else
        if (excess_b > 0)
            iter = 0
            while (excess_b > 0.0 && iter<5000)
                iter = iter + 1
                
                
                nom_i_guess -= 5E-3

                # calc excess demand
                rf_vec_temp = nom_i_guess ./ infl_nxt
                share_temp = copy(share_vec)
                
                
                share_temp, k_temp_vec, b_temp_vec, excess_b , constraint_binding_new = calc_excess_bond_nom(aggr_wealth,v_temp,mc_temp,rf_vec_temp,rk_vec,
                q_current, w_choice, l_aggr, k_nxt[1]*q_current, q_l_nxt,sss,
                big_weight_vec, c_temp, share_temp, constraint_binding, param)


            end
            if (iter > 4999)
                error("no convergence in finding nom_i_low")
            end
            nom_i_B = nom_i_guess
            excess_b_B = excess_b
        else 
            iter = 0
            while (excess_b < 0.0 && iter<5000)
                iter = iter + 1
                
                nom_i_guess += 5E-3

                # calc excess demand
                rf_vec_temp = nom_i_guess ./ infl_nxt
                share_temp = copy(share_vec)

                share_temp, k_temp_vec, b_temp_vec, excess_b, constraint_binding_new  = calc_excess_bond_nom(aggr_wealth,v_temp,mc_temp,rf_vec_temp,rk_vec,
                q_current, w_choice, l_aggr, k_nxt[1]*q_current, q_l_nxt,sss,
                big_weight_vec, c_temp, share_temp, constraint_binding , param)
            end
            if (iter > 4999)
                error("no convergence in finding nom_i_low")
            end
            nom_i_B = nom_i_guess
            excess_b_B = excess_b
        end


        if (excess_b_A * excess_b_B > 0)
            error("initial bracket does not contian root")
        end
        # B should be closer to 0 than A. if not, swap
        if (abs(excess_b_A)<abs(excess_b_B))
            nom_i_A, nom_i_B = switch_two_var(nom_i_A,nom_i_B)
            excess_b_A, excess_b_B = switch_two_var(excess_b_A, excess_b_B)  
        end
    

        # (3) Do Brent Update
        nom_i_C = copy(nom_i_A); excess_b_C = copy(excess_b_A); excess_b_S = copy(excess_b_B)
        nom_i_D = Inf  # set a starting value for first iteration, will update it to be a finite number
        mflag = true; nom_i_S = 0;
        iter = 0
        while ( (excess_b_S != 0.0 ) && (abs(nom_i_A-nom_i_B)>1E-15) )
            iter = iter + 1
            #println("iter is $iter")
            if iter > 1000
                error("fails to find bond clearing rate in $iter iterations")
            end
            nom_i_S, mflag = brent_new_guess(nom_i_A, nom_i_B, nom_i_C, 
                                            excess_b_A, excess_b_B,excess_b_C, mflag; x_d = nom_i_D)
            rf_vec_temp = nom_i_S ./ infl_nxt
            share_temp = copy(share_vec)

            share_temp, k_temp_vec,b_temp_vec, excess_b_S , constraint_binding_new = calc_excess_bond_nom(aggr_wealth,v_temp,mc_temp,rf_vec_temp,rk_vec,
                q_current, w_choice, l_aggr, k_nxt[1]*q_current, q_l_nxt,sss,
                big_weight_vec, c_temp, share_temp,  constraint_binding, param)
            

            # Update positions
            nom_i_D = copy(nom_i_C)
            nom_i_C = copy(nom_i_B) ## 2022/09/06
            excess_b_C = copy(excess_b_B)

            if ((excess_b_A * excess_b_S)<0) 
                nom_i_B= copy(nom_i_S)
                excess_b_B = copy(excess_b_S)
            else
                nom_i_A = copy(nom_i_S)
                excess_b_A = copy(excess_b_S)
            end
            # swap a and B
            if (abs(excess_b_A) < abs(excess_b_B))
                nom_i_A, nom_i_B = switch_two_var(nom_i_A, nom_i_B)
                excess_b_A, excess_b_B = switch_two_var(excess_b_A, excess_b_B)
            end
        end
        nom_i_guess = nom_i_S
        nom_i_new   = nom_i_guess
    end

    # potential failure
    if (excess_b_S > 1E-2)
        error("excess aggregate bond holdings")
    end

    # store results of this block
    # nom_i is exported
    share_vec = copy(share_temp)
    nom_i_in = copy(nom_i)
    nom_i = copy(nom_i_guess)
    rf_vec = nom_i ./ infl_nxt

    # basic idiosyncratic risk expansion
    if (use_idio_risk == 1)
        for kkk in 1:n_idio
            for iii in 1:n_I
                v_idio[ ((kkk-1)*n_quad+1):(kkk*n_quad),iii].=v_temp[:,iii]
                mc_idio[ ((kkk-1)*n_quad+1):(kkk*n_quad),iii].=mc_temp[:,iii]
                q_l_nxt_idio[ ((kkk-1)*n_quad+1):(kkk*n_quad),iii] .= q_l_nxt[:,iii]*exp(idio_shock_grid[kkk,iii])
            end
            dz_vec_idio[((kkk-1)*n_quad+1):(kkk*n_quad)] = copy(dz_vec)
            rf_vec_idio[((kkk-1)*n_quad+1):(kkk*n_quad)] = copy(rf_vec)
            big_weight_vec_idio[((kkk-1)*n_quad+1):(kkk*n_quad)] .= big_weight_vec .* idio_weight_vec[kkk]
        end
    end

    #  =============================================================================
    # (4) Store value function, next-period capital, next-period wealth shares
    # ==============================================================================
        
    c_temp    = copy(c_vec)
    savings_vec = zeros(n_I); survivor_wealth_vec = zeros(n_quad,n_I); v_new = zeros(n_I); mc_new = zeros(n_I)
    s_nxt = zeros(n_quad, 3); labor_frac_nxt = zeros(n_quad, n_I)

    for iii in 1:n_I
        # Read parameters for agent
        ies   = ies_vec[iii]
        bbeta = bbeta_vec[iii]
        gma   = gma_vec[iii]
        thtbar  = thtbar_vec[iii]

        wealth = aggr_wealth .* wealth_share_grid[iii,sss]
        labor  = l_aggr .* labor_alloc_vec[iii]

        # Enforce minimal spending and mininal saving
        
        consumption =  minimum([maximum([c_temp[iii], min_cons_sav]),
                        wealth.+ w_choice.*labor.-min_cons_sav])
        savings = wealth .+ w_choice .* labor .- consumption
        savings_vec[iii] = savings
    
        # Portfolio return and wealth transitions
        r_alpha_vec = share_vec[iii] .* rf_vec .+ (1.0 .- share_vec[iii]) .* rk_vec
        survivor_wealth_vec[:,iii] .= lmbd_vec[iii] .* savings .* (r_alpha_vec)
        

        if (use_idio_risk == 1)
            # ------------------------------------------------------------ !
            # With Idiosyncratic Risk
            # ------------------------------------------------------------ !
            # Expand portfolio return
            for kkk in 1:n_idio
                r_alpha_vec_idio[((kkk-1)*n_quad+1):(kkk*n_quad)] = @. share_vec[iii]*rf_vec + (1.0 - share_vec[iii])*rk_vec*exp(idio_shock_grid[kkk,iii])
            end

            # Get Next Period Wealth Share
            next_period_share_idio = (savings .* r_alpha_vec_idio .+ q_l_nxt_idio[:,iii] ) ./ (tot_wealth_ss_vec[iii] )
            v_vec_twisted_idio = (v_idio[:,iii] .* next_period_share_idio) .^ (1.0 .- gma)

            # Calculate Expected Value Function
            EV = sum(big_weight_vec_idio .* v_vec_twisted_idio)
        else
            # ------------------------------------------------------------ !
            # No Idiosyncratic Risk
            # ------------------------------------------------------------ !
            # Get Next Period Wealth Share
            next_period_share = (savings .* r_alpha_vec .+ q_l_nxt[:,iii] ) ./ (tot_wealth_ss_vec[iii] )
            v_vec_twisted = (v_temp[:,iii] .* next_period_share) .^ (1.0 .- gma)

            # Calculate Expected Value Function
            EV = sum(big_weight_vec .* v_vec_twisted)
        end
        

        # Value Function Update
        util, util_c_deriv, labor_part= util_fun(consumption, labor, iii, thtbar_vec, ies_vec, tht)
        objf = ( (1.0-bbeta) .* util .+ bbeta .*  (EV .^( (1.0 - (1.0 ./ ies)) ./ (1.0 .- gma)))) .^ (1.0 ./ (1.0-1.0 ./ ies))
        v_new[iii] = objf.* tot_wealth_ss_vec[iii]/(wealth .+ q_l_current[iii]) 
        mc_new[iii] = util_c_deriv .* (tot_wealth_ss_vec[iii]/(wealth + q_l_current[iii])).^(-1.0/ies)
    end

    # (4.2) Next period wealth share
    temp_nxt_aggr_wealth = vec(sum(survivor_wealth_vec, dims = 2)) # original code wrong
    # sum by the second dimension, then convert to vector (original output is a matrix)


    for iii in 1:n_I
        # Next period theta, accounting for survivors and newborns
        # EXPORTED
        s_nxt[:,iii] .= ((1.0-xi) .* survivor_wealth_vec[:,iii] .+ xi .* s_bar_vec[iii] .* temp_nxt_aggr_wealth) ./ (temp_nxt_aggr_wealth)
        labor_frac_nxt[:,iii] .= ((1.0-xi) .* (survivor_wealth_vec[:,iii] .+ lmbd_vec[iii] .* q_l_nxt[:,iii])) ./
            ((1.0-xi) .* survivor_wealth_vec[:,iii] .+ xi .* s_bar_vec[iii] .* temp_nxt_aggr_wealth .+ lmbd_vec[iii] .* q_l_nxt[:,iii])
    end

    # (4.3) Next period capital
    # EXPORTED
    k_next_new = sum(lmbd_vec .* savings_vec .* (1.0 .- share_vec) ./q_current)

    # =============================================================================
    # (6) Update Consumption Choice
    # =============================================================================
    c_temp = c_vec
    for iii in 1:n_I
        savings = 0.0
        # (6.1) Read parameters and quantities (same as in part (4) above)
        ies   = ies_vec[iii]
        bbeta = bbeta_vec[iii]
        gma   = gma_vec[iii]
        thtbar  = thtbar_vec[iii]

        wealth = aggr_wealth .* wealth_share_grid[iii,sss]
        labor  = l_aggr .* labor_alloc_vec[iii]
        
        consumption =  minimum([maximum([c_temp[iii], min_cons_sav]),
                          wealth .+ w_choice .* labor .- min_cons_sav])
        r_alpha_vec = share_vec[iii] .* rf_vec .+ (1.0 .- share_vec[iii]) .* rk_vec
        labor_part =  1.0 .+ (1.0 ./ ies .- 1.0) .* thtbar .* tht ./ (1.0 .+ tht) .* labor .^((1.0 .+ tht) ./ tht)
        # (6.2) Loop and update consumption
        diff_inner = 1.0
        iter_inner = 0
        while (diff_inner > 1e-08)
            iter_inner = iter_inner + 1

            # Based on equation (48/54) of equation doc
            savings = wealth .+ w_choice .* labor .- consumption


            if (use_idio_risk == 1)
                # ------------------------------------------------------------ !
                # With Idiosyncratic Risk
                # ------------------------------------------------------------ !
                # Expand portfolio return
                for kkk in 1:n_idio
                    r_alpha_vec_idio[((kkk-1)*n_quad+1):(kkk*n_quad)] = @. share_vec[iii]*rf_vec + (1.0 - share_vec[iii])*rk_vec*exp(idio_shock_grid[kkk,iii])
                end

            #Same question here as before

                # Other Operations
                next_period_share_idio = (savings .* r_alpha_vec_idio .+ q_l_nxt_idio[:,iii]) ./ (tot_wealth_ss_vec[iii])
                v_vec_twisted_idio = (v_idio[:,iii] .* next_period_share_idio) .^ (1.0 .- gma)
                EV = sum(big_weight_vec_idio .* v_vec_twisted_idio)

                deriv_helper_idio = bbeta.*(v_idio[:,iii].^((1.0/ies) - gma )) .* (next_period_share_idio.^(-gma)) .*(EV .^( (gma - (1.0/ies))/(1.0-gma) )) .* mc_idio[:,iii]

                if (constr_agents[iii] == 1 && constraint_binding[iii] == 1)
                    temp = sum(deriv_helper_idio .* rf_vec_idio .* big_weight_vec_idio)
                else
                    temp = sum(deriv_helper_idio .* r_alpha_vec_idio .* big_weight_vec_idio)
                end

            else
                # ------------------------------------------------------------ !
                # No Idiosyncratic Risk
                # ------------------------------------------------------------ !
                next_period_share = ( savings .* r_alpha_vec .+ q_l_nxt[:,iii]) ./ (tot_wealth_ss_vec[iii])
                v_vec_twisted = (v_temp[:,iii] .* next_period_share) .^(1.0 .- gma)
                EV = sum(big_weight_vec .* v_vec_twisted)

                deriv_helper = bbeta.*(v_temp[:,iii].^((1.0/ies) - gma )) .* (next_period_share.^(-gma)) .*(EV .^( (gma - (1.0/ies))/(1.0-gma) )) .* mc_temp[:,iii]

                if (constr_agents[iii] == 1 && constraint_binding[iii] == 1)
                    temp = sum(deriv_helper .* rf_vec .* big_weight_vec)
                else
                    temp = sum(deriv_helper .* r_alpha_vec .* big_weight_vec)
                end
            end

            cons_update = labor_part .* temp .^ (-ies)
            # Update
            cons_update = min(max(cons_update, min_cons_sav), (wealth .+ w_choice .*labor .-min_cons_sav) )
            diff_inner = abs(consumption .- cons_update)
            consumption = consumption .+ 0.5 .* (cons_update .- consumption)
            if (iter_inner > 1000)
                if (diff_inner > 0.1)
                    error("NO CONVERGENCE 1.")
                else
                    error("'WARNING: NO CONVERGENCE 1.")
                    diff_inner = 0.0
                end
            end
        end

        # EXPORT
        c_vec_new[iii] = consumption
        savings_vec[iii] = savings

    end


    # =============================================================================
    # (7) Update marginal utility and pricing kernel
    # =============================================================================
    M_vec_mat =zeros(n_quad,n_I)
    marg_util_c_vec = zeros(n_I)

    for iii in 1:n_I
        # (7.1) Read parameters and quantities (same as in part (4) above)
        ies   = ies_vec[iii]
        bbeta = bbeta_vec[iii]
        gma   = gma_vec[iii]
        thtbar  = thtbar_vec[iii]

        wealth = aggr_wealth .* wealth_share_grid[iii,sss]
        labor  = l_aggr .* labor_alloc_vec[iii]

        consumption =  min(max(c_vec_new[iii], min_cons_sav),
                            (wealth + w_choice *labor -min_cons_sav) )
        r_alpha_vec = share_vec[iii] .* rf_vec .+ (1.0 .- share_vec[iii]) .* rk_vec
        savings     = wealth .+ w_choice .* labor .- consumption

        # (7.2) Calculate (Realized) Pricing Kernel
        if (use_idio_risk == 1)
            # ------------------------------------------------------------ !
            # With Idiosyncratic Risk
            # ------------------------------------------------------------ !
            # Expand portfolio return
            for kkk in 1 : n_idio
                r_alpha_vec_idio[((kkk-1)*n_quad+1):(kkk*n_quad)] = @. share_vec[iii]*rf_vec + (1.0 - share_vec[iii])*rk_vec*exp(idio_shock_grid[kkk,iii])
            end
            # Same question as before

            # Other Operations
            next_period_share_idio = (savings .* r_alpha_vec_idio .+ q_l_nxt_idio[:,iii]) ./ (tot_wealth_ss_vec[iii])
            v_vec_twisted_idio = (v_idio[:,iii] .* next_period_share_idio) .^ (1.0 .- gma)
            EV = sum(big_weight_vec_idio .* v_vec_twisted_idio)

            util, util_c_deriv, labor_part = util_fun(consumption, labor, iii, thtbar_vec, ies_vec, tht  )

            objf = ((1.0 - bbeta) .* util .+ bbeta .* (EV .^( (1.0 .- (1.0 ./ ies)) ./ (1.0 .- gma)))) .^ (1.0 ./ (1.0-1.0 ./ ies))

            marg_util_c_vec[iii] = objf .^ (1.0 ./ ies) .*(1.0 - bbeta).* util_c_deriv

            # SDF: first, the full form
            M_idio = bbeta* ( next_period_share_idio .^ (-gma) ) .* (EV .^( (gma - (1.0/ies))/(1.0-gma) )) .*(v_idio[:,iii].^((1.0/ies) -gma )) .* mc_idio[:,iii] ./ util_c_deriv 

            M_idio_mat[:,iii] = M_idio

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

            util, util_c_deriv, labor_part = util_fun(consumption, labor, iii, thtbar_vec, ies_vec, tht  )

            objf = ( (1.0-bbeta) .* util .+  bbeta .*  (EV .^( (1.0 .- (1.0 ./ ies))/(1.0 .- gma)))) .^( 1.0 ./ (1.0 .- 1.0 ./ies))

            marg_util_c_vec[iii] = objf .^ (1.0 ./ ies) .*(1.0-bbeta) .* util_c_deriv

            M_vec_mat[:,iii] .=   bbeta* ( next_period_share .^ (-gma) ) .* (EV .^( (gma - (1.0/ies))/(1.0-gma) )) .* (v_temp[:,iii].^((1.0/ies) -gma )) .* mc_temp[:,iii] ./ util_c_deriv

        end


    end

       if (any(isnan.(M_vec_mat)))
               error("M_vec_mat is NaN.")
       end
    
       # =============================================================================
       # (8) Update Pricing of Time Endowment
       # =============================================================================
       for  iii in 1:n_I

           #  EXPORTED, CHANGED
           # Is changed a command in Julia?
           if (labor_alloc_vec[iii] > sqrt(eps()))
                if (use_idio_risk == 1)
                    q_l_new[iii] = w_choice*l_aggr*labor_alloc_vec[iii] + sum(big_weight_vec_idio.*M_idio_mat[:,iii].*q_l_nxt_idio[:,iii]) 
                else
                    q_l_new[iii] =  l_aggr*labor_alloc_vec[iii]*w_choice +  sum(big_weight_vec.*M_vec_mat[:,iii].*q_l_nxt[:,iii])              
                end
           else
               q_l_new .= 0.0
           end
       end

    # =============================================================================
    # (9) Update labor, wages
    # =============================================================================
    # Combined Weights
    weight_vec = lmbd_vec .* marg_util_c_vec .* labor_alloc_vec

    if (chiW >= -1.0)
        # Wage rigidity
        l_aggr_temp   = copy(l_aggr)
        w_choice_temp = copy(w_choice)

        # Find wage implied by Union optimization
        w_temp = copy(w_choice_temp)
        w_diff = 1.0
        iter = 0
    
        while ( (w_diff > 1e-08 && iter<100) || iter < 10 )
            iter = iter + 1

            w_current_next = w_temp .* ones(n_quad)
            w_current_next[n_quad] = w_temp .* exp.(dz_vec[n_quad])


            w_temp_new = 0.0

            for iii = 1:n_I
                # Need to make sure this equation is correct & understand how the union works
                if (labor_alloc_vec[iii] > sqrt(eps()))

                w_temp_new = w_temp_new .+  weight_vec[iii] ./ sum(weight_vec)  .* (vareps_w ./ ((1.0 .-vareps_w) .*  (1.0 .-tau_w))  .* ( -  1.0 ./ ies_vec[iii] .* c_vec[iii] .* thtbar_vec[iii] .* (l_aggr_temp .^ (1.0 ./ tht)) ./
                                (1.0 .+ (1.0 ./ ies_vec[iii] - 1.0) .* thtbar_vec[iii] .* tht ./(1.0 .+ tht) .*
                                l_aggr_temp .^((1.0 .+tht) ./tht))) + w_temp .*( 1.0 ./labor_alloc_vec[iii]) .*chiW ./ ( (1.0-vareps_w) .*  (1.0 .- tau_w) ) .* ( (
                                w_temp ./w_current .*infl .- 1.0  ) .*w_temp ./w_current .*infl  - sum( M_vec_mat[:,iii] .* big_weight_vec .* labor_frac_nxt[:,iii] .*
                                ( w_next_choice ./w_current_next .*infl_nxt .- 1.0 ) .*
                                (( w_next_choice ./w_current_next) .*( w_next_choice ./w_temp) ) .* infl_nxt .* l_aggr_nxt ./l_aggr_temp)))
                end
            end

            w_diff = abs(w_temp_new-w_temp)
            w_temp = w_temp + 0.005*(w_temp_new-w_temp)

            if (isnan.(w_temp))
                error("error w isnan, $weight_vec, $M_vec_mat")
            elseif (w_temp<=0)
                error("error w<0")
            end
            l_aggr_temp = k_aggr .* (((1.0 .- aalpha) ./ w_temp) .^ (1.0 ./ aalpha))
            w_choice_temp = copy(w_temp)
            
        end

        # EXPORTED
        l_aggr_new = copy(l_aggr_temp)

    else
        # No wage ridigity
        l_aggr_new = sum(lmbd_vec .* l_temp_vec)
    end

    # EXPORTED
    w_choice_new = w_choice

    # =============================================================================
    # (10) Update Inflation Rate
    # =============================================================================
    infl_new = ( exp(-tayl_ic) * exp(-state_grid[idx_m,sss] - monetary_shock)*nom_i_in  )^(1.0/phi)


    constraint_binding = copy(constraint_binding_new) 

    return share_vec, nom_i,q_new, q_l_new, k_next_new*ones(n_quad), c_vec_new, infl_new, v_new, mc_new, l_aggr_new, w_choice_new, s_nxt, constraint_binding_new

end
