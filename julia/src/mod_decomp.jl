#This file translates the code in mod_decomp.f90



function monetary_decomposition(param::mpr_economy,paramid::Int64)
    @unpack_mpr_economy param # unpack all the values


	# ============================================================ !
	# Step 1: Basic Loading and Processing
	# ============================================================ !

	n_nxt = 7 + 2*n_I;  ## 07/15/2022


	next_state_mat = zeros(n_quad, n_active_dims, n_states)
	next_state_monetary_mat = zeros(n_quad, n_active_dims, n_states)
	next_state_montransition_mat = zeros(1, n_active_dims, n_states)

	interp_input_mat = zeros(n_states, n_interp)
	interp_input_monetary_mat = zeros(n_states, n_interp)
	state_variables = zeros(smolyak_d)
	other_variables = zeros(n_interp)
	state_variables_star = zeros(smolyak_d)
	other_variables_star = zeros(n_interp)

	current_state = zeros(1,n_active_dims)
	next_state = zeros(1,n_active_dims)
	stochsteady_state = zeros(1,n_active_dims)

	state_series_0_mat=zeros(n_sample, smolyak_d)
	other_vars_series_0_mat=zeros(n_sample, n_interp)
	state_series_star_mat=zeros(n_sample, smolyak_d)
	other_vars_series_star_mat=zeros(n_sample, n_interp)
	state_series_1_mat=zeros(n_sample, smolyak_d)
	other_vars_series_1_mat=zeros(n_sample, n_interp)
	starting_states = zeros(n_sample,n_active_dims)
	S_0 = zeros(1,n_active_dims)
	S_Star = zeros(1,n_active_dims)
	S_Trans = zeros(1,n_active_dims)
	S_1 = zeros(1,n_active_dims)
    smol_coeffs_use=zeros(n_states,n_nxt)
	nxt_mat = zeros(n_quad, n_nxt)
	nxt_mat_2 = zeros(n_quad)
	interp_use = zeros(n_states, n_nxt)
	b4 = zeros(n_states, n_nxt)

	c_vec = zeros(n_I)
	share_vec = zeros(n_I)
	results_a = zeros(4)
	results_b = zeros(4)
	results_c = zeros(4)
	cf_policy_matA = zeros(n_sample,n_I, 4)
	cf_policy_matB = zeros(n_sample,n_I, 4)

	# ------------------------------------------------------------ !
	# Read Solution Files
	# ------------------------------------------------------------ !

	println("")
	println("READING RESULTS")

	next_state_mat = load(results_folder*"results.jld","next_state_mat")
	interp_input_mat = load(results_folder*"results.jld","results_mat")

	println(interp_input_mat[1,3:5])
	# error("asd")

	# ------------------------------------------------------------ !
	# Find Stochastic Steady State
	# ------------------------------------------------------------ !
	println("")
	println("FIND STOCHASTIC STEADY STATE")

	# (1.1) Polynomial coefficients for transition matrix
	nrhs = n_active_dims
	b=zeros(n_states,nrhs)
	smol_montransition_coeffs=zeros(1,n_states,nrhs)
	smol_monetary_coeffs=zeros(n_quad,n_states,nrhs)
	smol_coeffs=zeros(n_quad,n_states,nrhs)

	for qqq = 1 : n_quad
	    b = transpose(next_state_mat[qqq,:,:])
		smol_coeffs[qqq,:,:] =  smol_polynom\b
	end

	# (1.2) Polynomial coefficients for non-state variables
	nrhs2 = n_interp
	interp_coeffs = smol_polynom\interp_input_mat
	interp_monetary_coeffs=zeros(n_states,nrhs2)

	# (1.3) Polynomial coefficients for state variables
	b3 = transpose(state_grid)
	state_coeffs = smol_polynom\b3

	# (2.1) Start from some interior grid point and move forward
	# Use no_shock_idx: the quadrature point where there's no shock happening
	# current_state(1,:) = smol_grid(1,:)
	current_state[1,:] .= 0
	diff = 1.0

	while (diff > sqrt(eps()))
		polyn_points = Smolyak_Polynomial(current_state, n_active_dims,
							max_smol_level, smol_elem_ani)

		next_state = polyn_points*smol_coeffs[no_shock_idx,:,:]

		next_state[next_state.>1.0] .= 1.0
		next_state[next_state .< -1.0] .= -1.0

		diff = maximum(abs.(current_state .- next_state))
		current_state = next_state
	end
	stochsteady_state = current_state


	starting_states = load(results_folder*"starting_states.jld","starting_states")

	next_state_monetary_mat,next_state_montransition_mat, interp_input_monetary_mat = load(results_folder*"/result_irf_2.jld","next_state_monetary_mat", "next_state_montransition_mat", "results_monetary_mat")
	# ============================================================ !
	# Step 1: Find S* and n^i(S*)
	# For now, just start from the SSS, so S* is S0
	# But write the code for any transition
	# ============================================================ !
	# Start from S0 (for now stochastic steady state)
    cf_policy_matA .= 0.0
    cf_policy_matB .= 0.0

    for iii = 1:n_sample
            
        if (mod(iii, 50) == 0)
            println(iii)
        end

		S_0[1,:] = starting_states[iii,:]
		polyn_points = Smolyak_Polynomial(S_0, n_active_dims, max_smol_level, smol_elem_ani)

		# Policies (calc_bond_price format) and states at state S_0

		other_variables = polyn_points*interp_coeffs
		
		if other_variables[67] > 0.5 
			other_variables[36] = 0.0
			other_variables[45] = 0.0
			other_variables[54] = 0.0
			other_variables[57] = lmbd_vec[3]*kbar_constr[3]
		end		 

		state_variables = polyn_points*state_coeffs

		state_series_0_mat[iii,:]      = state_variables
		other_vars_series_0_mat[iii,:] = other_variables

		S_Star = polyn_points*smol_coeffs[no_shock_idx,:,:]

		polyn_points = Smolyak_Polynomial(S_Star, n_active_dims, max_smol_level, smol_elem_ani) 

		other_variables_star = polyn_points*interp_coeffs



		if other_variables_star[67] > 0.5 
			other_variables_star[36] = 0.0
			other_variables_star[45] = 0.0
			other_variables_star[54] = 0.0
			other_variables_star[57] = lmbd_vec[3]*kbar_constr[3]
		end		
		state_variables_star = polyn_points*state_coeffs

		state_series_star_mat[iii,:]      = state_variables_star 
		other_vars_series_star_mat[iii,:] = other_variables_star

		# Extract per-person wealth by type
		wealth_a = other_variables_star[15] * state_variables_star[2] / lmbd_vec[1]
		wealth_b = other_variables_star[15] * state_variables_star[3] / lmbd_vec[2]
		wealth_c = other_variables_star[15] * (1.0 - state_variables_star[2] - state_variables_star[3]) / lmbd_vec[3]


		# ============================================================ !
		# Step 2: Extract relevant information about the shocked state S1
		# ============================================================ !
		# Read Monetary Transition Files

		for qqq = 1:n_quad
			b = transpose(next_state_monetary_mat[qqq,:,:])
			smol_monetary_coeffs[qqq,:,:] = smol_polynom\b
		end

		smol_montransition_coeffs[1,:,:] = smol_polynom\transpose(next_state_montransition_mat[1,:,:])

		interp_monetary_coeffs = smol_polynom\interp_input_monetary_mat

		polyn_points = Smolyak_Polynomial(S_0, n_active_dims, max_smol_level, smol_elem_ani) 

		S_1 = polyn_points*smol_montransition_coeffs[1,:,:]

		polyn_points = Smolyak_Polynomial(S_1, n_active_dims, max_smol_level, smol_elem_ani) 

		other_variables = polyn_points*interp_monetary_coeffs

		if other_variables[67] > 0.5 
			other_variables[36] = 0.0
			other_variables[45] = 0.0
			other_variables[54] = 0.0
			other_variables[57] = lmbd_vec[3]*kbar_constr[3]
		end		

		state_variables = polyn_points*state_coeffs

		state_series_1_mat[iii,:]      = state_variables
		other_vars_series_1_mat[iii,:] = other_variables

		# Extract Relevant Information
		p_dis = exp(state_variables[6])
		k_aggr = state_variables[1]
		l_aggr = other_variables[2]
		q_current = other_variables[10]
		infl = other_variables[18]
		nom_i = other_variables[19]
		c_vec[1] = other_variables[3]
		c_vec[2] = other_variables[4]
		c_vec[3] = other_variables[5]

		share_vec[1] = other_variables[11]
		share_vec[2] = other_variables[12]
		share_vec[3] = other_variables[13]

		# nxt_mat, n_nxt, nxt_mat_2
		# Need the relevant interp_input_mat
		# Now we are standing in state S1 and there is no more shock, so use interp_input_mat directly
		interp_use[:,1] .= interp_input_mat[:,22]        # Value a
		interp_use[:,2] .= interp_input_mat[:,23]        # Value b
		interp_use[:,3] .= interp_input_mat[:,24]        # Value c
		interp_use[:,4] .= interp_input_mat[:,68]        # Value a
		interp_use[:,5] .= interp_input_mat[:,69]        # Value b
		interp_use[:,6] .= interp_input_mat[:,70]        # Value c
		interp_use[:,2*n_I + 1]  .= interp_input_mat[:,10] # q
		interp_use[:,2*n_I + 2]  .= interp_input_mat[:,2] # l
		interp_use[:,2*n_I + 3]  .= interp_input_mat[:,18] # Inflation
		interp_use[:,2*n_I + 4]  .= interp_input_mat[:,25] # q_l a
		interp_use[:,2*n_I + 5]  .= interp_input_mat[:,26] # q_l b
		interp_use[:,2*n_I + 6]  .= interp_input_mat[:,27] # q_l c

		# Coefficients
		b4 = copy(interp_use)
		smol_coeffs_use = smol_polynom \ b4
		# Create nxt_mat
		nxt_mat .= 0.0
		nxt_mat_2 .= 0.0

		for qqq = 1 : n_quad
			polyn_points = Smolyak_Polynomial(S_1, n_active_dims, max_smol_level, smol_elem_ani)

			# S_Trans is next state when n_quad qqq happens
			S_Trans = polyn_points*smol_monetary_coeffs[qqq,:,:]


			polyn_points = Smolyak_Polynomial(S_Trans, n_active_dims, max_smol_level, smol_elem_ani)
			
			nxt_mat[qqq,:] = polyn_points*smol_coeffs_use
		end

		k_next = other_variables[20] ./ other_variables[10]
		nxt_mat_2 = k_next ./ exp.(dz_vec_adj) .* exp.(dz_vec)

		# ============================================================ !
		# Step 3A: Compute policies functions like k^i(n^i(S^*), S_1)
		# Everything in step 1, except wealth
		# ============================================================ !
		# Current Variables
		results_a=calc_counterfactual_policies(1, wealth_a, p_dis, k_aggr, l_aggr, q_current, infl, nom_i, round(other_variables[65]),nxt_mat, n_nxt, nxt_mat_2, c_vec, share_vec,param,iii)

		results_b=calc_counterfactual_policies(2, wealth_b, p_dis, k_aggr, l_aggr, q_current, infl, nom_i, round(other_variables[66]),nxt_mat, n_nxt, nxt_mat_2, c_vec, share_vec,param,iii)

		results_c=calc_counterfactual_policies(3, wealth_c, p_dis, k_aggr, l_aggr, q_current, infl, nom_i, round(other_variables[67]),nxt_mat, n_nxt, nxt_mat_2, c_vec, share_vec,param,iii)

		# Counterfactuals
		# cf_policy_mat .= 0.0
		cf_policy_matA[iii,1,:] .= results_a
		cf_policy_matA[iii,2,:] .= results_b
		cf_policy_matA[iii,3,:] .= results_c

		
		# ============================================================ !
		# Step 3B: Hold Expected Returns Constant
		# ============================================================ !

		q_use = q_current .* other_variables[9] ./ other_variables_star[9]
		i_use = nom_i .* exp(other_variables_star[8]) ./ exp(other_variables[8])


		results_a=calc_counterfactual_policies(1, wealth_a, p_dis, k_aggr, l_aggr, q_use, infl, i_use, round(other_variables[65]), nxt_mat, n_nxt, nxt_mat_2, c_vec, share_vec,param,iii)
		results_b=calc_counterfactual_policies(2, wealth_b, p_dis, k_aggr, l_aggr, q_use, infl, i_use, round(other_variables[66]), nxt_mat, n_nxt, nxt_mat_2, c_vec, share_vec,param,iii)
		results_c=calc_counterfactual_policies(3, wealth_c, p_dis, k_aggr, l_aggr, q_use, infl, i_use, round(other_variables[67]), nxt_mat, n_nxt, nxt_mat_2, c_vec, share_vec,param,iii)

		cf_policy_matB[iii,1,:] .= results_a
		cf_policy_matB[iii,2,:] .= results_b
		cf_policy_matB[iii,3,:] .= results_c
	end

	
    # ============================================================ !
    # Step 4: Collect results for further processing
    # ============================================================ !

	save(results_folder*"Counterfactual_Policies_A.jld","cf_policy_mat",sum(cf_policy_matA,dims = 1)[1,:,:]./n_sample)
	save(results_folder*"Counterfactual_Policies_B.jld","cf_policy_mat",sum(cf_policy_matB,dims = 1)[1,:,:]./n_sample)

	save(results_folder*"Policies_S_0.jld","other_variables",sum(other_vars_series_0_mat, dims = 1)[1,:,:]./n_sample)
	save(results_folder*"States_S_0.jld","state_variables",sum(state_series_0_mat, dims = 1)[1,:,:]./n_sample)

	save(results_folder*"Policies_S_Star.jld","other_variables",sum(other_vars_series_star_mat, dims = 1)[1,:,:]./n_sample)
	save(results_folder*"States_S_Star.jld","state_variables",sum(state_series_star_mat, dims = 1)[1,:,:]./n_sample)


	save(results_folder*"Policies_S_1.jld","other_variables",sum(other_vars_series_1_mat, dims = 1)[1,:,:]./n_sample)
	save(results_folder*"States_S_1.jld","state_variables",sum(state_series_1_mat, dims = 1)[1,:,:]./n_sample)

end



function calc_counterfactual_policies(iii, wealth, p_dis, k_aggr, l_aggr, q_current, infl, nom_i, constraint_binding,
										nxt_mat, n_nxt, nxt_mat_2, c_vec, share_vec, param::mpr_economy,sample_id)

	@unpack_mpr_economy param
	
    results_vec  = zeros(4) # Consumption, Share, Bond, Capital
    consumption = 0.0
	share_update = 0.0
	savings = 0.0
	q_nxt = zeros(n_quad)
	rk_vec = zeros(n_quad)
    v_temp = zeros(n_quad,n_I)
	mc_temp = zeros(n_quad,n_I)
	w_next_choice = zeros(n_quad)
	pi_nxt = zeros(n_quad)
	y_next = zeros(n_quad)
	infl_nxt = zeros(n_quad)
	k_nxt = zeros(n_quad)
	rf_vec = zeros(n_quad)
	min_cons_sav = 1E-08
 
    r_alpha_vec = zeros(n_quad)
    q_l_nxt = zeros(n_quad, n_I)
	v_vec_twisted = zeros(n_quad)
    l_aggr_nxt = zeros(n_quad)
    big_weight_vec = zeros(n_quad)
	deriv_helper = zeros(n_quad)


    temp_vec = zeros(n_quad)

    v_idio = zeros(n_quad*n_idio, n_I)
	mc_idio = zeros(n_quad*n_idio, n_I)
	q_l_nxt_idio = zeros(n_quad*n_idio, n_I)

    dz_vec_idio = zeros(n_quad*n_idio)
	big_weight_vec_idio = zeros(n_quad*n_idio)
    r_alpha_vec_idio = zeros(n_quad*n_idio)
	next_period_share_idio = zeros(n_quad*n_idio)
    v_vec_twisted_idio = zeros(n_quad*n_idio)
	rf_vec_idio = zeros(n_quad*n_idio)
    deriv_helper_idio = zeros(n_quad*n_idio)


    # ============================================================ !
    # Preparations: Current Period
    # ============================================================ !
    big_weight_vec[1:n_quad-1] .= quad_weight_vec .* (1.0-p_dis)
    big_weight_vec[n_quad]     = p_dis
    pi_current = aalpha .* (l_aggr ./ k_aggr) .^ (1.0-aalpha)
    aggr_wealth = k_aggr .* (pi_current .+ (1.0-ddelta)*q_current)
    y_current = k_aggr .^ aalpha .* l_aggr .^ (1.0-aalpha)
    w_choice = (1.0-aalpha) .* y_current ./ l_aggr

    # ============================================================ !
    # Preparations: Net Period
    # ============================================================ !
    k_nxt = copy(nxt_mat_2)
    v_temp[:,1] .= nxt_mat[:,1] 
    v_temp[:,2] .= nxt_mat[:,2] 
    v_temp[:,3] .= nxt_mat[:,3]
	mc_temp[:,1] .= nxt_mat[:,4] 
    mc_temp[:,2] .= nxt_mat[:,5] 
    mc_temp[:,3] .= nxt_mat[:,6]
    q_nxt        = nxt_mat[:, 2*n_I+1]
    l_aggr_nxt   = nxt_mat[:, 2*n_I+2]
    infl_nxt     = nxt_mat[:, 2*n_I+3]
	# println(infl_nxt)
    q_l_nxt[:,1] .= nxt_mat[:, 2*n_I+4] .* exp.(dz_vec)
    q_l_nxt[:,2] .= nxt_mat[:, 2*n_I+5] .* exp.(dz_vec)
    q_l_nxt[:,3] .= nxt_mat[:, 2*n_I+6] .* exp.(dz_vec)

    # Output
    y_next = @. exp((dz_vec)*(1.0-aalpha))*l_aggr_nxt* ((k_nxt/l_aggr_nxt)^aalpha)

    # Wage
    w_next_choice = @. (1-aalpha) .*y_next ./ l_aggr_nxt

    pi_nxt = @. aalpha* exp( (dz_vec)* (1.0- aalpha))* (l_aggr_nxt^ (1.0- aalpha))*k_nxt^ (aalpha-1.0)

    # Gross Capital Return (Adjust in disaster state)
    rk_vec = @. ( (1.0 .- ddelta) .* q_nxt .+ pi_nxt ) ./ q_current
    rk_vec[n_quad] = @. exp(dz_vec[n_quad]) * rk_vec[n_quad]

    # Gross Risk-Free Return
    rf_vec = nom_i ./ infl_nxt


    # Basic Idiosyncratic Risk Expansion
    if (use_idio_risk == 1)
        for kkk = 1:n_idio
            v_idio[((kkk-1)*n_quad+1):(kkk*n_quad),1] .= v_temp[:,1]
            v_idio[((kkk-1)*n_quad+1):(kkk*n_quad),2] .= v_temp[:,2]
            v_idio[((kkk-1)*n_quad+1):(kkk*n_quad),3] .= v_temp[:,3]
			mc_idio[((kkk-1)*n_quad+1):(kkk*n_quad),1] .= mc_temp[:,1]
            mc_idio[((kkk-1)*n_quad+1):(kkk*n_quad),2] .= mc_temp[:,2]
            mc_idio[((kkk-1)*n_quad+1):(kkk*n_quad),3] .= mc_temp[:,3]
            q_l_nxt_idio[((kkk-1)*n_quad+1):(kkk*n_quad),1] = q_l_nxt[:,1]
            q_l_nxt_idio[((kkk-1)*n_quad+1):(kkk*n_quad),2] = q_l_nxt[:,2]
            q_l_nxt_idio[((kkk-1)*n_quad+1):(kkk*n_quad),3] = q_l_nxt[:,3]
            dz_vec_idio[((kkk-1)*n_quad+1):(kkk*n_quad)] .= dz_vec
            rf_vec_idio[((kkk-1)*n_quad+1):(kkk*n_quad)] .= rf_vec
            big_weight_vec_idio[((kkk-1)*n_quad+1):(kkk*n_quad)] .= big_weight_vec .* idio_weight_vec[kkk]
        end
    end

    # ============================================================ !
    # Calculation
    # ============================================================ !
    ies   = ies_vec[iii]
    bbeta = bbeta_vec[iii]
    gma   = gma_vec[iii]
    thtbar  = thtbar_vec[iii]
    labor  = l_aggr*labor_alloc_vec[iii]
    labor_part =  1.0 .+ (1.0 ./ies .- 1.0) .* thtbar .* tht/(1.0 .+ tht) .* labor .^( (1.0 .+ tht) ./ tht )
    cons_temp = c_vec[iii]
    share_temp = share_vec[iii]

	temp_vec = -rk_vec./(rf_vec.-rk_vec .+ sqrt(eps()))  ## Yinjie 20220103 sqrt_eps = 1e-16
	temp_vec[temp_vec.>=0.0].=-10000.0
	share_low = maximum([low_guess_fixed, maximum(temp_vec)+sqrt(eps())])

	#             ! find natural leverage limit to prevent negative payoffs
	temp_vec = -rk_vec./(rf_vec .-rk_vec .+ sqrt(eps()))
	temp_vec[temp_vec.<=0.0].=10000.0

	share_high = minimum([high_guess_fixed, minimum(temp_vec) - sqrt(eps()) ])

	#             ! start from share guess but stay in limits
	share_temp = maximum([minimum([share_vec[iii], share_high]), share_low]) 


    diff_inner = 1.0;
counter = 0;
    while (diff_inner > 1E-05)
        # ============================================================ !
        # Update Consumption Based on Portfolio Choice
        # ============================================================ !


        curr_share = share_temp
        consumption = cons_temp

        r_alpha_vec = curr_share .* rf_vec .+ (1.0 .- curr_share) .* rk_vec

        savings = wealth .+ w_choice .* labor .- consumption


        if (use_idio_risk == 1)
            # ------------------------------------------------------------ !
            # With Idiosyncratic Risk
            # ------------------------------------------------------------ !
            # Expand portfolio return
            for kkk = 1:n_idio
                r_alpha_vec_idio[((kkk-1)*n_quad+1):(kkk*n_quad)] .= curr_share.*rf_vec .+ (1.0 .- curr_share).*rk_vec.*exp.(idio_shock_grid[kkk,iii])
            end

            # Other Operations
            next_period_share_idio = (savings .* r_alpha_vec_idio .+ q_l_nxt_idio[:,iii] ) ./ (tot_wealth_ss_vec[iii])
            	v_vec_twisted_idio = (v_idio[:,iii] .* next_period_share_idio) .^ (1.0-gma)
            	EV = sum(big_weight_vec_idio .* v_vec_twisted_idio)


				deriv_helper_idio = bbeta.*(v_idio[:,iii].^((1.0/ies) - gma )) .* (next_period_share_idio.^(-gma)) .*(EV .^( (gma - (1.0/ies))/(1.0-gma) )) .* mc_idio[:,iii]

				if (constr_agents[iii] == 1 && constraint_binding == 1)
					temp = sum(deriv_helper_idio .* rf_vec_idio .* big_weight_vec_idio)
				else
					temp = sum(deriv_helper_idio .* r_alpha_vec_idio .* big_weight_vec_idio)
				end
        else
            # ------------------------------------------------------------ !
            # No Idiosyncratic risk
            # ------------------------------------------------------------ !
            next_period_share = (savings .* r_alpha_vec .+ q_l_nxt[:,iii] ) ./ tot_wealth_ss_vec[iii]
			v_vec_twisted = (v_temp[:,iii] .* next_period_share) .^ (1.0-gma)
			EV = sum(big_weight_vec .* v_vec_twisted)

			deriv_helper = bbeta.*(v_temp[:,iii].^((1.0/ies) - gma )) .* (next_period_share.^(-gma)) .*(EV .^( (gma - (1.0/ies))/(1.0-gma) )) .* mc_temp[:,iii]
			

			if (constr_agents[iii] == 1 && constraint_binding == 1)
				temp = sum(deriv_helper .* rf_vec .* big_weight_vec)
			else
				temp = sum(deriv_helper .* r_alpha_vec .* big_weight_vec)	
			end
        end

		if temp>0
        	cons_update = labor_part .* (temp ) .^ (-ies)
        	cons_update = minimum([maximum([cons_update, min_cons_sav]), (wealth .+ w_choice .* labor .- min_cons_sav)])
		else
			cons_update = min_cons_sav
		end

        # ============================================================ !
        # Update Portfolio Choice
        # ============================================================ !
        if (constr_agents[iii] == 1 && constraint_binding == 1)
            share_update = 1.0 .- kbar_constr[iii] .* q_current ./ savings
        else

            # find natural leverage limit to prevent negative payoffs
            temp_vec = -rk_vec ./ (rf_vec .- rk_vec .+ sqrt(eps()))
            temp_vec[temp_vec .>= 0.0] .= -10000.0

            share_low = maximum([low_guess_fixed, maximum(temp_vec) .+ sqrt(eps())])

            # find natural leverage limit to prevent negative payoffs
            temp_vec = -rk_vec ./ (rf_vec .- rk_vec .+ sqrt(eps()))
            temp_vec[temp_vec .<= 0.0] .= 10000.0
            share_high = minimum([high_guess_fixed, minimum(temp_vec) .- sqrt(eps()) ])

            # Calculate Portfolio Share, Given Updated Consumption and Savings
            share_update = copy(curr_share)
            share_update = calc_portfolio_share(1.0, share_low, share_high, iii, big_weight_vec, savings, labor, cons_update,  0.0, q_l_nxt[:,iii], v_temp[:,iii], mc_temp[:,iii], rf_vec, rk_vec, q_current, share_update,param)
            # Here I am not very sure about the output, needs to be revisited
		end

        # ============================================================ !
        # Update Values
        # ============================================================ !
        diff_inner = abs.(consumption .- cons_update)
        cons_temp = consumption .+ 0.2 .* (cons_update .- consumption)
        share_temp = curr_share .+ 0.2 .* (share_update .- curr_share)

    end

    # ============================================================ !
    # Collect policies
    # ============================================================ !
    results_vec[1] = consumption
    results_vec[2] = share_update
    results_vec[3] = savings .* share_update
    results_vec[4] = savings .* (1.0 .- share_update) ./ q_current
	return results_vec
end
