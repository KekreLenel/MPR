function create_results(param::mpr_economy,paramid::Int64)
    
    @unpack_mpr_economy param

    LSEED = 1
    SEED = Array{Int64,1}(undef,LSEED)
    n_spread = 20;
    state_spread_vec = zeros(n_spread)
    next_state_monetary_mat = zeros(n_quad, n_active_dims, n_states)
    next_state_montransition_mat = zeros(1, n_active_dims, n_states)
    interp_input_monetary_mat = zeros(n_states, n_interp)
    state_variables = zeros(smolyak_d)
    other_variables = zeros(n_interp)
    current_state = zeros(1,n_active_dims)
    next_state = zeros(1,n_active_dims)
    stochsteady_state = zeros(1,n_active_dims)
    big_weight_vec = zeros(n_quad)
    polyn_points = zeros(n_quad,n_states)

    println("")
    println("READING RESULTS")

    # ------------------------------------------------------------ !
    # Store Basic Grid Sizes, etc.
    # ------------------------------------------------------------ !

    next_state_mat = load(results_folder*"results.jld","next_state_mat")
    interp_input_mat = load(results_folder*"results.jld","results_mat")

    save(results_folder*"smol_grid.jld","smol_grid",smol_grid)
    save(results_folder*"num_params.jld", "n_I",n_I,"n_states", n_states,"n_active_dims",n_active_dims,
    "n_interp",n_interp,"n_shocks",n_shocks,"n_spread", n_spread,"smolyak_d",smolyak_d, "n_prm", n_prm,
    "n_sim_periods", n_sim_periods, "n_irf_periods", n_irf_periods)

    save(results_folder*"grid_locs.jld","k_grid_mean", k_grid_mean,"k_grid_dev",k_grid_dev,"s_a_grid_mean",s_a_grid_mean,"s_a_grid_dev",s_a_grid_dev,"s_c_grid_mean", s_c_grid_mean,
    "s_c_grid_dev",s_c_grid_dev,"w_grid_mean",w_grid_mean,"w_grid_dev",w_grid_dev,"dis_grid_mean",dis_grid_mean,
    "dis_grid_dev",dis_grid_dev, "m_grid_mean",m_grid_mean,"m_grid_dev", m_grid_dev)

    # ------------------------------------------------------------ !
    # Find Stochastic Steady State
    # ------------------------------------------------------------ !
    println("")
    println("Find STOCHASTIC STEADY STATE")

    # (1.1) Polynomial coefficients for transition matrix
    nrhs = n_active_dims

    smol_montransition_coeffs=zeros(1,n_states,nrhs)
    smol_monetary_coeffs=zeros(n_quad,n_states,nrhs)
    smol_coeffs=zeros(n_quad,n_states,nrhs)

    for qqq = 1:n_quad
        smol_coeffs[qqq,:,:] = smol_polynom\transpose(next_state_mat[qqq,:,:])
    end

    # (1.2) Polynomial coefficients for non-state variables
    nrhs2 = n_interp
    interp_coeffs = smol_polynom \interp_input_mat
    interp_monetary_coeffs = zeros(n_states,nrhs2)

    # (1.3) Polynomial coefficients for state variables
    state_coeffs = smol_polynom \transpose(state_grid)

    # (2.1) Start from some interior grid point and move forward
    # Use no_shock_idx: the quadrature point where there's no shock happening
    current_state[1,:] .= 0
    diff = 1.0


    while (diff > sqrt(eps()))
        polyn_points = Smolyak_Polynomial(current_state, n_active_dims,max_smol_level, smol_elem_ani)
        next_state = polyn_points*smol_coeffs[no_shock_idx,:,:]
        next_state[next_state.>1.0] .= 1.0
        next_state[next_state.<-1.0] .= -1.0

        diff = maximum(abs.(current_state .- next_state))
        current_state = copy(next_state)

        
    end

    stochsteady_state = copy(current_state)


    other_variables = polyn_points*interp_coeffs
    state_variables = polyn_points*state_coeffs


    # (2.2) Print the Stochastic SS
    println("Stochastic steady state -------------------------------------")
    println("-------------------------------------------------------------")

    println("Capital   : $(state_variables[idx_k])")
    println("Wealth   a: $(state_variables[idx_sa])")
    println("Wealth   c: $(state_variables[idx_sc])")
    println("m         : $(state_variables[idx_m])")
    println("w         : $(state_variables[idx_w])")
    println("Disaster p: $(state_variables[idx_dis])")

    # ------------------------------------------------------------ !
    # Policy and Price ~ State
    # ------------------------------------------------------------ !
    println("")
    println("State dependencies ------------------------------------------")
    println("-------------------------------------------------------------")

    state_spread_series=zeros(n_spread,n_interp,smolyak_d)
    state_spread_vec = collect(LinRange(-1.0, 1.0, n_spread))
    counter = 0
    for qqq = 1:smolyak_d
        if (vector_mus_dimensions[qqq] > 0)
        counter = counter + 1
        for ttt = 1:n_spread
            for aaa = 1:n_active_dims
                if (aaa == counter)
                    current_state[1,aaa] = state_spread_vec[ttt]
                else
                    current_state[1,aaa] = stochsteady_state[1,aaa]
                end
            end

            polyn_points = Smolyak_Polynomial(current_state, n_active_dims, max_smol_level, smol_elem_ani)
            state_spread_series[ttt,:,qqq] = polyn_points*interp_coeffs

        end
        end
    end

    save(results_folder*"state_spread_series.jld","state_spread_series",state_spread_series)


    # ------------------------------------------------------------ !
    # Business Cycle Moments, No Disaster
    # ------------------------------------------------------------ !
    println("")
    println("Business cycle moments  no disaster -------------------------")
    println("-------------------------------------------------------------")


    n_periods = n_sim_periods + n_burn
    simshock_series = zeros(n_periods)
    shock_series = zeros(n_periods, n_shocks)
    state_series = zeros(n_periods, smolyak_d)
    other_vars_series = zeros(n_periods, n_interp)
    state_collection = zeros(n_periods, n_active_dims)

    # seed random number generator
    Random.seed!(123)
    SEED    = 712
    # generate random numbers:
    ifail = 0
    shock_generator = Uniform()
    simshock_series =  rand(shock_generator,n_periods)

    current_state = stochsteady_state

    shock_series[1,:] .= shock_grid[no_shock_idx[1],:]

    # loop through
    for ttt = 1:n_periods

        if mod(ttt,5000)==0
            println(ttt," periods finished." )
        end 

        polyn_points = Smolyak_Polynomial(current_state, n_active_dims, max_smol_level, smol_elem_ani)

        other_vars_series[ttt,:] = polyn_points*interp_coeffs
        state_series[ttt,:] = polyn_points*state_coeffs

        checker = 0.0
        q_nxt   = 0
        while (checker <= simshock_series[ttt] && q_nxt < n_quad-1)
            q_nxt   = q_nxt + 1
            checker = checker .+ quad_weight_vec[q_nxt]
        end

        next_state = polyn_points*smol_coeffs[q_nxt,:,:]


        if (ttt<n_periods)
            shock_series[ttt+1,:] = shock_grid[q_nxt,:]
        end

        next_state[next_state.>1.0] .= 1.0
        next_state[next_state.<-1.0] .= -1.0

        state_collection[ttt,:]   .= current_state[1,:]
        current_state             = copy(next_state)

    end

    for idx = 1:length(other_vars_series[:,67])
        if other_vars_series[idx,67] > 0.5 
            other_vars_series[idx,36] = 0.0
            other_vars_series[idx,45] = 0.0
            other_vars_series[idx,54] = 0.0
            other_vars_series[idx,57] = lmbd_vec[3]*kbar_constr[3]
        end
    end

    save(results_folder*"sim_series.jld", "state_series",  state_series[n_burn+1:n_periods,:],
    "other_vars_series", other_vars_series[n_burn+1:n_periods,:], 
    "shock_series",shock_series[n_burn+1:n_periods,:])





    # ------------------------------------------------------------ !
    # Sample n_sample IRF starting points
    # ------------------------------------------------------------ !
    println("")
    println("Sampling IRF Starting Points -------------------------")
    println("-------------------------------------------------------------")

    # Generate random integers to sample starting points
    ifail = 0
    starting_locs=zeros(n_sample)
    starting_states=zeros(n_sample, n_active_dims)

    sample_length = n_sim_periods - n_burn
    sample_vector = zeros(sample_length)
    sample_vector = [ttt for ttt = n_burn+1:n_sim_periods]

    starting_locs = sample(Random.GLOBAL_RNG,sample_vector, n_sample; replace=false, ordered=true)

    # Obtain starting states
    if (n_sample > 1)
    starting_states = state_collection[starting_locs,:]
    else
    starting_states = copy(stochsteady_state)
    end

    save(results_folder*"starting_states.jld", "starting_states",  starting_states)  ## 0729
    

    # ------------------------------------------------------------ !
    # Business Cycle Moments, With Disaster
    # ------------------------------------------------------------ !
    println("")
    println("Business cycle moments  w. disaster -------------------------")
    println("-------------------------------------------------------------")

    # generate random numbers:
    ifail = 0
    simshock_series = rand(shock_generator,n_periods)

    current_state = copy(stochsteady_state)

    shock_series[1,:] .= shock_grid[no_shock_idx[1],:]

    # loop through
    for ttt = 1:n_periods
        if mod(ttt,5000)==0
            println(ttt," periods finished." )
        end 

        polyn_points = Smolyak_Polynomial(current_state, n_active_dims, max_smol_level, smol_elem_ani) 
        other_vars_series[ttt,:] = polyn_points*interp_coeffs
        state_series[ttt,:] = polyn_points*state_coeffs

        p_dis = exp(state_series[ttt,idx_dis])

        # current weight vec
        # CHANGEd: Only last state is disaster
        big_weight_vec[1:n_quad-1] = quad_weight_vec * (1.0-p_dis)
        big_weight_vec[n_quad]    =  p_dis

        checker = 0.0
        q_nxt   = 0
        while (checker <= simshock_series[ttt] && q_nxt < n_quad)
            q_nxt   = q_nxt .+ 1
            checker = checker .+ big_weight_vec[q_nxt]
        end

        next_state = polyn_points*smol_coeffs[q_nxt,:,:]


        if (ttt<n_periods)
        shock_series[ttt+1,:] .= shock_grid[q_nxt,:]
        end

        next_state[next_state.>1.0] .= 1.0
        next_state[next_state.<-1.0] .= -1.0

        current_state        = copy(next_state)

    end

    for idx = 1:length(other_vars_series[:,67])
        if other_vars_series[idx,67] > 0.5 
            other_vars_series[idx,36] = 0.0
            other_vars_series[idx,45] = 0.0
            other_vars_series[idx,54] = 0.0
            other_vars_series[idx,57] = lmbd_vec[3]*kbar_constr[3]
        end
    end

    save(results_folder*"simdis_series.jld", "state_series",  state_series[n_burn+1:n_periods,:],
    "other_vars_series", other_vars_series[n_burn+1:n_periods,:], 
    "shock_series",shock_series[n_burn+1:n_periods,:])


    # allocate various matrices

    n_periods = n_irf_periods

    state_series_N=zeros(n_sample, n_periods, smolyak_d)
    other_vars_series_N=zeros(n_sample, n_periods, n_interp)
    shock_series_N=zeros(n_sample, n_periods, n_shocks)
    state_series_N_helper=zeros(n_sample, n_periods, smolyak_d)
    other_vars_series_N_helper=zeros(n_sample, n_periods, n_interp)
    shock_series_N_helper=zeros(n_sample, n_periods, n_shocks)
    irf_factor=zeros(n_sample)
    state_series_M=zeros(n_sample, n_periods, smolyak_d)
    other_vars_series_M=zeros(n_sample, n_periods, n_interp)
    shock_series_M=zeros(n_sample, n_periods, n_shocks)


    state_series=zeros(n_periods, smolyak_d)
    other_vars_series=zeros(n_periods, n_interp)
    shock_series=zeros(n_periods, n_shocks)

    println("Calc IRF paths absent any shock")
    for iii = 1:n_sample

        if (mod(iii, 50) == 0)
            println(iii)
        end

        # Starting State
        current_state[1,:] .= starting_states[iii,:]

        # loop through
        for ttt = 1:n_periods

            shock_series[ttt,:] .= shock_grid[no_shock_idx[1],:]

            polyn_points = Smolyak_Polynomial(current_state, n_active_dims,   max_smol_level, smol_elem_ani)


            next_state = polyn_points*smol_coeffs[no_shock_idx,:,:]
            other_vars_series[ttt,:] = polyn_points*interp_coeffs
            state_series[ttt,:] = polyn_points*state_coeffs
            current_state        = copy(next_state)

        end

        for idx = 1:length(other_vars_series[:,67])
            if other_vars_series[idx,67] > 0.5 
                other_vars_series[idx,36] = 0.0
                other_vars_series[idx,45] = 0.0
                other_vars_series[idx,54] = 0.0
                other_vars_series[idx,57] = lmbd_vec[3]*kbar_constr[3]
            end
        end



        # Store
        state_series_N[iii,:,:] .= state_series
        other_vars_series_N[iii,:,:] .= other_vars_series
        shock_series_N[iii,:,:] .= shock_series

        for ttt = 1 : n_periods
        state_series_N_helper[iii,ttt,:] .= state_series[1,:]
        other_vars_series_N_helper[iii,ttt,:] .= other_vars_series[1,:]
        shock_series_N_helper[iii,ttt,:] .= shock_series[1,:]
        end

    end

    state_series      = (sum(state_series_N,dims = 1)[1,:,:]) ./ n_sample
    other_vars_series = (sum(other_vars_series_N,dims = 1)[1,:,:]) ./ n_sample
    shock_series      = (sum(shock_series_N,dims = 1)[1,:,:]) ./ n_sample

    save(results_folder*"none_irf_series.jld","state_series",state_series,"other_vars_series",other_vars_series,"shock_series",shock_series)

    println("Done")
    println("")


    println("Calc IRF paths with mit shock")
    for fff = 1 : size(irf_indices,1)

        # ------------------------------------------------------------ !
        # Loop Each IRF Calculation
        # ------------------------------------------------------------ !


        println("")
        print("Computing Impulse responses -------------------------")
        println("shock_char",irf_indices[fff])
        println("-------------------------------------------------------------")

        shock_char=string(irf_indices[fff])
        # Get all solution files
        next_state_monetary_mat = load(results_folder*"result_irf_"*shock_char*".jld","next_state_monetary_mat")
        next_state_montransition_mat = load(results_folder*"result_irf_"*shock_char*".jld","next_state_montransition_mat")
        interp_input_monetary_mat = load(results_folder*"result_irf_"*shock_char*".jld","results_monetary_mat")

        for qqq = 1:n_quad
            b = transpose(next_state_monetary_mat[qqq,:,:])
            smol_monetary_coeffs[qqq,:,:] = smol_polynom\b
        end
        

        b =  transpose(next_state_montransition_mat[1,:,:])
        smol_montransition_coeffs[1,:,:] = smol_polynom\b



        b2 = interp_input_monetary_mat
        interp_monetary_coeffs = smol_polynom\b2

        for iii = 1:n_sample

            if (mod(iii, 500) == 0)
                println(iii)
            end

            # ************************************************************ !
            # IRF With Shock
            # ************************************************************ !

            # Starting State
            current_state[1,:] .= starting_states[iii,:]

            polyn_points = Smolyak_Polynomial(current_state, n_active_dims, max_smol_level, smol_elem_ani) 
            other_vars_series[1,:] = polyn_points*interp_coeffs
            state_series[1,:] = polyn_points*state_coeffs
            next_state = polyn_points*smol_coeffs[no_shock_idx,:,:]
            current_state = copy(next_state)
            
            polyn_points = Smolyak_Polynomial(current_state, n_active_dims, max_smol_level, smol_elem_ani) 

            other_vars_series[2,:] = polyn_points*interp_coeffs
            state_series[2,:] = polyn_points*state_coeffs
            next_state = polyn_points*smol_montransition_coeffs[1,:,:]

            # in that state find state and variables
            current_state = copy(next_state)
            polyn_points = Smolyak_Polynomial(current_state, n_active_dims, max_smol_level, smol_elem_ani) 
            other_vars_series[3,:] = polyn_points*interp_monetary_coeffs
            state_series[3,:] = polyn_points*state_coeffs
            next_state = polyn_points*smol_monetary_coeffs[no_shock_idx,:,:]

            shock_series[1,:] = shock_grid[no_shock_idx[1],:]
            shock_series[2,:] = shock_grid[no_shock_idx[1],:]
            shock_series[3,:] = shock_grid[no_shock_idx[1],:]

            if (irf_indices[fff] == sidx_z)
                shock_series[3,sidx_z] = irf_shock_sizes[fff]
            end

            current_state = copy(next_state)

            # loop through
            for ttt = 4 : n_periods

                if (ttt<=n_periods)
                    shock_series[ttt,:] .= shock_grid[no_shock_idx[1],:]
                end

                polyn_points = Smolyak_Polynomial(current_state, n_active_dims, max_smol_level, smol_elem_ani) 
                next_state = polyn_points*smol_coeffs[no_shock_idx,:,:]
                other_vars_series[ttt,:] = polyn_points*interp_coeffs
                state_series[ttt,:] = polyn_points*state_coeffs
                current_state        = copy(next_state)

            end

            if (sidx_m == irf_indices[fff]) 
                irf_factor[iii] =  0.002/((log(other_vars_series[3,jx_1y]) - log(other_vars_series_N[iii,3,jx_1y])) - (log(other_vars_series[2,jx_1y]) - log(other_vars_series_N[iii,2,jx_1y])))
            else 
                irf_factor[iii]  = 1.0
            end
 
            if (sidx_m == irf_indices[fff] && fff == 2)  
                save(results_folder*"irf_factor.jld","irf_factor",sum(irf_factor)/n_sample)
            end
         
            # interpolation of MPKs of constrained agent does not work well due to sparse grid
            # assign heuristically
            for idx = 1:length(other_vars_series[:,67])
                if other_vars_series[idx,67] > 0.5 
                    other_vars_series[idx,36] = 0.0
                    other_vars_series[idx,45] = 0.0
                    other_vars_series[idx,54] = 0.0
                    other_vars_series[idx,57] = lmbd_vec[3]*kbar_constr[3]
                end
            end
             
         
            # Store
            state_series_M[iii,:,:]      = (state_series .- state_series_N[iii,:,:]).*irf_factor[iii]          .+ state_series_N_helper[iii,:,:]
            other_vars_series_M[iii,:,:] = (other_vars_series .- other_vars_series_N[iii,:,:]).*irf_factor[iii] .+ other_vars_series_N_helper[iii,:,:]
            shock_series_M[iii,:,:]      = (shock_series .- shock_series_N[iii,:,:]).*irf_factor[iii]          .+ shock_series_N_helper[iii,:,:]
        end

        state_series      = (sum(state_series_M,dims=1)[1,:,:]) ./ n_sample
        other_vars_series = (sum(other_vars_series_M,dims = 1)[1,:,:]) ./ n_sample
        shock_series      = (sum(shock_series_M,dims = 1)[1,:,:]) ./ n_sample

        if (irf_indices[fff] == sidx_m && fff == 2) 

            shock_name = "m"

        elseif (irf_indices[fff] == sidx_z) 

            shock_name = "g"

            println(other_vars_series[1:20,1]')
            println(shock_series[1:20,1])

        elseif (irf_indices[fff] == sidx_dis)  

            shock_name = "p"

        end

        save(results_folder*shock_name*"_irf_series.jld","state_series",state_series,"other_vars_series",other_vars_series,"shock_series",shock_series)


    end

end