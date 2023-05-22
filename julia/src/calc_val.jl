
function load_vars(folder,idx)
    next_state_mat = load(folder*"results.jld","next_state_mat")
    if idx == 0
    nxt_mat, k_next_mat, nom_i_mat, q_mat, c_mat, share_mat, infl_mat, l_aggr_mat, q_l_mat, w_choice_mat, v_mat, mc_mat, constraint_binding_mat,weights,smol_coeffs = load(folder*"variables.jld","nxt_mat","k_next_mat","nom_i_mat","q_mat","c_mat","share_mat","infl_mat","l_aggr_mat",
    "q_l_mat","w_choice_mat","v_mat","mc_mat","constraint_binding_mat","weights","smol_coeffs");
    return  next_state_mat, nxt_mat, k_next_mat ,nom_i_mat,q_mat ,c_mat ,share_mat ,infl_mat ,l_aggr_mat ,q_l_mat ,w_choice_mat, v_mat , mc_mat, constraint_binding_mat,weights,smol_coeffs
    else 
        nxt_mat, k_next_mat, nom_i_mat, q_mat, c_mat, share_mat, infl_mat, l_aggr_mat, q_l_mat, w_choice_mat, v_mat, mc_mat, constraint_binding_mat,weights,smol_coeffs,q_bond_mat = load(folder*"variables.jld","nxt_mat","k_next_mat","nom_i_mat","q_mat","c_mat","share_mat","infl_mat","l_aggr_mat",
    "q_l_mat","w_choice_mat","v_mat","mc_mat","constraint_binding_mat","weights","smol_coeffs","q_bond_mat");
    return  next_state_mat, nxt_mat, k_next_mat ,nom_i_mat,q_mat ,c_mat ,share_mat ,infl_mat ,l_aggr_mat ,q_l_mat ,w_choice_mat, v_mat , mc_mat, constraint_binding_mat,weights,smol_coeffs,q_bond_mat
    end
end 

function calc_val(mpr_vars::mpr_economy,mode_toirf::Int64,mode_tobond::Int64,param_id::Int64)
    # Mode==0 means to solve; Mode==1 means not to solve (skip iteration for equilibrium and bond prices)
    prog = ProgressUnknown("Diff is:"; showspeed=true)

    @unpack_mpr_economy mpr_vars
    n_nxt = 6 + 2*n_I;
    irf_ll_length = 20;
    maxiter = 5000;
    nrhs = n_nxt;
    ## Make sure to get precision small unless it's the debug mode
    conv_crit = 1e-6;
    constraint_binding = Array{Int64,1}(undef,n_I);
    constraint_binding_mat = Array{Int64,2}(undef,n_I, n_states);
    constraint_binding_monetary_mat = Array{Int64,2}(undef,n_I, n_states);
    constraint_binding .= 0;
    constraint_binding_mat .= 0;
    constraint_binding_monetary_mat .= 0;

    v_mat = Array{Float64, 2}(undef,n_I, n_states); 
    k_next_mat= Array{Float64, 2}(undef,n_quad, n_states);
    q_mat=Array{Float64, 1}(undef,n_states);
    l_aggr_mat=Array{Float64, 1}(undef,n_states);
    
    weights = zeros(n_states,n_nxt);
    smol_coeffs = zeros(n_states,nrhs);

    next_state_mat=zeros(n_quad, n_active_dims, n_states);
    next_state_mat_new=zeros(n_quad, n_active_dims, n_states);
    nxt_mat=zeros(n_quad, n_nxt,n_states);
    nxt_mat_trans=zeros(1, n_nxt,n_states);
    interp_input_mat=zeros(n_states, n_nxt);
    nom_i_mat=zeros(n_states);
    # nom_i = 0.0;
    w_choice_mat=zeros(n_states);
    infl_mat = zeros(n_states);
    c_mat=zeros(n_I, n_states);
    mc_mat = zeros(n_I, n_states);
    mc_mat_new = zeros(n_I, n_states);
    # mc_new = zeros(n_I);
    share_mat=zeros(n_I, n_states);
    q_l_mat=zeros(n_I, n_states);
    
    
    s_nxt=zeros(n_quad,3,n_states);
    v_mat_new=zeros(n_I, n_states);
    l_aggr_mat_new=zeros(n_states);
    q_mat_new=zeros(n_states);
    nom_i_mat_new=zeros(n_states);
    c_mat_new=zeros(n_I, n_states);
    w_choice_mat_new=zeros(n_states);
    k_next_mat_new=zeros(n_quad, n_states);
    q_l_mat_new=zeros(n_I, n_states);
    mc_mat_new = zeros(n_I, n_states);
    infl_mat_new=zeros(n_states);
    k_mat=zeros(n_states);
    share_mat_new=zeros(n_I, n_states);
    
    
    
    q_bond_mat=zeros(n_states, n_bond+1);
    results_vec=zeros( n_interp- n_bond*2,n_states);
    results_mat=zeros(n_states,  n_interp);
    nxt_bond_prices=zeros( n_quad,n_states);
    who_prices_bond= zeros(n_states);
    
    monetary_shock = 0.0;
    next_state_monetary_mat=zeros(n_quad, n_active_dims, n_states); 
    nom_i_monetary_mat=zeros(n_states);
    q_monetary_mat=zeros(n_states);
    c_monetary_mat=zeros(n_I,n_states);
    share_monetary_mat=zeros(n_I,n_states);
    infl_monetary_mat=zeros(n_states);
    l_aggr_monetary_mat=zeros(n_states);
    next_state_monetary_mat_new=zeros(n_quad, n_active_dims, n_states);
    q_l_monetary_mat=zeros(n_I,n_states);
    share_monetary_mat_new=zeros(n_I,n_states);
    w_choice_monetary_mat_new=zeros(n_states);
    infl_monetary_mat_new=zeros(n_states);
    v_monetary_mat=zeros(n_I,n_states);
    q_monetary_mat_new=zeros(n_states);
    l_aggr_monetary_mat_new=zeros(n_states);
    v_monetary_mat_new=zeros(n_I,n_states);
    nom_i_monetary_mat_new=zeros(n_states);
    c_monetary_mat_new=zeros(n_I,n_states);
    k_next_monetary_mat_new=zeros(n_quad,n_states);
    q_l_monetary_mat_new=zeros(n_I,n_states);
    next_state_montransition_mat_new=zeros(1,n_active_dims, n_states);
    next_state_montransition_mat=zeros(1,n_active_dims, n_states);
    k_next_monetary_mat=zeros(n_quad, n_states);
    w_choice_monetary_mat=zeros(n_states);
    trans_shock_vec=zeros(n_shocks);
    results_monetary_mat=zeros(n_states, n_interp);
    q_bond_monetary_mat=zeros(n_states, n_bond+1);
    mc_monetary_mat = zeros(n_I, n_states);
    mc_monetary_mat_new = zeros(n_I, n_states);
    
    

    shock_char = ""
    diff = 0.0;
    outer_iter =1;
    counter = Array{Int64,1}(undef,n_states)
    
    
    l_aggr_mat      .= l_ss
    q_mat           .= 1.0
    nom_i_mat       .= 1.0 + rf_ss
    infl_mat        .= 1.0
    share_mat       .= 0.0 
    q_l_mat[1,:]    .= q_l_ss_vec[1];
    q_l_mat[2,:]    .= q_l_ss_vec[2];
    q_l_mat[3,:]    .= q_l_ss_vec[3];
    
    
    if mode_toirf == 0
        if mode_tobond ==0
            counter = 0;
            # Value Function 
            for iii = 1:n_I
                v_mat[iii,:] .= v_ss_vec[iii]
                mc_mat[iii,:] .= mc_ss_vec[iii]  ## 0616
            end
            # Capital 
            k_mat = state_grid[idx_k,:]
            # Wage 
            w_choice_mat = (1.0-aalpha)* k_mat.^(aalpha) * l_ss^(-aalpha)
            # Consumption 
            for iii = 1:n_I
                c_mat[iii,:] = w_choice_mat .* l_ss .+ k_mat .* wealth_share_grid[iii,:] .* rf_ss
            end
            
            
            
            # (2) Initialize Next Period State Transition Matrix
            # Note that next_state_mat is on Smolyak grid (-1 to 1)
            next_state_mat .= 0.0
            for sss = 1:n_states
                counter = 0
            
            #   Initialize next period capital = current period capital * disaster shock
                k_next_mat[:, sss] .= state_grid[idx_k,sss] .* exp.(dz_vec) ./ exp.(dz_vec_adj)
            
                if (vector_mus_dimensions[1] > 0)
                    counter += 1
                    next_state_mat[:,counter,sss] .= (state_grid[idx_k,sss]./(exp.(dz_vec_adj)) .- k_grid_mean)./k_grid_dev
                end
            
                if (vector_mus_dimensions[2] > 0)
                    counter += 1
                    next_state_mat[:,counter,sss] .= 0.0
                end
            
                if (vector_mus_dimensions[3] > 0)
                    counter += 1
                    next_state_mat[:,counter,sss] .= 0.0
                end
            
                if (vector_mus_dimensions[4] > 0) 
                    counter += 1
                    next_state_mat[:,counter,sss] .= (next_m_mat[:,sss] .- m_grid_mean) ./m_grid_dev
                end
            
                if (vector_mus_dimensions[5] > 0) 
                    counter += 1
                    next_state_mat[:,counter,sss] .= (w_choice_mat[sss] .- w_grid_mean) ./w_grid_dev
                end
            
                if (vector_mus_dimensions[6] > 0) 
                    counter += 1
                    next_state_mat[:,counter,sss] .= (next_dis_mat[:,sss] .- dis_grid_mean) ./dis_grid_dev
                end
            
            
            end
            
            # Limit the grid points between -1 and 1
            
            next_state_mat[next_state_mat.>1.0] .= 1.0;
            next_state_mat[next_state_mat.<-1.0] .= -1.0;
            
            next_state_mat_new = copy(next_state_mat); # Pre-allocation
            
            
            outer_iter = 0;
            diff = 1.0;
            
            
            while (diff > conv_crit && outer_iter < maxiter)

                outer_iter += 1
            #     timer_start = omp_get_wtime()
            
            #         ************************************************************ !
            #         (1) Collect objects to be interpolated, and calculate Smolyak coefficients 
            #         ************************************************************ !
            #         9 objects needs interpolation: 3 value functions, q, labor, inflation
            #         3 time endowment prices
            #         ************************************************************ !
            #         Collect Objects
                for iii = 1:n_I
                    interp_input_mat[:,iii] .= v_mat[iii,:]
                    interp_input_mat[:,n_I+iii] .= mc_mat[iii,:]
                end
                interp_input_mat[:,2*n_I + 1]  = q_mat
                interp_input_mat[:,2*n_I + 2]  = l_aggr_mat
                interp_input_mat[:,2*n_I + 3]  = infl_mat
                interp_input_mat[:,2*n_I + 4]  = q_l_mat[1,:]
                interp_input_mat[:,2*n_I + 5]  = q_l_mat[2,:]
                interp_input_mat[:,2*n_I + 6]  = q_l_mat[3,:]


                smol_coeffs = smol_polynom \ interp_input_mat
                
            
            #         ! ************************************************************ !
            #         ! (2) Main updates state by state (parallelized)
            #         ! ************************************************************ !
            counter = Array{Int64,1}(undef,n_states)
            counter .= 0;
                
                    Threads.@threads for sss = 1:n_states
                    # for sss = 1:n_states

            #             ! ============================================================ !
            #             ! [1] Interpolate next period values from Smolyak Coefficients
            #             ! ============================================================ !
            #             ! Output is nxt_mat (n_quad * 9)
            #             ! ============================================================ !
            polyn_points = Smolyak_Polynomial2(next_state_mat[:,:,sss], n_active_dims, 
                                               n_quad, n_states, max_smol_level, smol_elem_ani)
                    
            nxt_mat[:,:,sss]  = polyn_points * smol_coeffs
            
            #             ============================================================ 
            #             [2] Calculate Equilibrium and Update
            #             ============================================================ 
            #             Load current states 

                share_mat_new[:,sss] , nom_i_mat_new[sss], q_mat_new[sss], q_l_mat_new[:,sss], k_next_mat_new[:,sss] , c_mat_new[:,sss], infl_mat_new[sss], v_mat_new[:,sss],
                mc_mat_new[:,sss] , l_aggr_mat_new[sss] , w_choice_mat_new[sss], s_nxt[:,:,sss], constraint_binding_mat[:,sss] = calc_equilibrium_and_update(nxt_mat[:,:,sss], k_next_mat[:,sss], sss,
                c_mat[:,sss], l_aggr_mat[sss], q_mat[sss], q_l_mat[:,sss], infl_mat[sss], nom_i_mat[sss], share_mat[:,sss], monetary_shock,
                constraint_binding_mat[:,sss], mpr_vars,outer_iter = outer_iter)
            
    
            
            #             ============================================================
            #             [3] Update next period state (next_state_mat_new)
            #             Note that it is on the Smolyak grid 
            #             ============================================================
            
                        if (vector_mus_dimensions[1] > 0) 
                            counter[sss] = counter[sss] + 1
                            next_state_mat_new[:, counter[sss], sss]   = (k_next_mat_new[:,sss] ./exp.(dz_vec_adj) .-k_grid_mean )./k_grid_dev
                        end
                        k_next_mat_new[:,sss]      = k_next_mat_new[:,sss]./exp.(dz_vec_adj).*exp.(dz_vec)
            
                        if (vector_mus_dimensions[2] > 0) 
                            counter[sss] = counter[sss] + 1
                            next_state_mat_new[:, counter[sss], sss] = (s_nxt[:,1,sss] .- s_a_grid_mean)./s_a_grid_dev
                        end
            
                        if (vector_mus_dimensions[3] > 0) 
                            counter[sss] = counter[sss] + 1
                            next_state_mat_new[:, counter[sss], sss] = (s_nxt[:,3,sss] .- s_c_grid_mean)./s_c_grid_dev
                        end
            
                        if (vector_mus_dimensions[4] > 0) 
                            counter[sss] = counter[sss] + 1
                            next_state_monetary_mat_new[:, counter[sss], sss] = (next_m_mat[:,sss] .+ rho_m.*monetary_shock   .- m_grid_mean) ./m_grid_dev
                        end
            
                        if (vector_mus_dimensions[5] > 0) 
                            counter[sss] = counter[sss] + 1
                            next_state_mat_new[:, counter[sss], sss]  = (w_choice_mat_new[sss]./exp.(dz_vec_adj) .-w_grid_mean) ./w_grid_dev
                        end
            
                        if (vector_mus_dimensions[6] > 0) 
                            counter[sss] = counter[sss] + 1
                        end
                        
                    end
                    
            
            #         ! ************************************************************ !
            #         ! (3) Update the matrices 
            #         ! ************************************************************ !
            #         ! Restrict Smolyak Grid
                    next_state_mat_new[next_state_mat_new.>1.0] .= 1.0;
                    next_state_mat_new[next_state_mat_new.<-1.0] .= -1.0;
                
            #         ! Calculate Differences
                    diff =  maximum([ maximum(abs.(log.(v_mat_new./v_mat))), 
                            maximum(abs.(log.(mc_mat_new./mc_mat))), 
                            maximum(abs.((q_l_mat_new.-q_l_mat))), 
                            maximum(abs.(log.(c_mat_new./c_mat))), 
                            maximum(abs.((share_mat_new.-share_mat))), 
                            maximum(abs.(log.(k_next_mat_new./k_next_mat))), 
                            maximum(abs.(log.(l_aggr_mat_new./l_aggr_mat))), 
                            maximum(abs.(log.(nom_i_mat_new./nom_i_mat))), 
                            maximum(abs.(log.(w_choice_mat_new./w_choice_mat))), 
                            maximum(abs.(log.(q_mat_new./q_mat))), 
                            maximum(abs.(log.(infl_mat_new./infl_mat))), 
                            maximum(abs.(next_state_mat .- next_state_mat_new)) ])
                    # sleep(0.5)
                    # ProgressMeter.next!(prog; showvalues = [(:outer_iter,outer_iter), (:diff,diff)])
                
                if (mod(outer_iter, 10) == 0)
            
                    println("")
                    println("Paramset $param_id, Iteration $outer_iter")
                    println("")
        
                    println("Changes (max diff: $diff )")
                    println("------------------------------------------------------------------------")
                    println("v         = ",maximum(abs.(log.(v_mat_new[1,:]./v_mat[1,:]))), " ",maximum(abs.(log.(v_mat_new[2,:]./v_mat[2,:]))),  
                    " ",maximum(abs.(log.(v_mat_new[3,:]./v_mat[3,:]))))
                    println("mc         = ",maximum(abs.(log.(mc_mat_new[1,:]./mc_mat[1,:]))), " ",maximum(abs.(log.(mc_mat_new[2,:]./mc_mat[2,:]))),  
                    " ",maximum(abs.(log.(mc_mat_new[3,:]./mc_mat[3,:]))))
                    println("q_l         = ",maximum(abs.((q_l_mat_new[1,:].-q_l_mat[1,:])))," ", maximum(abs.((q_l_mat_new[2,:].-q_l_mat[2,:]))),  
                    " ",maximum(abs.((q_l_mat_new[3,:].-q_l_mat[3,:]))))
                    println("c_spend   = ", maximum(abs.(log.(c_mat_new[1,:]./c_mat[1,:]))), " ", maximum(abs.(log.(c_mat_new[2,:]./c_mat[2,:]))),  
                    " ", maximum(abs.(log.(c_mat_new[3,:]./c_mat[3,:]))))
                    println("share_mat = " , maximum(abs.((share_mat_new[1,:].-share_mat[1,:])))," ", maximum(abs.((share_mat_new[2,:].-share_mat[2,:]))), 
                    " ",maximum(abs.((share_mat_new[3,:].-share_mat[3,:]))))
                    println("k         = " , maximum(abs.(log.(k_next_mat_new./k_next_mat))))
                    for sss = 1:n_active_dims 
                        println("next_state= ", maximum(abs.(next_state_mat_new[:,sss,:].-next_state_mat[:,sss,:])))
                    end
                    println("q         = ",  maximum(abs.(log.(q_mat_new./q_mat))))
                    println("nom_i     = ",  maximum(abs.(log.(nom_i_mat_new./nom_i_mat))))
                    println("l         = ", maximum(abs.(log.(l_aggr_mat_new./l_aggr_mat))))
                    println(" w         = ", maximum(abs.(log.(w_choice_mat_new./w_choice_mat))))
                    println("infl      = ", maximum(abs.(log.(infl_mat_new./infl_mat))))
        
                end
            
         
                
                    v_mat          = v_mat            .+ 1.0*(v_mat_new.-v_mat)
                    mc_mat         = mc_mat           .+ 1.0*(mc_mat_new-mc_mat)
                    next_state_mat = next_state_mat   .+  0.25*(next_state_mat_new.-next_state_mat);
                    q_mat          = q_mat            .+ 0.25*(q_mat_new.-q_mat);
                    l_aggr_mat     = l_aggr_mat       .+ 0.25*(l_aggr_mat_new.-l_aggr_mat); 
                    w_choice_mat   = w_choice_mat     .*(1 .+ 0.25*log.(w_choice_mat_new./w_choice_mat)); 
                    c_mat          = c_mat            .+ 0.5*(c_mat_new .- c_mat);
                    q_l_mat        = q_l_mat          .+ 0.3*(q_l_mat_new.-q_l_mat);
                    nom_i_mat      = nom_i_mat        .+ 1.0.*(nom_i_mat_new .- nom_i_mat);
                    infl_mat       = infl_mat         .+ 0.25.*(infl_mat_new .- infl_mat);
                    k_next_mat     = k_next_mat       .+ 0.25.*(k_next_mat_new .- k_next_mat);
                    share_mat      = share_mat        .+ 1.0.*(share_mat_new .- share_mat) ;
                
            end
            # ProgressMeter.finish!(prog)


            save(results_folder*"results.jld", "next_state_mat", next_state_mat)
            save(results_folder*"variables.jld","nxt_mat",nxt_mat,"k_next_mat",k_next_mat,"nom_i_mat",nom_i_mat,"q_mat",q_mat,
                    "c_mat",c_mat,"share_mat",share_mat,"infl_mat",infl_mat,"l_aggr_mat", l_aggr_mat,"q_l_mat",q_l_mat,
                    "w_choice_mat",w_choice_mat,"v_mat",v_mat,"mc_mat",mc_mat,"constraint_binding_mat",constraint_binding_mat,"weights",weights,"smol_coeffs",smol_coeffs)
        else
            next_state_mat, nxt_mat ,k_next_mat ,nom_i_mat,q_mat ,c_mat ,share_mat ,infl_mat ,l_aggr_mat ,q_l_mat ,w_choice_mat,
            v_mat ,mc_mat,constraint_binding_mat,weights,smol_coeffs = load_vars(results_folder,0)
        end

        q_bond_mat .= 1.0

        for  bbb = 1:n_bond
               
            bond_coeffs = smol_polynom \ q_bond_mat[:,bbb]
            # ! Get coefficients for interpolating bond prices
            

            Threads.@threads for sss = 1:n_states

                # Interpolation 
                polyn_points = Smolyak_Polynomial2(next_state_mat[:,:,sss], n_active_dims,
                                                   n_quad, n_states, max_smol_level, smol_elem_ani)

                nxt_mat[:,:,sss] = polyn_points * smol_coeffs

                # Bond Prices 
                nxt_bond_prices[:,sss] = polyn_points * bond_coeffs

                # ! Calculate Bond Prices 
                who_prices_bond[sss], q_bond_mat[sss,bbb+1], results_vec[:,sss]= calc_bond_prices(nxt_mat[:,:,sss], n_nxt, k_next_mat[:,sss], nxt_bond_prices[:,sss], sss, 
                c_mat[:,sss], l_aggr_mat[sss],q_mat[sss], q_l_mat[:,sss], infl_mat[sss], 
                nom_i_mat[sss], share_mat[:,sss],constraint_binding_mat[:,sss],mpr_vars)
                
                # ! Store in results matrix 
                results_mat[sss, 1:(n_interp-n_bond*2)]  = results_vec[:,sss]
                results_mat[sss, n_interp-n_bond*2+bbb]  = q_bond_mat[sss,bbb+1]
                results_mat[sss, n_interp-n_bond+bbb]    = who_prices_bond[sss]
            end
        end

        # ! ------------------------------------------------------------ !
        # ! Done with basic model solution here, store objects 
        # ! ------------------------------------------------------------ !
        println("DONE WITH MODEL SOLUTION - NOW STORE")
        # flush(6)
        save(results_folder*"results.jld", "next_state_mat", next_state_mat, "results_mat", results_mat)
        save(results_folder*"variables.jld","nxt_mat",nxt_mat,"k_next_mat",k_next_mat,"nom_i_mat",nom_i_mat,"q_mat",q_mat,
                "c_mat",c_mat,"share_mat",share_mat,"infl_mat",infl_mat,"l_aggr_mat", l_aggr_mat,"q_l_mat",q_l_mat, "q_bond_mat",q_bond_mat,
                "w_choice_mat",w_choice_mat,"v_mat",v_mat,"mc_mat",mc_mat,"constraint_binding_mat",constraint_binding_mat,"weights",weights,"smol_coeffs",smol_coeffs)
    else
    # ! print all results
        next_state_mat, nxt_mat ,k_next_mat ,nom_i_mat,q_mat ,c_mat ,share_mat ,infl_mat ,l_aggr_mat ,q_l_mat ,w_choice_mat,
            v_mat ,mc_mat, constraint_binding_mat,weights,smol_coeffs,q_bond_mat = load_vars(results_folder,1)
    end
    
    for fff = 1:size(irf_indices,1)
        if (irf_indices[fff] == sidx_m) 
            monetary_shock = irf_shock_sizes[irf_indices[fff]]
        else
            monetary_shock = 0.0
        end

        next_state_monetary_mat     = copy(next_state_mat)
        next_state_monetary_mat_new = copy(next_state_mat)
        k_next_monetary_mat         = copy(k_next_mat)
        nom_i_monetary_mat          = copy(nom_i_mat) 
        q_monetary_mat              = copy(q_mat)
        c_monetary_mat              = copy(c_mat )
        share_monetary_mat          = copy(share_mat)
        infl_monetary_mat           = copy(infl_mat )
        l_aggr_monetary_mat         = copy(l_aggr_mat) 
        q_l_monetary_mat            = copy(q_l_mat )
        w_choice_monetary_mat       = copy(w_choice_mat )
        v_monetary_mat              = copy(v_mat )
        mc_monetary_mat             = copy(mc_mat) 
        constraint_binding_monetary_mat = copy(constraint_binding_mat) 
        
        for iii = 1:n_I
            interp_input_mat[:,iii] = v_mat[iii,:]
            interp_input_mat[:,iii+n_I] = mc_mat[iii,:]
        end
        interp_input_mat[:,2*n_I + 1]  = q_mat
        interp_input_mat[:,2*n_I + 2]  = l_aggr_mat
        interp_input_mat[:,2*n_I + 3]  = infl_mat
        interp_input_mat[:,2*n_I + 4]  = q_l_mat[1,:]
        interp_input_mat[:,2*n_I + 5]  = q_l_mat[2,:]
        interp_input_mat[:,2*n_I + 6]  = q_l_mat[3,:]

        # ! Solve for Smolyak Coefficients (smol_coeffs)

        smol_coeffs = smol_polynom \ interp_input_mat

        outer_iter = 0
        diff = 1.0

        while (diff > conv_crit && outer_iter < maxiter) 
            outer_iter = outer_iter + 1
        
            counter .= 0;

            

            Threads.@threads for sss = 1:n_states
            #     ! ============================================================ !
            #     ! [1] Interpolate next period values from Smolyak Coefficients
            #     ! ============================================================ !
            #     ! Output is nxt_mat (n_quad * 9)
            #     ! ============================================================ !

                polyn_points = Smolyak_Polynomial2(next_state_monetary_mat[:,:,sss], n_active_dims, n_quad, n_states, max_smol_level, smol_elem_ani)
                nxt_mat[:,:,sss] = polyn_points * smol_coeffs
                nxt_mat_2 = k_next_monetary_mat[:,sss]
                # ! ============================================================ !
                # ! [2] Calculate Equilibrium and Update
                # ! ============================================================ !
                # ! Load current states 

                share_monetary_mat_new[:,sss] , nom_i_monetary_mat_new[sss],q_monetary_mat_new[sss], q_l_monetary_mat_new[:,sss],
                k_next_monetary_mat_new[:,sss] , c_monetary_mat_new[:,sss], infl_monetary_mat_new[sss], v_monetary_mat_new[:,sss] ,
                mc_monetary_mat_new[:,sss] , l_aggr_monetary_mat_new[sss] , w_choice_monetary_mat_new[sss],
                s_nxt[:,:,sss] = calc_equilibrium_and_update(nxt_mat[:,:,sss], k_next_monetary_mat[:,sss], sss,
                c_monetary_mat[:,sss], l_aggr_monetary_mat[sss], q_monetary_mat[sss], q_l_monetary_mat[:,sss],
                infl_monetary_mat[sss], nom_i_monetary_mat[sss], share_monetary_mat[:,sss],
                monetary_shock, constraint_binding_monetary_mat[:,sss], mpr_vars)

        

                # ! Find updated variables
                # ! ============================================================ !
                # ! [3] Update next period state (next_state_mat_new)
                # ! Note that it is on the Smolyak grid 
                # ! ============================================================ !

                # ! Example: maps k_next_new calculated above to Smolyak grid 
                if (vector_mus_dimensions[1] > 0) 
                    counter[sss] = counter[sss] + 1
                    next_state_monetary_mat_new[:, counter[sss], sss]   = (k_next_monetary_mat_new[:,sss]./exp.(dz_vec_adj) .-k_grid_mean )./k_grid_dev
                end

                if (vector_mus_dimensions[2] > 0) 
                    counter[sss] = counter[sss] + 1
                    next_state_monetary_mat_new[:, counter[sss], sss]  = (s_nxt[:,1,sss] .- s_a_grid_mean)./s_a_grid_dev
                end

                if (vector_mus_dimensions[3] > 0) 
                    counter[sss] = counter[sss] + 1
                    next_state_monetary_mat_new[:, counter[sss], sss]  = (s_nxt[:,3,sss] .- s_c_grid_mean)./s_c_grid_dev
                end

                if (vector_mus_dimensions[4] > 0) 
                    counter[sss] = counter[sss] + 1
                    next_state_monetary_mat_new[:, counter[sss], sss]  = (next_m_mat[:,sss] .+ rho_m .*monetary_shock   .- m_grid_mean) ./m_grid_dev
                end

                if (vector_mus_dimensions[5] > 0) 
                    counter[sss]= counter[sss] + 1
                    next_state_monetary_mat_new[:, counter[sss], sss]   = (w_choice_monetary_mat_new[sss]./exp.(dz_vec_adj) .-w_grid_mean) ./w_grid_dev
                end

                if (vector_mus_dimensions[6] > 0) 
                    counter[sss] = counter[sss] + 1
                end

                k_next_monetary_mat_new[:,sss]      = k_next_monetary_mat_new[:,sss]./exp.(dz_vec_adj).*exp.(dz_vec)

            end

            # !$OMP END DO 
            # !$OMP END PARALLEL


            next_state_monetary_mat_new[next_state_monetary_mat_new .> 1.0]   .= 1.0 
            next_state_monetary_mat_new[next_state_monetary_mat_new.<-1.0]   .= -1.0

            # ! check convergence
            if (outer_iter < 10) 
            diff = 1.0
            else
            # ! check convergence
            diff =  maximum([ maximum(abs.(log.(v_monetary_mat_new./v_monetary_mat))),
            maximum(abs.((mc_monetary_mat_new.-mc_monetary_mat))), 
            maximum(abs.((q_l_monetary_mat_new.-q_l_monetary_mat))), 
            maximum(abs.(log.(c_monetary_mat_new./c_monetary_mat))), 
            maximum(abs.((share_monetary_mat_new.-share_monetary_mat))), 
            maximum(abs.(log.(k_next_monetary_mat_new./k_next_monetary_mat))), 
            maximum(abs.(log.(l_aggr_monetary_mat_new./l_aggr_monetary_mat))), 
            maximum(abs.(log.(nom_i_monetary_mat_new./nom_i_monetary_mat))), 
            maximum(abs.(log.(w_choice_monetary_mat_new./w_choice_monetary_mat))), 
            maximum(abs.(log.(q_monetary_mat_new./q_monetary_mat))), 
            maximum(abs.(log.(infl_monetary_mat_new./infl_monetary_mat))), 
            maximum(abs.(next_state_monetary_mat .- next_state_monetary_mat_new)) ])
        
            end

            println( " ")
            println( "  Monetary shock iteration ", outer_iter)
            println( " Changes (max diff: ", diff , ")")
            println( " ")

            # ! Damping updates
            v_monetary_mat           = @. v_monetary_mat + 1.0*(v_monetary_mat_new-v_monetary_mat )
            mc_monetary_mat           = @. mc_monetary_mat + 1.0*(mc_monetary_mat_new-mc_monetary_mat )

            for sss = 1:n_states
                for i = 1:n_quad
                    for j = 1:n_active_dims
                        if abs((next_state_monetary_mat[i,j,sss]-next_state_monetary_mat_new[i,j,sss])) > 0.05
                            next_state_monetary_mat[i,j,sss] = next_state_monetary_mat[i,j,sss] + 0.01*sign(next_state_monetary_mat_new[i,j,sss] - next_state_monetary_mat[i,j,sss])
                        else
                            next_state_monetary_mat[i,j,sss] = next_state_monetary_mat[i,j,sss] + 0.2*(next_state_monetary_mat_new[i,j,sss] - next_state_monetary_mat[i,j,sss])
                        end
                    end
                end
            end

            q_monetary_mat = @. q_monetary_mat + 0.2*(q_monetary_mat_new-q_monetary_mat)
            l_aggr_monetary_mat = @.  l_aggr_monetary_mat + 0.2*(l_aggr_monetary_mat_new-l_aggr_monetary_mat)     
            w_choice_monetary_mat = @.  w_choice_monetary_mat*(1 + 0.2*log(w_choice_monetary_mat_new/w_choice_monetary_mat))
            q_l_monetary_mat = @.  q_l_monetary_mat  + 0.5*(q_l_monetary_mat_new-q_l_monetary_mat)
            c_monetary_mat = @.  c_monetary_mat      + 0.5*(c_monetary_mat_new - c_monetary_mat)
            nom_i_monetary_mat      = @.  nom_i_monetary_mat           + 0.2*(nom_i_monetary_mat_new - nom_i_monetary_mat)
            infl_monetary_mat      =  @.  infl_monetary_mat           + 0.2*(infl_monetary_mat_new - infl_monetary_mat)
            k_next_monetary_mat     = @.  k_next_monetary_mat          + 0.2*(k_next_monetary_mat_new - k_next_monetary_mat)
            share_monetary_mat      = @.  share_monetary_mat           + 1.0*(share_monetary_mat_new - share_monetary_mat) 

        end

        


        # ! Obtain Results 
        q_bond_monetary_mat .= 1.0
        # !DIR$ NOUNROLL
        for bbb = 1:n_bond

            bond_coeffs = smol_polynom \ q_bond_mat[:,bbb]

            Threads.@threads for sss = 1:n_states

                polyn_points = Smolyak_Polynomial2(next_state_monetary_mat[:,:,sss], n_active_dims, 
                                                   n_quad, n_states, max_smol_level, smol_elem_ani)

                nxt_mat[:,:,sss] = polyn_points * smol_coeffs

                nxt_mat_2 = k_next_monetary_mat[:,sss]

                nxt_bond_prices[:,sss] = polyn_points * bond_coeffs
                # println(nxt_bond_prices[:,sss])


                who_prices_bond[sss], q_bond_monetary_mat[sss,bbb+1], results_vec[:,sss]= calc_bond_prices(nxt_mat[:,:,sss], n_nxt, k_next_monetary_mat[:,sss], nxt_bond_prices[:,sss], sss, 
                c_monetary_mat[:,sss], l_aggr_monetary_mat[sss],q_monetary_mat[sss], q_l_monetary_mat[:,sss], infl_monetary_mat[sss], 
                nom_i_monetary_mat[sss], share_monetary_mat[:,sss],constraint_binding_monetary_mat[:,sss],mpr_vars)                
                

                results_monetary_mat[sss, 1:(n_interp-n_bond*2)]  = results_vec[:,sss]
                results_monetary_mat[sss, n_interp-n_bond*2+bbb]  = q_bond_monetary_mat[sss,bbb+1]
                results_monetary_mat[sss, n_interp-n_bond+bbb]    = who_prices_bond[sss]
            end
        end

    

        trans_shock_vec .= 0.0
        if (irf_indices[fff] == sidx_z) 
            trans_shock_vec[irf_indices[fff]] = irf_shock_sizes[fff]
        elseif (irf_indices[fff] == sidx_dis)
            trans_shock_vec[irf_indices[fff]] = irf_shock_sizes[fff]
        elseif (irf_indices[fff] == sidx_m) 
            trans_shock_vec[irf_indices[fff]] = 0*irf_shock_sizes[fff]
        end

        next_state_montransition_mat[1,:,:]     = next_state_mat[no_shock_idx,:,:]
        next_state_montransition_mat_new[1,:,:] = next_state_mat[no_shock_idx,:,:]

        for iii = 1:n_I
            interp_input_mat[:,iii] = v_monetary_mat[iii,:]
            interp_input_mat[:,iii+n_I] = mc_monetary_mat[iii,:]
        end
        interp_input_mat[:,2*n_I + 1]  = q_monetary_mat
        interp_input_mat[:,2*n_I + 2]  = l_aggr_monetary_mat
        interp_input_mat[:,2*n_I + 3]  = infl_monetary_mat
        interp_input_mat[:,2*n_I + 4]  = q_l_monetary_mat[1,:]
        interp_input_mat[:,2*n_I + 5]  = q_l_monetary_mat[2,:]
        interp_input_mat[:,2*n_I + 6]  = q_l_monetary_mat[3,:]

        # ! solve for polynomial coefficients
        
        diff = 1.0
        outer_iter = 0

        smol_coeffs = smol_polynom \ interp_input_mat


        while (diff > conv_crit && outer_iter < maxiter)
            outer_iter = outer_iter + 1
            counter .= 0

            Threads.@threads for sss = 1:n_states

                polyn_points_trans = Smolyak_Polynomial2(next_state_montransition_mat[:,:,sss], n_active_dims, 
                                                   1, n_states, max_smol_level, smol_elem_ani)

                nxt_mat_trans[:,:,sss] = polyn_points_trans * smol_coeffs

                # ! ============================================================ !
                # ! [1] Interpolate next period values from Smolyak Coefficients
                # ! ============================================================ !
                # ! Output is nxt_mat (n_quad * 9)
                # ! ============================================================ !
                

                # ! ============================================================ !
                # ! [2] Calculate Equilibrium and Update
                # ! ============================================================ !
                # ! Load current states 

                # ! Find updated variables
                s_nxt[1,:,sss] = calc_unexpected_transition(nxt_mat_trans[:,:,sss], n_nxt, k_next_mat[1,sss], sss, 
                                c_mat[:,sss], l_aggr_mat[sss],q_mat[sss], q_l_mat[:,sss], 
                                infl_mat[sss], nom_i_mat[sss], share_mat[:,sss], trans_shock_vec, s_nxt[1,:,sss],mpr_vars)

                if (vector_mus_dimensions[1] > 0) 
                    counter[sss] = counter[sss] + 1
                    next_state_montransition_mat_new[1, counter[sss], sss]   = (k_next_mat[1,sss]./exp.(trans_shock_vec[sidx_z]) .-k_grid_mean )./k_grid_dev
                end

                if (vector_mus_dimensions[2] > 0) 
                    counter[sss] = counter[sss] + 1
                    next_state_montransition_mat_new[1, counter[sss], sss] = (s_nxt[1,1,sss] .- s_a_grid_mean)./s_a_grid_dev
                end

                if (vector_mus_dimensions[3] > 0) 
                    counter[sss] = counter[sss] + 1
                    next_state_montransition_mat_new[1, counter[sss], sss] = (s_nxt[1,3,sss] .- s_c_grid_mean)./s_c_grid_dev
                end

                if (vector_mus_dimensions[4] > 0) 
                    counter[sss] = counter[sss] + 1

                    next_state_montransition_mat_new[1, counter[sss], sss] = (next_m_mat[no_shock_idx,sss] .+ trans_shock_vec[sidx_m]   .- m_grid_mean) ./m_grid_dev
                end

                if (vector_mus_dimensions[5] > 0) 
                    counter[sss] = counter[sss] + 1
                    next_state_montransition_mat_new[1, counter[sss], sss]  = (w_choice_mat[sss]./exp.(trans_shock_vec[sidx_z]) .-w_grid_mean) ./w_grid_dev
                end

                if (vector_mus_dimensions[6] > 0) 
                    counter[sss] = counter[sss] + 1
                    if (dis_grid_dev < sqrt(eps()))
                        next_state_montransition_mat_new[1, counter[sss], sss] = 0.0  .+ trans_shock_vec[sidx_dis]./dis_grid_dev
                    else
                        next_state_montransition_mat_new[1, counter[sss], sss] = (next_dis_mat[no_shock_idx,sss] .+ trans_shock_vec[sidx_dis]   .- dis_grid_mean) ./dis_grid_dev
                    end
                end

            end

 
            next_state_montransition_mat_new[next_state_montransition_mat_new .> 1.0] .= 1.0
            next_state_montransition_mat_new[next_state_montransition_mat_new .< -1.0] .= -1.0


            if (outer_iter < 10) 
                diff = 1.0
            else
                diff =  maximum(abs.(next_state_montransition_mat .- next_state_montransition_mat_new)) 
            end

            next_state_montransition_mat = @. next_state_montransition_mat +  0.2*(next_state_montransition_mat_new-next_state_montransition_mat)

        end

        println("DONE WITH IRF CALCULATION - NOW STORE")
        # flush(6)

        shock_char=string(irf_indices[fff])

        save(results_folder*"/result_irf_"*shock_char*".jld","next_state_monetary_mat",next_state_monetary_mat,
        "next_state_montransition_mat",next_state_montransition_mat,"results_monetary_mat",results_monetary_mat)
    end
end


