function calc_excess_bond_nom(aggr_wealth::Float64, v_temp::Array{Float64,2},mc_temp::Array{Float64,2},
    rf::Array{Float64,1}, rk::Array{Float64,1}, q_current::Float64,w_choice::Float64, 
    l::Float64, next_k_value::Float64, q_l_nxt::Array{Float64,2}, sss::Int64, 
    big_weight_vec::Array{Float64,1}, c_temp::Array{Float64,1}, share_guess_vec::Array{Float64, 1},
    constraint_binding::Array{Int64,1},param::mpr_economy)

    @unpack_mpr_economy param

    #     ! Calculate bond demand for each agent 
    constraint_binding_new = Array{Int64,1}(undef,n_I);
    constraint_binding_new .= 0;
    k_temp_vec = zeros(n_I);
    b_temp_vec = zeros(n_I);
    for iii = 1:n_I
        #         ! Basic information for the agent 
        wealth = aggr_wealth*wealth_share_grid[iii,sss]
        savings =  wealth + w_choice*l*labor_alloc_vec[iii] - c_temp[iii]
        if savings <0.0
            println(savings, sss)
            error("found negative savings here")
        end

        #             ! find natural leverage limit to prevent negative payoffs
        temp = -rk./(rf.-rk .+ sqrt(eps()))  ## Yinjie 20220103 sqrt_eps = 1e-16
        temp[temp.>=0.0].=-10000.0
        share_low = maximum([low_guess_fixed, maximum(temp)+sqrt(eps())])

        #             ! find natural leverage limit to prevent negative payoffs
        temp = -rk./(rf .-rk .+ sqrt(eps()))
        temp[temp.<=0.0].=10000.0

        share_high = minimum([high_guess_fixed, minimum(temp) - sqrt(eps()) ])

        #             ! start from share guess but stay in limits
        share_guess = maximum([minimum([share_guess_vec[iii], share_high]), share_low]) 


        share_guess = calc_portfolio_share(1.0, share_low, share_high, iii, big_weight_vec, savings, l, c_temp[iii],  
                next_k_value, q_l_nxt[:,iii], v_temp[:,iii],mc_temp[:,iii], rf, rk, q_current, 
                share_guess,param)



        #             ! Store calculated portfolio shares and bond demand 
        share_guess_vec[iii]= share_guess

        if (constr_agents[iii] == 1) # occasionally binding constraint
          share_guess = 1.0 - kbar_constr[iii]*q_current/savings
          if (share_guess < share_guess_vec[iii])  
              constraint_binding_new[iii] = 1
          end
          if (constraint_binding[iii] == 1)
              share_guess_vec[iii] = share_guess 
          end
        end

        # Store calculated portfolio shares and bond demand 
        k_temp_vec[iii] = savings*(1.0-share_guess_vec[iii])/q_current
        b_temp_vec[iii] = share_guess_vec[iii]*savings
    end

    excess_b =  sum(lmbd_vec.*b_temp_vec) 
    return share_guess_vec, k_temp_vec, b_temp_vec, excess_b, constraint_binding_new
end
