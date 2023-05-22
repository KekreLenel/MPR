function calc_portfolio_share(time_choice::Float64, share_low::Float64, share_high::Float64, iii::Int64, big_weight_vec::Array{Float64,1}, 
    savings::Float64, labor::Float64, consumption::Float64, next_k_value::Float64, q_l_nxt::Array{Float64,1}, 
    v::Array{Float64,1}, mc::Array{Float64,1}, rf::Array{Float64,1}, rk::Array{Float64,1}, q_current::Float64, share_guess::Float64, param::mpr_economy)
    @unpack_mpr_economy param
    
    v_use=zeros(n_quad*n_idio);
    mc_use=zeros(n_quad*n_idio);
    q_l_nxt_use=zeros(n_quad*n_idio);
    rf_use=zeros(n_quad*n_idio);
    rk_use=zeros(n_quad*n_idio);
    big_weight_vec_use=zeros(n_quad*n_idio);
    r_alpha=zeros(n_quad*n_idio);
    next_period_share=zeros(n_quad*n_idio);
    v_vec_twisted=zeros(n_quad*n_idio);
    dz_vec_use=zeros(n_quad*n_idio);
    share_guess_D = 0;
    share_guess_S = 0;

    gma = gma_vec[iii]
    ies = ies_vec[iii]
    bbeta = bbeta_vec[iii]

    #     ! ------------------------------------------------------------ !
    #     ! Initialize: expand grid
    #     ! ------------------------------------------------------------ !
    if use_idio_risk == 1
        #         ! Expand the grid
        for kkk = 1:n_idio
            #             ! Basic 
            v_use[((kkk-1)*n_quad+1):(kkk*n_quad)] .= v
            mc_use[((kkk-1)*n_quad+1):(kkk*n_quad)] .= mc
            rf_use[((kkk-1)*n_quad+1):(kkk*n_quad)] .= rf
            dz_vec_use[((kkk-1)*n_quad+1):(kkk*n_quad)] .= dz_vec

            #             ! risky assets
            q_l_nxt_use[((kkk-1)*n_quad+1):(kkk*n_quad)] .= q_l_nxt*exp(idio_shock_grid[kkk,iii])
            rk_use[((kkk-1)*n_quad+1):(kkk*n_quad)] .= rk*exp(idio_shock_grid[kkk,iii]) 

            #             ! Weights 
            big_weight_vec_use[((kkk-1)*n_quad+1):(kkk*n_quad)] .= big_weight_vec .* idio_weight_vec[kkk]
        end
    else
        #         ! No Idiosyncratic risk
        v_use = v 
        mc_use = mc 
        q_l_nxt_use = q_l_nxt 
        rf_use = rf 
        rk_use = rk 
        big_weight_vec_use = big_weight_vec
        dz_vec_use = dz_vec
    end


    # Start Solving
    share_temp = share_low
    r_alpha = share_temp.*rf_use .+ (1.0 .-share_temp).*rk_use
    next_period_share = (savings.*r_alpha  .+ q_l_nxt_use.*time_choice) ./ (tot_wealth_ss_vec[iii])
    v_vec_twisted = v_use.^(1.0/ies -gma).*mc_use
    E_MR = sum(v_vec_twisted.*(next_period_share.^(-gma)).*(rf_use .-  rk_use).*big_weight_vec_use)./abs.(sum(big_weight_vec_use.*v_vec_twisted))

    share_FOC_low = E_MR  

    #     ! share_FOC<0 means optimum lies to the left
    #     ! share_FOC>0 means optimum lies to the right 
    if share_FOC_low <= 0.0
        #         ! No solution --> just use lower value 
        share_guess = share_low
    else
        #         ! If not satisfied, use brent method to get port. share 

        #         ! ================================================== !
        #         ! (1) Evaluate FOC at high end 
        #         ! ================================================== !
        share_temp = share_high
        r_alpha = share_temp.*rf_use .+ (1.0 .-share_temp).*rk_use
        next_period_share = (savings.*r_alpha  .+ q_l_nxt_use.*time_choice) ./ (tot_wealth_ss_vec[iii])
        v_vec_twisted = v_use.^(1.0/ies .-gma).*mc_use
        E_MR = sum(v_vec_twisted.*(next_period_share.^(-gma)).*(rf_use .-  rk_use).*big_weight_vec_use)./abs.(sum(big_weight_vec_use.*v_vec_twisted))

        share_FOC_hi = E_MR    

        #         ! ================================================== !
        #         ! (2) Use Brent Method 
        #         ! ================================================== !

        if share_FOC_hi >= 0.0
            #             ! No solution --> just use higher value 
            share_guess = share_high
            # constraint_binding = -1
        else
            #! Evaluate FOC at current guess 
            share_temp = share_guess
            r_alpha = share_temp.*rf_use .+ (1.0 .-share_temp).*rk_use
            next_period_share = (savings.*r_alpha  .+ q_l_nxt_use.*time_choice) ./ (tot_wealth_ss_vec[iii])
            v_vec_twisted = v_use.^(1.0/ies .-gma).*mc_use
            E_MR = sum(v_vec_twisted.*(next_period_share.^(-gma)).*(rf_use .-  rk_use).*big_weight_vec_use)./abs.(sum(big_weight_vec_use.*v_vec_twisted))

            share_FOC = E_MR   

            #             ! Update using Brent Method 
            if share_FOC != 0.0
                share_guess_B = share_guess # ! make this highest value
                share_FOC_B   = share_FOC
                    

                if share_FOC <= 0.0 # ! share_FOC<0 means optimum lies to the left
                    share_guess_A = share_low
                    share_temp = share_guess_A
                    r_alpha = share_temp.*rf_use .+ (1.0 .-share_temp).*rk_use
                    next_period_share = (savings.*r_alpha  .+ q_l_nxt_use.*time_choice) ./ (tot_wealth_ss_vec[iii])
                    v_vec_twisted = v_use.^(1.0/ies .-gma).*mc_use
                    E_MR = sum(v_vec_twisted.*(next_period_share.^(-gma)).*(rf_use .-  rk_use).*big_weight_vec_use)./abs.(sum(big_weight_vec_use.*v_vec_twisted))
                    
                    share_FOC_A = E_MR   
                else # ! share_FOC > 0 so optimum lies to the right
                    share_guess_A = share_high
                    share_temp = share_guess_A
                    r_alpha = share_temp.*rf_use .+ (1.0 .-share_temp).*rk_use
                    next_period_share = (savings.*r_alpha  .+ q_l_nxt_use.*time_choice) ./ (tot_wealth_ss_vec[iii])
                    v_vec_twisted = v_use.^(1.0/ies .-gma).*mc_use
                    E_MR = sum(v_vec_twisted.*(next_period_share.^(-gma)).*(rf_use .-  rk_use).*big_weight_vec_use)./abs.(sum(big_weight_vec_use.*v_vec_twisted))

                    share_FOC_A = E_MR   
                end

                    #                 ! check that root is bracketed
                if (share_FOC_A*share_FOC_B) > 0
                    error("ERROR: Initial bracket does not contain root.")
                end

                #                 ! swap a and b as needed 
                if  abs(share_FOC_A) < abs(share_FOC_B) 
                    temp = share_guess_A
                    share_guess_A = share_guess_B
                    share_guess_B = temp

                    temp = share_FOC_A
                    share_FOC_A = share_FOC_B
                    share_FOC_B = temp
                end
                #                 ! Brent Solution 
                share_guess_C    = share_guess_A
                share_FOC_C      = share_FOC_A
                share_FOC_S      = share_FOC_B
                mflag = 1 # ! set mflag
                brent_delta = 5e-16
                while ( share_FOC_S != 0) && (abs(share_guess_A - share_guess_B) > 1e-14)
                    #   ! ************************************************************ !
                    #   ! Obtain updated position
                    #   ! ************************************************************ !
                    if ( (share_FOC_A != share_FOC_C) && (share_FOC_C != share_FOC_B) ) #  ! inverse quadratic interpolation
                        share_guess_S = share_guess_A * share_FOC_B * share_FOC_C / ( (share_FOC_A - share_FOC_B) * (share_FOC_A - share_FOC_C) ) +
                            share_guess_B * share_FOC_A * share_FOC_C / ( (share_FOC_B - share_FOC_A) * (share_FOC_B - share_FOC_C) ) + 
                            share_guess_C * share_FOC_A * share_FOC_B / ( (share_FOC_C - share_FOC_A) * (share_FOC_C - share_FOC_B) )
                    else #! secant method
                        share_guess_S = share_guess_B - share_FOC_B * (share_guess_B - share_guess_A) /   (share_FOC_B - share_FOC_A)
                    end

                    if ( ( (  ( share_guess_S > ( 3*share_guess_A + share_guess_B )/4 && share_guess_S < share_guess_B) || 
                    ( share_guess_S < ( 3*share_guess_A + share_guess_B )/4 && share_guess_S > share_guess_B)  ) == false ) || 
                    (mflag == 1 && abs(share_guess_S - share_guess_B) >= abs(share_guess_B - share_guess_C)/2  )             || 
                    (mflag == 0 && abs(share_guess_S - share_guess_B) >= abs(share_guess_B - share_guess_D)/2  )             ||
                    (mflag == 1 && abs(share_guess_B - share_guess_C) <  abs(brent_delta)  )                     ||
                    (mflag == 0 && abs(share_guess_B - share_guess_D) <  abs(brent_delta)  )   ) 

                        share_guess_S = (share_guess_A + share_guess_B )/ 2
                        mflag = 1
                    else
                        mflag = 0
                    end

                    #                     ! ************************************************************ !
                    #                     ! Evaluate FOC at share_guess_S
                    #                     ! ************************************************************ !
                    share_temp = share_guess_S
                    r_alpha = share_temp.*rf_use .+ (1.0 .-share_temp).*rk_use
                    next_period_share = (savings.*r_alpha  .+ q_l_nxt_use.*time_choice) ./ (tot_wealth_ss_vec[iii])
                    v_vec_twisted = v_use.^(1.0/ies .-gma).*mc_use
                    E_MR = sum(v_vec_twisted.*(next_period_share.^(-gma)).*(rf_use .-  rk_use).*big_weight_vec_use)./abs.(sum(big_weight_vec_use.*v_vec_twisted))
                        
                    share_FOC_S = E_MR

                    #                     ! ************************************************************ !
                    #                     ! Update Boundaries
                    #                     ! ************************************************************ !
                    share_guess_D = share_guess_C
                    share_guess_C = share_guess_B
                    share_FOC_C = share_FOC_B

                    if ((share_FOC_A*share_FOC_S) < 0) 
                    share_guess_B = share_guess_S
                    share_FOC_B = share_FOC_S
                    else
                    share_guess_A = share_guess_S
                    share_FOC_A = share_FOC_S
                    end

                    # ! swap a and b as needed 
                    if ( abs(share_FOC_A) < abs(share_FOC_B) ) 
                        temp = share_guess_A
                        share_guess_A = share_guess_B
                        share_guess_B = temp

                        temp = share_FOC_A
                        share_FOC_A = share_FOC_B
                        share_FOC_B = temp
                    end                  
                end

                    #! EXPORT 
                share_guess = share_guess_S
            end
        end
    end

    return share_guess
end 


