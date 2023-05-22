# ----------------------------------------------------------------- #
# 1. Read in parameters
# ----------------------------------------------------------------- #
function init(param::prm)
    @unpack_prm param

    smolyak_d =  6  
    n_I =  3  
    n_shocks =  3                     # number of states  number of agents  number of shocks
    n_uni_quad_g =  3  
    n_uni_quad_m =  3  
    n_uni_quad_dis =  3   # number of quadrature nodes in each dimension
    n_GH_points =  n_uni_quad_g*n_uni_quad_m*n_uni_quad_dis   # total number of quadrature shocks
    n_quad =  n_GH_points + 1  


    # States: 1. Capital          2. Wealth share a  3. Wealth share c  
    #         4. monetary policy  5. Wage            6. Disaster
    idx_k= 1;  idx_sa= 2;  idx_sc= 3;  
    idx_m= 4; idx_w= 5;  idx_dis= 6;

    # Shocks: 1. Growth  2. Monetary  3. Disaster 
    sidx_z =  1;  sidx_m =  2;   sidx_dis =  3;
    irf_indices = [1 ,2 , 3];

    lmbd_vec         = [lmbd_a, lmbd_b, 1.0-lmbd_a-lmbd_b] ;
    bbeta_vec        = [bbeta_a, bbeta_b,bbeta_c];
    gma_vec          = [gma_a , gma_b, gma_c];
    ies_vec          = [ies_a,  ies_b, ies_c];
    s_bar_vec           = [s_bar_a, 1.0 - s_bar_a - s_bar_c,s_bar_c];
    s_trgt_vec          = [s_trgt_a,1.0 - s_trgt_a - s_trgt_c,s_trgt_c];
    labor_alloc_vec     = [labor_alloc_a, labor_alloc_b,labor_alloc_c];
    kbar_constr         = [kbar,kbar,kbar];
    constr_agents       = [Int(constrained_a), Int(constrained_b),Int(constrained_c)];
    std_idio_risk_a     = idio_risk_a;
    std_idio_risk_b     = idio_risk_b;
    std_idio_risk_c     = idio_risk_c;
    s_a_grid_dev        = s_a_dev;
    s_c_grid_dev        = s_c_dev;
    ddelta = delta;
    tht = theta;
    thtbar_vec = zeros(n_I);
    v_normalization_vec = zeros(n_I);
    v_ss_vec = zeros(n_I);
    mc_ss_vec = zeros(n_I);
    tot_wealth_ss_vec = zeros(n_I);
    q_l_ss_vec = zeros(n_I);

    vector_mus_dimensions = [3, 3, 3, 3, 3, 3];
    # wage state variable not needed without stickiness
    if (chiW < 1E-3)
        vector_mus_dimensions[idx_w] = 0;
    end

    # monetary state variable needed for solving the model with monetary shock persistence 
    if (rho_m > 1E-6)
        vector_mus_dimensions[idx_m] = 3;
    end

    # disaster shock standard deviation
    sig_dis = disast_std*sqrt(1.0-rho_p^2)

    # min/max leverage constraints for numerical stability
    low_guess_fixed  = -20.0
    high_guess_fixed =   3.0

    println("---------------------------------------------")
    println("PARAMETRIZATION")
    println("---------------------------------------------")
    println("")
    println(" bbeta        = " , bbeta_vec)
    println(" gma          = " , gma_vec)
    println(" ies          = " , ies_vec)
    println(" s_bar        = " , s_bar_vec) 
    println(" s_target     = " , s_trgt_vec)
    println(" lmbd_vec     = " , lmbd_vec)
    println(" tht          = " , theta)
    println(" ddelta       = " , delta)
    println(" aalpha       = " , aalpha)
    println(" sig_g        = " , sig_z ) 
    println(" phi          = " , phi  )
    println(" tayl_ic      = " , tayl_ic)
    println(" sig_m        = " , sig_m  )
    println(" disast_p     = " , disast_p ) 
    println(" varphi_p     = " , varphi_p  )
    println(" rho_p        = " , rho_p  )
    println(" disast_std   = " , disast_std  )
    println(" vareps_w     = " , vareps_w  )
    println(" tau_w        = " , tau_w  )
    println(" chiW         = " , chiW  )
    println(" chiX         = " , chiX   ) 
    println(" xi           = " , xi    )       
    println(" l_target     = " , l_target   )  
    println("")
    println("---------------------------------------------")

    # ----------------------------------------------------------------- #
    # 2. Construct Quadratures
    # Ordering: g, m, dis
    # ----------------------------------------------------------------- #
    # (1) Variance-Covariance Matrix and its Cholesky Decomposition
    cov_mat = zeros(n_shocks,n_shocks);
    cov_mat[1,1] = sig_z^2.0   
    cov_mat[2,2] = sig_m^2.0
    cov_mat[3,3] = sig_dis^2.0
    cov_mat[:] = cholesky(cov_mat).U';  

    if (n_uni_quad_g > 1) 
        quad_vec_g, uni_weight_vec_g = get_quadrature_points(n_uni_quad_g);
    else
        quad_vec_g = 0.0
        uni_weight_vec_g = 1.0
    end

    if (n_uni_quad_m > 1) 
        quad_vec_m, uni_weight_vec_m = get_quadrature_points(n_uni_quad_m);
    else
        quad_vec_m = 0.0
        uni_weight_vec_m = 1.0
    end

    if (n_uni_quad_dis > 1) 
        quad_vec_dis, uni_weight_vec_dis = get_quadrature_points(n_uni_quad_dis);
    else
        quad_vec_dis = 0.0
        uni_weight_vec_dis = 1.0
    end
    counter = 0
    shock_grid = zeros(n_quad,n_shocks)
    quad_weight_vec = zeros(n_GH_points)
    for ggg = 1:n_uni_quad_g
        for mmm = 1:n_uni_quad_m
            for ddd = 1:n_uni_quad_dis
                counter = counter + 1
                quad_vec_temp  = [ quad_vec_g[ggg] quad_vec_m[mmm] quad_vec_dis[ddd]]
                
                for iii = 1:n_shocks
                    shock_grid[counter, iii] = sum(quad_vec_temp*cov_mat[iii,:]) 
                    # println(shock_grid[counter, iii])
                end

                ### weights of each realization of shock
                quad_weight_vec[counter] =  uni_weight_vec_g[ggg]*uni_weight_vec_m[mmm]*uni_weight_vec_dis[ddd]
            end
        end
    end

    # (4) Extend the shock grid to include disaster
    shock_grid[n_quad, 1] = -varphi_p
    shock_grid[n_quad, 2:n_shocks] .= 0.0

    # Growth vector dz_vec 
    dz_vec             = shock_grid[:,1] 

    # adjusted growth vector: dz_vec without disaster shock
    # needed for rescaling wage and capital for stationary solution 
    dz_vec_adj         = copy(dz_vec)
    dz_vec_adj[n_quad] = 0.0 

    # find shock index with no shock for stochastic steady state calculations
    no_shock_idx = sum(argmin(abs.(shock_grid[:,1])+abs.(shock_grid[:,2])+abs.(shock_grid[:,3]),dims = 1))
    # impulse response shock sizes
    irf_shock_sizes = [IRF_g*sig_z, IRF_m*sig_m, IRF_dis*sig_dis]

    n_idio = use_idio_risk + 1
    quad_vec_idio = zeros(n_idio)
    idio_weight_vec = zeros(n_idio)
    idio_shock_grid = zeros(n_idio, n_I)

    if (use_idio_risk == 1)
        quad_vec_idio, idio_weight_vec = get_quadrature_points(n_idio);
    else
        quad_vec_idio .= 0.0;
        idio_weight_vec .= 1.0;
    end


    idio_shock_grid[:,1] = std_idio_risk_a .* quad_vec_idio .- log(sum(idio_weight_vec.*exp.(quad_vec_idio.*std_idio_risk_a)))
    idio_shock_grid[:,2] = std_idio_risk_b .* quad_vec_idio .- log(sum(idio_weight_vec.*exp.(quad_vec_idio.*std_idio_risk_b)))
    idio_shock_grid[:,3] = std_idio_risk_c .* quad_vec_idio .- log(sum(idio_weight_vec.*exp.(quad_vec_idio.*std_idio_risk_c)))
    
    #---------------------------------------------------------#
    # 3. Create Smolyak Grid
    #---------------------------------------------------------#
    n_active_dims = count(vector_mus_dimensions .> 0);
    mus_dimensions_redux = zeros(n_active_dims);
    counter = 0;
    for ddd = 1:smolyak_d 
        if (vector_mus_dimensions[ddd] > 0)
            counter = counter + 1
            mus_dimensions_redux[counter] = vector_mus_dimensions[ddd]
        end
    end

    smolyak_elem_iso = Smolyak_Elem_Isotrop(n_active_dims, Int(maximum(mus_dimensions_redux)))
    smol_elem_ani    = Smolyak_Elem_Anisotrop(smolyak_elem_iso, n_active_dims, mus_dimensions_redux)
    smol_grid        = Smolyak_Grid(n_active_dims,Int(maximum(mus_dimensions_redux)), smol_elem_ani)
    smol_polynom     = Smolyak_Polynomial(smol_grid,n_active_dims,Int(maximum(mus_dimensions_redux)), smol_elem_ani); 
    n_states  = size(smol_grid,1)


    wealth_share_grid=zeros(n_I,n_states);
    state_grid = zeros(smolyak_d, n_states);
    next_m_mat = zeros(n_quad, n_states);
    next_dis_mat = zeros(n_quad, n_states);

    # (4) Set mean and standard deviations for most states 
    #     Standard Deviations
    m_grid_dev     = 2.5*sig_m/sqrt(1.0-rho_m^2)
    dis_grid_dev   = 2.5*disast_std
    #  Means
    s_a_grid_mean = s_trgt_vec[1]
    s_c_grid_mean = s_trgt_vec[3]
    m_grid_mean       = 0.0
    dis_grid_mean     = disast_p


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

    # labor endowment price
    q_l_ss_vec = w_ss*l_ss*labor_alloc_vec ./ (1.0 .-bbeta_vec )

    # total individual wealth including endowment claim
    tot_wealth_ss_vec = q_l_ss_vec .+ per_person_wealth*k_ss*(1.0+rk_ss)

    # ----------------------------------------------------------------- #
    # Value function; labor endowment prices
    # ----------------------------------------------------------------- #
    for iii = 1:n_I
        util, mc_ss_vec[iii],~ =  util_fun(c_ss_vec[iii], l_ss*labor_alloc_vec[iii], iii,thtbar_vec, ies_vec, tht)
        v_ss_vec[iii]  =  util^( 1.0/(1.0-1.0/ies_vec[iii]) )
    end  

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
    println( " mc_ss =            " , mc_ss_vec[1], " ", mc_ss_vec[2], " ", mc_ss_vec[3])
    println("---------------------------------------------")
    println("")

    # Export calibrated thtbar for result files
    save(results_folder*"extra_data.jld","thtbar_vec", thtbar_vec);

    counter = 0
        
    if (vector_mus_dimensions[1]>0) 
    counter = counter + 1
    state_grid[1,:] = @. smol_grid[:,counter]*k_grid_dev + k_grid_mean
    else 
    state_grid[1,:] .= k_grid_mean
    end

    if (vector_mus_dimensions[2]>0)
    counter = counter + 1
    state_grid[2,:] = @. smol_grid[:,counter]*s_a_grid_dev + s_a_grid_mean
    else 
    state_grid[2,:] .= s_a_grid_mean
    end

    if (vector_mus_dimensions[3]>0) 
    counter = counter + 1
    state_grid[3,:] = @. smol_grid[:,counter]*s_c_grid_dev + s_c_grid_mean
    else 
    state_grid[3,:] .= s_c_grid_mean
    end

    if (vector_mus_dimensions[4]>0)
    counter = counter + 1
    state_grid[4,:] = @. smol_grid[:,counter]*m_grid_dev + m_grid_mean
    else 
    state_grid[4,:] .= m_grid_mean
    end

    if (vector_mus_dimensions[5]>0)
    counter = counter + 1
    state_grid[5,:] = @. smol_grid[:,counter]*w_grid_dev + w_grid_mean
    else 
    state_grid[5,:] .= w_grid_mean
    end

    if (vector_mus_dimensions[6]>0) 
    counter = counter + 1
    state_grid[6,:] = @. smol_grid[:,counter]*dis_grid_dev + dis_grid_mean
    else 
    state_grid[6,:] .= dis_grid_mean
    end


    # ----------------------------------------------------------------- #
    # Derive wealth share grids from state grids
    # ----------------------------------------------------------------- #
    wealth_share_grid[1,:] = @. state_grid[idx_sa,:] / lmbd_vec[1]
    wealth_share_grid[2,:] = @. (1.0-state_grid[idx_sa,:]-state_grid[idx_sc,:]) / lmbd_vec[2]
    wealth_share_grid[3,:] = @. state_grid[idx_sc,:] / lmbd_vec[3]

    # ----------------------------------------------------------------- #
    # Next period position matrix 
    # Can be constructed directly for the 2 exogenous states 
    # ----------------------------------------------------------------- #
    for sss = 1:n_states
        next_m_mat[:,sss]   = @. (1.0 - rho_m)*0.0  + rho_m*state_grid[idx_m,sss]   + shock_grid[:,sidx_m]
        next_dis_mat[:,sss]  = @. (1.0 - rho_p)*disast_p + rho_p*state_grid[idx_dis,sss]  + shock_grid[:,sidx_dis]
    end


    return mpr_economy{Int64,Float64}(
        no_shock_idx = no_shock_idx,
        lmbd_vec         = lmbd_vec ,
        bbeta_vec        = bbeta_vec,
        gma_vec          = gma_vec,
        ies_vec          = ies_vec,
        tht                 = theta,
        ddelta              = ddelta,
        aalpha              = aalpha,
        sigma_z             = sig_z,
        chiX                = chiX,
        phi                 = phi,
        tayl_ic             = tayl_ic,
        sig_m               = sig_m,
        rho_m               = rho_m,
        disast_p            = disast_p,
        varphi_p            = varphi_p,
        rho_p               = rho_p,
        disast_std          = disast_std,
        vareps_w            = vareps_w,
        tau_w               = tau_w,
        chiW                = chiW,
        s_bar_vec           = s_bar_vec,
        s_trgt_vec          = s_trgt_vec,
        xi                  = xi,
        l_target            = l_target,
        labor_alloc_vec     = labor_alloc_vec,
        kbar_constr         = kbar_constr,
        gov_debt            = gov_debt,
        k_grid_adj          = k_grid_adj,
        w_grid_adj          = w_grid_adj,
        s_a_grid_dev        = s_a_dev,
        s_c_grid_dev        = s_c_dev,
        k_dev_param         = k_dev_param,
        w_dev_param         = w_dev_param,
        IRF_z               = IRF_g,
        IRF_m               = IRF_m,
        IRF_dis             = IRF_dis,
        irf_indices         = irf_indices,
        irf_shock_sizes     = irf_shock_sizes,
        constr_agents       = constr_agents,
        use_idio_risk       = use_idio_risk,
        std_idio_risk_a     = idio_risk_a,
        std_idio_risk_b     = idio_risk_b,
        std_idio_risk_c     = idio_risk_c,
        thtbar_vec          = thtbar_vec,
        v_normalization_vec =  v_normalization_vec,
        vector_mus_dimensions = vector_mus_dimensions,
        low_guess_fixed     = low_guess_fixed,
        high_guess_fixed    = high_guess_fixed,
        smol_grid           = smol_grid,
        state_grid          = state_grid,
        smol_polynom        = smol_polynom,
        smolyak_elem_iso    = smolyak_elem_iso,
        smol_elem_ani       = smol_elem_ani,
        wealth_share_grid   = wealth_share_grid,
        s_a_grid_mean       = s_a_grid_mean,
        s_c_grid_mean       = s_c_grid_mean,
        m_grid_mean         = m_grid_mean ,
        dis_grid_mean       = dis_grid_mean,
        n_active_dims       = n_active_dims,
        n_states            = n_states,
        n_idio              = n_idio,
        sig_dis             = sig_dis,
        k_grid_dev          = k_grid_dev,
        m_grid_dev          = m_grid_dev,
        w_grid_dev          = w_grid_dev,
        dis_grid_dev        = dis_grid_dev,
        k_grid_mean         = k_grid_mean,
        w_grid_mean         = w_grid_mean,
        shock_grid          = shock_grid,
        quad_weight_vec     = quad_weight_vec,
        dz_vec              = dz_vec,
        dz_vec_adj          = dz_vec_adj,
        next_m_mat          = next_m_mat,
        next_dis_mat        = next_dis_mat,
        quad_vec_idio       = quad_vec_idio,
        idio_weight_vec     = idio_weight_vec,
        idio_shock_grid     = idio_shock_grid,
        k_ss = k_ss,
        y_ss = y_ss,
        l_ss = l_ss,
        w_ss = w_ss,
        rf_ss = rf_ss,
        rk_ss = rk_ss,
        v_ss_vec = v_ss_vec,
        mc_ss_vec = mc_ss_vec,
        c_ss_vec = c_ss_vec,
        q_l_ss_vec = q_l_ss_vec,
        tot_wealth_ss_vec = tot_wealth_ss_vec
        );
end