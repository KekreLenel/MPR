function get_dis(;disast_p= 0.005,disast_std= 0.0025,i=1)
    function xfvalflag!(F,x)
        F[1]=1000*(disast_p - exp(x[1] + x[2]/2))
        F[2]=1000*(disast_std^2 - (exp(x[2])-1)*exp(2*x[1] + x[2]))
    end
    X = nlsolve(xfvalflag!,[-1.0; 1.0],autodiff =:forward);
    return X.zero[i];
end

@with_kw struct prm{F<:Int64,R<:Real}
    lmbd_a          ::R = 0.04 # pop share a
    lmbd_b          ::R = 0.36 # pop share b
    lmbd_c          ::R = 0.6
    debt_to_equity  ::R = 0.5

    # Preference Parameters
    bbeta_a         ::R = 0.98
    bbeta_b         ::R = 0.98
    bbeta_c         ::R = 0.98
    gma_a           ::R = 10.0 
    gma_b           ::R = 25.5
    gma_c           ::R = round((lmbd_a/(lmbd_a + lmbd_b)/gma_a + lmbd_b/(lmbd_a + lmbd_b)/gma_b)^(-1),digits=6)
    ies_a           ::R = 0.8
    ies_b           ::R = 0.8
    ies_c           ::R = 0.8
    theta           ::R = 1.0  # frisch     

    # Technology Parameters
    delta           ::R = 0.025
    aalpha          ::R = 0.33
    sig_z           ::R = 0.0055   
    chiX            ::R = 3.5      # capital adjustment cost

    # Monetary Policy
    phi             ::R = 1.5      # Taylor rule coefficient
    tayl_ic         ::R = 0.003    # Taylor rule intercept: approx. 0# inflation 
    sig_m           ::R = 0.0025/4 # Monetary policy shock
    rho_m           ::R = 0.0

    # Disaster Shock
    disast_p        ::R = round(get_dis(i = 1),digits=6)       # average log(p)
    varphi_p        ::R = 0.15       # disaster depth
    rho_p           ::R = 0.8        # disaster autocorrelation
    disast_std      ::R = round(sqrt(get_dis(i = 2)),digits=6) # unconditional std. log(p)

    # Rotemberg wage adjustment
    vareps_w        ::R = 10.0
    tau_w           ::R = round(1/(1-vareps_w),digits=6)
    chiW            ::R = 150.0

    # Newborn wealth allocation 
    s_bar_a         ::R = 0.0
    s_bar_c         ::R = -0.25 

    # Target wealth shares (a and c wealth are state variables)
    s_trgt_a        ::R = 0.18
    s_trgt_c        ::R = 0.23

    # Mortality rate
    xi              ::R = 0.01

    # Average labor supply 
    l_target        ::R = 1.0

    # Labor share: phi_i
    labor_alloc_a   ::R = round(0.03/lmbd_a ,digits=6)
    labor_alloc_b   ::R = round(0.14/lmbd_b ,digits=6)
    labor_alloc_c   ::R = round(0.83/(1.0 - lmbd_a - lmbd_b) ,digits=6)

    # lower bound capital holdings
    kbar            ::R =  10.0

    # Government debt 
    gov_debt        ::R = 0.10

    # Grid scaling parameters: width and center relative to steady state
    k_grid_adj      ::R = 1.25
    w_grid_adj      ::R = 1.06
    s_a_dev     ::R = 0.1
    s_c_dev     ::R = 0.1
    k_dev_param     ::R = 0.1
    w_dev_param     ::R = 0.05 

    # IRF shock sizes 
    IRF_g           ::R =  2.0
    IRF_m           ::R = -10.0
    IRF_dis         ::R =  2.0

    # for efficiency kbar constrained only activated for group c, not binding for a anddkj b 
    constrained_a   ::F = 0
    constrained_b   ::F = 0
    constrained_c   ::F = 1

    # idiosyncratic risk
    use_idio_risk   ::F = 0
    idio_risk_a     ::R = 0.00
    idio_risk_b     ::R = 0.00
    idio_risk_c     ::R = 0.00

end


@with_kw struct mpr_economy{F<:Integer,R<:Real}
    smolyak_d  ::F =  6 ; 
    n_I  ::F =  3;  
    n_shocks  ::F =  3;                     # number of states  number of agents  number of shocks
    n_uni_quad_g  ::F =  3;  
    n_uni_quad_m  ::F =  3;  
    n_uni_quad_dis  ::F =  3;   # number of quadrature nodes in each dimension
    n_GH_points  ::F =  n_uni_quad_g*n_uni_quad_m*n_uni_quad_dis;   # total number of quadrature shocks
    n_quad  ::F =  n_GH_points + 1 ; 
    max_smol_level  ::F =  3 ;            # number of shocks w. disaster shock  smolyak density
    n_bond  ::F =  4  ;
    n_interp ::F = 70+n_bond*2   ;                      # number of priced bonds  number of result variables
    n_sample ::F = 1000  ;  # 1000/100
    n_sim_periods  ::F =  50000 ;      # 50000/1000              # sample size for IRFs  number of simulation periods
    n_irf_periods  ::F =  200  ;     # 200/40
    n_burn  ::F = 5000 ;   # 5000/100
    n_prm  ::F =  54;
    jx_1y ::F = 74;

    # States: 1. Capital          2. Wealth share a  3. Wealth share c  
    #         4. monetary policy  5. Wage            6. Disaster
    idx_k ::F = 1;  idx_sa ::F = 2;  idx_sc ::F = 3;  
    idx_m ::F = 4; idx_w ::F = 5;  idx_dis ::F = 6;

    # Shocks: 1. Growth  2. Monetary  3. Disaster 
    sidx_z  ::F =  1;  sidx_m  ::F =  2;   sidx_dis  ::F =  3;

    vector_mus_dimensions ::Array{F,1};
    n_active_dims ::F;
    n_states ::F;

    # IRF Related
    irf_indices ::Array{F,1};
    irf_shock_sizes ::Array{R,1};
    no_shock_idx ::F = 14 ;

    # Parameters: will be assigned values from external parameter file
    constr_agents ::Array{F,1};
    use_idio_risk  ::F;
    n_idio::F = 1;

    tht::R; aalpha::R;  ddelta::R;  sigma_z::R;  chiX::R;  phi::R;  tayl_ic::R ;        
    rho_m::R;  disast_p::R;  varphi_p::R; rho_p::R;  disast_std::R;      
    vareps_w::R;  tau_w::R;  chiW::R;  xi::R;  l_target::R;  sig_dis::R;  
    bbeta_vec::Array{R,1};  gma_vec::Array{R,1};
    ies_vec::Array{R,1};  lmbd_vec::Array{R,1};            
    s_trgt_vec::Array{R,1};  s_bar_vec::Array{R,1};
    kbar_constr::Array{R,1};          
    k_grid_adj::R;  w_grid_adj::R;  k_dev_param::R;               
    w_dev_param::R;  thtbar_vec::Array{R,1};  v_normalization_vec::Array{R,1};        
    low_guess_fixed::R;  IRF_z::R;  IRF_m::R;  IRF_dis::R;  sig_m::R;  high_guess_fixed::R;   
    labor_alloc_vec::Array{R,1} = zeros(3);  std_idio_risk_a::R;            
    std_idio_risk_b::R;  std_idio_risk_c::R;  gov_debt::R;

    k_ss::R;y_ss::R;l_ss::R;w_ss::R;rf_ss::R;rk_ss::R;

    v_ss_vec::Array{R,1}; c_ss_vec::Array{R,1}; mc_ss_vec::Array{R,1}; q_l_ss_vec::Array{R,1};
    tot_wealth_ss_vec::Array{R,1};


    # Grid and basis functions
    smol_grid::Array{R,2};  state_grid::Array{R,2}; smol_polynom::Array{R,2};  wealth_share_grid::Array{R,2};
    smolyak_elem_iso::Array{R,2}; smol_elem_ani::Array{F,2};

    k_grid_dev::R;  s_a_grid_dev::R;  s_c_grid_dev::R;         
                  m_grid_dev::R;  w_grid_dev::R;        dis_grid_dev::R;            
                  k_grid_mean::R;  s_a_grid_mean::R;  s_c_grid_mean::R;      
                  m_grid_mean::R;  w_grid_mean::R;  dis_grid_mean::R;

    # Quadrature shocks and corresponding exogenous state transitions
    shock_grid::Array{R,2};  quad_weight_vec::Array{R,1};  
    dz_vec::Array{R,1};  dz_vec_adj::Array{R,1};
    next_m_mat::Array{R,2};  next_dis_mat::Array{R,2};  
    quad_vec_idio::Array{R,1};  idio_weight_vec::Array{R,1};  idio_shock_grid::Array{R,2};
end


# ======================================================

# number of parameters
n_params = size(fieldnames(prm),1);

# how many alternative calibrations
n_comp              = 9;

# how many variables varied
n_comp_prms = 33;

prm_list = ["rho_m", "chiX", "chiW", "use_idio_risk", "idio_risk_a","idio_risk_b", 
    "idio_risk_c", "bbeta_a", "bbeta_b", "bbeta_c", "gma_a", "gma_b", "gma_c", 
    "constrained_c", "xi", "s_bar_a", "s_bar_c", "lmbd_a", "lmbd_b", 
    "s_trgt_a", "s_trgt_c", "labor_alloc_a", "labor_alloc_b", 
    "labor_alloc_c", "k_grid_adj", "w_grid_adj", "s_a_dev", "s_c_dev",  
    "tayl_ic", "disast_p", "varphi_p", "rho_p", "disast_std"];

prm_vec = Array{prm}(undef, n_comp);

prm_vec[1] = prm{Int64,Float64}();
gamma = (prm_vec[1].s_trgt_a/prm_vec[1].gma_a  + (1 - prm_vec[1].s_trgt_a - prm_vec[1].s_trgt_c)/prm_vec[1].gma_b+ prm_vec[1].s_trgt_c/prm_vec[1].gma_c)^(-1)
prm_vec[2] = prm{Int64,Float64}(gma_a = gamma, gma_b = gamma, gma_c = gamma, constrained_c = 0.0, s_bar_a = 1/3,
        s_bar_c = 1/3, lmbd_a = 1/3, lmbd_b = 1/3, lmbd_c = 1/3, s_trgt_a = 1/3, s_trgt_c = 1/3, labor_alloc_a = 1.0,
        labor_alloc_b = 1.0, labor_alloc_c = 1.0, k_grid_adj = 1.05, w_grid_adj = 1.02, tayl_ic = 0.005);
prm_vec[3] = prm{Int64,Float64}(rho_m = 0.75);
prm_vec[4] = prm{Int64,Float64}(chiX = 0.0);
prm_vec[5] = prm{Int64,Float64}(chiW = 0.0);
prm_vec[6] = prm{Int64,Float64}(use_idio_risk = 1.0, idio_risk_a = 0.0, idio_risk_b = 0.036,
        idio_risk_c = sqrt(prm_vec[1].lmbd_a/(prm_vec[1].lmbd_a + prm_vec[1].lmbd_b)*0.0 + prm_vec[1].lmbd_b/(prm_vec[1].lmbd_a + prm_vec[1].lmbd_b)*(0.036^2)),
        bbeta_a = 0.98, bbeta_b = 0.98, bbeta_c = 0.98, gma_a = 11.0, gma_b = 11.0, gma_c = 11.0, xi = 0.01,
        s_bar_a = -0.02, s_bar_c = -0.15, k_grid_adj = 1.20, w_grid_adj = 1.05, tayl_ic = 0.0035,
        disast_p = get_dis(disast_p = 0.005, disast_std = 0.009, i = 1),
        disast_std = sqrt(get_dis(disast_p = 0.005, disast_std = 0.009, i = 2)));
risk = (prm_vec[1].s_trgt_a*prm_vec[6].idio_risk_a +  (1-prm_vec[1].s_trgt_a - prm_vec[1].s_trgt_c)*prm_vec[6].idio_risk_b+ prm_vec[1].s_trgt_c*prm_vec[6].idio_risk_c);          
risk = round(risk,digits=6)
prm_vec[7] = prm{Int64,Float64}(use_idio_risk = 1.0, idio_risk_a = risk, idio_risk_b = risk, idio_risk_c = risk,
        bbeta_a = prm_vec[6].bbeta_a, bbeta_b = prm_vec[6].bbeta_a, bbeta_c = prm_vec[6].bbeta_a,
        gma_a = prm_vec[6].gma_a, gma_b = prm_vec[6].gma_a, gma_c = prm_vec[6].gma_a,
        constrained_c = 0.0, xi = prm_vec[6].xi, s_bar_a = 1/3,
        s_bar_c = 1/3, lmbd_a = 1/3, lmbd_b = 1/3,lmbd_c = 1/3, s_trgt_a = 1/3, s_trgt_c = 1/3, labor_alloc_a = 1.0,
        labor_alloc_b = 1.0, labor_alloc_c = 1.0, k_grid_adj = 1.05, w_grid_adj = 1.02, tayl_ic = 0.006,
        disast_p = prm_vec[6].disast_p, disast_std = prm_vec[6].disast_std);
prm_vec[8] = prm{Int64,Float64}(bbeta_a = 0.984, bbeta_b = 0.984, bbeta_c = 0.984,
        gma_a = 2.5, gma_b = 21.0, gma_c = (0.004/0.4/2.5 + 0.396/0.4/21)^(-1),
        xi = 0.04, s_bar_a = 0.0,
        s_bar_c = 0.009/0.04, lmbd_a = 0.004, lmbd_b = 0.396, s_trgt_a = 0.02, s_trgt_c = 0.23, labor_alloc_a = 0.0,
        labor_alloc_b = 0.17/0.396, labor_alloc_c = 0.83/0.6, k_grid_adj = 1.08, w_grid_adj = 1.02, 
        s_a_dev = 0.019,tayl_ic = 0.003);
gamma = (prm_vec[8].s_trgt_a/prm_vec[8].gma_a + (1-prm_vec[8].s_trgt_a-prm_vec[8].s_trgt_c)/prm_vec[8].gma_b+prm_vec[8].s_trgt_c/prm_vec[8].gma_c)^(-1);
gamma = round(gamma,digits=6)
prm_vec[9] = prm{Int64,Float64}(bbeta_a = 0.984, bbeta_b = 0.984, bbeta_c = 0.984,
        gma_a = gamma , gma_b = gamma , gma_c = gamma,
        constrained_c = 0.0, xi =  prm_vec[8].xi, s_bar_a = 1/3,
        s_bar_c = 1/3, lmbd_a = 1/3, lmbd_b = 1/3,lmbd_c = 1/3, s_trgt_a = 1/3, s_trgt_c = 1/3, labor_alloc_a = 1.0,
        labor_alloc_b = 1.0, labor_alloc_c = 1.0, k_grid_adj = 1.05, w_grid_adj = 1.02, tayl_ic = 0.006);
println("ALL PARAMETERS GENERATED \n")





