! -------------------------------------------------------------------------
! mod_param.f90: parameter and setup module
! initializes parameters, variables, smolyak grid and basis functions 
! public subroutines:
! init_setup: loads parameter file and prepares smolyak grid and basis functions
! grid_setup: after steady state calculation, scales state grid
! -------------------------------------------------------------------------
! authors:         Rohan Kekre and Moritz Lenel
! for updates see: https://github.com/KekreLenel/MPR
! -------------------------------------------------------------------------
module mod_param    

use base_lib        ! load helper functions

implicit none
private

public :: init_setup, grid_setup, bbeta_vec, ies_vec, gma_vec, tht, aalpha, ddelta, sigma_z,      &
          chiX, phi, tayl_ic, sig_m, rho_m, disast_p, varphi_p, rho_p, disast_std,   & 
          vareps_w, tau_w, chiW, lmbd_vec, s_trgt_vec, l_target, k_grid_adj, w_grid_adj,          & 
          low_guess_fixed, high_guess_fixed, sig_dis, n_I, thtbar_vec, jx_1y,          &
          k_grid_dev, w_grid_dev, k_grid_mean, w_grid_mean, v_normalization_vec, k_dev_param,           &
          w_dev_param, quad_weight_vec, shock_grid, wealth_share_grid, next_m_mat, next_dis_mat,        & 
          state_grid, n_quad, n_states, n_active_dims, vector_mus_dimensions, n_shocks, idx_k, idx_sa, & 
          idx_sc, idx_m, idx_w, idx_dis, smolyak_d, smol_polynom, max_smol_level, smol_elem_ani,& 
          s_a_grid_mean, s_a_grid_dev, s_c_grid_mean, s_c_grid_dev, m_grid_mean,        &
          m_grid_dev, dis_grid_mean, dis_grid_dev, dz_vec, dz_vec_adj, results_folder, irf_indices,       & 
          irf_shock_sizes, sidx_m, sidx_z, sidx_dis, xi, s_bar_vec, n_interp, n_bond, &
          smol_grid, no_shock_idx, kbar_constr, constr_agents, labor_alloc_vec, use_idio_risk, &
          idio_weight_vec, idio_shock_grid, n_idio, n_uni_quad_dis, n_uni_quad_m, n_uni_quad_g, &
          command_input, n_sample, n_sim_periods, n_burn, n_irf_periods, n_prm

! -------------------------------------------------------------------------------- !
! Declare module variables
! -------------------------------------------------------------------------------- !

! numerical parameters: 
integer, parameter  :: smolyak_d = 6, n_I = 3, n_shocks = 3,                   & ! number of states, number of agents, number of shocks
                       n_uni_quad_g = 3, n_uni_quad_m = 3, n_uni_quad_dis = 3, & ! number of quadrature nodes in each dimension
                       n_GH_points = n_uni_quad_g*n_uni_quad_m*n_uni_quad_dis, & ! total number of quadrature shocks
                       n_quad = n_GH_points + 1, max_smol_level = 3,           & ! number of shocks w. disaster shock, smolyak density
                       n_bond = 4, n_interp=70+n_bond*2, jx_1y = 74,           & ! number of priced bonds, number of result variables, 1y yield index
                       n_sample=1000, n_sim_periods = 50000 ,                  & ! sample size for IRFs, number of simulation periods
                       n_irf_periods = 200, n_burn = 5000, n_prm = 54            ! length of impulse response, burn period, number of parameters

! States: 1. Capital,         2. Wealth share a, 3. Wealth share c, 
!         4. monetary policy, 5. Wage,           6. Disaster
integer, parameter  :: idx_k=1, idx_sa=2, idx_sc=3, idx_m=4, idx_w=5, idx_dis=6

! Shocks: 1. Growth, 2. Monetary, 3. Disaster 
integer, parameter  :: sidx_z = 1, sidx_m = 2, sidx_dis = 3

integer, save       :: vector_mus_dimensions(smolyak_d), n_active_dims, n_states
character           :: command_input*2 

! IRF Related
integer, parameter :: irf_indices(3) = [1, 2, 3]
real(dp), save     :: irf_shock_sizes(3)
integer, save      :: no_shock_idx(1)

! Parameters: will be assigned values from external parameter file
integer, save      :: constr_agents(n_I), use_idio_risk, n_idio
real(dp), save     :: tht, aalpha, ddelta, sigma_z, chiX, phi, tayl_ic,         &
                      rho_m, disast_p, varphi_p, rho_p, disast_std,             &
                      vareps_w, tau_w, chiW, xi, l_target, sig_dis,             & 
                      bbeta_vec(n_I), gma_vec(n_I), ies_vec(n_I), lmbd_vec(n_I),& 
                      s_trgt_vec(n_I), s_bar_vec(n_I),kbar_constr(n_I),         & 
                      k_grid_adj, w_grid_adj, k_dev_param,  high_guess_fixed,   & 
                      w_dev_param, thtbar_vec(n_I), v_normalization_vec(n_I),   & 
                      low_guess_fixed, IRF_z, IRF_m, IRF_dis, sig_m,            & 
                      labor_alloc_vec(3), std_idio_risk_a, std_idio_risk_b,     &
                      std_idio_risk_c, gov_debt

! Grid and basis functions
real(dp), allocatable, save :: smol_grid(:,:), state_grid(:,:), smol_polynom(:,:), wealth_share_grid(:,:)
integer, allocatable, save  ::  smolyak_elem_iso(:,:), smol_elem_ani(:,:)

real(dp), save :: k_grid_dev, s_a_grid_dev, s_c_grid_dev,       & 
                  m_grid_dev, w_grid_dev, dis_grid_dev,         &
                  k_grid_mean, s_a_grid_mean, s_c_grid_mean,    & 
                  m_grid_mean, w_grid_mean, dis_grid_mean

! Quadrature shocks and corresponding exogenous state transitions
real(dp), save :: shock_grid(n_quad,n_shocks), quad_weight_vec(n_GH_points), dz_vec(n_quad), dz_vec_adj(n_quad)
real(dp), allocatable, save :: next_m_mat(:,:), next_dis_mat(:,:), quad_vec_idio(:), idio_weight_vec(:), idio_shock_grid(:, :)

! string for name of results folder
character, save :: results_folder * 100

contains

! ---------------------------------------------------------------------- !
! Initial Setup 
! ---------------------------------------------------------------------- !
subroutine init_setup
    use mod_smolyak
    use base_lib, only: get_quadrature_points 

    ! ------------------------------------------------------------ !
    ! Declare Local Variables 
    ! ------------------------------------------------------------ !
    integer            :: command_len, narg
    character          :: param_file*100, n_prm_str*2
    character(12), dimension(n_prm) :: paramnames
    real(dp)           ::  param_input_vec(n_prm), cov_mat(n_shocks, n_shocks), &
                           quad_vec_temp(n_shocks), weight_tmp
    integer, allocatable :: mus_dimensions_redux(:)
    integer              :: counter, info, ggg, mmm, ddd, iii

    !! Quadrature internal variables (one-dimensional grids and weights)
    real(dp) :: uni_weight_vec_g(n_uni_quad_g), & 
                uni_weight_vec_m(n_uni_quad_m), & 
                uni_weight_vec_dis(n_uni_quad_dis)
    real(dp) :: quad_vec_g(n_uni_quad_g), &
                quad_vec_m(n_uni_quad_m), &
                quad_vec_dis(n_uni_quad_dis)

    ! ----------------------------------------------------------------- !
    ! 1. Read in parameters
    ! ----------------------------------------------------------------- !
    
    ! Read command line input: which parameterization to run 
    narg=command_argument_count()
    if (narg>0) then
        call get_command_argument(1,command_input,command_len)
        param_file = '../src/params/param_file_' // command_input(1:command_len) // '.csv'
        results_folder = '../output/tmp/res_' // command_input(1:command_len) // '/'
    else
        write(*,*) 'ERROR: No run index specified.'
        stop
    endif

    ! Open parameter file
    open (unit = 10, file = param_file, ACTION="read",  &
            & FORM="formatted", ACCESS="sequential")

    ! Read parameters
    read(10,*) paramnames
    write(n_prm_str,'(i2)') n_prm
    read(10,'(' // n_prm_str // 'f100.0)') param_input_vec
    close(10)

    ! Assign Parameters
    lmbd_vec(1)         = param_input_vec(1)
    lmbd_vec(2)         = param_input_vec(2)
    lmbd_vec(3)         = 1.0_dp - lmbd_vec(1) - lmbd_vec(2) 
    bbeta_vec(1)        = param_input_vec(3)
    bbeta_vec(2)        = param_input_vec(4)
    bbeta_vec(3)        = param_input_vec(5)
    gma_vec(1)          = param_input_vec(6)
    gma_vec(2)          = param_input_vec(7)
    gma_vec(3)          = param_input_vec(8)
    ies_vec(1)          = param_input_vec(9)
    ies_vec(2)          = param_input_vec(10)
    ies_vec(3)          = param_input_vec(11)
    tht                 = param_input_vec(12)
    ddelta              = param_input_vec(13)
    aalpha              = param_input_vec(14)
    sigma_z             = param_input_vec(15)
    chiX                = param_input_vec(16)
    phi                 = param_input_vec(17)
    tayl_ic             = param_input_vec(18)
    sig_m               = param_input_vec(19)
    rho_m               = param_input_vec(20)
    disast_p            = param_input_vec(21)
    varphi_p            = param_input_vec(22)
    rho_p               = param_input_vec(23)
    disast_std          = param_input_vec(24)
    vareps_w            = param_input_vec(25)
    tau_w               = param_input_vec(26)
    chiW                = param_input_vec(27)
    s_bar_vec(1)        = param_input_vec(28)
    s_bar_vec(3)        = param_input_vec(29)
    s_bar_vec(2)        = 1.0_dp - s_bar_vec(1) - s_bar_vec(3)
    s_trgt_vec(1)       = param_input_vec(30)
    s_trgt_vec(3)       = param_input_vec(31)
    s_trgt_vec(2)       = 1.0_dp - s_trgt_vec(1) - s_trgt_vec(3)
    xi                  = param_input_vec(32)
    l_target            = param_input_vec(33)
    labor_alloc_vec(1)  = param_input_vec(34)
    labor_alloc_vec(2)  = param_input_vec(35)
    labor_alloc_vec(3)  = param_input_vec(36)
    kbar_constr(1)      = param_input_vec(37)
    kbar_constr(2)      = param_input_vec(37)
    kbar_constr(3)      = param_input_vec(37)
    gov_debt            = param_input_vec(38)
    k_grid_adj          = param_input_vec(39)
    w_grid_adj          = param_input_vec(40)
    s_a_grid_dev        = param_input_vec(41)
    s_c_grid_dev        = param_input_vec(42)
    k_dev_param         = param_input_vec(43)
    w_dev_param         = param_input_vec(44)
    IRF_z               = param_input_vec(45)
    IRF_m               = param_input_vec(46)
    IRF_dis             = param_input_vec(47)
    constr_agents(1)    = int(param_input_vec(48))
    constr_agents(2)    = int(param_input_vec(49))
    constr_agents(3)    = int(param_input_vec(50))
    use_idio_risk       = int(param_input_vec(51))
    std_idio_risk_a     = param_input_vec(52)
    std_idio_risk_b     = param_input_vec(53)
    std_idio_risk_c     = param_input_vec(54)

    ! Grid density
    vector_mus_dimensions = [3, 3, 3, 3, 3, 3]
    
    ! wage state variable not needed without stickiness
    if (chiW < 1E-3) then
        vector_mus_dimensions(idx_w) = 0
    endif

    ! monetary state variable needed for solving the model with monetary shock persistence 
    if (rho_m > 1E-6) then
        vector_mus_dimensions(idx_m) = 3
    endif
    
    ! disaster shock standard deviation
    sig_dis = disast_std*sqrt(1.0-rho_p**2)

    ! min/max leverage constraints for numerical stability
    low_guess_fixed  = -20.0_dp
    high_guess_fixed =   3.0_dp

    ! Print parameters to console
    call print_parameters()

    ! ----------------------------------------------------------------- !
    ! 2. Construct Quadratures
    ! Ordering: g, m, dis
    ! ----------------------------------------------------------------- !
    ! (1) Variance-Covariance Matrix and its Cholesky Decomposition
    cov_mat = 0.0_dp
    cov_mat(1,1) = sigma_z**2   
    cov_mat(2,2) = sig_m**2
    cov_mat(3,3) = sig_dis**2 

    ! Cholesky decomposition
    call F07FDF('U', n_shocks, cov_mat, 3, info)
    if (info /= 0) then 
        write(*,*) 'ERROR IN FACTOR DECOMPOSITION', info
        stop
    endif

    ! (2) Get One-dimensional Quadrature Nodes & Weights 
    ! TFP
    if (n_uni_quad_g > 1) then
        call get_quadrature_points(0.0_dp, 1.0_dp, n_uni_quad_g, quad_vec_g, uni_weight_vec_g)
    else
        quad_vec_g = 0.0_dp
        uni_weight_vec_g = 1.0_dp
    endif

    ! Monetary Policy
    if (n_uni_quad_m > 1) then
        call get_quadrature_points(0.0_dp, 1.0_dp, n_uni_quad_m, quad_vec_m, uni_weight_vec_m)
    else
        quad_vec_m       = 0.0_dp
        uni_weight_vec_m = 1.0_dp
    endif

    ! Disaster
    if (n_uni_quad_dis > 1) then
        call get_quadrature_points(0.0_dp, 1.0_dp, n_uni_quad_dis, quad_vec_dis, uni_weight_vec_dis)
    else
        quad_vec_dis       = 0.0_dp
        uni_weight_vec_dis = 1.0_dp
    endif

    ! (3) Combine 
    counter = 0
    !DIR$ NOUNROLL
    do ggg = 1, n_uni_quad_g
    do mmm = 1, n_uni_quad_m
    do ddd = 1, n_uni_quad_dis

        counter = counter + 1

        !! Combine 1D nodes
        quad_vec_temp  = [ quad_vec_g(ggg), quad_vec_m(mmm), quad_vec_dis(ddd)]
        
        !! Scale by square root of variance-covariance matrix
        do iii = 1,n_shocks
            shock_grid(counter, iii) = sum(cov_mat(iii,:)*quad_vec_temp) 
        enddo

        !! Multiply 1D weights
        quad_weight_vec(counter) =  uni_weight_vec_g(ggg)*uni_weight_vec_m(mmm)*  & 
                                    uni_weight_vec_dis(ddd)

    end do
    end do
    end do

    ! (4) Extend the shock grid to include disaster
    shock_grid(n_quad, 1) = -varphi_p
    shock_grid(n_quad, 2:n_shocks) = 0.0_dp

    ! Growth vector dz_vec 
    dz_vec             = shock_grid(:,1) 

    ! adjusted growth vector: dz_vec without disaster shock
    ! needed for rescaling wage and capital for stationary solution 
    dz_vec_adj         = dz_vec
    dz_vec_adj(n_quad) = 0.0_dp 

    ! find shock index with no shock for stochastic steady state calculations
    no_shock_idx =  minloc( abs(shock_grid(:,1)) + &
                            abs(shock_grid(:,2)) + &
                            abs(shock_grid(:,3)))
    
    ! impulse response shock sizes
    irf_shock_sizes = [IRF_z*sigma_z, IRF_m*sig_m, IRF_dis*sig_dis]
    
    ! ----------------------------------------------------------------- !
    ! 2A. Construct idiosyncratic risk quadrature
    ! ----------------------------------------------------------------- !

    ! Quadrature Information
    if (use_idio_risk == 1) then
        n_idio = 2
        allocate(quad_vec_idio(n_idio), idio_weight_vec(n_idio), idio_shock_grid(n_idio, n_I))
        call get_quadrature_points(0.0_dp, 1.0_dp, n_idio, quad_vec_idio, idio_weight_vec)
    else
        n_idio = 1
        allocate(quad_vec_idio(n_idio), idio_weight_vec(n_idio), idio_shock_grid(n_idio, n_I))
        quad_vec_idio   = 0.0_dp
        idio_weight_vec = 1.0_dp
    endif
    
    ! Map to individual - adjust for jensen term
    idio_shock_grid(:,1) = std_idio_risk_a * quad_vec_idio - log(sum(idio_weight_vec*exp(quad_vec_idio*std_idio_risk_a)))
    idio_shock_grid(:,2) = std_idio_risk_b * quad_vec_idio - log(sum(idio_weight_vec*exp(quad_vec_idio*std_idio_risk_b)))
    idio_shock_grid(:,3) = std_idio_risk_c * quad_vec_idio - log(sum(idio_weight_vec*exp(quad_vec_idio*std_idio_risk_c)))

    ! ----------------------------------------------------------------- !
    ! 3. Create Smolyak Grid
    ! ----------------------------------------------------------------- !   
    ! (1) Get Active Dimensions 
    n_active_dims = count(vector_mus_dimensions > 0)
    allocate(mus_dimensions_redux(n_active_dims))
    counter = 0
    do ddd = 1,smolyak_d 
    if (vector_mus_dimensions(ddd) > 0) then
        counter = counter + 1
        mus_dimensions_redux(counter) = vector_mus_dimensions(ddd)
    endif
    enddo
    write(*,"(A19,I2,A1,I2)") ' Active dimensions:',n_active_dims,'/',smolyak_d

    ! (2) Create Smolyak Grid
    smolyak_elem_iso = Smolyak_Elem_Isotrop(n_active_dims, maxval(mus_dimensions_redux))
    smol_elem_ani    = Smolyak_Elem_Anisotrop(smolyak_elem_iso, n_active_dims, mus_dimensions_redux)
    smol_grid        = Smolyak_Grid(n_active_dims,maxval(mus_dimensions_redux), Smol_elem_ani)

    smol_polynom     = Smolyak_Polynomial(smol_grid,n_active_dims,maxval(mus_dimensions_redux), smol_elem_ani); 
    n_states         = size(smol_grid,1)
    write(*,"(A19,I5)") ' Number of states:  ', n_states

    ! (3) Initialize State Grid and next period matrices
    allocate(wealth_share_grid(n_I,n_states), state_grid(smolyak_d, n_states))
    allocate(next_m_mat(n_quad, n_states), next_dis_mat(n_quad, n_states))

    ! (4) Set mean and standard deviations for most states 

    !! grid width for exogenous states 
    m_grid_dev     = 2.5*sig_m/sqrt(1.0-rho_m**2)
    dis_grid_dev   = 2.5*disast_std

    !! grid means
    s_a_grid_mean = s_trgt_vec(1)
    s_c_grid_mean = s_trgt_vec(3)
    m_grid_mean       = 0.0
    dis_grid_mean     = disast_p

end subroutine init_setup

! ---------------------------------------------------------------------- !
! Grid Setup 
! ---------------------------------------------------------------------- !
subroutine grid_setup()

    integer :: counter, sss 

    ! ----------------------------------------------------------------- !
    ! Create state grid from Smolyak grid 
    ! ----------------------------------------------------------------- !
    counter = 0
    
    if (vector_mus_dimensions(1)>0) then 
    counter = counter + 1
    state_grid(1,:) = smol_grid(:,counter)*k_grid_dev + k_grid_mean
    else 
    state_grid(1,:) = k_grid_mean
    endif

    if (vector_mus_dimensions(2)>0) then 
    counter = counter + 1
    state_grid(2,:) = smol_grid(:,counter)*s_a_grid_dev + s_a_grid_mean
    else 
    state_grid(2,:) = s_a_grid_mean
    endif

    if (vector_mus_dimensions(3)>0) then 
    counter = counter + 1
    state_grid(3,:) = smol_grid(:,counter)*s_c_grid_dev + s_c_grid_mean
    else 
    state_grid(3,:) = s_c_grid_mean
    endif

    if (vector_mus_dimensions(4)>0) then 
    counter = counter + 1
    state_grid(4,:) = smol_grid(:,counter)*m_grid_dev + m_grid_mean
    else 
    state_grid(4,:) = m_grid_mean
    endif

    if (vector_mus_dimensions(5)>0) then 
    counter = counter + 1
    state_grid(5,:) = smol_grid(:,counter)*w_grid_dev + w_grid_mean
    else 
    state_grid(5,:) = w_grid_mean
    endif

    if (vector_mus_dimensions(6)>0) then 
    counter = counter + 1
    state_grid(6,:) = smol_grid(:,counter)*dis_grid_dev + dis_grid_mean
    else 
    state_grid(6,:) = dis_grid_mean
    endif


    ! ----------------------------------------------------------------- !
    ! Derive wealth share grids from state grids
    ! ----------------------------------------------------------------- !
    wealth_share_grid(1,:) = state_grid(idx_sa,:) / lmbd_vec(1)
    wealth_share_grid(2,:) = (1.0-state_grid(idx_sa,:)-state_grid(idx_sc,:)) / lmbd_vec(2)
    wealth_share_grid(3,:) = state_grid(idx_sc,:) / lmbd_vec(3)

    ! ----------------------------------------------------------------- !
    ! Next period position matrix 
    ! Can be constructed directly for the 2 exogenous states 
    ! ----------------------------------------------------------------- !
    do sss = 1,n_states
        next_m_mat(:,sss)   = (1.0 - rho_m)*0.0_dp        + rho_m*state_grid(idx_m,sss) &
                               + shock_grid(:,sidx_m)
        next_dis_mat(:,sss)  = (1.0 - rho_p)*disast_p + rho_p*state_grid(idx_dis,sss) & 
                               + shock_grid(:,sidx_dis)
    enddo


end subroutine grid_setup


subroutine print_parameters()
    ! prints subset of parameters to console 
    write(*,*)
    write(*,*) '---------------------------------------------'
    write(*,*) 'PARAMETRIZATION'
    write(*,*) '---------------------------------------------'
    write(*,*)
    write(*,"(A21, F10.4, F10.4, F10.4, F10.4)") ' bbeta        = ' , bbeta_vec
    write(*,"(A21, F10.4, F10.4, F10.4, F10.4)") ' gma          = ' , gma_vec
    write(*,"(A21, F10.4, F10.4, F10.4, F10.4)") ' ies          = ' , ies_vec
    write(*,"(A21, F10.4, F10.4, F10.4, F10.4)") ' s_bar        = ' , s_bar_vec  
    write(*,"(A21, F10.4, F10.4, F10.4, F10.4)") ' s_target     = ' , s_trgt_vec
    write(*,"(A21, F10.4, F10.4, F10.4, F10.4)") ' lmbd_vec     = ' , lmbd_vec
    write(*,"(A21, F10.4)")                      ' tht          = ' , tht
    write(*,"(A21, F10.4)")                      ' ddelta       = ' , ddelta
    write(*,"(A21, F10.4)")                      ' aalpha       = ' , aalpha
    write(*,"(A21, F10.4)")                      ' sig_g        = ' , sigma_z  
    write(*,"(A21, F10.4)")                      ' phi          = ' , phi  
    write(*,"(A21, F10.4)")                      ' tayl_ic      = ' , tayl_ic
    write(*,"(A21, F10.4)")                      ' sig_m        = ' , sig_m  
    write(*,"(A21, F10.4)")                      ' disast_p     = ' , disast_p  
    write(*,"(A21, F10.4)")                      ' varphi_p     = ' , varphi_p  
    write(*,"(A21, F10.4)")                      ' rho_p        = ' , rho_p  
    write(*,"(A21, F10.4)")                      ' disast_std   = ' , disast_std  
    write(*,"(A21, F10.4)")                      ' vareps_w     = ' , vareps_w  
    write(*,"(A21, F10.4)")                      ' tau_w        = ' , tau_w  
    write(*,"(A21, F10.4)")                      ' chiW         = ' , chiW  
    write(*,"(A21, F10.4)")                      ' chiX         = ' , chiX    
    write(*,"(A21, F10.4)")                      ' xi           = ' , xi           
    write(*,"(A21, F10.4)")                      ' l_target     = ' , l_target     
    write(*,*)
    write(*,*) '---------------------------------------------'
    write(*,*)
end subroutine print_parameters

end module mod_param












