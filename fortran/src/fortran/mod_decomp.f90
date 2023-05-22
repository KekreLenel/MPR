! -------------------------------------------------------------------------
! mod_decomp.f90: decomposition module: 
! creates inputs for tables VIII and IX in main text 
! -------------------------------------------------------------------------
! authors:         Rohan Kekre and Moritz Lenel
! for updates see: https://github.com/KekreLenel/MPR
! -------------------------------------------------------------------------
module mod_decomp

use omp_lib
use base_lib, only: dp, sqrt_eps, Fill_linspace_dp

implicit none
private

public :: monetary_decomposition

contains

subroutine monetary_decomposition()

    ! ============================================================ !
    ! Step 1: Basic Loading and Processing
    ! ============================================================ !
    use mod_smolyak, only : Smolyak_Polynomial
    use mod_param, only : n_states, &
                         ddelta, state_grid, n_quad, aalpha, smolyak_d, n_active_dims, &
                         next_m_mat, next_dis_mat, smol_polynom, max_smol_level, smol_elem_ani, &
                         n_shocks, dz_vec, dz_vec_adj, &
                         bbeta_vec, wealth_share_grid, n_sample, lmbd_vec, kbar_constr, &
                         vector_mus_dimensions, irf_indices, irf_shock_sizes, &
                         n_interp, n_bond, n_active_dims, n_I, &
                         smol_grid, no_shock_idx, shock_grid, quad_weight_vec, lmbd_vec


    character :: command_input*2, results_file*100, results_folder*100, filename*100
    integer :: narg, command_len, sss, fff, sample_idx, ppp, ddd, counter
    integer  ::  nrhs, nrhs2, nrhs3, lda, ldaf, ldb, ldx
    integer  ::  ipiv(n_states), iwork(n_states)
    integer  ::  info

    integer, parameter  :: n_nxt = 7 + 2*n_I

    INTEGER, parameter :: GENID = 3, SUBID = 1, LSEED = 1
    INTEGER            :: SEED(LSEED), STATE(633), LSTATE, IFAIL
    integer, parameter :: n_spread = 20
    real(dp) :: state_spread_vec(n_spread) 
    real(dp), allocatable :: state_series(:,:), other_vars_series(:,:), shock_series(:,:), simshock_series(:), &
                             state_spread_series(:,:,:)

    real(dp), allocatable :: b(:,:), interp_coeffs(:,:), & 
                             state_coeffs(:,:), smol_coeffs(:,:,:), & 
                             ferr(:), berr(:)
    real(dp), allocatable :: b2(:,:), ferr2(:), berr2(:), b3(:,:), ferr3(:), berr3(:)

    real(dp) ::  a(n_states,n_states), af(n_states,n_states), r(n_states), & 
             c(n_states),  work(4*n_states), polyn_points(1,n_states)

    real(dp) :: v_mat(n_states, n_I), k_next_mat(n_quad, n_states),  & 
                tht_next_mat(n_quad, n_states), q_mat(n_states), &
                l_aggr_mat(n_states), &
                next_state_mat(n_quad, n_active_dims, n_states), &
                next_state_monetary_mat(n_quad, n_active_dims, n_states), & 
                next_state_montransition_mat(1, n_active_dims, n_states)

    real(dp) :: start_vec(n_active_dims), diff
    integer  :: qqq, n_periods, n_burn, ttt, aaa
    real(dp) :: interp_input_mat(n_states, n_interp), interp_input_monetary_mat(n_states, n_interp), & 
                state_variables(smolyak_d), other_variables(n_interp), &
                state_variables_star(smolyak_d), other_variables_star(n_interp), &
                state_series_0_mat(n_sample, smolyak_d), other_vars_series_0_mat(n_sample, n_interp), &
                state_series_star_mat(n_sample, smolyak_d), other_vars_series_star_mat(n_sample, n_interp), &
                state_series_1_mat(n_sample, smolyak_d), other_vars_series_1_mat(n_sample, n_interp), starting_states(n_sample, n_active_dims)

    real(dp) :: wealth_a, wealth_b, wealth_c, q_use, i_use

    real(dp) ::  rcond
    character(1) :: equed
    character :: shock_char*1

    real(dp) :: current_state(1,n_active_dims), next_state(1,n_active_dims), stochsteady_state(1,n_active_dims), checker , &
                big_weight_vec(n_quad), next_state_lo(1,n_active_dims), next_state_hi(1,n_active_dims), &
                S_0(1,n_active_dims), S_Star(1,n_active_dims), S_Trans(1,n_active_dims), S_1(1,n_active_dims)
    integer  :: q_nxt, iii

    real(dp), allocatable :: interp_monetary_coeffs(:,:), & 
                         smol_monetary_coeffs(:,:,:), & 
                         smol_montransition_coeffs(:,:,:)

    real(dp) :: smol_coeffs_use(n_states,n_nxt), nxt_mat(n_quad, n_nxt), ferr4(n_nxt), berr4(n_nxt), & 
                nxt_mat_2(n_quad), k_next, interp_use(n_states, n_nxt), b4(n_states, n_nxt)

    real(dp) :: c_vec(n_I), share_vec(n_I), k_aggr, l_aggr, q_current, infl, nom_i, p_dis
    real(dp) :: results_a(4), results_b(4), results_c(4),cf_policy_matB(n_sample,n_I, 4), cf_policy_matA(n_sample,n_I, 4), cf_policy_mat(n_I, 4), irf_factor

    lda  = n_states; ldaf = n_states 
    ldb  = n_states; ldx  = n_states                    

    ! ------------------------------------------------------------ !
    ! Read Solution Files
    ! ------------------------------------------------------------ !
    narg=command_argument_count()
    if (narg>0) then
        call get_command_argument(1,command_input,command_len)
        results_file =   '../output/tmp/res_' // command_input(1:command_len) // '/results.txt'
        results_folder = '../output/tmp/res_' // command_input(1:command_len) // '/'
    else
        write(*,*) 'ERROR: No run index specified.'
        stop
    endif

    write(*,*)
    write(*,*) 'READING RESULTS'

    open (unit = 10, file = trim(results_folder) // 'next_state_mat.dat', ACTION="read", STATUS="replace", &
            & FORM="unformatted", ACCESS="STREAM")
    read(10) next_state_mat; close(10)

    open (unit = 10, file = trim(results_folder) // 'results_mat.dat', ACTION="read", STATUS="replace", &
            & FORM="unformatted", ACCESS="STREAM")
    read(10) interp_input_mat; close(10)

    ! ------------------------------------------------------------ !
    ! Find Stochastic Steady State
    ! ------------------------------------------------------------ !
    write(*,*)
    write(*,*) 'FIND STOCHASTIC STEADY STATE'

    ! (1.1) Polynomial coefficients for transition matrix
    nrhs = n_active_dims
    allocate(b(n_states,nrhs), smol_montransition_coeffs(1,n_states,nrhs), &
            smol_monetary_coeffs(n_quad,n_states,nrhs), smol_coeffs(n_quad,n_states,nrhs), ferr(nrhs), &  
             berr(nrhs))
    do qqq = 1,n_quad
        b = transpose(next_state_mat(qqq,:,:))
        call F07ABF( 'E','N',n_states, nrhs, smol_polynom, lda,af,ldaf,ipiv,  &
                     equed, r, c, b, ldb, smol_coeffs(qqq,:,:), & 
                     ldx, rcond, ferr, berr, work, iwork, info)
    enddo

    ! (1.2) Polynomial coefficients for non-state variables
    nrhs2 = n_interp
    allocate(b2(n_states,nrhs2), interp_coeffs(n_states,nrhs2), interp_monetary_coeffs(n_states,nrhs2), ferr2(nrhs2), berr2(nrhs2))
    b2 = interp_input_mat
    call F07ABF( 'E','N',n_states, nrhs2, smol_polynom, lda,af,ldaf,ipiv,  &
             equed, r, c, b2, ldb, interp_coeffs, & 
             ldx, rcond, ferr2, berr2, work, iwork, info )

    ! (1.3) Polynomial coefficients for state variables
    nrhs3 = smolyak_d
    allocate(b3(n_states,nrhs3), state_coeffs(n_states,nrhs3), ferr3(nrhs3), berr3(nrhs3))
    b3 = transpose(state_grid)
    call F07ABF( 'E','N',n_states, nrhs3, smol_polynom, lda,af,ldaf,ipiv,      &
                 equed, r, c, b3, ldb, state_coeffs, & 
                 ldx, rcond, ferr3, berr3, work, iwork, info )

    ! (2.1) Start from some interior grid point and move forward
    ! Use no_shock_idx: the quadrature point where there's no shock happening
    ! current_state(1,:) = smol_grid(1,:)
    current_state(1,:) = 0
    diff = 1.0_dp 

    do while (diff > sqrt_eps) 
        polyn_points = Smolyak_Polynomial(current_state, n_active_dims, &
                            max_smol_level, smol_elem_ani) 
        
        CALL DGEMM('N','N', 1, nrhs, n_states, 1.0_dp, polyn_points, 1, & 
                           smol_coeffs(no_shock_idx,:,:), n_states, 0.0_dp, & 
                           next_state, 1)   
        
        where (next_state > 1.0_dp)
            next_state = 1.0_dp 
        elsewhere (next_state < -1.0_dp)
            next_state = -1.0_dp 
        endwhere

        ! write(*,*) next_state(1,1)

        diff = maxval(abs(current_state - next_state))
        current_state = next_state
    enddo
    stochsteady_state = current_state


    ! ============================================================ !
    ! Step 1: Find S* and n^i(S*)
    ! For now, just start from the SSS, so S* is S0
    ! But write the code for any transition
    ! ============================================================ !
    open (unit = 10, file = trim(results_folder) // 'starting_states.dat', ACTION="read", STATUS="replace", &
    & FORM="unformatted", ACCESS="STREAM")
    read(10) starting_states; close(10)

    ! Read Monetary Transition Files
    open (unit = 10, file = trim(results_folder) // 'next_state_shock_mat_2.dat', ACTION="read", STATUS="replace", &
                  & FORM="unformatted", ACCESS="STREAM")
          read(10) next_state_monetary_mat; close(10)

    open (unit = 10, file = trim(results_folder) // 'next_state_shocktrans_mat_2.dat', ACTION="read", STATUS="replace", &
          & FORM="unformatted", ACCESS="STREAM")
    read(10) next_state_montransition_mat; close(10)

    open (unit = 10, file = trim(results_folder) // 'results_shock_mat_2.dat', ACTION="read", STATUS="replace", &
          & FORM="unformatted", ACCESS="STREAM")
    read(10) interp_input_monetary_mat; close(10)
    
    ! Start from S0 (for now stochastic steady state)
    cf_policy_matA = 0.0_dp 
    cf_policy_matB = 0.0_dp 
    do iii = 1,n_sample
            
        if (mod(iii, 50) == 0) then
            write(*,*) iii
        endif
    
    S_0(1,:) = starting_states(iii,:)
    polyn_points = Smolyak_Polynomial(S_0, n_active_dims, max_smol_level, smol_elem_ani) 


    ! Policies (calc_bond_price format) and states at state S_0
    CALL DGEMM('N','N', 1, nrhs2, n_states, 1.0_dp, polyn_points, 1, & 
                           interp_coeffs, n_states, 0.0_dp, & 
                           other_variables, 1)   
            ! interpolation of MPKs of constrained agent does not work well due to sparse grid
            ! adjust manually
            if (other_variables(67) > 0.5) then
                other_variables(36) = 0.0_dp
                other_variables(45) = 0.0_dp
                other_variables(54) = 0.0_dp
                other_variables(57) = lmbd_vec(3)*kbar_constr(3)
            endif 
    
    CALL DGEMM('N','N', 1, nrhs3, n_states, 1.0_dp, polyn_points, 1, & 
                       state_coeffs, n_states, 0.0_dp, & state_variables, 1)   

    state_series_0_mat(iii,:)      = state_variables
    other_vars_series_0_mat(iii,:) = other_variables



    ! Transition one period without shock 
    CALL DGEMM('N','N', 1, nrhs, n_states, 1.0_dp, polyn_points, 1, & 
                           smol_coeffs(no_shock_idx,:,:), n_states, 0.0_dp, & 
                           S_Star, 1)  
    polyn_points = Smolyak_Polynomial(S_Star, n_active_dims, max_smol_level, smol_elem_ani) 

    ! Obtain variables at S_Star
    ! Dimension n_interp
    CALL DGEMM('N','N', 1, nrhs2, n_states, 1.0_dp, polyn_points, 1, & 
                           interp_coeffs, n_states, 0.0_dp, & 
                           other_variables_star, 1)   
            ! interpolation of MPKs of constrained agent does not work well due to sparse grid
            ! adjust manually
            if (other_variables_star(67) > 0.5) then
                other_variables_star(36) = 0.0_dp
                other_variables_star(45) = 0.0_dp
                other_variables_star(54) = 0.0_dp
                other_variables_star(57) = lmbd_vec(3)*kbar_constr(3)
            endif 

    ! State variables at S_Star
    ! Dimension smolyak_d
    CALL DGEMM('N','N', 1, nrhs3, n_states, 1.0_dp, polyn_points, 1, & 
                       state_coeffs, n_states, 0.0_dp, & 
                       state_variables_star, 1)   

    state_series_star_mat(iii,:)      = state_variables_star 
    other_vars_series_star_mat(iii,:) = other_variables_star

    ! Extract per-person wealth by type 
    wealth_a = other_variables_star(15) * state_variables_star(2) / lmbd_vec(1)
    wealth_b = other_variables_star(15) * state_variables_star(3) / lmbd_vec(2)
    wealth_c = other_variables_star(15) * (1.0_dp - state_variables_star(2) - state_variables_star(3)) / lmbd_vec(3)
    
    ! ============================================================ !
    ! Step 2: Extract relevant information about the shocked state S1
    ! ============================================================ !
    

    ! Next State 
    do qqq = 1,n_quad
        b = transpose(next_state_monetary_mat(qqq,:,:))
        call F07ABF( 'E','N',n_states, nrhs, smol_polynom, lda,af,ldaf,ipiv,      &
                     equed, r, c, b, ldb, smol_monetary_coeffs(qqq,:,:), & 
                     ldx, rcond, ferr, berr, work, iwork, info)
    enddo

    ! Transition Coefficients
    b = transpose(next_state_montransition_mat(1,:,:))
    call F07ABF( 'E','N',n_states, nrhs, smol_polynom, lda,af,ldaf,ipiv,      &
                 equed, r, c, b, ldb, smol_montransition_coeffs(1,:,:), & 
                 ldx, rcond, ferr, berr, work, iwork, info)

    ! Relevant variable coefficients 
    b2 = interp_input_monetary_mat
    call F07ABF( 'E','N',n_states, nrhs2, smol_polynom, lda,af,ldaf,ipiv,      &
                 equed, r, c, b2, ldb, interp_monetary_coeffs, & 
                 ldx, rcond, ferr2, berr2, work, iwork, info )

    ! Find Shocked State S_1
    polyn_points = Smolyak_Polynomial(S_0, n_active_dims, max_smol_level, smol_elem_ani) 
    CALL DGEMM('N','N', 1, nrhs, n_states, 1.0_dp, polyn_points, 1, & 
                               smol_montransition_coeffs(1,:,:), n_states, 0.0_dp, & 
                               S_1, 1)  


    ! Find other / state variables in state S_1
    polyn_points = Smolyak_Polynomial(S_1, n_active_dims, max_smol_level, smol_elem_ani) 
    ! polyn_points = Smolyak_Polynomial(S_0, n_active_dims, max_smol_level, smol_elem_ani) 
    CALL DGEMM('N','N', 1, nrhs2, n_states, 1.0_dp, polyn_points, 1, & 
                           interp_monetary_coeffs, n_states, 0.0_dp, & 
                           other_variables, 1)   
            ! interpolation of MPKs of constrained agent does not work well due to sparse grid
            ! adjust manually
            if (other_variables(67) > 0.5) then
                other_variables(36) = 0.0_dp
                other_variables(45) = 0.0_dp
                other_variables(54) = 0.0_dp
                other_variables(57) = lmbd_vec(3)*kbar_constr(3)
            endif    

    CALL DGEMM('N','N', 1, nrhs3, n_states, 1.0_dp, polyn_points, 1, & 
                       state_coeffs, n_states, 0.0_dp, & 
                       state_variables, 1)   
    
    state_series_1_mat(iii,:)      = state_variables
    other_vars_series_1_mat(iii,:) = other_variables
    
    ! Extract Relevant Information 
    p_dis = exp(state_variables(6))
    k_aggr = state_variables(1)
    l_aggr = other_variables(2)
    q_current = other_variables(10)
    infl = other_variables(18)
    nom_i = other_variables(19)
    c_vec(1) = other_variables(3)
    c_vec(2) = other_variables(4)
    c_vec(3) = other_variables(5)
    share_vec(1) = other_variables(11)
    share_vec(2) = other_variables(12)
    share_vec(3) = other_variables(13)
    
    ! nxt_mat, n_nxt, nxt_mat_2
    ! Need the relevant interp_input_mat 
    ! Now we are standing in state S1 and there is no more shock, so use interp_input_mat directly
    interp_use(:,1) = interp_input_mat(:,22)        ! Value a 
    interp_use(:,2) = interp_input_mat(:,23)        ! Value b
    interp_use(:,3) = interp_input_mat(:,24)        ! Value c
    interp_use(:,4) = interp_input_mat(:,68)        ! mc a 
    interp_use(:,5) = interp_input_mat(:,69)        ! mc b
    interp_use(:,6) = interp_input_mat(:,70)        ! mc c
    interp_use(:,2*n_I + 1)  = interp_input_mat(:,10) ! q
    interp_use(:,2*n_I + 2)  = interp_input_mat(:,2)  ! l
    interp_use(:,2*n_I + 3)  = interp_input_mat(:,18) ! Inflation
    interp_use(:,2*n_I + 4)  = interp_input_mat(:,25) ! q_l a 
    interp_use(:,2*n_I + 5)  = interp_input_mat(:,26) ! q_l b
    interp_use(:,2*n_I + 6)  = interp_input_mat(:,27) ! q_l c

    ! Coefficients 
    b4 = interp_use
    call F07ABF( 'Equilibration','No transpose',n_states, n_nxt, smol_polynom, lda,af,ldaf,ipiv,    &
                 equed,r,c,b4,ldb,smol_coeffs_use,ldx,rcond,ferr4,berr4,work,iwork,info )

    ! Create nxt_mat 
    nxt_mat = 0.0_dp 
    nxt_mat_2 = 0.0_dp 

    do qqq = 1,n_quad
        ! transition out of S1 going forward
        polyn_points = Smolyak_Polynomial(S_1, n_active_dims, max_smol_level, smol_elem_ani) 

        ! S_Trans is next state when n_quad qqq happens 
        CALL DGEMM('N','N', 1, nrhs, n_states, 1.0_dp, polyn_points, 1, & 
                   smol_monetary_coeffs(qqq,:,:), n_states, 0.0_dp, & 
                   S_Trans, 1) 
        polyn_points = Smolyak_Polynomial(S_Trans, n_active_dims, max_smol_level, smol_elem_ani) 

        ! Fill in nxt_mat 
        CALL DGEMM('N','N', 1, n_nxt, n_states, 1.0_dp, polyn_points, 1, & 
                           smol_coeffs_use, n_states, 0.0_dp, nxt_mat(qqq,:), 1)  
    enddo

    ! nxt_mat_2: mostly just the current capital holding, except disaster state
    k_next = other_variables(20) / other_variables(10)
    nxt_mat_2 = k_next/exp(dz_vec_adj)*exp(dz_vec)



    ! ============================================================ !
    ! Step 3A: Compute policies functions like k^i(n^i(S^*), S_1)
    ! Everything in step 1, except wealth
    ! ============================================================ !
    ! Current Variables 
    call calc_counterfactual_policies(1, wealth_a, p_dis, k_aggr, l_aggr, q_current, infl, nom_i, nint(other_variables(65)), &
                                      nxt_mat, n_nxt, nxt_mat_2, c_vec, share_vec, results_a)
    call calc_counterfactual_policies(2, wealth_b, p_dis, k_aggr, l_aggr, q_current, infl, nom_i, nint(other_variables(66)), &
                                      nxt_mat, n_nxt, nxt_mat_2, c_vec, share_vec, results_b)
    call calc_counterfactual_policies(3, wealth_c, p_dis, k_aggr, l_aggr, q_current, infl, nom_i, nint(other_variables(67)), &
                                      nxt_mat, n_nxt, nxt_mat_2, c_vec, share_vec, results_c)
    ! Counterfactuals
    cf_policy_matA(iii,1,:) = results_a 
    cf_policy_matA(iii,2,:) = results_b
    cf_policy_matA(iii,3,:) = results_c 


    ! ============================================================ !
    ! Step 3B: Hold Expected Returns Constant
    ! ============================================================ !
    q_use = q_current * other_variables(9) / other_variables_star(9)
    i_use = nom_i * other_variables_star(8) / other_variables(8)
    call calc_counterfactual_policies(1, wealth_a, p_dis, k_aggr, l_aggr, q_use, infl, i_use, nint(other_variables(65)), &
                                      nxt_mat, n_nxt, nxt_mat_2, c_vec, share_vec, results_a)
    call calc_counterfactual_policies(2, wealth_b, p_dis, k_aggr, l_aggr, q_use, infl, i_use, nint(other_variables(66)),&
                                      nxt_mat, n_nxt, nxt_mat_2, c_vec, share_vec, results_b)
    call calc_counterfactual_policies(3, wealth_c, p_dis, k_aggr, l_aggr, q_use, infl, i_use, nint(other_variables(67)),&
                                      nxt_mat, n_nxt, nxt_mat_2, c_vec, share_vec, results_c)

    ! Counterfactuals
    cf_policy_matB(iii,1,:) = results_a 
    cf_policy_matB(iii,2,:) = results_b
    cf_policy_matB(iii,3,:) = results_c 
    
    enddo
    open (unit = 10, file = trim(results_folder) // 'Counterfactual_Policies_A.txt', ACTION="write", STATUS="replace", &
            & FORM="formatted", ACCESS="STREAM")
    write(10,*) sum(cf_policy_matA     ,1)/n_sample; close(10)
    
    open (unit = 10, file = trim(results_folder) // 'Counterfactual_Policies_B.txt', ACTION="write", STATUS="replace", &
            & FORM="formatted", ACCESS="STREAM")
    write(10,*) sum(cf_policy_matB     ,1)/n_sample; close(10)

    ! ============================================================ !
    ! Step 4: Collect results for further processing
    ! ============================================================ !
    open (unit = 10, file = trim(results_folder) // 'Policies_S_0.txt', ACTION="write", STATUS="replace", &
            & FORM="formatted", ACCESS="STREAM")
    write(10,*) sum(other_vars_series_0_mat,1)/n_sample; close(10)
    open (unit = 10, file = trim(results_folder) // 'States_S_0.txt', ACTION="write", STATUS="replace", &
            & FORM="formatted", ACCESS="STREAM")
    write(10,*) sum(state_series_0_mat,1)/n_sample; close(10)

    ! Policies (calc_bond_price format) and states at state S_Star
    open (unit = 10, file = trim(results_folder) // 'Policies_S_Star.txt', ACTION="write", STATUS="replace", &
            & FORM="formatted", ACCESS="STREAM")
    write(10,*) sum(other_vars_series_star_mat,1)/n_sample; close(10)
    open (unit = 10, file = trim(results_folder) // 'States_S_Star.txt', ACTION="write", STATUS="replace", &
            & FORM="formatted", ACCESS="STREAM")
    write(10,*) sum(state_series_star_mat,1)/n_sample; close(10)

    ! Policies (calc_bond_price format) and states at state S_1
    open (unit = 10, file = trim(results_folder) // 'Policies_S_1.txt', ACTION="write", STATUS="replace", &
            & FORM="formatted", ACCESS="STREAM")
    write(10,*) sum(other_vars_series_1_mat,1)/n_sample; close(10)

    open (unit = 10, file = trim(results_folder) // 'States_S_1.txt', ACTION="write", STATUS="replace", &
            & FORM="formatted", ACCESS="STREAM")
    write(10,*) sum(state_series_1_mat,1)/n_sample; close(10)



end subroutine monetary_decomposition



subroutine calc_counterfactual_policies(iii, wealth, p_dis, k_aggr, l_aggr, q_current, infl, nom_i, constraint_binding, &
										nxt_mat, n_nxt, nxt_mat_2, c_vec, share_vec, results_vec)

    ! ============================================================ !
    ! Declarations
    ! ============================================================ !
    use mod_param, only: n_I, n_quad, &
                         aalpha, smolyak_d, quad_weight_vec, dz_vec, ddelta, state_grid, &
                         chiX, xi, ies_vec, gma_vec, dz_vec_adj, bbeta_vec, thtbar_vec, tht, &
                         wealth_share_grid, n_bond, n_interp, lmbd_vec, s_bar_vec, &
                         vareps_w, tau_w, chiW, tayl_ic, phi, disast_p, &
                         disast_std, s_trgt_vec, &
                         constr_agents, kbar_constr, low_guess_fixed, high_guess_fixed, &
                         labor_alloc_vec, use_idio_risk, idio_weight_vec, idio_shock_grid, n_idio
    
    use mod_calc, only: calc_portfolio_share, tot_wealth_ss_vec

    integer, intent(in) :: iii, n_nxt, constraint_binding
    real(dp), intent(in) :: p_dis, wealth, nxt_mat_2(n_quad), l_aggr, &
                            nxt_mat(n_quad, n_nxt), q_current, c_vec(n_I), infl, &
                            share_vec(n_I), nom_i
    real(dp), intent(out) :: results_vec(4) ! Consumption, Share, Bond, Capital

    integer :: xxx, ggg, ifail, ccc, nnn

    real(dp) :: k_aggr, & 
                q_nxt(n_quad), rk_vec(n_quad), &
                aggr_wealth, v_temp(n_quad, n_I),  &
                w_next_choice(n_quad), & 
                pi_nxt(n_quad), y_next(n_quad), &
                infl_nxt(n_quad), mc_temp(n_quad, n_I), &
                k_nxt(n_quad), &
                rf_vec(n_quad), pi_current, y_current, &
                excess_b, k_temp_vec(n_I), b_temp_vec(n_I), c_temp(n_I), & 
                thtbar, savings_vec(n_I)

    real(dp) :: E_rk_new, E_rf_new, next_k

    real(dp), parameter :: brent_delta = 5E-16_dp, min_cons_sav = 1E-08_dp
 
    real(dp) :: w_choice, consumption, r_alpha_vec(n_quad), &
             diff_inner, savings, objf, cons_update, IES, bbeta, &
             next_period_share(n_quad), v_vec(n_quad), EV, &
             M_vec_mat(n_quad, n_I), &
             survivor_wealth(n_quad), &
             gma, labor, q_l_nxt(n_quad, n_I), & 
             per_person_wealth_vec(n_I), v_vec_twisted(n_quad), & 
             l_aggr_nxt(n_quad), survivor_wealth_vec(n_quad,n_I), temp_nxt_aggr_wealth(n_quad), &
             diff, w_current_next(n_quad), &
             w_current, w_choice_new, util, util_c_deriv, labor_part, &
             l_aggr_vec_temp, l_aggr_temp, l_aggr_temp_new, w_choice_temp, w_choice_temp_new, &
             big_weight_vec(n_quad), deriv_helper(n_quad), & 
             cons_temp, share_temp
 
    integer :: iter, mflag, ix1, ix2, iter_inner, kkk

    real(dp) :: share_update, temp_vec(n_quad), share_high, share_low, &
                curr_share, temp, deriv, &
                MPC_mat(n_I), MPS_mat(n_I), MPK_mat(n_I), &
                MPC_mat2(n_I), MPS_mat2(n_I), MPK_mat2(n_I), MPL_mat2(n_I), &
                time_choice, time_update, q_l, RHS, Time_Deviation_mat(n_I)

    real(dp) :: v_idio(n_quad*n_idio, n_I), q_l_nxt_idio(n_quad*n_idio, n_I), &
                dz_vec_idio(n_quad*n_idio), big_weight_vec_idio(n_quad*n_idio), &
                r_alpha_vec_idio(n_quad*n_idio), next_period_share_idio(n_quad*n_idio), &
                v_vec_twisted_idio(n_quad*n_idio), rf_vec_idio(n_quad*n_idio), &
                deriv_helper_idio(n_quad*n_idio), M_idio(n_quad*n_idio), mc_idio(n_quad*n_idio, n_I)


    ! ============================================================ !
    ! Preparations: Current Period
    ! ============================================================ !
    big_weight_vec(1:n_quad-1) = quad_weight_vec*(1.0-p_dis)
    big_weight_vec(n_quad)     = p_dis

    pi_current = aalpha * (l_aggr/k_aggr)**(1.0-aalpha)
    aggr_wealth = k_aggr * (pi_current + (1.0-ddelta)*q_current)
    y_current = k_aggr**aalpha * l_aggr**(1.0-aalpha)
    w_choice = (1.0-aalpha)*y_current/l_aggr 

    ! ============================================================ !
    ! Preparations: Net Period
    ! ============================================================ !
    k_nxt = nxt_mat_2 
    v_temp(:,1) = nxt_mat(:,1) 
    v_temp(:,2) = nxt_mat(:,2)
    v_temp(:,3) = nxt_mat(:,3)
    mc_temp(:,1) = nxt_mat(:,4) 
    mc_temp(:,2) = nxt_mat(:,5)
    mc_temp(:,3) = nxt_mat(:,6)
    q_nxt        = nxt_mat(:, 2*n_I+1)
    l_aggr_nxt   = nxt_mat(:, 2*n_I+2)
    infl_nxt     = nxt_mat(:, 2*n_I+3)
    q_l_nxt(:,1) = nxt_mat(:, 2*n_I+4) * exp(dz_vec)
    q_l_nxt(:,2) = nxt_mat(:, 2*n_I+5) * exp(dz_vec)
    q_l_nxt(:,3) = nxt_mat(:, 2*n_I+6) * exp(dz_vec)


    ! Output
    y_next = exp((dz_vec)*(1.0_dp-aalpha))*l_aggr_nxt* & 
                  ( (k_nxt/l_aggr_nxt)**aalpha)

    ! Wage
    w_next_choice = (1-aalpha)*y_next/l_aggr_nxt

    pi_nxt = aalpha*exp( (dz_vec)*(1.0_dp-aalpha)) &
           *(l_aggr_nxt**(1.0_dp-aalpha))*k_nxt**(aalpha-1.0)

    ! Gross Capital Return (Adjust in disaster state)
    rk_vec = ( (1.0_dp-ddelta)*q_nxt + pi_nxt ) / q_current 
    rk_vec(n_quad) = exp(dz_vec(n_quad)) * rk_vec(n_quad)

    ! Gross Risk-Free Return
    rf_vec = nom_i / infl_nxt

    ! Basic Idiosyncratic Risk Expansion
    if (use_idio_risk == 1) then       
        do kkk = 1,n_idio
            v_idio(((kkk-1)*n_quad+1):(kkk*n_quad),1) = v_temp(:,1)
            v_idio(((kkk-1)*n_quad+1):(kkk*n_quad),2) = v_temp(:,2)
            v_idio(((kkk-1)*n_quad+1):(kkk*n_quad),3) = v_temp(:,3)
            mc_idio(((kkk-1)*n_quad+1):(kkk*n_quad),1) = mc_temp(:,1)
            mc_idio(((kkk-1)*n_quad+1):(kkk*n_quad),2) = mc_temp(:,2)
            mc_idio(((kkk-1)*n_quad+1):(kkk*n_quad),3) = mc_temp(:,3)
            q_l_nxt_idio(((kkk-1)*n_quad+1):(kkk*n_quad),1) = q_l_nxt(:,1)
            q_l_nxt_idio(((kkk-1)*n_quad+1):(kkk*n_quad),2) = q_l_nxt(:,2)
            q_l_nxt_idio(((kkk-1)*n_quad+1):(kkk*n_quad),3) = q_l_nxt(:,3)
            dz_vec_idio(((kkk-1)*n_quad+1):(kkk*n_quad)) = dz_vec
            rf_vec_idio(((kkk-1)*n_quad+1):(kkk*n_quad)) = rf_vec
            big_weight_vec_idio(((kkk-1)*n_quad+1):(kkk*n_quad)) = big_weight_vec * idio_weight_vec(kkk)
        enddo
    endif

    ! ============================================================ !
    ! Calculation
    ! ============================================================ !
    IES   = IES_vec(iii)
    bbeta = bbeta_vec(iii)
    gma   = gma_vec(iii)
    thtbar  = thtbar_vec(iii)
    labor  = (l_aggr*labor_alloc_vec(iii))
    labor_part =  1.0_dp + (1.0_dp/IES - 1.0_dp) * thtbar*tht/(1.0_dp + tht) * labor**( (1.0_dp+tht)/tht )
    cons_temp = c_vec(iii)
    share_temp = share_vec(iii)
    
    diff_inner = 1.0_dp 
    do while(diff_inner > 1E-05_dp)
        ! ============================================================ !
        ! Update Consumption Based on Portfolio Choice 
        ! ============================================================ !
        curr_share = share_temp 
        consumption = cons_temp
        r_alpha_vec = curr_share*rf_vec + (1.0_dp - curr_share)*rk_vec
        savings = wealth + w_choice*labor - consumption

        if (use_idio_risk == 1) then
            ! ------------------------------------------------------------ !
            ! With Idiosyncratic Risk
            ! ------------------------------------------------------------ !
            ! Expand portfolio return 
            do kkk = 1,n_idio
                r_alpha_vec_idio(((kkk-1)*n_quad+1):(kkk*n_quad)) = curr_share*rf_vec + (1.0_dp - curr_share)*rk_vec*exp(idio_shock_grid(kkk,iii))
            enddo

            ! Other Operations
            next_period_share_idio = (savings*r_alpha_vec_idio + q_l_nxt_idio(:,iii) )/tot_wealth_ss_vec(iii)
            v_vec_twisted_idio = (v_idio(:,iii)*next_period_share_idio) ** (1.0-gma)
            EV = sum(big_weight_vec_idio*v_vec_twisted_idio)

            deriv_helper_idio = bbeta * (v_idio(:,iii)**((1.0_dp/IES) - gma )) * (next_period_share_idio**(-gma)) & 
                              * (EV **( (gma - (1.0_dp/IES))/(1.0_dp-gma) )) * mc_idio(:,iii)

            if (constr_agents(iii) == 1 .and. constraint_binding == 1) then
                temp = sum(deriv_helper_idio*rf_vec_idio*big_weight_vec_idio)
            else
                temp = sum(deriv_helper_idio*r_alpha_vec_idio*big_weight_vec_idio)
            endif       

        else
            ! ------------------------------------------------------------ !
            ! No Idiosyncratic risk
            ! ------------------------------------------------------------ !
            next_period_share = (savings*r_alpha_vec + q_l_nxt(:,iii) )/tot_wealth_ss_vec(iii)
            v_vec_twisted = (v_temp(:,iii)*next_period_share) ** (1.0-gma)
            EV = sum(big_weight_vec*v_vec_twisted)

            deriv_helper = bbeta * (v_temp(:,iii)**((1.0_dp/IES) - gma )) * (next_period_share**(-gma))  & 
                           * (EV **( (gma - (1.0_dp/IES))/(1.0_dp-gma) )) * mc_temp(:,iii)
            if (constr_agents(iii) == 1 .and. constraint_binding == 1) then
                temp = sum(deriv_helper*rf_vec*big_weight_vec)
            else
                temp = sum(deriv_helper*r_alpha_vec*big_weight_vec)
            endif

        endif

        cons_update = labor_part * ( temp**(-IES) ) 
        cons_update = minval([maxval([cons_update, min_cons_sav]), (wealth + w_choice*labor-min_cons_sav)])                

        ! ============================================================ !
        ! Update Portfolio Choice 
        ! ============================================================ !
        if (constr_agents(iii) == 1 .and. constraint_binding == 1) then
            share_update = 1.0_dp - kbar_constr(iii)*q_current/savings
        else

            ! find natural leverage limit to prevent negative payoffs
            temp_vec = -rk_vec/(rf_vec-rk_vec + sqrt_eps)
            where (temp_vec >= 0.0_dp)
                temp_vec = -10000.0_dp
            end where
            share_low = maxval([low_guess_fixed, maxval(temp_vec)+ sqrt_eps])

            ! find natural leverage limit to prevent negative payoffs
            temp_vec = -rk_vec/(rf_vec-rk_vec + sqrt_eps)
            where (temp_vec <= 0.0_dp)
                temp_vec = 10000.0_dp
            end where
            share_high = minval([high_guess_fixed, minval(temp_vec) - sqrt_eps ])

            ! Calculate Portfolio Share, Given Updated Consumption and Savings
            share_update = curr_share
            call calc_portfolio_share(1.0_dp, share_low, share_high, iii, big_weight_vec, savings, labor, cons_update,  &
                              0.0_dp, q_l_nxt(:,iii), v_temp(:,iii), mc_temp(:,iii), rf_vec, rk_vec, q_current, &
                              share_update)

        endif

        ! ============================================================ !
        ! Update Values
        ! ============================================================ !
        diff_inner = abs(consumption-cons_update)   
        cons_temp = consumption + 0.2*(cons_update-consumption)
        share_temp = curr_share + 0.2*(share_update-curr_share)

    enddo
    

    ! ============================================================ !
    ! Collect policies
    ! ============================================================ !
    results_vec(1) = consumption
    results_vec(2) = share_update
    results_vec(3) = savings * share_update 
    results_vec(4) = savings * (1.0_dp - share_update) / q_current


end subroutine calc_counterfactual_policies



end module mod_decomp

