! -------------------------------------------------------------------------
! mod_calc.f90: solution module, public subroutines:
! - calc_steady: solves for steady state to center smolyak grid 
! - calc_sol: solves for stationary equilibrium and creates impulse respones
! - calc_portfolio_share: solves for portfolio share between safe and risky assets
! -------------------------------------------------------------------------
! authors:         Rohan Kekre and Moritz Lenel
! for updates see: https://github.com/KekreLenel/MPR
! -------------------------------------------------------------------------
module mod_calc

use omp_lib
use base_lib,  only: dp, sqrt_eps, eps, Fill_linspace_dp
use mod_param, only: n_I

implicit none
private

public :: calc_steady , calc_sol, calc_portfolio_share, tot_wealth_ss_vec

! global steady state variables 
real(dp) :: v_ss_vec(n_I), mc_ss_vec(n_I), rf_ss, l_ss, tot_wealth_ss_vec(n_I), q_l_ss_vec(n_I)

contains

subroutine calc_steady()
! steady state calculation used for grid scaling and initial guesses 

    ! load parameters
    use mod_param, only :   bbeta_vec, ies_vec, ddelta, aalpha, s_trgt_vec, l_target,    &
                            tht, thtbar_vec, k_grid_dev, w_grid_dev, k_grid_mean,        & 
                            w_grid_mean, lmbd_vec, k_dev_param, w_dev_param, w_grid_adj, &
                            k_grid_adj, gma_vec, results_folder, labor_alloc_vec

    ! Local Variables   
    real(dp) :: bbeta_avrg, gma_avrg, rk_ss, per_person_wealth(n_I), & 
                c_ss_vec(n_I),  s_ss_vec(n_I), diff, y_ss, util, labor_util, & 
                thtbar, thtbar_new, c_ss, ies_avrg, k_ss, w_ss
    integer  :: iter, iii  

    ! average preference parameters 
    bbeta_avrg = sum(s_trgt_vec*bbeta_vec) 
    gma_avrg   = sum(s_trgt_vec*gma_vec)
    ies_avrg   = sum(s_trgt_vec*ies_vec)

    ! corresponding returns
    rk_ss = 1.0_dp/bbeta_avrg  - 1.0_dp 
    rf_ss = 1.0_dp/bbeta_avrg  - 1.0_dp 
    
    ! steady state labor
    l_ss = l_target 

    ! Capital, wage, output in stationary economy
    k_ss = l_ss * ((rk_ss + ddelta) / aalpha)**(1.0_dp / (aalpha - 1.0_dp))
    w_ss = (1.0_dp - aalpha) * (k_ss**aalpha) * (l_ss**(-aalpha))
    y_ss = (k_ss**aalpha) * l_ss**(1-aalpha)

    ! Wealth shares, consumption
    s_ss_vec = s_trgt_vec
    do iii = 1,n_I
        per_person_wealth(iii) = s_ss_vec(iii)/lmbd_vec(iii)
        c_ss_vec(iii) = w_ss*l_target*labor_alloc_vec(iii) + per_person_wealth(iii) * rk_ss * k_ss
    enddo

    ! ----------------------------------------------------------------- !
    ! Find disutility parameters thtbar to match labor target 
    ! ----------------------------------------------------------------- !
    do iii = 1,n_I
    if (labor_alloc_vec(iii) > sqrt_eps) then 
        thtbar = 0.5_dp 
        diff = 1.0_dp 
        iter = 0
        do while (diff > 1E-12_dp)

            iter = iter + 1
            thtbar_new = ies_vec(iii)*w_ss/c_ss_vec(iii) * (l_ss*labor_alloc_vec(iii))**(-1.0_dp/tht) * &
                       (1.0_dp + (1.0_dp/ies_vec(iii)-1.0_dp) * thtbar * ((l_ss*labor_alloc_vec(iii))**(1.0_dp+1.0/tht)) / (1.0_dp+1.0/tht))
            diff = abs(thtbar-thtbar_new)
            thtbar = thtbar + 0.5*(thtbar_new-thtbar)

            if (iter > 1000) then
                write(*,*) 'ERROR: No convergence in finding thtbar.'
                stop
            endif

        enddo
        thtbar_vec(iii) = thtbar 
    else
        thtbar_vec(iii) = 0.0_dp
    endif
    enddo

    ! ----------------------------------------------------------------- !
    ! Value function; labor endowment prices
    ! ----------------------------------------------------------------- !

    ! labor endowment price
    q_l_ss_vec = w_ss*l_ss*labor_alloc_vec / (1.0-bbeta_vec )

    ! total individual wealth including endowment claim
    tot_wealth_ss_vec = q_l_ss_vec + per_person_wealth*k_ss*(1.0+rk_ss)

    do iii = 1,n_I

        
        call  util_fun(c_ss_vec(iii), l_ss*labor_alloc_vec(iii), iii, util, mc_ss_vec(iii), labor_util)
        
        ! value function
        v_ss_vec(iii)  =  util**( 1.0/(1.0-1.0/ies_vec(iii)) )


    enddo  
    

    ! ----------------------------------------------------------------- !
    ! Scaling parameters of wage and capital grid 
    ! ----------------------------------------------------------------- !
    k_grid_mean   = k_grid_adj*k_ss
    w_grid_mean   = w_grid_adj*w_ss

    k_grid_dev   = k_dev_param*k_grid_mean
    w_grid_dev   = w_dev_param*w_grid_mean

    ! ----------------------------------------------------------------- !
    ! Print SS
    ! ----------------------------------------------------------------- !
    write(*,*) '---------------------------------------------'
    write(*,*) 'STEADY STATE'
    write(*,*) '---------------------------------------------'
    write(*,"(A21, F10.4)")                       ' k_ss         =      ' , k_ss
    write(*,"(A21, F10.4)")                       ' k_grid_mean  =      ' , k_grid_mean
    write(*,"(A21, F10.4)")                       ' y_ss         =      ' , y_ss
    write(*,"(A21, F10.4)")                       ' l_ss         =      ' , l_ss
    write(*,*)
    write(*,*) '---------------------------------------------'
    write(*,"(A21, F10.4, A1, F10.4, A1, F10.4)") ' q_l_ss       =      ' , q_l_ss_vec(1), ' ', q_l_ss_vec(2), ' ', q_l_ss_vec(3)
    write(*,"(A21, F10.4, A1, F10.4, A1, F10.4)") ' thtbar       =      ' , thtbar_vec(1), ' ', thtbar_vec(2), ' ', thtbar_vec(3)
    write(*,"(A21, F10.4, A1, F10.4, A1, F10.4)") ' c_ss         =      ' , c_ss_vec(1),   ' ', c_ss_vec(2), ' ', c_ss_vec(3)
    write(*,"(A21, F10.4, A1, F10.4, A1, F10.4)") ' v_ss         =      ' , v_ss_vec(1),   ' ', v_ss_vec(2), ' ', v_ss_vec(3)
    write(*,"(A21, F10.4, A1, F10.4, A1, F10.4)") ' mc_ss        =      ' , mc_ss_vec(1),  ' ', mc_ss_vec(2), ' ', mc_ss_vec(3)
    write(*,*) '---------------------------------------------'
    write(*,*)

    ! Export calibrated thtbar for result files
    open (unit = 10, file = trim(results_folder) // 'extra_data.csv', ACTION="write",  &
        & FORM="formatted", ACCESS="sequential")
    write(10,'(3F10.4)') thtbar_vec(1), thtbar_vec(2), thtbar_vec(3)
    close(10)


end subroutine calc_steady

subroutine calc_sol()
    ! ------------------------------------------------------------ !
    ! Load Parameters 
    ! ------------------------------------------------------------ !
    use mod_smolyak, only: Smolyak_Polynomial2
    use mod_param, only: n_states, idx_k, idx_m, idx_w, idx_dis,              &
                         state_grid, n_quad, aalpha, smolyak_d, n_active_dims,                  &
                         next_m_mat, next_dis_mat, smol_polynom, max_smol_level, smol_elem_ani, &
                         k_grid_mean, k_grid_dev, s_a_grid_mean, s_a_grid_dev,          &
                         s_c_grid_mean, s_c_grid_dev, m_grid_mean, m_grid_dev,          &
                         w_grid_mean, w_grid_dev, dis_grid_mean, dis_grid_dev,                  &
                         n_shocks, dz_vec, dz_vec_adj, no_shock_idx, smol_grid, command_input,    &
                         wealth_share_grid, results_folder, labor_alloc_vec,                    &
                         vector_mus_dimensions, irf_indices, irf_shock_sizes,                   &
                         sidx_z, sidx_m, sidx_dis, n_interp, n_bond, n_shocks, rho_m
                         

    ! ------------------------------------------------------------ !
    ! Declare Variables 
    ! ------------------------------------------------------------ !
    integer, parameter  :: n_nxt = 6 + 2*n_I, maxiter = 10000, irf_ll_length = 20
    real(dp), parameter :: conv_crit = 1E-06_dp

    real(dp) :: v_mat(n_I, n_states), k_next_mat(n_quad, n_states), q_mat(n_states), l_aggr_mat(n_states), &
                next_state_mat(n_quad, n_active_dims, n_states), next_state_mat_new(n_quad, n_active_dims, n_states), &
                nxt_mat(n_quad, n_nxt), nxt_mat_2(n_quad), nxt_mat_trans(1, n_nxt), nxt_mat_2_trans, &
                interp_input_mat(n_states, n_nxt), nom_i_mat(n_states), current_state_vec(smolyak_d),  & 
                nom_i, w_choice_mat(n_states), q_current, infl_mat(n_states), c_vec(n_I), c_mat(n_I, n_states), &
                share_mat(n_I, n_states), q_l_mat(n_I, n_states), share_vec(n_I), infl, w_choice, l_aggr, q_new, & 
                k_next_new, c_vec_new(n_I), infl_new, w_choice_new, v_new(n_I), s_nxt(n_quad,3), l_aggr_new, &
                v_mat_new(n_I, n_states), l_aggr_mat_new(n_states), q_mat_new(n_states), nom_i_mat_new(n_states), &
                c_mat_new(n_I, n_states), w_choice_mat_new(n_states), k_next_mat_new(n_quad, n_states), &
                q_l_current(n_I), q_l_new(n_I), q_l_mat_new(n_I, n_states), infl_mat_new(n_states), &  
                k_mat(n_states), share_mat_new(n_I, n_states), q_bond_mat(n_states, n_bond+1), mc_new(n_I),  & 
                results_vec(n_interp-n_bond*2), results_mat(n_states, n_interp), nxt_bond_prices(n_quad), who_prices_bond, &
                monetary_shock, next_state_monetary_mat(n_quad, n_active_dims, n_states), mc_mat_new(n_I, n_states), &  
                nom_i_monetary_mat(n_states), q_monetary_mat(n_states), c_monetary_mat(n_I,n_states), & 
                share_monetary_mat(n_I,n_states), infl_monetary_mat(n_states), l_aggr_monetary_mat(n_states), & 
                next_state_monetary_mat_new(n_quad, n_active_dims, n_states), q_l_monetary_mat(n_I,n_states), &
                share_monetary_mat_new(n_I,n_states), w_choice_monetary_mat_new(n_states), mc_monetary_mat(n_I,n_states), & 
                infl_monetary_mat_new(n_states), v_monetary_mat(n_I,n_states), q_monetary_mat_new(n_states),  &
                l_aggr_monetary_mat_new(n_states), v_monetary_mat_new(n_I,n_states), nom_i_monetary_mat_new(n_states), & 
                c_monetary_mat_new(n_I,n_states), k_next_monetary_mat_new(n_quad,n_states),  mc_monetary_mat_new(n_I,n_states), & 
                q_l_monetary_mat_new(n_I,n_states), next_state_montransition_mat_new(1,n_active_dims, n_states), & 
                next_state_montransition_mat(1,n_active_dims, n_states), k_next_monetary_mat(n_quad, n_states), & 
                w_choice_monetary_mat(n_states), trans_shock_vec(n_shocks), results_monetary_mat(n_states, n_interp), & 
                q_bond_monetary_mat(n_states, n_bond+1), mc_mat(n_I, n_states), bond_coeffs(n_states) 
                
    integer  :: constraint_binding(n_I), constraint_binding_mat(n_I, n_states), constraint_binding_monetary_mat(n_I, n_states)


    ! Smolyak interpolation variables
    integer, parameter :: nrhs = n_nxt
    integer  ::  lda, ldaf, ldb, ldx
    integer  ::  ipiv(n_states), iwork(n_states)
    integer  ::  info, ifail
    real(dp) ::  a(n_states,n_states), af(n_states,n_states), r(n_states), & 
                 c(n_states), smol_coeffs(n_states,nrhs), work(3*n_states), ferr(nrhs), &  
                 berr(nrhs), polyn_points(n_quad,n_states), polyn_points_trans(1,n_states)
    real(dp) ::  rcond
    character(1) :: equed

    character :: shock_char*1

    integer, parameter :: nrhs2 = 1
    real(dp) ::  ferr2, berr2

    real(dp) :: diff, timer_start, timer_end
    integer  :: outer_iter, sss, iii, counter, bbb, fff

    ! ------------------------------------------------------------ !
    ! 1. Initialization
    ! ------------------------------------------------------------ !
    lda  = n_states
    ldaf = n_states 
    ldb  = n_states 
    ldx  = n_states

    ! (1) Current Period Initial Guesses
    l_aggr_mat      = l_ss
    q_mat           = 1.0_dp
    nom_i_mat       = 1.0 + rf_ss
    infl_mat        = 1.0_dp
    share_mat       = 0.0_dp 
    q_l_mat(1,:)    = q_l_ss_vec(1)
    q_l_mat(2,:)    = q_l_ss_vec(2)
    q_l_mat(3,:)    = q_l_ss_vec(3)

    !! Value Function and marginal utility of consumption
    do iii = 1,n_I
         v_mat(iii,:) =  v_ss_vec(iii)
        mc_mat(iii,:) = mc_ss_vec(iii)
    enddo

    !! Capital 
    k_mat = state_grid(idx_k,:)
    !! Wage 
    w_choice_mat = (1.0-aalpha)* k_mat**aalpha * l_ss**(-aalpha)
    !! Consumption 
    do iii = 1,n_I
        c_mat(iii,:) = w_choice_mat * l_ss + k_mat * wealth_share_grid(iii,:) * rf_ss
    enddo

    ! (2) Initialize Next Period State Transition Matrix
    ! next_state_mat is defined on unscaled Smolyak grid (-1 to 1)
    next_state_mat = 0.0_dp 
    do sss = 1,n_states
        counter = 0

        !! Initialize next period capital = current period capital * disaster shock
        k_next_mat(:, sss) = state_grid(idx_k,sss) * exp(dz_vec) / exp(dz_vec_adj)

        if (vector_mus_dimensions(1) > 0) then
            counter = counter + 1
            next_state_mat(:,counter,sss) = (state_grid(idx_k,sss)/exp(dz_vec_adj) - k_grid_mean ) /k_grid_dev
        endif

        if (vector_mus_dimensions(2) > 0) then
            counter = counter + 1
            next_state_mat(:,counter,sss) = 0.0_dp
        endif

        if (vector_mus_dimensions(3) > 0) then
            counter = counter + 1
            next_state_mat(:,counter,sss) = 0.0_dp 
        endif

        if (vector_mus_dimensions(4) > 0) then
            counter = counter + 1
            next_state_mat(:,counter,sss) = (next_m_mat(:,sss) - m_grid_mean) /m_grid_dev
        endif

        if (vector_mus_dimensions(5) > 0) then
            counter = counter + 1
            next_state_mat(:,counter,sss) = (w_choice_mat(sss) - w_grid_mean) /w_grid_dev
        endif

        if (vector_mus_dimensions(6) > 0) then
            counter = counter + 1
            next_state_mat(:,counter,sss) = (next_dis_mat(:,sss) - dis_grid_mean) /dis_grid_dev
        endif

    enddo
    
    constraint_binding_mat = 0 

    !! Limit transition matrix to remain on grid boundaries
    where (next_state_mat > 1.0_dp)
        next_state_mat = 1.0_dp 
    elsewhere (next_state_mat < -1.0_dp)
        next_state_mat = -1.0_dp 
    endwhere

    next_state_mat_new = next_state_mat !! Pre-allocation

    outer_iter = 0
    diff = 1.0_dp
    do while (diff > conv_crit .and. outer_iter < maxiter)

        outer_iter = outer_iter + 1
        timer_start = omp_get_wtime()

        ! ************************************************************ !
        ! (1) Collect variables for which expectations need to be calculated 
        ! ************************************************************ !
        do iii = 1,n_I
            interp_input_mat(:,iii)       =  v_mat(iii,:)
            interp_input_mat(:,n_I + iii) = mc_mat(iii,:)
        enddo
        interp_input_mat(:,2*n_I + 1)  = q_mat
        interp_input_mat(:,2*n_I + 2)  = l_aggr_mat
        interp_input_mat(:,2*n_I + 3)  = infl_mat
        interp_input_mat(:,2*n_I + 4)  = q_l_mat(1,:)
        interp_input_mat(:,2*n_I + 5)  = q_l_mat(2,:)
        interp_input_mat(:,2*n_I + 6)  = q_l_mat(3,:)

        ! Solve for Smolyak Coefficients (smol_coeffs) for those variables
        ! smol_coeffs has dimension n_states * 9
        call F07ABF( 'Equilibration','No transpose',n_states, nrhs, smol_polynom, lda,af,ldaf,ipiv,    &
                     equed,r,c,interp_input_mat,ldb,smol_coeffs,ldx,rcond,ferr,berr,work,iwork,info )
            
        ! ************************************************************ !
        ! (2) Main updates state by state (parallelized)
        ! ************************************************************ !
        !$OMP PARALLEL DEFAULT(NONE) &
        !$OMP SHARED( & 
        !$OMP next_state_mat, smol_elem_ani, smol_coeffs, state_grid, &
        !$OMP nom_i_mat, w_choice_mat, n_states, dz_vec, & 
        !$OMP k_grid_dev, k_grid_mean, q_mat, next_state_mat_new, &
        !$OMP c_mat, share_mat, infl_mat, l_aggr_mat, q_l_mat, &
        !$OMP w_grid_mean, w_grid_dev, n_active_dims, share_mat_new, &
        !$OMP s_a_grid_mean, s_a_grid_dev, s_c_grid_mean, s_c_grid_dev, &
        !$OMP w_choice_mat_new, infl_mat_new, constraint_binding_mat, &
        !$OMP q_mat_new, l_aggr_mat_new, v_mat_new, nom_i_mat_new, &
        !$OMP c_mat_new, k_next_mat_new, dz_vec_adj, q_l_mat_new, &
        !$OMP k_next_mat, vector_mus_dimensions, mc_mat_new) &
        !$OMP PRIVATE( &
        !$OMP polyn_points, nxt_mat, current_state_vec,    &
        !$OMP nom_i, w_choice, nxt_mat_2, counter, mc_new, &
        !$OMP q_current, c_vec, share_vec, q_l_new,        &
        !$OMP infl , l_aggr, q_l_current, constraint_binding,  &
        !$OMP q_new, k_next_new, c_vec_new, infl_new,      &
        !$OMP sss, v_new, l_aggr_new, w_choice_new, s_nxt)
        !$OMP DO SCHEDULE(static)
        do sss = 1, n_states
            ! ============================================================ !
            ! [1] Interpolate next period values from Smolyak Coefficients
            ! ============================================================ !
            ! Output is nxt_mat (n_quad * 9)
            ! ============================================================ !
            
            ! Get the Smolyak polynomials for potential next states
            ! next_state_mat(:,:,sss) has dimension n_quad * n_active_dims
            ! Its (i,j) element is the next period value of state j if today's state is sss and quadrature point i realizes
            ! polyn_points is the evaluated polynomial matrix, which has dimension n_quad * n_states 
            polyn_points = Smolyak_Polynomial2(next_state_mat(:,:,sss), n_active_dims, &
                                               n_quad, n_states, max_smol_level, smol_elem_ani)

            ! Pre-Multiply polyn_points to smol_coeffs 
            ! The product is assigned to nxt_mat 
            CALL DGEMM('N','N', n_quad, nrhs, n_states, 1.0_dp, polyn_points, n_quad, & 
                       smol_coeffs, n_states, 0.0_dp, nxt_mat, n_quad)   

            ! k doesn't need interpolation
            ! k_next_mat is defined above 
            nxt_mat_2 = k_next_mat(:,sss)

            ! ============================================================ !
            ! [2] Calculate Equilibrium and Update
            ! ============================================================ !
            ! Load current states 
            nom_i              = nom_i_mat(sss)
            q_current          = q_mat(sss)
            c_vec              = c_mat(:,sss)
            share_vec          = share_mat(:,sss)
            infl               = infl_mat(sss)
            l_aggr             = l_aggr_mat(sss)
            q_l_current        = q_l_mat(:,sss)
            constraint_binding = constraint_binding_mat(:,sss)

            ! Update guess for equilibrium solution
            call calc_equilibrium_and_update(nxt_mat, n_nxt, nxt_mat_2, sss, c_vec, l_aggr, q_current, q_l_current, &
                                             infl, nom_i, share_vec, 0.0_dp, q_new, q_l_new, k_next_new, c_vec_new, &
                                             infl_new, v_new, mc_new, l_aggr_new, w_choice_new, s_nxt, constraint_binding)

            ! Update value functions, policy functions and prices
            v_mat_new(:,sss)           = v_new 
            mc_mat_new(:,sss)          = mc_new 
            c_mat_new(:,sss)           = c_vec_new
            k_next_mat_new(:,sss)      = k_next_new/exp(dz_vec_adj)*exp(dz_vec)
            q_mat_new(sss)             = q_new 
            l_aggr_mat_new(sss)        = l_aggr_new
            nom_i_mat_new(sss)         = nom_i
            w_choice_mat_new(sss)      = w_choice_new
            q_l_mat_new(:,sss)         = q_l_new
            share_mat_new(:,sss)       = share_vec          
            infl_mat_new(sss)          = infl_new
            constraint_binding_mat(:,sss) = constraint_binding 

            ! ============================================================ !
            ! [3] Update next period state (next_state_mat_new)
            ! ============================================================ !
            counter = 0

            if (vector_mus_dimensions(1) > 0) then
                counter = counter + 1
                next_state_mat_new(:, counter, sss)   = (k_next_new/exp(dz_vec_adj) -k_grid_mean )/k_grid_dev
            endif

            if (vector_mus_dimensions(2) > 0) then
                counter = counter + 1
                next_state_mat_new(:, counter, sss) = (s_nxt(:,1) - s_a_grid_mean)/s_a_grid_dev
            endif

            if (vector_mus_dimensions(3) > 0) then
                counter = counter + 1
                next_state_mat_new(:, counter, sss) = (s_nxt(:,3) - s_c_grid_mean)/s_c_grid_dev
            endif

            if (vector_mus_dimensions(4) > 0) then
                counter = counter + 1
            endif

            if (vector_mus_dimensions(5) > 0) then
                counter = counter + 1
                next_state_mat_new(:, counter, sss)  = (w_choice_new/exp(dz_vec_adj) - w_grid_mean) /w_grid_dev
            endif

            if (vector_mus_dimensions(6) > 0) then
                counter = counter + 1
            endif

        enddo
        !$OMP END DO
        !$OMP END PARALLEL

        ! ************************************************************ !
        ! (3) Update the matrices 
        ! ************************************************************ !
        ! Restrict Smolyak Grid
        where(next_state_mat_new > 1.0_dp) 
            next_state_mat_new   = 1.0_dp 
        endwhere
        where(next_state_mat_new < -1.0_dp)
            next_state_mat_new   = -1.0_dp
        endwhere

        ! Calculate Differences
        diff =  maxval([ maxval(abs(log(v_mat_new/v_mat))), &
                 maxval(abs(log(mc_mat_new/mc_mat))), &
                 maxval(abs((q_l_mat_new-q_l_mat))), &
                 maxval(abs(log(c_mat_new/c_mat))), &
                 maxval(abs((share_mat_new-share_mat))), &
                 maxval(abs(log(k_next_mat_new/k_next_mat))), &
                 maxval(abs(log(l_aggr_mat_new/l_aggr_mat))), &
                 maxval(abs(log(nom_i_mat_new/nom_i_mat))), &
                 maxval(abs(log(w_choice_mat_new/w_choice_mat))), &
                 maxval(abs(log(q_mat_new/q_mat))), &
                 maxval(abs(log(infl_mat_new/infl_mat))), &
                 maxval(abs(next_state_mat - next_state_mat_new)) ])

        timer_end = omp_get_wtime()

        ! Write differences between new and current values 
        if (mod(outer_iter, 10) == 0) then

            write(*,*) ''
            write(*,"(A11, i5, i5)") ' Iteration ', outer_iter
            write(*,"(A11, f8.2)") ' Calc. time', timer_end-timer_start
            write(*,*) ''

            write(*,"(A19,e12.4,A1)") 'Changes (max diff: ', diff , ')'
            write(*,*) '------------------------------------------------------------------------'
            write(*,"(A13, e12.4, e12.4, e12.4)")     ' v         = ' , &
            maxval(abs(log(v_mat_new(1,:)/v_mat(1,:)))), & 
            maxval(abs(log(v_mat_new(2,:)/v_mat(2,:)))), & 
            maxval(abs(log(v_mat_new(3,:)/v_mat(3,:))))
            write(*,"(A13, e12.4, e12.4, e12.4)")     ' mc         = ' , &
            maxval(abs(log(mc_mat_new(1,:)/mc_mat(1,:)))), & 
            maxval(abs(log(mc_mat_new(2,:)/mc_mat(2,:)))), & 
            maxval(abs(log(mc_mat_new(3,:)/mc_mat(3,:))))
            write(*,"(A13, e12.4, e12.4, e12.4)")     ' q_l       = ' , &
            maxval(abs((q_l_mat_new(1,:)-q_l_mat(1,:)))), & 
            maxval(abs((q_l_mat_new(2,:)-q_l_mat(2,:)))), & 
            maxval(abs((q_l_mat_new(3,:)-q_l_mat(3,:))))
            write(*,"(A13, e12.4, e12.4, e12.4)")     ' c_spend   = ' , &
            maxval(abs(log(c_mat_new(1,:)/c_mat(1,:)))), & 
            maxval(abs(log(c_mat_new(2,:)/c_mat(2,:)))), & 
            maxval(abs(log(c_mat_new(3,:)/c_mat(3,:))))
            write(*,"(A13, e12.4, e12.4, e12.4)")     ' share_mat = ' , &
            maxval(abs((share_mat_new(1,:)-share_mat(1,:)))), & 
            maxval(abs((share_mat_new(2,:)-share_mat(2,:)))), & 
            maxval(abs((share_mat_new(3,:)-share_mat(3,:))))
            write(*,"(A13, e12.4)")     ' k         = ' , &
            maxval(abs(log(k_next_mat_new/k_next_mat)))
            do sss = 1,n_active_dims 
            write(*,"(A13, e12.4)")     ' next_state= ' , &
            maxval(abs((next_state_mat_new(:,sss,:)-next_state_mat(:,sss,:))))
            enddo
            write(*,"(A13, e12.4)")     ' q         = ' , &
            maxval(abs(log(q_mat_new/q_mat)))
            write(*,"(A13, e12.4)")     ' nom_i     = ' , &
            maxval(abs(log(nom_i_mat_new/nom_i_mat)))
            write(*,"(A13, e12.4)")     ' l         = ' , &
            maxval(abs(log(l_aggr_mat_new/l_aggr_mat)))
            write(*,"(A13, e12.4)")     ' w         = ' , &
            maxval(abs(log(w_choice_mat_new/w_choice_mat)))
            write(*,"(A13, e12.4)")     ' infl      = ' , &
            maxval(abs(log(infl_mat_new/infl_mat)))

        endif

        ! Dampened updates
        v_mat           = v_mat          + 1.0*(v_mat_new-v_mat )
        mc_mat          = mc_mat         + 1.0*(mc_mat_new-mc_mat )
        next_state_mat  = next_state_mat + 0.25_dp*(next_state_mat_new-next_state_mat)
        q_mat           = q_mat          + 0.25_dp*(q_mat_new-q_mat)
        l_aggr_mat      = l_aggr_mat     + 0.25_dp*(l_aggr_mat_new-l_aggr_mat)     
        w_choice_mat    = w_choice_mat*(1 + 0.25_dp*log(w_choice_mat_new/w_choice_mat))
        c_mat           = c_mat          + 0.5*(c_mat_new - c_mat)
        q_l_mat         = q_l_mat        + 0.3_dp*(q_l_mat_new-q_l_mat)
        nom_i_mat       = nom_i_mat      + 1.0*(nom_i_mat_new - nom_i_mat)
        infl_mat        = infl_mat       + 0.25*(infl_mat_new - infl_mat)
        k_next_mat      = k_next_mat     + 0.25*(k_next_mat_new - k_next_mat)
        share_mat       = share_mat      + 1.0*(share_mat_new - share_mat) 

    enddo

    ! ------------------------------------------------------------ !
    ! 3. Calculate Bond prices and collect results 
    ! ------------------------------------------------------------ !
    q_bond_mat = 1.0_dp
    !DIR$ NOUNROLL
    do bbb = 1,n_bond
        ! Get coefficients for interpolating bond prices
        call F07ABF( 'Equilibration','No transpose',n_states, nrhs2, smol_polynom, lda,af,ldaf,ipiv,    &
                  equed,r,c, q_bond_mat(:,bbb),ldb,bond_coeffs,ldx,rcond,ferr2,berr2,work,iwork,info)

        !$OMP PARALLEL DEFAULT(NONE) &
        !$OMP SHARED( & 
        !$OMP next_state_mat, smol_elem_ani, smol_coeffs, nom_i_mat, &
        !$OMP n_states, q_mat, c_mat, share_mat, infl_mat, l_aggr_mat, &
        !$OMP n_active_dims, constraint_binding_mat, q_l_mat, &
        !$OMP k_next_mat, results_mat, bbb, bond_coeffs, q_bond_mat) &
        !$OMP PRIVATE( &
        !$OMP polyn_points, nxt_mat, current_state_vec, &
        !$OMP nom_i, w_choice, nxt_mat_2, counter, constraint_binding, &
        !$OMP q_current, c_vec, share_vec, q_l_new, infl, l_aggr, q_l_current, &
        !$OMP q_new, k_next_new, c_vec_new, infl_new, results_vec, &
        !$OMP sss, v_new, l_aggr_new, w_choice_new, s_nxt, &
        !$OMP who_prices_bond, nxt_bond_prices)
        !$OMP DO SCHEDULE(static)
        ! Loop over states
        do sss = 1, n_states

            ! Interpolation 
            polyn_points = Smolyak_Polynomial2(next_state_mat(:,:,sss), n_active_dims, &
                                               n_quad, n_states, max_smol_level, smol_elem_ani)
            ! Standard variables
            CALL DGEMM('N','N', n_quad, nrhs, n_states, 1.0_dp, polyn_points, n_quad, & 
                           smol_coeffs, n_states, 0.0_dp, & 
                           nxt_mat, n_quad)   
            nxt_mat_2 = k_next_mat(:,sss) 
            ! Bond Prices 
            CALL DGEMM('N','N', n_quad, nrhs2, n_states, 1.0_dp, polyn_points, n_quad, & 
                               bond_coeffs, n_states, 0.0_dp, & 
                               nxt_bond_prices, n_quad)   

            ! Current State 
            nom_i = nom_i_mat(sss)
            q_current = q_mat(sss)
            c_vec = c_mat(:,sss)
            share_vec = share_mat(:,sss)
            infl = infl_mat(sss)
            l_aggr = l_aggr_mat(sss)
            q_l_current = q_l_mat(:,sss)
            constraint_binding = constraint_binding_mat(:,sss)

            ! Calculate Bond Prices 
            call calc_bond_prices(nxt_mat, n_nxt, nxt_mat_2, nxt_bond_prices, sss, &
                            c_vec, l_aggr, q_current, q_l_current, infl, &
                            nom_i, share_vec, constraint_binding, who_prices_bond, q_bond_mat(sss,bbb+1), results_vec)
            
            ! Store in results matrix 
            results_mat(sss, 1:(n_interp-n_bond*2)) = results_vec
            results_mat(sss, n_interp-n_bond*2+1)   = q_bond_mat(sss,2)
            results_mat(sss, n_interp-n_bond*2+2)   = q_bond_mat(sss,3)
            results_mat(sss, n_interp-n_bond*2+3)   = q_bond_mat(sss,4)
            results_mat(sss, n_interp-n_bond*2+4)   = q_bond_mat(sss,5)
            results_mat(sss, n_interp-n_bond+bbb)   = who_prices_bond

        enddo

        !$OMP END DO 
        !$OMP END PARALLEL


    enddo

    ! ------------------------------------------------------------ !
    ! Done with basic model solution here, store objects 
    ! ------------------------------------------------------------ !
    write(*,*) 'DONE WITH MODEL SOLUTION - NOW STORE'
    flush(6)

    ! print all results
    open (unit = 10, file = trim(results_folder) // 'next_state_mat.dat', ACTION="write", STATUS="replace", &
            & FORM="unformatted", ACCESS="STREAM")
    write(10) next_state_mat; close(10)

    open (unit = 10, file = trim(results_folder) // 'next_state_mat.txt', ACTION="write", STATUS="replace", &
            & FORM="formatted", ACCESS="STREAM")
    write(10,*) next_state_mat; close(10)

    ! print all results
    open (unit = 10, file = trim(results_folder) // 'results_mat.dat', ACTION="write", STATUS="replace", &
            & FORM="unformatted", ACCESS="STREAM")
    write(10) results_mat; close(10)

    ! print all results
    open (unit = 10, file = trim(results_folder) // 'results_mat.txt', ACTION="write", STATUS="replace", &
            & FORM="formatted", ACCESS="STREAM")
    write(10,*) results_mat; close(10)

    ! ==================================================================================================
    ! IRF SHOCKS
    ! ==================================================================================================
    do fff = 1,size(irf_indices,1)
        if (irf_indices(fff) == sidx_m) then
            monetary_shock = irf_shock_sizes(irf_indices(fff))
        else
            monetary_shock = 0.0_dp
        endif

        next_state_monetary_mat     = next_state_mat
        next_state_monetary_mat_new = next_state_mat
        k_next_monetary_mat         = k_next_mat
        nom_i_monetary_mat          = nom_i_mat 
        q_monetary_mat              = q_mat
        c_monetary_mat              = c_mat 
        share_monetary_mat          = share_mat 
        infl_monetary_mat           = infl_mat 
        l_aggr_monetary_mat         = l_aggr_mat 
        q_l_monetary_mat            = q_l_mat 
        w_choice_monetary_mat       = w_choice_mat 
        v_monetary_mat              = v_mat 
        mc_monetary_mat             = mc_mat 
        constraint_binding_monetary_mat = constraint_binding_mat 

        do iii = 1,n_I
            interp_input_mat(:,iii)       = v_mat(iii,:)
            interp_input_mat(:,iii + n_I) = mc_mat(iii,:)
        enddo
        interp_input_mat(:,2*n_I + 1)  = q_mat
        interp_input_mat(:,2*n_I + 2)  = l_aggr_mat
        interp_input_mat(:,2*n_I + 3)  = infl_mat
        interp_input_mat(:,2*n_I + 4)  = q_l_mat(1,:)
        interp_input_mat(:,2*n_I + 5)  = q_l_mat(2,:)
        interp_input_mat(:,2*n_I + 6)  = q_l_mat(3,:)

        ! Solve for Smolyak Coefficients (smol_coeffs)
        ! smol_coeffs has dimension n_states * 9
        call F07ABF( 'Equilibration','No transpose',n_states, nrhs, smol_polynom, lda,af,ldaf,ipiv,    &
                     equed,r,c,interp_input_mat,ldb,smol_coeffs,ldx,rcond,ferr,berr,work,iwork,info )

        outer_iter = 0
        diff = 1.0_dp
        do while (diff > conv_crit .and. outer_iter < maxiter) 
            outer_iter = outer_iter + 1

            !$OMP PARALLEL DEFAULT(NONE) &
            !$OMP SHARED( & 
            !$OMP next_state_monetary_mat, smol_elem_ani, smol_coeffs, state_grid, &
            !$OMP nom_i_monetary_mat, w_choice_monetary_mat, n_states, dz_vec, & 
            !$OMP k_grid_dev, k_grid_mean, q_monetary_mat, rho_m,  next_m_mat, m_grid_dev, &
            !$OMP c_monetary_mat, share_monetary_mat, infl_monetary_mat, m_grid_mean, & 
            !$OMP l_aggr_monetary_mat, next_state_monetary_mat_new, q_l_monetary_mat, &
            !$OMP w_grid_mean, w_grid_dev, n_active_dims, constraint_binding_monetary_mat, &
            !$OMP s_a_grid_mean, s_a_grid_dev, s_c_grid_mean, s_c_grid_dev, mc_monetary_mat_new, &
            !$OMP w_choice_monetary_mat_new, infl_monetary_mat_new, k_next_monetary_mat, &
            !$OMP q_monetary_mat_new, l_aggr_monetary_mat_new, v_monetary_mat_new, nom_i_monetary_mat_new, &
            !$OMP c_monetary_mat_new, k_next_monetary_mat_new, dz_vec_adj, q_l_monetary_mat_new, &
            !$OMP share_monetary_mat_new, monetary_shock, vector_mus_dimensions) &
            !$OMP PRIVATE( &
            !$OMP polyn_points, nxt_mat, current_state_vec, nom_i, mc_new, w_choice, nxt_mat_2, counter, &
            !$OMP q_current, c_vec, share_vec, q_l_new, infl , l_aggr, q_l_current, constraint_binding,  &
            !$OMP q_new, k_next_new, c_vec_new, infl_new, sss, v_new, l_aggr_new, w_choice_new, s_nxt)
            !$OMP DO SCHEDULE(static)
            do sss = 1, n_states
                ! ============================================================ !
                ! [1] Interpolate next period values from Smolyak Coefficients
                ! ============================================================ !
                ! Output is nxt_mat (n_quad * 9)
                ! ============================================================ !
                
                ! Get the Smolyak polynomials for potential next states
                ! next_state_mat(:,:,sss) has dimension n_quad * n_active_dims
                ! Its (i,j) element is the next period value of state j if today's state is sss and quadrature point i realizes
                ! polyn_points is the evaluated polynomial matrix, which has dimension n_quad * n_states 
                polyn_points = Smolyak_Polynomial2(next_state_monetary_mat(:,:,sss), n_active_dims, &
                                                   n_quad, n_states, max_smol_level, smol_elem_ani)

                ! Pre-Multiply polyn_points to smol_coeffs 
                ! The product is assigned to nxt_mat 
                CALL DGEMM('N','N', n_quad, nrhs, n_states, 1.0_dp, polyn_points, n_quad, & 
                           smol_coeffs, n_states, 0.0_dp, nxt_mat, n_quad)   

                ! k doesn't need interpolation
                ! k_next_mat is defined above 
                nxt_mat_2 = k_next_monetary_mat(:,sss)

                ! ============================================================ !
                ! [2] Calculate Equilibrium and Update
                ! ============================================================ !
                ! Load current states 
                nom_i             = nom_i_monetary_mat(sss)
                q_current         = q_monetary_mat(sss)
                c_vec             = c_monetary_mat(:,sss)
                share_vec         = share_monetary_mat(:,sss)
                infl              = infl_monetary_mat(sss)
                l_aggr            = l_aggr_monetary_mat(sss)
                q_l_current       = q_l_monetary_mat(:,sss)
                constraint_binding = constraint_binding_monetary_mat(:,sss) 

                ! Find updated variables
                call calc_equilibrium_and_update(nxt_mat, n_nxt, nxt_mat_2, sss, &
                                                 c_vec, l_aggr, q_current, q_l_current, &
                                                 infl, nom_i, share_vec, monetary_shock, &
                                                 q_new, q_l_new, k_next_new, c_vec_new, &
                                                 infl_new, v_new, mc_new, l_aggr_new, w_choice_new, s_nxt, constraint_binding)


                ! Update 
                v_monetary_mat_new(:,sss)           = v_new 
                mc_monetary_mat_new(:,sss)          = mc_new 
                c_monetary_mat_new(:,sss)           = c_vec_new
                k_next_monetary_mat_new(:,sss)      = k_next_new/exp(dz_vec_adj)*exp(dz_vec)
                q_monetary_mat_new(sss)             = q_new 
                l_aggr_monetary_mat_new(sss)        = l_aggr_new
                nom_i_monetary_mat_new(sss)         = nom_i
                w_choice_monetary_mat_new(sss)      = w_choice_new
                q_l_monetary_mat_new(:,sss)         = q_l_new
                share_monetary_mat_new(:,sss)       = share_vec          
                infl_monetary_mat_new(sss)          = infl_new
                constraint_binding_monetary_mat(:,sss) = constraint_binding

                ! ============================================================ !
                ! [3] Update next period state (next_state_mat_new)
                ! Note that it is on the Smolyak grid 
                ! ============================================================ !
                counter = 0

                ! Example: maps k_next_new calculated above to Smolyak grid 
                if (vector_mus_dimensions(1) > 0) then
                    counter = counter + 1
                    next_state_monetary_mat_new(:, counter, sss)   = (k_next_new/exp(dz_vec_adj) -k_grid_mean )/k_grid_dev
                endif

                if (vector_mus_dimensions(2) > 0) then
                    counter = counter + 1
                    next_state_monetary_mat_new(:, counter, sss) = (s_nxt(:,1) - s_a_grid_mean)/s_a_grid_dev
                endif

                if (vector_mus_dimensions(3) > 0) then
                    counter = counter + 1
                    next_state_monetary_mat_new(:, counter, sss) = (s_nxt(:,3) - s_c_grid_mean)/s_c_grid_dev
                endif

                if (vector_mus_dimensions(4) > 0) then
                    counter = counter + 1
                    next_state_monetary_mat_new(:,counter,sss) = (next_m_mat(:,sss) + rho_m*monetary_shock   - m_grid_mean) /m_grid_dev
                endif

                if (vector_mus_dimensions(5) > 0) then
                    counter = counter + 1
                    next_state_monetary_mat_new(:, counter, sss)  = (w_choice_new/exp(dz_vec_adj) -w_grid_mean) /w_grid_dev
                endif

                if (vector_mus_dimensions(6) > 0) then
                    counter = counter + 1
                endif

            enddo

            !$OMP END DO 
            !$OMP END PARALLEL

            where(next_state_monetary_mat_new > 1.0_dp) 
                next_state_monetary_mat_new   = 1.0_dp 
            endwhere
            where(next_state_monetary_mat_new < -1.0_dp)
                next_state_monetary_mat_new   = -1.0_dp
            endwhere

            ! check convergence
            if (outer_iter < 10) then
            diff = 1.0_dp
            else
            ! check convergence
            diff =  maxval([ maxval(abs(log(v_monetary_mat_new/v_monetary_mat))), &
                 maxval(abs(log(mc_monetary_mat_new/mc_monetary_mat))), &
                 maxval(abs((q_l_monetary_mat_new-q_l_monetary_mat))), &
                 maxval(abs(log(c_monetary_mat_new/c_monetary_mat))), &
                 maxval(abs((share_monetary_mat_new-share_monetary_mat))), &
                 maxval(abs(log(k_next_monetary_mat_new/k_next_monetary_mat))), &
                 maxval(abs(log(l_aggr_monetary_mat_new/l_aggr_monetary_mat))), &
                 maxval(abs(log(nom_i_monetary_mat_new/nom_i_monetary_mat))), &
                 maxval(abs(log(w_choice_monetary_mat_new/w_choice_monetary_mat))), &
                 maxval(abs(log(q_monetary_mat_new/q_monetary_mat))), &
                 maxval(abs(log(infl_monetary_mat_new/infl_monetary_mat))), &
                 maxval(abs(next_state_monetary_mat - next_state_monetary_mat_new)) ])
        
            endif

            timer_end = omp_get_wtime()

            write(*,*) ''
            write(*,"(A11, i5, i5)") ' Monetary shock iteration ', outer_iter
            write(*,"(A11, f8.2)") ' Calc. time', timer_end-timer_start
            write(*,"(A19,e12.4,A1)") 'Changes (max diff: ', diff , ')'
            write(*,*)

            ! Damping updates
            v_monetary_mat              = v_monetary_mat + 1.0*(v_monetary_mat_new-v_monetary_mat )
            mc_monetary_mat             = mc_monetary_mat + 1.0*(mc_monetary_mat_new-mc_monetary_mat )
            where (abs((next_state_monetary_mat-next_state_monetary_mat_new)) > 0.05)
                next_state_monetary_mat = next_state_monetary_mat + 0.01_dp*sign(1.0_dp, next_state_monetary_mat_new - next_state_monetary_mat)
            elsewhere
                next_state_monetary_mat = next_state_monetary_mat +  0.2_dp*(next_state_monetary_mat_new-next_state_monetary_mat)
            endwhere
            q_monetary_mat          = q_monetary_mat        + 0.2_dp*(q_monetary_mat_new-q_monetary_mat)
            l_aggr_monetary_mat     = l_aggr_monetary_mat   + 0.2_dp*(l_aggr_monetary_mat_new-l_aggr_monetary_mat)     
            w_choice_monetary_mat   = w_choice_monetary_mat*(1 + 0.2_dp*log(w_choice_monetary_mat_new/w_choice_monetary_mat))
            q_l_monetary_mat        = q_l_monetary_mat      + 0.5_dp*(q_l_monetary_mat_new-q_l_monetary_mat)
            c_monetary_mat          = c_monetary_mat        + 0.5*(c_monetary_mat_new - c_monetary_mat)
            nom_i_monetary_mat      = nom_i_monetary_mat    + 0.2*(nom_i_monetary_mat_new - nom_i_monetary_mat)
            infl_monetary_mat       =  infl_monetary_mat    + 0.2*(infl_monetary_mat_new - infl_monetary_mat)
            k_next_monetary_mat     = k_next_monetary_mat   + 0.2*(k_next_monetary_mat_new - k_next_monetary_mat)
            share_monetary_mat      = share_monetary_mat    + 1.0*(share_monetary_mat_new - share_monetary_mat) 

        enddo

        ! Obtain Results 
        q_bond_monetary_mat = 1.0
        !DIR$ NOUNROLL
        do bbb = 1,n_bond
            call F07ABF( 'Equilibration','No transpose',n_states, nrhs2, smol_polynom, lda,af,ldaf,ipiv,    &
                      equed,r,c, q_bond_mat(:,bbb),ldb,bond_coeffs,ldx,rcond,ferr2,berr2,work,iwork,info)
        
            !$OMP PARALLEL DEFAULT(NONE) &
            !$OMP SHARED( & 
            !$OMP next_state_monetary_mat, smol_elem_ani, smol_coeffs, &
            !$OMP nom_i_monetary_mat, n_states, & 
            !$OMP q_monetary_mat, k_next_monetary_mat, &
            !$OMP c_monetary_mat, share_monetary_mat, infl_monetary_mat, & 
            !$OMP l_aggr_monetary_mat, q_l_monetary_mat, &
            !$OMP n_active_dims, constraint_binding_monetary_mat, &
            !$OMP results_monetary_mat, bbb, bond_coeffs, q_bond_monetary_mat) &
            !$OMP PRIVATE( &
            !$OMP polyn_points, nxt_mat, current_state_vec, &
            !$OMP nom_i, w_choice, nxt_mat_2, counter, &
            !$OMP q_current, c_vec, share_vec, q_l_new, &
            !$OMP infl , l_aggr, q_l_current, constraint_binding, &
            !$OMP q_new, k_next_new, c_vec_new, infl_new, &
            !$OMP sss, v_new, l_aggr_new, w_choice_new, s_nxt, results_vec, &
            !$OMP who_prices_bond, nxt_bond_prices)
            !$OMP DO SCHEDULE(static)
            ! Loop over states
            do sss = 1, n_states

                ! Interpolation 
                polyn_points = Smolyak_Polynomial2(next_state_monetary_mat(:,:,sss), n_active_dims, &
                                                   n_quad, n_states, max_smol_level, smol_elem_ani)
                !! Standard variables
                CALL DGEMM('N','N', n_quad, nrhs, n_states, 1.0_dp, polyn_points, n_quad, & 
                               smol_coeffs, n_states, 0.0_dp, & 
                               nxt_mat, n_quad)   
                nxt_mat_2 = k_next_monetary_mat(:,sss) 
                !! Bond Prices 
                CALL DGEMM('N','N', n_quad, nrhs2, n_states, 1.0_dp, polyn_points, n_quad, & 
                                   bond_coeffs, n_states, 0.0_dp, & 
                                   nxt_bond_prices, n_quad)   

                ! Current State 
                nom_i = nom_i_monetary_mat(sss)
                q_current = q_monetary_mat(sss)
                c_vec = c_monetary_mat(:,sss)
                share_vec = share_monetary_mat(:,sss)
                infl = infl_monetary_mat(sss)
                l_aggr = l_aggr_monetary_mat(sss)
                q_l_current = q_l_monetary_mat(:,sss)
                constraint_binding = constraint_binding_monetary_mat(:,sss)

                ! Calculate Bond Prices 
                call calc_bond_prices(nxt_mat, n_nxt, nxt_mat_2, nxt_bond_prices, sss, &
                                c_vec, l_aggr, q_current, q_l_current, infl, &
                                nom_i, share_vec, constraint_binding, who_prices_bond, q_bond_monetary_mat(sss,bbb+1), results_vec)

                
                ! Store in results matrix 
                results_monetary_mat(sss, 1:n_interp-n_bond*2)                     = results_vec
                results_monetary_mat(sss, n_interp-n_bond*2+1)                     = q_bond_monetary_mat(sss,2)
                results_monetary_mat(sss, n_interp-n_bond*2+2)                     = q_bond_monetary_mat(sss,3)
                results_monetary_mat(sss, n_interp-n_bond*2+3)                     = q_bond_monetary_mat(sss,4)
                results_monetary_mat(sss, n_interp-n_bond*2+4)                     = q_bond_monetary_mat(sss,5)
                results_monetary_mat(sss, n_interp-n_bond+bbb)                     = who_prices_bond

            enddo

            !$OMP END DO 
            !$OMP END PARALLEL
        enddo

        trans_shock_vec = 0.0_dp
        if (irf_indices(fff) == sidx_z) then 
            trans_shock_vec(irf_indices(fff)) = irf_shock_sizes(fff)
        elseif (irf_indices(fff) == sidx_dis) then 
            trans_shock_vec(irf_indices(fff)) = irf_shock_sizes(fff)
        elseif (irf_indices(fff) == sidx_m) then 
            trans_shock_vec(irf_indices(fff)) = 0*irf_shock_sizes(fff) ! not this shock, would be double counting
        endif

        next_state_montransition_mat     = next_state_mat(no_shock_idx,:,:)
        next_state_montransition_mat_new = next_state_mat(no_shock_idx,:,:)

        do iii = 1,n_I
            interp_input_mat(:,iii)       = v_monetary_mat(iii,:)
            interp_input_mat(:,iii + n_I) = mc_monetary_mat(iii,:)
        enddo
        interp_input_mat(:,2*n_I + 1)  = q_monetary_mat
        interp_input_mat(:,2*n_I + 2)  = l_aggr_monetary_mat
        interp_input_mat(:,2*n_I + 3)  = infl_monetary_mat
        interp_input_mat(:,2*n_I + 4)  = q_l_monetary_mat(1,:)
        interp_input_mat(:,2*n_I + 5)  = q_l_monetary_mat(2,:)
        interp_input_mat(:,2*n_I + 6)  = q_l_monetary_mat(3,:)

        ! solve for polynomial coefficients
        ! call omp_set_num_threads(4)
        call F07ABF( 'Equilibration','No transpose',n_states, nrhs, smol_polynom, lda,af,ldaf,ipiv,    &
                      equed,r,c,interp_input_mat,ldb,smol_coeffs,ldx,rcond,ferr,berr,work,iwork,info)
        
        diff = 1.0_dp
        outer_iter = 0
        do while (diff > conv_crit .and. outer_iter < maxiter)
            outer_iter = outer_iter + 1

            ! ************************************************************ !
            !$OMP PARALLEL DEFAULT(NONE) &
            !$OMP SHARED( & 
            !$OMP fff, next_state_montransition_mat, smol_elem_ani, smol_coeffs, state_grid, &
            !$OMP nom_i_mat, w_choice_mat, n_states, dz_vec, & 
            !$OMP k_grid_dev, k_grid_mean, q_mat, &
            !$OMP c_mat, share_mat, infl_mat, next_state_mat, & 
            !$OMP l_aggr_mat, next_state_mat_new, q_l_mat, &
            !$OMP w_grid_mean, w_grid_dev, n_active_dims, &
            !$OMP s_a_grid_mean, s_a_grid_dev, s_c_grid_mean, s_c_grid_dev, &
            !$OMP w_choice_mat_new, infl_mat_new, &
            !$OMP q_mat_new, l_aggr_mat_new, v_mat_new, nom_i_mat_new, &
            !$OMP c_mat_new, k_next_mat_new, dz_vec_adj, q_l_mat_new, &
            !$OMP k_next_mat, share_mat_new, next_state_montransition_mat_new, &
            !$OMP next_m_mat, m_grid_mean, m_grid_dev, &
            !$OMP trans_shock_vec, dis_grid_dev, next_dis_mat, no_shock_idx, dis_grid_mean, vector_mus_dimensions) &
            !$OMP PRIVATE( &
            !$OMP polyn_points_trans, nxt_mat_trans, current_state_vec, &
            !$OMP nom_i, w_choice, nxt_mat_2_trans, counter, &
            !$OMP q_current, c_vec, share_vec, q_l_new, &
            !$OMP infl , l_aggr, q_l_current,  &
            !$OMP q_new, k_next_new, c_vec_new, infl_new, &
            !$OMP sss, v_new, l_aggr_new, w_choice_new, s_nxt)
            !$OMP DO SCHEDULE(static)
            do sss = 1, n_states
                ! ============================================================ !
                ! [1] Interpolate next period values from Smolyak Coefficients
                ! ============================================================ !
                ! Output is nxt_mat (n_quad * 9)
                ! ============================================================ !
                
                ! Get the Smolyak polynomials for potential next states
                ! next_state_mat(:,:,sss) has dimension n_quad * n_active_dims
                ! Its (i,j) element is the next period value of state j if today's state is sss and quadrature point i realizes
                ! polyn_points is the evaluated polynomial matrix, which has dimension n_quad * n_states 
                polyn_points_trans = Smolyak_Polynomial2(next_state_montransition_mat(:,:,sss), n_active_dims, &
                                                   1, n_states, max_smol_level, smol_elem_ani)

                ! Pre-Multiply polyn_points to smol_coeffs 
                ! The product is assigned to nxt_mat 
                CALL DGEMM('N','N', 1, nrhs, n_states, 1.0_dp, polyn_points_trans, 1, & 
                           smol_coeffs, n_states, 0.0_dp, nxt_mat_trans, 1)   

                ! k doesn't need interpolation
                ! k_next_mat is defined above 
                nxt_mat_2_trans = k_next_mat(1,sss)

                ! ============================================================ !
                ! [2] Calculate Equilibrium and Update
                ! ============================================================ !
                ! Load current states 
                nom_i             = nom_i_mat(sss)
                q_current         = q_mat(sss)
                c_vec             = c_mat(:,sss)
                share_vec         = share_mat(:,sss)
                infl              = infl_mat(sss)
                l_aggr            = l_aggr_mat(sss)
                q_l_current       = q_l_mat(:,sss)

                ! Find updated variables
                call calc_unexpected_transition(nxt_mat_trans, n_nxt, nxt_mat_2_trans, sss, &
                                                 c_vec, l_aggr, q_current, q_l_current, &
                                                 infl, nom_i, share_vec, &
                                                 trans_shock_vec, s_nxt(1,:))


                w_choice_new = w_choice_mat(sss)

                counter = 0

                if (vector_mus_dimensions(1) > 0) then
                    counter = counter + 1
                    next_state_montransition_mat_new(1, counter, sss)   = (k_next_mat(1,sss)/exp(trans_shock_vec(sidx_z)) -k_grid_mean )/k_grid_dev
                endif

                if (vector_mus_dimensions(2) > 0) then
                    counter = counter + 1
                    next_state_montransition_mat_new(:, counter, sss) = (s_nxt(1,1) - s_a_grid_mean)/s_a_grid_dev
                endif

                if (vector_mus_dimensions(3) > 0) then
                    counter = counter + 1
                    next_state_montransition_mat_new(:, counter, sss) = (s_nxt(1,3) - s_c_grid_mean)/s_c_grid_dev
                endif

                if (vector_mus_dimensions(4) > 0) then
                    counter = counter + 1
                    next_state_montransition_mat_new(:,counter,sss) = (next_m_mat(no_shock_idx,sss) + trans_shock_vec(sidx_m)   - m_grid_mean) /m_grid_dev
                endif

                if (vector_mus_dimensions(5) > 0) then
                    counter = counter + 1
                    next_state_montransition_mat_new(:, counter, sss)  = (w_choice_new/exp(trans_shock_vec(sidx_z)) -w_grid_mean) /w_grid_dev
                endif

                if (vector_mus_dimensions(6) > 0) then
                    counter = counter + 1
                    if (dis_grid_dev < sqrt_eps) then
                        next_state_montransition_mat_new(1,counter,sss) = 0.0_dp  + trans_shock_vec(sidx_dis)/dis_grid_dev
                    else
                        next_state_montransition_mat_new(:,counter,sss) = (next_dis_mat(no_shock_idx,sss) + trans_shock_vec(sidx_dis)   - dis_grid_mean) /dis_grid_dev
                    endif
                endif


            enddo
            !$OMP END DO
            !$OMP END PARALLEL
            

            where(next_state_montransition_mat_new> 1.0_dp) 
                next_state_montransition_mat_new= 1.0_dp 
            endwhere
            where(next_state_montransition_mat_new < -1.0_dp)
                next_state_montransition_mat_new   = -1.0_dp
            endwhere


            if (outer_iter < 10) then
                diff = 1.0_dp
            else
                diff =  maxval(abs(next_state_montransition_mat - next_state_montransition_mat_new)) 
            endif

            next_state_montransition_mat = next_state_montransition_mat +  0.2_dp*(next_state_montransition_mat_new-next_state_montransition_mat)

            

        enddo

        write(*,*) 'DONE WITH IRF CALCULATION - NOW STORE'
        flush(6)

        write(shock_char,'(i1)') fff
            
        ! print all results
        open (unit = 10, file = trim(results_folder) // 'next_state_shock_mat_' // shock_char // '.dat', ACTION="write", STATUS="replace", &
                & FORM="unformatted", ACCESS="STREAM")
        write(10) next_state_monetary_mat; close(10)

        ! print all results
        open (unit = 10, file = trim(results_folder) // 'next_state_shocktrans_mat_' // shock_char // '.dat', ACTION="write", STATUS="replace", &
                & FORM="unformatted", ACCESS="STREAM")
        write(10) next_state_montransition_mat; close(10)

        ! print all results
        open (unit = 10, file = trim(results_folder) // 'next_state_shocktrans_mat_' // shock_char // '.txt', ACTION="write", STATUS="replace", &
                & FORM="formatted", ACCESS="STREAM")
        write(10,*) next_state_montransition_mat; close(10)

        ! print all results
        open (unit = 10, file = trim(results_folder) // 'results_shock_mat_' // shock_char // '.dat', ACTION="write", STATUS="replace", &
                & FORM="unformatted", ACCESS="STREAM")
        write(10) results_monetary_mat; close(10)

        open (unit = 10, file = trim(results_folder) // 'results_shock_mat_' // shock_char // '.txt', ACTION="write", STATUS="replace", &
                & FORM="formatted", ACCESS="STREAM")
        write(10,*) results_monetary_mat; close(10)

    enddo


end subroutine calc_sol


! for impulse responses: find wealth transition in response to unexpected shock
subroutine calc_unexpected_transition( nxt_mat, n_nxt, nxt_mat_2, sss, c_vec, l_aggr, q_current, & 
                                       q_l_current, infl, nom_i, share_vec, trans_shock_vec, s_nxt)

     ! Load Parameters 
    use mod_param, only: n_I, n_quad, idx_k, aalpha, smolyak_d, quad_weight_vec, &
                         dz_vec, ddelta, state_grid, idx_dis, & 
                         chiX, xi, ies_vec, gma_vec, dz_vec_adj, bbeta_vec,      & 
                         thtbar_vec, tht, wealth_share_grid, &
                         lmbd_vec, s_bar_vec, vareps_w, tau_w, idx_w, & 
                         idx_m, phi, tayl_ic, s_trgt_vec, n_shocks, &
                         sidx_z, sidx_m, sidx_dis, labor_alloc_vec

    ! Declare input & output variables 
    integer, intent(in)  :: sss, n_nxt
    real(dp), intent(in) :: nxt_mat_2, l_aggr, nxt_mat(1, n_nxt), q_l_current(n_I), &
                            q_current, c_vec(n_I), infl, trans_shock_vec(n_shocks)
    real(dp), intent(in) :: share_vec(n_I), nom_i
    real(dp), intent(out)   :: s_nxt(1, 3)


    real(dp), parameter :: brent_delta = 5E-16_dp, min_cons_sav = 1E-08_dp
    real(dp) :: z_temp

    real(dp) :: k_aggr, p_dis, w_current, pi_current, mc_temp(1, n_I), &
                aggr_wealth, y_current, w_choice, k_nxt(1), next_k, v_temp(1, n_I), &
                q_nxt(1), l_aggr_nxt(1), infl_nxt(1), q_l_nxt(1, n_I), y_next(1), &
                w_next_choice(1), pi_nxt(1), rk_vec(1), rf_vec(1), c_temp(n_I), &
                savings, temp_nxt_aggr_wealth(1), survivor_wealth_vec(1, n_I), &
                savings_vec(n_I), ies, bbeta, gma, thtbar, wealth, labor, consumption, r_alpha_vec(1)


    integer :: iii, xxx, ggg, ifail, ccc, nnn, mflag

    ! =============================================================================
    ! CURRENT STATE ===============================================================
    ! =============================================================================
    k_aggr = state_grid(idx_k, sss)
    p_dis   = exp(state_grid(idx_dis, sss))                             
    w_current = state_grid(idx_w, sss)

    ! =============================================================================
    ! CURRENT PERIOD GUESSES ======================================================
    ! =============================================================================
    pi_current   = aalpha * (l_aggr/k_aggr)**(1.0-aalpha)    
    aggr_wealth  = k_aggr*((1-ddelta)*q_current + pi_current)   
    y_current    = k_aggr**aalpha * l_aggr**(1.0-aalpha)
    w_choice = (1.0-aalpha)*y_current/l_aggr

    ! =============================================================================
    ! NEXT PERIOD =================================================================
    ! =============================================================================
    z_temp = trans_shock_vec(sidx_z)

    k_nxt      = nxt_mat_2

    do iii = 1,n_I
        v_temp(:,iii)  = nxt_mat(:,      iii) 
        mc_temp(:,iii) = nxt_mat(:,n_I + iii) 
    enddo
    q_nxt        = nxt_mat(:, 2*n_I+1)
    l_aggr_nxt   = nxt_mat(:, 2*n_I+2)
    infl_nxt     = nxt_mat(:, 2*n_I+3)
    q_l_nxt(:,1) = nxt_mat(:, 2*n_I+4) * exp(z_temp)
    q_l_nxt(:,2) = nxt_mat(:, 2*n_I+5) * exp(z_temp)
    q_l_nxt(:,3) = nxt_mat(:, 2*n_I+6) * exp(z_temp)

    
    ! profit next period per unit of capital
    pi_nxt = aalpha*exp(z_temp*(1.0_dp-aalpha)) &
           *(l_aggr_nxt**(1.0_dp-aalpha))*k_nxt**(aalpha-1.0)
   
    ! capital return, adjustment in disaster state      
    rk_vec = ( (1.0_dp-ddelta)*q_nxt + pi_nxt)/q_current
    rf_vec = nom_i/infl_nxt 

    !! =============================================================================
    !! STORE VALUE FUNCTION, NEXT CAPITAL, NEXT WEALTH SHARES ======================
    !! =============================================================================
    c_temp    = c_vec
    do iii = 1,n_I
        ! Read parameters for agent 
        IES   = IES_vec(iii)
        bbeta = bbeta_vec(iii)
        gma   = gma_vec(iii)
        thtbar  = thtbar_vec(iii)

        wealth = aggr_wealth*wealth_share_grid(iii,sss) 
        labor  = l_aggr*labor_alloc_vec(iii)

        ! Enforce minimal spending and mininal saving 
        consumption =  minval( [maxval([c_temp(iii), min_cons_sav]), & 
                          wealth + w_choice*labor-min_cons_sav] )

        ! Get Next Period Wealth Share 
        r_alpha_vec = share_vec(iii)*rf_vec + (1.0_dp - share_vec(iii))*rk_vec
        savings = wealth + w_choice * labor - consumption

        ! update wealth transitions accordingly
        survivor_wealth_vec(:,iii) = lmbd_vec(iii)*savings*r_alpha_vec
        savings_vec(iii) = savings
    enddo

    ! (4.2) Next period wealth share 
    temp_nxt_aggr_wealth = sum(survivor_wealth_vec,2)
    ! write(*,"(A22, F10.4)") 'temp_nxt_aggr_wealth: ', temp_nxt_aggr_wealth
    do iii = 1,n_I
            ! Next period wealth, accounting for survivors and newborns
            ! EXPORTED
            s_nxt(:,iii) = ((1.0-xi)*survivor_wealth_vec(:,iii) + xi*s_bar_vec(iii)*temp_nxt_aggr_wealth) & 
                                  /(temp_nxt_aggr_wealth)
    enddo

end subroutine calc_unexpected_transition


! find equilbrium and polies for current state
subroutine calc_equilibrium_and_update(nxt_mat, n_nxt, nxt_mat_2, sss, &
                                       c_vec, l_aggr, q_current, q_l_current, &
                                       infl, nom_i, share_vec, monetary_shock, &
                                       q_new, q_l_new, k_next_new, c_vec_new, &
                                       infl_new, v_new, mc_new, l_aggr_new, w_choice_new, s_nxt, constraint_binding)
    ! Load Parameters 
    use mod_param, only: n_I, n_quad, idx_k, aalpha, smolyak_d, quad_weight_vec, &
                         dz_vec, ddelta, state_grid, idx_dis, &
                         chiX, xi, ies_vec, gma_vec, dz_vec_adj, &
                         bbeta_vec, thtbar_vec, tht, wealth_share_grid, &
                         lmbd_vec, s_bar_vec, vareps_w, &
                         tau_w, chiW, idx_w, idx_m, phi, tayl_ic, & 
                         s_trgt_vec, constr_agents, &
                         labor_alloc_vec, use_idio_risk, idio_weight_vec, idio_shock_grid, n_idio

    ! Declare input & output variables 
    integer, intent(in)  :: sss, n_nxt
    real(dp), intent(in) :: nxt_mat_2(n_quad), l_aggr, monetary_shock, &
                            nxt_mat(n_quad, n_nxt), q_l_current(n_I), &
                            q_current, c_vec(n_I), infl
    real(dp), intent(inout) :: share_vec(n_I), nom_i
    real(dp), intent(out)   :: q_new, q_l_new(n_I), k_next_new, c_vec_new(n_I), infl_new, &
                               v_new(n_I), mc_new(n_I), l_aggr_new, w_choice_new, s_nxt(n_quad, 3)
    integer, intent(inout) :: constraint_binding(n_I)

    real(dp), parameter :: brent_delta = 5E-16_dp, min_cons_sav = 1E-08_dp

    real(dp) :: k_aggr, p_dis, w_current, big_weight_vec(n_quad), pi_current, &
                aggr_wealth, y_current, w_choice, k_nxt(n_quad), next_k, v_temp(n_quad, n_I), &
                q_nxt(n_quad), l_aggr_nxt(n_quad), infl_nxt(n_quad), q_l_nxt(n_quad, n_I), &
                w_next_choice(n_quad), pi_nxt(n_quad), rk_vec(n_quad), rf_vec(n_quad), &
                c_temp(n_I), rf_vec_temp(n_quad), share_temp(n_I), k_temp_vec(n_I),& 
                b_temp_vec(n_I), excess_b, nom_i_in, nom_i_guess, nom_i_A, &
                excess_b_A, nom_i_new, excess_b_S, nom_i_B, excess_b_B,  mc_temp(n_quad, n_I), &
                temp, nom_i_C, excess_b_C, nom_i_S, nom_i_D, savings, EV, objf, &
                temp_nxt_aggr_wealth(n_quad), survivor_wealth_vec(n_quad, n_I), &
                next_period_share(n_quad), labor_frac_nxt(n_quad, n_I), &
                diff_inner, labor_part, deriv_helper(n_quad), &
                cons_update, labor_new, l_temp_vec(n_I), savings_vec(n_I), &
                ies, bbeta, gma, thtbar, wealth, labor, &
                consumption, r_alpha_vec(n_quad), v_vec_twisted(n_quad), util, &
                util_c_deriv, M_vec_mat(n_quad, n_I), marg_util_c_vec(n_I), &
                weight_vec(n_I), l_aggr_temp, w_choice_temp, &
                w_temp, w_diff, w_current_next(n_quad), w_temp_new, y_next(n_quad)


    integer :: iii, xxx, ggg, ifail, ccc, nnn, kkk
    integer :: iter, mflag, iter_inner, constraint_binding_new(n_I)

    real(dp) ::v_idio(n_quad*n_idio, n_I), mc_idio(n_quad*n_idio, n_I), q_l_nxt_idio(n_quad*n_idio, n_I), &
            dz_vec_idio(n_quad*n_idio), big_weight_vec_idio(n_quad*n_idio), &
            r_alpha_vec_idio(n_quad*n_idio), next_period_share_idio(n_quad*n_idio), &
            v_vec_twisted_idio(n_quad*n_idio), rf_vec_idio(n_quad*n_idio), &
            deriv_helper_idio(n_quad*n_idio), M_idio(n_quad*n_idio), M_idio_mat(n_quad*n_idio, n_I)


    ! ============================================================ !
    ! (1.1) Load Current State 
    ! ============================================================ !
    k_aggr = state_grid(idx_k, sss) ! Current Capital Stock 
    p_dis   = exp(state_grid(idx_dis, sss)) ! Current disaster prob
    w_current = state_grid(idx_w, sss) ! Current wage 

    ! Weighting vector across n_quad quadrature points, account for disaster 
    big_weight_vec(1:n_quad-1) = quad_weight_vec*(1.0-p_dis)
    big_weight_vec(n_quad)     = p_dis

    ! ============================================================ !
    ! (1.2) Current Guesses 
    ! ============================================================ !
    ! Calculate profit \pi, aggrgate wealth, output, wage
    pi_current = aalpha * (l_aggr/k_aggr)**(1.0-aalpha)
    aggr_wealth = k_aggr * (pi_current + (1.0-ddelta)*q_current)
    y_current = k_aggr**aalpha * l_aggr**(1.0-aalpha)
    w_choice = (1.0-aalpha)*y_current/l_aggr 



    ! =============================================================================
    ! (2) Basic Next-Period Updates 
    ! =============================================================================
    ! Load next-period values: those all have dimensions n_quad 
    k_nxt = nxt_mat_2 
    do iii = 1,n_I
        v_temp(:,iii)  = nxt_mat(:,      iii) 
        mc_temp(:,iii) = nxt_mat(:,n_I + iii) 
    enddo
    q_nxt        = nxt_mat(:, 2*n_I+1)
    l_aggr_nxt   = nxt_mat(:, 2*n_I+2)
    infl_nxt     = nxt_mat(:, 2*n_I+3)
    q_l_nxt(:,1) = nxt_mat(:, 2*n_I+4) * exp(dz_vec)
    q_l_nxt(:,2) = nxt_mat(:, 2*n_I+5) * exp(dz_vec)
    q_l_nxt(:,3) = nxt_mat(:, 2*n_I+6) * exp(dz_vec)


    ! Price of capital (EXPORTED)
    q_new = (k_nxt(1)/k_aggr)**chiX

    ! Output
    y_next = exp((dz_vec)*(1.0_dp-aalpha))*l_aggr_nxt* & 
                  ( (k_nxt/l_aggr_nxt)**aalpha)

    ! Wage
    w_next_choice = (1-aalpha)*y_next/l_aggr_nxt
    
    ! Profit
    pi_nxt = aalpha*exp( (dz_vec)*(1.0_dp-aalpha)) *(l_aggr_nxt**(1.0_dp-aalpha))*k_nxt**(aalpha-1.0)

    ! Gross Capital Return (Adjust in disaster state)
    rk_vec = ( (1.0_dp-ddelta)*q_nxt + pi_nxt ) / q_current 
    rk_vec(n_quad) = exp(dz_vec(n_quad)) * rk_vec(n_quad) ! adjust for disaster depreciation

    ! Gross Risk-Free Return
    rf_vec = nom_i / infl_nxt

    ! =============================================================================
    ! (3) Solve for Nominal Rate 
    ! =============================================================================
    ! at this point inflation and consumption choices are held fixed
    ! excess demand is solved for aggregate bond portfolios

    ! Current consumption and rate 
    c_temp    = c_vec
    rf_vec_temp = rf_vec
    share_temp = share_vec

    ! Get Excess Bond 
    call calc_excess_bond_nom(aggr_wealth, v_temp, mc_temp, rf_vec_temp, rk_vec,  &  
                              q_current, w_choice, l_aggr, k_nxt(1)*q_current, q_l_nxt, sss, &
                              big_weight_vec, c_temp, share_temp, &
                              k_temp_vec, b_temp_vec, excess_b, constraint_binding, constraint_binding_new)

    ! Initial Guess 
    nom_i_guess = nom_i
    nom_i_A    = nom_i_guess
    excess_b_A = excess_b

    ! If no excess bond demand, just return the guess 
    if (excess_b == 0) then 
        nom_i_new = nom_i_guess
        excess_b_S = excess_b
    ! If not......
    else
        ! (1) Positive excess bond 
        ! If initial rate guess is too high, reduce the rate until excess demand <= 0
        ! This rate is then the lower bond
        if (excess_b > 0.0_dp) then 
            iter = 0
            do while (excess_b > 0.0_dp .and. iter <= 5000)
                iter = iter + 1

                ! Reduce nominal rate 
                nom_i_guess = nom_i_guess - 5E-03_dp

                ! Calculate excess demand
                rf_vec_temp = nom_i_guess/infl_nxt
                share_temp = share_vec

                call calc_excess_bond_nom(aggr_wealth, v_temp, mc_temp, rf_vec_temp, rk_vec,  &  
                                          q_current, w_choice, l_aggr, k_nxt(1)*q_current, q_l_nxt, sss, &
                                          big_weight_vec, c_temp, share_temp, &
                                          k_temp_vec, b_temp_vec, excess_b, constraint_binding, constraint_binding_new)
            enddo

            if (iter >= 5000) then
                write(*,*) 'No convergence in finding nom_i_low'
                stop
            endif

            nom_i_B = nom_i_guess
            excess_b_B = excess_b

        ! (2) Negative excess demand
        ! If initial rate guess is too low, increase the rate until excess demand >= 0
        ! This rate is then the upper bond
        else
            iter = 0
            do while (excess_b < 0.0_dp .and. iter <= 5000)
                iter = iter + 1

                ! Reduce nominal rate 
                nom_i_guess = nom_i_guess + 5E-03_dp

                ! Calculate excess demand
                rf_vec_temp = nom_i_guess/infl_nxt
                share_temp = share_vec

                call calc_excess_bond_nom(aggr_wealth, v_temp, mc_temp, rf_vec_temp, rk_vec,  &  
                                          q_current, w_choice, l_aggr, k_nxt(1)*q_current, q_l_nxt, sss, &
                                          big_weight_vec, c_temp, share_temp, &
                                          k_temp_vec, b_temp_vec, excess_b, constraint_binding, constraint_binding_new)                
            enddo
            if (iter >= 5000) then
                write(*,*) 'No convergence in finding nom_i_high'
                stop
            endif

            nom_i_B = nom_i_guess
            excess_b_B = excess_b

        endif

        ! If the initial excess demand excess_b_A and the demand at the other end of bracket has the same sign 
        ! Then this signals problem: true rate is not in bracket 
        ! (Though given previous procedure, I'm not sure why this would ever occur?)
        if ((excess_b_A*excess_b_B) > 0) then
           write(*,*) 'ERROR: Initial bracket does not contain root.'
           stop
        endif

        ! B should be closer to 0 than A. If not, swap 
        if ( abs(excess_b_A) < abs(excess_b_B) )then
           temp = nom_i_A
           nom_i_A = nom_i_B
           nom_i_B = temp

           temp = excess_b_A
           excess_b_A = excess_b_B
           excess_b_B = temp
        endif

        ! ************************************************************ !
        ! (3) Do BRENT update 
        ! NOTE: Haven't looked into how this works yet
        ! ************************************************************ !
        nom_i_C    = nom_i_A
        excess_b_C = excess_b_A
        excess_b_S = excess_b_B
        mflag = 1

        ! Iterate until excess bond demand is zero:
        do while (excess_b_S /= 0 .and. abs(nom_i_A - nom_i_B) > 1E-15_dp)
            ! ************************************************************ !
            ! Obtain updated position
            ! ************************************************************ !
            if ( (excess_b_A /= excess_b_C) .and. (excess_b_C /= excess_b_B) ) then ! inverse quadratic interpolation

                nom_i_S = nom_i_A * excess_b_B * excess_b_C / ( (excess_b_A - excess_b_B) * (excess_b_A - excess_b_C) ) + &
                          nom_i_B * excess_b_A * excess_b_C / ( (excess_b_B - excess_b_A) * (excess_b_B - excess_b_C) ) + &
                          nom_i_C * excess_b_A * excess_b_B / ( (excess_b_C - excess_b_A) * (excess_b_C - excess_b_B) )

            else ! secant method

                nom_i_S = nom_i_B - excess_b_B * (nom_i_B - nom_i_A) /   (excess_b_B - excess_b_A)

            endif

            if ( ( (  ( nom_i_S > ( 3*nom_i_A + nom_i_B )/4 .and. nom_i_S < nom_i_B) .or. &
                        ( nom_i_S < ( 3*nom_i_A + nom_i_B )/4 .and. nom_i_S > nom_i_B)  ) == .FALSE. ) .or. &
                   (mflag == 1 .and. abs(nom_i_S - nom_i_B) >= abs(nom_i_B - nom_i_C)/2  )             .or. &
                   (mflag == 0 .and. abs(nom_i_S - nom_i_B) >= abs(nom_i_B - nom_i_D)/2  )             .or. &
                   (mflag == 1 .and. abs(nom_i_B - nom_i_C) <  abs(brent_delta)  )                     .or. &
                   (mflag == 0 .and. abs(nom_i_B - nom_i_D) <  abs(brent_delta)  )                          &
                 ) then

                nom_i_S = (nom_i_A + nom_i_B )/ 2
                mflag = 1
            else
                mflag = 0
            endif

            ! ************************************************************ !
            ! check whether it generates excess demand or not
            ! ************************************************************ !
            rf_vec_temp = nom_i_S/infl_nxt
            share_temp = share_vec
            call calc_excess_bond_nom(aggr_wealth, v_temp, mc_temp, rf_vec_temp, rk_vec,  &  
                                      q_current, w_choice, l_aggr, k_nxt(1)*q_current, q_l_nxt, sss, &
                                      big_weight_vec, c_temp, share_temp, &
                                      k_temp_vec, b_temp_vec, excess_b_S, constraint_binding, constraint_binding_new)

            ! ************************************************************ !
            ! Update Positions
            ! ************************************************************ !
            nom_i_D = nom_i_C
            nom_i_C = nom_i_B
            excess_b_C = excess_b_B

            if ((excess_b_A*excess_b_S) < 0) then
            nom_i_B = nom_i_S
            excess_b_B = excess_b_S
            else
            nom_i_A = nom_i_S
            excess_b_A = excess_b_S
            endif

            ! swap a and b
            if ( abs(excess_b_A) < abs(excess_b_B) )then
             temp    = nom_i_A
             nom_i_A = nom_i_B
             nom_i_B = temp

             temp       = excess_b_A
             excess_b_A = excess_b_B
             excess_b_B = temp
            endif

        enddo

        ! Found updated nominal rate that clears bond market 
        nom_i_guess = nom_i_S
        nom_i_new   = nom_i_guess

    endif
    
    ! Potential failure
    if (excess_b_S>1E-02_dp) then
        write(*,*) 'WARNING: EXCESS AGGREGATE BONDS HOLDINGS.'
        stop
    endif

    ! Store results of this block 
    ! nom_i is exported 
    share_vec = share_temp
    nom_i_in = nom_i
    nom_i = nom_i_guess 
    rf_vec = nom_i / infl_nxt


    ! Idiosyncratic Risk Expansion
    if (use_idio_risk == 1) then      
        do kkk = 1,n_idio
            do iii = 1,n_I
                v_idio(((kkk-1)*n_quad+1):(kkk*n_quad),iii)  = v_temp(:,iii)
                mc_idio(((kkk-1)*n_quad+1):(kkk*n_quad),iii) = mc_temp(:,iii)
                q_l_nxt_idio(((kkk-1)*n_quad+1):(kkk*n_quad),iii) = q_l_nxt(:,iii)*exp(idio_shock_grid(kkk,iii))
            enddo
            dz_vec_idio(((kkk-1)*n_quad+1):(kkk*n_quad)) = dz_vec
            rf_vec_idio(((kkk-1)*n_quad+1):(kkk*n_quad)) = rf_vec
            big_weight_vec_idio(((kkk-1)*n_quad+1):(kkk*n_quad)) = big_weight_vec * idio_weight_vec(kkk)
        enddo
    endif

    ! =============================================================================
    ! (4) Store value function, next-period capital, next-period wealth shares 
    ! =============================================================================
    c_temp    = c_vec
    do iii = 1,n_I
        ! Read parameters for agent 
        IES   = IES_vec(iii)
        bbeta = bbeta_vec(iii)
        gma   = gma_vec(iii)
        thtbar  = thtbar_vec(iii)

        wealth = aggr_wealth*wealth_share_grid(iii,sss) 
        labor  = l_aggr*labor_alloc_vec(iii)

        ! Enforce minimal spending and mininal saving 
        consumption =  minval( [maxval([c_temp(iii), min_cons_sav]), & 
                          wealth + w_choice*labor-min_cons_sav] )
        savings = wealth + w_choice * labor - consumption
        savings_vec(iii) = savings

        ! Portfolio return and wealth transitions
        r_alpha_vec = share_vec(iii)*rf_vec + (1.0_dp - share_vec(iii))*rk_vec
        survivor_wealth_vec(:,iii) = lmbd_vec(iii)*savings*(r_alpha_vec)
        
        if (use_idio_risk == 1) then
            ! ------------------------------------------------------------ !
            ! With Idiosyncratic Risk
            ! ------------------------------------------------------------ !
            ! Expand portfolio return 
            do kkk = 1,n_idio
                r_alpha_vec_idio(((kkk-1)*n_quad+1):(kkk*n_quad)) = share_vec(iii)*rf_vec + (1.0_dp - share_vec(iii))*rk_vec*exp(idio_shock_grid(kkk,iii))
            enddo

            ! Get Next Period Wealth Share 
            next_period_share_idio = (savings*r_alpha_vec_idio + q_l_nxt_idio(:,iii) )/tot_wealth_ss_vec(iii)
            v_vec_twisted_idio = (v_idio(:,iii)*next_period_share_idio) ** (1.0-gma)

            ! Calculate Expected Value Function
            EV = sum(big_weight_vec_idio*v_vec_twisted_idio)
        else
            ! ------------------------------------------------------------ !
            ! No Idiosyncratic Risk
            ! ------------------------------------------------------------ !
            ! Get Next Period Wealth Share 
            next_period_share = (savings*r_alpha_vec + q_l_nxt(:,iii) )/tot_wealth_ss_vec(iii)
            v_vec_twisted = (v_temp(:,iii)*next_period_share) ** (1.0-gma)

            ! Calculate Expected Value Function
            EV = sum(big_weight_vec*v_vec_twisted)
        endif

        ! Value Function Update
        call util_fun(consumption, labor, iii, util, util_c_deriv, labor_part)
        objf = ( (1.0-bbeta)*util +  bbeta*  (EV**( (1.0_dp - (1.0_dp/IES))/(1.0_dp-gma) )) )**(1.0_dp/(1.0_dp-1.0_dp/IES))
        v_new(iii) = objf*tot_wealth_ss_vec(iii)/(wealth + q_l_current(iii))
        mc_new(iii) = util_c_deriv * (tot_wealth_ss_vec(iii)/(wealth + q_l_current(iii)))**(-1.0_dp/IES)
    enddo

    ! (4.2) Next period wealth share 
    temp_nxt_aggr_wealth = sum(survivor_wealth_vec,2)
    do iii = 1,n_I
        ! Next period wealth, accounting for survivors and newborns
        s_nxt(:,iii) = ((1.0-xi)*survivor_wealth_vec(:,iii) + xi*s_bar_vec(iii)*temp_nxt_aggr_wealth) & 
                              /(temp_nxt_aggr_wealth)
        labor_frac_nxt(:,iii)   = ((1.0-xi)*( survivor_wealth_vec(:,iii) + lmbd_vec(iii)*q_l_nxt(:,iii) ))/ &
                           ((1.0-xi)*survivor_wealth_vec(:,iii) + xi*s_bar_vec(iii)*temp_nxt_aggr_wealth + lmbd_vec(iii)*q_l_nxt(:,iii) )
    enddo

    ! (4.3) Next period capital 
    ! EXPORTED
    k_next_new = sum(lmbd_vec*savings_vec*(1.0-share_vec)/q_current)

    ! =============================================================================
    ! (6) Update Consumption Choice
    ! =============================================================================
    c_temp = c_vec 
    do iii = 1,n_I

        ! (6.1) Read parameters and quantities (same as in part (4) above)
        IES   = IES_vec(iii)
        bbeta = bbeta_vec(iii)
        gma   = gma_vec(iii)
        thtbar  = thtbar_vec(iii)

        wealth = aggr_wealth*wealth_share_grid(iii,sss) 
        labor  = l_aggr*labor_alloc_vec(iii)

        consumption =  minval( [maxval([c_temp(iii), min_cons_sav]), & 
                          wealth + w_choice*labor-min_cons_sav] )
        r_alpha_vec = share_vec(iii)*rf_vec + (1.0_dp - share_vec(iii))*rk_vec
        labor_part =  1.0_dp + (1.0_dp/IES - 1.0_dp) * thtbar*tht/(1.0_dp + tht) * labor**( (1.0_dp+tht)/tht )

        ! (6.2) Loop and update consumption
        diff_inner = 1.0_dp 
        iter_inner = 0
        do while (diff_inner > 1E-08_dp)
            iter_inner = iter_inner + 1

            ! Based on equation (48/54) of equation doc 
            savings = wealth + w_choice*labor - consumption


            if (use_idio_risk == 1) then
                ! ------------------------------------------------------------ !
                ! With Idiosyncratic Risk 
                ! ------------------------------------------------------------ !
                ! Expand portfolio return 
                do kkk = 1,n_idio
                r_alpha_vec_idio(((kkk-1)*n_quad+1):(kkk*n_quad)) = share_vec(iii)*rf_vec + (1.0_dp - share_vec(iii))*rk_vec*exp(idio_shock_grid(kkk,iii))
                enddo

                ! Other Operations
                next_period_share_idio = (savings*r_alpha_vec_idio + q_l_nxt_idio(:,iii) )/tot_wealth_ss_vec(iii) 
                v_vec_twisted_idio = (v_idio(:,iii)*next_period_share_idio) ** (1.0-gma)
                EV = sum(big_weight_vec_idio*v_vec_twisted_idio)

                deriv_helper_idio = bbeta * (v_idio(:,iii)**((1.0_dp/IES) - gma )) * (next_period_share_idio**(-gma)) & 
                                  * (EV **( (gma - (1.0_dp/IES))/(1.0_dp-gma) )) * mc_idio(:,iii)

                if (constr_agents(iii) == 1 .and. constraint_binding(iii) == 1) then
                    temp = sum(deriv_helper_idio*rf_vec_idio*big_weight_vec_idio)
                else
                    temp = sum(deriv_helper_idio*r_alpha_vec_idio*big_weight_vec_idio)
                endif       

            else
                ! ------------------------------------------------------------ !
                ! No Idiosyncratic Risk
                ! ------------------------------------------------------------ !
                next_period_share = (savings*r_alpha_vec + q_l_nxt(:,iii) )/tot_wealth_ss_vec(iii) 
                v_vec_twisted = (v_temp(:,iii)*next_period_share) ** (1.0-gma)
                EV = sum(big_weight_vec*v_vec_twisted)

                deriv_helper = bbeta * (v_temp(:,iii)**((1.0_dp/IES) - gma )) * (next_period_share**(-gma))  & 
                               * (EV **( (gma - (1.0_dp/IES))/(1.0_dp-gma) )) * mc_temp(:,iii)
                if (constr_agents(iii) == 1 .and. constraint_binding(iii) == 1) then
                    temp = sum(deriv_helper*rf_vec*big_weight_vec)
                else
                    temp = sum(deriv_helper*r_alpha_vec*big_weight_vec)
                endif
            endif


            cons_update = labor_part * ( temp**(-IES) ) 

            ! Update 
            cons_update = minval([maxval([cons_update, min_cons_sav]), (wealth + w_choice*labor-min_cons_sav)])                
            diff_inner = abs(consumption-cons_update)                    
            consumption = consumption + 0.5*(cons_update-consumption)

            if (iter_inner > 1000) then
                if (diff_inner > 0.1_dp) then
                    write(*,*) 'NO CONVERGENCE 1.'
                    write(*,*) sss, iii 
                    stop
                else
                    write(*,*) 'WARNING: NO CONVERGENCE 1.'
                    diff_inner = 0.0_dp
                endif
            endif
        enddo
        ! EXPORT
        c_vec_new(iii) = consumption 
        savings_vec(iii) = savings 

    enddo

    ! =============================================================================
    ! (7) Update marginal utility and pricing kernel 
    ! =============================================================================
    M_vec_mat = 0.0_dp
    marg_util_c_vec = 0.0_dp

    do iii = 1,n_I
        ! (7.1) Read parameters and quantities (same as in part (4) above)
        IES   = IES_vec(iii)
        bbeta = bbeta_vec(iii)
        gma   = gma_vec(iii)
        thtbar  = thtbar_vec(iii)

        wealth = aggr_wealth*wealth_share_grid(iii,sss) 
        labor  = l_aggr*labor_alloc_vec(iii)

        consumption =  minval( [maxval([c_vec_new(iii), min_cons_sav]), & 
                          wealth + w_choice*labor-min_cons_sav] )
        r_alpha_vec = share_vec(iii)*rf_vec + (1.0_dp - share_vec(iii))*rk_vec
        savings     = wealth + w_choice*labor - consumption

        ! (7.2) Calculate (Realized) Pricing Kernel
        if (use_idio_risk == 1) then
            ! ------------------------------------------------------------ !
            ! With Idiosyncratic Risk
            ! ------------------------------------------------------------ !
            ! Expand portfolio return 
            do kkk = 1,n_idio
                r_alpha_vec_idio(((kkk-1)*n_quad+1):(kkk*n_quad)) = share_vec(iii)*rf_vec + (1.0_dp - share_vec(iii))*rk_vec*exp(idio_shock_grid(kkk,iii))
            enddo

            ! Other Operations
            next_period_share_idio = (savings*r_alpha_vec_idio + q_l_nxt_idio(:,iii) )/tot_wealth_ss_vec(iii)
            v_vec_twisted_idio = (v_idio(:,iii)*next_period_share_idio) ** (1.0-gma)
            EV = sum(big_weight_vec_idio*v_vec_twisted_idio)

            call util_fun(consumption, labor, iii, util, util_c_deriv, labor_part)

            objf = ( (1.0 - bbeta)*util +  bbeta*  (EV**( (1.0_dp - (1.0_dp/IES))/(1.0_dp-gma) )) )**(1.0_dp/(1.0_dp-1.0_dp/IES))

            marg_util_c_vec(iii) = objf**(1.0_dp/IES)*(1.0-bbeta)*util_c_deriv 

            ! SDF: first, the full form 
            M_idio = bbeta* ( next_period_share_idio ** (-gma) ) * (EV **( (gma - (1.0_dp/IES))/(1.0_dp-gma) ))  & 
                     * (v_idio(:,iii)**((1.0_dp/IES) -gma )) * mc_idio(:,iii) / util_c_deriv 

            M_idio_mat(:,iii) = M_idio
            ! Then for each of the n_quad aggregate nodes, average across the n_idio Idiosyncratic nodes
            M_vec_mat(:,iii)  = 0.0_dp 
            do kkk = 1,n_idio
                M_vec_mat(:,iii) = M_vec_mat(:,iii) + M_idio(((kkk-1)*n_quad+1):(kkk*n_quad)) * idio_weight_vec(kkk)
            enddo

        else
            ! ------------------------------------------------------------ !
            ! No Idiosyncratic Risk
            ! ------------------------------------------------------------ !
            next_period_share = (savings*r_alpha_vec + q_l_nxt(:,iii) )/tot_wealth_ss_vec(iii)
            v_vec_twisted = (v_temp(:,iii)*next_period_share)**(1.0-gma)
            EV = sum(big_weight_vec*v_vec_twisted)

            call util_fun(consumption, labor, iii, util, util_c_deriv, labor_part)

            objf = ( (1.0-bbeta)*util +  bbeta*  (EV**( (1.0_dp - (1.0_dp/IES))/(1.0_dp-gma) )) )**(1.0_dp/(1.0_dp-1.0_dp/IES))

            marg_util_c_vec(iii) = objf**(1.0_dp/IES)*(1.0-bbeta)*util_c_deriv 

            M_vec_mat(:,iii) = bbeta* ( next_period_share ** (-gma) ) * (EV **( (gma - (1.0_dp/IES))/(1.0_dp-gma) ))  & 
                     * (v_temp(:,iii)**((1.0_dp/IES) -gma )) * mc_temp(:,iii) / util_c_deriv 
        endif

        
    enddo

    if (any(isnan(M_vec_mat))) then
            write(*,*) 'M_vec_mat is NaN.'
            stop
    endif

    ! =============================================================================
    ! (8) Update Pricing of Time Endowment 
    ! =============================================================================
    do iii = 1,n_I 
        
        
        ! EXPORTED, CHANGED
        if (labor_alloc_vec(iii) > sqrt_eps) then  
         if (use_idio_risk == 1) then
         q_l_new(iii) =  w_choice*l_aggr*labor_alloc_vec(iii) + sum(big_weight_vec_idio*M_idio_mat(:,iii)*q_l_nxt_idio(:,iii))             
         else
         q_l_new(iii) =  l_aggr*labor_alloc_vec(iii)*w_choice +  sum(big_weight_vec*M_vec_mat(:,iii)*q_l_nxt(:,iii))             
         endif
        else
            q_l_new = 0.0_dp
        endif
    enddo

    ! =============================================================================
    ! (9) Update labor, wages
    ! =============================================================================
    ! Combined Weights
    weight_vec = lmbd_vec*marg_util_c_vec*labor_alloc_vec

    if (chiW >= -1.0_dp) then
        ! Wage rigidity
        l_aggr_temp   = l_aggr
        w_choice_temp = w_choice

        ! Find wage implied by Union optimization 
        w_temp = w_choice_temp
        w_diff = 1.0_dp 
        iter = 0

        do while ( (w_diff > 1E-08_dp .and. iter<100) .or. iter < 10 )
            iter = iter + 1

            w_current_next = w_temp 
            w_current_next(n_quad) = w_temp*exp(dz_vec(n_quad))

            w_temp_new = 0.0_dp

            do iii = 1,n_I
                ! Need to make sure this equation is correct & understand how the union works
                if (labor_alloc_vec(iii) > sqrt_eps) then  
                w_temp_new = w_temp_new +  weight_vec(iii)/sum(weight_vec)  *( &
                                vareps_w/ ( (1.0_dp-vareps_w) *  (1.0_dp-tau_w) )   &
                                * ( -  1.0_dp/IES_vec(iii) * c_vec(iii) * thtbar_vec(iii) * (l_aggr_temp**(1.0/tht)) / &
                                   (1.0_dp + (1.0_dp/IES_vec(iii) - 1.0_dp) * thtbar_vec(iii)*tht/(1.0_dp + tht) * &
                               l_aggr_temp**( (1.0_dp+tht)/tht ) ) )    &
                            + w_temp*(1.0/labor_alloc_vec(iii))*chiW/( (1.0_dp-vareps_w) *  (1.0-tau_w) ) * ( ( &
                                    w_temp/w_current*infl - 1.0 )   &
                                   *w_temp/w_current*infl  &
                                - sum( M_vec_mat(:,iii) * big_weight_vec * labor_frac_nxt(:,iii)* & 
                                  ( w_next_choice/w_current_next*infl_nxt - 1.0 ) * &
                                  (( w_next_choice/w_current_next)*( w_next_choice/w_temp)) * infl_nxt* l_aggr_nxt/l_aggr_temp ) &
                                 ) & 
                                 )                
                endif
            enddo

            w_diff = abs(w_temp_new-w_temp)
            w_temp = w_temp + 0.005_dp*(w_temp_new-w_temp)

            if (isnan(w_temp)) then
                write(*,*) 'error w isnan', weight_vec, M_vec_mat
                stop
            elseif (w_temp<=0) then
                write(*,*) 'error w<0'
                stop
            endif

            l_aggr_temp = k_aggr*(((1.0_dp-aalpha)/w_temp)**(1.0_dp/aalpha))
            w_choice_temp = w_temp

        enddo

        ! EXPORTED 
        l_aggr_new = l_aggr_temp

    else
        ! No wage ridigity
        l_aggr_new = sum(lmbd_vec*l_temp_vec)
    endif

    ! EXPORTED
    w_choice_new = w_choice

    ! =============================================================================
    ! (10) Update Inflation Rate  
    ! =============================================================================
    infl_new = ( exp(-tayl_ic) * exp(-state_grid(idx_m,sss) - monetary_shock)*nom_i_in  )**(1.0_dp/phi)

    constraint_binding = constraint_binding_new 

end subroutine calc_equilibrium_and_update



! ******************************************************************************** !
! Subroutine: Calculate Bond Prices
! ******************************************************************************** !
subroutine calc_bond_prices(nxt_mat, n_nxt, nxt_mat_2, nxt_bond_prices, sss, &
                            c_vec, l_aggr, q_current, q_l_current, infl, &
                            nom_i, share_vec, constraint_binding, who_prices_bond, bond_prices_out, results_vec)
    ! ============================================================ !
    ! Declarations
    ! ============================================================ !
    use mod_param, only: n_I, n_quad, idx_k, idx_m, idx_w, idx_dis, aalpha, smolyak_d,          &
                         quad_weight_vec, dz_vec, ddelta, state_grid, chiX, xi, ies_vec,        &
                         gma_vec, dz_vec_adj, bbeta_vec, thtbar_vec, tht, wealth_share_grid,    & 
                         n_bond, n_interp, lmbd_vec, s_bar_vec, vareps_w, tau_w, tayl_ic, phi,  &
                         s_trgt_vec, constr_agents, kbar_constr, low_guess_fixed, n_idio,       & 
                         high_guess_fixed,  labor_alloc_vec,  use_idio_risk, idio_weight_vec, idio_shock_grid
    real(dp), intent(in) :: nxt_mat_2(n_quad), l_aggr, nxt_bond_prices(n_quad), &
                            nxt_mat(n_quad, n_nxt), q_l_current(n_I), q_current, c_vec(n_I), infl 
    real(dp), intent(in) :: share_vec(n_I), nom_i
    real(dp), intent(out) :: who_prices_bond, bond_prices_out, results_vec(n_interp - n_bond*2)

    integer, intent(in)  :: sss, n_nxt, constraint_binding(n_I)

    integer :: iii, xxx, ggg, ifail, ccc, nnn

    real(dp) :: k_aggr, q_nxt(n_quad), rk_vec(n_quad), aggr_wealth, v_temp(n_quad, n_I),  &
                w_next_choice(n_quad), pi_nxt(n_quad), y_next(n_quad),  mc_temp(n_quad, n_I), &
                infl_nxt(n_quad), k_nxt(n_quad), rf_vec(n_quad), pi_current, y_current, &
                excess_b, k_temp_vec(n_I), b_temp_vec(n_I), c_temp(n_I), w_current, & 
                share_temp(n_I), thtbar, savings_vec(n_I), l_aggr_nxt(n_quad)

    real(dp) :: E_rk_new, E_rf_new, next_k

    real(dp), parameter :: brent_delta = 5E-16_dp, min_cons_sav = 1E-08_dp
 
    real(dp) :: w_choice, wealth, consumption, r_alpha_vec(n_quad), diff_inner, & 
                savings, objf, cons_update, IES, bbeta, next_period_share(n_quad), & 
                v_vec(n_quad), EV, M_vec_mat(n_quad, n_I), survivor_wealth(n_quad), & 
                mc_store(n_I), gma, labor, q_l_nxt(n_quad, n_I),  v_vec_twisted(n_quad), & 
                diff, util, util_c_deriv, labor_part, p_dis, big_weight_vec(n_quad),  &
                deriv_helper(n_quad), E_rf, E_rk, v_store(n_I)
 
    integer :: iter, iter_inner, kkk

    ! Individual Wealth Grid and Things
    integer, parameter :: n_wealth = 11

    real(dp) :: ratio_grid(n_wealth), cons_grid(n_wealth), &
                saving_grid(n_wealth), share_grid(n_wealth), wealth_grid(n_wealth), &
                capital_grid(n_wealth), deriv_temp(n_wealth), time_grid(n_wealth), labor_part_grid(n_wealth), &
                q_grid(n_wealth), nom_i_grid(n_wealth), q_inc(n_wealth), nom_i_inc(n_wealth), &
                E_rk_grid(n_wealth), E_rf_grid(n_wealth)

    real(dp) :: share_update, temp_vec(n_quad), share_high, share_low, &
                curr_share, temp, deriv, MPC_mat(n_I), MPS_mat(n_I), MPK_mat(n_I), &
                MPC_mat2(n_I), MPS_mat2(n_I), MPK_mat2(n_I), MPL_mat2(n_I), &
                time_choice, time_update, q_l, RHS, Time_Deviation_mat(n_I)

    real(dp) :: v_idio(n_quad*n_idio, n_I), mc_idio(n_quad*n_idio, n_I), q_l_nxt_idio(n_quad*n_idio, n_I), &
                dz_vec_idio(n_quad*n_idio), big_weight_vec_idio(n_quad*n_idio), &
                r_alpha_vec_idio(n_quad*n_idio), next_period_share_idio(n_quad*n_idio), &
                v_vec_twisted_idio(n_quad*n_idio), rf_vec_idio(n_quad*n_idio), &
                deriv_helper_idio(n_quad*n_idio), M_idio(n_quad*n_idio)

    real(dp) :: curr_q, dc_dErk_mat(n_I), db_dErk_mat(n_I), dk_dErk_mat(n_I), bond_grid(n_wealth), &
                rk_vec_use(n_quad), rf_vec_use(n_quad), curr_nom_i, dc_dErf_mat(n_I), & 
                db_dErf_mat(n_I), dk_dErf_mat(n_I)

    
    ! Prep grids for marginal response calculation 
    call Fill_linspace_dp(-1.0_dp, 1.0_dp, ratio_grid)
    call Fill_linspace_dp(0.0025_dp, -0.0025_dp, q_inc)
    call Fill_linspace_dp(-0.0025_dp, 0.0025_dp, nom_i_inc)

    ! ============================================================ !
    ! Current State and variables 
    ! ============================================================ !
    ! Same as in equilibrium computation 
    k_aggr = state_grid(idx_k, sss) ! Current Capital Stock 
    p_dis   = exp(state_grid(idx_dis, sss)) ! Current disaster prob
    w_current = state_grid(idx_w, sss) ! Current wage 
    big_weight_vec(1:n_quad-1) = quad_weight_vec*(1.0-p_dis)
    big_weight_vec(n_quad)     = p_dis
    pi_current = aalpha * (l_aggr/k_aggr)**(1.0-aalpha)
    aggr_wealth = k_aggr * (pi_current + (1.0-ddelta)*q_current)
    y_current = k_aggr**aalpha * l_aggr**(1.0-aalpha)
    w_choice = (1.0-aalpha)*y_current/l_aggr 

    ! ============================================================ !
    ! Quantities and Prices Next Period 
    ! ============================================================ !
    k_nxt = nxt_mat_2 
    do iii = 1,n_I
        v_temp(:,iii)  = nxt_mat(:,      iii) 
        mc_temp(:,iii) = nxt_mat(:,n_I + iii) 
    enddo
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
            do iii = 1,n_I
                v_idio(((kkk-1)*n_quad+1):(kkk*n_quad),iii)  = v_temp(:,iii)
                mc_idio(((kkk-1)*n_quad+1):(kkk*n_quad),iii) = mc_temp(:,iii)
                q_l_nxt_idio(((kkk-1)*n_quad+1):(kkk*n_quad),iii) = q_l_nxt(:,iii)*exp(idio_shock_grid(kkk,iii))
            enddo
            dz_vec_idio(((kkk-1)*n_quad+1):(kkk*n_quad)) = dz_vec
            rf_vec_idio(((kkk-1)*n_quad+1):(kkk*n_quad)) = rf_vec
            big_weight_vec_idio(((kkk-1)*n_quad+1):(kkk*n_quad)) = big_weight_vec * idio_weight_vec(kkk)
        enddo
    endif



    ! ============================================================ !
    ! Calculate Pricing Kernel 
    ! ============================================================ !
    ! Same as in part (7) of equilbrium solution 
    M_vec_mat = 0.0_dp

    do iii = 1,n_I
        ! (7.1) Read parameters and quantities (same as in part (4) above)
        IES   = IES_vec(iii)
        bbeta = bbeta_vec(iii)
        gma   = gma_vec(iii)
        thtbar  = thtbar_vec(iii)

        wealth = aggr_wealth*wealth_share_grid(iii,sss) 
        labor  = l_aggr*labor_alloc_vec(iii)

        consumption =  minval( [maxval([c_vec(iii), min_cons_sav]), & 
                          wealth + w_choice*labor-min_cons_sav] )
        r_alpha_vec = share_vec(iii)*rf_vec + (1.0_dp - share_vec(iii))*rk_vec
        savings     = wealth + w_choice*labor - consumption
        savings_vec(iii) = savings

        ! (7.2) Calculate (Realized) Pricing Kernel
        if (use_idio_risk == 1) then
            ! ------------------------------------------------------------ !
            ! With Idiosyncratic Risk
            ! ------------------------------------------------------------ !
            ! Expand portfolio return 
            do kkk = 1,n_idio
                r_alpha_vec_idio(((kkk-1)*n_quad+1):(kkk*n_quad)) = share_vec(iii)*rf_vec + (1.0_dp - share_vec(iii))*rk_vec*exp(idio_shock_grid(kkk,iii))
            enddo

            ! Other Operations
            next_period_share_idio = (savings*r_alpha_vec_idio + q_l_nxt_idio(:,iii) )/tot_wealth_ss_vec(iii) 
            v_vec_twisted_idio = (v_idio(:,iii)*next_period_share_idio) ** (1.0-gma)
            EV = sum(big_weight_vec_idio*v_vec_twisted_idio)

            call util_fun(consumption, labor, iii, util, util_c_deriv, labor_part)

            objf = ( (1.0-bbeta)*util +  bbeta*  (EV**( (1.0_dp - (1.0_dp/IES))/(1.0_dp-gma) )) )**(1.0_dp/(1.0_dp-1.0_dp/IES))

            mc_store(iii) = util_c_deriv
            v_store(iii) = objf*( (wealth + q_l_current(iii))/tot_wealth_ss_vec(iii))**(-1)

            ! SDF: first, the full form 
            M_idio = bbeta* ( next_period_share_idio ** (-gma) ) * (EV **( (gma - (1.0_dp/IES))/(1.0_dp-gma) ))  & 
                     * (v_idio(:,iii)**((1.0_dp/IES) -gma )) * mc_idio(:,iii) / util_c_deriv 

            ! Then for each of the n_quad aggregate nodes, average across the n_idio Idiosyncratic nodes
            M_vec_mat(:,iii)  = 0.0_dp 
            do kkk = 1,n_idio
                M_vec_mat(:,iii) = M_vec_mat(:,iii) + M_idio(((kkk-1)*n_quad+1):(kkk*n_quad)) * idio_weight_vec(kkk)
            enddo
        else
            ! ------------------------------------------------------------ !
            ! No Idiosyncratic risk
            ! ------------------------------------------------------------ !
            next_period_share = (savings*r_alpha_vec + q_l_nxt(:,iii) )/tot_wealth_ss_vec(iii)
            v_vec_twisted = (v_temp(:,iii)*next_period_share)**(1.0-gma)
            EV = sum(big_weight_vec*v_vec_twisted)

            call util_fun(consumption, labor, iii, util, util_c_deriv, labor_part)

            objf = ( (1.0-bbeta)*util +  bbeta*  (EV**( (1.0_dp - (1.0_dp/IES))/(1.0_dp-gma) )) )**(1.0_dp/(1.0_dp-1.0_dp/IES))
            
            mc_store(iii) = util_c_deriv
            v_store(iii) = objf*( (wealth + q_l_current(iii))/tot_wealth_ss_vec(iii))**(-1)

            M_vec_mat(:,iii) = bbeta* ( next_period_share ** (-gma) ) * (EV **( (gma - (1.0_dp/IES))/(1.0_dp-gma) ))  & 
                     * (v_temp(:,iii)**((1.0_dp/IES) -gma )) * mc_temp(:,iii) / util_c_deriv 

        endif


    enddo

    ! ============================================================ !
    ! Price Bonds
    ! ============================================================ !
    bond_prices_out = maxval([sum(big_weight_vec*M_vec_mat(:,1)*(nxt_bond_prices/infl_nxt)), &
                              sum(big_weight_vec*M_vec_mat(:,2)*(nxt_bond_prices/infl_nxt)), & 
                              sum(big_weight_vec*M_vec_mat(:,3)*(nxt_bond_prices/infl_nxt))])

    who_prices_bond = real(maxloc([sum(big_weight_vec*M_vec_mat(:,1)*(nxt_bond_prices/infl_nxt)), &
                              sum(big_weight_vec*M_vec_mat(:,2)*(nxt_bond_prices/infl_nxt)), & 
                              sum(big_weight_vec*M_vec_mat(:,3)*(nxt_bond_prices/infl_nxt))], 1, .true.), dp)

    
    ! ============================================================ !
    ! Calculate Expected Returns 
    ! ============================================================ !    
    E_rf = sum(big_weight_vec * rf_vec)    
    E_rk = sum(big_weight_vec * rk_vec)    

    ! ============================================================ !
    ! Calculate MPC and MPR: No Labor Endowment Trades
    ! ============================================================ !
    do iii=1,n_I
        IES   = IES_vec(iii)
        bbeta = bbeta_vec(iii)
        gma   = gma_vec(iii)
        thtbar  = thtbar_vec(iii)

        wealth = aggr_wealth*wealth_share_grid(iii,sss) 
        labor  = (l_aggr*labor_alloc_vec(iii))
        labor_part =  1.0_dp + (1.0_dp/IES - 1.0_dp) * thtbar*tht/(1.0_dp + tht) * labor**( (1.0_dp+tht)/tht )

        ! Initialize
        wealth_grid = ratio_grid + wealth
        cons_grid = c_vec(iii)
        share_grid = share_vec(iii)

        ! Loop Update
        do nnn = 1,n_wealth
            diff_inner = 1.0_dp 
            do while(diff_inner > 1E-05_dp)
                ! ============================================================ !
                ! Update Consumption Based on Portfolio Choice 
                ! ============================================================ !
                wealth = wealth_grid(nnn)
                curr_share = share_grid(nnn)
                r_alpha_vec = curr_share*rf_vec + (1.0_dp - curr_share)*rk_vec
                consumption = cons_grid(nnn)
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

                    if (constr_agents(iii) == 1 .and. constraint_binding(iii) == 1) then
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
                    if (constr_agents(iii) == 1 .and. constraint_binding(iii) == 1) then
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
                if (constr_agents(iii) == 1 .and. constraint_binding(iii) == 1) then
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
                cons_grid(nnn) = consumption + 0.5*(cons_update-consumption)
                share_grid(nnn) = share_update

            enddo

        enddo

        ! Compute MPC, MPR at the center point
        saving_grid = wealth_grid + w_choice*labor - cons_grid
        wealth = aggr_wealth*wealth_share_grid(iii,sss) 

        ! MPC
        call E01BEF( n_wealth, wealth_grid, cons_grid, deriv_temp, ifail)
        call E01BGF( n_wealth, wealth_grid, cons_grid, deriv_temp, 1 , wealth, temp, deriv, ifail)
        MPC_mat(iii) = deriv

        ! MPS (Savings)
        call E01BEF( n_wealth, wealth_grid, saving_grid, deriv_temp, ifail)
        call E01BGF( n_wealth, wealth_grid, saving_grid, deriv_temp, 1 , wealth, temp, deriv, ifail)
        MPS_mat(iii) = deriv

        ! MPK (Price * Capital)
        capital_grid = (1.0 - share_grid) * saving_grid
        call E01BEF( n_wealth, wealth_grid, capital_grid, deriv_temp, ifail)
        call E01BGF( n_wealth, wealth_grid, capital_grid, deriv_temp, 1 , wealth, temp, deriv, ifail)
        MPK_mat(iii) = deriv

        if (constr_agents(iii) == 1 .and. constraint_binding(iii) == 1) then
            MPK_mat(iii) = 0.0_dp
        endif
    enddo


    ! ============================================================ !
    ! Derivatives to Price of Capital 
    ! ============================================================ !
    do iii=1,n_I
        IES   = IES_vec(iii)
        bbeta = bbeta_vec(iii)
        gma   = gma_vec(iii)
        thtbar  = thtbar_vec(iii)

        wealth = aggr_wealth*wealth_share_grid(iii,sss) 
        labor  = (l_aggr*labor_alloc_vec(iii))
        labor_part =  1.0_dp + (1.0_dp/IES - 1.0_dp) * thtbar*tht/(1.0_dp + tht) * labor**( (1.0_dp+tht)/tht )

        ! Initialize
        q_grid = q_current + q_inc
        cons_grid = c_vec(iii)
        share_grid = share_vec(iii)

        ! Loop Update
        do nnn = 1,n_wealth
            diff_inner = 1.0_dp 
            do while(diff_inner > 1E-05_dp)
                ! ============================================================ !
                ! Update Consumption Based on Portfolio Choice 
                ! ============================================================ !
                curr_q = q_grid(nnn)
                curr_share = share_grid(nnn)

                rk_vec_use = ( (1.0_dp-ddelta)*q_nxt + pi_nxt ) / curr_q 
                rk_vec_use(n_quad) = exp(dz_vec(n_quad)) * rk_vec_use(n_quad)

                r_alpha_vec = curr_share*rf_vec + (1.0_dp - curr_share)*rk_vec_use
                consumption = cons_grid(nnn)
                savings = wealth + w_choice*labor - consumption

                if (use_idio_risk == 1) then
                    ! ------------------------------------------------------------ !
                    ! With Idiosyncratic Risk
                    ! ------------------------------------------------------------ !
                    ! Expand portfolio return 
                    do kkk = 1,n_idio
                        r_alpha_vec_idio(((kkk-1)*n_quad+1):(kkk*n_quad)) = curr_share*rf_vec + (1.0_dp - curr_share)*rk_vec_use*exp(idio_shock_grid(kkk,iii))
                    enddo

                    ! Other Operations
                    next_period_share_idio = (savings*r_alpha_vec_idio + q_l_nxt_idio(:,iii) )/tot_wealth_ss_vec(iii)
                    v_vec_twisted_idio = (v_idio(:,iii)*next_period_share_idio) ** (1.0-gma)
                    EV = sum(big_weight_vec_idio*v_vec_twisted_idio)

                    deriv_helper_idio = bbeta * (v_idio(:,iii)**((1.0_dp/IES) - gma )) * (next_period_share_idio**(-gma)) & 
                                      * (EV **( (gma - (1.0_dp/IES))/(1.0_dp-gma) )) * mc_idio(:,iii)

                    if (constr_agents(iii) == 1 .and. constraint_binding(iii) == 1) then
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
                    if (constr_agents(iii) == 1 .and. constraint_binding(iii) == 1) then
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
                if (constr_agents(iii) == 1 .and. constraint_binding(iii) == 1) then
                    share_update = 1.0_dp - kbar_constr(iii)*curr_q/savings
                else

                    ! find natural leverage limit to prevent negative payoffs
                    temp_vec = -rk_vec_use/(rf_vec-rk_vec_use + sqrt_eps)
                    where (temp_vec >= 0.0_dp)
                        temp_vec = -10000.0_dp
                    end where
                    share_low = maxval([low_guess_fixed, maxval(temp_vec)+ sqrt_eps])

                    ! find natural leverage limit to prevent negative payoffs
                    temp_vec = -rk_vec_use/(rf_vec-rk_vec_use + sqrt_eps)
                    where (temp_vec <= 0.0_dp)
                        temp_vec = 10000.0_dp
                    end where
                    share_high = minval([high_guess_fixed, minval(temp_vec) - sqrt_eps ])

                    ! Calculate Portfolio Share, Given Updated Consumption and Savings
                    share_update = curr_share
                    call calc_portfolio_share(1.0_dp, share_low, share_high, iii, big_weight_vec, savings, labor, cons_update,  &
                                      0.0_dp, q_l_nxt(:,iii), v_temp(:,iii), mc_temp(:,iii), rf_vec, rk_vec_use, q_current, &
                                      share_update)

                endif

                ! ============================================================ !
                ! Update Values
                ! ============================================================ !
                diff_inner = abs(consumption-cons_update)   
                cons_grid(nnn) = consumption + 0.5*(cons_update-consumption)
                share_grid(nnn) = share_update

            enddo

        enddo

        ! Compute MPC, MPR at the center point
        saving_grid = wealth + w_choice*labor - cons_grid
        E_rk_grid = E_rk * q_current / q_grid

        ! MPC
        call E01BEF( n_wealth, E_rk_grid, cons_grid, deriv_temp, ifail)
        call E01BGF( n_wealth, E_rk_grid, cons_grid, deriv_temp, 1 , E_rk, temp, deriv, ifail)
        dc_dErk_mat(iii) = deriv

        ! MPB Actual (Bond)
        bond_grid = share_grid * saving_grid
        call E01BEF( n_wealth, E_rk_grid, bond_grid, deriv_temp, ifail)
        call E01BGF( n_wealth, E_rk_grid, bond_grid, deriv_temp, 1 , E_rk, temp, deriv, ifail)
        db_dErk_mat(iii) = deriv

        ! MPK Actual (Capital)
        capital_grid = (1-share_grid) * saving_grid / q_grid
        call E01BEF( n_wealth, E_rk_grid, capital_grid, deriv_temp, ifail)
        call E01BGF( n_wealth, E_rk_grid, capital_grid, deriv_temp, 1 , E_rk, temp, deriv, ifail)
        dk_dErk_mat(iii) = deriv

        if (constr_agents(iii) == 1 .and. constraint_binding(iii) == 1) then
            dk_dErk_mat(iii) = 0.0_dp
        endif

    enddo

    ! ============================================================ !
    ! Derivatives to Nominal Rate
    ! ============================================================ !
    do iii=1,n_I
        IES   = IES_vec(iii)
        bbeta = bbeta_vec(iii)
        gma   = gma_vec(iii)
        thtbar  = thtbar_vec(iii)

        wealth = aggr_wealth*wealth_share_grid(iii,sss) 
        labor  = (l_aggr*labor_alloc_vec(iii))
        labor_part =  1.0_dp + (1.0_dp/IES - 1.0_dp) * thtbar*tht/(1.0_dp + tht) * labor**( (1.0_dp+tht)/tht )

        ! Initialize
        nom_i_grid = nom_i + nom_i_inc
        cons_grid = c_vec(iii)
        share_grid = share_vec(iii)

        ! Loop Update
        do nnn = 1,n_wealth
            diff_inner = 1.0_dp 
            do while(diff_inner > 1E-05_dp)
                ! ============================================================ !
                ! Update Consumption Based on Portfolio Choice 
                ! ============================================================ !
                curr_nom_i = nom_i_grid(nnn)
                curr_share = share_grid(nnn)

                rf_vec_use = curr_nom_i / infl_nxt

                r_alpha_vec = curr_share*rf_vec_use + (1.0_dp - curr_share)*rk_vec
                consumption = cons_grid(nnn)
                savings = wealth + w_choice*labor - consumption

                if (use_idio_risk == 1) then
                    ! ------------------------------------------------------------ !
                    ! With Idiosyncratic Risk
                    ! ------------------------------------------------------------ !
                    ! Expand portfolio return 
                    do kkk = 1,n_idio
                        r_alpha_vec_idio(((kkk-1)*n_quad+1):(kkk*n_quad)) = curr_share*rf_vec_use + (1.0_dp - curr_share)*rk_vec*exp(idio_shock_grid(kkk,iii))
                        rf_vec_idio(((kkk-1)*n_quad+1):(kkk*n_quad)) = rf_vec_use 
                    enddo

                    ! Other Operations
                    next_period_share_idio = (savings*r_alpha_vec_idio + q_l_nxt_idio(:,iii) )/tot_wealth_ss_vec(iii)
                    v_vec_twisted_idio = (v_idio(:,iii)*next_period_share_idio) ** (1.0-gma)
                    EV = sum(big_weight_vec_idio*v_vec_twisted_idio)

                    deriv_helper_idio = bbeta * (v_idio(:,iii)**((1.0_dp/IES) - gma )) * (next_period_share_idio**(-gma)) & 
                                      * (EV **( (gma - (1.0_dp/IES))/(1.0_dp-gma) )) * mc_idio(:,iii)

                    if (constr_agents(iii) == 1 .and. constraint_binding(iii) == 1) then
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
                    if (constr_agents(iii) == 1 .and. constraint_binding(iii) == 1) then
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
                if (constr_agents(iii) == 1 .and. constraint_binding(iii) == 1) then
                    share_update = 1.0_dp - kbar_constr(iii)*q_current/savings
                else

                    ! find natural leverage limit to prevent negative payoffs
                    temp_vec = -rk_vec/(rf_vec_use-rk_vec + sqrt_eps)
                    where (temp_vec >= 0.0_dp)
                        temp_vec = -10000.0_dp
                    end where
                    share_low = maxval([low_guess_fixed, maxval(temp_vec)+ sqrt_eps])

                    ! find natural leverage limit to prevent negative payoffs
                    temp_vec = -rk_vec/(rf_vec_use-rk_vec + sqrt_eps)
                    where (temp_vec <= 0.0_dp)
                        temp_vec = 10000.0_dp
                    end where
                    share_high = minval([high_guess_fixed, minval(temp_vec) - sqrt_eps ])

                    ! Calculate Portfolio Share, Given Updated Consumption and Savings
                    share_update = curr_share
                    call calc_portfolio_share(1.0_dp, share_low, share_high, iii, big_weight_vec, savings, labor, cons_update,  &
                                      0.0_dp, q_l_nxt(:,iii), v_temp(:,iii), mc_temp(:,iii), rf_vec_use, rk_vec, q_current, &
                                      share_update)

                endif

                ! ============================================================ !
                ! Update Values
                ! ============================================================ !
                diff_inner = abs(consumption-cons_update)   
                cons_grid(nnn) = consumption + 0.5*(cons_update-consumption)
                share_grid(nnn) = share_update

            enddo

        enddo

        ! Compute MPC, MPR at the center point
        saving_grid = wealth + w_choice*labor - cons_grid
        E_rf_grid = E_rf * nom_i_grid / nom_i

        ! MPC
        call E01BEF( n_wealth, E_rf_grid, cons_grid, deriv_temp, ifail)
        call E01BGF( n_wealth, E_rf_grid, cons_grid, deriv_temp, 1 , E_rf, temp, deriv, ifail)
        dc_dErf_mat(iii) = deriv

        ! MPB Actual (Bond)
        bond_grid = share_grid * saving_grid
        call E01BEF( n_wealth, E_rf_grid, bond_grid, deriv_temp, ifail)
        call E01BGF( n_wealth, E_rf_grid, bond_grid, deriv_temp, 1 , E_rf, temp, deriv, ifail)
        db_dErf_mat(iii) = deriv

        ! MPK Actual (Capital)
        capital_grid = (1-share_grid) * saving_grid / q_current
        call E01BEF( n_wealth, E_rf_grid, capital_grid, deriv_temp, ifail)
        call E01BGF( n_wealth, E_rf_grid, capital_grid, deriv_temp, 1 , E_rf, temp, deriv, ifail)
        dk_dErf_mat(iii) = deriv

        if (constr_agents(iii) == 1 .and. constraint_binding(iii) == 1) then
            dk_dErf_mat(iii) = 0.0_dp
        endif
    enddo
    
    ! also adjust to log returns above
    E_rf = sum(big_weight_vec * log(rf_vec) )   
    E_rk = sum(big_weight_vec * log(rk_vec) )   

    ! ============================================================ !
    ! Store all relevant variables
    ! ============================================================ !
    results_vec = 0.0_dp
    
    results_vec(1) = y_current 
    results_vec(2) = l_aggr
    results_vec(3) = c_vec(1)
    results_vec(4) = c_vec(2)
    results_vec(5) = c_vec(3)
    results_vec(6) = k_nxt(1) - (1.0-ddelta)*k_aggr
    results_vec(7) = pi_current
    results_vec(8) = E_rf 
    results_vec(9) = E_rk 
    results_vec(10) = q_current 
    results_vec(11) = share_vec(1)
    results_vec(12) = share_vec(2)
    results_vec(13) = share_vec(3)
    results_vec(14) = sum(lmbd_vec*savings_vec)
    results_vec(15) = aggr_wealth
    results_vec(16) = w_choice
    results_vec(17) = E_rk - E_rf 
    results_vec(18) = infl 
    results_vec(19) = nom_i 
    results_vec(20) = sum(lmbd_vec*savings_vec*(1.0-share_vec))
    results_vec(21) = k_aggr 
    results_vec(22) = v_store(1)
    results_vec(23) = v_store(2)
    results_vec(24) = v_store(3)
    results_vec(25) = q_l_current(1)
    results_vec(26) = q_l_current(2)
    results_vec(27) = q_l_current(3)
    results_vec(28) = MPC_mat(1)
    results_vec(29) = MPC_mat(2)
    results_vec(30) = MPC_mat(3)
    results_vec(31) = MPS_mat(1)
    results_vec(32) = MPS_mat(2)
    results_vec(33) = MPS_mat(3)
    results_vec(34) = MPK_mat(1)
    results_vec(35) = MPK_mat(2)
    results_vec(36) = MPK_mat(3)
    results_vec(37) = dc_dErk_mat(1)
    results_vec(38) = dc_dErk_mat(2)
    results_vec(39) = dc_dErk_mat(3)
    results_vec(40) = db_dErk_mat(1)
    results_vec(41) = db_dErk_mat(2)
    results_vec(42) = db_dErk_mat(3)
    results_vec(43) = dk_dErk_mat(1)
    results_vec(44) = dk_dErk_mat(2)
    results_vec(45) = dk_dErk_mat(3)
    results_vec(46) = dc_dErf_mat(1)
    results_vec(47) = dc_dErf_mat(2)
    results_vec(48) = dc_dErf_mat(3)
    results_vec(49) = db_dErf_mat(1)
    results_vec(50) = db_dErf_mat(2)
    results_vec(51) = db_dErf_mat(3)
    results_vec(52) = dk_dErf_mat(1)
    results_vec(53) = dk_dErf_mat(2)
    results_vec(54) = dk_dErf_mat(3)
    results_vec(55) = lmbd_vec(1)*savings_vec(1)*(1.0-share_vec(1))/q_current
    results_vec(56) = lmbd_vec(2)*savings_vec(2)*(1.0-share_vec(2))/q_current
    results_vec(57) = lmbd_vec(3)*savings_vec(3)*(1.0-share_vec(3))/q_current
    results_vec(58) = lmbd_vec(1)*savings_vec(1)
    results_vec(59) = lmbd_vec(2)*savings_vec(2)
    results_vec(60) = lmbd_vec(3)*savings_vec(3)
    results_vec(61) = lmbd_vec(1)*aggr_wealth*wealth_share_grid(1,sss)
    results_vec(62) = lmbd_vec(2)*aggr_wealth*wealth_share_grid(2,sss)
    results_vec(63) = lmbd_vec(3)*aggr_wealth*wealth_share_grid(3,sss)
    results_vec(64) = sum(big_weight_vec*M_vec_mat(:,3)*rk_vec)
    results_vec(65) = constraint_binding(1) 
    results_vec(66) = constraint_binding(2)
    results_vec(67) = constraint_binding(3)
    results_vec(68) = mc_store(1) 
    results_vec(69) = mc_store(2)
    results_vec(70) = mc_store(3)

end subroutine calc_bond_prices



! ******************************************************************************** !
! Subroutine: Calculate Excess Bond Demand 
! ******************************************************************************** !
subroutine calc_excess_bond_nom(aggr_wealth, v_temp, mc_temp, rf, rk,  q_current, w_choice, &
                                l, next_k_value, q_l_nxt, sss, big_weight_vec, c_temp, share_guess_vec, & 
                                k_temp_vec, b_temp_vec, excess_b, constraint_binding, constraint_binding_new)
    
    ! Declarations
    use mod_param, only: n_quad, n_I, low_guess_fixed, high_guess_fixed, &
                         n_states, state_grid, wealth_share_grid, lmbd_vec, &
                         kbar_constr, constr_agents, lmbd_vec, labor_alloc_vec

    integer, intent(in)  :: sss, constraint_binding(n_I)
    integer, intent(out) :: constraint_binding_new(n_I)
    real(dp), intent(in) :: aggr_wealth, v_temp(n_quad, n_I), mc_temp(n_quad, n_I), q_l_nxt(n_quad, n_I), & 
                            rf(n_quad), rk(n_quad), big_weight_vec(n_quad),  q_current, w_choice, l, next_k_value, c_temp(n_I)
    real(dp), intent(inout) :: share_guess_vec(n_I)
    real(dp), intent(out)   :: k_temp_vec(n_I), b_temp_vec(n_I), excess_b

    ! Local Variables
    real(dp) :: share_guess, share_low, share_high, wealth, temp(n_quad), savings, objf, labor
    integer  :: iii
    
    constraint_binding_new = 0
    ! Calculate bond demand for each agent 
    do iii = 1,n_I
        ! Basic information for the agent 
        wealth = aggr_wealth*wealth_share_grid(iii,sss)
        savings =  wealth + w_choice*l*labor_alloc_vec(iii) - c_temp(iii)
        if (savings <0.0_dp) then
            write(*,*) 'found negative savings here'
            write(*,*) savings, sss
            stop
        endif


        ! find natural leverage limit to prevent negative payoffs
        temp = -rk/(rf-rk + sqrt_eps)
        where (temp >= 0.0_dp)
            temp = -10000.0_dp
        end where
        share_low = maxval([low_guess_fixed, maxval(temp)+ sqrt_eps])

        ! find natural leverage limit to prevent negative payoffs
        temp = -rk/(rf-rk + sqrt_eps)
        where (temp <= 0.0_dp)
            temp = 10000.0_dp
        end where
        share_high = minval([high_guess_fixed, minval(temp) - sqrt_eps ])

        ! start from share guess but stay in limits
        share_guess = maxval([minval([share_guess_vec(iii), share_high]), share_low]) 

        ! calculate the portfolio share given the initial guess 
        call calc_portfolio_share(1.0_dp, share_low, share_high, iii, big_weight_vec, savings, l, c_temp(iii),  &
                                  next_k_value, q_l_nxt(:,iii), v_temp(:,iii), mc_temp(:,iii), rf, rk, q_current, share_guess)
        
        share_guess_vec(iii) = share_guess

        if (constr_agents(iii) == 1) then ! occasionally binding constraint
            share_guess = 1.0_dp - kbar_constr(iii)*q_current/savings
            if (share_guess < share_guess_vec(iii)) then  
                constraint_binding_new(iii) = 1
            endif
            if (constraint_binding(iii) == 1) then  
                share_guess_vec(iii) = share_guess 
            endif
        endif

        ! Store calculated portfolio shares and bond demand 
        k_temp_vec(iii) = savings*(1.0_dp-share_guess_vec(iii))/q_current
        b_temp_vec(iii) = share_guess_vec(iii)*savings


    enddo

    ! Excess Demand 
    excess_b =  sum(lmbd_vec*b_temp_vec) 
    
end subroutine calc_excess_bond_nom

! ******************************************************************************** !
! Subroutine: Calculate optimal portfolio given a total savings decision
! ******************************************************************************** !
subroutine calc_portfolio_share(time_choice, share_low, share_high, iii, big_weight_vec, &
                                savings, labor, consumption, next_k_value, q_l_nxt, &
                                v, mc, rf, rk, q_current, share_guess)

    ! ------------------------------------------------------------ !
    ! Declarations
    ! ------------------------------------------------------------ !
    use mod_param, only: gma_vec, ies_vec, bbeta_vec, n_quad, dz_vec, n_idio, &
                         use_idio_risk, idio_weight_vec, idio_shock_grid 
                     

    real(dp), intent(in)    :: time_choice, v(n_quad), share_low, share_high, q_l_nxt(n_quad),&
                               rf(n_quad), rk(n_quad), savings, next_k_value, labor, q_current, &
                               big_weight_vec(n_quad), consumption, mc(n_quad)
    integer, intent(in)     :: iii
    real(dp), intent(inout) :: share_guess

    real(dp) :: IES, share_FOC_low, share_FOC_hi, share_guess_A, share_guess_B, share_guess_C, share_guess_S,   & 
                share_guess_D, share_FOC, share_FOC_A, share_FOC_B, share_FOC_C, share_FOC_S, temp_vec(n_quad), & 
                bbeta, temp, brent_delta, share_temp, gma, E_MR

    integer :: mflag, kkk

    real(dp) :: v_use(n_quad*n_idio), q_l_nxt_use(n_quad*n_idio), rf_use(n_quad*n_idio), rk_use(n_quad*n_idio), &
                big_weight_vec_use(n_quad*n_idio), mc_use(n_quad*n_idio), r_alpha(n_quad*n_idio), & 
                next_period_share(n_quad*n_idio), v_vec_twisted(n_quad*n_idio), dz_vec_use(n_quad*n_idio)

    ! ------------------------------------------------------------ !
    ! solve portfolio choice through E[M * (rk -rf)] = 0
    ! ------------------------------------------------------------ !

    ! preference parameters
    gma             = gma_vec(iii)
    ies             = ies_vec(iii)
    bbeta           = bbeta_vec(iii)

    ! ------------------------------------------------------------ !
    ! Initialize: expand for idiosyncratic risk
    ! ------------------------------------------------------------ !
    if (use_idio_risk == 1) then
        ! expand for idiosyncratic risk
        do kkk = 1,n_idio
            
            v_use(((kkk-1)*n_quad+1):(kkk*n_quad))      = v
            mc_use(((kkk-1)*n_quad+1):(kkk*n_quad))     = mc
            rf_use(((kkk-1)*n_quad+1):(kkk*n_quad))     = rf
            dz_vec_use(((kkk-1)*n_quad+1):(kkk*n_quad)) = dz_vec

            ! risky assets
            q_l_nxt_use(((kkk-1)*n_quad+1):(kkk*n_quad)) = q_l_nxt*exp(idio_shock_grid(kkk,iii))
            rk_use(((kkk-1)*n_quad+1):(kkk*n_quad)) = rk*exp(idio_shock_grid(kkk,iii)) 

            ! weights 
            big_weight_vec_use(((kkk-1)*n_quad+1):(kkk*n_quad)) = big_weight_vec * idio_weight_vec(kkk)
        enddo
    else
        v_use = v 
        mc_use = mc 
        q_l_nxt_use = q_l_nxt 
        rf_use = rf 
        rk_use = rk 
        big_weight_vec_use = big_weight_vec
        dz_vec_use = dz_vec
    endif


    ! ------------------------------------------------------------ !
    ! Start Solving
    ! ------------------------------------------------------------ !
    share_temp = share_low
    r_alpha = share_temp*rf_use + (1.0_dp -share_temp)*rk_use
    next_period_share = (savings*r_alpha  + q_l_nxt_use*time_choice) / tot_wealth_ss_vec(iii)
    v_vec_twisted = v_use**(1.0/ies-gma) * mc_use
    E_MR = sum( v_vec_twisted * (next_period_share**(-gma)) *(rf_use - rk_use) * big_weight_vec_use) &
            /abs(sum(big_weight_vec_use*v_vec_twisted))

    share_FOC_low = E_MR   

    ! share_FOC<0 means optimum lies to the left
    ! share_FOC>0 means optimum lies to the right 
    if (share_FOC_low <= 0.0_dp) then
        ! No interior solution --> just use low value 
        share_guess = share_low
    else
        ! If not satisfied, use brent method to get port. share 

        ! ================================================== !
        ! (1) Evaluate FOC at high end 
        ! ================================================== !
        share_temp = share_high
        r_alpha = share_temp*rf_use + (1.0_dp -share_temp)*rk_use
        next_period_share = (savings*r_alpha  + q_l_nxt_use*time_choice) / tot_wealth_ss_vec(iii)
        v_vec_twisted = v_use**(1.0/ies-gma) * mc_use
        E_MR = sum( v_vec_twisted * (next_period_share**(-gma)) *(rf_use - rk_use) * big_weight_vec_use) &
                /abs(sum(big_weight_vec_use*v_vec_twisted))

        share_FOC_hi = E_MR    

        ! ================================================== !
        ! (2) Use Brent Method 
        ! ================================================== !
        if (share_FOC_hi >= 0.0_dp) then 
            ! No interior solution --> just use high value 
            share_guess = share_high
        else
            ! Evaluate FOC at current guess 
            share_temp = share_guess
            r_alpha = share_temp*rf_use + (1.0_dp -share_temp)*rk_use
            next_period_share = (savings*r_alpha  + q_l_nxt_use*time_choice) / tot_wealth_ss_vec(iii)
            v_vec_twisted = v_use**(1.0/ies-gma) * mc_use
            E_MR = sum( v_vec_twisted * (next_period_share**(-gma)) *(rf_use - rk_use) * big_weight_vec_use) &
                    /abs(sum(big_weight_vec_use*v_vec_twisted))

            share_FOC = E_MR   

            ! Update using Brent Method 
            if (share_FOC /= 0.0_dp) then
                share_guess_B = share_guess ! make this highest value
                share_FOC_B   = share_FOC

                if (share_FOC <= 0.0_dp) then ! share_FOC<0 means optimum lies to the left
                    share_guess_A = share_low
                    share_temp = share_guess_A
                    r_alpha = share_temp*rf_use + (1.0_dp -share_temp)*rk_use
                    next_period_share = (savings*r_alpha  + q_l_nxt_use*time_choice) / tot_wealth_ss_vec(iii)
                    v_vec_twisted = v_use**(1.0/ies-gma) * mc_use
                    E_MR = sum( v_vec_twisted * (next_period_share**(-gma)) *(rf_use - rk_use) * big_weight_vec_use) &
                            /abs(sum(big_weight_vec_use*v_vec_twisted))

                    share_FOC_A = E_MR   
                else ! share_FOC > 0 so optimum lies to the right
                    share_guess_A = share_high
                    share_temp = share_guess_A
                    r_alpha = share_temp*rf_use + (1.0_dp -share_temp)*rk_use
                    next_period_share = (savings*r_alpha  + q_l_nxt_use*time_choice) / tot_wealth_ss_vec(iii)
                    v_vec_twisted = v_use**(1.0/ies-gma) * mc_use
                    E_MR = sum( v_vec_twisted * (next_period_share**(-gma)) *(rf_use - rk_use) * big_weight_vec_use) &
                            /abs(sum(big_weight_vec_use*v_vec_twisted))

                    share_FOC_A = E_MR   
                endif

                ! check that root is bracketed
                if ((share_FOC_A*share_FOC_B) > 0) then
                    write(*,*) 'ERROR: Initial bracket does not contain root.'
                    stop
                endif

                ! swap a and b as needed 
                if ( abs(share_FOC_A) < abs(share_FOC_B) )then
                    temp = share_guess_A
                    share_guess_A = share_guess_B
                    share_guess_B = temp

                    temp = share_FOC_A
                    share_FOC_A = share_FOC_B
                    share_FOC_B = temp
                endif

                ! Brent Solution 
                share_guess_C    = share_guess_A
                share_FOC_C      = share_FOC_A
                share_FOC_S      = share_FOC_B
                mflag = 1 ! set mflag
                brent_delta = 5E-16_dp
                do while ( share_FOC_S /= 0 .and. abs(share_guess_A - share_guess_B) > 1E-14_dp)
                    ! ************************************************************ !
                    ! Obtain updated position
                    ! ************************************************************ !
                    if ( (share_FOC_A /= share_FOC_C) .and. (share_FOC_C /= share_FOC_B) ) then ! inverse quadratic interpolation
                        share_guess_S = share_guess_A * share_FOC_B * share_FOC_C / ( (share_FOC_A - share_FOC_B) * (share_FOC_A - share_FOC_C) ) + &
                                 share_guess_B * share_FOC_A * share_FOC_C / ( (share_FOC_B - share_FOC_A) * (share_FOC_B - share_FOC_C) ) + &
                                 share_guess_C * share_FOC_A * share_FOC_B / ( (share_FOC_C - share_FOC_A) * (share_FOC_C - share_FOC_B) )
                    else ! secant method
                        share_guess_S = share_guess_B - share_FOC_B * (share_guess_B - share_guess_A) /   (share_FOC_B - share_FOC_A)
                    endif

                    if ( ( (  ( share_guess_S > ( 3*share_guess_A + share_guess_B )/4 .and. share_guess_S < share_guess_B) .or. &
                           ( share_guess_S < ( 3*share_guess_A + share_guess_B )/4 .and. share_guess_S > share_guess_B)  ) == .FALSE. ) .or. &
                      (mflag == 1 .and. abs(share_guess_S - share_guess_B) >= abs(share_guess_B - share_guess_C)/2  )             .or. &
                      (mflag == 0 .and. abs(share_guess_S - share_guess_B) >= abs(share_guess_B - share_guess_D)/2  )             .or. &
                      (mflag == 1 .and. abs(share_guess_B - share_guess_C) <  abs(brent_delta)  )                     .or. &
                      (mflag == 0 .and. abs(share_guess_B - share_guess_D) <  abs(brent_delta)  )                          &
                    ) then

                        share_guess_S = (share_guess_A + share_guess_B )/ 2
                        mflag = 1
                    else
                        mflag = 0
                    endif

                    ! ************************************************************ !
                    ! Evaluate FOC at share_guess_S
                    ! ************************************************************ !
                    share_temp = share_guess_S
                    r_alpha = share_temp*rf_use + (1.0_dp -share_temp)*rk_use
                    next_period_share = (savings*r_alpha  + q_l_nxt_use*time_choice) / tot_wealth_ss_vec(iii)
                    v_vec_twisted = v_use**(1.0/ies-gma) * mc_use
                    E_MR = sum( v_vec_twisted * (next_period_share**(-gma)) *(rf_use - rk_use) * big_weight_vec_use) &
                            /abs(sum(big_weight_vec_use*v_vec_twisted))

                    share_FOC_S = E_MR

                    ! ************************************************************ !
                    ! Update Boundaries
                    ! ************************************************************ !
                    share_guess_D = share_guess_C
                    share_guess_C = share_guess_B
                    share_FOC_C = share_FOC_B

                    if ((share_FOC_A*share_FOC_S) < 0) then
                        share_guess_B = share_guess_S
                        share_FOC_B = share_FOC_S
                    else
                        share_guess_A = share_guess_S
                        share_FOC_A = share_FOC_S
                    endif

                    ! swap a and b as needed 
                    if ( abs(share_FOC_A) < abs(share_FOC_B) ) then
                        temp = share_guess_A
                        share_guess_A = share_guess_B
                        share_guess_B = temp

                        temp = share_FOC_A
                        share_FOC_A = share_FOC_B
                        share_FOC_B = temp
                    endif                    
                enddo

                ! EXPORT 
                share_guess = share_guess_S
            endif
        endif
    endif


end subroutine calc_portfolio_share


! ******************************************************************************** !
! Subroutine: instantaenous utility and its derivative
! ******************************************************************************** !
subroutine util_fun(consumption, labor, iii, util, util_c_deriv, labor_part) 
    use mod_param, only: thtbar_vec, ies_vec, tht

    integer, intent(in)   :: iii
    real(dp), intent(in)  :: consumption, labor
    real(dp), intent(out) :: util, util_c_deriv, labor_part

    real(dp) :: thtbar, ies 

    thtbar = thtbar_vec(iii)
    ies  = ies_vec(iii)

    if (IES > 1.0_dp-sqrt_eps .and. IES < 1.0_dp + sqrt_eps) then 
        labor_part = 1.0_dp
        util = log(consumption) -  thtbar*tht/(1.0_dp + tht) * labor**( (1.0_dp+tht)/tht ) 
        util_c_deriv = 1.0_dp/consumption    
    else 
        labor_part =  (1.0_dp + (1.0_dp/IES - 1.0_dp) * thtbar*tht/(1.0_dp + tht) * labor**( (1.0_dp+tht)/tht ) )**(1.0_dp/IES)
        util = (consumption**(1.0_dp-1.0_dp/IES)) * labor_part
        util_c_deriv = util/consumption    
    endif

end subroutine util_fun

end module mod_calc
