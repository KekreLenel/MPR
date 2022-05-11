! -------------------------------------------------------------------------
! mod_results.f90: results module, public subroutines:
! - create_results: given model solution, creates and stores simulations and IRFs 
! -------------------------------------------------------------------------
! authors:         Rohan Kekre and Moritz Lenel
! for updates see: https://github.com/KekreLenel/MPR
! -------------------------------------------------------------------------
module mod_results

use omp_lib
use base_lib, only: dp, sqrt_eps, Fill_linspace_dp

implicit none
private

public :: create_results

contains

subroutine create_results()
    use mod_smolyak, only : Smolyak_Polynomial
    use mod_param, only : n_states, idx_k, idx_sa, idx_sc, idx_m, idx_w, idx_dis, &
                          state_grid, smol_polynom, max_smol_level, smol_elem_ani, &
                          k_grid_mean, k_grid_dev, s_a_grid_mean, s_a_grid_dev, &
                          s_c_grid_mean, s_c_grid_dev, m_grid_mean, m_grid_dev, &
                          w_grid_mean, w_grid_dev, dis_grid_mean, dis_grid_dev, lmbd_vec, &
                          n_shocks, kbar_constr, results_folder, n_sim_periods, & 
                          n_prm, n_burn, n_irf_periods, n_sample, n_quad, smolyak_d, &
                          vector_mus_dimensions, irf_indices, irf_shock_sizes, &
                          sidx_z, sidx_m, sidx_dis, n_interp, n_bond, n_active_dims, n_I, &
                          smol_grid, no_shock_idx, shock_grid, quad_weight_vec, jx_1y
    
    character :: command_input*2, filename*100
    integer   :: narg, command_len, sss, fff, sample_idx, ppp, ddd, counter
    integer   ::  nrhs, nrhs2, nrhs3, lda, ldaf, ldb, ldx
    integer   ::  ipiv(n_states), iwork(n_states), info

    INTEGER, parameter :: GENID = 3, SUBID = 1, LSEED = 1
    INTEGER            :: SEED(LSEED), STATE(633), LSTATE, IFAIL
    integer, parameter :: n_spread = 20
    real(dp) :: state_spread_vec(n_spread) 
    real(dp), allocatable :: state_series(:,:), other_vars_series(:,:), shock_series(:,:), simshock_series(:), &
                             state_spread_series(:,:,:), shock_grid_temp(:,:), shock_tmp(:), other_vars_series_tmp(:,:)

    real(dp), allocatable :: b(:,:), interp_coeffs(:,:), state_coeffs(:,:), smol_coeffs(:,:,:), & 
                             ferr(:), berr(:), b2(:,:), ferr2(:), berr2(:), b3(:,:), ferr3(:), berr3(:)

    real(dp) ::  a(n_states,n_states), af(n_states,n_states), r(n_states), & 
                 c(n_states),  work(4*n_states), polyn_points(1,n_states)

    real(dp) :: next_state_mat(n_quad, n_active_dims, n_states), &
                next_state_monetary_mat(n_quad, n_active_dims, n_states), & 
                next_state_montransition_mat(1, n_active_dims, n_states)

    real(dp) :: start_vec(n_active_dims), diff
    integer  :: qqq, n_periods, ttt, aaa
    real(dp) :: interp_input_mat(n_states, n_interp), interp_input_monetary_mat(n_states, n_interp), & 
                state_variables(smolyak_d), other_variables(n_interp)

    real(dp) ::  rcond
    character(1) :: equed
    character :: shock_char*1, iter_char*4 

    real(dp) :: current_state(1,n_active_dims), next_state_tmp(n_quad,n_active_dims), next_state(1,n_active_dims), & 
                stochsteady_state(1,n_active_dims), checker , p_dis, big_weight_vec(n_quad)
    integer  :: q_nxt, sample_length, iii, iter


    integer, allocatable  :: starting_locs(:), sample_vector(:)
    real(dp), allocatable :: starting_states(:,:), state_collection(:,:), &
                             interp_monetary_coeffs(:,:), smol_monetary_coeffs(:,:,:), smol_montransition_coeffs(:,:,:), &
                             state_series_N(:,:,:), other_vars_series_N(:,:,:), shock_series_N(:,:,:), irf_factor(:), &
                             state_series_N_helper(:,:,:), other_vars_series_N_helper(:,:,:), shock_series_N_helper(:,:,:), &
                             state_series_M(:,:,:), other_vars_series_M(:,:,:), shock_series_M(:,:,:)


    narg=command_argument_count()
    if (narg>0) then
        call get_command_argument(1,command_input,command_len)
    endif
 
    lda  = n_states; ldaf = n_states 
    ldb  = n_states; ldx  = n_states
    
    write(*,*)
    write(*,*) 'READING RESULTS'

    open (unit = 10, file = trim(results_folder) // 'next_state_mat.dat', ACTION="read", STATUS="replace", &
            & FORM="unformatted", ACCESS="STREAM")
    read(10) next_state_mat; close(10)

    open (unit = 10, file = trim(results_folder) // 'results_mat.dat', ACTION="read", STATUS="replace", &
            & FORM="unformatted", ACCESS="STREAM")
    read(10) interp_input_mat; close(10)

    open (unit = 10, file = trim(results_folder) // 'next_state.txt', ACTION="write", STATUS="replace", &
            & FORM="formatted", ACCESS="STREAM")
    write(10,*) next_state_mat; close(10)

    open (unit = 10, file = trim(results_folder) // 'smol_grid.txt', ACTION="write", STATUS="replace", &
            & FORM="formatted", ACCESS="STREAM")
    write(10,*) smol_grid; close(10)

    open (unit = 10, file = trim(results_folder) // 'interp_input.txt', ACTION="write", STATUS="replace", &
            & FORM="formatted", ACCESS="STREAM")
    write(10,*) interp_input_mat; close(10)



    ! ------------------------------------------------------------ !
    ! Store Basic Grid Sizes, etc.
    ! ------------------------------------------------------------ !
    open (unit = 10, file = trim(results_folder) // 'num_params.csv', ACTION="write",  &
            & FORM="formatted", ACCESS="sequential")
    write(10,'(11i6)')  n_I, n_states, n_active_dims, n_interp, n_shocks, n_spread,  & 
                       smolyak_d, n_prm, n_sim_periods, n_irf_periods
    close(10)

    open (unit = 10, file = trim(results_folder) // 'grid_locs.csv', ACTION="write",  &
            & FORM="formatted", ACCESS="sequential")
    write(10,'(12f10.4)')  k_grid_mean, k_grid_dev, s_a_grid_mean, s_a_grid_dev, &
                          s_c_grid_mean, s_c_grid_dev, w_grid_mean, w_grid_dev, & 
                          dis_grid_mean, dis_grid_dev, m_grid_mean, m_grid_dev
    close(10)

    ! ------------------------------------------------------------ !
    ! Find Stochastic Steady State
    ! ------------------------------------------------------------ !

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

    write(*,*)
    write(*,*) 'FIND STOCHASTIC STEADY STATE'
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

    CALL DGEMM('N','N', 1, nrhs2, n_states, 1.0_dp, polyn_points, 1, & 
                           interp_coeffs, n_states, 0.0_dp, & 
                           other_variables, 1)   

    CALL DGEMM('N','N', 1, nrhs3, n_states, 1.0_dp, polyn_points, 1, & 
                           state_coeffs, n_states, 0.0_dp, & 
                           state_variables, 1)   

    ! (2.2) Print the Stochastic SS 
    write(*,*) 'Stochastic steady state -------------------------------------' 
    write(*,*) '-------------------------------------------------------------' 

    write(*,"(A12,F10.4)") 'Capital   :', state_variables(idx_k)
    write(*,"(A12,F10.4)") 'Wealth   a:', state_variables(idx_sa) 
    write(*,"(A12,F10.4)") 'Wealth   c:', state_variables(idx_sc) 
    write(*,"(A12,F10.4)") 'm         :', state_variables(idx_m) 
    write(*,"(A12,F10.4)") 'w         :', state_variables(idx_w) 
    write(*,"(A12,F10.4)") 'Disaster p:', state_variables(idx_dis)

    ! ------------------------------------------------------------ !
    ! Policy and Price ~ State
    ! ------------------------------------------------------------ !
    write(*,*) 
    write(*,*) 'State dependencies ------------------------------------------' 
    write(*,*) '-------------------------------------------------------------' 

    allocate(state_spread_series(n_spread,n_interp,smolyak_d))
    call Fill_linspace_dp(-1.0_dp, 1.0_dp, state_spread_vec)

    counter = 0
    do qqq = 1,smolyak_d
        if (vector_mus_dimensions(qqq) > 0) then 
        counter = counter + 1
        do ttt = 1,n_spread
            do aaa = 1,n_active_dims
                if (aaa == counter) then
                    current_state(1,aaa) = state_spread_vec(ttt)
                else
                    current_state(1,aaa) = stochsteady_state(1,aaa)
                endif
            enddo

            polyn_points = Smolyak_Polynomial(current_state, n_active_dims, &
                                max_smol_level, smol_elem_ani) 
             
            CALL DGEMM('N','N', 1, nrhs2, n_states, 1.0_dp, polyn_points, 1, & 
                               interp_coeffs, n_states, 0.0_dp, & 
                               state_spread_series(ttt,:,qqq), 1)
        enddo
        endif
    enddo

    open (unit = 10, file = trim(results_folder) // 'state_spread_series.txt', ACTION="write", STATUS="replace", &
            & FORM="formatted", ACCESS="STREAM")
    write(10,*) state_spread_series; close(10)


    n_periods = n_sim_periods + n_burn
    allocate( simshock_series(n_periods), shock_series(n_periods, n_shocks), state_series(n_periods, smolyak_d), other_vars_series(n_periods, n_interp) )
    allocate(state_collection(n_periods, n_active_dims))
    
    ! ------------------------------------------------------------ !
    ! Business Cycle Moments, No Disaster
    ! ------------------------------------------------------------ !
    write(*,*)
    write(*,*) 'Business cycle moments  no disaster -------------------------' 
    write(*,*) '-------------------------------------------------------------' 



    ! seed random number generator
    LSTATE  = 633
    SEED    = 712
    ifail   = 0
    call G05KFF( GENID, SUBID, SEED, LSEED, STATE, LSTATE, IFAIL)

    ! generate random numbers:
    ifail = 0
    call g05saf( n_periods, STATE, simshock_series, IFAIL)
    
    current_state = stochsteady_state 

    shock_series(1,:) = shock_grid(no_shock_idx(1),:)

    ! loop through
    do ttt = 1,n_periods
        
        polyn_points = Smolyak_Polynomial(current_state, n_active_dims, &
                            max_smol_level, smol_elem_ani) 

        CALL DGEMM('N','N', 1, nrhs2, n_states, 1.0_dp, polyn_points, 1, & 
                           interp_coeffs, n_states, 0.0_dp, & 
                           other_vars_series(ttt,:), 1)   
        


        CALL DGEMM('N','N', 1, nrhs3, n_states, 1.0_dp, polyn_points, 1, & 
                           state_coeffs, n_states, 0.0_dp, & 
                           state_series(ttt,:), 1)   

        checker = 0.0_dp
        q_nxt   = 0
        do while (checker <= simshock_series(ttt) .and. q_nxt < n_quad-1) 
            q_nxt   = q_nxt + 1
            checker = checker + quad_weight_vec(q_nxt)
        enddo
        
        
        CALL DGEMM('N','N', 1, nrhs, n_states, 1.0_dp, polyn_points, 1, & 
                           smol_coeffs(q_nxt,:,:), n_states, 0.0_dp, & 
                           next_state, 1)   
            
       
        if (ttt<n_periods) then 
        shock_series(ttt+1,:) = shock_grid(q_nxt,:)
        endif
        
        where (next_state > 1.0_dp)
            next_state = 1.0_dp 
        elsewhere (next_state < -1.0_dp)
            next_state = -1.0_dp 
        endwhere

        state_collection(ttt,:)   = current_state(1,:)
        current_state             = next_state
        
    enddo
    
    ! interpolation of MPKs of constrained agent does not work well due to sparse grid
    ! adjust manually
    where (other_vars_series(:,67) > 0.5) 
        other_vars_series(:,36) = 0.0_dp
        other_vars_series(:,45) = 0.0_dp
        other_vars_series(:,54) = 0.0_dp
        other_vars_series(:,57) = lmbd_vec(3)*kbar_constr(3)
    endwhere 

    open (unit = 10, file = trim(results_folder) // 'sim_state_series.txt', ACTION="write", STATUS="replace", &
            & FORM="formatted", ACCESS="STREAM")
    write(10,*) state_series(n_burn+1:n_periods,:); close(10)
   !   
    open (unit = 10, file = trim(results_folder) // 'sim_vars_series.txt', ACTION="write", STATUS="replace", &
            & FORM="formatted", ACCESS="STREAM")
    write(10,*) other_vars_series(n_burn+1:n_periods,:); close(10)

    open (unit = 10, file = trim(results_folder) // 'sim_shock_series.txt', ACTION="write", STATUS="replace", &
    & FORM="formatted", ACCESS="STREAM")
    write(10,*) shock_series(n_burn+1:n_periods,:); close(10)

    ! ------------------------------------------------------------ !
    ! Sample n_sample IRF starting points
    ! ------------------------------------------------------------ !
    write(*,*)
    write(*,*) 'Sampling IRF Starting Points -------------------------' 
    write(*,*) '-------------------------------------------------------------' 

    ! Generate random integers to sample starting points
    ifail = 0
    allocate(starting_locs(n_sample), starting_states(n_sample, n_active_dims))

    sample_length = n_sim_periods - n_burn
    allocate(sample_vector(sample_length))
    sample_vector = (/(ttt, ttt=n_burn+1,n_sim_periods, 1)/) 
    call g05ndf(sample_vector, sample_length, starting_locs, n_sample, STATE, IFAIL)
    ! Obtain starting states
    if (n_sample > 1) then  
    starting_states = state_collection(starting_locs,:)
    else    
    starting_states = stochsteady_state 
    endif
    
    open (unit = 10, file = trim(results_folder) // 'starting_states.dat', ACTION="write", STATUS="replace", &
    & FORM="unformatted", ACCESS="STREAM")
    write(10) starting_states; close(10)

    ! ------------------------------------------------------------ !
    ! Business Cycle Moments, With Disaster
    ! ------------------------------------------------------------ !
    write(*,*)
    write(*,*) 'Business cycle moments  w. disaster -------------------------' 
    write(*,*) '-------------------------------------------------------------' 

    ! generate random numbers:
    ifail = 0
    call g05saf( n_periods, STATE, simshock_series, IFAIL)
    
    current_state = stochsteady_state 

    shock_series(1,:) = shock_grid(no_shock_idx(1),:)

    ! loop through
    do ttt = 1,n_periods
        
        polyn_points = Smolyak_Polynomial(current_state, n_active_dims, &
                            max_smol_level, smol_elem_ani) 

        CALL DGEMM('N','N', 1, nrhs2, n_states, 1.0_dp, polyn_points, 1, & 
                           interp_coeffs, n_states, 0.0_dp, & 
                           other_vars_series(ttt,:), 1)   
        
        CALL DGEMM('N','N', 1, nrhs3, n_states, 1.0_dp, polyn_points, 1, & 
                           state_coeffs, n_states, 0.0_dp, & 
                           state_series(ttt,:), 1)   
        
        p_dis = exp(state_series(ttt,idx_dis)) 

        ! current weight vec
        ! CHANGED: Only last state is disaster
        big_weight_vec(1:n_quad-1) = quad_weight_vec*(1.0-p_dis)
        big_weight_vec(n_quad)     = p_dis

        checker = 0.0_dp
        q_nxt   = 0
        do while (checker <= simshock_series(ttt) .and. q_nxt < n_quad) 
            q_nxt   = q_nxt + 1
            checker = checker + big_weight_vec(q_nxt) 
        enddo
        
        
        CALL DGEMM('N','N', 1, nrhs, n_states, 1.0_dp, polyn_points, 1, & 
                           smol_coeffs(q_nxt,:,:), n_states, 0.0_dp, & 
                           next_state, 1)   
            
       
        if (ttt<n_periods) then 
        shock_series(ttt+1,:) = shock_grid(q_nxt,:)
        endif
        
        where (next_state > 1.0_dp)
            next_state = 1.0_dp 
        elsewhere (next_state < -1.0_dp)
            next_state = -1.0_dp 
        endwhere

        current_state        = next_state
        
    enddo
    
    ! interpolation of MPKs of constrained agent does not work well due to sparse grid
    ! adjust manually
    where (other_vars_series(:,67) > 0.5) 
        other_vars_series(:,36) = 0.0_dp
        other_vars_series(:,45) = 0.0_dp
        other_vars_series(:,54) = 0.0_dp
        other_vars_series(:,57) = lmbd_vec(3)*kbar_constr(3)
    endwhere 

    open (unit = 10, file = trim(results_folder) // 'simDis_state_series.txt', ACTION="write", STATUS="replace", &
            & FORM="formatted", ACCESS="STREAM")
    write(10,*) state_series(n_burn+1:n_periods,:); close(10)
   !   
    open (unit = 10, file = trim(results_folder) // 'simDis_vars_series.txt', ACTION="write", STATUS="replace", &
            & FORM="formatted", ACCESS="STREAM")
    write(10,*) other_vars_series(n_burn+1:n_periods,:); close(10)

    open (unit = 10, file = trim(results_folder) // 'simDis_shock_series.txt', ACTION="write", STATUS="replace", &
    & FORM="formatted", ACCESS="STREAM")
    write(10,*) shock_series(n_burn+1:n_periods,:); close(10)




    ! allocate various matrices

    n_periods = n_irf_periods

    allocate(state_series_N(n_sample, n_periods, smolyak_d), other_vars_series_N(n_sample, n_periods, n_interp), shock_series_N(n_sample, n_periods, n_shocks) )
    allocate(state_series_N_helper(n_sample, n_periods, smolyak_d), other_vars_series_N_helper(n_sample, n_periods, n_interp), shock_series_N_helper(n_sample, n_periods, n_shocks), irf_factor(n_sample))
    allocate(state_series_M(n_sample, n_periods, smolyak_d), other_vars_series_M(n_sample, n_periods, n_interp), shock_series_M(n_sample, n_periods, n_shocks) )
    deallocate(state_series, other_vars_series, shock_series )
    allocate( state_series(n_periods, smolyak_d), other_vars_series(n_periods, n_interp), shock_series(n_periods, n_shocks), other_vars_series_tmp(n_periods, n_interp) )

    write(*,*) "Calc IRF paths absent any shock"
    do iii = 1,n_sample
            
        if (mod(iii, 50) == 0) then
            write(*,*) iii
        endif

        ! Starting State
        current_state(1,:) = starting_states(iii,:) 
        
        ! loop through
        do ttt = 1,n_periods
              
             shock_series(ttt,:) = shock_grid(no_shock_idx(1),:)

             polyn_points = Smolyak_Polynomial(current_state, n_active_dims, &
                                 max_smol_level, smol_elem_ani) 
             
             CALL DGEMM('N','N', 1, nrhs, n_states, 1.0_dp, polyn_points, 1, & 
                                smol_coeffs(no_shock_idx,:,:), n_states, 0.0_dp, & 
                                next_state, 1)   
                 
             CALL DGEMM('N','N', 1, nrhs2, n_states, 1.0_dp, polyn_points, 1, & 
                                interp_coeffs, n_states, 0.0_dp, & 
                                other_vars_series(ttt,:), 1)   
             
             CALL DGEMM('N','N', 1, nrhs3, n_states, 1.0_dp, polyn_points, 1, & 
                                state_coeffs, n_states, 0.0_dp, & 
                                state_series(ttt,:), 1)   
             
             current_state        = next_state
             
        enddo
    
        ! interpolation of MPKs of constrained agent does not work well due to sparse grid
        ! adjust manually
        where (other_vars_series(:,67) > 0.5) 
            other_vars_series(:,36) = 0.0_dp
            other_vars_series(:,45) = 0.0_dp
            other_vars_series(:,54) = 0.0_dp
            other_vars_series(:,57) = lmbd_vec(3)*kbar_constr(3)
        endwhere 
        
        ! Store
        state_series_N(iii,:,:) = state_series 
        other_vars_series_N(iii,:,:) = other_vars_series
        shock_series_N(iii,:,:) = shock_series

        do ttt = 1,n_periods
        state_series_N_helper(iii,ttt,:) = state_series(1,:) 
        other_vars_series_N_helper(iii,ttt,:) = other_vars_series(1,:)
        shock_series_N_helper(iii,ttt,:) = shock_series(1,:)
        enddo

    enddo

    state_series      = sum(state_series_N     ,1)/n_sample 
    other_vars_series = sum(other_vars_series_N,1)/n_sample 
    shock_series      = sum(shock_series_N     ,1)/n_sample 
   

    open (unit = 10, file = trim(results_folder) // 'none_irf_state_series.txt', ACTION="write", STATUS="replace", &
            & FORM="formatted", ACCESS="STREAM")
    write(10,*) state_series; close(10)
      
    open (unit = 10, file = trim(results_folder) // 'none_irf_vars_series.txt', ACTION="write", STATUS="replace", &
            & FORM="formatted", ACCESS="STREAM")
    write(10,*) other_vars_series; close(10)

    open (unit = 10, file = trim(results_folder) // 'none_irf_shock_series.txt', ACTION="write", STATUS="replace", &
            & FORM="formatted", ACCESS="STREAM")
    write(10,*) shock_series; close(10)
    write(*,*) "DONE"
    write(*,*) ""


    write(*,*) "Calc IRF paths with mit shock"
    do fff = 1,size(irf_indices,1)

        ! ------------------------------------------------------------ !
        ! Loop Each IRF Calculation
        ! ------------------------------------------------------------ !
        write(shock_char,'(i1)') fff

        write(*,*)
        write(*,*) 'Computing Impulse responses -------------------------' 
        write(*,*) shock_char
        write(*,*) '-------------------------------------------------------------' 


        ! Get all solution files
        open (unit = 10, file = trim(results_folder) // 'next_state_shock_mat_' // shock_char // '.dat', ACTION="read", STATUS="replace", &
              & FORM="unformatted", ACCESS="STREAM")
        read(10) next_state_monetary_mat; close(10)

        ! Get all solution files
        open (unit = 10, file = trim(results_folder) // 'next_state_shocktrans_mat_' // shock_char // '.dat', ACTION="read", STATUS="replace", &
              & FORM="unformatted", ACCESS="STREAM")
        read(10) next_state_montransition_mat; close(10)

        ! print all results
        open (unit = 10, file = trim(results_folder) // 'results_shock_mat_' // shock_char // '.dat', ACTION="read", STATUS="replace", &
              & FORM="unformatted", ACCESS="STREAM")
        read(10) interp_input_monetary_mat; close(10)

        ! polynomial coefficients for transition matrix
        do qqq = 1,n_quad

            b = transpose(next_state_monetary_mat(qqq,:,:))

            call F07ABF( 'E','N',n_states, nrhs, smol_polynom, lda,af,ldaf,ipiv,      &
                         equed, r, c, b, ldb, smol_monetary_coeffs(qqq,:,:), & 
                         ldx, rcond, ferr, berr, work, iwork, info)

        enddo

        b = transpose(next_state_montransition_mat(1,:,:))

        call F07ABF( 'E','N',n_states, nrhs, smol_polynom, lda,af,ldaf,ipiv,      &
                     equed, r, c, b, ldb, smol_montransition_coeffs(1,:,:), & 
                     ldx, rcond, ferr, berr, work, iwork, info)

     
        b2 = interp_input_monetary_mat
         
        call F07ABF( 'E','N',n_states, nrhs2, smol_polynom, lda,af,ldaf,ipiv,      &
                     equed, r, c, b2, ldb, interp_monetary_coeffs, & 
                     ldx, rcond, ferr2, berr2, work, iwork, info )

        do iii = 1,n_sample
            
            if (mod(iii, 50) == 0) then
                write(*,*) iii
            endif

            ! ************************************************************ !
            ! IRF With Shock
            ! ************************************************************ !

            ! Starting State
            current_state(1,:) = starting_states(iii,:) 
            polyn_points = Smolyak_Polynomial(current_state, n_active_dims, &
                                    max_smol_level, smol_elem_ani) 

            ! state and variables in stochastic steady state
            CALL DGEMM('N','N', 1, nrhs2, n_states, 1.0_dp, polyn_points, 1, & 
                                   interp_coeffs, n_states, 0.0_dp, & 
                                   other_vars_series(1,:), 1)   

            CALL DGEMM('N','N', 1, nrhs3, n_states, 1.0_dp, polyn_points, 1, & 
                                   state_coeffs, n_states, 0.0_dp, & 
                                   state_series(1,:), 1)   
            
            CALL DGEMM('N','N', 1, nrhs, n_states, 1.0_dp, polyn_points, 1, & 
                                   smol_coeffs(no_shock_idx,:,:), n_states, 0.0_dp, & 
                                   next_state, 1)   
            
            current_state = next_state
            
            polyn_points = Smolyak_Polynomial(current_state, n_active_dims, &
                                    max_smol_level, smol_elem_ani) 

            ! state and variables in stochastic steady state
            CALL DGEMM('N','N', 1, nrhs2, n_states, 1.0_dp, polyn_points, 1, & 
                                   interp_coeffs, n_states, 0.0_dp, & 
                                   other_vars_series(2,:), 1)   

            CALL DGEMM('N','N', 1, nrhs3, n_states, 1.0_dp, polyn_points, 1, & 
                                   state_coeffs, n_states, 0.0_dp, & 
                                   state_series(2,:), 1)   

            ! transition to shocked state
            CALL DGEMM('N','N', 1, nrhs, n_states, 1.0_dp, polyn_points, 1, & 
                                   smol_montransition_coeffs(1,:,:), n_states, 0.0_dp, & 
                                   next_state, 1)   
            
            ! in that state find state and variables
            current_state = next_state 
            polyn_points = Smolyak_Polynomial(current_state, n_active_dims, &
                                    max_smol_level, smol_elem_ani) 
            
            CALL DGEMM('N','N', 1, nrhs2, n_states, 1.0_dp, polyn_points, 1, & 
                                   interp_monetary_coeffs, n_states, 0.0_dp, & 
                                   other_vars_series(3,:), 1)   
            
            CALL DGEMM('N','N', 1, nrhs3, n_states, 1.0_dp, polyn_points, 1, & 
                                   state_coeffs, n_states, 0.0_dp, & 
                                   state_series(3,:), 1)   

            ! transition out of this state
            CALL DGEMM('N','N', 1, nrhs, n_states, 1.0_dp, polyn_points, 1, & 
                                   smol_monetary_coeffs(no_shock_idx,:,:), n_states, 0.0_dp, & 
                                   next_state, 1)   
            
            shock_series(1,:) = shock_grid(no_shock_idx(1),:)
            shock_series(2,:) = shock_grid(no_shock_idx(1),:)
            shock_series(3,:) = shock_grid(no_shock_idx(1),:)

            if (irf_indices(fff) == sidx_z) then
                shock_series(3,sidx_z) = irf_shock_sizes(fff)
            endif

            current_state = next_state

            ! loop through
            do ttt = 4,n_periods
                 
                if (ttt<=n_periods) then 
                    shock_series(ttt,:) = shock_grid(no_shock_idx(1),:)
                endif

                polyn_points = Smolyak_Polynomial(current_state, n_active_dims, &
                                    max_smol_level, smol_elem_ani) 
                
                CALL DGEMM('N','N', 1, nrhs, n_states, 1.0_dp, polyn_points, 1, & 
                                   smol_coeffs(no_shock_idx,:,:), n_states, 0.0_dp, & 
                                   next_state, 1)   
                    
                CALL DGEMM('N','N', 1, nrhs2, n_states, 1.0_dp, polyn_points, 1, & 
                                   interp_coeffs, n_states, 0.0_dp, & 
                                   other_vars_series(ttt,:), 1)   
                
                CALL DGEMM('N','N', 1, nrhs3, n_states, 1.0_dp, polyn_points, 1, & 
                                   state_coeffs, n_states, 0.0_dp, & 
                                   state_series(ttt,:), 1)   
                
                current_state        = next_state
                
           enddo

         
           if (sidx_m == irf_indices(fff)) then
               irf_factor(iii) =  0.002/((log(other_vars_series(3,jx_1y)) - log(other_vars_series_N(iii,3,jx_1y))) - (log(other_vars_series(2,jx_1y)) - log(other_vars_series_N(iii,2,jx_1y))))
           else 
                irf_factor(iii)  = 1.0_dp
           endif

            if (sidx_m == irf_indices(fff) .and. fff == 2) then  
                open (unit = 10, file = trim(results_folder) // 'irf_factor.txt', ACTION="write", STATUS="replace", &
                        & FORM="formatted", ACCESS="STREAM")
                write(10,*) sum(irf_factor)/n_sample; close(10)
            endif
        
            ! interpolation of MPKs of constrained agent does not work well due to sparse grid
            ! assign heuristically
            where (other_vars_series(:,67) > 0.5) 
                other_vars_series(:,36) = 0.0_dp
                other_vars_series(:,45) = 0.0_dp
                other_vars_series(:,54) = 0.0_dp
                other_vars_series(:,57) = lmbd_vec(3)*kbar_constr(3)
            endwhere 
        
            ! Store
            state_series_M(iii,:,:)      = (state_series - state_series_N(iii,:,:))*irf_factor(iii)          + state_series_N_helper(iii,:,:)
            other_vars_series_M(iii,:,:) = (other_vars_series- other_vars_series_N(iii,:,:))*irf_factor(iii) + other_vars_series_N_helper(iii,:,:)
            shock_series_M(iii,:,:)      = (shock_series - shock_series_N(iii,:,:))*irf_factor(iii)          + shock_series_N_helper(iii,:,:)

        enddo

        ! Take Averages across sample points, Shocked
        state_series      = sum(state_series_M      ,1)/n_sample 
        other_vars_series = sum(other_vars_series_M ,1)/n_sample 
        shock_series      = sum(shock_series_M      ,1)/n_sample 

        if (irf_indices(fff) == sidx_m .and. fff == 2) then 
            open (unit = 10, file = trim(results_folder) // 'm_irf_state_series.txt', ACTION="write", STATUS="replace", &
                    & FORM="formatted", ACCESS="STREAM")
            write(10,*) state_series; close(10)
            
            open (unit = 10, file = trim(results_folder) // 'm_irf_vars_series.txt', ACTION="write", STATUS="replace", &
                    & FORM="formatted", ACCESS="STREAM")
            write(10,*) other_vars_series; close(10)

            open (unit = 10, file = trim(results_folder) // 'm_irf_shock_series.txt', ACTION="write", STATUS="replace", &
                    & FORM="formatted", ACCESS="STREAM")
            write(10,*) shock_series; close(10)
            
        elseif (irf_indices(fff) == sidx_z) then 

            open (unit = 10, file = trim(results_folder) // 'g_irf_state_series.txt', ACTION="write", STATUS="replace", &
                    & FORM="formatted", ACCESS="STREAM")
            write(10,*) state_series; close(10)
              
            open (unit = 10, file = trim(results_folder) // 'g_irf_vars_series.txt', ACTION="write", STATUS="replace", &
                    & FORM="formatted", ACCESS="STREAM")
            write(10,*) other_vars_series; close(10)

            open (unit = 10, file = trim(results_folder) // 'g_irf_shock_series.txt', ACTION="write", STATUS="replace", &
                    & FORM="formatted", ACCESS="STREAM")
            write(10,*) shock_series; close(10)

        elseif (irf_indices(fff) == sidx_dis) then 
            open (unit = 10, file = trim(results_folder) // 'p_irf_state_series.txt', ACTION="write", STATUS="replace", &
                    & FORM="formatted", ACCESS="STREAM")
            write(10,*) state_series; close(10)
              
            open (unit = 10, file = trim(results_folder) // 'p_irf_vars_series.txt', ACTION="write", STATUS="replace", &
                    & FORM="formatted", ACCESS="STREAM")
            write(10,*) other_vars_series; close(10)

            open (unit = 10, file = trim(results_folder) // 'p_irf_shock_series.txt', ACTION="write", STATUS="replace", &
                    & FORM="formatted", ACCESS="STREAM")
            write(10,*) shock_series; close(10)
        endif

    enddo



end subroutine create_results



end module mod_results


