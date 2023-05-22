! -------------------------------------------------------------------------
! main.f90: main program 
! -------------------------------------------------------------------------
! authors:         Rohan Kekre and Moritz Lenel
! for updates see: https://github.com/KekreLenel/MPR
! -------------------------------------------------------------------------
program main
    
    use mod_param,   only: init_setup , grid_setup  
    use mod_calc,    only: calc_steady, calc_sol
    use mod_results, only: create_results
    use mod_decomp,  only: monetary_decomposition

    implicit none

    call init_setup()               ! load parameters 
    call calc_steady()              ! calculate steady state
    call grid_setup()               ! prepare Smolyak grid 
    call calc_sol()                 ! solve model 
    call create_results()           ! create result files
    call monetary_decomposition()   ! additional decomposition results

end program main
