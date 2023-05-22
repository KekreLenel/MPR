! This module provides a Fortran translation of the codes
! written by Kenneth L. Judd, Lilia Maliar, Serguei Maliar 
! and Rafael Valero for their paper (2014),  "Smolyak method 
! for solving dynamic economic models: Lagrange interpolation,  
! anisotropic grid and adaptive domain" Journal of Economic Dynamics  
! and Control 44, 92 123 

! Copyright 2014 by Lilia Maliar, Serguei Maliar and Rafael Valero. All 
! rights reserved. The code may be used, modified and redistributed under  
! the terms provided in the license agreement at the bottom of this file as
! provided by Kenneth L. Judd, Lilia Maliar, Serguei Maliar and Rafael Valero
! See: https://lmaliar.ws.gc.cuny.edu/codes/

! this translation to Fortran: Rohan Kekre and Moritz Lenel, March 2022
! for updates see: https://github.com/KekreLenel/MPR

module mod_smolyak

use base_lib

implicit none
private

public  :: Smolyak_Elem_Isotrop, Smolyak_Elem_Anisotrop, Smolyak_Grid, Smolyak_Polynomial, Smolyak_Polynomial2

contains

! Smolyak_Elem_Isotrop is a routine that constructs the subindices of the
! Smolyak elements (grid points and basis functions) for the isotropic case;
! see "Smolyak method for solving dynamic economic models: Lagrange interpo-
! lation, anisotropic grid and adaptive domain" by Kenneth L. Judd, Lilia Maliar, 
! Serguei Maliar and Rafael Valero, (2014), Journal of Economic Dynamics and 
! Control 44, 92 123 (henceforth, JMMV (2014)), Section 2.2.3
!
! This version: Novenber 5, 2014. First version: December 17, 2012.
! This translation: February 2020
! -------------------------------------------------------------------------
! Input:   "d" is the number of dimensions (the number of state  variables) 
!          "mu" is the level of approximation
!  
! Output:  "Smolyak_elem_iso" is the vector of subindices of unidimensional 
!           elements (Smolyak grid points or polynomial basis function) that 
!           constitute a multidimensional element (Smolyak grid point or 
!           polynomial basis function) for the isotropic case

function Smolyak_Elem_Isotrop(d, mu) result(Smolyak_elem_iso)

integer, intent(in) 	:: d, mu
integer, allocatable 	:: Smolyak_elem_iso(:,:), Smol_rule(:,:), incr_Smol_rule(:,:), &
                           prev_incr(:,:), aux(:,:), aux_new(:,:), augmented(:,:), &
                           new_Smol_rule(:,:), one_comb(:), Smolyak_elem_iso_new(:,:), &
                           incr_indices(:,:), prev_indices(:,:), z(:,:), indices_elem_jd(:,:), & 
                           a(:,:), b(:,:), ones_vec(:,:)

integer                 :: i, jd, j, m, id, current_size, current_size2, n_comb, &
                           a_columns, b_columns, a_rows, b_rows, k
                           
                        

! 1. Identify the indices of disjoint sets A's, i1,...,id, that satisfy the 
! Smolyak rule, d<=i1+i2+...+id<=|i|; see equation (1) in JMMV(2014)
! -------------------------------------------------------------------------
   
!     Smol_rule = [];     ! This will be a matrix of subindices i1,...,id of 
!                         ! disjoint sets A_i's such that their sum across all  
!                         ! dimensions, |i|=i1+i2+...+id,  satisfies  the 
!                         ! Smolyak rule, d<=|i|<=d+mu; initially, the matrix 
!                         ! is empty; at the intermediate steps j=0,...,mu, 
!                         ! the subindices satisfy d<=|i|=d+j
                       
!     incr_Smol_rule = [];! The matrix of subindices i1,...,id of disjoint   
!                         ! sets A_i's, such that their sum across all dimensions,  
!                         ! |i|=i1+i2+...+id, is exactly equal to d+j, i.e., 
!                         ! |i|=d+j; initially, the matrix is empty; when j 
!                         ! increases, this matrix is concatinated to matrix  
!                         ! "Smol_rule" obtained for the previous value of j,
!                         ! i.e., j-1
do j = 0,mu  
       
    if (j > 0) then   
    prev_incr = incr_Smol_rule
    endif
    ! For the previous value of j, call "incr_Smol_rule"  
    ! as "prev_incr"
     
    ! Identify new subindices of unidimensional sets A's i1,i2,...,id that 
    ! jointly satisfy the Smolyak rule as the sum of subindices increases 
    ! from j-1 to j; see JMMV (2014), Section 2.2 
    ! ---------------------------------------------------------------------                
    if (j == 0) then

        allocate(incr_Smol_rule(1:1,1:d))
        incr_Smol_rule(1:1,1:d) = 1  
                                    ! The matrix of subindices for j=0 is a 
                                    ! 1-by-d vector of ones, (1,...,1)
    else  
        m = size(prev_incr,1)
                                    ! Compute the number of columns in the
                                    ! previous matrix of subindices that
                                    ! satisfy the Smolyak rule "prev_incr"

        ! incr_Smol_rule = [];      ! Initially, the matrix of subindices is 
        
        if (allocated(aux)) then; deallocate(aux); endif
        allocate(aux(m,d))
        aux = 0                     ! Allocate memory to an initial auxiliary
                                    ! matrix that will be added to
                                    ! "prev_incr"

        if (allocated(incr_Smol_rule)) then; deallocate(incr_Smol_rule); endif
        allocate(incr_Smol_rule(m*d,d))
        do id = 1,d

            aux_new = aux         	! New auxiliary matrix is equal to the old 
                                    ! one
            aux_new(:,id) = 1       ! For a dimension i, set elements of this
                                    ! new auxiliary matrix to 1 instead of 0


                                  
            augmented = prev_incr + aux_new
                                    ! Increase the subinices of
                                    ! "prevoius_incr" by 1 in dimension id
            
            current_size = size(incr_Smol_rule,1)

            incr_Smol_rule(1+(id-1)*m:id*m,1:d) = augmented                    
                                    ! Concatenate "incr_Smol_rule" and
                                    ! "augmented": the first row of "augmented"
                                    ! goes after the last row of
                                    ! "incr_Smol_rule" 

        enddo
        
    endif
    
    incr_Smol_rule = remove_dups_2d(incr_Smol_rule)
                            ! Eliminate the repeated indices in the
                            ! matrix "incr_Smol_rule"   
                                   
    ! Concatenate the matrix of newly constructed indices to the previously 
    ! obtained matrix of indices satisfying the Smolyak rule
    ! ---------------------------------------------------------------------
    if (allocated(Smol_rule)) then 
        current_size  = size(Smol_rule,1)
        current_size2 = size(incr_Smol_rule,1)
        
        if (allocated(new_Smol_rule)) then; deallocate(new_Smol_rule); endif
        allocate(new_Smol_rule(current_size+current_size2,d))

        new_Smol_rule(1:current_size,1:d) = Smol_rule
        new_Smol_rule(current_size+1:current_size+current_size2,1:d) = incr_Smol_rule 
        Smol_rule = new_Smol_rule   ! E.g., for mu=1 and d=2,
    else                            ! Smol_rule=[1 1; 1 2; 2 1]
        Smol_rule = incr_Smol_rule
    endif

enddo
  
n_comb = size(Smol_rule,1);
            ! The total number of combinations of indices, i1,...,id, 
            ! of the unidimensional disjoint sets A's that satisfy the 
            ! Smolyak rule; e.g., for mu=1 and d=2, n_comb=3 such as (1)
            ! i1=1 and i2=1, (2) i1=1 and i2=2, (3) i1=2 and i2=1

! Smol_rule = remove_dups_2d(Smol_rule)              
n_comb = size(Smol_rule,1);

! 2. Construct the multidimensional indices of elements as a Cartesian product
! of unidimensional indices
! -------------------------------------------------------------------------

! Smolyak_elem_iso = []; 
!                 ! The matrix of multidimensional indices of elements (points)
!                 !  belonging to the disjoint subsets A's that satisfy the 
!                 ! Smolyak rule (every singl   e point is indexed); initially, 
!                 ! the matrix is empty
              
do i = 1,n_comb            ! For each combination of subindices i1,...,id 
                           ! that satisfies the Smolyak rule

    if (allocated(incr_indices)) then 
        deallocate(incr_indices)
    endif                   ! This is a matrix of multidimensional indices
                            ! of unidimensional grid points that will be
                            ! added to the final multidimensional matrix of 
                            ! indices "Smolyak_elem_iso", as n_comb increases;
                            ! initially, the matrix is empty

    one_comb = Smol_rule(i,:) 
                            ! Consider the i-th row of "Smol_rule" (i.e.,
                            ! one particular combination of the subindices 
                            ! of the disjoint sets satisfying the Smolyak
                            ! rule); e.g., for mu=1 and d=2, Smol_rule=
                            ! [1 1; 1 2; 2 1] and one_comb for i=1 is [1 1];
                            ! 1-by-d

    do jd = 1,d             ! For each dimension jd, ...

        if (allocated(prev_indices)) then; deallocate(prev_indices); endif
        if (allocated(incr_indices)) then
            prev_indices = incr_indices   
        endif

        
            
        !  Compute the indices of elements (points) of the unidimensional 
        !  set
        ! ----------------------------------------------------------------
        if (one_comb(jd) == 1) then
            ! Take a jd-th element of the row vector "one_comb" 
            ! that corresponds to dimension jd (this is the 
            ! subindex of the unidimensional disjoint set A_i 
            ! from which this element comes from; if an element
            ! (point) is from the disjoint set
            ! A_1,... 
            if (allocated(indices_elem_jd)) then; deallocate(indices_elem_jd); endif
            allocate(indices_elem_jd(1:1,1:1))             
            indices_elem_jd = 1                        
            ! A_1 contains one element; this element is indexed "1"
   
        elseif (one_comb(jd) == 2) then
                        ! If an element (point) is from the disjoint set
                        ! A_2,...
            if (allocated(indices_elem_jd)) then; deallocate(indices_elem_jd); endif
            allocate(indices_elem_jd(2,1:1))
            indices_elem_jd(:,1) = [2, 2**(one_comb(jd)-1)+1]            
                        ! A_2 contains two elements; these elements are
                        ! indexed "2" and "3", so that indices_elem_jd=[2;3]
        else
            if (allocated(indices_elem_jd)) then; deallocate(indices_elem_jd); endif
            allocate(indices_elem_jd(1:1 + 2**(one_comb(jd)-1)+1 - 2**(one_comb(jd)-2)-2, 1:1))
            do id = 1, 1 + 2**(one_comb(jd)-1)+1 - 2**(one_comb(jd)-2)-2
                indices_elem_jd(id, 1) = id + 2**(one_comb(jd)-2)+2 - 1
            enddo
                        ! The subsequent disjoint sets contain the elements 
                        ! from m(one_comb(jd)-1)+1=2^(one_comb(jd)-2)+2 to 
                        ! m(one_comb(jd))=2^(one_comb(jd)-1)+1; e.g., the 
                        ! disjoint set A_4 contains the elements indexed 
                        ! 6, 7, 8 and 9, so that indices_elem_jd=[6,7,8,9]
        endif                      
        
        ! Create a Cartesian product of two sets, "prev_indices" and 
        ! "indices_elem_jd"
        !-----------------------------------------------------------
        

        if (allocated(a)) then; deallocate(a); endif
        if (allocated(prev_indices)) then  
        a = prev_indices         ! Call one set "a" 
        endif

        b = indices_elem_jd      ! Call the other set "b"
        
        if (.not.(allocated(a))) then ! If "a" is empty, ...
            z = b                ! The Cartesian product "z" is given by "b"   
        else
            if (allocated(z)) then; deallocate(z); endif
            ! For the case of non-empty sets, ...
            ! Initially, the Cartesian product is empty
            ! rows and column notatation switched due to legacy code
            a_columns = size(a,1)  
            b_columns = size(b,1)
            a_rows    = size(a,2)  
            b_rows    = size(b,2)
            allocate(z(a_columns*b_columns, a_rows+b_rows))

            do k = 1,b_columns
                z(1+a_columns*(k-1):a_columns*k,1:a_rows) = a 
 
                if (allocated(ones_vec)) then; deallocate(ones_vec); endif
                allocate(ones_vec(1:a_columns,1:1)); ones_vec = 1
                
                z(1+a_columns*(k-1):a_columns*k,a_rows+1:b_rows+a_rows) = matmul(ones_vec, b(k:k,1:b_rows))


            enddo
        endif
           
        incr_indices = z 

    enddo

    if (allocated(Smolyak_elem_iso_new)) then; deallocate(Smolyak_elem_iso_new); endif
    
    if (i == 1) then 
        Smolyak_elem_iso = incr_indices
    else 
        allocate(Smolyak_elem_iso_new(size(Smolyak_elem_iso,1)+size(incr_indices,1),size(Smolyak_elem_iso,2))) 
        Smolyak_elem_iso_new(1:size(Smolyak_elem_iso,1),:) = Smolyak_elem_iso
        Smolyak_elem_iso_new(size(Smolyak_elem_iso,1)+1:size(Smolyak_elem_iso,1)+size(incr_indices,1),:) = incr_indices 
                               ! Construct the matrix of multidimensional indices 
                               ! of elements by concatenating "incr_indices" to 
                               ! the previous "Smolyak_elem_iso" obtained for the 
                               ! previous combination "n_comb"

        Smolyak_elem_iso = Smolyak_elem_iso_new
    endif

enddo

end function Smolyak_Elem_Isotrop


! Smolyak_Elem_Anisotrop is a routine that selects a subset of the subindices 
! of Smolyak elements corresponding to the given anisotropic case from a set of 
! subindices of the Smolyak isotropic elements; see "Smolyak method for solving  
! dynamic economic models: Lagrange interpolation, anisotropic grid and   
! adaptive domain" by Kenneth L. Judd, Lilia Maliar, Serguei Maliar and  
! Rafael Valero, (2014), Journal of Economic Dynamics and Control 44, 92�123 
! (henceforth, JMMV (2014)), Section 2.2.3 
!
! This version: Novenber 5, 2014. First version: December 17, 2012.
! This translation: February 2020
! -------------------------------------------------------------------------
! Input:   "smol_elem_iso"         is the matrix of subindices of unidi- 
!                                  mensional elements that constitute multi-
!                                  dimensional elements for the isotropic case
!          "vector_mus_dimensions" is the vector of the levels of
!                                  approximations in all dimensions  
! Output:  "smol_elem_ani"         is the matrix of subindices of unidimen-
!                                  sional elements that constitute multidi-
!                                  mensional elements for the given anisotropic 
!                                  case
! 
! -------------------------------------------------------------------------
! Copyright 2014 by Lilia Maliar, Serguei Maliar and Rafael Valero. All 
! rights reserved. The code may be used, modified and redistributed under  
! the terms provided in the file "License_Agreement.txt".
! -------------------------------------------------------------------------


function  Smolyak_Elem_Anisotrop(smol_elem_iso, smolyak_d, vector_mus_dimensions) result(smol_elem_ani)

 integer, intent(in)    :: smol_elem_iso(:,:), smolyak_d, vector_mus_dimensions(smolyak_d)
 integer, allocatable   :: smol_elem_ani(:,:), smol_elem_ani_temp(:,:)

 integer                :: points_dimensions(smolyak_d), aux, i, counter


points_dimensions = 0   ! This vector will tell how many
                        ! unidimensional elements in
                        ! each dimension we consider

do i = 1,smolyak_d
    aux = vector_mus_dimensions(i)
    if (aux == 0) then                     ! If the approximation level in
                                           ! the i-th dimension is 0, ...
        points_dimensions(i) = 1           ! The number of unidimensional 
                                           ! elements is 1
    else                                   ! If the approximation level in
                                           ! the i-th dimension is not 0,...
        points_dimensions(i) = 2**(aux)+1  ! Compute the number of unidimensional
                                           ! elements using the formula
    endif
enddo

smol_elem_ani = smol_elem_iso

do i = 1,smolyak_d
    where (smol_elem_ani(:,i) > points_dimensions(i))
        ! If a subindex (i.e., number of elements) of the isotropic case is 
        ! larger than that of the anisotropic case, set an auxiliary variable
        ! "aux1" to 1
        smol_elem_ani(:,i) = 0
    endwhere
enddo

smol_elem_ani_temp = 0*smol_elem_iso

counter = 0
do i = 1,size(smol_elem_ani,1)
    if (.not.any(smol_elem_ani(i,:) == 0) ) then 
        counter = counter + 1
        smol_elem_ani_temp(counter,:) = smol_elem_ani(i,:)
    endif
enddo

smol_elem_ani = smol_elem_ani_temp(1:counter,:)

end function Smolyak_Elem_Anisotrop

! Smolyak_Grid is a routine that constructs a multidimensional Smolyak  
! grid in the hypercube [-1,1]^d; see "Smolyak method for solving dynamic 
! economic models: Lagrange interpolation, anisotropic grid and adaptive 
! domain" by Kenneth L. Judd, Lilia Maliar, Serguei Maliar and Rafael Valero,
! (2014), Journal of Economic Dynamics and Control 44, 92�123 (henceforth, 
! JMMV (2014)), Section 2.2.3 
!
! This version: Novenber 5, 2014. First version: December 17, 2012.
! -------------------------------------------------------------------------
! Input:   "d"         is the number of dimensions (the number of state 
!                      variables)
!          "mu"        is the level of approximation (in the anisotropic
!                      case, this is maximum level of approximation)
!          "smol_elem" is the matrix of the subindices of the Smolyak
!                      unidimensional elements; these elements can be either 
!                      isotropic (produced by Smolyak_Elem_Isotrop.m) or 
!                      anisotropic (produced by Smolyak_Elem_Anisotrop.m); 
!                      in the former case, the indices i1,...,id that jointly 
!                      satisfy the Smolyak rule, d<=|i|<=d+mu, where 
!                      |i|=i1+i2+...+id; see JMMV (2014), Section 3.2.3; 
!                      in the later case, they are a subset of the above 
!                      indices 
!
! Output:  "smol_grid" is the multidimensional Smolyak grid 
! -------------------------------------------------------------------------
! Copyright � 2014 by Lilia Maliar, Serguei Maliar and Rafael Valero. All 
! rights reserved. The code may be used, modified and redistributed under  
! the terms provided in the file "License_Agreement.txt".
! -------------------------------------------------------------------------


function  Smolyak_Grid(d,mu,smol_elem) result(smol_grid)

integer, intent(in)     :: d, mu, smol_elem(:,:)
real(dp), allocatable    :: smol_grid(:,:)

integer              :: i_max, i, m_i, j, numb_points_md, jp, index_row(d), jd

real(dp), allocatable :: extrem_Cheb_1d(:), points_1d(:), points_1d_temp(:)

! 1. Compute the vector of extrema of Chebyshev polynomials corresponding 
! to the given level of Smolyak approximation mu
! -----------------------------------------------------------------------

! These points will be ordered as in Section 2.2.1 of JMMV(2014); e.g., for
! mu=1, the set of points is {0,-1,1}
         
!points_1d = [];                   ! Initially, the set of unidimensional 
                                   ! points "points_1d" is empty; see JMMV  
                                   ! (2014), Section 2.2.1
i_max = mu+1                       ! The maximum subindex of unidimensional
                                   ! set A_i whose points are used to
                                   ! construct Smolyak grid of the given mu; 
                                   !  e.g., for mu=1, we consider up to 
                                   ! A_i_max={-1,1} where i_max=1+1=2
do i = 1,i_max                     ! A subindex of a unidimensional set of
                                   ! points                                

    ! Compute the number of elements, m(i),(using m(i)=2^(i-1)+1) in the  
    ! i-th unidimensional set of points; see Section 2.2.1 in JMMV (2014)
    !---------------------------------------------------------------------
    if (i==1) then 
        m_i = 1          ! If i=1, then m(i)=1
    else
        m_i = 2**(i-1)+1 ! If i>1, then m(i) = 2^(i-1)+1
    endif
     

     ! Construct the extrema of Chebyshev polynomials used as unidimensional 
     ! grid points in the Smolyak method
     !---------------------------------------------------------------------
     if (m_i==1) then
        extrem_Cheb_1d = [0]
     else 
        
        if (allocated(extrem_Cheb_1d)) then; deallocate(extrem_Cheb_1d); endif
        allocate(extrem_Cheb_1d(m_i))
        do j = 1,m_i
                                        ! For j=1,...,m_i,...
            extrem_Cheb_1d(j) = -cos(pi*real((j-1),dp)/real(m_i-1,dp))    
                                        ! Chebyshev polynomials are defined in 
                                        ! the interval [-1,1]

            if (abs(extrem_Cheb_1d(j))<1E-12_dp) then
                extrem_Cheb_1d(j) = 0.0 ! Round "extrem_Cheb_1d" to 0 if its  
                                        ! absolute value is smaller than 1d-12
            elseif (1.0-extrem_Cheb_1d(j)<1E-12_dp) then
                extrem_Cheb_1d(j) = 1.0  ! Round "extrem_Cheb_1d" to 1 if    
                                         ! 1-extrem_Cheb_1d is smaller than 1d-12
            elseif (1.0+extrem_Cheb_1d(j)<1E-12_dp) then
                extrem_Cheb_1d(j) = -1.0 ! Round "extrem_Cheb_1d" to -1 if   
                                         ! 1+extrem_Cheb_1d is smaller than 1d-12
            endif

        enddo

     endif 


    if (i>1) then
        points_1d_temp = points_1d
        deallocate(points_1d)
        allocate(points_1d(1:size(points_1d_temp) + size(extrem_Cheb_1d)))
        points_1d(1:size(points_1d_temp)) = points_1d_temp
        points_1d(size(points_1d_temp)+1:size(points_1d)) = extrem_Cheb_1d

    else 
        points_1d = extrem_Cheb_1d
    endif
    
    


                                   ! Add to the previous set "points_1d" new 
                                   ! points (given by extrema of unidimensional  
                                   ! Chebyshev polynomials) as i increases
    
    points_1d = remove_dups_1d_real(points_1d)
                                ! Choose the unrepeated points and order 
                                ! them as in Section 2.2.1 of JMMV (2014);
                                ! no reordering! 
                                   
    
enddo              


 ! 2. Construct the matrix multidimensional Smolyak grid points for the   
 ! required level of Smolyak approximation, mu; see JMMV (2014), Sections 2.2.3 
 ! for examples
 ! -------------------------------------------------------------------------
 allocate(smol_grid(size(smol_elem,1),size(smol_elem,2)))
 smol_grid = 0.0_dp     ! Initialize the matrix of multidimensional
                        ! Smolyak grid points
                               
numb_points_md = size(smol_grid,1);
                    ! Compute the number of multidimensional 
                    ! Smolyak grid points                                 
 
do jp = 1,numb_points_md     ! For each multidimensional grid point, ...
                                
    index_row = smol_elem(jp,:)
                                ! Identify the subindex of the unidimensional
                                ! grid point jp; this is a jp-th row of matrix 
                                ! "smol_elem"; 1-by-d
    
    do jd = 1,d                 ! For each dimension (state variable), ...
                                ! A subindex of a unidimensional grid point 
                                ! in a dimension jd is denoted n
        
         smol_grid(jp,jd) = points_1d(index_row(jd))
                                ! Find the corresponding unidimensional grid
                                ! point in the vector "points_1d"
     enddo

enddo

end function Smolyak_Grid


! Smolyak_Polynomial.m is a routine that constructs the multidimensional 
! basis functions of Smolyak polynomial of the approximation level that  
! corresponds to the previously constructed Smolyak (multidimensional) grid 
! points; see "Smolyak method for solving dynamic economic models: Lagrange 
! interpolation, anisotropic grid and adaptive domain" by Kenneth L. Judd, 
! Lilia Maliar, Serguei Maliar and Rafael, (2014). Journal of Economic 
! Dynamics and Control 44, 92�123 (henceforth, JMMV (2014)). 
!
! This version: Novenber 5, 2014. First version: May 30, 2011.
! -------------------------------------------------------------------------
! Inputs:  "points"    is the matrix of points in which the polynomial basis 
!                      functions must be evaluated; numb_pts-by-d
!          "d"         is the number of dimensions (state variables)
!          "smol_elem" is the matrix of the subindices of the Smolyak
!                      unidimensional elements; these elements can be either 
!                      isotropic (produced by Smolyak_Elem_Isotrop.m) or 
!                      anisotropic (produced by Smolyak_Elem_Anisotrop.m); 
!                      in the former case, the indices i1,...,id that jointly 
!                      satisfy the Smolyak rule, d<=|i|<=d+mu, where 
!                      |i|=i1+i2+...+id; see JMMV (2014), Section 3.2.3; 
!                      in the later case, they are a subset of the above 
!                      indices 
!
! Output:  "Smol_bases" is the matrix of multidimensional basis functions of 
!                      Smolyak polynomial of the given level of approximation, 
!                      evaluated in data matrix "points"
! -------------------------------------------------------------------------
! Copyright 2014 by Lilia Maliar, Serguei Maliar and Rafael Valero. All 
! rights reserved. The code may be used, modified and redistributed under  
! the terms provided in the file "License_Agreement.txt".
! -------------------------------------------------------------------------


function  Smolyak_Polynomial(points,d,mu,smol_elem) result(smol_bases)
 
integer, intent(in)  :: d, mu, smol_elem(:,:)
real(dp), intent(in) :: points(:,:)
real(dp), allocatable :: smol_bases(:,:), smol_bases_temp(:,:)

integer  :: i_max, numb_pts, numb_terms, jt, jd, n, m_i_max, j
real(dp), allocatable :: phi(:,:,:), pproduct(:)
integer, allocatable :: index_row(:)

! Smolyak polynomial is given by the sum of multidimensional basis functions, 
! multiplied by the coefficients; see formula (15) in JMMV (2014). By 
! convention, the first basis function is given by 1 (unity). 

! Unidimensional basis functions are given by Chebyshev polynomial bases; 
! in JMMV (2014), a unidimensional Chebyshev polynomial basis function of   
! degree n-1 is denoted by "phi_n", i.e., has a subindex n and we follow  
! this notation here

! 1. Construct the unidimensional basis functions and evaluate them in
! all the points of matrix "points"
! -------------------------------------------------------------------------
i_max = mu+1                    ! The maximum subindex of unidimensional
                                ! set S_i whose points are used to
                                ! construct Smolyak grid; e.g., for mu=1, 
                                ! we consider up to S_i_max={0,-1,1} where
                                ! i_max=mu+1=1+1=2
                                
! ! Compute the number of elements in the i_max-th unidimensional set of 
! ! elements S_i_max; this coincides with a maximum subindex of elements 
! ! (unidimensional grid point or unidimensional basis function); 
! ! see Section 2.2.1 in JMMV (2014)

if (i_max == 1) then 
    m_i_max = 1         ! If i_max=1, then m(i_max)=1, i.e., set  
                        ! S_1={0} and the maximum subindex is 1
else
    m_i_max =  2**(i_max-1) + 1
endif                   ! If i_max>1, then m(i_max)= 2^(i_max-1)+1;
                        ! e.g., for S_2={0,-1,1}, the maximum
                        ! subindex is 3

                                 
numb_pts = size(points,1)       ! Compute the number of points (rows),   
                                ! "numb_pts" in the matrix of points 
                                ! "points", in which the polynomial bases  
                                ! must be evaluated              
allocate(phi(1:numb_pts,1:d,1:m_i_max))
phi = 1
                                ! Allocate memory to a matrix of the
                                ! unidimensional polynomial bases "phi_n",  
                                ! evaluated in all the points 
                     
! For a polynomial bases "phi_n" with n=1, we have phi_n(x)=1 for all x; 
! our phi(:,:,1) is a matrix of ones of size numb_pts-by-d by the above 
! construction
                               
phi(:,:,2) = points;            ! For a polynomial bases "phi_n" with n=2, 
                                ! we have phi_n(x) is x; evaluating it in 
                                ! all the points gives us matrix "points"; 
                                ! numb_pts-by-d                    
do j = 3,m_i_max                ! For polynomial bases "phi_n", from n=3 to
                                ! n=m_i_max, ...
    phi(:,:,j) = 2*phi(:,:,2)*phi(:,:,j-1) - phi(:,:,j-2);
                                ! Use the recurrence formula to compute the 
                                ! Chebyshev polynomial basis functions of 
                                ! the degrees from 2 to m_i_max-1
enddo 

 
! 2. Form the multidimensional polynomial bases of Smolyak polynomial of the 
! required level of Smolyak approximation; see JMMV (2014), Sections 3.3.3 
! and 3.4.2 for examples
! ----------------------------------------------------------------------
! Smol_bases = [];             ! Initially, the matrix of multidimensional 
                               ! polynomial bases is empty
                               
numb_terms = size(smol_elem,1);
                              ! Compute the number of terms (i.e., multi-
                              ! dimensional polynomial bases) in Smolyak 
                              ! polynomial                                

allocate(pproduct(1:numb_pts))
 
do jt = 1,numb_terms          ! For each term of Smolyak polynomial, ...
                                
     index_row = smol_elem(jt,:)
                               ! Identify the subindices of the unidimensional
                               ! basis function that constitute a jt-th multi-
                               ! dimensional basis function; this is a jt-th
                               ! row of matrix "smol_elem"; 1-by-d
     pproduct = 1.0_dp
                               ! Initialize a vector of pproducts of unidimen-
                               ! sional basis functions; numb_pts-by-1
     do jd = 1,d              ! For each dimension (state variable), ...
         n = index_row(jd)    ! A subindex of a unidimensional basis function 
                               ! phi_n in a dimension jd is denoted n
         if (.not.(n == 1)) then             ! If the subindex of unidimensional basis 
                               ! function is not equal to unity, ...
         pproduct(1:numb_pts) = pproduct(1:numb_pts)*phi(:,jd,n);
                               ! Compute the pproduct of basis functions
                               
         ! Otherwise, i.e., if n = 1, there is no need to compute the 
         ! pproduct of unidimensional basis functions, as it's equal to unity 
        endif
     enddo


     if (allocated(Smol_bases)) then
     Smol_bases_temp = Smol_bases
     deallocate(smol_bases)
     allocate(Smol_bases(1:numb_pts,size(Smol_bases_temp,2)+1))
     Smol_bases(:,1:size(Smol_bases_temp,2)) = Smol_bases_temp
     Smol_bases(:,size(Smol_bases_temp,2)+1) = pproduct
     else 
     allocate(Smol_bases(1:numb_pts,1:1))
     Smol_bases(1:numb_pts,1) = pproduct
     endif
                               ! Attach to the previously obtained matrix of 
                               ! multidimensional basis functions a new
                               ! product of unidimensional basis functions;
                               ! e.g., for mu=1 and d=2, basis_bs is of
                               ! size numb_pts-by-5
enddo

end function Smolyak_Polynomial

function  Smolyak_Polynomial2(points,d,numb_pts,numb_terms,mu,smol_elem) result(smol_bases)
 
integer, intent(in)  :: d, numb_pts, numb_terms,  mu, smol_elem(numb_terms,d)
real(dp), intent(in) :: points(numb_pts,d)
real(dp) :: smol_bases(numb_pts,numb_terms)

integer  :: i_max,  jt, jd, n, m_i_max, j
real(dp) :: phi(1:numb_pts,1:d,1:2**mu + 1), pproduct(1:numb_pts)
integer :: index_row(d)

! numb_pts = size(points,1)       ! Compute the number of points (rows),   
                                ! "numb_pts" in the matrix of points 
                                ! "points", in which the polynomial bases  
                                ! must be evaluated              
! allocate(phi(1:numb_pts,1:d,1:m_i_max))
! Smolyak polynomial is given by the sum of multidimensional basis functions, 
! multiplied by the coefficients; see formula (15) in JMMV (2014). By 
! convention, the first basis function is given by 1 (unity). 

! Unidimensional basis functions are given by Chebyshev polynomial bases; 
! in JMMV (2014), a unidimensional Chebyshev polynomial basis function of   
! degree n-1 is denoted by "phi_n", i.e., has a subindex n and we follow  
! this notation here

! 1. Construct the unidimensional basis functions and evaluate them in
! all the points of matrix "points"
! -------------------------------------------------------------------------
i_max = mu+1                    ! The maximum subindex of unidimensional
                                ! set S_i whose points are used to
                                ! construct Smolyak grid; e.g., for mu=1, 
                                ! we consider up to S_i_max={0,-1,1} where
                                ! i_max=mu+1=1+1=2
                                
! ! Compute the number of elements in the i_max-th unidimensional set of 
! ! elements S_i_max; this coincides with a maximum subindex of elements 
! ! (unidimensional grid point or unidimensional basis function); 
! ! see Section 2.2.1 in JMMV (2014)

if (i_max == 1) then 
    m_i_max = 1         ! If i_max=1, then m(i_max)=1, i.e., set  
                        ! S_1={0} and the maximum subindex is 1
else
    m_i_max =  2**(i_max-1) + 1
endif                   ! If i_max>1, then m(i_max)= 2^(i_max-1)+1;
                        ! e.g., for S_2={0,-1,1}, the maximum
                        ! subindex is 3

                                 
phi = 1.0
                                ! Allocate memory to a matrix of the
                                ! unidimensional polynomial bases "phi_n",  
                                ! evaluated in all the points 
                     
! For a polynomial bases "phi_n" with n=1, we have phi_n(x)=1 for all x; 
! our phi(:,:,1) is a matrix of ones of size numb_pts-by-d by the above 
! construction
                               
phi(:,:,2) = points;            ! For a polynomial bases "phi_n" with n=2, 
                                ! we have phi_n(x) is x; evaluating it in 
                                ! all the points gives us matrix "points"; 
                                ! numb_pts-by-d                    
do j = 3,m_i_max                ! For polynomial bases "phi_n", from n=3 to
                                ! n=m_i_max, ...
    phi(:,:,j) = 2*phi(:,:,2)*phi(:,:,j-1) - phi(:,:,j-2);
                                ! Use the recurrence formula to compute the 
                                ! Chebyshev polynomial basis functions of 
                                ! the degrees from 2 to m_i_max-1
enddo 

 
! 2. Form the multidimensional polynomial bases of Smolyak polynomial of the 
! required level of Smolyak approximation; see JMMV (2014), Sections 3.3.3 
! and 3.4.2 for examples
! ----------------------------------------------------------------------
! Smol_bases = [];             ! Initially, the matrix of multidimensional 
                               ! polynomial bases is empty
                               
! numb_terms = size(smol_elem,1);
                              ! Compute the number of terms (i.e., multi-
                              ! dimensional polynomial bases) in Smolyak 
                              ! polynomial                                

! allocate(pproduct(1:numb_pts))
 
do jt = 1,numb_terms          ! For each term of Smolyak polynomial, ...
                                
     index_row = smol_elem(jt,:)
                               ! Identify the subindices of the unidimensional
                               ! basis function that constitute a jt-th multi-
                               ! dimensional basis function; this is a jt-th
                               ! row of matrix "smol_elem"; 1-by-d
     pproduct = 1.0_dp
                               ! Initialize a vector of pproducts of unidimen-
                               ! sional basis functions; numb_pts-by-1
     do jd = 1,d              ! For each dimension (state variable), ...
         n = index_row(jd)    ! A subindex of a unidimensional basis function 
                               ! phi_n in a dimension jd is denoted n
         if (.not.(n == 1)) then             ! If the subindex of unidimensional basis 
                               ! function is not equal to unity, ...
         pproduct(1:numb_pts) = pproduct(1:numb_pts)*phi(:,jd,n);
                               ! Compute the pproduct of basis functions
                               
         ! Otherwise, i.e., if n = 1, there is no need to compute the 
         ! pproduct of unidimensional basis functions, as it's equal to unity 
        endif
     enddo

     smol_bases(:,jt) = pproduct

   !  if (allocated(Smol_bases)) then
   !  Smol_bases_temp = Smol_bases
   !  deallocate(smol_bases)
   !  allocate(Smol_bases(1:numb_pts,size(Smol_bases_temp,2)+1))
   !  Smol_bases(:,1:size(Smol_bases_temp,2)) = Smol_bases_temp
   !  Smol_bases(:,size(Smol_bases_temp,2)+1) = pproduct
   !  else 
   !  allocate(Smol_bases(1:numb_pts,1:1))
   !  Smol_bases(1:numb_pts,1) = pproduct
   !  endif
                               ! Attach to the previously obtained matrix of 
                               ! multidimensional basis functions a new
                               ! product of unidimensional basis functions;
                               ! e.g., for mu=1 and d=2, basis_bs is of
                               ! size numb_pts-by-5
enddo

end function Smolyak_Polynomial2

end module mod_smolyak

! LICENSE AGREEMENT 
!
! Lilia Maliar, Serguei Maliar and Rafael Valero agree to make the Software for 
! Smolyak algorithm accompanying the article "Smolyak method for solving dynamic
! economic models: Lagrange interpolation, anisotropic grid and adaptive
! domain" by Kenneth L. Judd, Lilia Maliar, Serguei Maliar and Rafael Valero 
! (2014). Journal of Economic Dynamics and Control. Volume 44, July 2014, 
! Pages 92-123 available to you under the following conditions.
!
! 1. Any work that contains some results derived from the Software or a
!    modified version including translations to other languages must: 
! 
!    a. provide a prominent description of the use of the Software in the text; 
! 
!    b. cite the article "Smolyak method for solving dynamic economic models: 
!    Lagrange interpolation, anisotropic grid and adaptive domain" by Kenneth L. 
!    Judd, Lilia Maliar, Serguei Maliar and Rafael Valero (2014). Journal of 
!    Economic Dynamics and Control. Volume 44, July 2014, Pages 92-123. 
! 
! 2. The Software and any modified version may be redistributed under the 
!    terms of the present license agreement. Redistributions must include 
!    the original Software, this license agreement and the source code of 
!    the modified version with a prominent notice stating that the Software 
!    was changed and the date of the change.
! 
! 3. The Software comes as is and is used at your own risk. Any modification 
!    of the Software is your own responsibility.  The Authors assume no 
!    obligation to provide assistance with using the Software, make no 
!    warranties concerning the performance of the Software and are not 
!    liable to you or any entity for any damage related to the use of the 
!    Software. 
! 
! 4. The use of the Software is restricted to noncommercial research and 
!    educational purposes. If the Software or a modified version leads to the 
!    development of materials (including software, teaching materials, patents), 
!    which may be used for commercial purposes, a specific agreement must be 
!    negotiated between the authors and the beneficiary. 
! 
! 5. The Software is protected by copyright and other applicable laws. 
