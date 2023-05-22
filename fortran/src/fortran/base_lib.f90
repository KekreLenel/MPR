! -------------------------------------------------------------------------
! base_lib.f90: this module contains various auxiliary functions, 
! most of which were copied from other sources 
! -------------------------------------------------------------------------
! authors:         Rohan Kekre and Moritz Lenel
! for updates see: https://github.com/KekreLenel/MPR
! -------------------------------------------------------------------------
module base_lib
use omp_lib
use nag_library, only: nag_wp
implicit none
private
public :: Fill_linspace_dp, m_choose_r, pi, dp, eps, sqrt_eps, choldc, QsortC, calc_var, calc_cov, & 
          remove_dups_2d, remove_dups_1d_real, get_quadrature_points

integer, parameter   :: dp       = nag_wp

real(dp), parameter  :: pi       = 4*atan(1.0_dp), &
                     &  eps      = epsilon(eps),   &
                     &  sqrt_eps = sqrt(eps)

contains


    ! Presumable copied from the excellent library of John Burkardt
    ! https://people.sc.fsu.edu/~jburkardt/f_src/
    subroutine get_quadrature_points(mu, sigma, n, points, weights)
    
    real(kind=dp), intent(IN) :: mu, sigma
    integer, intent(IN) :: n
    real(kind=dp), intent(OUT) :: points(n), weights(n)
    integer :: i, info
    real(kind=dp) ::   workspace(4*n)
    real(kind=dp) :: diag(n), off_diag(n-1), J_mat(n,n)
    
        diag(:) = 0.0_dp;
    
    
        do i = 1,n-1
            off_diag(i) = sqrt(1.0_dp*i)
        end do
        call F08JEF('I', n, diag, off_diag,  J_mat, n,  workspace, INFO)
    
        points = diag
        points = mu+sigma*points
        weights = J_mat(1,1:n)**2
    
    
    end subroutine get_quadrature_points

   ! Copied from https://github.com/astrofrog/fortranlib/blob/master/src/lib_array.f90
   subroutine Fill_linspace_dp(xmin,xmax,x)
      implicit none
      real(dp),intent(in) :: xmin,xmax
      real(dp),intent(inout) :: x(:)
      integer :: i,n
      n = size(x)
      if (n == 1) then
         if(xmin /= xmax) then
            write(0,'("ERROR: Cannot call linspace with n=1 and xmin /= xmax")')
            stop
         else
            x = xmin
         end if
      else
         x(1) = xmin
         do i=2,n-1
             x(i) = (xmax-xmin) * real(i-1,dp) / real(n-1,dp) + xmin
         end do
         x(n) = xmax
      end if
   end subroutine Fill_linspace_dp


    RECURSIVE FUNCTION m_choose_r(m,r_temp) RESULT(res)
    INTEGER :: m, r_temp, r, res
    IF (m < r_temp) THEN
         write(*,*) 'Warning M_CHOOSE_R'
         stop
    END IF

    r = min(m-r_temp, r_temp)
    IF (m == r) THEN
         res = 1
    ELSE IF (r == 0) THEN
         res = 1
    ELSE IF (r == 1) THEN
         res = m
    ELSE
         res = m*m_choose_r_second(m-1,r-1)/r
    END IF
   END FUNCTION m_choose_r

    RECURSIVE FUNCTION m_choose_r_second(m,r) RESULT(res)
    INTEGER :: m, r, res

    IF (r == 1) THEN
         res = m
    ELSE
         res = m*m_choose_r_second(m-1,r-1)/r
    END IF
    END FUNCTION m_choose_r_second


   SUBROUTINE choldc(a,n,p,ifail)
   INTEGER, intent(IN) :: n
   REAL(kind=dp), intent(INOUT) :: a(n,n)
   REAL(kind=dp), intent(OUT) :: p(n)
   INTEGER, intent(INOUT) :: ifail
   INTEGER :: i,j,k
   REAL(kind=dp) :: sum
   do i=1,n
      do j=i,n
         sum=a(i,j)
         do k=i-1,1,-1
            sum=sum-a(i,k)*a(j,k)
         enddo
         if(i.eq.j)then
            if(sum.le.0.) then
               ifail = 1
               sum = sqrt(sum**2)
            end if
            p(i)=sqrt(sum)
         else
            a(j,i)=sum/p(i)
         endif
      enddo
   enddo
   return
   END

recursive subroutine QsortC(A)
  real(dp), intent(in out), dimension(:) :: A
  integer :: iq

  if(size(A) > 1) then
     call Partition(A, iq)
     call QsortC(A(:iq-1))
     call QsortC(A(iq:))
  endif
end subroutine QsortC

subroutine Partition(A, marker)
  real(dp), intent(in out), dimension(:) :: A
  integer, intent(out) :: marker
  integer :: i, j
  real(dp) :: temp
  real(dp) :: x      ! pivot point
  x = A(1)
  i= 0
  j= size(A) + 1

  do
     j = j-1
     do
        if (A(j) <= x) exit
        j = j-1
     end do
     i = i+1
     do
        if (A(i) >= x) exit
        i = i+1
     end do
     if (i < j) then
        ! exchange A(i) and A(j)
        temp = A(i)
        A(i) = A(j)
        A(j) = temp
     elseif (i == j) then
        marker = i+1
        return
     else
        marker = i
        return
     endif
  end do

end subroutine Partition


function calc_mean(vec, n_vec) result(mean)
    integer, intent(in)   :: n_vec
    real(dp), intent(in)  :: vec(n_vec)
    real(dp) :: mean

    mean = sum(vec)/real(n_vec,dp)

end function calc_mean

function calc_var(vec, n_vec) result(var)
    integer, intent(in)   :: n_vec
    real(dp), intent(in)  :: vec(n_vec)
    real(dp) :: var

    real(dp) :: mean

    mean = calc_mean(vec, n_vec)
    var  = sum((vec-mean)**2)/real(n_vec-1,dp)

end function calc_var

function calc_cov(vec1, vec2, n_vec) result(cov)
    integer, intent(in)   :: n_vec
    real(dp), intent(in)  :: vec1(n_vec), vec2(n_vec)
    real(dp) :: cov

    real(dp) :: mean1, mean2

    mean1 = calc_mean(vec1, n_vec)
    mean2 = calc_mean(vec2, n_vec)
    cov  = sum((vec1-mean1)*(vec2-mean2))/real(n_vec-1,dp)

end function calc_cov

! adapted from https://rosettacode.org for 2d matrices along rows
function remove_dups_2d(input_mat) result(output_mat)

    implicit none
    integer, intent(in)  :: input_mat(:,:)  ! The input
    integer, allocatable :: output_mat(:,:) ! The output
    integer :: k                            ! The number of unique elements
    integer :: i, j
 
    k = 1
    output_mat = input_mat
    output_mat = -9999
    output_mat(1,:) = input_mat(1,:)
    outer: do i=2,size(input_mat,1)
        do j=1,k
            if (all(output_mat(j,:).eq.input_mat(i,:))) then
                ! Found a match so start looking again
                cycle outer
            end if
        end do
        ! No match found so add it to the output
        k = k + 1
        output_mat(k,:) = input_mat(i,:)
    end do outer
  
  output_mat = output_mat(1:k,:)

end function remove_dups_2d


! adapted from https://rosettacode.org for 2d matrices along rows
function remove_dups_1d_real(input_mat) result(output_mat)

    implicit none
    real(dp), intent(in)  :: input_mat(:)  ! The input
    real(dp), allocatable :: output_mat(:) ! The output
    integer :: k                          ! The number of unique elements
    integer :: i, j
 
    k = 1
    output_mat = input_mat
    output_mat = -9999
    output_mat(1) = input_mat(1)
    outer: do i=2,size(input_mat)
        do j=1,k
            if ( output_mat(j) == input_mat(i) ) then
                ! Found a match so start looking again
                cycle outer
            end if
        end do
        ! No match found so add it to the output
        k = k + 1
        output_mat(k) = input_mat(i)
    end do outer
  
    output_mat = output_mat(1:k)

end function remove_dups_1d_real

end module base_lib
