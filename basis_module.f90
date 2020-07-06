module basis_change
  !  This module will contain all the procedures to do a change of basis.
  !
  !  globals is another modules that has the global variables and derived
  !  types defined.
  !
  !  There are three subroutines defined indes this module.  "compute_factors"
  !   computes some factors used in determining the transformation matrix,
  !  "compute_delta" actutally computes the general transformation matrix,
  !   and "change_basis" changes the basis of a wave function.
  !
  !  The global variable delta is the transformation matrix
  use precision
  use globals
  use operators
  
  implicit none
  complex(cdp), dimension(:,:), allocatable, private :: delta
contains
  subroutine compute_delta(b1,b2)
    
    !  This computes the transformation matrix delta between two basis
    !
    !  b1 is the basis to change from
    !  b2 is the basis to change to
    !  facd and faccd are two factors that are summed over to get delta
    !  i is the variable used to loop over basis b2
    !  j is the variable used to loop over basis b1
    !  l is a loop variable
    !  minij contains the minimum of i,j at various stages in the program
    !  d is the difference between the two basis
    !  cd is the conjugate of d
    !  fac a tempory variable
    !  dcd is d*cd
    
    type(basis), intent(in) :: b1, b2
    complex(cdp), dimension(:,:), allocatable :: facd, faccd
    integer(idp) :: i,j,l,minij
    complex(cdp) :: cd, d
    real(fdp) :: fac, dcd
    if (allocated(delta)) then
       deallocate(delta)
    end if
    allocate(facd(0:b2%n,0:b2%n), faccd(0:b1%n,0:b1%n), delta(0:b2%n,0:b1%n))
    delta=cmplx(0)
    facd=cmplx(0)
    faccd=cmplx(0)
    d=b1%alpha - b2%alpha
    cd=conjg(d)
    dcd=real(d*cd)
    !
    !  Sets up facd and faccd for iterative computation
    !
    
    do i=0,b2%n
       facd(i,i)=cmplx(1)
    end do
    do j=0,b1%n
       faccd(j,j)=cmplx(1)
    end do
    !
    !  Computes (iteratively) facd and faccd
    !
    do i=1,b2%n
       do l=i-1,0,-1
          facd(i,l)=sqrt(dble(l+1))/dble(i-l)*d*facd(i,l+1)
       end do
    end do
    do j=1,b1%n
       do l=j-1,0,-1
          faccd(j,l)=-sqrt(dble(l+1))/dble(j-l)*cd*faccd(j,l+1)
          
       end do
    end do
    
    !
    !  Computes delta from facd and faccd
    !
    
    do i = 0,b2%n
       do j = 0,b1%n
          minij=min(i,j)
          do l = 0,minij-1,2
             fac=1-(i-l)*(j-l)/(dcd*(l+1))
             delta(i,j)=delta(i,j) + facd(i,l)*faccd(j,l)*fac
          end do
          if ((minij/2)*2 == minij) then
             delta(i,j)=delta(i,j) + facd(i,minij)*faccd(j,minij)
          end if
       end do
    end do
    deallocate(facd, faccd)
    
  end subroutine compute_delta
  subroutine change_basis(y1,y2)
    !  This subroutine changes the wave function y1 to y2, depending on 
    !  the basis these wave functions are defined in.
    type(wave_fn), intent(in) :: y1
    type(wave_fn), intent(inout) :: y2
    call compute_delta(y1%in_basis,y2%in_basis)
    
    y2%co_eff = matmul(delta,y1%co_eff)
    
    y2%co_eff = y2%co_eff / normalisation(y2)
    
  end subroutine change_basis
end module basis_change
