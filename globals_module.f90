module globals
  !  This module contains the global variables and derived types used
  !  throughout the program
  use precision
  public
  !
  !  cut_off is the maximum size of the wave function at the last basis vector
  !
  
  real(fdp), parameter :: cut_off=1.0d-5
  !
  !  pi is the parameter value of pi
  !
  real(fdp), parameter :: pi=3.14159265359
  !
  !  hbar is the parameter value of the physical constant hbar
  !
  
  real(fdp), parameter :: hbar=1.01d-34
  !
  !  w is the angular frequency of the harmonic oscillator
  !
  
  real(fdp), parameter :: w=1.0d0
  !
  !  I is the complex variable i
  !
  
  complex(cdp), parameter :: I=(0,1)
  !
  !  Various derived types
  !
  
  type basis
     !  This derived type uses a displaced number state to represent the basis
     !  The alpha is the displacement, and n is the number of basis states
     !  being used.
     complex(cdp) :: alpha
     integer(isp) :: n
  end type basis
  type wave_fn
     !  This derived type defines a wave function as a set of co-efficients
     !  and the basis that the wave function is represented on.
     complex(cdp), dimension(:), pointer :: co_eff
     type (basis) :: in_basis
  end type wave_fn
  type gauss
     
     !  This derived type defines the data structure of the wiener
     !  random numbers
     real(fdp), dimension(:,:), pointer :: num
  end type gauss
end module globals
