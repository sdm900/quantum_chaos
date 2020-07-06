program test
  !  This test program sets up the variables to be integrated
  use integrate_sde
  use operators
  use globals
  use basis_change
  integer(idp) :: n=30
  type(wave_fn) :: y
  !
  !  This is where the initial wave function is set up.  This defines
  !  a coherent state as the initial condition.  In the displaced basis
  !  states, a coherent state is simply a 1 in the |alpha,0> state
  !
  y%in_basis%n=n
  y%in_basis%alpha=cmplx(10.0d0,10.0d0)
  allocate(y%co_eff(0:n))
  y%co_eff=cmplx(0.0d0,0.0d0)
  y%co_eff(0)=cmplx(1.0d0)
  y%co_eff = y%co_eff / normalisation(y)
  open(unit=1,file="wave.dat")
  write(1,'(e15.7,e15.7)') y%co_eff
  close(1)
  
  !
  !  Perform the integration
  !
  y=integrate(y,0.0d0,100.0d0,10000000)
  deallocate(y%co_eff)
end program test
