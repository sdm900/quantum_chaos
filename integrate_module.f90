module integrate_sde
  !  This module will provide the routines necessary to 
  !  integrate a stochastic differential equation.
  use precision
  use globals
  use rand_num_gen
  use basis_change
  use operators
  implicit none
  
  
contains
  function gauss_rand_num(num_vars,nstep) result(dW)
    
    !  This subroutine specifies the gaussian random numbers for use
    !  in the stochastic integration process.  The random numbers have
    !  variance 1, and average 0
    !
    !  This routine was motivated by the routine given in Numerical Recipes
    !
    !  num_vars is the number of variables in the stochastic integration
    !  nstep is the number of integration steps requiring random numbers
    !  dW are the random numbers in the required form
    !  loop1 and loop2 loop over num_vars and nstep respectively
    !  tmp_ran1 and tmp_ran2 are two temporary random numbers
    
    integer(idp), intent(in) :: nstep
    integer(isp), intent(in) :: num_vars
    type(gauss) :: dW
    integer(idp) :: loop1, loop2
    real(fdp) :: tmp_ran1, tmp_ran2
    allocate(dW%num(0:num_vars,nstep))       
    
    dW%num=dble(0)
    do loop1=0,num_vars
       do loop2=1,nstep-1,2
          tmp_ran1=random_num()
          tmp_ran2=random_num()
          
          dW%num(loop1,loop2)=sqrt(-2*log(tmp_ran1))*cos(2*pi*tmp_ran2)
          dW%num(loop1,loop2+1)=sqrt(-2*log(tmp_ran1))*sin(2*pi*tmp_ran2)
       end do
       if ((nstep/2)*2 /= nstep) then
          dW%num(loop1,nstep)=sqrt(-2*log(random_num()))*cos(2*pi*random_num())
       end if
    end do
  end function gauss_rand_num
  
  
  
  function integrate(initial, ti, tf, nstep) result(y)
    
    !  This function carries out the integration of the function f,
    !  given its derivatives as defined in the function func_derivs,
    !  from ti to tf, with nstep number of integration steps.
    !
    !  initial is the initial condition for integration
    !  ti and tf define the integration period
    !  nstep is the number of steps in the integration
    !  y is the result of the integration
    !  tmpy is a temporary wave function used for changing basis
    !  det_component is the deterministic component of integration
    !  stoc_component is the stochastic component of integration
    !  dW are the random numbers for integration of stochastic component
    !  tmp_det_comp is the co-efficients of the deterministic component
    !  tmp_stoc_comp is the co-efficients of the stochastic component
    !  loop is the current integration step
    !  dt is the time step
    !  ave is a variable of the average of the tail of the wave function
    !  ea is the expectation value <y|a|y>
    
    integer(idp), intent(in) :: nstep
    real(fdp), intent(in) :: ti, tf
    type(wave_fn), intent(in) :: initial
    type(wave_fn) :: y, tmpy, det_component, stoc_component
    type(gauss) :: dW
    complex(cdp), dimension(:), allocatable :: tmp_det_comp, tmp_stoc_comp
    integer(idp) :: loop
    real(fdp) :: dt, ave
    complex(cdp) :: ea
    
    allocate(y%co_eff(0:initial%in_basis%n))
    y%co_eff=initial%co_eff
    y%in_basis=initial%in_basis
    
    dt=(tf-ti)/dble(nstep)
    loop=0
    ea=y%in_basis%alpha
    open(unit=1,file="y.dat")
    
    !
    !  The integration process
    !
    do while (loop < nstep)
    
       !
       !  Integration loop while the wave function is localised
       !
    
       do while ((loop < nstep) .and. (abs(y%in_basis%alpha-ea) < 1) .and. (ave < 1.0d-6))
          dW=gauss_rand_num(y%in_basis%n,1)
          !
          !  Compute the derivatives
          !
          det_component = det_deriv(y)
          stoc_component = stoc_deriv(y)
          allocate(tmp_det_comp(0:y%in_basis%n), tmp_stoc_comp(0:y%in_basis%n))
          !
          !  Integrate
          !
          tmp_det_comp = -I*dt* det_component%co_eff
          
          tmp_stoc_comp = stoc_component%co_eff *dW%num(:,1)*sqrt(dt)
          y%co_eff=y%co_eff + tmp_det_comp
          y%co_eff = y%co_eff / normalisation(y)
          loop=loop+1
          ea=y.a.y
          if ((loop/100)*100==loop) then
             write(1,'(e15.7,e15.7)') ea
          
          end if
          deallocate(dW%num, tmp_det_comp, tmp_stoc_comp, det_component%co_eff, stoc_component%co_eff)
          !
          !  Compute the average of the last five co-efficients
          !
          ave=sum(abs(y%co_eff(y%in_basis%n-5:y%in_basis%n)))/5
          
       end do
       
       !
       !  Move basis to be centered about the wave function
       !
       tmpy%in_basis%n=y%in_basis%n
       tmpy%in_basis%alpha=ea
       allocate(tmpy%co_eff(0:tmpy%in_basis%n))
       tmpy%co_eff=cmplx(0)
       
       write(*,*) 'Before change ',y.a.y
       call change_basis(y,tmpy)
       write(*,*) 'Changed Basis at step ', loop
       write(*,*) 'After change ',tmpy.a.tmpy
       
       deallocate(y%co_eff)
       y%in_basis=tmpy%in_basis
       
       allocate(y%co_eff(0:y%in_basis%n))
       y%co_eff=tmpy%co_eff
       deallocate(tmpy%co_eff)
    end do
    close(1)
  end function integrate
  
     
  
  function det_deriv(y) result(tmpy)
    
    !  This function returns the deterministic derivatives used in the 
    !  integration process
    !
    !  y is the wave function to take the derivative of
    !  tmpy is the derivative of y
    !  loop is  a dummy variable to loop over
    !  n is the size of the basis
    !  alpha is the centre of the basis
    !  calpha is the conjugate of alpha
    !  aalpha is alpha*calpha
    
    type(wave_fn), intent(in) :: y
    type(wave_fn) :: tmpy
    integer(idp) :: loop, n
    complex(cdp) :: alpha, calpha
    real(fdp) :: aalpha
    tmpy%in_basis = y%in_basis
    alpha=y%in_basis%alpha
    calpha=conjg(alpha)
    aalpha=real(calpha*alpha)
    n=tmpy%in_basis%n
    allocate(tmpy%co_eff(0:tmpy%in_basis%n))
    !
    !  Computes the derivative according to the Hamiltonian
    !
    do loop=1,n-1
       tmpy%co_eff(loop)=(loop+0.5d0+aalpha)*y%co_eff(loop) + calpha*sqrt(dble(loop+1))*y%co_eff(loop+1) + alpha*sqrt(dble(loop))*y%co_eff(loop-1)
    end do
    tmpy%co_eff(0)=(0.5d0+aalpha)*y%co_eff(0) + calpha*y%co_eff(1)
    tmpy%co_eff(n)=(0.5d0+aalpha)*y%co_eff(n) + alpha*sqrt(dble(n))*y%co_eff(n-1) 
  end function det_deriv
  
  
  
  function stoc_deriv(y) result(tmpy)
    
    !  This function returns the stochastic component used in the 
    !  integration process
    !
    !  y is the wave function to take the derivative of
    !  tmpy is the derivative of y
    
    type(wave_fn), intent(in) :: y
    type(wave_fn) :: tmpy
    
    allocate(tmpy%co_eff(0:y%in_basis%n))
    tmpy%in_basis = y%in_basis
    tmpy%co_eff = cmplx(0)
  end function stoc_deriv
  
  
  
end module integrate_sde
