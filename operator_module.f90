module operators
  !  This module contains functions that evaluates the effect of operators on
  !  the wave functions
  use precision
  use globals
  interface operator (.a.)
     !  This operator (a) is the annihilation operator
     !
     !  expectation_a is the function that caclulates <y|a|y>
     !  oper_a is the function that calculates a|y>
     module procedure expectation_a
     module procedure oper_a
  end interface
  interface operator (.ad.)
     !  This operator (ad) is the creation operator
     !
     !  expectation_ad is the function that caclulates <y|ad|y>
     !  oper_ad is the function that calculates ad|y>
     module procedure expectation_ad
     module procedure oper_ad
  end interface
  interface operator (+)
     !  This overloads the operator + to apply to wave functions.
     !
     !  add_wave_fn adds co-efficients of the wave functions together
     module procedure add_wave_fn
  end interface
  interface operator (*)
     !  This overload the operator * to apply to wave functions.
     !
     !  scalar_mult_wave_fn multiplies the co-efficients of a wave function
     !    by a scalar
     !  wave_fn_mult_scalar multiplies the co-efficients of a wave function
     !    by a scalar
     module procedure scalar_mult_wave_fn
     module procedure wave_fn_mult_scalar
  end interface
  
  
contains
  
  function add_wave_fn(y1,y2) result(tmpy)
    !  This function add the co-efficients of wave functions together
    !
    !  y1 and y2 are the input wave functions to add together
    !  tmpy is the result of the operation
    type(wave_fn), intent(in) :: y1, y2
    type(wave_fn) :: tmpy
    if (y1%in_basis%alpha /= y2%in_basis%alpha .or. y1%in_basis%n /= y2%in_basis%n) then
       
       write(*,*) 'Can not add two wave functions that are dissimilar'
    else
       
       tmpy%in_basis=y1%in_basis
       allocate(tmpy%co_eff(tmpy%in_basis%n))
       tmpy%co_eff=y1%co_eff + y2%co_eff
    end if
  end function add_wave_fn
  function scalar_mult_wave_fn(a,y) result(tmpy)
    !  This multiplies the co-efficients of a wave function by a scalar
    !
    !  a is the scalar to multiply the wave function by
    !  y is the wave function to be multiplied by the scalar
    !  tmpy is the result of the operation
    type(wave_fn), intent(in) :: y
    complex(cdp), intent(in) :: a
    type(wave_fn) :: tmpy
    tmpy%in_basis=y%in_basis
    allocate(tmpy%co_eff(0:tmpy%in_basis%n))
    tmpy%co_eff=a*y%co_eff
  end function scalar_mult_wave_fn
  function wave_fn_mult_scalar(y,a) result(tmpy)
    
    !  This multiplies the co-efficients of a wave function by a scalar
    !
    !  a is the scalar to multiply the wave function by
    !  y is the wave function to be multiplied by the scalar
    !  tmpy is the result of the operation
    
    type(wave_fn), intent(in) :: y
    complex(cdp), intent(in) :: a
    type(wave_fn) :: tmpy
    
    
    tmpy%in_basis=y%in_basis
    
    allocate(tmpy%co_eff(0:tmpy%in_basis%n))
    
    tmpy%co_eff=a*y%co_eff
  end function wave_fn_mult_scalar
  
  function oper_a(y) result(tmpy)
    !  This func calculates a|y> and returns the appropriate value
    !
    !  y is the wave function is the operation a|y>
    !  tmpy is the result of the operation
    !  loop is used to loop over all the co-efficients in the wave function
    !  alpha is the point that the basis is centered at
    type(wave_fn), intent(in) :: y
    type(wave_fn) :: tmpy
    integer(isp) :: loop
    complex(cdp) :: alpha
    tmpy%in_basis = y%in_basis
    alpha = y%in_basis%alpha
    
    allocate(tmpy%co_eff(0:tmpy%in_basis%n))
    tmpy%co_eff=cmplx(0)
    do loop=0,tmpy%in_basis%n-1
     
       tmpy%co_eff(loop) = sqrt(dble(loop+1))*y%co_eff(loop+1) + alpha*y%co_eff(loop)
    end do
    tmpy%co_eff(tmpy%in_basis%n) = alpha*y%co_eff(y%in_basis%n)
  end function oper_a
  function expectation_a(y1,y2) result(alpha)
    
    !  This func calculates the expectation value of <y|a|y>
    !
    !  y1 and y2 are the wave functions for the operation <y1|a|y2>
    !  alpha is the point that the basis of y2 is centered upon
    !  loop is used to loop over the co-efficients of the wave functions
    !  minn is the minimum number of basis states of the two wave functions
    
    type(wave_fn), intent(in) :: y1, y2
    complex(cdp) :: alpha
    integer(isp) :: loop, minn
    minn=min(y1%in_basis%n,y2%in_basis%n)
    
    alpha=y2%in_basis%alpha
    
    do loop=0,minn-1
       
       alpha=alpha + sqrt(dble(loop+1))*conjg(y1%co_eff(loop))*y2%co_eff(loop+1)
       
    end do
    
  end function expectation_a
  
  function oper_ad(y) result(tmpy)
    !  This func calculates ad|y> and returns the appropriate value
    !
    !  y is the wave function is the operation ad|y>
    !  tmpy is the result of the operation
    !  loop is used to loop over all the co-efficients in the wave function
    !  alpha is the point that the basis is centered at
    type(wave_fn), intent(in) :: y
    type(wave_fn) :: tmpy
    integer(isp) :: loop
    complex(cdp) :: alpha
    tmpy%in_basis = y%in_basis
    alpha=y%in_basis%alpha
    
    allocate(tmpy%co_eff(0:tmpy%in_basis%n))
    tmpy%co_eff=cmplx(0)
    do loop=1,tmpy%in_basis%n
     
       tmpy%co_eff(loop) = sqrt(dble(loop))*y%co_eff(loop-1) + conjg(alpha)*y%co_eff(loop)
    end do
    tmpy%co_eff(0) = conjg(alpha)*y%co_eff(0)
  end function oper_ad
  function expectation_ad(y1,y2) result(alpha)
    
    !  This func calculates the expectation value of <y|ad|y>
    !
    !  y1 and y2 are the wave functions for the operation <y1|ad|y2>
    !  alpha is the conjugate of the point the basis of y2 is centered upon
    !  loop is used to loop over the co-efficients of the wave functions
    !  minn is the minimum number of basis states of the two wave functions
    
    type(wave_fn), intent(in) :: y1, y2
    complex(cdp) :: alpha
    integer(isp) :: loop,minn
    minn=min(y1%in_basis%n,y2%in_basis%n)
    
    alpha=conjg(y2%in_basis%alpha)
    
    do loop=1,minn
       
       alpha=alpha + sqrt(dble(loop))*conjg(y1%co_eff(loop))*y2%co_eff(loop-1)
       
    end do
    
  end function expectation_ad
  function normalisation(y) result(norm)
    !  This function will return a constant to normalise a wave function
    !
    !  y is the wave function that is needed to be normalised
    !  norm is the constant that will normalise the wave function y
    type(wave_fn), intent(in) :: y
    real(fdp) :: norm
    norm = real(sqrt(sum(y%co_eff*conjg(y%co_eff))))
  end function normalisation
    
end module operators
