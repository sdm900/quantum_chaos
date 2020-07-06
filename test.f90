program test
  
  use precision
  use globals
  use operators
  type(wave_fn) :: y, y2
  y%in_basis%alpha=cmplx(0)
  y%in_basis%n=20
  allocate(y%co_eff(0:y%in_basis%n))
  y%co_eff=cmplx(0)
  y%co_eff(2)=cmplx(1,1)
  y%co_eff=y%co_eff/normalisation(y)
  write(*,'(f15.7,f15.7)') y%co_eff
  write(*,*) '=========='
  y2=.ad.y
  write(*,'(f15.7,f15.7)') y2%co_eff
end program test
