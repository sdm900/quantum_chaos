program test
  !  This program test the change of basis code.
  
  use precision
  use globals
  use basis_change
  use operators
  
  implicit none
  type(wave_fn) :: y1, y2, y1back
  integer(idp) ::  loop
  complex(cdp) :: alpha
  
  alpha=cmplx(2,2)
  y1%in_basis%n=500
  y2%in_basis%n=500
  y1back%in_basis%n=500
  allocate(y1%co_eff(0:y1%in_basis%n), y2%co_eff(0:y2%in_basis%n), y1back%co_eff(0:y1%in_basis%n))
  
  y1%in_basis%alpha=cmplx(0)
  y1%co_eff=cmplx(1.0d0)
  y1back%in_basis=y1%in_basis
  y1back%co_eff=cmplx(0)
  y2%co_eff=cmplx(0.0d0)
  
  do loop=1,y1%in_basis%n
     
     y1%co_eff(loop)=alpha/sqrt(dble(loop))*y1%co_eff(loop-1)
     
  end do
  y1%co_eff=y1%co_eff / normalisation(y1)
  y2%in_basis%alpha=y1.a.y1
  write(*,*) 'expectation y1 ',y1.a.y1
  call change_basis(y1,y2)
  write(*,*) 'expectation y2 ',y2.a.y2
  call change_basis(y2,y1back)
  write(*,*) 'expectation y1back ',y1back.a.y1back
  write(*,*) 'error in change ',sum(abs(y1%co_eff-y1back%co_eff))/sum(abs(y1%co_eff))
  open(unit=1, file="y1.dat")
  write(1,'(e15.7,e15.7)') y1%co_eff
  close(1)
  open(unit=1, file="y2.dat")
  write(1,'(e15.7,e15.7)') y2%co_eff
  close(1)
  open(unit=1, file="y1back.dat")
  write(1,'(e15.7,e15.7)') y1back%co_eff
  close(1)
  deallocate(y1%co_eff, y2%co_eff, y1back%co_eff)
  
end program test
