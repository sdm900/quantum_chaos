program test
  use integrate_sde
  use operators
  use globals
  use basis_change
  integer(idp) :: n=200
  type(wave_fn) :: y, y1, y2
  complex(cdp) :: alpha
  y%in_basis%n=n
  y%in_basis%alpha=cmplx(0)
  allocate(y%co_eff(0:n), y2%co_eff(0:n))
  y%co_eff=cmplx(1.0d0,0.0d0)
  write(*,*) 'Enter alpha '
  read(*,*) alpha
  do j=1,y%in_basis%n
     
     y%co_eff(j)=alpha/sqrt(dble(j))*y%co_eff(j-1)
  end do
  y%co_eff = y%co_eff / normalisation(y)
  y1=.ad.y
  y2%co_eff=(y.ad.y)*y%co_eff
  open(unit=1, file="y1.dat")
  write(1,'(e20.8, e20.8)') y1%co_eff
  close(1)
  open(unit=1, file="y2.dat")
  write(1,'(e20.8, e20.8)') y2%co_eff
  close(1)
  write(*,*) y.ad.y,(sum(abs(y1%co_eff-y2%co_eff)))/sum(abs(y2%co_eff))
  deallocate(y%co_eff)
end program test
