program test
  use rand_num_gen
  use precision
  integer(idp) :: loop1,loop2
  real(fdp) :: ave, var, inv
  complex(cdp) :: tmp
  real(fdp), parameter :: pi=3.14159265359
  open(unit=1,file="data/var.good.dat")
  open(unit=2,file="data/ave.good.dat")
  do loop1=0,8
     var=dble(0)
     ave=dble(0)
     do loop2=1,10**loop1
        tmp=sqrt(-2*log(random_num()))*exp(cmplx(0,2*pi*random_num()))
        
        inv=1/dble(loop2)
        ave=ave - (ave - real(tmp))*inv
        var=var - (var - real(tmp)**2)*inv
        
     end do
     
     write(1,'(e15.7)') var
     write(2,'(e15.7)') ave
     write(*,*) 'var ',var
     write(*,*) 'ave ',ave
     write(*,*) ' '
  end do
  close(1)
  close(2)
end program test
