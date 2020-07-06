Module rand_num_gen
  !  This module will contain a routine to generate random numbers
  !
  !  This is taken from numerical recipes for F90
  !
  !  It has been modified slightly to make it easily used.  The basic workings
  !  of it have not been changed in any way.
  !
  !  The variable names don't seem to have a lot of meaning (I simply copied
  !  them from Numerical Recipes)
  use precision
  
  implicit none
  integer(idp), save, private :: idum=-1, ix=-1, iy=-1, c
  real(fdp), save, private :: am
  
contains
  function random_num(init) result(number)
    !  This function returns a rantom number.  The optional input init, tells
    !  the random number generator to re-initialise, if the user desires.
    integer(idp), optional, intent(in) :: init
    real(fdp) :: number
    integer(idp), parameter :: ia=16807, im=2147483647, iq=127773, ir=2836
    if (present(init)) then
       
       idum=init
    end if
    if (idum<=0 .or. iy < 0) then
       
       am=nearest(1.0,-1.0)/im
       iy=ior(ieor(888889999,abs(idum)),1)
       ix=ieor(777755555,abs(idum))
       idum=abs(idum)+1
    end if
    ix=ieor(ix,ishft(ix,13))
    ix=ieor(ix,ishft(ix,-17))
    ix=ieor(ix,ishft(ix,5))
    c=iy/iq
    iy=ia*(iy-c*iq)-ir*c
    if (iy < 0) iy=iy+im
    number=am*ior(iand(im,ieor(ix,iy)),1)
  end function random_num
end module rand_num_gen
