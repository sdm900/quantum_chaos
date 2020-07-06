module precision
  !  This module will assign fdp (floating point double precision)
  !  isp (integer single precision), idp (integer double precision)
  !  csp (complex single precision), cdp (complex double precision)
  !  fqp (floating point quadratic precision)
  implicit none
  integer, parameter :: fqp=selected_real_kind(20,200)
  integer, parameter :: fdp=selected_real_kind(10,200)
  integer, parameter :: fsp=selected_real_kind(5,50)
  integer, parameter :: isp=selected_int_kind(6)
  integer, parameter :: idp=selected_int_kind(9)
  integer, parameter :: cdp=fdp
  integer, parameter :: csp=fsp
end module precision
