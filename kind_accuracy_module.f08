!*****************************************************************************************
!> author: Kyriacos Skoufaris
!
!  Define the numeric kinds.
module KIND_ACCURACY

  implicit none

  !!! for integer values that are larger than selected_int_kind(6) (+-2147483647) you must declare it explicitly ex. 10000000000_iac
  !!! for real values that are more accurate than selected_int_kind(6) (8 desimals) you must declare it explicitly ex. 1._rac = 1.000000000000000

  integer, parameter :: iac = selected_int_kind(15)
  integer, parameter :: rac = selected_real_kind(15, 307)
  !integer, parameter :: rac = selected_real_kind(33, 4931)

  !!! SELECTED_INT_KIND(R) return the kind value of the smallest integer type that can represent all values ranging from -10^R (exclusive) to 10^R (exclusive)
  !integer, parameter :: k8 = selected_int_kind(8)

  !!! SELECTED_REAL_KIND(P,R) returns the kind value of a real data type with decimal precision of at least P digits, exponent range of at least R
  !!! precision of 32-, 64-, and 128-bit reals
  !integer, parameter :: sp = selected_real_kind(6, 37)
  !integer, parameter :: dp = selected_real_kind(15, 307)
  !integer, parameter :: qp = selected_real_kind(33, 4931)

  ! precision up to that of the machine-compiler-specifics,
  ! ensures that the double and quad types are actually twice and four times the precision of a single
  !!!integer, parameter ::                            &
  !!!  sp = kind(1.0),                                &
  !!!  dp = selected_real_kind(2*precision(1.0_sp)),  &
  !!!  qp = selected_real_kind(2*precision(1.0_dp))

end module KIND_ACCURACY
!*****************************************************************************************
