!*****************************************************************************************
!> author: Kyriacos Skoufaris
!
!  Defines some numeric parameters.

    module NUMBERS_CONSTANTS

    use KIND_ACCURACY

    private

    real (kind=rac), parameter, public :: zero       = 0.0_rac
    real (kind=rac), parameter, public :: one        = 1.0_rac
    real (kind=rac), parameter, public :: two        = 2.0_rac
    real (kind=rac), parameter, public :: three      = 3.0_rac
    real (kind=rac), parameter, public :: four       = 4.0_rac
    real (kind=rac), parameter, public :: five       = 5.0_rac
    real (kind=rac), parameter, public :: six        = 6.0_rac
    real (kind=rac), parameter, public :: seven      = 7.0_rac
    real (kind=rac), parameter, public :: eight      = 8.0_rac
    real (kind=rac), parameter, public :: nine       = 9.0_rac
    real (kind=rac), parameter, public :: ten        = 10.0_rac

    real (kind=rac), parameter, public :: pi         = atan2(one,zero)*two ! acos(-one)*two
    real (kind=rac), parameter, public :: twopi      = two*pi
    real (kind=rac), parameter, public :: fourpi     = four*pi
    real (kind=rac), parameter, public :: clight = 299792458._rac ![m/s]
    real (kind=rac), parameter, public :: proton_mass = 0.9382720813_rac ![GeV/c^2]

    real (kind=rac), parameter, public :: universal_grav_constant = 6.67408e-20_rac !! CODATA-recommended universal gravitational
                                                                          !! constant \( km^3/kg-s^2  \)

    !> 3x3 identity matrix:
    real (kind=rac), dimension(3,3), parameter, public :: identity_3x3 = reshape(&
                                                         [[one,zero,zero],&
                                                          [zero,one,zero],&
                                                          [zero,zero,one]],[3,3])

    !> 6x6 identity matrix:
    real (kind=rac), dimension(6,6), parameter, public :: identity_6x6= reshape(&
                                                         [[one,zero,zero,zero,zero,zero],&
                                                          [zero,one,zero,zero,zero,zero],&
                                                          [zero,zero,one,zero,zero,zero],&
                                                          [zero,zero,zero,one,zero,zero],&
                                                          [zero,zero,zero,zero,one,zero],&
                                                          [zero,zero,zero,zero,zero,one] ],[6,6])

    end module NUMBERS_CONSTANTS
!*****************************************************************************************
