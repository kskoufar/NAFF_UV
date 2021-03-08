!*****************************************************************************************
!> author: Jacob Williams
!> modified by: Kyriacos Skoufaris
!
!  Brent algorithms for minimization and root solving without derivatives.
!
!# See also
!
!  1. R. Brent, "Algorithms for Minimization Without Derivatives",
!     Prentice-Hall, Inc., 1973.

    module BRENT_MNBRAK_METHODS

    use KIND_ACCURACY
    use NUMBERS_CONSTANTS

    implicit none

    private

    type,public :: brent_mnbrak_class
        !! the main class
        procedure(func),pointer :: f => null()  !! function to be minimized
    contains
        procedure :: set_function
        procedure :: braket_min => mnbrak
        procedure :: global_min => fglomin
        procedure :: local_min => flocmin
        procedure :: find_zero => zeroin
    end type brent_mnbrak_class

    abstract interface
        function func(me,x) result(f)    !! Interface to the function to be minimized.
                                         !! It should evaluate f(x) for any x in the interval (ax,bx)
            import :: rac,brent_mnbrak_class
            implicit none
            class(brent_mnbrak_class),intent(inout) :: me
            real(rac),intent(in) :: x
            real(rac) :: f
        end function func
    end interface

    !unit test routine:
    public :: brent_test

    contains
!*****************************************************************************************

!*****************************************************************************************
!> author: Jacob Williams
!  date: 7/19/2014
!> modified by: Kyriacos Skoufaris
!
!  Set the function to be minimized.

    subroutine set_function(me,f)

    implicit none

    class (brent_mnbrak_class), intent(inout) :: me
    procedure (func) :: f

    me%f => f

    end subroutine set_function
!*****************************************************************************************

!*****************************************************************************************
    SUBROUTINE fglomin(me,ax,bx,cx,tol,xmin,f_at_xmin)

    IMPLICIT NONE
    CLASS (brent_mnbrak_class), INTENT(INOUT) :: me
    REAL (kind=rac), INTENT(IN) :: ax,bx,cx,tol
    REAL (kind=rac), INTENT(OUT) :: xmin, f_at_xmin
    INTEGER (kind=iac), PARAMETER :: ITMAX=1000
    REAL (kind=rac), PARAMETER :: CGOLD=0.5_rac*(3.0_rac-sqrt(5.0_rac)),ZEPS=1.0e-3_rac*epsilon(ax)
    INTEGER (kind=iac) :: iter
    REAL (kind=rac) :: a,b,d,e,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm

    a=min(ax,cx)
    b=max(ax,cx)
    v=bx
    w=v
    x=v
    e=0.0_rac
    fx=me%f(x)
    fv=fx
    fw=fx
    do iter=1,ITMAX
      xm=0.5_rac*(a+b)
      tol1=tol*abs(x)+ZEPS
      tol2=2.0_rac*tol1
      if (abs(x-xm) <= (tol2-0.5_rac*(b-a))) then
        xmin=x
        f_at_xmin=fx
        RETURN
      end if
      if (abs(e) > tol1) then
        r=(x-w)*(fx-fv)
        q=(x-v)*(fx-fw)
        p=(x-v)*q-(x-w)*r
        q=2.0_rac*(q-r)
        if (q > 0.0_rac) p=-p
        q=abs(q)
        etemp=e
        e=d
        if (abs(p) >= abs(0.5_rac*q*etemp) .or. &
          p <= q*(a-x) .or. p >= q*(b-x)) then
          e=merge(a-x,b-x, x >= xm )
          d=CGOLD*e
        else
          d=p/q
          u=x+d
          if (u-a < tol2 .or. b-u < tol2) d=sign(tol1,xm-x)
        end if
      else
        e=merge(a-x,b-x, x >= xm )
        d=CGOLD*e
      end if
      u=merge(x+d,x+sign(tol1,d), abs(d) >= tol1 )
      fu=me%f(u)
      if (fu <= fx) then
        if (u >= x) then
          a=x
        else
          b=x
        end if
        call shft(v,w,x,u)
        call shft(fv,fw,fx,fu)
      else
        if (u < x) then
          a=u
        else
          b=u
        end if
        if (fu <= fw .or. w == x) then
          v=w
          fv=fw
          w=u
          fw=fu
        else if (fu <= fv .or. v == x .or. v == w) then
          v=u
          fv=fu
        end if
      end if
    end do

    print*, "ERROR - FUNCTION fglomin !!!"
    print*, 'brent (fglomin): exceeded maximum iterations before it finds the global minimum'

    CONTAINS
!BL
    SUBROUTINE shft(a,b,c,d)
    REAL (kind=rac), INTENT(OUT) :: a
    REAL (kind=rac), INTENT(INOUT) :: b,c
    REAL (kind=rac), INTENT(IN) :: d
    a=b
    b=c
    c=d
    END SUBROUTINE shft

    END SUBROUTINE fglomin
!*****************************************************************************************


!*****************************************************************************************
!>
!  An approximation x to the point where f attains a local minimum on
!  the interval (ax,bx) is determined.
!
!  the method used is a combination of golden section search and
!  successive parabolic interpolation. convergence is never much slower
!  than that for a fibonacci search. if f has a continuous second
!  derivative which is positive at the minimum (which is not at ax or
!  bx), then convergence is superlinear, and usually of the order of
!  about 1.324.
!
!  the function f is never evaluated at two points closer together
!  than eps*abs(flocmin) + (tol/3), where eps is approximately the square
!  root of the relative machine precision. if f is a unimodal
!  function and the computed values of f are always unimodal when
!  separated by at least eps*abs(x) + (tol/3), then flocmin approximates
!  the abcissa of the global minimum of f on the interval ax,bx with
!  an error less than 3*eps*abs(flocmin) + tol. if f is not unimodal,
!  then flocmin may approximate a local, but perhaps non-global, minimum to
!  the same accuracy.
!
!  this function subprogram is a slightly modified version of the
!  algol 60 procedure localmin given in richard brent, algorithms for
!  minimization without derivatives, prentice - hall, inc. (1973).
!
!# See also
!  [1] http://www.netlib.no/netlib/fmm/fmin.f fmin.f->flocmin

    SUBROUTINE flocmin(me,ax,bx,tol,xmin,fxmin)

    implicit none

    class(brent_mnbrak_class),intent(inout) :: me
    real (kind=rac), intent(in)  :: ax    !! left endpoint of initial interval
    real (kind=rac), intent(in)  :: bx    !! right endpoint of initial interval
    real (kind=rac), intent(in)  :: tol   !! desired length of the interval of uncertainty of the final result (>=0)
    real (kind=rac), intent(out) :: xmin  !! abcissa approximating the point where f attains a minimum
    real (kind=rac), intent(out) :: fxmin !! f at xmin
    
    real (kind=rac) :: a,b,d,e,xm,p,q,r,tol1,tol2,u,v,w
    real (kind=rac) :: fu,fv,fw,fx,x
    real (kind=rac) :: abs,sqrt,sign

    real (kind=rac), parameter :: half = 0.5_rac
    real (kind=rac), parameter :: c = (three-sqrt(five))*half    !! squared inverse of golden ratio
    real (kind=rac), parameter :: eps = epsilon(one) !sqrt(epsilon(one))

    !initialization

    a = ax
    b = bx
    v = a + c*(b - a)
    w = v
    x = v
    e = zero
    fx = me%f(x)
    fv = fx
    fw = fx

    do    !  main loop starts here

        xm = half*(a + b)
        tol1 = eps*abs(x) + tol/three
        tol2 = two*tol1
        

        !  check stopping criterion

        if (abs(x - xm) <= (tol2 - half*(b - a))) exit

        ! is golden-section necessary

        if (abs(e) <= tol1) go to 40

        !  fit parabola

        r = (x - w)*(fx - fv)
        q = (x - v)*(fx - fw)
        p = (x - v)*q - (x - w)*r
        q = two*(q - r)
        if (q > zero) p = -p
        q =  abs(q)
        r = e
        e = d

        !  is parabola acceptable

        if (abs(p) >= abs(half*q*r)) go to 40
        if (p <= q*(a - x)) go to 40
        if (p >= q*(b - x)) go to 40

        !  a parabolic interpolation step

        d = p/q
        u = x + d

        !  f must not be evaluated too close to ax or bx

        if ((u - a) < tol2) d = sign(tol1, xm - x)
        if ((b - u) < tol2) d = sign(tol1, xm - x)
        go to 50

        !  a golden-section step

    40  if (x >= xm) e = a - x
        if (x < xm) e = b - x
        d = c*e

        !  f must not be evaluated too close to x

    50  if (abs(d) >= tol1) u = x + d
        if (abs(d) < tol1) u = x + sign(tol1, d)
        fu = me%f(u)

        !  update a, b, v, w, and x

        if (fu <= fx) then
            if (u >= x) a = x
            if (u < x) b = x
            v = w
            fv = fw
            w = x
            fw = fx
            x = u
            fx = fu
            cycle
        end if

        if (u < x) a = u
        if (u >= x) b = u

        if (fu <= fw .or. w == x) then
            v = w
            fv = fw
            w = u
            fw = fu
            cycle
        end if

        if (fu <= fv .or. v == x .or. v == w ) then
            v = u
            fv = fu
            cycle
        end if

    end do    !  end of main loop

    xmin = x
    fxmin = fx

    END SUBROUTINE flocmin
!*****************************************************************************************

!*****************************************************************************************
!>
!  Find a zero of the function \( f(x) \) in the given interval
!  \( [a_x,b_x] \) to within a tolerance \( 4 \epsilon |x| + tol \),
!  where \( \epsilon \) is the relative machine precision defined as
!  the smallest representable number such that \( 1.0 + \epsilon > 1.0 \).
!
!  It is assumed that \( f(a_x) \) and \( f(b_x) \) have opposite signs.
!
!#References
!  * R. P. Brent, "[An algorithm with guaranteed convergence for
!    finding a zero of a function](http://maths-people.anu.edu.au/~brent/pd/rpb005.pdf)",
!    The Computer Journal, Vol 14, No. 4., 1971.
!  * R. P. Brent, "[Algorithms for minimization without derivatives](http://maths-people.anu.edu.au/~brent/pub/pub011.html)",
!    Prentice-Hall, Inc., 1973.
!
!# See also
!  1. [zeroin.f](http://www.netlib.org/go/zeroin.f) from Netlib

    subroutine zeroin(me,ax,bx,tol,xzero,fzero,iflag,fax,fbx)

    use iso_fortran_env, only: error_unit

    implicit none

    class(brent_mnbrak_class),intent(inout) :: me
    real (kind=rac), intent(in) :: ax      !! left endpoint of initial interval
    real (kind=rac), intent(in) :: bx      !! right endpoint of initial interval
    real (kind=rac), intent(in) :: tol     !! desired length of the interval of uncertainty of the final result (>=0)
    real (kind=rac), intent(out) :: xzero   !! abscissa approximating a zero of `f` in the interval `ax`,`bx`
    real (kind=rac), intent(out) :: fzero   !! value of `f` at the root (`f(xzero)`)
    integer (kind=iac), intent(out)      :: iflag   !! status flag (`-1`=error, `0`=root found)
    real (kind=rac), intent(in),optional :: fax     !! if `f(ax)` is already known, it can be input here
    real (kind=rac), intent(in),optional :: fbx     !! if `f(ax)` is already known, it can be input here

    real (kind=rac), parameter :: eps   = epsilon(one)  !! original code had d1mach(4)
    real (kind=rac) :: a,b,c,d,e,fa,fb,fc,tol1,xm,p,q,r,s

    tol1 = eps+one

    a=ax
    b=bx

    if (present(fax)) then
        fa = fax
    else
        fa=me%f(a)
    end if
    if (present(fbx)) then
        fb = fbx
    else
        fb=me%f(b)
    end if

    !check trivial cases first:
    if (fa==zero) then

        iflag = 0
        xzero = a
        fzero = fa

    elseif (fb==zero) then

        iflag = 0
        xzero = b
        fzero = fb

    elseif (fa*(fb/abs(fb))<zero) then  ! check that f(ax) and f(bx) have different signs

        c=a
        fc=fa
        d=b-a
        e=d

        do

            if (abs(fc)<abs(fb)) then
                a=b
                b=c
                c=a
                fa=fb
                fb=fc
                fc=fa
            end if

            tol1=two*eps*abs(b)+0.5_rac*tol
            xm = 0.5_rac*(c-b)
            if ((abs(xm)<=tol1).or.(fb==zero)) exit

            ! see if a bisection is forced
            if ((abs(e)>=tol1).and.(abs(fa)>abs(fb))) then
                s=fb/fa
                if (a/=c) then
                    ! inverse quadratic interpolation
                    q=fa/fc
                    r=fb/fc
                    p=s*(two*xm*q*(q-r)-(b-a)*(r-one))
                    q=(q-one)*(r-one)*(s-one)
                else
                    ! linear interpolation
                    p=two*xm*s
                    q=one-s
                end if
                if (p<=zero) then
                    p=-p
                else
                    q=-q
                end if
                s=e
                e=d
                if (((two*p)>=(three*xm*q-abs(tol1*q))) .or. &
                    (p>=abs(0.5_rac*s*q))) then
                    d=xm
                    e=d
                else
                    d=p/q
                end if
            else
                d=xm
                e=d
            end if

            a=b
            fa=fb
            if (abs(d)<=tol1) then
                if (xm<=zero) then
                    b=b-tol1
                else
                    b=b+tol1
                end if
            else
                b=b+d
            end if
            fb=me%f(b)
            if ((fb*(fc/abs(fc)))>zero) then
                c=a
                fc=fa
                d=b-a
                e=d
            end if

        end do

        iflag = 0
        xzero = b
        fzero = fb

    else

        iflag = -1
        write(error_unit,'(A)')&
            'Error in zeroin: f(ax) and f(bx) do not have different signs.'

    end if

    end subroutine zeroin
!*****************************************************************************************

!*****************************************************************************************
!> author: Jacob Williams
!  date: 7/16/2014
!> modified by: Kyriacos Skoufaris
!
!  Test of the flocmin and zeroin functions.

    subroutine brent_test()

    implicit none

    real(rac) :: r,ff_fx,fzero
    integer(iac) :: iflag

    real(rac),parameter :: ax = zero
    real(rac),parameter :: bx = two*pi
    real(rac),parameter :: tol = 1.0e-6_rac

    type,extends(brent_mnbrak_class) :: myfunc_type
        integer :: i = 0    !! function counter
    end type myfunc_type
    type(myfunc_type) :: myfunc

    write(*,*) ''
    write(*,*) '---------------'
    write(*,*) ' brent_test'
    write(*,*) '---------------'
    write(*,*) ''

    call myfunc%set_function(sin_func)    !set the function

    !call flocmin:
    ! [the minimum is at 270 deg]
    myfunc%i = 0
    !r = myfunc%local_min(ax,bx,tol)
    call myfunc%local_min(ax,bx,tol,r,ff_fx)
    write(*,*) 'minimum of sin(x) at: ', r*180.0_rac/pi,' deg'
    write(*,*) 'number of function calls: ', myfunc%i

    !call zeroin:
    ! [the root is at pi]
    myfunc%i = 0
    call myfunc%find_zero(ax+0.0001_rac,bx/two+0.0002,tol,r,fzero,iflag)
    write(*,*) 'root of sin(x) at: ', r*180.0_rac/pi,' deg'
    write(*,*) 'number of function calls: ', myfunc%i

    contains

        function sin_func(me,x) result(f)
        !! Example function to minimize: sin(x)

        implicit none

        class(brent_mnbrak_class),intent(inout) :: me
        real(rac),intent(in) :: x
        real(rac) :: f

        f = sin(x)

        select type (me)
        class is (myfunc_type)
            me%i = me%i + 1 !number of function calls
        end select

        end function sin_func

    end subroutine brent_test
 !*****************************************************************************************



!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! http://jean-pierre.moreau.pagesperso-orange.fr/Fortran/mnbrak_f90.txt
! with modifications from Kyriacos Skoufaris
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!Given a function me -> F(X), and given distinct initial points AX and
!BX, this routine searches in the downhill direction (defined by the
!function as evaluated at the initial points) and returns new points
!AX, BX, CX which bracket a minimum of the function. Also returned
!are the function values at the three points, FA, FB and FC.

    SUBROUTINE mnbrak(me,ax,bx,cx,fa,fb,fc)
    CLASS (brent_mnbrak_class), INTENT(INOUT) :: me
    REAL ( kind = rac ), INTENT(INOUT) :: ax,bx
    REAL ( kind = rac ), INTENT(OUT) :: cx,fa,fb,fc
    REAL ( kind = rac ), PARAMETER :: GOLD=0.5_rac*(1.0_rac+sqrt(5.0_rac)), GLIMIT=100.0_rac,TINY=1.0e-20_rac !1.0e-3_rac*epsilon(ax)
    REAL ( kind = rac ) :: dum,fu,q,r,u,ulim
    LOGICAL DONE
    fa=me%f(ax)
    fb=me%f(bx)
    if(fb>fa) then
      dum=ax
      ax=bx
      bx=dum
      dum=fb
      fb=fa
      fa=dum
    endif
    cx=bx+GOLD*(bx-ax)
    fc=me%f(cx)
    do
      if(fb<fc) exit
      done=.TRUE.
      r=(bx-ax)*(fb-fc)
      q=(bx-cx)*(fb-fa)
      u=bx-((bx-cx)*q-(bx-ax)*r)/(2.0_rac * sign(max(abs(q-r),TINY),q-r))
      ulim=bx+GLIMIT*(cx-bx)
      if((bx-u)*(u-cx)>0.0_rac) then
        fu=me%f(u)
        if(fu<fc) then
          ax=bx
          fa=fb
          bx=u
          fb=fu
          return
        else if(fu>fb) then
          cx=u
          fc=fu
          return
        endif
        u=cx+GOLD*(cx-bx)
        fu=me%f(u)
      else if((cx-u)*(u-ulim)>0.0_rac) then
        fu=me%f(u)
        if(fu<fc) then
          bx=cx
          cx=u
          u=cx+GOLD*(cx-bx)
          fb=fc
          fc=fu
          fu=me%f(u)
        endif
      else if((u-ulim)*(ulim-cx)>=0.0_rac) then
        u=ulim
        fu=me%f(u)
      else
        u=cx+GOLD*(cx-bx)
        fu=me%f(u)
      endif
      if (done) then
        ax=bx
        bx=cx
        cx=u
        fa=fb
        fb=fc
        fc=fu
      else
        done=.FALSE.
      end if
      if (.not.done) exit
    end do
    END SUBROUTINE mnbrak

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!*****************************************************************************************
    end module BRENT_MNBRAK_METHODS
!*****************************************************************************************
