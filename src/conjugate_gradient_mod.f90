module conjugate_gradient_mod
  USE nrtype
  implicit none
  !******* nr variables ********!
  !*****************************!
  type, abstract, public :: conjugate_gradient
    real(8), allocatable    :: ftol
    real(8), dimension(:), pointer       :: pcom, xicom
    integer                              :: ncom
  contains
    procedure :: conjugate_gradient_init
    procedure(abstract_get_func), deferred :: get_func
    procedure(abstract_get_dfunc), deferred :: get_dfunc
    procedure :: brent
    procedure :: mnbrak
    procedure :: linmin
    procedure :: f1dim
    procedure :: frprmn
    procedure :: shft
  end type conjugate_gradient
  abstract interface
    function abstract_get_func(this, val) result(r)
      import :: conjugate_gradient
      class(conjugate_gradient), intent(in)  :: this
      real(8), dimension(:), intent(in)      :: val
      real(8)                                :: r
    end function abstract_get_func
    function abstract_get_dfunc(this, val) result(r)
      import :: conjugate_gradient
      class(conjugate_gradient), intent(in)  :: this
      real(8), dimension(:), intent(in)      :: val
      real(8), dimension(size(val))          :: r
    end function abstract_get_dfunc
  end interface
  !*****************************!
contains
  !*** nr module subroutines ***!
  FUNCTION brent(this, ax,bx,cx,tol,xmin)
    USE nrtype; USE nrutil, ONLY : nrerror
    IMPLICIT NONE
    class(conjugate_gradient), intent(inout)  :: this
    REAL(DP), INTENT(IN) :: ax,bx,cx,tol
    REAL(DP), INTENT(OUT) :: xmin
    REAL(DP) :: brent
    INTEGER(I4B), PARAMETER :: ITMAX=100
    REAL(DP), PARAMETER :: CGOLD=0.3819660_dp, ZEPS=1.0e-3_dp*epsilon(ax)
    INTEGER(I4B) :: iter
    REAL(DP) :: a,b,d,e,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm
    a=min(ax,cx)
    b=max(ax,cx)
    v=bx
    w=v
    x=v
    e=0.0
    fx=this%f1dim(x)
    fv=fx
    fw=fx
    do iter=1,ITMAX
      xm=0.5_dp*(a+b)
      tol1=tol*abs(x)+ZEPS
      tol2=2.0_dp*tol1
      if (abs(x-xm) <= (tol2-0.5_dp*(b-a))) then
        xmin=x
        brent=fx
        RETURN
      end if
      if (abs(e) > tol1) then
        r=(x-w)*(fx-fv)
        q=(x-v)*(fx-fw)
        p=(x-v)*q-(x-w)*r
        q=2.0_dp*(q-r)
        if (q > 0.0) p=-p
        q=abs(q)
        etemp=e
        e=d
        if (abs(p) >= abs(0.5_dp*q*etemp) .or. &
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
      fu=this%f1dim(u)
      if (fu <= fx) then
        if (u >= x) then
          a=x
        else
          b=x
        end if
        call this%shft(v,w,x,u)
        call this%shft(fv,fw,fx,fu)
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
    call nrerror('brent: exceed maximum iterations')
  END FUNCTION brent
  SUBROUTINE shft(this,a,b,c,d)
    class(conjugate_gradient), intent(inout)  :: this
    REAL(DP), INTENT(OUT) :: a
    REAL(DP), INTENT(INOUT) :: b,c
    REAL(DP), INTENT(IN) :: d
    a=b
    b=c
    c=d
  END SUBROUTINE shft
  SUBROUTINE mnbrak(this,ax,bx,cx,fa,fb,fc)
    USE nrtype; USE nrutil, ONLY : swap
    class(conjugate_gradient), intent(inout)  :: this
    REAL(DP), INTENT(INOUT) :: ax,bx
    REAL(DP), INTENT(OUT) :: cx,fa,fb,fc
    REAL(DP), PARAMETER :: GOLD=1.618034_dp,GLIMIT=100.0_dp,TINY=1.0e-20_dp
    REAL(DP) :: fu,q,r,u,ulim
    fa=this%f1dim(ax)
    fb=this%f1dim(bx)
    if (fb > fa) then
      call swap(ax,bx)
      call swap(fa,fb)
    end if
    cx=bx+GOLD*(bx-ax)
    fc=this%f1dim(cx)
    do
!      write(6, *) ax,bx,cx,fa,fb,fc !
      if (fb < fc) RETURN
      r=(bx-ax)*(fb-fc)
      q=(bx-cx)*(fb-fa)
      u=bx-((bx-cx)*q-(bx-ax)*r)/(2.0_dp*sign(max(abs(q-r),TINY),q-r))
      ulim=bx+GLIMIT*(cx-bx)
      if ((bx-u)*(u-cx) > 0.0) then
        fu=this%f1dim(u)
        if (fu < fc) then
          ax=bx
          fa=fb
          bx=u
          fb=fu
 !         write(6, *) ax,bx,cx,fa,fb,fc !
          RETURN
        else if (fu > fb) then
          cx=u
          fc=fu
!          write(6, *) ax,bx,cx,fa,fb,fc !
          RETURN
        end if
        u=cx+GOLD*(cx-bx)
        fu=this%f1dim(u)
      else if ((cx-u)*(u-ulim) > 0.0) then
        fu=this%f1dim(u)
        if (fu < fc) then
          bx=cx
          cx=u
          u=cx+GOLD*(cx-bx)
          call this%shft(fb,fc,fu,this%f1dim(u))
        end if
      else if ((u-ulim)*(ulim-cx) >= 0.0) then
        u=ulim
        fu=this%f1dim(u)
      else
        u=cx+GOLD*(cx-bx)
        fu=this%f1dim(u)
      end if
      call this%shft(ax,bx,cx,u)
      call this%shft(fa,fb,fc,fu)
    end do
  END SUBROUTINE mnbrak
  FUNCTION f1dim(this,x)
    class(conjugate_gradient), intent(in)  :: this
    REAL(DP), INTENT(IN) :: x
    REAL(DP) :: f1dim
    REAL(DP), DIMENSION(:), ALLOCATABLE :: xt
    allocate(xt(this%ncom))
    xt(:)=this%pcom(:)+x*this%xicom(:)
    f1dim=this%get_func(xt)
    deallocate(xt)
  END FUNCTION f1dim
  SUBROUTINE linmin(this,p,xi,fret)
    USE nrtype; USE nrutil, ONLY : assert_eq
    class(conjugate_gradient), intent(inout)  :: this
    REAL(DP), INTENT(OUT) :: fret
    REAL(DP), DIMENSION(:), TARGET, INTENT(INOUT) :: p,xi
    REAL(DP), PARAMETER :: TOL=1.0e-4_dp
    REAL(DP) :: ax,bx,fa,fb,fx,xmin,xx
    this%ncom=assert_eq(size(p),size(xi),'linmin')
    this%pcom=>p
    this%xicom=>xi
    ax=0.0
    xx=1.0
    call this%mnbrak(ax,xx,bx,fa,fx,fb)
    fret=this%brent(ax,xx,bx,TOL,xmin)
    xi=xmin*xi
    p=p+xi
  END SUBROUTINE linmin
  SUBROUTINE frprmn(this,p,ftol,iter,fret,optimized)
    USE nrtype; USE nrutil, ONLY : nrerror
    class(conjugate_gradient), intent(inout)  :: this
    INTEGER(I4B), INTENT(OUT) :: iter
    REAL(DP), INTENT(IN) :: ftol
    REAL(DP), INTENT(OUT) :: fret
    REAL(DP), DIMENSION(:), INTENT(INOUT) :: p
    INTEGER(I4B), PARAMETER :: ITMAX=100 !<- changed
    REAL(DP), PARAMETER :: EPS=1.0e-10_dp
    INTEGER(I4B) :: its
    REAL(DP) :: dgg,fp,gam,gg
    REAL(DP), DIMENSION(size(p)) :: g,h,xi
    integer :: i !
    logical, intent(out), optional     :: optimized
    fp=this%get_func(p)
    xi=this%get_dfunc(p)
    g=-xi
    h=g
    xi=h
    do its=1,ITMAX
      iter=its
      !write(*,*)'xi1',p
      !write(*,*)'func',fp
      !write(*,*)'dfunc',xi
      !write(*,*)''
      call this%linmin(p,xi,fret)
      if (2.0_dp*abs(fret-fp) <= ftol*(abs(fret)+abs(fp)+EPS)) then
        if(present(optimized)) then
          optimized = .true.
        end if
        RETURN
      endif
      fp=fret
      xi=this%get_dfunc(p)
      gg=dot_product(g,g)
      dgg=dot_product(xi+g,xi)
      if (gg == 0.0) RETURN
      gam=dgg/gg
      g=-xi
      h=g+gam*h
      xi=h
    end do
    if(present(optimized)) then 
      optimized = .false.
    end if
    !call nrerror('frprmn: maximum iterations exceeded')
  END SUBROUTINE frprmn
  !*****************************!
  ! module subprograms
  ! public type-bound procedures
  subroutine conjugate_gradient_init(this)
    class(conjugate_gradient),intent(inout)  :: this
    if (allocated(this%ftol))deallocate(this%ftol)
    allocate(this%ftol)
    this%ftol = 1.0d-16
  end subroutine conjugate_gradient_init

end module conjugate_gradient_mod
