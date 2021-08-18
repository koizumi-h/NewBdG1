module math_mod
  implicit none
  real(8)   , parameter, public :: pi = 2.0d0 * acos(0.0d0)
  complex(8), parameter, public :: ui = (0.0d0, 1.0d0)
  real(8), dimension(3), parameter :: ex = [1.d0, 0.d0, 0.d0]
  real(8), dimension(3), parameter :: ey = [0.d0, 1.d0, 0.d0]
  real(8), dimension(3), parameter :: ez = [0.d0, 0.d0, 1.d0]
  type, public :: vector
    real(8), dimension(3) :: value
  end type vector
!  real(8), parameter, public :: alpha = 0.0158d0, beta = 0.0910d0
!  real(8), parameter, public :: alpha = 0.0158d0, beta = 0.02884d0
  real(8), parameter :: alpha = 0.00791d0, beta = 0.02884d0
  real(8), parameter :: mf_scale = 2.431d-4
  interface cross_product
    module procedure cross_product1, cross_product2
  end interface
  interface sign2
    module procedure dsign2, isign2
  end interface sign2
  interface swap
    module procedure iswap, iaswap, rswap
  end interface
  interface inv
    module procedure dinv
    module procedure zinv
  end interface
  interface det
    module procedure zdet
  end interface
  interface vecmatmul
    module procedure dvecmatmul
    module procedure zvecmatmul
  end interface
  interface matvecmul
    module procedure dmatvecmul
    module procedure zmatvecmul
  end interface
  interface g_dot_product
    module procedure zg_dot_product
  end interface
  interface lt
    module procedure zlt
  end interface
  interface avg
    module procedure davg, iavg
  end interface
contains
  pure function ue(i) result(r)
    integer, intent(in)   :: i
    integer               :: j
    real(8), dimension(3) :: r
    j = mod(i, 3)
    select case(j)
      case(1); r = ex
      case(2); r = ey
      case(0); r = ez
    end select
  end function ue
  pure function cross_product1(u, v) result(r)
    real(8), dimension(3), intent(in) :: u, v
    real(8), dimension(3)             :: r
    r = [u(2)*v(3)-u(3)*v(2), u(3)*v(1)-u(1)*v(3), u(1)*v(2)-u(2)*v(1)]
  end function cross_product1
  real(8) pure function cross_product2(u, v, i) result(r)
    real(8), dimension(3), intent(in) :: u, v
    integer, intent(in)               :: i
    r = u(mod(i, 3)+1)*v(mod(i+1, 3)+1)-u(mod(i+1, 3)+1)*v(mod(i, 3)+1)
  end function cross_product2
  real(8) pure function vabs(u) result(r)
    real(8), dimension(3), intent(in) :: u
    r = sqrt(dot_product(u, u))
  end function vabs
  pure subroutine iswap(a, b)
    integer, intent(inout) :: a, b
    integer                :: temp
    !if (a > b) then
      temp = a; a = b; b = temp
    !end if
  end subroutine iswap
  subroutine iaswap(int_array)
    integer, dimension(2) :: int_array
    integer               :: temp
    !if (int_array(1) > int_array(2)) then
      temp = int_array(1)
      int_array(1) = int_array(2)
      int_array(2) = temp
    !end if
  end subroutine iaswap
  pure subroutine rswap(a, b)
    real(8), intent(inout) :: a, b
    real(8)                :: temp
    temp = a; a = b; b = temp
  end subroutine rswap
  pure real(8) function dsign2(x, m) result(r)
    real(8), intent(in) :: x
    integer, intent(in) :: m
    select case(m)
      case(:-1); r = -x
      case(0);   r = 0.d0
      case(1:);  r = x
    end select
  end function dsign2
  pure integer function isign2(i, m) result(r)
    integer, intent(in) :: i, m
    select case(m)
      case(:-1); r = -i
      case(0);   r = 0
      case(1:);  r = i
    end select
  end function isign2
  real(8) pure function gauss(x, sigma, xmat, gmat) result(r)
    real(8), intent(in)                        :: x, sigma
    real(8), dimension(:), intent(in)          :: xmat
    real(8), dimension(size(xmat)), intent(in) :: gmat
    integer :: i
    r = sum([(gmat(i)*exp(-(x-xmat(i))**2 / (2.0d0*sigma**2)), &
    & i = 1, size(xmat))]) / (sqrt(2.0d0*pi) * sigma)
  end function gauss
  pure elemental real(8) function principal(x) result(r)
    real(8), intent(in) :: x
    r = modulo(x - pi, 2.0d0*pi) - pi
  end function principal
  pure elemental real(8) function principal2(x) result(r)
    real(8), intent(in) :: x
    r = modulo(x - 0.5d0*pi, pi) - 0.5d0*pi
  end function principal2
  elemental real(8) function principal3(x) result(r)
    real(8), intent(in) :: x
    r = modulo(x, pi*2.d0)
  end function principal3
  subroutine my_zheev(uplo, m, e, c)
  ! solves eigenvalue problem MC=EC using zheev.
    character(1), intent(in)                    :: uplo
    complex(8), dimension(:, :), intent(in)     :: m
    real(8), dimension(size(m, 1)), intent(out) :: e
    complex(8), dimension(size(m, 1), size(m, 2)), &
    & intent(out), optional                     :: c
    complex(8), dimension(2*size(m, 1)-1)       :: work
    real(8), dimension(3*size(m, 1)-2)          :: rwork
    integer                                     :: info
    complex(8), dimension(size(m, 1), size(m, 2)) &
    &                                           :: temp
    if (uplo /= 'u' .and. uplo /= 'l') then
      write(6, *) "uplo should be 'u' or 'l'."
      stop
    end if
    if (size(m, 1) /= size(m, 2)) then
      write(6, *) 'M must be square matrix.'
    else
      temp(:, :) = m(:, :)
      if (present(c)) then
        call zheev('v', uplo, size(temp, 1), temp, size(temp, 1), &
        & e, work, size(work), rwork, info)
        c(:, :) = temp(:, :)
      else
        call zheev('n', uplo, size(temp, 1), temp, size(temp, 1), &
        & e, work, size(work), rwork, info)
      end if
      if (info /= 0) write(6, *) 'info = ', info
    end if
  end subroutine my_zheev
  subroutine my_zhegv(uplo, m, s, e, c, lu)
  ! solves generalized eigenvalue problem MC=ESC using zhegv.
    character(1), intent(in)                                   :: uplo
    complex(8), dimension(:, :), intent(in)                    :: m
    complex(8), dimension(size(m, 1), size(m, 2)), intent(in)  :: s
    real(8), dimension(size(m, 1)), intent(out)                :: e
    complex(8), dimension(size(m, 1), size(m, 2)), intent(out) :: c, lu
    complex(8), dimension(2*size(m, 1)-1)                      :: work
    real(8), dimension(3*size(m, 1)-2)                         :: rwork
    integer                                                    :: info
    if (uplo /= 'u' .and. uplo /= 'l') then
      write(6, *) "uplo should be 'u' or 'l'."
      stop
    end if
    if (size(m, 1) /= size(m, 2)) then
      write(6, *) 'M must be square matrix.'
    else
      e(:) = 0.d0
      c(:, :) = m(:, :)
      lu(:, :) = s(:, :)
      call zhegv(1, 'v', 'l', size(m, 1), c, size(m, 1), lu, &
      & size(m, 1), e, work, size(work), rwork, info)
      if (info /= 0) write(6, *) 'info = ', info
      select case(uplo)
        case('u'); lu(:, :) = lu(:, :)
        case('l'); lu(:, :) = lt(lu)
      end select
    end if
  end subroutine my_zhegv
  subroutine my_dsysv(m, c, w)
  ! solves simultaneous equations problem MC=W (M is symmetric matrix)
    real(8), dimension(:, :), intent(in) :: m
    real(8), dimension(size(m,1)), intent(inout)       :: w
    real(8), dimension(size(m,1)), intent(out)         :: c
    real(8), dimension(size(m,1))                      :: temp
    real(8), dimension(size(m,1), size(m,1))             :: mm
    integer, dimension(size(m,1))                      :: ipiv
    real(8), dimension(size(m,1)**2)                   :: work
    integer                                          :: info, lwork
    lwork = size(m,1)**2
    temp = w
    mm = m
    call dsysv('u', size(m,1), 1, mm, size(m,1), ipiv, temp, &
    & size(m,1), work, lwork, info)
    c = temp
  end subroutine my_dsysv
  subroutine my_dgetrs(A, b, x)
    real(8), dimension(:, :), intent(in) :: A
    real(8), dimension(size(A,1)), intent(in)     :: b
    real(8), dimension(size(A,1)), intent(out)    :: x
    real(8), dimension(size(A,1), size(A,1))      :: At
    integer, dimension(size(A,1))                 :: ipiv
    integer                                       :: info
    real(8), dimension(size(A,1))                 :: tmp
    integer                                       :: n
    At = A
    tmp = b
    n = size(A,1)
    call dgetrf(n, n, At, n, ipiv, info)
    call dgetrs('N', n, 1, At, n, ipiv, tmp, n, info)
    x = tmp
  end subroutine my_dgetrs
  subroutine my_dgetrs2(A, b, x)
    real(8), dimension(:, :), intent(in) :: A
    real(8), dimension(:, :), intent(in)     :: b
    real(8), dimension(size(A,1), size(b,2)), intent(out)    :: x
    real(8), dimension(size(A,1), size(A,1))      :: At
    integer, dimension(size(A,1))                 :: ipiv
    integer                                       :: info
    real(8), dimension(size(A,1), size(b,2))      :: tmp
    integer                                       :: n, m
    At = A
    tmp = b
    n = size(A,1)
    m = size(b,2)
    if(size(A,1) /= size(A,2) .or. size(A,1) /= size(B,1)) then
      stop "ERROR : incorrect dimensions for martix"
    end if
    call dgetrf(n, n, At, n, ipiv, info)
    call dgetrs('N', n, m, At, n, ipiv, tmp, n, info)
    x = tmp
  end subroutine my_dgetrs2
  subroutine my_zgetrs(A, b, x)
    complex(8), dimension(:, :), intent(in) :: A
    complex(8), dimension(size(A,1)), intent(in)     :: b
    complex(8), dimension(size(A,1)), intent(out)    :: x
    complex(8), dimension(size(A,1), size(A,1))      :: At
    integer, dimension(size(A,1))                 :: ipiv
    integer                                       :: info
    complex(8), dimension(size(A,1))                 :: tmp
    integer                                       :: n
    At = A
    tmp = b
    n = size(A,1)
    call zgetrf(n, n, At, n, ipiv, info)
    call zgetrs('N', n, 1, At, n, ipiv, tmp, n, info)
    x = tmp
  end subroutine my_zgetrs
  subroutine my_zgetrs2(A, b, x)
    complex(8), dimension(:, :), intent(in) :: A
    complex(8), dimension(:, :), intent(in)     :: b
    complex(8), dimension(size(A,1), size(b,2)), intent(out)    :: x
    complex(8), dimension(size(A,1), size(A,1))      :: At
    integer, dimension(size(A,1))                 :: ipiv
    integer                                       :: info
    complex(8), dimension(size(A,1), size(b,2))      :: tmp
    integer                                       :: n, m
    At = A
    tmp = b
    n = size(A,1)
    m = size(b,2)
    if(size(A,1) /= size(A,2) .or. size(A,1) /= size(B,1)) then
      stop "ERROR : incorrect dimensions for martix"
    end if
    call zgetrf(n, n, At, n, ipiv, info)
    call zgetrs('N', n, m, At, n, ipiv, tmp, n, info)
    x = tmp
  end subroutine my_zgetrs2
  function dinv(m) result(r)
  ! returns the inverse matrix of the real matrix M.
    real(8), dimension(:, :), intent(in)       :: m
    real(8), dimension(size(m, 1), size(m, 2)) :: r
    real(8), dimension(size(m, 1))             :: work
    integer, dimension(size(m, 1))                :: ipiv
    integer                                       :: info
    if (size(m, 1) /= size(m, 2)) then
      write(6, *) 'A must be square matrix.'
    else
      r(:, :) = m(:, :)
      call dgetrf(size(m, 1), size(m, 1), r, size(m, 1), ipiv, info)
      call dgetri(size(m, 1), r, size(m, 1), ipiv, work, size(work), info)
      if (info /= 0) write(6, *) 'info = ', info
    end if
  end function dinv
  function zinv(m) result(r)
  ! returns the inverse matrix of the complex matrix M.
    complex(8), dimension(:, :), intent(in)       :: m
    complex(8), dimension(size(m, 1), size(m, 2)) :: r
    complex(8), dimension(size(m, 1))             :: work
    integer, dimension(size(m, 1))                :: ipiv
    integer                                       :: info
    if (size(m, 1) /= size(m, 2)) then
      write(6, *) 'A must be square matrix.'
    else
      r(:, :) = m(:, :)
      call zgetrf(size(m, 1), size(m, 1), r, size(m, 1), ipiv, info)
      call zgetri(size(m, 1), r, size(m, 1), ipiv, work, size(work), info)
      if (info /= 0) write(6, *) 'info = ', info
    end if
  end function zinv
  complex(8) function zdet(m) result(r)
  ! returns the determinant of the complex matrix M.
    complex(8), dimension(:, :), intent(in)       :: m
    complex(8), dimension(size(m, 1), size(m, 2)) :: temp
    integer, dimension(size(m, 1))                :: ipiv
    integer                                       :: info, i
    if (size(m, 1) /= size(m, 2)) then
      write(6, *) 'A must be square matrix.'
    else
      temp(:, :) = m(:, :)
      call zgetrf(size(m, 1), size(m, 1), temp, size(m, 1), ipiv, info)
      if (info /= 0) write(6, *) 'info = ', info
      r = product([(temp(i, i), i = 1, size(m, 1))])
    end if
  end function zdet
  function dmatvecmul(m, v) result(r)
  ! returns M * V**H.
    real(8), dimension(:, :), intent(in)       :: m
    real(8), dimension(size(m, 2)), intent(in) :: v
    real(8), dimension(size(m, 1))             :: r
    integer :: i, j
    do i = 1, size(m, 1)
      r(i) = 0.d0
      do j = 1, size(m, 2)
        r(i) = r(i) + m(i, j) * v(j)
      end do
    end do
  end function dmatvecmul
  function zmatvecmul(m, v) result(r)
  ! returns M * V**H.
    complex(8), dimension(:, :), intent(in)       :: m
    complex(8), dimension(size(m, 2)), intent(in) :: v
    complex(8), dimension(size(m, 1))             :: r
    integer :: i, j
    do i = 1, size(m, 1)
      r(i) = (0.d0, 0.d0)
      do j = 1, size(m, 2)
        r(i) = r(i) + m(i, j) * v(j)
      end do
    end do
  end function zmatvecmul
  function dvecmatmul(v, m) result(r)
  ! returns V**H * M.
    real(8), dimension(:), intent(in)    :: v
    real(8), dimension(size(v), size(v)) :: m
    real(8), dimension(size(v))          :: r
    integer                                 :: i, j
    do i = 1, size(v)
      r(i) = 0.d0
      do j = 1, size(v)
        r(i) = r(i) + v(j) * m(j, i)
      end do
    end do
    r(:) = r(:)
  end function dvecmatmul
  function zvecmatmul(v, m) result(r)
  ! returns V**H * M.
    complex(8), dimension(:), intent(in)    :: v
    complex(8), dimension(size(v), size(v)) :: m
    complex(8), dimension(size(v))          :: r
    integer                                 :: i, j
    do i = 1, size(v)
      r(i) = (0.d0, 0.d0)
      do j = 1, size(v)
        r(i) = r(i) + conjg(v(j)) * m(j, i)
      end do
    end do
    r(:) = conjg(r(:))
  end function zvecmatmul
  real(8) function zg_dot_product(v1, v2, s) result(r)
  ! returns V1**H * V2 where S is overlap matrix.
    complex(8), dimension(:), intent(in)                  :: v1
    complex(8), dimension(size(v1)), intent(in)           :: v2
    complex(8), dimension(size(v1), size(v1)), intent(in) :: s
    complex(8)                                            :: temp
    integer                                               :: i, j
    temp = (0.d0, 0.d0)
    do i = 1, size(v1)
      do j = 1, size(v1)
        temp = temp + conjg(v1(i)) * v2(j) * s(i, j)
      end do
    end do
    r = real(temp, 8)
  end function zg_dot_product
  function hc(m) result(r)
  ! returns the Hermitian conjugate of the complex matrix M.
    complex(8), dimension(:, :), intent(in)       :: m
    complex(8), dimension(size(m, 1), size(m, 2)) :: r
    if (size(m, 1) /= size(m, 2)) then
      write(6, *) 'A must be square matrix.'
    else
      r(:, :) = conjg(reshape(m, shape(m), order = [2, 1]))
    end if
  end function hc
  function zlt(m) result(r)
  ! returns the lower triangular matrix of the complex matrix M.
    complex(8), dimension(:, :), intent(in)       :: m
    complex(8), dimension(size(m, 1), size(m, 2)) :: r
    integer                                       :: i
    if (size(m, 1) /= size(m, 2)) then
      write(6, *) 'A must be square matrix.'
    else
      r(:, :) = m(:, :)
      forall(i=1:size(m, 1)-1) r(i, i+1:size(m, 1)) = 0.d0
    end if
  end function zlt
  function davg(v) result(r)
    real(8), dimension(:), intent(in)       :: v
    real(8)                                 :: r
    r = sum(v)/real(size(v),kind(0d0))
  end function davg
  function iavg(v) result(r)
    integer(8), dimension(:), intent(in)       :: v
    real(8)                                 :: r
    r = sum(v)/real(size(v),kind(0d0))
  end function iavg
  function serror(x, avgx) result(r)
    real(8), dimension(:), intent(in)       :: x
    real(8), intent(in), optional           :: avgx
    real(8)                                 :: temp, r
    if (present(avgx)) then
      temp = avgx
    else
      temp = davg(x)
    end if
    r = sum((x-temp)**2)/real(size(x)-1,kind(0d0))
  end function serror
  function get_ixiy(nx, site_pos) result(r)
    integer, intent(in) :: nx, site_pos
    integer, dimension(2)  :: r
    r(1) = mod(site_pos, nx) + 1
    r(2) = (site_pos - r(1)) / nx
  end function get_ixiy
  function get_site_pos(nx, ixiy) result(r)
    integer, intent(in)               :: nx
    integer, intent(in), dimension(2) ::  ixiy
    integer                           :: r
    r = ixiy(1) + ixiy(2)*nx
  end function get_site_pos

  pure elemental function fermi(e,kbT) result(f)
    implicit none
    real(8) f
    real(8),intent(in) ::e,kbT
    if((e)/kbT > 50) then
      f = 0.d0
    else 
      f = 1.d0/(exp((e)/kbT)+1.d0)
    end if
  end function fermi

  function dfermi(e,kbT) result(f)
    implicit none
    real(8) f
    real(8),intent(in) ::e,kbT
      !f = 1/(kbt*(cosh(e/(2*kbt)))**2)
      if( abs(e/kbt) > 700d0 ) then  ! more about 700, cosh occur buffer-overflow
        f = 0d0
      else
        f = -0.25d0/(kbt*(cosh(0.5d0*e/kbt))**2d0)
      end if
  end function dfermi
end module math_mod
