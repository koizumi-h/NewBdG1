module diagonalize_mod
  use math_mod, only : ui, pi,fermi,dfermi, my_zheev, principal
  use path_list_mod, only : path, path_list
  use site_list_mod
  use parameter_mod, only : site_max
  use base_mod
  implicit none
  
  type, extends(base) :: diag
    real(8), dimension(:), allocatable :: eg
    ! complex(8),dimension(:,:), allocatable   :: wf ! already defined in base
    contains
    procedure :: diag_init
    procedure :: diag_surface => diagonalize_surface_hamiltonian
    procedure :: diag_bulk => diagonalize_bulk_hamiltonian
    procedure :: diag_bulk_xi => diagonalize_bulk_xi_hamiltonian
    procedure :: diag_sb_2lay => diagonalize_sb_2lay_hamiltonian
    procedure :: diag_sb_2lay_xi => diagonalize_sb_2lay_xi_hamiltonian
    !procedure :: diag_sb_2lay_xichi => diagonalize_sb_2lay_xichi_hamiltonian
    procedure :: calc_ne => diag_mod_calc_ne
    procedure :: calc_spin => diag_mod_calc_spin
    !procedure :: calc_surface_delta => diag_mod_calc_surface_delta
    procedure :: calc_delta => diag_mod_calc_delta_sb_2lay
    !procedure :: calc_delta_xi => diag_mod_calc_delta_sb_2lay_xi
    procedure :: calc_tdelta_xi => diag_mod_calc_tdelta_sb_2lay_xi
    procedure :: calc_spin_xi => diag_mod_calc_spin_xi
    procedure :: calc_xi => diag_mod_calc_xi
    !procedure :: calc_eg_sb_2lay_xi
  end type diag
  
  contains

  subroutine diag_init(this, pshape, hole, t, u, jd, lm)
    class(diag), intent(inout)        :: this
    integer, dimension(:, :), intent(in)    :: pshape
    integer, dimension(:, :), intent(in) :: hole
    real(8), dimension(:), intent(in)    :: t
    real(8), intent(in)                  :: u, jd, lm
    real(8), dimension(:), allocatable :: dummy_ne, dummy_xi, dummy_chi
    integer :: n, nz
    n = site_max(pshape)
    nz = size(pshape(1,:))
    allocate(dummy_ne(nz), source = 0d0)
    allocate(dummy_xi(n), source = 0d0)
    allocate(dummy_chi(n), source = 0d0)
    call this%base_init(pshape, hole, dummy_ne, t, u, jd, lm, dummy_xi, dummy_chi)
    if(allocated(this%eg))deallocate(this%eg)
    allocate(this%eg(4*this%n_acc))
    if(allocated(this%wf))deallocate(this%wf)
    allocate(this%wf(4*this%n_acc, 4*this%n_acc))
  end subroutine diag_init

  subroutine diagonalize_surface_hamiltonian(this, mu, delta, ne_u, ne_d)
    implicit none
    class(diag), intent(inout)  :: this
    real(8),intent(in)    :: mu
    complex(8), dimension(this%pl_1st%length()), intent(in) :: delta
    real(8), dimension(this%n), intent(in)   :: ne_u, ne_d
    complex(8),DIMENSION(4*this%n,4*this%n)   :: ham
    integer :: i, j, k

    !! error cheack
    if(this%nz /= 1 .or. this%nh /= 0) then
      print *, "nz : ", this%nz
      print *, "nh : ", this%nh
      stop "error : wrong data in diagonalize_surface_hamiltonian"
    end if

    ham = 0d0
    do i = 1, this%n
      ham(4*i-3, 4*i-3) = -mu
      ham(4*i-2, 4*i-2) = -mu
      ham(4*i-1, 4*i-1) =  mu
      ham(4*i  , 4*i  ) =  mu
    end do
    do i = 1, this%pl_1st%length()
      j = this%pl_1st%value(i)%i
      k = this%pl_1st%value(i)%f
      ham(4*j-3, 4*k-3) = -this%t(1)*(1d0-ne_d(j))*(1d0-ne_d(k))
      ham(4*j-2, 4*k-2) = -this%t(1)*(1d0-ne_u(j))*(1d0-ne_u(k))
      ham(4*j-1, 4*k-1) =  this%t(1)*(1d0-ne_d(j))*(1d0-ne_d(k))
      ham(4*j  , 4*k  ) =  this%t(1)*(1d0-ne_u(j))*(1d0-ne_u(k))
      ham(4*k-3, 4*j-3) = -this%t(1)*(1d0-ne_d(j))*(1d0-ne_d(k))
      ham(4*k-2, 4*j-2) = -this%t(1)*(1d0-ne_u(j))*(1d0-ne_u(k))
      ham(4*k-1, 4*j-1) =  this%t(1)*(1d0-ne_d(j))*(1d0-ne_d(k))
      ham(4*k  , 4*j  ) =  this%t(1)*(1d0-ne_u(j))*(1d0-ne_u(k))
    end do
    do i = 1, this%pl_2nd%length()
      j = this%pl_2nd%value(i)%i
      k = this%pl_2nd%value(i)%f
      ham(4*j-3, 4*k-3) = -this%t(2)*(1d0-ne_d(j))*(1d0-ne_d(k))
      ham(4*j-2, 4*k-2) = -this%t(2)*(1d0-ne_u(j))*(1d0-ne_u(k))
      ham(4*j-1, 4*k-1) =  this%t(2)*(1d0-ne_d(j))*(1d0-ne_d(k))
      ham(4*j  , 4*k  ) =  this%t(2)*(1d0-ne_u(j))*(1d0-ne_u(k))
      ham(4*k-3, 4*j-3) = -this%t(2)*(1d0-ne_d(j))*(1d0-ne_d(k))
      ham(4*k-2, 4*j-2) = -this%t(2)*(1d0-ne_u(j))*(1d0-ne_u(k))
      ham(4*k-1, 4*j-1) =  this%t(2)*(1d0-ne_d(j))*(1d0-ne_d(k))
      ham(4*k  , 4*j  ) =  this%t(2)*(1d0-ne_u(j))*(1d0-ne_u(k))
    end do
    do i = 1, this%pl_1st%length()
      j = this%pl_1st%value(i)%i
      k = this%pl_1st%value(i)%f
      ham(4*j-3, 4*k  ) = delta(i)
      ham(4*j-2, 4*k-1) = delta(i)
      ham(4*j-1, 4*k-2) = conjg(delta(i))
      ham(4*j  , 4*k-3) = conjg(delta(i))
      ham(4*k-3, 4*j  ) = delta(i)
      ham(4*k-2, 4*j-1) = delta(i)
      ham(4*k-1, 4*j-2) = conjg(delta(i))
      ham(4*k  , 4*j-3) = conjg(delta(i))
    end do
  
    call my_zheev("l", ham, this%eg, this%wf)

  end subroutine diagonalize_surface_hamiltonian

  !subroutine diagonalize_bulk_hamiltonian(this, mu, xi, chi)
  subroutine diagonalize_bulk_hamiltonian(this, mu, Sx, Sy, Sz)
    implicit none
    class(diag), intent(inout)   :: this
    !REAL(8), DIMENSION(4*this%n_acc), intent(out) :: eg
    !complex(8),DIMENSION(4*this%n_acc, 4*this%n_acc), intent(out)   :: wf
    real(8),intent(in)    :: mu
    !real(8), dimension(this%n), intent(in) :: xi, chi
    real(8), dimension(this%n), intent(in) :: Sx, Sy, Sz
    complex(8),DIMENSION(4*this%n_acc, 4*this%n_acc)   :: ham
    !real(8), dimension(this%n) :: ne_u, ne_d
    integer :: i, j, k, is, js, ks

    !! error cheack
    if(this%nz /= 1) then
      print *, "nz : ", this%nz
      stop "error : wrong data in diagonalize_bulk_hamiltonian"
    end if

    !do i = 1 ,this%n
    !  Sx(i) = 0.5d0*cos(xi(i))
    !  Sy(i) = 0.5d0*sin(xi(i))
    !end do
    !Sz = 0d0

    ham = 0d0
    do i = 1, this%n_acc
      is = this%ets(i)
      ham(4*i-3, 4*i-3) = this%U*(-2d0/3d0*Sz(is)+0.5d0)
      ham(4*i-2, 4*i-2) = this%U*( 2d0/3d0*Sz(is)+0.5d0)
      ham(4*i-3, 4*i-2) = - this%U*2d0/3d0*(Sx(is)-ui*Sy(is))
      ham(4*i-2, 4*i-3) = - this%U*2d0/3d0*(Sx(is)+ui*Sy(is))
      ham(4*i-1, 4*i-1) = - this%U*(-2d0/3d0*Sz(is)+0.5d0) ! - conjg
      ham(4*i  , 4*i  ) = - this%U*( 2d0/3d0*Sz(is)+0.5d0) ! - conjg
      ham(4*i-1, 4*i  ) = conjg(- this%U*2d0/3d0*(Sx(is)-ui*Sy(is))) ! + conjg
      ham(4*i  , 4*i-1) = conjg(- this%U*2d0/3d0*(Sx(is)+ui*Sy(is))) ! + conjg
    end do
    do i = 1, this%n_acc
      ham(4*i-3, 4*i-3) = ham(4*i-3, 4*i-3) - mu
      ham(4*i-2, 4*i-2) = ham(4*i-2, 4*i-2) - mu
      ham(4*i-1, 4*i-1) = ham(4*i-1, 4*i-1) + mu
      ham(4*i  , 4*i  ) = ham(4*i  , 4*i  ) + mu
    end do
    do i = 1, this%pl_1st%length()
      js = this%pl_1st%value(i)%i
      ks = this%pl_1st%value(i)%f
      j = this%ste(js)
      k = this%ste(ks)
      ham(4*j-3, 4*k-3) = -this%t(1)
      ham(4*j-2, 4*k-2) = -this%t(1)
      ham(4*j-1, 4*k-1) =  this%t(1)
      ham(4*j  , 4*k  ) =  this%t(1)
      ham(4*k-3, 4*j-3) = -this%t(1)
      ham(4*k-2, 4*j-2) = -this%t(1)
      ham(4*k-1, 4*j-1) =  this%t(1)
      ham(4*k  , 4*j  ) =  this%t(1)
    end do
    do i = 1, this%pl_2nd%length()
      js = this%pl_2nd%value(i)%i
      ks = this%pl_2nd%value(i)%f
      j = this%ste(js)
      k = this%ste(ks)
      ham(4*j-3, 4*k-3) = -this%t(2)
      ham(4*j-2, 4*k-2) = -this%t(2)
      ham(4*j-1, 4*k-1) =  this%t(2)
      ham(4*j  , 4*k  ) =  this%t(2)
      ham(4*k-3, 4*j-3) = -this%t(2)
      ham(4*k-2, 4*j-2) = -this%t(2)
      ham(4*k-1, 4*j-1) =  this%t(2)
      ham(4*k  , 4*j  ) =  this%t(2)
    end do
    do i = 1, this%pl_h%length()
      js = this%pl_h%value(i)%i
      ks = this%pl_h%value(i)%f
      j = this%ste(js)
      k = this%ste(ks)
      ham(4*j-3, 4*j-3) = ham(4*j-3, 4*j-3) + 0.5d0*this%jd*Sz(ks)
      ham(4*j-2, 4*j-2) = ham(4*j-2, 4*j-2) - 0.5d0*this%jd*Sz(ks)
      ham(4*j-3, 4*j-2) = ham(4*j-3, 4*j-2) + 0.5d0*this%jd*(Sx(ks) - ui*Sy(ks))
      ham(4*j-2, 4*j-3) = ham(4*j-2, 4*j-3) + 0.5d0*this%jd*(Sx(ks) + ui*Sy(ks))
      ham(4*j-1, 4*j-1) = ham(4*j-1, 4*j-1) - 0.5d0*this%jd*Sz(ks) ! -conjg
      ham(4*j  , 4*j  ) = ham(4*j  , 4*j  ) + 0.5d0*this%jd*Sz(ks) ! -conjg
      ham(4*j-1, 4*j  ) = ham(4*j-1, 4*j  ) + conjg( 0.5d0*this%jd*(Sx(ks) - ui*Sy(ks)) ) ! + conjg
      ham(4*j  , 4*j-1) = ham(4*j  , 4*j-1) + conjg( 0.5d0*this%jd*(Sx(ks) + ui*Sy(ks)) ) ! + conjg
      
      ham(4*k-3, 4*k-3) = ham(4*k-3, 4*k-3) + 0.5d0*this%jd*Sz(js)
      ham(4*k-2, 4*k-2) = ham(4*k-2, 4*k-2) - 0.5d0*this%jd*Sz(js)
      ham(4*k-3, 4*k-2) = ham(4*k-3, 4*k-2) + 0.5d0*this%jd*(Sx(js) - ui*Sy(js))
      ham(4*k-2, 4*k-3) = ham(4*k-2, 4*k-3) + 0.5d0*this%jd*(Sx(js) + ui*Sy(js))
      ham(4*k-1, 4*k-1) = ham(4*k-1, 4*k-1) - 0.5d0*this%jd*Sz(js) ! -conjg
      ham(4*k  , 4*k  ) = ham(4*k  , 4*k  ) + 0.5d0*this%jd*Sz(js) ! -conjg
      ham(4*k-1, 4*k  ) = ham(4*k-1, 4*k  ) + conjg( 0.5d0*this%jd*(Sx(js) - ui*Sy(js)) ) ! + conjg
      ham(4*k  , 4*k-1) = ham(4*k  , 4*k-1) + conjg( 0.5d0*this%jd*(Sx(js) + ui*Sy(js)) ) ! + conjg
    end do
  
    call my_zheev("l", ham, this%eg, this%wf)

  end subroutine diagonalize_bulk_hamiltonian
  
  subroutine diagonalize_bulk_xi_hamiltonian(this, mu, Sx, Sy, Sz, xi)
    implicit none
    class(diag), intent(inout)   :: this
    !REAL(8), DIMENSION(4*this%n_acc), intent(out) :: eg
    !complex(8),DIMENSION(4*this%n_acc, 4*this%n_acc), intent(out)   :: wf
    real(8),intent(in)    :: mu
    !real(8), dimension(this%n), intent(in) :: xi, chi
    real(8), dimension(this%n), intent(in) :: Sx, Sy, Sz, xi
    complex(8),DIMENSION(4*this%n_acc, 4*this%n_acc)   :: ham
    !real(8), dimension(this%n) :: ne_u, ne_d
    integer :: i, j, k, is, js, ks
    complex(8) :: exi,cexi

    !! error cheack
    if(this%nz /= 1) then
      print *, "nz : ", this%nz
      stop "error : wrong data in diagonalize_bulk_hamiltonian"
    end if

    !do i = 1 ,this%n
    !  Sx(i) = 0.5d0*cos(xi(i))
    !  Sy(i) = 0.5d0*sin(xi(i))
    !end do
    !Sz = 0d0

    ham = 0d0
    do i = 1, this%n_acc
      is = this%ets(i)
      exi  = exp(ui*xi(is))
      cexi = conjg(exi)
      ham(4*i-3, 4*i-3) = this%U*(-2d0/3d0*Sz(is)+0.5d0)
      ham(4*i-2, 4*i-2) = this%U*( 2d0/3d0*Sz(is)+0.5d0)
      ham(4*i-3, 4*i-2) = -  exi*this%U*2d0/3d0*(Sx(is)-ui*Sy(is))
      ham(4*i-2, 4*i-3) = - cexi*this%U*2d0/3d0*(Sx(is)+ui*Sy(is))
      ham(4*i-1, 4*i-1) = - this%U*(-2d0/3d0*Sz(is)+0.5d0) ! - conjg
      ham(4*i  , 4*i  ) = - this%U*( 2d0/3d0*Sz(is)+0.5d0) ! - conjg
      ham(4*i-1, 4*i  ) = conjg(- exi*this%U*2d0/3d0*(Sx(is)-ui*Sy(is))) ! + conjg
      ham(4*i  , 4*i-1) = conjg(- cexi*this%U*2d0/3d0*(Sx(is)+ui*Sy(is))) ! + conjg
    end do
    do i = 1, this%n_acc
      ham(4*i-3, 4*i-3) = ham(4*i-3, 4*i-3) - mu
      ham(4*i-2, 4*i-2) = ham(4*i-2, 4*i-2) - mu
      ham(4*i-1, 4*i-1) = ham(4*i-1, 4*i-1) + mu
      ham(4*i  , 4*i  ) = ham(4*i  , 4*i  ) + mu
    end do
    do i = 1, this%pl_1st%length()
      js = this%pl_1st%value(i)%i
      ks = this%pl_1st%value(i)%f
      j = this%ste(js)
      k = this%ste(ks)
      exi  = exp(0.5d0*ui*modulo(xi(js)-xi(ks), 2d0*pi))
      !exi  = exp(0.5d0*ui*modphase(xi(js)-xi(ks)))
      !print *,modulo(xi(js)-xi(ks), 2d0*pi), modphase(xi(js)-xi(ks))
      cexi = conjg(exi)
      ham(4*j-3, 4*k-3) = - exi*this%t(1)
      ham(4*j-2, 4*k-2) = -cexi*this%t(1)
      ham(4*j-1, 4*k-1) =  cexi*this%t(1)
      ham(4*j  , 4*k  ) =   exi*this%t(1)
      ham(4*k-3, 4*j-3) = -cexi*this%t(1)
      ham(4*k-2, 4*j-2) = - exi*this%t(1)
      ham(4*k-1, 4*j-1) =   exi*this%t(1)
      ham(4*k  , 4*j  ) =  cexi*this%t(1)
    end do
    !stop
    do i = 1, this%pl_2nd%length()
      js = this%pl_2nd%value(i)%i
      ks = this%pl_2nd%value(i)%f
      j = this%ste(js)
      k = this%ste(ks)
      exi  = exp(0.5d0*ui*modulo(xi(js)-xi(ks), 2d0*pi))
      !exi  = exp(0.5d0*ui*modphase(xi(js)-xi(ks)))
      cexi = conjg(exi)
      ham(4*j-3, 4*k-3) = - exi*this%t(2)
      ham(4*j-2, 4*k-2) = -cexi*this%t(2)
      ham(4*j-1, 4*k-1) =  cexi*this%t(2)
      ham(4*j  , 4*k  ) =   exi*this%t(2)
      ham(4*k-3, 4*j-3) = -cexi*this%t(2)
      ham(4*k-2, 4*j-2) = - exi*this%t(2)
      ham(4*k-1, 4*j-1) =   exi*this%t(2)
      ham(4*k  , 4*j  ) =  cexi*this%t(2)
    end do
    do i = 1, this%pl_h%length()
      js = this%pl_h%value(i)%i
      ks = this%pl_h%value(i)%f
      j = this%ste(js)
      k = this%ste(ks)
      exi  = exp(ui*xi(js))
      cexi = conjg(exi)
      ham(4*j-3, 4*j-3) = ham(4*j-3, 4*j-3) + 0.5d0*this%jd*Sz(ks)
      ham(4*j-2, 4*j-2) = ham(4*j-2, 4*j-2) - 0.5d0*this%jd*Sz(ks)
      ham(4*j-3, 4*j-2) = ham(4*j-3, 4*j-2) + 0.5d0*this%jd*(Sx(ks) - ui*Sy(ks))*exi
      ham(4*j-2, 4*j-3) = ham(4*j-2, 4*j-3) + 0.5d0*this%jd*(Sx(ks) + ui*Sy(ks))*cexi
      ham(4*j-1, 4*j-1) = ham(4*j-1, 4*j-1) - 0.5d0*this%jd*Sz(ks) ! -conjg
      ham(4*j  , 4*j  ) = ham(4*j  , 4*j  ) + 0.5d0*this%jd*Sz(ks) ! -conjg
      ham(4*j-1, 4*j  ) = ham(4*j-1, 4*j  ) + conjg( 0.5d0*this%jd*(Sx(ks) - ui*Sy(ks))*exi ) ! + conjg
      ham(4*j  , 4*j-1) = ham(4*j  , 4*j-1) + conjg( 0.5d0*this%jd*(Sx(ks) + ui*Sy(ks))*cexi ) ! + conjg
      
      exi  = exp(ui*xi(ks))
      cexi = conjg(exi)
      ham(4*k-3, 4*k-3) = ham(4*k-3, 4*k-3) + 0.5d0*this%jd*Sz(js)
      ham(4*k-2, 4*k-2) = ham(4*k-2, 4*k-2) - 0.5d0*this%jd*Sz(js)
      ham(4*k-3, 4*k-2) = ham(4*k-3, 4*k-2) + 0.5d0*this%jd*(Sx(js) - ui*Sy(js))*exi
      ham(4*k-2, 4*k-3) = ham(4*k-2, 4*k-3) + 0.5d0*this%jd*(Sx(js) + ui*Sy(js))*cexi
      ham(4*k-1, 4*k-1) = ham(4*k-1, 4*k-1) - 0.5d0*this%jd*Sz(js) ! -conjg
      ham(4*k  , 4*k  ) = ham(4*k  , 4*k  ) + 0.5d0*this%jd*Sz(js) ! -conjg
      ham(4*k-1, 4*k  ) = ham(4*k-1, 4*k  ) + conjg( 0.5d0*this%jd*(Sx(js) - ui*Sy(js))*exi ) ! + conjg
      ham(4*k  , 4*k-1) = ham(4*k  , 4*k-1) + conjg( 0.5d0*this%jd*(Sx(js) + ui*Sy(js))*cexi ) ! + conjg
    end do
  
    call my_zheev("l", ham, this%eg, this%wf)

  end subroutine diagonalize_bulk_xi_hamiltonian

  subroutine diagonalize_sb_2lay_hamiltonian(this, mu, delta, ne_u, ne_d, Sx, Sy, Sz)
    ! z =1 : surface
    ! z =2 : bulk
    implicit none
    class(diag), intent(inout)   :: this
    !real(8), dimension(4*this%n_acc), intent(out) :: eg
    !complex(8),dimension(4*this%n_acc, 4*this%n_acc), intent(out)   :: wf
    real(8),intent(in)    :: mu(2)
    complex(8), dimension(this%pl_1st_l(1)%length()), intent(in) :: delta
    real(8), dimension(this%n), intent(in)   :: ne_u, ne_d
    !real(8), dimension(this%n), intent(in) :: xi, chi
    real(8), dimension(this%n), intent(in) :: Sx, Sy, Sz
    complex(8),DIMENSION(4*this%n_acc, 4*this%n_acc)   :: ham
    !real(8), dimension(this%n) :: ne_u, ne_d
    integer :: i, j, k, is, js, ks
    real(8) :: t3U

    !! error cheack
    if(this%nz /= 2 .or. this%nh_l(1) /= 0) then
      print *, "nz : ", this%nz
      print *, "nh(1) : ", this%nh_l(1)
      stop "error : wrong data in diagonalize_sb_2lay_hamiltonian"
      ! this subroutine is available for only 2 layer"
    endif

    t3U = 2d0*this%t(3)**2d0/this%U

    !ne_u = 0.5d0
    !ne_d = 0.5d0

    ham = 0d0
    
    !!!!!!!!!!!!!!!!!!    surface    !!!!!!!!!!!!!!!!!!
    do i = this%n_acc_start(1), this%n_acc_end(1)
      ham(4*i-3, 4*i-3) = -mu(1)
      ham(4*i-2, 4*i-2) = -mu(1)
      ham(4*i-1, 4*i-1) =  mu(1)
      ham(4*i  , 4*i  ) =  mu(1)
    end do
    do i = 1, this%pl_1st_l(1)%length()
      js = this%pl_1st_l(1)%value(i)%i
      ks = this%pl_1st_l(1)%value(i)%f
      j = this%ste(js) ! in the 1st layer (surface), js = j, ks = k
      k = this%ste(ks)
      ham(4*j-3, 4*k-3) = -this%t(1)*(1d0-ne_d(js))*(1d0-ne_d(ks))
      ham(4*j-2, 4*k-2) = -this%t(1)*(1d0-ne_u(js))*(1d0-ne_u(ks))
      ham(4*j-1, 4*k-1) =  this%t(1)*(1d0-ne_d(js))*(1d0-ne_d(ks))
      ham(4*j  , 4*k  ) =  this%t(1)*(1d0-ne_u(js))*(1d0-ne_u(ks))
      ham(4*k-3, 4*j-3) = -this%t(1)*(1d0-ne_d(js))*(1d0-ne_d(ks))
      ham(4*k-2, 4*j-2) = -this%t(1)*(1d0-ne_u(js))*(1d0-ne_u(ks))
      ham(4*k-1, 4*j-1) =  this%t(1)*(1d0-ne_d(js))*(1d0-ne_d(ks))
      ham(4*k  , 4*j  ) =  this%t(1)*(1d0-ne_u(js))*(1d0-ne_u(ks))
    end do
    do i = 1, this%pl_2nd_l(1)%length()
      js = this%pl_2nd_l(1)%value(i)%i
      ks = this%pl_2nd_l(1)%value(i)%f
      j = this%ste(js) ! in the 1st layer (surface), js = j, ks = k
      k = this%ste(ks)
      ham(4*j-3, 4*k-3) = -this%t(2)*(1d0-ne_d(js))*(1d0-ne_d(ks))
      ham(4*j-2, 4*k-2) = -this%t(2)*(1d0-ne_u(js))*(1d0-ne_u(ks))
      ham(4*j-1, 4*k-1) =  this%t(2)*(1d0-ne_d(js))*(1d0-ne_d(ks))
      ham(4*j  , 4*k  ) =  this%t(2)*(1d0-ne_u(js))*(1d0-ne_u(ks))
      ham(4*k-3, 4*j-3) = -this%t(2)*(1d0-ne_d(js))*(1d0-ne_d(ks))
      ham(4*k-2, 4*j-2) = -this%t(2)*(1d0-ne_u(js))*(1d0-ne_u(ks))
      ham(4*k-1, 4*j-1) =  this%t(2)*(1d0-ne_d(js))*(1d0-ne_d(ks))
      ham(4*k  , 4*j  ) =  this%t(2)*(1d0-ne_u(js))*(1d0-ne_u(ks))
    end do
    do i = 1, this%pl_1st_l(1)%length()
      ! here, delta's index corresponds to pl_1st_l(1)'s index
      js = this%pl_1st_l(1)%value(i)%i
      ks = this%pl_1st_l(1)%value(i)%f
      j = this%ste(js) ! in the 1st layer ( surface ), js = j, ks = k
      k = this%ste(ks)
      ham(4*j-3, 4*k  ) = delta(i)
      ham(4*j-2, 4*k-1) = delta(i)
      ham(4*j-1, 4*k-2) = conjg(delta(i))
      ham(4*j  , 4*k-3) = conjg(delta(i))
      ham(4*k-3, 4*j  ) = delta(i)
      ham(4*k-2, 4*j-1) = delta(i)
      ham(4*k-1, 4*j-2) = conjg(delta(i))
      ham(4*k  , 4*j-3) = conjg(delta(i))
    end do
    !!!!!!!!!!!!!!    end surface    !!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!    bulk    !!!!!!!!!!!!!!!!
    do i = this%n_acc_start(2), this%n_acc_end(2)
      is = this%ets(i)
      ham(4*i-3, 4*i-3) = this%U*(-2d0/3d0*Sz(is)+0.5d0)
      ham(4*i-2, 4*i-2) = this%U*( 2d0/3d0*Sz(is)+0.5d0)
      ham(4*i-3, 4*i-2) = - this%U*2d0/3d0*(Sx(is)-ui*Sy(is))
      ham(4*i-2, 4*i-3) = - this%U*2d0/3d0*(Sx(is)+ui*Sy(is))
      ham(4*i-1, 4*i-1) = - this%U*(-2d0/3d0*Sz(is)+0.5d0) ! - conjg
      ham(4*i  , 4*i  ) = - this%U*( 2d0/3d0*Sz(is)+0.5d0) ! - conjg
      ham(4*i-1, 4*i  ) = conjg(- this%U*2d0/3d0*(Sx(is)-ui*Sy(is))) ! + conjg
      ham(4*i  , 4*i-1) = conjg(- this%U*2d0/3d0*(Sx(is)+ui*Sy(is))) ! + conjg
    end do
    do i = this%n_acc_start(2), this%n_acc_end(2)
      ham(4*i-3, 4*i-3) = ham(4*i-3, 4*i-3) - mu(2)
      ham(4*i-2, 4*i-2) = ham(4*i-2, 4*i-2) - mu(2)
      ham(4*i-1, 4*i-1) = ham(4*i-1, 4*i-1) + mu(2)
      ham(4*i  , 4*i  ) = ham(4*i  , 4*i  ) + mu(2)
    end do
    do i = 1, this%pl_1st_l(2)%length()
      js = this%pl_1st_l(2)%value(i)%i
      ks = this%pl_1st_l(2)%value(i)%f
      j = this%ste(js)
      k = this%ste(ks)
      ham(4*j-3, 4*k-3) = -this%t(1)
      ham(4*j-2, 4*k-2) = -this%t(1)
      ham(4*j-1, 4*k-1) =  this%t(1)
      ham(4*j  , 4*k  ) =  this%t(1)
      ham(4*k-3, 4*j-3) = -this%t(1)
      ham(4*k-2, 4*j-2) = -this%t(1)
      ham(4*k-1, 4*j-1) =  this%t(1)
      ham(4*k  , 4*j  ) =  this%t(1)
    end do
    do i = 1, this%pl_2nd_l(2)%length()
      js = this%pl_2nd_l(2)%value(i)%i
      ks = this%pl_2nd_l(2)%value(i)%f
      j = this%ste(js)
      k = this%ste(ks)
      ham(4*j-3, 4*k-3) = -this%t(2)
      ham(4*j-2, 4*k-2) = -this%t(2)
      ham(4*j-1, 4*k-1) =  this%t(2)
      ham(4*j  , 4*k  ) =  this%t(2)
      ham(4*k-3, 4*j-3) = -this%t(2)
      ham(4*k-2, 4*j-2) = -this%t(2)
      ham(4*k-1, 4*j-1) =  this%t(2)
      ham(4*k  , 4*j  ) =  this%t(2)
    end do
    do i = 1, this%pl_h%length()
      js = this%pl_h%value(i)%i
      ks = this%pl_h%value(i)%f
      j = this%ste(js)
      k = this%ste(ks)
      ham(4*j-3, 4*j-3) = ham(4*j-3, 4*j-3) + 0.5d0*this%jd*Sz(ks)
      ham(4*j-2, 4*j-2) = ham(4*j-2, 4*j-2) - 0.5d0*this%jd*Sz(ks)
      ham(4*j-3, 4*j-2) = ham(4*j-3, 4*j-2) + 0.5d0*this%jd*(Sx(ks) - ui*Sy(ks))
      ham(4*j-2, 4*j-3) = ham(4*j-2, 4*j-3) + 0.5d0*this%jd*(Sx(ks) + ui*Sy(ks))
      ham(4*j-1, 4*j-1) = ham(4*j-1, 4*j-1) - 0.5d0*this%jd*Sz(ks) ! -conjg
      ham(4*j  , 4*j  ) = ham(4*j  , 4*j  ) + 0.5d0*this%jd*Sz(ks) ! -conjg
      ham(4*j-1, 4*j  ) = ham(4*j-1, 4*j  ) + conjg( 0.5d0*this%jd*(Sx(ks) - ui*Sy(ks)) ) ! + conjg
      ham(4*j  , 4*j-1) = ham(4*j  , 4*j-1) + conjg( 0.5d0*this%jd*(Sx(ks) + ui*Sy(ks)) ) ! + conjg
      
      ham(4*k-3, 4*k-3) = ham(4*k-3, 4*k-3) + 0.5d0*this%jd*Sz(js)
      ham(4*k-2, 4*k-2) = ham(4*k-2, 4*k-2) - 0.5d0*this%jd*Sz(js)
      ham(4*k-3, 4*k-2) = ham(4*k-3, 4*k-2) + 0.5d0*this%jd*(Sx(js) - ui*Sy(js))
      ham(4*k-2, 4*k-3) = ham(4*k-2, 4*k-3) + 0.5d0*this%jd*(Sx(js) + ui*Sy(js))
      ham(4*k-1, 4*k-1) = ham(4*k-1, 4*k-1) - 0.5d0*this%jd*Sz(js) ! -conjg
      ham(4*k  , 4*k  ) = ham(4*k  , 4*k  ) + 0.5d0*this%jd*Sz(js) ! -conjg
      ham(4*k-1, 4*k  ) = ham(4*k-1, 4*k  ) + conjg( 0.5d0*this%jd*(Sx(js) - ui*Sy(js)) ) ! + conjg
      ham(4*k  , 4*k-1) = ham(4*k  , 4*k-1) + conjg( 0.5d0*this%jd*(Sx(js) + ui*Sy(js)) ) ! + conjg
    end do
    !!!!!!!!!!!!!!    end bulk    !!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!    interlayer    !!!!!!!!!!!!!!!!
    ! pl_z_l(i) is z-dir path list between ith layer and (i+1)th layer.
    ! here, pl_z_l(1) is interlayer path ( between 1st(surfaca) and 2nd(bulk) ).
    do i = 1, this%pl_z_l(1)%length()
      js = this%pl_z_l(1)%value(i)%i
      ks = this%pl_z_l(1)%value(i)%f
      j = this%ste(js)
      k = this%ste(ks)
      ham(4*j-3, 4*k-3) = -this%t(3)*(1d0-ne_d(js))*(1d0-ne_d(ks))
      ham(4*j-2, 4*k-2) = -this%t(3)*(1d0-ne_u(js))*(1d0-ne_u(ks))
      ham(4*j-1, 4*k-1) =  this%t(3)*(1d0-ne_d(js))*(1d0-ne_d(ks))
      ham(4*j  , 4*k  ) =  this%t(3)*(1d0-ne_u(js))*(1d0-ne_u(ks))
      ham(4*k-3, 4*j-3) = -this%t(3)*(1d0-ne_d(js))*(1d0-ne_d(ks))
      ham(4*k-2, 4*j-2) = -this%t(3)*(1d0-ne_u(js))*(1d0-ne_u(ks))
      ham(4*k-1, 4*j-1) =  this%t(3)*(1d0-ne_d(js))*(1d0-ne_d(ks))
      ham(4*k  , 4*j  ) =  this%t(3)*(1d0-ne_u(js))*(1d0-ne_u(ks))
      
      ham(4*j-3, 4*j-3) = ham(4*j-3, 4*j-3) + t3U*(Sz(ks) - 0.5d0*(ne_u(ks) + ne_d(ks)))
      ham(4*j-2, 4*j-2) = ham(4*j-2, 4*j-2) - t3U*(Sz(ks) + 0.5d0*(ne_u(ks) + ne_d(ks)))
      ham(4*j-3, 4*j-2) = ham(4*j-3, 4*j-2) + t3U*(Sx(ks) - ui*Sy(ks))
      ham(4*j-2, 4*j-3) = ham(4*j-2, 4*j-3) + t3U*(Sx(ks) + ui*Sy(ks))
      ham(4*j-1, 4*j-1) = ham(4*j-1, 4*j-1) - t3U*(Sz(ks) - 0.5d0*(ne_u(ks) + ne_d(ks))) ! -conjg
      ham(4*j  , 4*j  ) = ham(4*j  , 4*j  ) + t3U*(Sz(ks) + 0.5d0*(ne_u(ks) + ne_d(ks))) ! -conjg
      ham(4*j-1, 4*j  ) = ham(4*j-1, 4*j  ) + conjg( t3U*(Sx(ks) - ui*Sy(ks)) ) ! + conjg
      ham(4*j  , 4*j-1) = ham(4*j  , 4*j-1) + conjg( t3U*(Sx(ks) + ui*Sy(ks)) ) ! + conjg
      
      ham(4*k-3, 4*k-3) = ham(4*k-3, 4*k-3) + t3U*(Sz(js) - 0.5d0*(ne_u(js) + ne_d(js)))
      ham(4*k-2, 4*k-2) = ham(4*k-2, 4*k-2) - t3U*(Sz(js) + 0.5d0*(ne_u(js) + ne_d(js)))
      ham(4*k-3, 4*k-2) = ham(4*k-3, 4*k-2) + t3U*(Sx(js) - ui*Sy(js))
      ham(4*k-2, 4*k-3) = ham(4*k-2, 4*k-3) + t3U*(Sx(js) + ui*Sy(js))
      ham(4*k-1, 4*k-1) = ham(4*k-1, 4*k-1) - t3U*(Sz(js) - 0.5d0*(ne_u(js) + ne_d(js))) ! -conjg
      ham(4*k  , 4*k  ) = ham(4*k  , 4*k  ) + t3U*(Sz(js) + 0.5d0*(ne_u(js) + ne_d(js))) ! -conjg
      ham(4*k-1, 4*k  ) = ham(4*k-1, 4*k  ) + conjg( t3U*(Sx(js) - ui*Sy(js)) ) ! + conjg
      ham(4*k  , 4*k-1) = ham(4*k  , 4*k-1) + conjg( t3U*(Sx(js) + ui*Sy(js)) ) ! + conjg
    end do
    !!!!!!!!!!!!!!   end  interlayer    !!!!!!!!!!!!!!!!
  
    call my_zheev("l", ham, this%eg, this%wf)

  end subroutine diagonalize_sb_2lay_hamiltonian

  subroutine diagonalize_sb_2lay_xi_hamiltonian(this, mu, delta, ne_u, ne_d, Sx, Sy, Sz, xi)
    ! z =1 : surface
    ! z =2 : bulk
    implicit none
    class(diag), intent(inout)   :: this
    !real(8), dimension(4*this%n_acc), intent(out) :: eg
    !complex(8),dimension(4*this%n_acc, 4*this%n_acc), intent(out)   :: wf
    real(8),intent(in)    :: mu(2)
    complex(8), dimension(this%pl_1st_l(1)%length()), intent(in) :: delta
    real(8), dimension(this%n), intent(in)   :: ne_u, ne_d
    !real(8), dimension(this%n), intent(in) :: xi, chi
    real(8), dimension(this%n), intent(in) :: Sx, Sy, Sz, xi
    complex(8),DIMENSION(4*this%n_acc, 4*this%n_acc)   :: ham
    !real(8), dimension(this%n) :: ne_u, ne_d
    integer :: i, j, k, is, js, ks
    real(8) :: t3U
    complex(8) :: exi,cexi

    !! error cheack
    if(this%nz /= 2 .or. this%nh_l(1) /= 0) then
      print *, "nz : ", this%nz
      print *, "nh(1) : ", this%nh_l(1)
      stop "error : wrong data in diagonalize_sb_2lay_xi_hamiltonian"
      ! this subroutine is available for only 2 layer"
    endif

    t3U = 2*this%t(3)**2/this%U

    !ne_u = 0.5d0
    !ne_d = 0.5d0

    ham = 0d0
    
    !!!!!!!!!!!!!!!!!!    surface    !!!!!!!!!!!!!!!!!!
    do i = this%n_acc_start(1), this%n_acc_end(1)
      ham(4*i-3, 4*i-3) = -mu(1)
      ham(4*i-2, 4*i-2) = -mu(1)
      ham(4*i-1, 4*i-1) =  mu(1)
      ham(4*i  , 4*i  ) =  mu(1)
    end do
    do i = 1, this%pl_1st_l(1)%length()
      js = this%pl_1st_l(1)%value(i)%f
      ks = this%pl_1st_l(1)%value(i)%i
      j = this%ste(js) ! in the 1st layer (surface), js = j, ks = k
      k = this%ste(ks)
      !exi  = exp(0.5d0*ui*principal(xi(js)-xi(ks)))
      exi  = exp(0.5d0*ui*modulo(xi(js)-xi(ks), 2d0*pi))
      !exi  = exp(-0.5d0*ui*modulo(xi(js)-xi(ks), 2d0*pi))
      !exi  = exp(0.5d0*ui*(xi(js)-xi(ks)))
      !exi  = exp(-0.5d0*ui*(xi(js)-xi(ks)))
      !exi  = 1d0
      cexi = conjg(exi)
      ham(4*j-3, 4*k-3) = - exi*this%t(1)*(1d0-ne_d(js))*(1d0-ne_d(ks))
      ham(4*j-2, 4*k-2) = -cexi*this%t(1)*(1d0-ne_u(js))*(1d0-ne_u(ks))
      ham(4*j-1, 4*k-1) =  cexi*this%t(1)*(1d0-ne_d(js))*(1d0-ne_d(ks))
      ham(4*j  , 4*k  ) =   exi*this%t(1)*(1d0-ne_u(js))*(1d0-ne_u(ks))
      ham(4*k-3, 4*j-3) = -cexi*this%t(1)*(1d0-ne_d(js))*(1d0-ne_d(ks))
      ham(4*k-2, 4*j-2) = - exi*this%t(1)*(1d0-ne_u(js))*(1d0-ne_u(ks))
      ham(4*k-1, 4*j-1) =   exi*this%t(1)*(1d0-ne_d(js))*(1d0-ne_d(ks))
      ham(4*k  , 4*j  ) =  cexi*this%t(1)*(1d0-ne_u(js))*(1d0-ne_u(ks))
    end do
    do i = 1, this%pl_2nd_l(1)%length()
      js = this%pl_2nd_l(1)%value(i)%f
      ks = this%pl_2nd_l(1)%value(i)%i
      j = this%ste(js) ! in the 1st layer (surface), js = j, ks = k
      k = this%ste(ks)
      !exi  = exp(0.5d0*ui*principal(xi(js)-xi(ks)))
      exi  = exp(0.5d0*ui*modulo(xi(js)-xi(ks), 2d0*pi))
      !exi  = exp(-0.5d0*ui*modulo(xi(js)-xi(ks), 2d0*pi))
      !exi  = exp(0.5d0*ui*(xi(js)-xi(ks)))
      !exi  = exp(-0.5d0*ui*(xi(js)-xi(ks)))
      !exi  = 1d0
      cexi = conjg(exi)
      ham(4*j-3, 4*k-3) = - exi*this%t(2)*(1d0-ne_d(js))*(1d0-ne_d(ks))
      ham(4*j-2, 4*k-2) = -cexi*this%t(2)*(1d0-ne_u(js))*(1d0-ne_u(ks))
      ham(4*j-1, 4*k-1) =  cexi*this%t(2)*(1d0-ne_d(js))*(1d0-ne_d(ks))
      ham(4*j  , 4*k  ) =   exi*this%t(2)*(1d0-ne_u(js))*(1d0-ne_u(ks))
      ham(4*k-3, 4*j-3) = -cexi*this%t(2)*(1d0-ne_d(js))*(1d0-ne_d(ks))
      ham(4*k-2, 4*j-2) = - exi*this%t(2)*(1d0-ne_u(js))*(1d0-ne_u(ks))
      ham(4*k-1, 4*j-1) =   exi*this%t(2)*(1d0-ne_d(js))*(1d0-ne_d(ks))
      ham(4*k  , 4*j  ) =  cexi*this%t(2)*(1d0-ne_u(js))*(1d0-ne_u(ks))
    end do
    do i = 1, this%pl_1st_l(1)%length()
      ! here, delta's index corresponds to pl_1st_l(1)'s index
      js = this%pl_1st_l(1)%value(i)%f
      ks = this%pl_1st_l(1)%value(i)%i
      j = this%ste(js) ! in the 1st layer ( surface ), js = j, ks = k
      k = this%ste(ks)
      !exi  = exp(0.5d0*ui*principal(xi(js)-xi(ks)))
      exi  = exp(0.5d0*ui*modulo(xi(js)-xi(ks), 2d0*pi))
      !exi  = exp(-0.5d0*ui*modulo(xi(js)-xi(ks), 2d0*pi))
      !exi  = exp(0.5d0*ui*(xi(js)-xi(ks)))
      !exi  = exp(-0.5d0*ui*(xi(js)-xi(ks)))
      !exi  = exp(-0.5d0*ui*modulo(xi(js)-xi(ks), 2d0*pi))
      !exi  = 1d0
      cexi = conjg(exi)
     
      !ham(4*j-3, 4*k  ) = exi*delta(i)
      !ham(4*j-2, 4*k-1) = cexi*delta(i)
      !ham(4*j-1, 4*k-2) = exi*conjg(delta(i))
      !ham(4*j  , 4*k-3) = cexi*conjg(delta(i))
      !ham(4*k-3, 4*j  ) = cexi*delta(i)
      !ham(4*k-2, 4*j-1) = exi*delta(i)
      !ham(4*k-1, 4*j-2) = cexi*conjg(delta(i))
      !ham(4*k  , 4*j-3) = exi*conjg(delta(i))
      
      ham(4*j-3, 4*k  ) = exi*delta(i)
      ham(4*j-2, 4*k-1) = cexi*delta(i)
      ham(4*j-1, 4*k-2) = cexi*conjg(delta(i))
      ham(4*j  , 4*k-3) = exi*conjg(delta(i))
      ham(4*k-3, 4*j  ) = cexi*delta(i)
      ham(4*k-2, 4*j-1) = exi*delta(i)
      ham(4*k-1, 4*j-2) = exi*conjg(delta(i))
      ham(4*k  , 4*j-3) = cexi*conjg(delta(i))
    end do
    !!!!!!!!!!!!!!    end surface    !!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!    bulk    !!!!!!!!!!!!!!!!
    do i = this%n_acc_start(2), this%n_acc_end(2)
      is = this%ets(i)
      exi  = exp(ui*xi(is))
      !exi  = exp(-ui*xi(is))
      !exi=1d0
      cexi = conjg(exi)
      ham(4*i-3, 4*i-3) = this%U*(-2d0/3d0*Sz(is)+0.5d0)
      ham(4*i-2, 4*i-2) = this%U*( 2d0/3d0*Sz(is)+0.5d0)
      ham(4*i-3, 4*i-2) = -  exi*this%U*2d0/3d0*(Sx(is)-ui*Sy(is))
      ham(4*i-2, 4*i-3) = - cexi*this%U*2d0/3d0*(Sx(is)+ui*Sy(is))
      ham(4*i-1, 4*i-1) = - this%U*(-2d0/3d0*Sz(is)+0.5d0) ! - conjg
      ham(4*i  , 4*i  ) = - this%U*( 2d0/3d0*Sz(is)+0.5d0) ! - conjg
      ham(4*i-1, 4*i  ) = conjg(-  exi*this%U*2d0/3d0*(Sx(is)-ui*Sy(is))) ! + conjg
      ham(4*i  , 4*i-1) = conjg(- cexi*this%U*2d0/3d0*(Sx(is)+ui*Sy(is))) ! + conjg
    end do
    do i = this%n_acc_start(2), this%n_acc_end(2)
      ham(4*i-3, 4*i-3) = ham(4*i-3, 4*i-3) - mu(2)
      ham(4*i-2, 4*i-2) = ham(4*i-2, 4*i-2) - mu(2)
      ham(4*i-1, 4*i-1) = ham(4*i-1, 4*i-1) + mu(2)
      ham(4*i  , 4*i  ) = ham(4*i  , 4*i  ) + mu(2)
    end do
    do i = 1, this%pl_1st_l(2)%length()
      js = this%pl_1st_l(2)%value(i)%f
      ks = this%pl_1st_l(2)%value(i)%i
      j = this%ste(js)
      k = this%ste(ks)
      !exi  = exp(0.5d0*ui*principal(xi(js)-xi(ks)))
      exi  = exp(0.5d0*ui*modulo(xi(js)-xi(ks), 2d0*pi))
      !exi  = exp(-0.5d0*ui*modulo(xi(js)-xi(ks), 2d0*pi))
      !exi  = exp(0.5d0*ui*(xi(js)-xi(ks)))
      !exi  = exp(-0.5d0*ui*(xi(js)-xi(ks)))
      !exi  = exp(-0.5d0*ui*modulo(xi(js)-xi(ks), 2d0*pi))
      cexi = conjg(exi)
      ham(4*j-3, 4*k-3) = - exi*this%t(1)
      ham(4*j-2, 4*k-2) = -cexi*this%t(1)
      ham(4*j-1, 4*k-1) =  cexi*this%t(1)
      ham(4*j  , 4*k  ) =   exi*this%t(1)
      ham(4*k-3, 4*j-3) = -cexi*this%t(1)
      ham(4*k-2, 4*j-2) = - exi*this%t(1)
      ham(4*k-1, 4*j-1) =   exi*this%t(1)
      ham(4*k  , 4*j  ) =  cexi*this%t(1)
    end do
    do i = 1, this%pl_2nd_l(2)%length()
      js = this%pl_2nd_l(2)%value(i)%f
      ks = this%pl_2nd_l(2)%value(i)%i
      j = this%ste(js)
      k = this%ste(ks)
      !exi  = exp(0.5d0*ui*principal(xi(js)-xi(ks)))
      exi  = exp(0.5d0*ui*modulo(xi(js)-xi(ks), 2d0*pi))
      !exi  = exp(-0.5d0*ui*modulo(xi(js)-xi(ks), 2d0*pi))
      !exi  = exp(0.5d0*ui*(xi(js)-xi(ks)))
      !exi  = exp(-0.5d0*ui*(xi(js)-xi(ks)))
      !exi  = exp(-0.5d0*ui*modulo(xi(js)-xi(ks), 2d0*pi))
      cexi = conjg(exi)
      ham(4*j-3, 4*k-3) = - exi*this%t(2)
      ham(4*j-2, 4*k-2) = -cexi*this%t(2)
      ham(4*j-1, 4*k-1) =  cexi*this%t(2)
      ham(4*j  , 4*k  ) =   exi*this%t(2)
      ham(4*k-3, 4*j-3) = -cexi*this%t(2)
      ham(4*k-2, 4*j-2) = - exi*this%t(2)
      ham(4*k-1, 4*j-1) =   exi*this%t(2)
      ham(4*k  , 4*j  ) =  cexi*this%t(2)
    end do
    do i = 1, this%pl_h%length()
      js = this%pl_h%value(i)%f
      ks = this%pl_h%value(i)%i
      j = this%ste(js)
      k = this%ste(ks)
      exi  = exp(ui*xi(js))
      !exi  = exp(-ui*xi(js))
      !exi  = exp(ui*xi(ks))
      cexi = conjg(exi)
      ham(4*j-3, 4*j-3) = ham(4*j-3, 4*j-3) + 0.5d0*this%jd*Sz(ks)
      ham(4*j-2, 4*j-2) = ham(4*j-2, 4*j-2) - 0.5d0*this%jd*Sz(ks)
      ham(4*j-3, 4*j-2) = ham(4*j-3, 4*j-2) + 0.5d0*this%jd*(Sx(ks) - ui*Sy(ks))*exi
      ham(4*j-2, 4*j-3) = ham(4*j-2, 4*j-3) + 0.5d0*this%jd*(Sx(ks) + ui*Sy(ks))*cexi
      ham(4*j-1, 4*j-1) = ham(4*j-1, 4*j-1) - 0.5d0*this%jd*Sz(ks) ! -conjg
      ham(4*j  , 4*j  ) = ham(4*j  , 4*j  ) + 0.5d0*this%jd*Sz(ks) ! -conjg
      ham(4*j-1, 4*j  ) = ham(4*j-1, 4*j  ) + conjg( 0.5d0*this%jd*(Sx(ks) - ui*Sy(ks))*exi ) ! + conjg
      ham(4*j  , 4*j-1) = ham(4*j  , 4*j-1) + conjg( 0.5d0*this%jd*(Sx(ks) + ui*Sy(ks))*cexi ) ! + conjg
      
      exi  = exp(ui*xi(ks))
      !exi  = exp(-ui*xi(ks))
      !exi  = exp(ui*xi(js))
      cexi = conjg(exi)
      ham(4*k-3, 4*k-3) = ham(4*k-3, 4*k-3) + 0.5d0*this%jd*Sz(js)
      ham(4*k-2, 4*k-2) = ham(4*k-2, 4*k-2) - 0.5d0*this%jd*Sz(js)
      ham(4*k-3, 4*k-2) = ham(4*k-3, 4*k-2) + 0.5d0*this%jd*(Sx(js) - ui*Sy(js))*exi
      ham(4*k-2, 4*k-3) = ham(4*k-2, 4*k-3) + 0.5d0*this%jd*(Sx(js) + ui*Sy(js))*cexi
      ham(4*k-1, 4*k-1) = ham(4*k-1, 4*k-1) - 0.5d0*this%jd*Sz(js) ! -conjg
      ham(4*k  , 4*k  ) = ham(4*k  , 4*k  ) + 0.5d0*this%jd*Sz(js) ! -conjg
      ham(4*k-1, 4*k  ) = ham(4*k-1, 4*k  ) + conjg( 0.5d0*this%jd*(Sx(js) - ui*Sy(js))*exi ) ! + conjg
      ham(4*k  , 4*k-1) = ham(4*k  , 4*k-1) + conjg( 0.5d0*this%jd*(Sx(js) + ui*Sy(js))*cexi ) ! + conjg
    end do
    !!!!!!!!!!!!!!    end bulk    !!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!    interlayer    !!!!!!!!!!!!!!!!
    ! pl_z_l(i) is z-dir path list between ith layer and (i+1)th layer.
    ! here, pl_z_l(1) is interlayer path ( between 1st(surfaca) and 2nd(bulk) ).
    do i = 1, this%pl_z_l(1)%length()
      js = this%pl_z_l(1)%value(i)%f
      ks = this%pl_z_l(1)%value(i)%i
      j = this%ste(js)
      k = this%ste(ks)
      !exi  = exp(0.5d0*ui*principal(xi(js)-xi(ks)))
      exi  = exp(0.5d0*ui*modulo(xi(js)-xi(ks), 2d0*pi))
      !exi  = exp(-0.5d0*ui*modulo(xi(js)-xi(ks), 2d0*pi))
      !exi  = exp(0.5d0*ui*(xi(js)-xi(ks)))
      !exi  = exp(-0.5d0*ui*(xi(js)-xi(ks)))
      !exi  = exp(-0.5d0*ui*xi(ks))
      !exi  = exp(0.5d0*ui*xi(ks))
      cexi = conjg(exi)
      ham(4*j-3, 4*k-3) = - exi*this%t(3)*(1d0-ne_d(js))*(1d0-ne_d(ks))
      ham(4*j-2, 4*k-2) = -cexi*this%t(3)*(1d0-ne_u(js))*(1d0-ne_u(ks))
      ham(4*j-1, 4*k-1) =  cexi*this%t(3)*(1d0-ne_d(js))*(1d0-ne_d(ks))
      ham(4*j  , 4*k  ) =   exi*this%t(3)*(1d0-ne_u(js))*(1d0-ne_u(ks))
      ham(4*k-3, 4*j-3) = -cexi*this%t(3)*(1d0-ne_d(js))*(1d0-ne_d(ks))
      ham(4*k-2, 4*j-2) = - exi*this%t(3)*(1d0-ne_u(js))*(1d0-ne_u(ks))
      ham(4*k-1, 4*j-1) =   exi*this%t(3)*(1d0-ne_d(js))*(1d0-ne_d(ks))
      ham(4*k  , 4*j  ) =  cexi*this%t(3)*(1d0-ne_u(js))*(1d0-ne_u(ks))
      
      exi  = exp(ui*xi(js))
      !exi  = exp(-ui*xi(js))
      !exi  = 1d0
      cexi = conjg(exi)
      ham(4*j-3, 4*j-3) = ham(4*j-3, 4*j-3) + t3U*(Sz(ks) - 0.5d0*(ne_u(ks) + ne_d(ks)))
      ham(4*j-2, 4*j-2) = ham(4*j-2, 4*j-2) - t3U*(Sz(ks) + 0.5d0*(ne_u(ks) + ne_d(ks)))
      ham(4*j-3, 4*j-2) = ham(4*j-3, 4*j-2) + t3U*(Sx(ks) - ui*Sy(ks))*exi
      ham(4*j-2, 4*j-3) = ham(4*j-2, 4*j-3) + t3U*(Sx(ks) + ui*Sy(ks))*cexi
      ham(4*j-1, 4*j-1) = ham(4*j-1, 4*j-1) - t3U*(Sz(ks) - 0.5d0*(ne_u(ks) + ne_d(ks))) ! -conjg
      ham(4*j  , 4*j  ) = ham(4*j  , 4*j  ) + t3U*(Sz(ks) + 0.5d0*(ne_u(ks) + ne_d(ks))) ! -conjg
      ham(4*j-1, 4*j  ) = ham(4*j-1, 4*j  ) + conjg( t3U*(Sx(ks) - ui*Sy(ks))*exi ) ! + conjg
      ham(4*j  , 4*j-1) = ham(4*j  , 4*j-1) + conjg( t3U*(Sx(ks) + ui*Sy(ks))*cexi ) ! + conjg
      
      exi  = exp(ui*xi(ks))
      !exi  = exp(-ui*xi(ks))
      cexi = conjg(exi)
      ham(4*k-3, 4*k-3) = ham(4*k-3, 4*k-3) + t3U*(Sz(js) - 0.5d0*(ne_u(js) + ne_d(js)))
      ham(4*k-2, 4*k-2) = ham(4*k-2, 4*k-2) - t3U*(Sz(js) + 0.5d0*(ne_u(js) + ne_d(js)))
      ham(4*k-3, 4*k-2) = ham(4*k-3, 4*k-2) + t3U*(Sx(js) - ui*Sy(js))*exi
      ham(4*k-2, 4*k-3) = ham(4*k-2, 4*k-3) + t3U*(Sx(js) + ui*Sy(js))*cexi
      ham(4*k-1, 4*k-1) = ham(4*k-1, 4*k-1) - t3U*(Sz(js) - 0.5d0*(ne_u(js) + ne_d(js))) ! -conjg
      ham(4*k  , 4*k  ) = ham(4*k  , 4*k  ) + t3U*(Sz(js) + 0.5d0*(ne_u(js) + ne_d(js))) ! -conjg
      ham(4*k-1, 4*k  ) = ham(4*k-1, 4*k  ) + conjg( t3U*(Sx(js) - ui*Sy(js))*exi ) ! + conjg
      ham(4*k  , 4*k-1) = ham(4*k  , 4*k-1) + conjg( t3U*(Sx(js) + ui*Sy(js))*cexi ) ! + conjg
    end do
    !!!!!!!!!!!!!!   end  interlayer    !!!!!!!!!!!!!!!!
  
    call my_zheev("l", ham, this%eg, this%wf)

  end subroutine diagonalize_sb_2lay_xi_hamiltonian
  
!  subroutine diagonalize_sb_2lay_xichi_hamiltonian(this, mu, delta, ne_u, ne_d, Sx, Sy, Sz, xi, chi)
!    ! z =1 : surface
!    ! z =2 : bulk
!    implicit none
!    class(diag), intent(inout)   :: this
!    !real(8), dimension(4*this%n_acc), intent(out) :: eg
!    !complex(8),dimension(4*this%n_acc, 4*this%n_acc), intent(out)   :: wf
!    real(8),intent(in)    :: mu(2)
!    complex(8), dimension(this%pl_1st_l(1)%length()), intent(in) :: delta
!    real(8), dimension(this%n), intent(in)   :: ne_u, ne_d
!    !real(8), dimension(this%n), intent(in) :: xi, chi
!    real(8), dimension(this%n), intent(in) :: Sx, Sy, Sz, xi, chi
!    complex(8),DIMENSION(4*this%n_acc, 4*this%n_acc)   :: ham
!    !real(8), dimension(this%n) :: ne_u, ne_d
!    integer :: i, j, k, is, js, ks
!    real(8) :: t3U
!    complex(8) :: exi, cexi, echi, cechi
!
!    !! error cheack
!    if(this%nz /= 2 .or. this%nh_l(1) /= 0) then
!      print *, "nz : ", this%nz
!      print *, "nh(1) : ", this%nh_l(1)
!      stop "error : wrong data in diagonalize_sb_2lay_xichi_hamiltonian"
!      ! this subroutine is available for only 2 layer"
!    endif
!
!    t3U = 2*this%t(3)**2/this%U
!
!    !ne_u = 0.5d0
!    !ne_d = 0.5d0
!
!    ham = 0d0
!    
!    !!!!!!!!!!!!!!!!!!    surface    !!!!!!!!!!!!!!!!!!
!    do i = this%n_acc_start(1), this%n_acc_end(1)
!      ham(4*i-3, 4*i-3) = -mu(1)
!      ham(4*i-2, 4*i-2) = -mu(1)
!      ham(4*i-1, 4*i-1) =  mu(1)
!      ham(4*i  , 4*i  ) =  mu(1)
!    end do
!    do i = 1, this%pl_1st_l(1)%length()
!      js = this%pl_1st_l(1)%value(i)%i
!      ks = this%pl_1st_l(1)%value(i)%f
!      j = this%ste(js) ! in the 1st layer (surface), js = j, ks = k
!      k = this%ste(ks)
!      !exi  = exp(0.5d0*ui*modulo(xi(js)-xi(ks), 2d0*pi))
!      !exi  = exp(-0.5d0*ui*modulo(xi(js)-xi(ks), 2d0*pi))
!      exi  = exp(-0.5d0*ui*(xi(js)-xi(ks)))
!      !exi  = exp(0.5d0*ui*(xi(js)-xi(ks)))
!      !exi  = exp(-0.5d0*ui*(xi(js)-xi(ks)))
!      !exi  = 1d0
!      cexi = conjg(exi)
!      echi  = exp(0.5d0*ui*(chi(js)-chi(ks)))
!      cechi = conjg(echi)
!      ham(4*j-3, 4*k-3) = -  echi* exi*this%t(1)*(1d0-ne_d(js))*(1d0-ne_d(ks))
!      ham(4*j-2, 4*k-2) = -  echi*cexi*this%t(1)*(1d0-ne_u(js))*(1d0-ne_u(ks))
!      ham(4*j-1, 4*k-1) =   cechi*cexi*this%t(1)*(1d0-ne_d(js))*(1d0-ne_d(ks))
!      ham(4*j  , 4*k  ) =   cechi* exi*this%t(1)*(1d0-ne_u(js))*(1d0-ne_u(ks))
!      ham(4*k-3, 4*j-3) = - cechi*cexi*this%t(1)*(1d0-ne_d(js))*(1d0-ne_d(ks))
!      ham(4*k-2, 4*j-2) = - cechi* exi*this%t(1)*(1d0-ne_u(js))*(1d0-ne_u(ks))
!      ham(4*k-1, 4*j-1) =    echi* exi*this%t(1)*(1d0-ne_d(js))*(1d0-ne_d(ks))
!      ham(4*k  , 4*j  ) =    echi*cexi*this%t(1)*(1d0-ne_u(js))*(1d0-ne_u(ks))
!    end do
!    do i = 1, this%pl_2nd_l(1)%length()
!      js = this%pl_2nd_l(1)%value(i)%i
!      ks = this%pl_2nd_l(1)%value(i)%f
!      j = this%ste(js) ! in the 1st layer (surface), js = j, ks = k
!      k = this%ste(ks)
!      !exi  = exp(0.5d0*ui*modulo(xi(js)-xi(ks), 2d0*pi))
!      !exi  = exp(-0.5d0*ui*modulo(xi(js)-xi(ks), 2d0*pi))
!      exi  = exp(-0.5d0*ui*(xi(js)-xi(ks)))
!      !exi  = exp(0.5d0*ui*(xi(js)-xi(ks)))
!      !exi  = exp(-0.5d0*ui*(xi(js)-xi(ks)))
!      !exi  = 1d0
!      cexi = conjg(exi)
!      echi  = exp(0.5d0*ui*(chi(js)-chi(ks)))
!      cechi = conjg(echi)
!      ham(4*j-3, 4*k-3) = -  echi* exi*this%t(2)*(1d0-ne_d(js))*(1d0-ne_d(ks))
!      ham(4*j-2, 4*k-2) = -  echi*cexi*this%t(2)*(1d0-ne_u(js))*(1d0-ne_u(ks))
!      ham(4*j-1, 4*k-1) =   cechi*cexi*this%t(2)*(1d0-ne_d(js))*(1d0-ne_d(ks))
!      ham(4*j  , 4*k  ) =   cechi* exi*this%t(2)*(1d0-ne_u(js))*(1d0-ne_u(ks))
!      ham(4*k-3, 4*j-3) = - cechi*cexi*this%t(2)*(1d0-ne_d(js))*(1d0-ne_d(ks))
!      ham(4*k-2, 4*j-2) = - cechi* exi*this%t(2)*(1d0-ne_u(js))*(1d0-ne_u(ks))
!      ham(4*k-1, 4*j-1) =    echi* exi*this%t(2)*(1d0-ne_d(js))*(1d0-ne_d(ks))
!      ham(4*k  , 4*j  ) =    echi*cexi*this%t(2)*(1d0-ne_u(js))*(1d0-ne_u(ks))
!    end do
!    do i = 1, this%pl_1st_l(1)%length()
!      ! here, delta's index corresponds to pl_1st_l(1)'s index
!      js = this%pl_1st_l(1)%value(i)%i
!      ks = this%pl_1st_l(1)%value(i)%f
!      j = this%ste(js) ! in the 1st layer ( surface ), js = j, ks = k
!      k = this%ste(ks)
!      !exi  = exp(0.5d0*ui*modulo(xi(js)-xi(ks), 2d0*pi))
!      !exi  = exp(-0.5d0*ui*modulo(xi(js)-xi(ks), 2d0*pi))
!      !exi  = exp(0.5d0*ui*(xi(js)-xi(ks)))
!      exi  = exp(-0.5d0*ui*(xi(js)-xi(ks)))
!      !exi  = exp(-0.5d0*ui*modulo(xi(js)-xi(ks), 2d0*pi))
!      !exi  = 1d0
!      cexi = conjg(exi)
!      echi  = exp(0.5d0*ui*(chi(js)+chi(ks)) )! + chi_j + chi_k
!      !echi  = 1d0
!      cechi = conjg(echi)
!      ham(4*j-3, 4*k  ) = echi*exi*delta(i)
!      ham(4*j-2, 4*k-1) = echi*cexi*delta(i)
!      ham(4*j-1, 4*k-2) = cechi*cexi*conjg(delta(i))
!      ham(4*j  , 4*k-3) = cechi*exi*conjg(delta(i))
!      ham(4*k-3, 4*j  ) = cechi*cexi*delta(i)
!      ham(4*k-2, 4*j-1) = cechi*exi*delta(i)
!      ham(4*k-1, 4*j-2) = echi*exi*conjg(delta(i))
!      ham(4*k  , 4*j-3) = echi*cexi*conjg(delta(i))
!    end do
!    !!!!!!!!!!!!!!    end surface    !!!!!!!!!!!!!!!!!!
!
!    !!!!!!!!!!!!!!!!!!    bulk    !!!!!!!!!!!!!!!!
!    do i = this%n_acc_start(2), this%n_acc_end(2)
!      is = this%ets(i)
!      !exi  = exp(ui*xi(is))
!      exi  = exp(-ui*xi(is))
!      !exi=1d0
!      cexi = conjg(exi)
!      ham(4*i-3, 4*i-3) = this%U*(-2d0/3d0*Sz(is)+0.5d0)
!      ham(4*i-2, 4*i-2) = this%U*( 2d0/3d0*Sz(is)+0.5d0)
!      ham(4*i-3, 4*i-2) = -  exi*this%U*2d0/3d0*(Sx(is)-ui*Sy(is))
!      ham(4*i-2, 4*i-3) = - cexi*this%U*2d0/3d0*(Sx(is)+ui*Sy(is))
!      ham(4*i-1, 4*i-1) = - this%U*(-2d0/3d0*Sz(is)+0.5d0) ! - conjg
!      ham(4*i  , 4*i  ) = - this%U*( 2d0/3d0*Sz(is)+0.5d0) ! - conjg
!      ham(4*i-1, 4*i  ) = conjg(-  exi*this%U*2d0/3d0*(Sx(is)-ui*Sy(is))) ! + conjg
!      ham(4*i  , 4*i-1) = conjg(- cexi*this%U*2d0/3d0*(Sx(is)+ui*Sy(is))) ! + conjg
!    end do
!    do i = this%n_acc_start(2), this%n_acc_end(2)
!      ham(4*i-3, 4*i-3) = ham(4*i-3, 4*i-3) - mu(2)
!      ham(4*i-2, 4*i-2) = ham(4*i-2, 4*i-2) - mu(2)
!      ham(4*i-1, 4*i-1) = ham(4*i-1, 4*i-1) + mu(2)
!      ham(4*i  , 4*i  ) = ham(4*i  , 4*i  ) + mu(2)
!    end do
!    do i = 1, this%pl_1st_l(2)%length()
!      js = this%pl_1st_l(2)%value(i)%i
!      ks = this%pl_1st_l(2)%value(i)%f
!      j = this%ste(js)
!      k = this%ste(ks)
!      !exi  = exp(0.5d0*ui*modulo(xi(js)-xi(ks), 2d0*pi))
!      !exi  = exp(-0.5d0*ui*modulo(xi(js)-xi(ks), 2d0*pi))
!      !exi  = exp(0.5d0*ui*(xi(js)-xi(ks)))
!      exi  = exp(-0.5d0*ui*(xi(js)-xi(ks)))
!      !exi  = exp(-0.5d0*ui*modulo(xi(js)-xi(ks), 2d0*pi))
!      cexi = conjg(exi)
!      echi  = exp(0.5d0*ui*(chi(js)-chi(ks)))
!      cechi = conjg(echi)
!      ham(4*j-3, 4*k-3) = -  echi* exi*this%t(1)
!      ham(4*j-2, 4*k-2) = -  echi*cexi*this%t(1)
!      ham(4*j-1, 4*k-1) =   cechi*cexi*this%t(1)
!      ham(4*j  , 4*k  ) =   cechi* exi*this%t(1)
!      ham(4*k-3, 4*j-3) = - cechi*cexi*this%t(1)
!      ham(4*k-2, 4*j-2) = - cechi* exi*this%t(1)
!      ham(4*k-1, 4*j-1) =    echi* exi*this%t(1)
!      ham(4*k  , 4*j  ) =    echi*cexi*this%t(1)
!    end do
!    do i = 1, this%pl_2nd_l(2)%length()
!      js = this%pl_2nd_l(2)%value(i)%i
!      ks = this%pl_2nd_l(2)%value(i)%f
!      j = this%ste(js)
!      k = this%ste(ks)
!      !exi  = exp(0.5d0*ui*modulo(xi(js)-xi(ks), 2d0*pi))
!      !exi  = exp(-0.5d0*ui*modulo(xi(js)-xi(ks), 2d0*pi))
!      !exi  = exp(0.5d0*ui*(xi(js)-xi(ks)))
!      exi  = exp(-0.5d0*ui*(xi(js)-xi(ks)))
!      !exi  = exp(-0.5d0*ui*modulo(xi(js)-xi(ks), 2d0*pi))
!      cexi = conjg(exi)
!      echi  = exp(0.5d0*ui*(chi(js)-chi(ks)))
!      cechi = conjg(echi)
!      ham(4*j-3, 4*k-3) = -  echi* exi*this%t(2)
!      ham(4*j-2, 4*k-2) = -  echi*cexi*this%t(2)
!      ham(4*j-1, 4*k-1) =   cechi*cexi*this%t(2)
!      ham(4*j  , 4*k  ) =   cechi* exi*this%t(2)
!      ham(4*k-3, 4*j-3) = - cechi*cexi*this%t(2)
!      ham(4*k-2, 4*j-2) = - cechi* exi*this%t(2)
!      ham(4*k-1, 4*j-1) =    echi* exi*this%t(2)
!      ham(4*k  , 4*j  ) =    echi*cexi*this%t(2)
!    end do
!    do i = 1, this%pl_h%length()
!      js = this%pl_h%value(i)%i
!      ks = this%pl_h%value(i)%f
!      j = this%ste(js)
!      k = this%ste(ks)
!      !exi  = exp(ui*xi(js))
!      exi  = exp(-ui*xi(js))
!      !exi  = exp(ui*xi(ks))
!      cexi = conjg(exi)
!      ham(4*j-3, 4*j-3) = ham(4*j-3, 4*j-3) + 0.5d0*this%jd*Sz(ks)
!      ham(4*j-2, 4*j-2) = ham(4*j-2, 4*j-2) - 0.5d0*this%jd*Sz(ks)
!      ham(4*j-3, 4*j-2) = ham(4*j-3, 4*j-2) + 0.5d0*this%jd*(Sx(ks) - ui*Sy(ks))*exi
!      ham(4*j-2, 4*j-3) = ham(4*j-2, 4*j-3) + 0.5d0*this%jd*(Sx(ks) + ui*Sy(ks))*cexi
!      ham(4*j-1, 4*j-1) = ham(4*j-1, 4*j-1) - 0.5d0*this%jd*Sz(ks) ! -conjg
!      ham(4*j  , 4*j  ) = ham(4*j  , 4*j  ) + 0.5d0*this%jd*Sz(ks) ! -conjg
!      ham(4*j-1, 4*j  ) = ham(4*j-1, 4*j  ) + conjg( 0.5d0*this%jd*(Sx(ks) - ui*Sy(ks))*exi ) ! + conjg
!      ham(4*j  , 4*j-1) = ham(4*j  , 4*j-1) + conjg( 0.5d0*this%jd*(Sx(ks) + ui*Sy(ks))*cexi ) ! + conjg
!      
!      !exi  = exp(ui*xi(ks))
!      exi  = exp(-ui*xi(ks))
!      !exi  = exp(ui*xi(js))
!      cexi = conjg(exi)
!      ham(4*k-3, 4*k-3) = ham(4*k-3, 4*k-3) + 0.5d0*this%jd*Sz(js)
!      ham(4*k-2, 4*k-2) = ham(4*k-2, 4*k-2) - 0.5d0*this%jd*Sz(js)
!      ham(4*k-3, 4*k-2) = ham(4*k-3, 4*k-2) + 0.5d0*this%jd*(Sx(js) - ui*Sy(js))*exi
!      ham(4*k-2, 4*k-3) = ham(4*k-2, 4*k-3) + 0.5d0*this%jd*(Sx(js) + ui*Sy(js))*cexi
!      ham(4*k-1, 4*k-1) = ham(4*k-1, 4*k-1) - 0.5d0*this%jd*Sz(js) ! -conjg
!      ham(4*k  , 4*k  ) = ham(4*k  , 4*k  ) + 0.5d0*this%jd*Sz(js) ! -conjg
!      ham(4*k-1, 4*k  ) = ham(4*k-1, 4*k  ) + conjg( 0.5d0*this%jd*(Sx(js) - ui*Sy(js))*exi ) ! + conjg
!      ham(4*k  , 4*k-1) = ham(4*k  , 4*k-1) + conjg( 0.5d0*this%jd*(Sx(js) + ui*Sy(js))*cexi ) ! + conjg
!    end do
!    do i = 1, this%pl_rh_1%length()  ! rashba (upper-right path)
!      js = this%pl_h%value(i)%i ! h-y   |   h-x
!      ks = this%pl_h%value(i)%f ! h+x   |   h+y
!      j = this%ste(js)
!      k = this%ste(ks)
!      exi  = exp(0.5d0*ui*(xi(js)+xi(ks)))
!      cexi = conjg(exi)
!      echi  = exp(0.5d0*ui*(chi(js)-chi(ks)))
!      cechi = conjg(echi)
!      ham(4*j-2, 4*k-3) =    echi* exi*this%lm
!      ham(4*j-3, 4*k-2) = -  echi*cexi*this%lm
!      ham(4*j  , 4*k-1) =   cechi*cexi*this%lm ! + conjg
!      ham(4*j-1, 4*k  ) = - cechi* exi*this%lm ! + conjg
!      ham(4*k-2, 4*j-3) =   cechi*cexi*this%lm
!      ham(4*k-3, 4*j-2) = - cechi* exi*this%lm
!      ham(4*k  , 4*j-1) =    echi* exi*this%lm ! + conjg
!      ham(4*k-1, 4*j  ) = -  echi*cexi*this%lm ! + conjg
!    end do
!    do i = 1, this%pl_rh_2%length()  ! rashba (upper-left path)
!      js = this%pl_h%value(i)%i ! h-y   |   h+x
!      ks = this%pl_h%value(i)%f ! h-x   |   h+y
!      j = this%ste(js)
!      k = this%ste(ks)
!      exi  = exp(0.5d0*ui*(xi(js)+xi(ks)))
!      cexi = conjg(exi)
!      echi  = exp(0.5d0*ui*(chi(js)-chi(ks)))
!      cechi = conjg(echi)
!      ham(4*j-2, 4*k-3) =   ui* echi* exi*this%lm
!      ham(4*j-3, 4*k-2) =   ui* echi*cexi*this%lm
!      ham(4*j  , 4*k-1) = - ui*cechi*cexi*this%lm ! + conjg
!      ham(4*j-1, 4*k  ) = - ui*cechi* exi*this%lm ! + conjg
!      ham(4*k-2, 4*j-3) = - ui*cechi*cexi*this%lm
!      ham(4*k-3, 4*j-2) = - ui*cechi* exi*this%lm
!      ham(4*k  , 4*j-1) =   ui* echi* exi*this%lm ! + conjg
!      ham(4*k-1, 4*j  ) =   ui* echi*cexi*this%lm ! + conjg
!    end do
!    !!!!!!!!!!!!!!    end bulk    !!!!!!!!!!!!!!!!!!
!
!    !!!!!!!!!!!!!!!!!!    interlayer    !!!!!!!!!!!!!!!!
!    ! pl_z_l(i) is z-dir path list between ith layer and (i+1)th layer.
!    ! here, pl_z_l(1) is interlayer path ( between 1st(surfaca) and 2nd(bulk) ).
!    do i = 1, this%pl_z_l(1)%length()
!      js = this%pl_z_l(1)%value(i)%i
!      ks = this%pl_z_l(1)%value(i)%f
!      j = this%ste(js)
!      k = this%ste(ks)
!      !exi  = exp(0.5d0*ui*modulo(xi(js)-xi(ks), 2d0*pi))
!      !exi  = exp(-0.5d0*ui*modulo(xi(js)-xi(ks), 2d0*pi))
!      !exi  = exp(0.5d0*ui*(xi(js)-xi(ks)))
!      exi  = exp(-0.5d0*ui*(xi(js)-xi(ks)))
!      !exi  = exp(-0.5d0*ui*xi(ks))
!      !exi  = exp(0.5d0*ui*xi(ks))
!      cexi = conjg(exi)
!      echi  = exp(0.5d0*ui*(chi(js)-chi(ks)))
!      cechi = conjg(echi)
!      ham(4*j-3, 4*k-3) = -  echi* exi*this%t(3)*(1d0-ne_d(js))*(1d0-ne_d(ks))
!      ham(4*j-2, 4*k-2) = -  echi*cexi*this%t(3)*(1d0-ne_u(js))*(1d0-ne_u(ks))
!      ham(4*j-1, 4*k-1) =   cechi*cexi*this%t(3)*(1d0-ne_d(js))*(1d0-ne_d(ks))
!      ham(4*j  , 4*k  ) =   cechi* exi*this%t(3)*(1d0-ne_u(js))*(1d0-ne_u(ks))
!      ham(4*k-3, 4*j-3) = - cechi*cexi*this%t(3)*(1d0-ne_d(js))*(1d0-ne_d(ks))
!      ham(4*k-2, 4*j-2) = - cechi* exi*this%t(3)*(1d0-ne_u(js))*(1d0-ne_u(ks))
!      ham(4*k-1, 4*j-1) =    echi* exi*this%t(3)*(1d0-ne_d(js))*(1d0-ne_d(ks))
!      ham(4*k  , 4*j  ) =    echi*cexi*this%t(3)*(1d0-ne_u(js))*(1d0-ne_u(ks))
!      
!      !exi  = exp(ui*xi(js))
!      exi  = exp(-ui*xi(js))
!      !exi  = 1d0
!      cexi = conjg(exi)
!      ham(4*j-3, 4*j-3) = ham(4*j-3, 4*j-3) - t3U*Sz(ks)
!      ham(4*j-2, 4*j-2) = ham(4*j-2, 4*j-2) + t3U*Sz(ks)
!      ham(4*j-3, 4*j-2) = ham(4*j-3, 4*j-2) - t3U*(Sx(ks) - ui*Sy(ks))*exi
!      ham(4*j-2, 4*j-3) = ham(4*j-2, 4*j-3) - t3U*(Sx(ks) + ui*Sy(ks))*cexi
!      ham(4*j-1, 4*j-1) = ham(4*j-1, 4*j-1) + t3U*Sz(ks) ! -conjg
!      ham(4*j  , 4*j  ) = ham(4*j  , 4*j  ) - t3U*Sz(ks) ! -conjg
!      ham(4*j-1, 4*j  ) = ham(4*j-1, 4*j  ) - conjg( t3U*(Sx(ks) - ui*Sy(ks))*exi ) ! + conjg
!      ham(4*j  , 4*j-1) = ham(4*j  , 4*j-1) - conjg( t3U*(Sx(ks) + ui*Sy(ks))*cexi ) ! + conjg
!      
!      !exi  = exp(ui*xi(ks))
!      exi  = exp(-ui*xi(ks))
!      cexi = conjg(exi)
!      ham(4*k-3, 4*k-3) = ham(4*k-3, 4*k-3) - t3U*Sz(js)
!      ham(4*k-2, 4*k-2) = ham(4*k-2, 4*k-2) + t3U*Sz(js)
!      ham(4*k-3, 4*k-2) = ham(4*k-3, 4*k-2) - t3U*(Sx(js) - ui*Sy(js))*exi
!      ham(4*k-2, 4*k-3) = ham(4*k-2, 4*k-3) - t3U*(Sx(js) + ui*Sy(js))*cexi
!      ham(4*k-1, 4*k-1) = ham(4*k-1, 4*k-1) + t3U*Sz(js) ! -conjg
!      ham(4*k  , 4*k  ) = ham(4*k  , 4*k  ) - t3U*Sz(js) ! -conjg
!      ham(4*k-1, 4*k  ) = ham(4*k-1, 4*k  ) - conjg( t3U*(Sx(js) - ui*Sy(js))*exi ) ! + conjg
!      ham(4*k  , 4*k-1) = ham(4*k  , 4*k-1) - conjg( t3U*(Sx(js) + ui*Sy(js))*cexi ) ! + conjg
!    end do
!    !!!!!!!!!!!!!!   end  interlayer    !!!!!!!!!!!!!!!!
!  
!    call my_zheev("l", ham, this%eg, this%wf)
!
!  end subroutine diagonalize_sb_2lay_xichi_hamiltonian

!  subroutine calc_eg_sb_2lay_xi(this, mu, delta, ne_u, ne_d, Sx, Sy, Sz, xi)
!    ! z =1 : surface
!    ! z =2 : bulk
!    implicit none
!    class(diag), intent(inout)   :: this
!    !real(8), dimension(4*this%n_acc), intent(out) :: eg
!    !complex(8),dimension(4*this%n_acc, 4*this%n_acc), intent(out)   :: wf
!    real(8),intent(in)    :: mu(2)
!    complex(8), dimension(this%pl_1st_l(1)%length()), intent(in) :: delta
!    real(8), dimension(this%n), intent(in)   :: ne_u, ne_d
!    !real(8), dimension(this%n), intent(in) :: xi, chi
!    real(8), dimension(this%n), intent(in) :: Sx, Sy, Sz, xi
!    complex(8),DIMENSION(4*this%n_acc, 4*this%n_acc)   :: twf
!    !real(8), dimension(this%n) :: ne_u, ne_d
!    integer :: i, j, k, is, js, ks, r
!    real(8) :: t3U
!    complex(8) :: exi,cexi
!
!    !! error cheack
!    if(this%nz /= 2 .or. this%nh_l(1) /= 0) then
!      print *, "nz : ", this%nz
!      print *, "nh(1) : ", this%nh_l(1)
!      stop "error : wrong data in calc_eg_sb_2lay_xi_hamiltonian"
!      ! this subroutine is available for only 2 layer"
!    endif
!
!    t3U = 2*this%t(3)**2/this%U
!
!    !ne_u = 0.5d0
!    !ne_d = 0.5d0
!
!    this%eg = 0d0
!    twf = this%wf
!
!   ! do i = 1, this%n_acc
!   !   is = this%ets(i)
!   !   twf(4*i-3, :) = exp( 0.5d0*ui*xi(is))*twf(4*i-3, :)
!   !   twf(4*i-2, :) = exp(-0.5d0*ui*xi(is))*twf(4*i-2, :)
!   !   twf(4*i-1, :) = exp(-0.5d0*ui*xi(is))*twf(4*i-1, :)
!   !   twf(4*i  , :) = exp( 0.5d0*ui*xi(is))*twf(4*i  , :)
!   !   
!   ! !  twf(4*i-3, :) = exp(-0.5d0*ui*xi(is))*twf(4*i-3, :)
!   ! !  twf(4*i-2, :) = exp( 0.5d0*ui*xi(is))*twf(4*i-2, :)
!   ! !  twf(4*i-1, :) = exp( 0.5d0*ui*xi(is))*twf(4*i-1, :)
!   ! !  twf(4*i  , :) = exp(-0.5d0*ui*xi(is))*twf(4*i  , :)
!   ! end do
!
!
!    do r = 1, 4*this%n_acc
!    !!!!!!!!!!!!!!!!!!    surface    !!!!!!!!!!!!!!!!!!
!    do i = this%n_acc_start(1), this%n_acc_end(1)
!      this%eg(r) = this%eg(r) - mu(1)*conjg(twf(4*i-3, r))*twf(4*i-3, r)
!      this%eg(r) = this%eg(r) - mu(1)*conjg(twf(4*i-2, r))*twf(4*i-2, r)
!      this%eg(r) = this%eg(r) + mu(1)*conjg(twf(4*i-1, r))*twf(4*i-1, r)
!      this%eg(r) = this%eg(r) + mu(1)*conjg(twf(4*i  , r))*twf(4*i  , r)
!    end do
!    do i = 1, this%pl_1st_l(1)%length()
!      js = this%pl_1st_l(1)%value(i)%i
!      ks = this%pl_1st_l(1)%value(i)%f
!      j = this%ste(js) ! in the 1st layer (surface), js = j, ks = k
!      k = this%ste(ks)
!      !exi  = exp(0.5d0*ui*modulo(xi(js)-xi(ks), 2d0*pi))
!      exi  = exp(-0.5d0*ui*modulo(xi(js)-xi(ks), 2d0*pi))
!      !exi  = 1d0
!      cexi = conjg(exi)
!      this%eg(r) = this%eg(r) - exi*this%t(1)*(1d0-ne_d(js))*(1d0-ne_d(ks))*conjg(twf(4*j-3, r))*twf(4*k-3, r)
!      this%eg(r) = this%eg(r) -cexi*this%t(1)*(1d0-ne_u(js))*(1d0-ne_u(ks))*conjg(twf(4*j-2, r))*twf(4*k-2, r)
!      this%eg(r) = this%eg(r) +cexi*this%t(1)*(1d0-ne_d(js))*(1d0-ne_d(ks))*conjg(twf(4*j-1, r))*twf(4*k-1, r)
!      this%eg(r) = this%eg(r) + exi*this%t(1)*(1d0-ne_u(js))*(1d0-ne_u(ks))*conjg(twf(4*j  , r))*twf(4*k  , r)
!      this%eg(r) = this%eg(r) -cexi*this%t(1)*(1d0-ne_d(js))*(1d0-ne_d(ks))*conjg(twf(4*k-3, r))*twf(4*j-3, r)
!      this%eg(r) = this%eg(r) - exi*this%t(1)*(1d0-ne_u(js))*(1d0-ne_u(ks))*conjg(twf(4*k-2, r))*twf(4*j-2, r)
!      this%eg(r) = this%eg(r) + exi*this%t(1)*(1d0-ne_d(js))*(1d0-ne_d(ks))*conjg(twf(4*k-1, r))*twf(4*j-1, r)
!      this%eg(r) = this%eg(r) +cexi*this%t(1)*(1d0-ne_u(js))*(1d0-ne_u(ks))*conjg(twf(4*k  , r))*twf(4*j  , r)
!      
!    end do
!    do i = 1, this%pl_2nd_l(1)%length()
!      js = this%pl_2nd_l(1)%value(i)%i
!      ks = this%pl_2nd_l(1)%value(i)%f
!      j = this%ste(js) ! in the 1st layer (surface), js = j, ks = k
!      k = this%ste(ks)
!      !exi  = exp(0.5d0*ui*modulo(xi(js)-xi(ks), 2d0*pi))
!      exi  = exp(-0.5d0*ui*modulo(xi(js)-xi(ks), 2d0*pi))
!      !exi  = 1d0
!      cexi = conjg(exi)
!      this%eg(r) = this%eg(r) - exi*this%t(2)*(1d0-ne_d(js))*(1d0-ne_d(ks))*conjg(twf(4*j-3, r))*twf(4*k-3, r)
!      this%eg(r) = this%eg(r) -cexi*this%t(2)*(1d0-ne_u(js))*(1d0-ne_u(ks))*conjg(twf(4*j-2, r))*twf(4*k-2, r)
!      this%eg(r) = this%eg(r) +cexi*this%t(2)*(1d0-ne_d(js))*(1d0-ne_d(ks))*conjg(twf(4*j-1, r))*twf(4*k-1, r)
!      this%eg(r) = this%eg(r) + exi*this%t(2)*(1d0-ne_u(js))*(1d0-ne_u(ks))*conjg(twf(4*j  , r))*twf(4*k  , r)
!      this%eg(r) = this%eg(r) -cexi*this%t(2)*(1d0-ne_d(js))*(1d0-ne_d(ks))*conjg(twf(4*k-3, r))*twf(4*j-3, r)
!      this%eg(r) = this%eg(r) - exi*this%t(2)*(1d0-ne_u(js))*(1d0-ne_u(ks))*conjg(twf(4*k-2, r))*twf(4*j-2, r)
!      this%eg(r) = this%eg(r) + exi*this%t(2)*(1d0-ne_d(js))*(1d0-ne_d(ks))*conjg(twf(4*k-1, r))*twf(4*j-1, r)
!      this%eg(r) = this%eg(r) +cexi*this%t(2)*(1d0-ne_u(js))*(1d0-ne_u(ks))*conjg(twf(4*k  , r))*twf(4*j  , r)
!    end do
!    do i = 1, this%pl_1st_l(1)%length()
!      ! here, delta's index corresponds to pl_1st_l(1)'s index
!      js = this%pl_1st_l(1)%value(i)%i
!      ks = this%pl_1st_l(1)%value(i)%f
!      j = this%ste(js) ! in the 1st layer ( surface ), js = j, ks = k
!      k = this%ste(ks)
!      !exi  = exp(0.5d0*ui*modulo(xi(js)-xi(ks), 2d0*pi))
!      exi  = exp(-0.5d0*ui*modulo(xi(js)-xi(ks), 2d0*pi))
!      !exi  = 1d0
!      cexi = conjg(exi)
!     ! ham(4*j-3, 4*k  ) = exi*delta(i)
!     ! ham(4*j-2, 4*k-1) = cexi*delta(i)
!     ! ham(4*j-1, 4*k-2) = cexi*conjg(delta(i))
!     ! ham(4*j  , 4*k-3) = exi*conjg(delta(i))
!     ! ham(4*k-3, 4*j  ) = cexi*delta(i)
!     ! ham(4*k-2, 4*j-1) = exi*delta(i)
!     ! ham(4*k-1, 4*j-2) = exi*conjg(delta(i))
!     ! ham(4*k  , 4*j-3) = cexi*conjg(delta(i))
!      this%eg(r) = this%eg(r) +        delta(i)* exi*conjg(twf(4*j-3, r))*twf(4*k  , r)
!      this%eg(r) = this%eg(r) +        delta(i)*cexi*conjg(twf(4*j-2, r))*twf(4*k-1, r)
!      this%eg(r) = this%eg(r) + conjg(delta(i))*cexi*conjg(twf(4*j-1, r))*twf(4*k-2, r)
!      this%eg(r) = this%eg(r) + conjg(delta(i))* exi*conjg(twf(4*j  , r))*twf(4*k-3, r)
!      this%eg(r) = this%eg(r) +        delta(i)*cexi*conjg(twf(4*k-3, r))*twf(4*j  , r)
!      this%eg(r) = this%eg(r) +        delta(i)* exi*conjg(twf(4*k-2, r))*twf(4*j-1, r)
!      this%eg(r) = this%eg(r) + conjg(delta(i))* exi*conjg(twf(4*k-1, r))*twf(4*j-2, r)
!      this%eg(r) = this%eg(r) + conjg(delta(i))*cexi*conjg(twf(4*k  , r))*twf(4*j-3, r)
!    end do
!    !!!!!!!!!!!!!!    end surface    !!!!!!!!!!!!!!!!!!
!
!    !!!!!!!!!!!!!!!!!!    bulk    !!!!!!!!!!!!!!!!
!    do i = this%n_acc_start(2), this%n_acc_end(2)
!      is = this%ets(i)
!      !exi  = exp(ui*xi(is))
!      exi  = exp(-ui*xi(is))
!      !exi=1d0
!      cexi = conjg(exi)
!      this%eg(r) = this%eg(r) + this%U*(-2d0/3d0*Sz(is)+0.5d0)*conjg(twf(4*i-3, r))*twf(4*i-3, r)
!      this%eg(r) = this%eg(r) + this%U*( 2d0/3d0*Sz(is)+0.5d0)*conjg(twf(4*i-2, r))*twf(4*i-2, r)
!      this%eg(r) = this%eg(r) - this%U*(-2d0/3d0*Sz(is)+0.5d0)*conjg(twf(4*i-1, r))*twf(4*i-1, r)
!      this%eg(r) = this%eg(r) - this%U*( 2d0/3d0*Sz(is)+0.5d0)*conjg(twf(4*i  , r))*twf(4*i  , r)
!      this%eg(r) = this%eg(r) - exi*this%U*2d0/3d0*(Sx(is)-ui*Sy(is))*conjg(twf(4*i-3, r))*twf(4*i-2, r)
!      this%eg(r) = this%eg(r) -cexi*this%U*2d0/3d0*(Sx(is)+ui*Sy(is))*conjg(twf(4*i-2, r))*twf(4*i-3, r)
!      this%eg(r) = this%eg(r) -cexi*this%U*2d0/3d0*(Sx(is)+ui*Sy(is))*conjg(twf(4*i-1, r))*twf(4*i  , r)
!      this%eg(r) = this%eg(r) - exi*this%U*2d0/3d0*(Sx(is)-ui*Sy(is))*conjg(twf(4*i  , r))*twf(4*i-1, r)
!      
!      !ham(4*i-3, 4*i-3) = this%U*(-2d0/3d0*Sz(is)+0.5d0)
!      !ham(4*i-2, 4*i-2) = this%U*( 2d0/3d0*Sz(is)+0.5d0)
!      !ham(4*i-3, 4*i-2) = -  exi*this%U*2d0/3d0*(Sx(is)-ui*Sy(is))
!      !ham(4*i-2, 4*i-3) = - cexi*this%U*2d0/3d0*(Sx(is)+ui*Sy(is))
!      !ham(4*i-1, 4*i-1) = - this%U*(-2d0/3d0*Sz(is)+0.5d0) ! - conjg
!      !ham(4*i  , 4*i  ) = - this%U*( 2d0/3d0*Sz(is)+0.5d0) ! - conjg
!      !ham(4*i-1, 4*i  ) = conjg(-  exi*this%U*2d0/3d0*(Sx(is)-ui*Sy(is))) ! + conjg
!      !ham(4*i  , 4*i-1) = conjg(- cexi*this%U*2d0/3d0*(Sx(is)+ui*Sy(is))) ! + conjg
!    end do
!    do i = this%n_acc_start(2), this%n_acc_end(2)
!      this%eg(r) = this%eg(r) - mu(2)*conjg(twf(4*i-3, r))*twf(4*i-3, r)
!      this%eg(r) = this%eg(r) - mu(2)*conjg(twf(4*i-2, r))*twf(4*i-2, r)
!      this%eg(r) = this%eg(r) + mu(2)*conjg(twf(4*i-1, r))*twf(4*i-1, r)
!      this%eg(r) = this%eg(r) + mu(2)*conjg(twf(4*i  , r))*twf(4*i  , r)
!!      ham(4*i-3, 4*i-3) = ham(4*i-3, 4*i-3) - mu(2)
!!      ham(4*i-2, 4*i-2) = ham(4*i-2, 4*i-2) - mu(2)
!!      ham(4*i-1, 4*i-1) = ham(4*i-1, 4*i-1) + mu(2)
!!      ham(4*i  , 4*i  ) = ham(4*i  , 4*i  ) + mu(2)
!    end do
!    do i = 1, this%pl_1st_l(2)%length()
!      js = this%pl_1st_l(2)%value(i)%i
!      ks = this%pl_1st_l(2)%value(i)%f
!      j = this%ste(js)
!      k = this%ste(ks)
!      !exi  = exp(0.5d0*ui*modulo(xi(js)-xi(ks), 2d0*pi))
!      exi  = exp(-0.5d0*ui*modulo(xi(js)-xi(ks), 2d0*pi))
!      !exi = 1d0
!      cexi = conjg(exi)
!      this%eg(r) = this%eg(r) - exi*this%t(1)*conjg(twf(4*j-3, r))*twf(4*k-3, r)
!      this%eg(r) = this%eg(r) -cexi*this%t(1)*conjg(twf(4*j-2, r))*twf(4*k-2, r)
!      this%eg(r) = this%eg(r) +cexi*this%t(1)*conjg(twf(4*j-1, r))*twf(4*k-1, r)
!      this%eg(r) = this%eg(r) + exi*this%t(1)*conjg(twf(4*j  , r))*twf(4*k  , r)
!      this%eg(r) = this%eg(r) -cexi*this%t(1)*conjg(twf(4*k-3, r))*twf(4*j-3, r)
!      this%eg(r) = this%eg(r) - exi*this%t(1)*conjg(twf(4*k-2, r))*twf(4*j-2, r)
!      this%eg(r) = this%eg(r) + exi*this%t(1)*conjg(twf(4*k-1, r))*twf(4*j-1, r)
!      this%eg(r) = this%eg(r) +cexi*this%t(1)*conjg(twf(4*k  , r))*twf(4*j  , r)
! !     ham(4*j-3, 4*k-3) = - exi*this%t(1)
! !     ham(4*j-2, 4*k-2) = -cexi*this%t(1)
! !     ham(4*j-1, 4*k-1) =  cexi*this%t(1)
! !     ham(4*j  , 4*k  ) =   exi*this%t(1)
! !     ham(4*k-3, 4*j-3) = -cexi*this%t(1)
! !     ham(4*k-2, 4*j-2) = - exi*this%t(1)
! !     ham(4*k-1, 4*j-1) =   exi*this%t(1)
! !     ham(4*k  , 4*j  ) =  cexi*this%t(1)
!    end do
!    do i = 1, this%pl_2nd_l(2)%length()
!      js = this%pl_2nd_l(2)%value(i)%i
!      ks = this%pl_2nd_l(2)%value(i)%f
!      j = this%ste(js)
!      k = this%ste(ks)
!      !exi  = exp(0.5d0*ui*modulo(xi(js)-xi(ks), 2d0*pi))
!      exi  = exp(-0.5d0*ui*modulo(xi(js)-xi(ks), 2d0*pi))
!      !exi = 1d0
!      cexi = conjg(exi)
!      this%eg(r) = this%eg(r) - exi*this%t(2)*conjg(twf(4*j-3, r))*twf(4*k-3, r)
!      this%eg(r) = this%eg(r) -cexi*this%t(2)*conjg(twf(4*j-2, r))*twf(4*k-2, r)
!      this%eg(r) = this%eg(r) +cexi*this%t(2)*conjg(twf(4*j-1, r))*twf(4*k-1, r)
!      this%eg(r) = this%eg(r) + exi*this%t(2)*conjg(twf(4*j  , r))*twf(4*k  , r)
!      this%eg(r) = this%eg(r) -cexi*this%t(2)*conjg(twf(4*k-3, r))*twf(4*j-3, r)
!      this%eg(r) = this%eg(r) - exi*this%t(2)*conjg(twf(4*k-2, r))*twf(4*j-2, r)
!      this%eg(r) = this%eg(r) + exi*this%t(2)*conjg(twf(4*k-1, r))*twf(4*j-1, r)
!      this%eg(r) = this%eg(r) +cexi*this%t(2)*conjg(twf(4*k  , r))*twf(4*j  , r)
!!      ham(4*j-3, 4*k-3) = - exi*this%t(2)
!!      ham(4*j-2, 4*k-2) = -cexi*this%t(2)
!!      ham(4*j-1, 4*k-1) =  cexi*this%t(2)
!!      ham(4*j  , 4*k  ) =   exi*this%t(2)
!!      ham(4*k-3, 4*j-3) = -cexi*this%t(2)
!!      ham(4*k-2, 4*j-2) = - exi*this%t(2)
!!      ham(4*k-1, 4*j-1) =   exi*this%t(2)
!!      ham(4*k  , 4*j  ) =  cexi*this%t(2)
!    end do
!    do i = 1, this%pl_h%length()
!      js = this%pl_h%value(i)%i
!      ks = this%pl_h%value(i)%f
!      j = this%ste(js)
!      k = this%ste(ks)
!      !exi  = exp(ui*xi(js))
!      exi  = exp(-ui*xi(js))
!      !exi  = exp(ui*xi(ks))
!      !exi = 1d0
!      cexi = conjg(exi)
!      
!      this%eg(r) = this%eg(r) + 0.5d0*this%jd*Sz(ks)*conjg(twf(4*j-3, r))*twf(4*j-3, r)
!      this%eg(r) = this%eg(r) - 0.5d0*this%jd*Sz(ks)*conjg(twf(4*j-2, r))*twf(4*j-2, r)
!      this%eg(r) = this%eg(r) - 0.5d0*this%jd*Sz(ks)*conjg(twf(4*j-1, r))*twf(4*j-1, r)
!      this%eg(r) = this%eg(r) + 0.5d0*this%jd*Sz(ks)*conjg(twf(4*j  , r))*twf(4*j  , r)
!      this%eg(r) = this%eg(r) + 0.5d0*this%jd*(Sx(ks) - ui*Sy(ks))*exi*conjg(twf(4*j-3, r))*twf(4*j-2, r)
!      this%eg(r) = this%eg(r) + 0.5d0*this%jd*(Sx(ks) + ui*Sy(ks))*cexi*conjg(twf(4*j-2, r))*twf(4*j-3, r)
!      this%eg(r) = this%eg(r) + 0.5d0*this%jd*(Sx(ks) + ui*Sy(ks))*cexi*conjg(twf(4*j-1, r))*twf(4*j  , r)
!      this%eg(r) = this%eg(r) + 0.5d0*this%jd*(Sx(ks) - ui*Sy(ks))*exi*conjg(twf(4*j  , r))*twf(4*j-1, r)
!!      ham(4*j-3, 4*j-3) = ham(4*j-3, 4*j-3) + 0.5d0*this%jd*Sz(ks)
!!      ham(4*j-2, 4*j-2) = ham(4*j-2, 4*j-2) - 0.5d0*this%jd*Sz(ks)
!!      ham(4*j-3, 4*j-2) = ham(4*j-3, 4*j-2) + 0.5d0*this%jd*(Sx(ks) - ui*Sy(ks))*exi
!!      ham(4*j-2, 4*j-3) = ham(4*j-2, 4*j-3) + 0.5d0*this%jd*(Sx(ks) + ui*Sy(ks))*cexi
!!      ham(4*j-1, 4*j-1) = ham(4*j-1, 4*j-1) - 0.5d0*this%jd*Sz(ks) ! -conjg
!!      ham(4*j  , 4*j  ) = ham(4*j  , 4*j  ) + 0.5d0*this%jd*Sz(ks) ! -conjg
!!      ham(4*j-1, 4*j  ) = ham(4*j-1, 4*j  ) + conjg( 0.5d0*this%jd*(Sx(ks) - ui*Sy(ks))*exi ) ! + conjg
!!      ham(4*j  , 4*j-1) = ham(4*j  , 4*j-1) + conjg( 0.5d0*this%jd*(Sx(ks) + ui*Sy(ks))*cexi ) ! + conjg
!!      
!      !exi  = exp(ui*xi(ks))
!      exi  = exp(-ui*xi(ks))
!      !exi  = exp(ui*xi(js))
!      !exi = 1d0
!      cexi = conjg(exi)
!      this%eg(r) = this%eg(r) + 0.5d0*this%jd*Sz(js)*conjg(twf(4*k-3, r))*twf(4*k-3, r)
!      this%eg(r) = this%eg(r) - 0.5d0*this%jd*Sz(js)*conjg(twf(4*k-2, r))*twf(4*k-2, r)
!      this%eg(r) = this%eg(r) - 0.5d0*this%jd*Sz(js)*conjg(twf(4*k-1, r))*twf(4*k-1, r)
!      this%eg(r) = this%eg(r) + 0.5d0*this%jd*Sz(js)*conjg(twf(4*k  , r))*twf(4*k  , r)
!      this%eg(r) = this%eg(r) + 0.5d0*this%jd*(Sx(js) - ui*Sy(js))*exi*conjg(twf(4*k-3, r))*twf(4*k-2, r)
!      this%eg(r) = this%eg(r) + 0.5d0*this%jd*(Sx(js) + ui*Sy(js))*cexi*conjg(twf(4*k-2, r))*twf(4*k-3, r)
!      this%eg(r) = this%eg(r) + 0.5d0*this%jd*(Sx(js) + ui*Sy(js))*cexi*conjg(twf(4*k-1, r))*twf(4*k  , r)
!      this%eg(r) = this%eg(r) + 0.5d0*this%jd*(Sx(js) - ui*Sy(js))*exi*conjg(twf(4*k  , r))*twf(4*k-1, r)
!!      ham(4*k-3, 4*k-3) = ham(4*k-3, 4*k-3) + 0.5d0*this%jd*Sz(js)
!!      ham(4*k-2, 4*k-2) = ham(4*k-2, 4*k-2) - 0.5d0*this%jd*Sz(js)
!!      ham(4*k-3, 4*k-2) = ham(4*k-3, 4*k-2) + 0.5d0*this%jd*(Sx(js) - ui*Sy(js))*exi
!!      ham(4*k-2, 4*k-3) = ham(4*k-2, 4*k-3) + 0.5d0*this%jd*(Sx(js) + ui*Sy(js))*cexi
!!      ham(4*k-1, 4*k-1) = ham(4*k-1, 4*k-1) - 0.5d0*this%jd*Sz(js) ! -conjg
!!      ham(4*k  , 4*k  ) = ham(4*k  , 4*k  ) + 0.5d0*this%jd*Sz(js) ! -conjg
!!      ham(4*k-1, 4*k  ) = ham(4*k-1, 4*k  ) + conjg( 0.5d0*this%jd*(Sx(js) - ui*Sy(js))*exi ) ! + conjg
!!      ham(4*k  , 4*k-1) = ham(4*k  , 4*k-1) + conjg( 0.5d0*this%jd*(Sx(js) + ui*Sy(js))*cexi ) ! + conjg
!    end do
!    !!!!!!!!!!!!!!    end bulk    !!!!!!!!!!!!!!!!!!
!!
!    !!!!!!!!!!!!!!!!!!    interlayer    !!!!!!!!!!!!!!!!
!    ! pl_z_l(i) is z-dir path list between ith layer and (i+1)th layer.
!    ! here, pl_z_l(1) is interlayer path ( between 1st(surfaca) and 2nd(bulk) ).
!    do i = 1, this%pl_z_l(1)%length()
!      js = this%pl_z_l(1)%value(i)%i
!      ks = this%pl_z_l(1)%value(i)%f
!      j = this%ste(js)
!      k = this%ste(ks)
!      !exi  = exp(0.5d0*ui*modulo(xi(js)-xi(ks), 2d0*pi))
!      exi  = exp(-0.5d0*ui*modulo(xi(js)-xi(ks), 2d0*pi))
!      !exi  = exp(-0.5d0*ui*xi(ks))
!      !exi  = exp(0.5d0*ui*xi(ks))
!      !exi = 1d0
!      cexi = conjg(exi)
!      this%eg(r) = this%eg(r) - exi*this%t(3)*(1d0-ne_d(js))*(1d0-ne_d(ks))*conjg(twf(4*j-3, r))*twf(4*k-3, r)
!      this%eg(r) = this%eg(r) -cexi*this%t(3)*(1d0-ne_u(js))*(1d0-ne_u(ks))*conjg(twf(4*j-2, r))*twf(4*k-2, r)
!      this%eg(r) = this%eg(r) +cexi*this%t(3)*(1d0-ne_d(js))*(1d0-ne_d(ks))*conjg(twf(4*j-1, r))*twf(4*k-1, r)
!      this%eg(r) = this%eg(r) + exi*this%t(3)*(1d0-ne_u(js))*(1d0-ne_u(ks))*conjg(twf(4*j  , r))*twf(4*k  , r)
!      this%eg(r) = this%eg(r) -cexi*this%t(3)*(1d0-ne_d(js))*(1d0-ne_d(ks))*conjg(twf(4*k-3, r))*twf(4*j-3, r)
!      this%eg(r) = this%eg(r) - exi*this%t(3)*(1d0-ne_u(js))*(1d0-ne_u(ks))*conjg(twf(4*k-2, r))*twf(4*j-2, r)
!      this%eg(r) = this%eg(r) + exi*this%t(3)*(1d0-ne_d(js))*(1d0-ne_d(ks))*conjg(twf(4*k-1, r))*twf(4*j-1, r)
!      this%eg(r) = this%eg(r) +cexi*this%t(3)*(1d0-ne_u(js))*(1d0-ne_u(ks))*conjg(twf(4*k  , r))*twf(4*j  , r)
!     ! ham(4*j-3, 4*k-3) = - exi*this%t(3)*(1d0-ne_d(js))*(1d0-ne_d(ks))
!     ! ham(4*j-2, 4*k-2) = -cexi*this%t(3)*(1d0-ne_u(js))*(1d0-ne_u(ks))
!     ! ham(4*j-1, 4*k-1) =  cexi*this%t(3)*(1d0-ne_d(js))*(1d0-ne_d(ks))
!     ! ham(4*j  , 4*k  ) =   exi*this%t(3)*(1d0-ne_u(js))*(1d0-ne_u(ks))
!     ! ham(4*k-3, 4*j-3) = -cexi*this%t(3)*(1d0-ne_d(js))*(1d0-ne_d(ks))
!     ! ham(4*k-2, 4*j-2) = - exi*this%t(3)*(1d0-ne_u(js))*(1d0-ne_u(ks))
!     ! ham(4*k-1, 4*j-1) =   exi*this%t(3)*(1d0-ne_d(js))*(1d0-ne_d(ks))
!     ! ham(4*k  , 4*j  ) =  cexi*this%t(3)*(1d0-ne_u(js))*(1d0-ne_u(ks))
!      
!      !exi  = exp(ui*xi(js))
!      exi  = exp(-ui*xi(js))
!      !exi  = 1d0
!      cexi = conjg(exi)
!      this%eg(r) = this%eg(r) - t3U*Sz(ks)*conjg(twf(4*j-3, r))*twf(4*j-3, r)
!      this%eg(r) = this%eg(r) + t3U*Sz(ks)*conjg(twf(4*j-2, r))*twf(4*j-2, r)
!      this%eg(r) = this%eg(r) + t3U*Sz(ks)*conjg(twf(4*j-1, r))*twf(4*j-1, r)
!      this%eg(r) = this%eg(r) - t3U*Sz(ks)*conjg(twf(4*j  , r))*twf(4*j  , r)
!      this%eg(r) = this%eg(r) - t3U*(Sx(ks) - ui*Sy(ks))*exi*conjg(twf(4*j-3, r))*twf(4*j-2, r)
!      this%eg(r) = this%eg(r) - t3U*(Sx(ks) + ui*Sy(ks))*cexi*conjg(twf(4*j-2, r))*twf(4*j-3, r)
!      this%eg(r) = this%eg(r) - t3U*(Sx(ks) + ui*Sy(ks))*cexi*conjg(twf(4*j-1, r))*twf(4*j  , r)
!      this%eg(r) = this%eg(r) - t3U*(Sx(ks) - ui*Sy(ks))*exi*conjg(twf(4*j  , r))*twf(4*j-1, r)
!      !ham(4*j-3, 4*j-3) = ham(4*j-3, 4*j-3) - t3U*Sz(ks)
!      !ham(4*j-2, 4*j-2) = ham(4*j-2, 4*j-2) + t3U*Sz(ks)
!      !ham(4*j-3, 4*j-2) = ham(4*j-3, 4*j-2) - t3U*(Sx(ks) - ui*Sy(ks))*exi
!      !ham(4*j-2, 4*j-3) = ham(4*j-2, 4*j-3) - t3U*(Sx(ks) + ui*Sy(ks))*cexi
!      !ham(4*j-1, 4*j-1) = ham(4*j-1, 4*j-1) + t3U*Sz(ks) ! -conjg
!      !ham(4*j  , 4*j  ) = ham(4*j  , 4*j  ) - t3U*Sz(ks) ! -conjg
!      !ham(4*j-1, 4*j  ) = ham(4*j-1, 4*j  ) - conjg( t3U*(Sx(ks) - ui*Sy(ks))*exi ) ! + conjg
!      !ham(4*j  , 4*j-1) = ham(4*j  , 4*j-1) - conjg( t3U*(Sx(ks) + ui*Sy(ks))*cexi ) ! + conjg
!      
!      !exi  = exp(ui*xi(ks))
!      exi  = exp(-ui*xi(ks))
!      !exi = 1d0
!      cexi = conjg(exi)
!      this%eg(r) = this%eg(r) - t3U*Sz(js)*conjg(twf(4*k-3, r))*twf(4*k-3, r)
!      this%eg(r) = this%eg(r) + t3U*Sz(js)*conjg(twf(4*k-2, r))*twf(4*k-2, r)
!      this%eg(r) = this%eg(r) + t3U*Sz(js)*conjg(twf(4*k-1, r))*twf(4*k-1, r)
!      this%eg(r) = this%eg(r) - t3U*Sz(js)*conjg(twf(4*k  , r))*twf(4*k  , r)
!      this%eg(r) = this%eg(r) - t3U*(Sx(js) - ui*Sy(js))*exi*conjg(twf(4*k-3, r))*twf(4*k-2, r)
!      this%eg(r) = this%eg(r) - t3U*(Sx(js) + ui*Sy(js))*cexi*conjg(twf(4*k-2, r))*twf(4*k-3, r)
!      this%eg(r) = this%eg(r) - t3U*(Sx(js) + ui*Sy(js))*cexi*conjg(twf(4*k-1, r))*twf(4*k  , r)
!      this%eg(r) = this%eg(r) - t3U*(Sx(js) - ui*Sy(js))*exi*conjg(twf(4*k  , r))*twf(4*k-1, r)
!     ! ham(4*k-3, 4*k-3) = ham(4*k-3, 4*k-3) - t3U*Sz(js)
!     ! ham(4*k-2, 4*k-2) = ham(4*k-2, 4*k-2) + t3U*Sz(js)
!     ! ham(4*k-3, 4*k-2) = ham(4*k-3, 4*k-2) - t3U*(Sx(js) - ui*Sy(js))*exi
!     ! ham(4*k-2, 4*k-3) = ham(4*k-2, 4*k-3) - t3U*(Sx(js) + ui*Sy(js))*cexi
!     ! ham(4*k-1, 4*k-1) = ham(4*k-1, 4*k-1) + t3U*Sz(js) ! -conjg
!     ! ham(4*k  , 4*k  ) = ham(4*k  , 4*k  ) - t3U*Sz(js) ! -conjg
!     ! ham(4*k-1, 4*k  ) = ham(4*k-1, 4*k  ) - conjg( t3U*(Sx(js) - ui*Sy(js))*exi ) ! + conjg
!     ! ham(4*k  , 4*k-1) = ham(4*k  , 4*k-1) - conjg( t3U*(Sx(js) + ui*Sy(js))*cexi ) ! + conjg
!    end do
!    !!!!!!!!!!!!!!   end  interlayer    !!!!!!!!!!!!!!!!
!!  
!!    call my_zheev("l", ham, this%eg, this%wf)
!
!    end do
!  end subroutine calc_eg_sb_2lay_xi

  subroutine diag_mod_calc_ne(this, ne_u, ne_d, kbt)
    class(diag), intent(in) :: this
    real(8), dimension(this%n), intent(out) :: ne_u, ne_d
    real(8), intent(in) :: kbt
    integer :: i, is, a
    ne_u = 0.d0
    ne_d = 0.d0
    do i = 1, this%n_acc
      is = this%ets(i)
      do a = 2*this%n_acc+1, 4*this%n_acc
        ne_u(is) = ne_u(is) + fermi( this%eg(a), kbT)*conjg(this%wf(4*i-3, a))*this%wf(4*i-3, a) &
                  &         + fermi(-this%eg(a), kbT)*conjg(this%wf(4*i-1, a))*this%wf(4*i-1, a)
        ne_d(is) = ne_d(is) + fermi( this%eg(a), kbT)*conjg(this%wf(4*i-2, a))*this%wf(4*i-2, a) &
                  &         + fermi(-this%eg(a), kbT)*conjg(this%wf(4*i  , a))*this%wf(4*i  , a) 
      end do!a
    end do!i
  end subroutine diag_mod_calc_ne
  
  subroutine diag_mod_calc_spin(this, Sx, Sy, Sz, kbt)
    class(diag), intent(in) :: this
    real(8), dimension(:), intent(out) :: Sx, Sy, Sz
    real(8), intent(in) :: kbt
    integer :: i, is, a
    complex(8) :: uu, ud, vu, vd
    Sx = 0d0
    Sy = 0d0
    Sz = 0d0
    do i = 1, this%n_acc
      is = this%ets(i)
      do a = 2*this%n_acc + 1, 4*this%n_acc 
        uu = this%wf(4*i-3, a)
        ud = this%wf(4*i-2, a)
        vu = this%wf(4*i-1, a)
        vd = this%wf(4*i  , a)

        !Sx(is) = Sx(is) - 0.5d0*(vu*conjg(vd) + conjg(vu)*vd)
        !Sy(is) = Sy(is) + 0.5d0*ui*(vu*conjg(vd) - conjg(vu)*vd)
        !Sz(is) = Sz(is) + 0.5d0*(vu*conjg(vu)-vd*conjg(vd))
        
        Sx(is) = Sx(is) + 0.5d0*(uu*conjg(ud) + conjg(uu)*ud)*fermi( this%eg(a), kbT) &
               &        - 0.5d0*(vu*conjg(vd) + conjg(vu)*vd)*fermi(-this%eg(a), kbT) 
        Sy(is) = Sy(is) - 0.5d0*ui*(uu*conjg(ud) - conjg(uu)*ud)*fermi( this%eg(a), kbT) &
               &        + 0.5d0*ui*(vu*conjg(vd) - conjg(vu)*vd)*fermi(-this%eg(a), kbT) 
        Sz(is) = Sz(is) + 0.5d0*(uu*conjg(uu) - ud*conjg(ud))*fermi( this%eg(a), kbT) &
               &        + 0.5d0*(vu*conjg(vu) - vd*conjg(vd))*fermi(-this%eg(a), kbT) 
      end do
    end do
  end subroutine diag_mod_calc_spin

  !subroutine diag_mod_calc_surface_delta(this, delta, kbt)
  !  class(diag), intent(in) :: this
  !  complex(8), dimension(:), intent(out) :: delta
  !  real(8), intent(in) :: kbt
  !  integer :: i, j, k, a
  !  !! error cheack
  !  if(this%nh /= 0) then
  !    print *, "nh : ", this%nh
  !    print *,"this subroutine can use only surface."
  !    stop "error : wrong data in calc_suface_delta"
  !  end if
  !  delta  = 0.d0
  !  do a = 2*this%n+1, 4*this%n
  !    do k = 1, this%pl_1st%length()
  !      i = this%pl_1st%value(k)%i
  !      j = this%pl_1st%value(k)%f
  !      ! jd : 2*t(1)**2 / U
  !      delta(k) = delta(k) + (this%jd*(this%wf(4*i-3, a)*conjg(this%wf(4*j  , a)) &
  !                &                   + this%wf(4*j-2, a)*conjg(this%wf(4*i-1, a)) &
  !                &                   + this%wf(4*i-2, a)*conjg(this%wf(4*j-1, a)) &
  !                &                   + this%wf(4*j-3, a)*conjg(this%wf(4*i  , a))))*tanh(this%eg(a)/kbT*0.5d0)
  !    end do!k
  !  end do!a
  !end subroutine diag_mod_calc_surface_delta
  
  subroutine diag_mod_calc_delta_sb_2lay(this, delta, kbt)
    class(diag), intent(in) :: this
    complex(8), dimension(:), allocatable, intent(out) :: delta
    real(8), intent(in) :: kbt
    integer :: i, j, k, a
    !! error cheack
    if(this%nz /= 2 .or. this%nh_l(1) /= 0) then
      print *, "nz : ", this%nz
      print *, "nh(1) : ", this%nh_l(1)
      stop "error : wrong data in calc_delta_sb_2lay"
      ! this subroutine is available for only 2 layer"
    endif
    if(allocated(delta)) deallocate(delta)
    allocate(delta(this%pl_1st_l(1)%length()))
    delta  = 0.d0
    do a = 2*this%n_acc+1, 4*this%n_acc
      do k = 1, this%pl_1st_l(1)%length()
        i = this%pl_1st_l(1)%value(k)%i
        j = this%pl_1st_l(1)%value(k)%f
        ! jd : 2*t(1)**2 / U
        delta(k) = delta(k) + (this%jd*(this%wf(4*i-3, a)*conjg(this%wf(4*j  , a)) &
                  &                   + this%wf(4*j-2, a)*conjg(this%wf(4*i-1, a)) &
                  &                   + this%wf(4*i-2, a)*conjg(this%wf(4*j-1, a)) &
                  &                   + this%wf(4*j-3, a)*conjg(this%wf(4*i  , a))))*tanh(this%eg(a)/kbT*0.5d0)
      end do!k
    end do!a
  end subroutine diag_mod_calc_delta_sb_2lay
  
 ! subroutine diag_mod_calc_delta_sb_2lay_xi(this, delta, xi, kbt)
 !   class(diag), intent(in) :: this
 !   complex(8), dimension(:), allocatable, intent(out) :: delta
 !   real(8), dimension(:), intent(in) :: xi
 !   real(8), intent(in) :: kbt
 !   integer :: i, j, k, a, is, js
 !   complex(8) :: exi, cexi
 !   !! error cheack
 !   if(this%nz /= 2 .or. this%nh_l(1) /= 0) then
 !     print *, "nz : ", this%nz
 !     print *, "nh(1) : ", this%nh_l(1)
 !     stop "error : wrong data in calc_delta_sb_2lay_xi"
 !     ! this subroutine is available for only 2 layer"
 !   endif
 !   if(allocated(delta)) deallocate(delta)
 !   allocate(delta(this%pl_1st_l(1)%length()))
 !   delta  = 0.d0
 !   do a = 2*this%n_acc+1, 4*this%n_acc
 !     do k = 1, this%pl_1st_l(1)%length()
 !       i = this%pl_1st_l(1)%value(k)%i
 !       j = this%pl_1st_l(1)%value(k)%f
 !       is = i
 !       js = j
 !       exi = exp(ui*(xi(is)-xi(js)))
 !       ! jd : 2*t(1)**2 / U
 !       delta(k) = delta(k) + (this%jd*(     this%wf(4*i-3, a)*conjg(this%wf(4*j  , a)) &
 !                 &                   +  exi*this%wf(4*j-2, a)*conjg(this%wf(4*i-1, a)) &
 !                 &                   +  exi*this%wf(4*i-2, a)*conjg(this%wf(4*j-1, a)) &
 !                 &                   +      this%wf(4*j-3, a)*conjg(this%wf(4*i  , a))))*tanh(this%eg(a)/kbT*0.5d0)
 !     !  exi  = exp(0.5d0*ui*modulo(xi(is)-xi(js), 2d0*pi))
 !     !  !exi  = exp(-0.5d0*ui*(xi(is)-xi(js)))
 !     !  cexi = conjg(exi)
 !     !  ! jd : 2*t(1)**2 / U
 !     !  delta(k) = delta(k) + (this%jd*( exi*this%wf(4*i-3, a)*conjg(this%wf(4*j  , a)) &
 !     !            &                   +  exi*this%wf(4*j-2, a)*conjg(this%wf(4*i-1, a)) &
 !     !            &                   + cexi*this%wf(4*i-2, a)*conjg(this%wf(4*j-1, a)) &
 !     !            &                   + cexi*this%wf(4*j-3, a)*conjg(this%wf(4*i  , a))))*tanh(this%eg(a)/kbT*0.5d0)
 !     end do!k
 !   end do!a
 ! end subroutine diag_mod_calc_delta_sb_2lay_xi
  
  subroutine diag_mod_calc_tdelta_sb_2lay_xi(this, tdelta_for, tdelta_rev, xi, kbt)
    ! tdelta_for : tilde{delta}_{ij} 
    ! tdelta_rev : tilde{delta}_{ji} 
    class(diag), intent(in) :: this
    !complex(8), dimension(:), allocatable, intent(out) :: delta
    complex(8), dimension(:), allocatable, intent(out) :: tdelta_for, tdelta_rev
    real(8), dimension(:), intent(in) :: xi
    real(8), intent(in) :: kbt
    integer :: i, j, k, a, is, js
    complex(8) :: exi, cexi
    !! error cheack
    if(this%nz /= 2 .or. this%nh_l(1) /= 0) then
      print *, "nz : ", this%nz
      print *, "nh(1) : ", this%nh_l(1)
      stop "error : wrong data in calc_delta_sb_2lay_xi"
      ! this subroutine is available for only 2 layer"
    endif
    !if(allocated(tdelta_for)) deallocate(tdelta_for)
    allocate(tdelta_for(this%pl_1st_l(1)%length()))
    allocate(tdelta_rev(this%pl_1st_l(1)%length()))
    tdelta_for  = 0.d0
    tdelta_rev  = 0.d0
    do a = 2*this%n_acc+1, 4*this%n_acc
      do k = 1, this%pl_1st_l(1)%length()
        i = this%pl_1st_l(1)%value(k)%i
        j = this%pl_1st_l(1)%value(k)%f
        is = i
        js = j
        exi = exp(-ui*(xi(is)-xi(js)))
        !exi = exp(ui*(xi(is)-xi(js)))
        cexi = conjg(exi)
        ! jd : 2*t(1)**2 / U
        tdelta_for(k) = tdelta_for(k) &
                  & + (this%jd*(     this%wf(4*i-3, a)*conjg(this%wf(4*j  , a)) &
                  &           +      this%wf(4*j-2, a)*conjg(this%wf(4*i-1, a)) &
                  &           +  exi*this%wf(4*i-2, a)*conjg(this%wf(4*j-1, a)) &
                  &           +  exi*this%wf(4*j-3, a)*conjg(this%wf(4*i  , a))))*tanh(this%eg(a)/kbT*0.5d0)
        tdelta_rev(k) = tdelta_rev(k) &
                  & + (this%jd*(cexi*this%wf(4*i-3, a)*conjg(this%wf(4*j  , a)) &
                  &           + cexi*this%wf(4*j-2, a)*conjg(this%wf(4*i-1, a)) &
                  &           +      this%wf(4*i-2, a)*conjg(this%wf(4*j-1, a)) &
                  &           +      this%wf(4*j-3, a)*conjg(this%wf(4*i  , a))))*tanh(this%eg(a)/kbT*0.5d0)
      end do!k
    end do!a
  end subroutine diag_mod_calc_tdelta_sb_2lay_xi
  
  subroutine diag_mod_calc_spin_xi(this, Sx, Sy, Sz, xi, kbt)
    class(diag), intent(in) :: this
    real(8), dimension(:), intent(out) :: Sx, Sy, Sz
    real(8), dimension(:), intent(in) :: xi
    real(8), intent(in) :: kbt
    integer :: i, is, a
    complex(8) :: uu, ud, vu, vd
    Sx = 0d0
    Sy = 0d0
    Sz = 0d0
    do i = 1, this%n_acc
      is = this%ets(i)
      do a = 2*this%n_acc + 1, 4*this%n_acc 
        uu = this%wf(4*i-3, a)*exp(-0.5d0*ui*xi(is))
        ud = this%wf(4*i-2, a)*exp( 0.5d0*ui*xi(is))
        vu = this%wf(4*i-1, a)*exp( 0.5d0*ui*xi(is))
        vd = this%wf(4*i  , a)*exp(-0.5d0*ui*xi(is))
        !uu = this%wf(4*i-3, a)*exp( 0.5d0*ui*xi(is))
        !ud = this%wf(4*i-2, a)*exp(-0.5d0*ui*xi(is))
        !vu = this%wf(4*i-1, a)*exp(-0.5d0*ui*xi(is))
        !vd = this%wf(4*i  , a)*exp( 0.5d0*ui*xi(is))

        !Sx(is) = Sx(is) - 0.5d0*(vu*conjg(vd) + conjg(vu)*vd)
        !Sy(is) = Sy(is) + 0.5d0*ui*(vu*conjg(vd) - conjg(vu)*vd)
        !Sz(is) = Sz(is) + 0.5d0*(vu*conjg(vu)-vd*conjg(vd))
        
        Sx(is) = Sx(is) + 0.5d0*(uu*conjg(ud) + conjg(uu)*ud)*fermi( this%eg(a), kbT) &
               &        - 0.5d0*(vu*conjg(vd) + conjg(vu)*vd)*fermi(-this%eg(a), kbT) 
        Sy(is) = Sy(is) - 0.5d0*ui*(uu*conjg(ud) - conjg(uu)*ud)*fermi( this%eg(a), kbT) &
               &        + 0.5d0*ui*(vu*conjg(vd) - conjg(vu)*vd)*fermi(-this%eg(a), kbT) 
        Sz(is) = Sz(is) + 0.5d0*(uu*conjg(uu) - ud*conjg(ud))*fermi( this%eg(a), kbT) &
               &        + 0.5d0*(vu*conjg(vu) - vd*conjg(vd))*fermi(-this%eg(a), kbT) 
        
        !Sx(is) = Sx(is) - real(vu*conjg(vd))
        !Sy(is) = Sy(is) - aimag(vu*conjg(vd))
        !Sz(is) = Sz(is) + real(vu*conjg(vu)-ud*conjg(vd))/2
      end do
    end do
  end subroutine diag_mod_calc_spin_xi

  subroutine diag_mod_calc_xi(this, xi, Sx, Sy)
    class(diag), intent(in) :: this
    real(8), dimension(:), intent(out) :: xi
    real(8), dimension(:), intent(in) :: Sx, Sy
    !real(8), intent(in) :: kbt
    integer :: i
    xi(:) = 0d0
    do i = 1, this%n
      if( Sx(1) /= 0d0 .and. Sy(2) /= 0d0) then
        xi(i) = atan2(Sy(i), Sx(i))
      end if
    end do
  end subroutine diag_mod_calc_xi
  
  !!!! public procedures

  function new_diag(pshape, hole, t, u, jd, lm) RESULT(r)
    type(diag)               :: r
    integer, dimension(:,:), intent(in)    :: pshape
    integer, dimension(:, :), intent(in) :: hole
    real(8), dimension(:), intent(in)    :: t
    real(8), intent(in)                  :: u, jd, lm
    CALL r%diag_init(pshape, hole, t, u, jd, lm)
  END FUNCTION new_diag

end module diagonalize_mod