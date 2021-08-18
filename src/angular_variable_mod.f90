module angular_variable_mod
  use base_mod,        only : base
  use path_list_mod,  only : path
  use math_mod,       only : pi, principal
  use hop_iterator_mod,   only : get_hopping_path
  use parameter_mod
  implicit none
  type, extends(base) :: angular_variable
  contains
    procedure :: angular_variable_init
    procedure :: angular_variable_3D_init
    procedure :: get_angular_variable    => angular_variable_get_angular_variable
    procedure :: get_af_angular_variable => angular_variable_get_af_angular_variable
    procedure :: rebuilt_angval_xi       => angular_variable_rebuilt_angval_xi
    procedure :: rebuilt_angval_chi      => angular_variable_rebuilt_angval_chi
    procedure :: get_and_rebuilt_chi    => angular_variable_get_and_rebuilt_chi
    procedure :: get_eta                => angular_variable_get_eta
    !procedure :: rebuild_eta            => angular_variable_rebuild_eta
    procedure :: d_angval1                => angular_variable_d_angval1
    procedure :: d_angval2                => angular_variable_d_angval2
    generic   :: d_angval_xi                => d_angval1, d_angval2
    procedure :: d_angval3                => angular_variable_d_angval3
    procedure :: d_angval4                => angular_variable_d_angval4
    generic   :: d_angval_chi                => d_angval3, d_angval4
    procedure :: get_rebuilt_phase       => angular_variable_get_rebuilt_phase
    procedure :: get_atan_phase_on_xy     => angular_variable_get_atan_phase_on_xy
  end type
contains
  subroutine angular_variable_init(this, pshape, hole)
    class(angular_variable), intent(inout)  :: this
    integer,intent(in)            :: pshape(2)
    integer,intent(in)            :: hole(:,:)
    call this%base_init_w_(pshape, hole)
  end subroutine angular_variable_init

  subroutine angular_variable_3D_init(this, pshape, hole)
    class(angular_variable), intent(inout)  :: this
    !integer,intent(in)            :: pshape(2)
    integer,intent(in)            :: pshape(:,:)
    integer,intent(in)            :: hole(:,:)
    !!!!! for 3D
    call this%base_init_w(pshape, hole)
    !!call this%base_init_w(pshape(:,1), hole)
    !! set parameters
    !!this%pshape = pshape(:,1)
    !!this%nx     = pshape(1,1)
    !!this%ny     = pshape(2,1)
    !this%n      = site_max(pshape)
    !if (allocated(this%nx)) deallocate(this%nx)
    !if (allocated(this%ny)) deallocate(this%ny)
    !allocate(this%nx(this%nz))
    !allocate(this%ny(this%nz))
    !this%nx = pshape(1,:)
    !this%ny = pshape(2,:)
    !this%nh     = size(hole,1)
    !this%ne     = layers(:)%ne
    !this%nz     = size(pshape(1,:))
    !this%pshape     = pshape
    !if (allocated(this%hole)) deallocate(this%hole)
    !allocate(this%hole(this%nh, 2), source = hole)
    !call this%set_site_to_e()
!   !set path
    !call this%cpl%alloc_value()
    !call this%pl_1st%alloc_value()
    !call this%pl_z%alloc_value()
    !call this%cpl%extend(get_hopping_path(pshape, [1,2,3,4,7], hole(:, 1)))
    !call this%pl_1st%extend(get_hopping_path(pshape, [1,2], hole(:, 1)))
    !call this%pl_z%extend(get_hopping_path(pshape, [7], hole(:, 1)))
    !this%np  = this%cpl%length()
    !this%np_1st  = this%pl_1st%length()
    !this%np_z  = this%pl_z%length()
    !!!!
  end subroutine angular_variable_3D_init

   function angular_variable_get_angular_variable (this) result(r)
    class(angular_variable), intent(in)           :: this
    real(8), dimension(:), allocatable          :: r
    integer                                     :: i,j,site,dy,dx
    if(allocated(r))deallocate(r)
    allocate(r(this%n))
    r = 0.d0
!   create init angular variable
    do j = 1, this%n_acc
      site = this%ets(j)
      do i = 1, this%nh
        if(this%z(site) == this%z(this%hole(i, 1)) ) then
          dy = this%y(site)-this%y(this%hole(i, 1))
          dx = this%x(site)-this%x(this%hole(i, 1))
          r(site) = r(site) + this%hole(i, 2) * atan2(real(dy,8),real(dx,8))
  !       write(*,*)this%hole(i,2)*atan2(real(dy,8),real(dx,8))
        end if
      enddo
 !     write(*,*)j
 !       if(j==2)stop
    enddo
 !   stop
  end function angular_variable_get_angular_variable

  function angular_variable_get_and_rebuilt_chi(this,hole_chi,val, fval) result(r)
   class(angular_variable), intent(in)           :: this
   real(8), dimension(:), allocatable          :: r
   integer                                     :: i,j,site,dy,dx
   integer, intent(in)                         :: hole_chi(this%nh, 2)
   real(8), intent(in), optional              :: fval
   real(8), dimension(this%n), intent(in)     :: val
   real(8), dimension(this%np)                :: dval
   if(allocated(r))deallocate(r)
   allocate(r(this%n))
   r = 0.d0
!   create init angular variable
   do j = 1, this%n_acc
     site = this%ets(j)
     do i = 1, this%nh
      if(this%z(site) == this%z(this%hole(i, 1)) ) then
        dy = this%y(site)-this%y(hole_chi(i, 1))
        dx = this%x(site)-this%x(hole_chi(i, 1))
        r(site) = r(site) + hole_chi(i, 2) * atan2(real(dy,8),real(dx,8))
!       write(*,*)hole_chi(i,2)*atan2(real(dy,8),real(dx,8))
      end if
    enddo
!     write(*,*)j
!       if(j==2)stop
   enddo
   dval = this%d_angval_chi(val)
   r = this%get_rebuilt_phase(dval, fval)
 end function angular_variable_get_and_rebuilt_chi

  function angular_variable_get_eta(this, cp_xi) result(eta)
    ! eta = xi - pi * (j_x + j_y) and -pi < eta <= pi
    ! We are going to make the rebuilt xi from the eta that is calculated from
    ! the xi of the Car-Parrinello method.
    class(angular_variable), intent(in)           :: this
    real(8), dimension(this%n)          :: eta
    integer                                     :: j,site
    real(8), dimension(this%n), intent(in)     :: cp_xi
    eta = 0.d0
    ! create init angular variable
    do j = 1, this%n_acc
      site = this%ets(j)
      eta(site) = cp_xi(site) - pi * real((this%x(site) + this%y(site) + this%z(site)), 8)
      eta(site) = modulo(eta(site), 2.d0 * pi)
      if(eta(site) < pi) eta(site) = eta(site) - 2.d0 * pi
    enddo
 end function angular_variable_get_eta

!  function angular_variable_rebuild_eta(this, eta) result(rebuilt_eta)
    ! add plus or minus pi to the eta in the upper right corner
    ! of the squares centered on the holes.
!    class(angular_variable), intent(in)           :: this
!    real(8), dimension(this%n)          :: rebuilt_eta
!    integer                                     :: i
!    real(8),  intent(in)            :: eta(this%n)
!    d_eta = 0.d0
!    rebuilt_eta = 0.d0
!    do i = 1,this%nh
!      
!    enddo
! end function angular_variable_rebuild_eta

   function angular_variable_get_af_angular_variable (this) result(r)
    class(angular_variable), intent(in)           :: this
    real(8), dimension(:), allocatable          :: r
    integer                                     :: i,j,site,dy,dx
    if(allocated(r))deallocate(r)
    allocate(r(this%n))
    r = 0.d0
!   create init angular variable
    do j = 1,this%n_acc
      site = this%ets(j)
      do i = 1, this%nh
        if(this%z(site) == this%z(this%hole(i, 1)) ) then
          dy = this%y(site)-this%y(this%hole(i, 1))
          dx = this%x(site)-this%x(this%hole(i, 1))
          r(site) = r(site) + this%hole(i, 2) * atan2(real(dy,8),real(dx,8))
        end if
      enddo
    enddo
!   add antiferro
    do i = 1,this%n
      r(i) = r(i) + (pi-1.e-2)*dble(this%x(i) + this%y(i) + this%z(i) - 1)
    enddo
  end function angular_variable_get_af_angular_variable

  function angular_variable_rebuilt_angval_xi(this, val, fval) result(r)
    class(angular_variable), intent(in)        :: this
    real(8), intent(in), optional              :: fval
    real(8), dimension(this%n), intent(in)     :: val
    real(8), dimension(this%n)                 :: r
    real(8), dimension(this%np)                :: dval
    dval = this%d_angval_xi(val)
    r = this%get_rebuilt_phase(dval, fval)
  end function angular_variable_rebuilt_angval_xi

  function angular_variable_rebuilt_angval_chi(this, val, fval) result(r)
    class(angular_variable), intent(in)        :: this
    real(8), intent(in), optional              :: fval
    real(8), dimension(this%n), intent(in)     :: val
    real(8), dimension(this%n)                 :: r
    real(8), dimension(this%np)                :: dval
    dval = this%d_angval_chi(val)
    r = this%get_rebuilt_phase(dval, fval)
  end function angular_variable_rebuilt_angval_chi

  function angular_variable_d_angval1(this, val) result(r)
    class(angular_variable), intent(in)        :: this
    real(8), dimension(this%n), intent(in)     :: val
    real(8), dimension(this%np)                :: r
    integer                                    :: i
    type(path)                                 :: p
    do i = 1, this%np
      p = this%cpl%value(i)
      r(i) = modulo(val(p%f) - val(p%i), 2.d0*pi)
    end do
  end function angular_variable_d_angval1

  function angular_variable_d_angval2(this, fphase, iphase) result(r)
    class(angular_variable), intent(in)        :: this
    real(8), intent(in)                        :: iphase, fphase
    real(8)                                    :: r
      r = modulo(fphase - iphase, 2.d0*pi)
  end function angular_variable_d_angval2

  function angular_variable_d_angval3(this, val) result(r)
    class(angular_variable), intent(in)        :: this
    real(8), dimension(this%n), intent(in)     :: val
    real(8), dimension(this%np)                :: r
    integer                                    :: i
    type(path)                                 :: p
    do i = 1, this%np
      p = this%cpl%value(i)
      r(i) = modulo(val(p%f) - val(p%i), 2.d0*pi)
     if (abs(r(i)) > pi) r(i)=r(i)-sign(1.d0,r(i))*2.d0*pi
    end do
  end function angular_variable_d_angval3

  function angular_variable_d_angval4(this, fphase, iphase) result(r)
    class(angular_variable), intent(in)        :: this
    real(8), intent(in)                        :: iphase, fphase
    real(8)                                    :: r
      r = modulo(fphase - iphase, 2.d0*pi)
  end function angular_variable_d_angval4

  function angular_variable_get_rebuilt_phase(this, dphase, fphase) result(r)
    class(angular_variable), intent(in)     :: this
    real(8), dimension(this%np), intent(in) :: dphase !difference between bonds
    real(8), intent(in), optional           :: fphase !first phase
    real(8), dimension(this%n)              :: r
    type(path)                              :: p
    integer                                 :: i
    logical, dimension(this%n)              :: flag
    integer                                 :: cnt, z, t
    logical,dimension(this%nz)              :: f_lay
    cnt = 1
    r(:) = 0.d0
    if (present(fphase)) then
      r(1) = fphase
    else
       r(1) = 0.d0!pi       pi/2.d0
    end if
    do i = 1, this%scpl%length()
      p = this%scpl%index(i)
      t = this%cpl%seek(p)
      r(p%f) = r(p%i) + dphase(t)
    end do

    !flag(:) = .false.
    !if (present(fphase)) then
    !  r(1) = fphase
    !else
    !   r(1) = 0.d0!pi       pi/2.d0
    !end if
    !flag(1) = .true.
    !f_lay(:) = .false.
    !f_lay(1) = .true.
    !do while(.not. all(f_lay))
    !  do i = 1, this%np_z
    !    p = this%pl_z%index(i)
    !    if (flag(p%i) .and. (.not. flag(p%f))) then
    !      t = this%cpl%seek(p)
    !      r(p%f) = r(p%i) + dphase(t)
    !      flag(p%f) = .true.
    !      z = site_to_z(p%f,this%pshape)
    !      f_lay(z) = .true.
    !      exit
    !    end if
    !  end do
    !end do
    !forall(i=1:this%nh) flag(this%hole(i, 1)) = .true.
    !do while(.not. all(flag))
    !!print *,flag(1:25)
    !!print *,flag(26:50)
    !  !do i = 1, this%np - this%np_z
    !  do i = 1, this%np_1st
    !    p = this%cpl%index(i)
    !    if (flag(p%i) .and. (.not. flag(p%f))) then
    !      r(p%f) = r(p%i) + dphase(i)
    !      flag(p%f) = .true.
    !      cycle !
    !    end if
    !  end do
    !  cnt = cnt + 1
    !end do
!    r(:) = modulo(r(:), 2.d0*pi)
  end function angular_variable_get_rebuilt_phase

  ! constructor
  function new_angular_variable (pshape, hole) result(r)
    type(angular_variable)        :: r
    integer,intent(in)            :: pshape(2)
    integer,intent(in)            :: hole(:, :)
    call r%angular_variable_init(pshape, hole)
  end function new_angular_variable
  ! constructor for 3D
  function new_angular_variable_3D (pshape, hole) result(r)
    type(angular_variable)        :: r
    integer,intent(in)            :: pshape(:, :)
    integer,intent(in)            :: hole(:, :)
    call r%angular_variable_3D_init(pshape, hole)
  end function new_angular_variable_3D

   !function angular_variable_atan_phase(this, pole, w, pole_xz, w_zx, pole_zy, w_zy) result(r)
   function angular_variable_get_atan_phase_on_xy(this, pole, w, z) result(r)
    class(angular_variable), intent(in)           :: this
    real(8), dimension(:,:), intent(in)           :: pole ! dimension(2, nh)
    integer, dimension(:), intent(in)             :: w
    integer, intent(in)                           :: z
    real(8), dimension(:), allocatable          :: r
    integer                                     :: i,j,site
    real(8) :: dy,dx
    if(allocated(r))deallocate(r)
    allocate(r(this%n))
    r = 0.d0
!   create init angular variable
    do j = 1, this%n_acc
      site = this%ets(j)
      do i = 1, size(pole(1,:))
        if(this%z(site) == z ) then
          dx = this%x(site)-pole(1, i)
          dy = this%y(site)-pole(2, i)
          r(site) = r(site) + w(i) * atan2(real(dy,8),real(dx,8))
        end if
      enddo
 !     write(*,*)j
 !       if(j==2)stop
    enddo
 !   stop
  end function angular_variable_get_atan_phase_on_xy
end module angular_variable_mod
