module base_mod
  use path_list_mod,      only : path, path_list
  use hop_iterator_mod
  use parameter_mod
  implicit none

  TYPE   ::  near_site
    INTEGER, ALLOCATABLE  :: s(:)
  END TYPE near_site
  
  type,abstract   :: base
    integer                              :: n, nz, nh, n_acc
    integer,dimension(:), allocatable    :: nl, nx, ny, nh_l
    type(cuo2_plane), dimension(:), allocatable :: layers
    integer                   :: np, np_1st, np_2nd, np_rh, np_z, np_h
    integer, dimension(:,:),allocatable      :: pshape
    integer, dimension(:),allocatable        :: ste, ets
    !integer, dimension(:),allocatable        :: site_to_acc, acc_to_site
    integer, dimension(:,:), allocatable    :: hole
    integer,dimension(:), allocatable    :: n_start, n_end, n_acc_start, n_acc_end
    real(8), dimension(:),allocatable        :: t, ne
    real(8)                                  :: u,jd,lm
    real(8),    dimension(:),    allocatable :: xi,chi
    complex(8), dimension(:, :), allocatable :: wf
    type(path_list)                          :: cpl
    TYPE(path_list)                       :: pl_1st, pl_2nd, pl_rh, pl_z, pl_h
    TYPE(path_list)                       :: pl_rh_1, pl_rh_2
    TYPE(path_list), dimension(:), allocatable :: pl_1st_l, pl_2nd_l, pl_z_l
    type(path_list)                          :: scpl ! simply connected path list
    TYPE(near_site), dimension(:),ALLOCATABLE :: ns_1st,ns_2nd, ns_z
  contains
    ! public type-bound procedures
    !procedure :: base_init => base_init_g
    !generic :: base_init => base_init_l,base_init_p
    generic :: base_init => base_init_p
    !procedure :: base_init_l => base_init_layer
    procedure :: base_init_p => base_init_pshape
    procedure :: base_modify_init
    procedure :: base_init_w_
    procedure :: base_init_w
    procedure :: base_init_lattice
    ! protected type-bound procedures
    procedure :: x             => base_x
    procedure :: y             => base_y
    procedure :: z             => base_z
    procedure :: r             => base_r
    !procedure :: site          => base_site
    procedure :: set_site_to_e => base_set_site_to_e
    procedure :: check_hole    => base_check_hole
    !procedure :: shift_site    => base_shift_site
    procedure :: shift_xyz_site    => base_shift_xyz_site
    procedure,private :: set_path    => base_set_path
    procedure,private :: set_scpl    => base_set_scpl
    procedure,private :: set_adjacent_site => base_set_adjacent_site
  end type base
contains
  ! public type-bound procedures
  subroutine base_init_pshape(this, pshape, hole, ne, t, u, jd, lm, xi, chi, wf)
    class(base), intent(inout)            :: this
    integer, intent(in)                   :: pshape(:,:)
    integer, intent(in)                   :: hole(:,:)
    real(8), intent(in)                   :: ne(:)
    real(8), intent(in)                   :: t(:)
    real(8), intent(in)                   :: u, jd, lm
    real(8), intent(in)                   :: xi(:), chi(:)
    complex(8), intent(in),optional       :: wf(:,:)
    integer i
    !type(cuo2_plane),intent(in),allocatable :: layers(:)
! set parameters
    this%pshape = pshape
    this%nz     = size(pshape(1,:))
    this%n      = site_max(pshape)
    if (allocated(this%nx)) deallocate(this%nx)
    if (allocated(this%ny)) deallocate(this%ny)
    if (allocated(this%nl)) deallocate(this%nl)
    if (allocated(this%ne)) deallocate(this%ne)
    allocate(this%nx(this%nz))
    allocate(this%ny(this%nz))
    allocate(this%nl(this%nz))
    allocate(this%nh_l(this%nz))
    allocate(this%ne(this%nz))
    this%nx(:) = pshape(1,:)
    this%ny(:) = pshape(2,:)
    this%nh = size(hole(:,1))
    do i = 1, this%nz
      this%nl(i) = this%nx(i) * this%ny(i)
    end do
    this%nh_l(:) = 0
    do i = 1, this%nh
      this%nh_l(this%z(hole(i,1))) = this%nh_l(this%z(hole(i,1))) + 1
    end do
    this%ne     = ne
    this%n_acc  = this%n - this%nh
    !this%layers     = layers
    allocate(this%xi(this%n), source = xi(1:this%n))
    allocate(this%chi(this%n), source = chi(1:this%n))
    allocate(this%wf(4*this%n_acc, 4*this%n_acc))
    if (present(wf)) this%wf = wf
    if (allocated(this%hole)) deallocate(this%hole)
    !allocate(this%hole(size(hole, 1), size(hole, 2)), source = hole)
    allocate(this%hole(this%nh,2))
    this%hole = hole
    if(.not. allocated(this%t)) allocate(this%t(size(t)), source = t)
    this%u  = u
    this%jd = jd
    this%lm = lm
    call this%set_site_to_e()
!   set path
    call this%cpl%alloc_value() ! allocate cpl%value
   ! call this%pl_1st%alloc_value()
   ! call this%pl_z%alloc_value()
   ! call this%pl_1st%extend(get_hopping_path(pshape, [1,2], hole(:, 1)))
   ! call this%pl_z%extend(get_hopping_path(pshape, [7], hole(:, 1)))
    call this%cpl%extend(get_hopping_path(pshape, [1,2,7], hole(:, 1)))
    call this%set_path()
    call this%set_adjacent_site()
    this%np  = this%cpl%length()
    this%np_1st  = this%pl_1st%length()
    this%np_2nd  = this%pl_2nd%length()
    this%np_z  = this%np - this%np_1st
    if (allocated(this%n_start)) deallocate(this%n_start)
    if (allocated(this%n_end)) deallocate(this%n_end)
    if (allocated(this%n_acc_start)) deallocate(this%n_acc_start)
    if (allocated(this%n_acc_end)) deallocate(this%n_acc_end)
    allocate(this%n_start(this%nz))
    allocate(this%n_end(this%nz))
    allocate(this%n_acc_start(this%nz))
    allocate(this%n_acc_end(this%nz))
    this%n_start(1) = 1
    this%n_end(1) = this%nx(1) * this%ny(1)
    do i = 2, this%nz
      this%n_start(i) = this%n_end(i-1) + 1
      this%n_end(i) = this%n_end(i-1) + this%nx(i) * this%ny(i)
    end do
    do i = 1, this%nz
      this%n_acc_start(i) = this%ste( this%n_start(i) )
      this%n_acc_end(i) = this%ste( this%n_end(i) )
    end do
    call this%set_scpl()
  end subroutine base_init_pshape

!  subroutine base_init_layer(this, layers, t, u, jd, lm, xi, chi, wf)
!    class(base), intent(inout)            :: this
!    !integer, intent(in)                   :: pshape(:,:)
!    !integer, intent(in)                   :: hole(:,:)
!    type(cuo2_plane),intent(in),allocatable :: layers(:)
!    real(8), intent(in)                   :: t(:)
!    real(8), intent(in)                   :: u, jd, lm
!    real(8), intent(in)                   :: xi(:), chi(:)
!    complex(8), intent(in),optional       :: wf(:,:)
!    integer i
!! set parameters
!    !this%pshape = pshape
!    this%pshape = make_pshape(layers)
!    !this%nx     = pshape(1,:)
!    !this%ny     = pshape(2,:)
!    !this%nl      = pshape(1,:) * pshape(2,:)
!    this%nz     = size(layers)
!    this%n      = layers(this%nz)%site_end
!    if (allocated(this%nx)) deallocate(this%nx)
!    if (allocated(this%ny)) deallocate(this%ny)
!    if (allocated(this%nl)) deallocate(this%nl)
!    if (allocated(this%ne)) deallocate(this%ne)
!    allocate(this%nx(this%nz))
!    allocate(this%ny(this%nz))
!    allocate(this%nl(this%nz))
!    allocate(this%nh_l(this%nz))
!    allocate(this%ne(this%nz))
!    this%nx(:) = layers(:)%nx
!    this%ny(:) = layers(:)%ny
!    this%nh = 0
!    do i = 1, this%nz
!      this%nl(i) = this%nx(i) * this%ny(i)
!      this%nh   =  this%nh + size(layers(i)%hole)
!      this%nh_l(i) = size(layers(i)%hole)
!    end do
!    this%ne     = layers(:)%ne
!    this%n_acc  = this%n - this%nh
!    this%layers     = layers
!    allocate(this%xi(this%n), source = xi)
!    allocate(this%chi(this%n), source = chi)
!    allocate(this%wf(4*this%n, 4*this%n))
!    if (present(wf)) this%wf = wf
!    if (allocated(this%hole)) deallocate(this%hole)
!    !allocate(this%hole(size(hole, 1), size(hole, 2)), source = hole)
!    allocate(this%hole(this%nh,2))
!    this%hole(:,1) = make_holearray(layers)
!    !this%hole(:,2) = ?
!    if(.not. allocated(this%t)) allocate(this%t(size(t)), source = t)
!    this%u  = u
!    this%jd = jd
!    this%lm = lm
!    call this%set_site_to_e()
!!   set path
!    call this%cpl%alloc_value() ! allocate cpl%value
!    call this%cpl%extend(get_hopping_path_l(this%layers, [1,2,7]))
!    !call this%pl_1st%alloc_value()
!    !call this%pl_z%alloc_value()
!    !call this%pl_1st%extend(get_hopping_path(this%pshape, [1,2], this%hole(:, 1)))
!    !call this%pl_z%extend(get_hopping_path(this%pshape, [7], this%hole(:, 1)))
!    call this%set_path()
!    call this%set_adjacent_site()
!    this%np  = this%cpl%length()
!    this%np_1st  = this%pl_1st%length()
!    this%np_z  = this%pl_z%length()
!    if (allocated(this%n_start)) deallocate(this%n_start)
!    if (allocated(this%n_end)) deallocate(this%n_end)
!    if (allocated(this%n_acc_start)) deallocate(this%n_acc_start)
!    if (allocated(this%n_acc_end)) deallocate(this%n_acc_end)
!    allocate(this%n_start(this%nz))
!    allocate(this%n_end(this%nz))
!    allocate(this%n_acc_start(this%nz))
!    allocate(this%n_acc_end(this%nz))
!    this%n_start(1) = 1
!    this%n_end(1) = this%nx(1) * this%ny(1)
!    do i = 2, this%nz
!      this%n_start(i) = this%n_end(i-1) + 1
!      this%n_end(i) = this%n_end(i-1) + this%nx(i) * this%ny(i)
!    end do
!    do i = 1, this%nz
!      this%n_acc_start(i) = this%ste( this%n_start(i) )
!      this%n_acc_end(i) = this%ste( this%n_end(i) )
!    end do
!    call this%set_scpl()
!  end subroutine base_init_layer
  subroutine base_modify_init(this,hole,xi,chi,wf)
    class(base), intent(inout)            :: this
    integer, intent(in)                             :: hole(this%nh, 2)
    real(8), intent(in)                             :: xi(this%n), chi(this%n)
    complex(8), intent(in), optional                :: wf(4*this%n_acc, 4*this%n_acc)
    this%xi(:) = xi(:)
    this%chi(:) = chi(:)
    this%hole(:,:) = hole(:,:)
    if(present(wf)) this%wf(:,:) = wf
  end subroutine base_modify_init

!  limited init
  subroutine base_init_w_(this, pshape, hole)
    class(base), intent(inout)            :: this
    integer, intent(in)                   :: pshape(2)
    integer, intent(in)                   :: hole(:,:)
    !type(cuo2_plane),intent(in),allocatable :: layers(:)

    call this%base_init_w(reshape(pshape,[2,1]),hole)
  end subroutine base_init_w_
  
  
  subroutine base_init_w(this, pshape, hole)
    class(base), intent(inout)            :: this
    integer, intent(in)                   :: pshape(:,:), hole(:,:)
    !type(cuo2_plane),intent(in),allocatable :: layers(:)
    integer i
    this%pshape = pshape
    this%nz     = size(pshape(1,:))
    if (allocated(this%ne)) deallocate(this%ne)
    allocate(this%ne(this%nz))
    if (allocated(this%nx)) deallocate(this%nx)
    allocate(this%nx(this%nz))
    if (allocated(this%ny)) deallocate(this%ny)
    allocate(this%ny(this%nz))
    do i = 1, this%nz
      this%nx(i) = pshape(1,i)
      this%ny(i) = pshape(2,i)
    end do
    this%n = 0
    do i = 1, this%nz
      this%n      = this%n + pshape(1,i) * pshape(2,i)
    end do
    this%nh     = size(hole(:,1))
    this%n_acc  = this%n - this%nh
    this%np_1st = this%pl_1st%length()
    this%np_2nd = this%pl_2nd%length()
    this%np_z = this%pl_z%length()
    if (allocated(this%hole)) deallocate(this%hole)
    allocate(this%hole(this%nh, 2), source = hole)
    call this%set_site_to_e()
!   set path
    call this%cpl%alloc_value()
    !call this%cpl%extend(get_hopping_path(pshape, [1,2]))
    call this%cpl%extend(get_hopping_path(pshape, [1,2,3,4,7],hole(:,1)))
    this%np  = this%cpl%length()

    call this%set_path()
    call this%set_adjacent_site()
    call this%set_scpl()
  end subroutine base_init_w

  subroutine base_init_lattice(this, pshape)
    class(base), intent(inout)            :: this
    integer, intent(in)                   :: pshape(:,:)
    this%pshape = pshape
    this%n = site_max(pshape)
  end subroutine base_init_lattice

  ! protected type-bound procedures
  !! return x-coord of <site_num>
  pure elemental integer function base_x(this, site_num) result(r)
    class(base), intent(in) :: this
    integer, intent(in)               :: site_num
    !r = mod(site_num-1, this%nx(1)) + 1
    r = site_to_x(site_num,this%pshape)
  end function base_x
  
  !! return y-coord of <site_num>
  pure elemental integer function base_y(this, site_num) result(r)
    class(base), intent(in) :: this
    integer, intent(in)               :: site_num
    !r = (site_num-1)/this%nx + 1
    r = site_to_y(site_num,this%pshape)
  end function base_y
  
  !! return z-coord of <site_num>
  pure elemental integer function base_z(this, site_num) result(r)
    class(base), intent(in) :: this
    integer, intent(in)               :: site_num
    !r = (site_num-1)/this%nx + 1
    r = site_to_z(site_num,this%pshape)
  end function base_z

  !! return r-coord of <site_num>
  pure function base_r(this, site_num) result(r)
    class(base), intent(in) :: this
    integer, intent(in)               :: site_num
    integer, dimension(3)       :: r
    !r = (site_num-1)/this%nx + 1
    r = site_to_r(site_num,this%pshape)
  end function base_r
  
  !! return site number of point (<x_coord>, <y_coord>)
  !pure integer function base_site(this, x_coord, y_coord, z_coord) result(r)
  !  class(base), intent(in) :: this
  !  integer, intent(in)               :: x_coord, y_coord, z_coord
  !  !r = x_coord + (y_coord-1)*this%nx
  !  r = xyz_to_site(x_coord,y_coord,z_coord,this%pshape)
  !end function base_site

! set electron site
  subroutine base_set_site_to_e(this)
    class(base), intent(inout)  :: this
    integer                                      :: i, j
    if(allocated(this%ste))deallocate(this%ste)
    allocate(this%ste(this%n))
    if(allocated(this%ets))deallocate(this%ets)
    !allocate(this%ets(this%n-this%nh))
    allocate(this%ets(this%n_acc))
    j = 0
    do i = 1, this%n
      if (this%check_hole(i)) then
        j = j+1
        this%ste(i) = 0
      else
        this%ste(i) = i-j
        this%ets(i-j) = i
      end if
    end do
  end subroutine base_set_site_to_e
  function base_check_hole(this, i) result(r)
    class(base), intent(in)  :: this
    integer, intent(in)                       :: i
    logical                                   :: r
    integer                                   :: j
    r = .false.
    do j = 1, this%nh
      if(i == this%hole(j, 1)) r = .true.
    end do
  end function base_check_hole

  !pure elemental integer function &
  !& base_shift_site(this, dir, i) result(r)
  !  class(base), intent(in) :: this
  !  integer, intent(in) :: dir, i
  !  r = 0
  !  select case(dir)
  !    case(1)
  !      if (.not. this%x(i)+1 > this%nx) then
  !        r = this%site(this%x(i)+1, this%y(i)  )
  !      end if
  !    case(2)
  !      if (.not. this%y(i)+1 > this%ny) then
  !        r = this%site(this%x(i)  , this%y(i)+1)
  !      end if
  !    case(3)
  !      if (.not. this%x(i)-1 < 1) then
  !        r = this%site(this%x(i)-1, this%y(i)  )
  !      end if
  !    case(4)
  !      if (.not. this%y(i)-1 < 1) then
  !        r = this%site(this%x(i)  , this%y(i)-1)
  !      end if
  !    case(5)
  !      if (.not. (this%x(i)+1 > this%nx .or. this%y(i)+1 > this%ny)) then
  !        r = this%site(this%x(i)+1, this%y(i)+1)
  !      end if
  !    case(6)
  !      if (.not. (this%x(i)-1 < 1       .or. this%y(i)+1 > this%ny)) then
  !        r = this%site(this%x(i)-1, this%y(i)+1)
  !      end if
  !    case(7)
  !      if (.not. (this%x(i)-1 < 1       .or. this%y(i)-1 < 1      )) then
  !        r = this%site(this%x(i)-1, this%y(i)-1)
  !      end if
  !    case(8)
  !      if (.not. (this%x(i) >= this%nx .or. this%y(i) <= 1      )) then
  !        r = this%site(this%x(i)+1, this%y(i)-1)
  !      end if
  !  end select
  !end function base_shift_site
  
  pure elemental integer function &
  & base_shift_xyz_site(this, dx, dy, dz, site) result(r)
    class(base), intent(in) :: this
    integer, intent(in) :: dx, dy, dz, site
    integer r_site(3)

    r_site = site_to_r(site,this%pshape)
    r_site = r_site + [dx, dy, dz]
    if( r_site(1) <= 0 .or. &
      & r_site(2) <= 0 .or. &
      & r_site(3) <= 0 .or. &
      & r_site(3) > this%nz) then
      !stop "Error"
      r = -1
    else if(r_site(1) > this%nx(r_site(3)) .or. r_site(2) > this%ny(r_site(3))) then
      r = -1
    else
      r = r_to_site(r_site, this%pshape)
    end if
    
  end function base_shift_xyz_site
  subroutine base_set_path(this)
    class(base), intent(inout)  :: this
    type(hop_iterator), allocatable           :: iter
    integer i, z
    call this%pl_1st%alloc_value()
    call this%pl_1st%extend(get_hopping_path(this%pshape, [1,2]))

    call this%pl_2nd%alloc_value()
    call this%pl_2nd%extend(get_hopping_path(this%pshape, [5,6]))

    call this%pl_rh%alloc_value()
    call this%pl_rh%extend(get_hopping_path(this%pshape, [3,4], this%hole(:,1)))
    
    call this%pl_rh_1%alloc_value()
    call this%pl_rh_1%extend(get_hopping_path(this%pshape, [3], this%hole(:,1)))
    call this%pl_rh_2%alloc_value()
    call this%pl_rh_2%extend(get_hopping_path(this%pshape, [4], this%hole(:,1)))
    
    !call this%pl_2nd%alloc_value()
    !call this%pl_2nd%extend(get_hopping_path(this%pshape, [5,6]))
    
    call this%pl_z%alloc_value()
    call this%pl_z%extend(get_hopping_path(this%pshape, [7]))
    
    call this%pl_h%alloc_value()
    allocate(iter, source = new_hop_iterator(this%pshape))
    call iter%set_path_hole(this%hole(:,1))
    call this%pl_h%extend(iter%pl)
    deallocate(iter)
    
    do i = 1, size(this%hole(:,1))
      call this%pl_1st%remove(this%hole(i,1))
      call this%pl_rh%remove(this%hole(i,1))
      call this%pl_2nd%remove(this%hole(i,1))
      call this%pl_z%remove(this%hole(i,1))
    end do

    ! for each layer path list
    if(allocated(this%pl_1st_l)) deallocate(this%pl_1st_l)
    allocate(this%pl_1st_l(this%nz))
    if(allocated(this%pl_2nd_l)) deallocate(this%pl_2nd_l)
    allocate(this%pl_2nd_l(this%nz))
    if(allocated(this%pl_z_l)) deallocate(this%pl_z_l)
    allocate(this%pl_z_l(this%nz-1))
    do z = 1, this%nz
      call this%pl_1st_l(z)%alloc_value()
      call this%pl_1st_l(z)%extend( this%pl_1st )
      call this%pl_2nd_l(z)%alloc_value()
      call this%pl_2nd_l(z)%extend( this%pl_2nd )
      !do i = this%xyz_to_site(1,1,z), this%xyz_to_site(this%nx(z),this%ny(z),z)
      do i = 1, this%n
        if(this%z(i) /= z) then
          call this%pl_1st_l(z)%remove(i)
          call this%pl_2nd_l(z)%remove(i)
        end if
      end do
    end do
    do z = 1, this%nz-1
      call this%pl_z_l(z)%alloc_value()
      call this%pl_z_l(z)%extend( this%pl_z )
      do i = 1, this%n
        if(this%z(i) /= z .and. this%z(i) /= z+1) then
          call this%pl_z_l(z)%remove(i)
        end if
      end do
    end do

    this%np_1st  = this%pl_1st%length()
    this%np_2nd  = this%pl_2nd%length()
    this%np_z  = this%pl_z%length()
    this%np_rh  = this%pl_rh%length()
    this%np_h  = this%pl_h%length()

  end subroutine base_set_path

  subroutine base_set_adjacent_site(this)
    class(base), intent(inout)  :: this
    integer i, is ,k
    type(path), allocatable :: p(:)

    if(allocated(this%ns_1st))deallocate(this%ns_1st)
    allocate(this%ns_1st(this%n))
    if(allocated(this%ns_2nd))deallocate(this%ns_2nd)
    allocate(this%ns_2nd(this%n))
    if(allocated(this%ns_z))deallocate(this%ns_z)
    allocate(this%ns_z(this%n))
    !do i = 1, this%n - this%nh
    do i = 1, this%n_acc
      is = this%ets(i)
      call this%pl_1st%seek3(is, p)
      allocate(this%ns_1st(is)%s(size(p)))
      do k = 1, size(p)
        if(p(k)%i == is) then
          this%ns_1st(is)%s(k) = p(k)%f
        elseif (p(k)%f == is) then
          this%ns_1st(is)%s(k) = p(k)%i
        else
          stop 'Error'
        end if
      end do
      call this%pl_2nd%seek3(is, p)
      allocate(this%ns_2nd(is)%s(size(p)))
      do k = 1, size(p)
        if(p(k)%i == is) then
          this%ns_2nd(is)%s(k) = p(k)%f
        elseif (p(k)%f == is) then
          this%ns_2nd(is)%s(k) = p(k)%i
        else
          stop 'Error'
        end if
      end do
      call this%pl_z%seek3(is, p)
      allocate(this%ns_z(is)%s(size(p)))
      do k = 1, size(p)
        if(p(k)%i == is) then
          this%ns_z(is)%s(k) = p(k)%f
        elseif (p(k)%f == is) then
          this%ns_z(is)%s(k) = p(k)%i
        else
          stop 'Error'
        end if
      end do
      deallocate(p)
    end do
  end subroutine base_set_adjacent_site

  subroutine base_set_scpl(this)
    ! scpl is simply connected path list.
    ! order of priority is x-dir path > y-dir path > z-dir path.
    ! z-dir path in scpl is only one each interlayer; scpl contains nz-1 paths.
    class(base), intent(inout)  :: this
    logical, dimension(this%n)  :: flag
    logical, dimension(this%nz)  :: f_lay
    type(path), allocatable :: p
    integer i,z, iter

    flag(:) = .false.
    flag(1) = .true.

    f_lay(:) = .false.
    f_lay(1) = .true.
    
    iter = 1
    call this%scpl%alloc_value()
    do while(.not. all(f_lay))
      do i = 1, this%np_z
        p = this%pl_z%index(i)
        if (flag(p%i) .and. (.not. flag(p%f))) then
          call this%scpl%append(p)
          flag(p%f) = .true.
          z = site_to_z(p%f,this%pshape)
          f_lay(z) = .true.
          exit
        end if
      end do
      iter = iter + 1
      if(iter > 10000) then
        stop "error : unable to creat simply connected path list, in base_set_scpl, base_mod"
      end if
    end do

    do i=1, this%nh
     flag(this%hole(i, 1)) = .true.
    end do
    
    do while(.not. all(flag))
      do i = 1, this%np_1st
        p = this%pl_1st%index(i)
        if (flag(p%i) .and. (.not. flag(p%f))) then
          call this%scpl%append(p)
          flag(p%f) = .true.
        end if
      end do
      iter = iter + 1
      if(iter > 10000) then
        stop "error : unable to creat simply connected path list, in base_set_scpl, base_mod"
      end if
    end do

  end subroutine base_set_scpl
end module base_mod
