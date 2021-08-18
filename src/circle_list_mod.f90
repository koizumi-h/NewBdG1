module circle_list_mod
  use base_mod,           only : base
  use math_mod,           only : pi, ui, principal, sign2, principal3
  use path_list_mod,      only : path, path_list
  use int_list_mod,       only : int_list, new_int_list
  !use cnode_iterator_mod, only : cnode_iterator, new_cnode_iterator
  !use cnode_creator_mod
  use parameter_mod
  use hop_iterator_mod,   only : get_hopping_path
  implicit none
  type, extends(base), public :: circle_list
    ! c1:small loop, c2:large loop, c:all loop
    ! ptoc:ptoc(bond_num) contain circle_num with + or -
    type(int_list), dimension(:), allocatable :: c1il, c2il, c3il, c4il, c5il, c6il, c7il, c8il,c9il,c10il,c11il, cil, ptoc
    integer                                   :: nc1, nc2, nc3, nc4, nc5 ,nc6, nc7, nc8, nc9, nc10, nc11, nc
    type(int_list), dimension(:), allocatable :: c_1st_il, c_rh_il, c_2nd_il, c_z_il
    !type(int_list), dimension(:), allocatable :: c_1st_il, c_rh_il, c_z_il
    integer                                   :: nc_1st, nc_rh, nc_2nd, nc_z
    !integer                                   :: nc_1st, nc_rh, nc_z
    !type(int_list), dimension(:,:), allocatable :: cpl_id !dim(dir, z)
  contains
    procedure :: circle_list_init
    procedure :: get_c1p           => circle_list_get_c1p
    procedure :: get_c2p           => circle_list_get_c2p
    procedure :: get_c3p           => circle_list_get_c3p
    procedure :: get_c4p           => circle_list_get_c4p
    procedure :: get_c5p           => circle_list_get_c5p
    procedure :: get_c6p           => circle_list_get_c6p
    procedure :: get_c7p           => circle_list_get_c7p
    procedure :: get_c8p           => circle_list_get_c8p
    procedure :: get_c9p           => circle_list_get_c9p
    procedure :: get_c10p           => circle_list_get_c10p
    procedure :: get_c11p           => circle_list_get_c11p
    procedure :: set_ptoc          => circle_list_set_ptoc
    procedure :: set_c1il          => circle_list_set_c1il
    !procedure :: set_c2il          => circle_list_set_c2il
    procedure :: set_c3il          => circle_list_set_c3il
    procedure :: set_c4il          => circle_list_set_c4il
    procedure :: set_c5il          => circle_list_set_c5il
    procedure :: set_c6il          => circle_list_set_c6il
    procedure :: set_c7il          => circle_list_set_c7il
    procedure :: set_c8il          => circle_list_set_c8il
    procedure :: set_c9il          => circle_list_set_c9il
    procedure :: set_c10il          => circle_list_set_c10il
    procedure :: set_c11il          => circle_list_set_c11il
    procedure :: set_cil           => circle_list_set_cil
    procedure :: get_delta_phase   => circle_list_get_delta_phase
    procedure :: check_cil         => circle_list_check_cil
    procedure :: get_circle_center         => circle_list_get_circle_center
    !procedure :: get_winding_number => circle_list_get_winding_number
    procedure :: judge_circle => circle_list_judge_circle
    procedure :: get_node_c_2nd => circle_list_get_node_c_2nd
    procedure :: get_node_c11 => circle_list_get_node_c11

  end type circle_list
contains
  ! public type-bound procedures
  subroutine circle_list_init(this, pshape, hole, path_2nd)
    class(circle_list), intent(inout)         :: this
    integer, dimension(:,:),intent(in)       :: pshape
    integer, dimension(:,:), intent(in)    :: hole
    logical, intent(in), optional   :: path_2nd
    !type(cuo2_plane), dimension(:), allocatable :: layers
    !character(1),optional                   :: rs
    integer i, j ! for debug
    logical :: flag

    flag = .true.
    if(present(path_2nd)) then
      flag = path_2nd
    end if

    ! if preset(rs) path=[1,2,3,4] ,else path=[1,2]
    !call this%base_init_w(layers,pshape,hole)
    call this%base_init_w(pshape,hole)
    
    call this%cpl%alloc_value()
    call this%pl_1st%alloc_value()
    call this%pl_rh%alloc_value()
    call this%pl_2nd%alloc_value()
    call this%pl_z%alloc_value()
    !call this%pl_h%alloc_value()
    
    if(flag) then
      call this%cpl%extend(get_hopping_path(pshape, [1,2,3,4,5,6,7], hole(:, 1)))
    else
      call this%cpl%extend(get_hopping_path(pshape, [1,2,3,4,7], hole(:, 1)))  !not include 2nd hop path
    endif

    call this%pl_1st%extend(get_hopping_path(pshape, [1,2], hole(:, 1)))
    call this%pl_rh%extend(get_hopping_path(pshape, [3,4], hole(:, 1)))
    if(flag) then
      call this%pl_2nd%extend(get_hopping_path(pshape, [5,6], hole(:, 1)))
    endif
    call this%pl_z%extend(get_hopping_path(pshape, [7], hole(:, 1)))
    this%np = this%cpl%length()
    this%np_1st = this%pl_1st%length()
    this%np_rh = this%pl_rh%length()
    this%np_2nd = this%pl_2nd%length()
    !this%np_2nd = this%pl_rh%length()
    this%np_z = this%pl_z%length()
    !this%np_h = this%pl_h%length()
    call this%set_c1il() ! x + y
    call this%set_c3il() ! around hole 1
    call this%set_c4il() ! around hole 2
    call this%set_c5il() ! around hole 3
    call this%set_c6il() ! around hole 4
    call this%set_c7il() ! around hole 5
    call this%set_c8il() ! 2nd direction(1,1) + x + y
    call this%set_c9il() ! 2nd direction(-1,1) + x + y
    call this%set_c10il() ! z + x
    call this%set_c11il() ! z + y
    call this%set_cil()
    this%nc1 = size(this%c1il)
    this%nc3 = size(this%c3il)
    this%nc4 = size(this%c4il)
    this%nc5 = size(this%c5il)
    this%nc6 = size(this%c6il)
    this%nc7 = size(this%c7il)
    this%nc8 = size(this%c8il)
    this%nc9 = size(this%c9il)
    this%nc10 = size(this%c10il)
    this%nc11 = size(this%c11il)
    !this%nc  = this%nc1 + this%nc3 + this%nc4 + this%nc5 + this%nc6 + this%nc7 + this%nc10 + this%nc11
    this%nc  = this%nc1 + this%nc3 + this%nc4 + this%nc5 + this%nc6 + this%nc7 &
          &  + this%nc8 + this%nc9 + this%nc10 + this%nc11
    !this%nc  = this%nc1 + this%nc3 + this%nc4 + this%nc5 + this%nc6 + this%nc7 &
    !      &  + this%nc10 + this%nc11
    !  call this%cpl%alloc_value()
    !  call this%cpl%extend(get_hopping_path(pshape, [1,2,3,4], hole(:, 1)))
    !  this%np = this%cpl%length()
    !  this%np_h = this%np - this%np_1st
      !print *,size(this%c3il)
      !print *,size(this%c4il)
      !print *,size(this%c5il)
      !print *,size(this%c6il)
      !print *,size(this%c7il)
    !  call this%set_c3il()
    !  call this%set_c4il()
    !  call this%set_c5il()
    !  call this%set_c6il()
    !  call this%set_c7il()
    !  this%nc3 = size(this%c3il)
    !  this%nc4 = size(this%c4il)
    !  this%nc5 = size(this%c5il)
    !  this%nc6 = size(this%c6il)
    !  this%nc7 = size(this%c7il)
    !  !this%nc  = this%nc1 + this%nc2 + this%nc3 + this%nc4 + this%nc5 + this%nc6
    !  this%nc  = this%nc1 + this%nc3 + this%nc4 + this%nc5 + this%nc6 + this%nc7
    !else
    !  this%np_h = 0
    !  allocate(this%c3il(0))
    !  allocate(this%c4il(0))
    !  allocate(this%c5il(0))
    !  allocate(this%c6il(0))
    !  allocate(this%c7il(0))
    !  this%nc3 = 0
    !  this%nc4 = 0
    !  this%nc5 = 0
    !  this%nc6 = 0
    !  this%nc7 = 0
    !endif
    !!if(present(hop2))then
    !if(.true.)then 
    !  call this%cpl%alloc_value()
    !  call this%cpl%extend(get_hopping_path(pshape, [1,2,3,4,5,6], hole(:, 1)))
    !  this%np = this%cpl%length()
    !  this%np_2nd = this%np - this%np_1st - this%np_h
    !  call this%set_c8il()
    !  call this%set_c9il()
    !  this%nc8 = size(this%c8il)
    !  this%nc9 = size(this%c9il)
    !  this%nc  = this%nc1 + this%nc3 + this%nc4 + this%nc5 + this%nc6 + this%nc7 + this%nc8 + this%nc9
    !else
    !  this%np_2nd = 0
    !  allocate(this%c8il(0))
    !  allocate(this%c9il(0))
    !  this%nc8 = 0
    !  this%nc9 = 0
    !end if


    !! z hop
    !if(.true.)then 
    !  call this%cpl%alloc_value()
    !  call this%cpl%extend(get_hopping_path(pshape, [1,2,3,4,5,6,7], hole(:, 1)))
    !  this%np = this%cpl%length()
    !  this%np_z = this%np - this%np_1st - this%np_h - this%np_2nd
    !  call this%set_c10il()
    !  call this%set_c11il()
    !  this%nc10 = size(this%c10il)
    !  this%nc11 = size(this%c11il)
    !  this%nc  = this%nc1 + this%nc3 + this%nc4 + this%nc5 + this%nc6  &
    !      &       + this%nc7 + this%nc8 + this%nc9 + this%nc10 + this%nc11 
    !else
    !  this%np_z = 0
    !  allocate(this%c10il(0))
    !  allocate(this%c11il(0))
    !  this%nc10 = 0
    !  this%nc11 = 0
    !end if
    !!call this%set_cil(rs)
    !call this%set_cil(rs,hop2)
      
      
      !!!!! debug
      !write(*,*) "nc is", this%nc
      !do i=1, this%nc
      !  write(*,*) "-------------",i,"----------------"
      !  do j=1, this%cil(i)%length()
      !    write(*,*) this%cil(i)%value(j), this%cpl%value(abs(this%cil(i)%value(j)))
      !  end do
      !  !call this%c8il(i)%display 
      !end do
      !stop
      ! end debug
    !!!!call this%set_ptoc()
  end subroutine circle_list_init

  ! path to circle init
!  subroutine circle_list_set_ptoc(this)
!    class(circle_list),intent(inout)  :: this
!    integer                              :: i, j, k, l, temp,ptemp(4)
!    integer, dimension(:), allocatable   :: a
!    if(allocated(this%ptoc)) deallocate(this%ptoc)
!    allocate(this%ptoc(this%np))
!    do i =1 ,this%np
!      l = 0
!      temp = 0
!      ptemp = 0
!     do j = 1, this%nc
!       temp = this%cil(j)%seek2(i)
!       if (temp/=0)then
!         l=l+1
!         ptemp(l)=sign(1,temp)*j
!       endif
!     enddo
!         write(*,*)i,ptemp,l
!     allocate(this%ptoc(i)%value(l))
!     do j = 1, l
!       this%ptoc(i)%value(j)=ptemp(j)
!     enddo
!   enddo
! end subroutine circle_list_set_ptoc

  ! path to circle init
  subroutine circle_list_set_ptoc(this)
    class(circle_list),intent(inout)  :: this
    integer                              :: i, j, k
    integer, dimension(:), allocatable   :: a, temp
    if(allocated(this%ptoc)) deallocate(this%ptoc)
    allocate(this%ptoc(this%np))
   do i = 1, this%nc
     allocate(temp(this%cil(i)%length()),source = this%cil(i)%value)
     do j = 1, size(temp) ! the numbers of path
       k = abs(temp(j)) ! path
       if (allocated(this%ptoc(k)%value))then
         allocate(a(this%ptoc(k)%length()), source = this%ptoc(k)%value)
         deallocate(this%ptoc(k)%value)
       else
         allocate(a(0))
       end if
        !write(*,*) 'f',i, j, k
       allocate(this%ptoc(k)%value(size(a)+1))
       this%ptoc(k)%value = [a, (sign(1,temp(j))*i)]
       deallocate(a)
        !write(*,*) 'g', j, k
     end do
     deallocate(temp)
   end do
 end subroutine circle_list_set_ptoc

  ! protected type-bound procedures
  ! start cil init
  function circle_list_get_c1p(this, i) result(r)
    class(circle_list), intent(in) :: this
    integer, intent(in)               :: i
    integer, dimension(4)             :: r
    integer site_r(3), n
    integer r1(3), r2(3)
    site_r = site_to_r(i,this%pshape)
    !r(1) = this%cpl%seek2(path(i, i+1))
    !r(2) = this%cpl%seek2(path(i+1, i+1+this%nx))
    !r(3) = this%cpl%seek2(path(i+1+this%nx, i+this%nx))
    !r(4) = this%cpl%seek2(path(i+this%nx, i))

    !r(1) = this%cpl%seek2(path(r_to_site(site_r + [0,0,0], this%pshape), r_to_site(site_r + [1,0,0], this%pshape)))
    !r(2) = this%cpl%seek2(path(r_to_site(site_r + [1,0,0], this%pshape), r_to_site(site_r + [1,1,0], this%pshape)))
    !r(3) = this%cpl%seek2(path(r_to_site(site_r + [1,1,0], this%pshape), r_to_site(site_r + [0,1,0], this%pshape)))
    !r(4) = this%cpl%seek2(path(r_to_site(site_r + [0,1,0], this%pshape), r_to_site(site_r + [0,0,0], this%pshape)))
    r1 = site_r
    r2 = site_r + [1,0,0]
    r(1) = this%cpl%seek2(path(r_to_site(r1, this%pshape), r_to_site(r2, this%pshape)))
    r1 = r2
    r2 = site_r + [1,1,0]
    r(2) = this%cpl%seek2(path(r_to_site(r1, this%pshape), r_to_site(r2, this%pshape)))
    r1 = r2
    r2 = site_r + [0,1,0]
    r(3) = this%cpl%seek2(path(r_to_site(r1, this%pshape), r_to_site(r2, this%pshape)))
    r1 = r2
    r2 = site_r
    r(4) = this%cpl%seek2(path(r_to_site(r1, this%pshape), r_to_site(r2, this%pshape)))
  end function circle_list_get_c1p

  function circle_list_get_c2p(this, i) result(r)
    class(circle_list), intent(in) :: this
    integer, intent(in)               :: i ! i is node
    integer, dimension(8)             :: r
    !r(1) = this%cpl%seek2(path(i            , i+1          ))
    !r(2) = this%cpl%seek2(path(i+1          , i+2          ))
    !r(3) = this%cpl%seek2(path(i+2          , i+2+  this%nx))
    !r(4) = this%cpl%seek2(path(i+2+  this%nx, i+2+2*this%nx))
    !r(5) = this%cpl%seek2(path(i+2+2*this%nx, i+1+2*this%nx))
    !r(6) = this%cpl%seek2(path(i+1+2*this%nx, i+  2*this%nx))
    !r(7) = this%cpl%seek2(path(i+  2*this%nx, i+this%nx    ))
    !r(8) = this%cpl%seek2(path(i+    this%nx, i            ))
  end function circle_list_get_c2p

  function circle_list_get_c3p(this, h) result(r)
    class(circle_list), intent(in) :: this
    !integer, intent(in)               :: i
    integer, intent(in)               :: h  ! h is hole site
    integer, dimension(3)             :: r
    integer site_r(3)
    integer r1(3), r2(3)
    site_r = site_to_r(h,this%pshape)
    !r(1) = this%cpl%seek2(path(i, i+1))
    !r(2) = this%cpl%seek2(path(i+1, i+this%nx))
    !r(3) = this%cpl%seek2(path(i+this%nx, i))
    !r(1) = this%cpl%seek2(path(r_to_site(site_r + [0,0,0], this%pshape), r_to_site(site_r + [1,0,0], this%pshape)))
    !r(1) = this%cpl%seek2(path(r_to_site(site_r + [1,0,0], this%pshape), r_to_site(site_r + [1,1,0], this%pshape)))
    !r(1) = this%cpl%seek2(path(r_to_site(site_r + [1,1,0], this%pshape), r_to_site(site_r + [0,0,0], this%pshape)))
    
    !r(1) = this%cpl%seek2(path(r_to_site(site_r + [-1,-1,0], this%pshape), r_to_site(site_r + [0,-1,0], this%pshape)))
    !r(2) = this%cpl%seek2(path(r_to_site(site_r + [0,-1,0], this%pshape), r_to_site(site_r + [-1,0,0], this%pshape)))
    !r(3) = this%cpl%seek2(path(r_to_site(site_r + [-1,0,0], this%pshape), r_to_site(site_r + [-1,-1,0], this%pshape)))
    r1 = site_r + [-1,-1,0]
    r2 = site_r + [0,-1,0]
    r(1) = this%cpl%seek2(path(r_to_site(r1, this%pshape), r_to_site(r2, this%pshape)))
    r1 = r2
    r2 = site_r + [-1,0,0]
    r(2) = this%cpl%seek2(path(r_to_site(r1, this%pshape), r_to_site(r2, this%pshape)))
    r1 = r2
    r2 = site_r + [-1,-1,0]
    r(3) = this%cpl%seek2(path(r_to_site(r1, this%pshape), r_to_site(r2, this%pshape)))
  end function circle_list_get_c3p

  function circle_list_get_c4p(this, h) result(r)
    class(circle_list), intent(in) :: this
    integer, intent(in)               :: h
    integer, dimension(3)             :: r
    integer site_r(3)
    integer r1(3), r2(3)
    site_r = site_to_r(h,this%pshape)
    !r(1) = this%cpl%seek2(path(r_to_site(site_r + [1,-1,0], this%pshape), r_to_site(site_r + [1,0,0], this%pshape)))
    !r(2) = this%cpl%seek2(path(r_to_site(site_r + [1,0,0], this%pshape), r_to_site(site_r + [0,-1,0], this%pshape)))
    !r(3) = this%cpl%seek2(path(r_to_site(site_r + [0,-1,0], this%pshape), r_to_site(site_r + [1,-1,0], this%pshape)))
    r1 = site_r + [1,-1,0]
    r2 = site_r + [1,0,0]
    r(1) = this%cpl%seek2(path(r_to_site(r1, this%pshape), r_to_site(r2, this%pshape)))
    r1 = r2
    r2 = site_r + [0,-1,0]
    r(2) = this%cpl%seek2(path(r_to_site(r1, this%pshape), r_to_site(r2, this%pshape)))
    r1 = r2
    r2 = site_r + [1,-1,0]
    r(3) = this%cpl%seek2(path(r_to_site(r1, this%pshape), r_to_site(r2, this%pshape)))
  end function circle_list_get_c4p

  function circle_list_get_c5p(this, h) result(r)
    class(circle_list), intent(in) :: this
    integer, intent(in)               :: h
    integer, dimension(3)             :: r
    integer site_r(3)
    integer r1(3), r2(3)
    site_r = site_to_r(h,this%pshape)
    !r(1) = this%cpl%seek2(path(r_to_site(site_r + [-1,1,0], this%pshape), r_to_site(site_r + [-1,0,0], this%pshape)))
    !r(2) = this%cpl%seek2(path(r_to_site(site_r + [-1,0,0], this%pshape), r_to_site(site_r + [0,1,0], this%pshape)))
    !r(3) = this%cpl%seek2(path(r_to_site(site_r + [0,1,0], this%pshape), r_to_site(site_r + [-1,1,0], this%pshape)))
    r1 = site_r + [-1,1,0]
    r2 = site_r + [-1,0,0]
    r(1) = this%cpl%seek2(path(r_to_site(r1, this%pshape), r_to_site(r2, this%pshape)))
    r1 = r2
    r2 = site_r + [0,1,0]
    r(2) = this%cpl%seek2(path(r_to_site(r1, this%pshape), r_to_site(r2, this%pshape)))
    r1 = r2
    r2 = site_r + [-1,1,0]
    r(3) = this%cpl%seek2(path(r_to_site(r1, this%pshape), r_to_site(r2, this%pshape)))
  end function circle_list_get_c5p

  function circle_list_get_c6p(this, h) result(r)
    class(circle_list), intent(in) :: this
    integer, intent(in)               :: h
    integer, dimension(3)             :: r
    integer site_r(3)
    integer r1(3), r2(3)
    site_r = site_to_r(h,this%pshape)
    !r(1) = this%cpl%seek2(path(r_to_site(site_r + [1,1,0], this%pshape), r_to_site(site_r + [0,1,0], this%pshape)))
    !r(2) = this%cpl%seek2(path(r_to_site(site_r + [0,1,0], this%pshape), r_to_site(site_r + [1,0,0], this%pshape)))
    !r(3) = this%cpl%seek2(path(r_to_site(site_r + [1,0,0], this%pshape), r_to_site(site_r + [1,1,0], this%pshape)))
    r1 = site_r + [1,1,0]
    r2 = site_r + [0,1,0]
    r(1) = this%cpl%seek2(path(r_to_site(r1, this%pshape), r_to_site(r2, this%pshape)))
    r1 = r2
    r2 = site_r + [1,0,0]
    r(2) = this%cpl%seek2(path(r_to_site(r1, this%pshape), r_to_site(r2, this%pshape)))
    r1 = r2
    r2 = site_r + [1,1,0]
    r(3) = this%cpl%seek2(path(r_to_site(r1, this%pshape), r_to_site(r2, this%pshape)))
  end function circle_list_get_c6p


  function circle_list_get_c7p(this, h) result(r)
    class(circle_list), intent(in) :: this
    integer, intent(in)               :: h
    integer, dimension(4)             :: r
    integer site_r(3)
    integer r1(3), r2(3)
    site_r = site_to_r(h,this%pshape)
    !r(1) = this%cpl%seek2(path(i, i+1+this%nx))
    !r(2) = this%cpl%seek2(path(i+1+this%nx, i+2*this%nx))
    !r(3) = this%cpl%seek2(path(i+2*this%nx, i-1+this%nx))
    !r(4) = this%cpl%seek2(path(i-1+this%nx, i))
    !r(1) = this%cpl%seek2(path(r_to_site(site_r + [0,-1,0], this%pshape), r_to_site(site_r + [1,0,0], this%pshape)))
    !r(2) = this%cpl%seek2(path(r_to_site(site_r + [1,0,0], this%pshape), r_to_site(site_r + [0,1,0], this%pshape)))
    !r(3) = this%cpl%seek2(path(r_to_site(site_r + [0,1,0], this%pshape), r_to_site(site_r + [-1,0,0], this%pshape)))
    !r(4) = this%cpl%seek2(path(r_to_site(site_r + [-1,0,0], this%pshape), r_to_site(site_r + [0,-1,0], this%pshape)))
    r1 = site_r + [0,-1,0]
    r2 = site_r + [1,0,0]
    r(1) = this%cpl%seek2(path(r_to_site(r1, this%pshape), r_to_site(r2, this%pshape)))
    r1 = r2
    r2 = site_r + [0,1,0]
    r(2) = this%cpl%seek2(path(r_to_site(r1, this%pshape), r_to_site(r2, this%pshape)))
    r1 = r2
    r2 = site_r + [-1,0,0]
    r(3) = this%cpl%seek2(path(r_to_site(r1, this%pshape), r_to_site(r2, this%pshape)))
    r1 = r2
    r2 = site_r + [0,-1,0]
    r(4) = this%cpl%seek2(path(r_to_site(r1, this%pshape), r_to_site(r2, this%pshape)))
  end function circle_list_get_c7p

  !     /|
  !   /  |
  ! /____| circle (i : point of left bottom)
  function circle_list_get_c8p(this, i) result(r)
    class(circle_list), intent(in) :: this
    integer, intent(in)               :: i
    integer, dimension(3)             :: r
    integer site_r(3)
    integer r1(3), r2(3)
    site_r = site_to_r(i,this%pshape)
    !r(1) = this%cpl%seek2(path(i, i+1))
    !r(2) = this%cpl%seek2(path(i+1, i+this%nx+1))
    !r(3) = this%cpl%seek2(path(i+this%nx+1, i))
    r1 = site_r
    r2 = site_r + [1,0,0]
    r(1) = this%cpl%seek2(path(r_to_site(r1, this%pshape), r_to_site(r2, this%pshape)))
    r1 = r2
    r2 = site_r + [1,1,0]
    r(2) = this%cpl%seek2(path(r_to_site(r1, this%pshape), r_to_site(r2, this%pshape)))
    r1 = r2
    r2 = site_r + [0,0,0]
    r(3) = this%cpl%seek2(path(r_to_site(r1, this%pshape), r_to_site(r2, this%pshape)))
  end function circle_list_get_c8p

  ! |\
  ! |  \
  ! |____\ circle (i : point of left bottom)
  function circle_list_get_c9p(this, i) result(r)
    class(circle_list), intent(in) :: this
    integer, intent(in)               :: i
    integer, dimension(3)             :: r
    integer site_r(3)
    integer r1(3), r2(3)
    site_r = site_to_r(i,this%pshape)
    !r(1) = this%cpl%seek2(path(i, i+1))
    !r(2) = this%cpl%seek2(path(i+1, i+this%nx))
    !r(3) = this%cpl%seek2(path(i+this%nx, i))
    r1 = site_r
    r2 = site_r + [1,0,0]
    r(1) = this%cpl%seek2(path(r_to_site(r1, this%pshape), r_to_site(r2, this%pshape)))
    r1 = r2
    r2 = site_r + [0,1,0]
    r(2) = this%cpl%seek2(path(r_to_site(r1, this%pshape), r_to_site(r2, this%pshape)))
    r1 = r2
    r2 = site_r + [0,0,0]
    r(3) = this%cpl%seek2(path(r_to_site(r1, this%pshape), r_to_site(r2, this%pshape)))
  end function circle_list_get_c9p

  ! circle of path_x + path_z
  function circle_list_get_c10p(this, i) result(r)
    class(circle_list), intent(in) :: this
    integer, intent(in)               :: i
    integer, dimension(4)             :: r
    integer site_r(3)
    integer r1(3), r2(3)
    site_r = site_to_r(i,this%pshape)
    !r(1) = this%cpl%seek2(path(r_to_site(site_r + [0,0,0], this%pshape), r_to_site(site_r + [1,0,0], this%pshape)))
    !r(2) = this%cpl%seek2(path(r_to_site(site_r + [1,0,0], this%pshape), r_to_site(site_r + [1,0,1], this%pshape)))
    !r(3) = this%cpl%seek2(path(r_to_site(site_r + [1,0,1], this%pshape), r_to_site(site_r + [0,0,1], this%pshape)))
    !r(4) = this%cpl%seek2(path(r_to_site(site_r + [0,0,1], this%pshape), r_to_site(site_r + [0,0,0], this%pshape)))
    r1 = site_r
    r2 = site_r + [1,0,0]
    r(1) = this%cpl%seek2(path(r_to_site(r1, this%pshape), r_to_site(r2, this%pshape)))
    r1 = r2
    r2 = site_r + [1,0,1]
    r(2) = this%cpl%seek2(path(r_to_site(r1, this%pshape), r_to_site(r2, this%pshape)))
    r1 = r2
    r2 = site_r + [0,0,1]
    r(3) = this%cpl%seek2(path(r_to_site(r1, this%pshape), r_to_site(r2, this%pshape)))
    r1 = r2
    r2 = site_r
    r(4) = this%cpl%seek2(path(r_to_site(r1, this%pshape), r_to_site(r2, this%pshape)))
  end function circle_list_get_c10p

  ! circle of path_x + path_z
  function circle_list_get_c11p(this, i) result(r)
    class(circle_list), intent(in) :: this
    integer, intent(in)               :: i
    integer, dimension(4)             :: r
    integer site_r(3)
    integer r1(3), r2(3)
    site_r = site_to_r(i,this%pshape)
    !r(1) = this%cpl%seek2(path(r_to_site(site_r + [0,0,0], this%pshape), r_to_site(site_r + [0,1,0], this%pshape)))
    !r(2) = this%cpl%seek2(path(r_to_site(site_r + [0,1,0], this%pshape), r_to_site(site_r + [0,1,1], this%pshape)))
    !r(3) = this%cpl%seek2(path(r_to_site(site_r + [0,1,1], this%pshape), r_to_site(site_r + [0,0,1], this%pshape)))
    !r(4) = this%cpl%seek2(path(r_to_site(site_r + [0,0,1], this%pshape), r_to_site(site_r + [0,0,0], this%pshape)))
    r1 = site_r
    r2 = site_r + [0,1,0]
    r(1) = this%cpl%seek2(path(r_to_site(r1, this%pshape), r_to_site(r2, this%pshape)))
    r1 = r2
    r2 = site_r + [0,1,1]
    r(2) = this%cpl%seek2(path(r_to_site(r1, this%pshape), r_to_site(r2, this%pshape)))
    r1 = r2
    r2 = site_r + [0,0,1]
    r(3) = this%cpl%seek2(path(r_to_site(r1, this%pshape), r_to_site(r2, this%pshape)))
    r1 = r2
    r2 = site_r
    r(4) = this%cpl%seek2(path(r_to_site(r1, this%pshape), r_to_site(r2, this%pshape)))
  end function circle_list_get_c11p

  subroutine circle_list_set_c1il(this)
    class(circle_list), intent(inout)      :: this
    !integer,dimension(:), allocatable         :: cnode
    integer                                   :: j
    type(int_list), dimension(:), allocatable :: temp
    integer, dimension(4) :: circle
    !logical  flag
    if (allocated(this%c1il)) deallocate(this%c1il)
    !cnode = get_c1_node(this%pshape)
    !do j = 1 , size(cnode)
    do j = 1 , this%n
      !circle = this%get_c1p(cnode(j))
      circle = this%get_c1p(j)
      !flag = .true.
      !do h = 1, this%nh
      !  !do k = 1, 4
      !  do k = 1, 4, 2
      !    if(this%hole(h,1) == this%cpl%value(abs(circle(k)))%i) flag = .false.
      !    if(this%hole(h,1) == this%cpl%value(abs(circle(k)))%f) flag = .false.
      !  end do
      !end do
      !print *, circle
      if(this%judge_circle(circle)) then
        allocate(temp(size(this%c1il)+1), &
        & source = [this%c1il, new_int_list(circle)])
        call move_alloc(temp, this%c1il)
      end if
    end do
  end subroutine circle_list_set_c1il

  !subroutine circle_list_set_c2il(this)
  !  class(circle_list), intent(inout)      :: this
  !  type(cnode_iterator), allocatable         :: iter
  !  integer                                   :: j
  !  type(int_list), dimension(:), allocatable :: temp
  !  if (allocated(this%c2il)) deallocate(this%c2il)
  !  allocate(iter, source = new_cnode_iterator(this%pshape))
  !  call iter%set_c2node(this%hole(:, 1))
  !  do while(iter%has_next())
  !    call iter%next(j)
  !    allocate(temp(size(this%c2il)+1), &
  !    & source = [this%c2il, new_int_list(this%get_c2p(j))])
  !    call move_alloc(temp, this%c2il)
  !  end do
  !end subroutine circle_list_set_c2il

  subroutine circle_list_set_c3il(this)
    class(circle_list), intent(inout)      :: this
    !integer,dimension(:), allocatable         :: cnode
    integer                                   :: j
    type(int_list), dimension(:), allocatable :: temp
    integer, dimension(3) :: circle
    if (allocated(this%c3il)) deallocate(this%c3il)
    do j = 1, this%nh
    circle = this%get_c3p(this%hole(j,1))
      allocate(temp(size(this%c3il)+1), &
      & source = [this%c3il, new_int_list(circle)])
      call move_alloc(temp, this%c3il)
    end do
  end subroutine circle_list_set_c3il

  subroutine circle_list_set_c4il(this)
    class(circle_list), intent(inout)      :: this
    !integer,dimension(:), allocatable         :: cnode
    integer                                   :: j
    type(int_list), dimension(:), allocatable :: temp
    integer, dimension(3) :: circle
    if (allocated(this%c4il)) deallocate(this%c4il)
    do j = 1, this%nh
    circle = this%get_c4p(this%hole(j,1))
      allocate(temp(size(this%c4il)+1), &
      & source = [this%c4il, new_int_list(circle)])
      call move_alloc(temp, this%c4il)
    end do
  end subroutine circle_list_set_c4il

  subroutine circle_list_set_c5il(this)
    class(circle_list), intent(inout)      :: this
    !integer,dimension(:), allocatable         :: cnode
    integer                                   :: j
    type(int_list), dimension(:), allocatable :: temp
    integer, dimension(3) :: circle
    if (allocated(this%c5il)) deallocate(this%c5il)
    do j = 1, this%nh
    circle = this%get_c5p(this%hole(j,1))
      allocate(temp(size(this%c5il)+1), &
      & source = [this%c5il, new_int_list(circle)])
      call move_alloc(temp, this%c5il)
    end do
  end subroutine circle_list_set_c5il

  subroutine circle_list_set_c6il(this)
    class(circle_list), intent(inout)      :: this
    !integer,dimension(:), allocatable         :: cnode
    integer                                   :: j
    type(int_list), dimension(:), allocatable :: temp
    integer, dimension(3) :: circle
    if (allocated(this%c6il)) deallocate(this%c6il)
    do j = 1, this%nh
    circle = this%get_c6p(this%hole(j,1))
      allocate(temp(size(this%c6il)+1), &
      & source = [this%c6il, new_int_list(circle)])
      call move_alloc(temp, this%c6il)
    end do
  end subroutine circle_list_set_c6il

  subroutine circle_list_set_c7il(this)
    class(circle_list), intent(inout)      :: this
    !integer,dimension(:), allocatable         :: cnode
    integer                                   :: j
    type(int_list), dimension(:), allocatable :: temp
    integer, dimension(4) :: circle
    if (allocated(this%c7il)) deallocate(this%c7il)
    do j = 1, this%nh
    circle = this%get_c7p(this%hole(j,1))
      allocate(temp(size(this%c7il)+1), &
      & source = [this%c7il, new_int_list(circle)])
      call move_alloc(temp, this%c7il)
    end do
  end subroutine circle_list_set_c7il

  subroutine circle_list_set_c8il(this)
    class(circle_list), intent(inout)      :: this
    integer                                   :: i, j
    type(int_list), dimension(:), allocatable :: temp
    type(int_list), allocatable :: node
    integer, dimension(3) :: circle
    if (allocated(this%c8il)) deallocate(this%c8il)
    node = this%get_node_c_2nd()
    do i = 1 , node%length()
      j = node%value(i)
      circle = this%get_c8p(j)
      if(this%judge_circle(circle)) then
        allocate(temp(size(this%c8il)+1), &
        & source = [this%c8il, new_int_list(circle)])
        call move_alloc(temp, this%c8il)
      end if
    end do
  end subroutine circle_list_set_c8il
  
  subroutine circle_list_set_c9il(this)
    class(circle_list), intent(inout)      :: this
    integer                                   :: i,j
    type(int_list), dimension(:), allocatable :: temp
    type(int_list), allocatable :: node
    integer, dimension(3) :: circle
    if (allocated(this%c9il)) deallocate(this%c9il)
    node = this%get_node_c_2nd()
    do i = 1 , node%length()
      j = node%value(i)
      circle = this%get_c9p(j)
      if(this%judge_circle(circle)) then
        allocate(temp(size(this%c9il)+1), &
        & source = [this%c9il, new_int_list(circle)])
        call move_alloc(temp, this%c9il)
      end if
    end do
  end subroutine circle_list_set_c9il
  
  subroutine circle_list_set_c10il(this)
    class(circle_list), intent(inout)      :: this
    integer,dimension(:), allocatable         :: cnode
    integer                                   :: j
    type(int_list), dimension(:), allocatable :: temp
    integer, dimension(4) :: circle
    logical  flag
    if (allocated(this%c10il)) deallocate(this%c10il)
    !cnode = get_c10_node(this%pshape)
    do j = 1 , this%n
      circle = this%get_c10p(j)

      if(this%judge_circle(circle)) then
        allocate(temp(size(this%c10il)+1), &
        & source = [this%c10il, new_int_list(circle)])
        call move_alloc(temp, this%c10il)
      end if
    end do
  end subroutine circle_list_set_c10il

  subroutine circle_list_set_c11il(this)
    class(circle_list), intent(inout)      :: this
    integer,dimension(:), allocatable         :: cnode
    integer                                   :: i,j
    type(int_list), dimension(:), allocatable :: temp
    type(int_list), allocatable :: node
    integer, dimension(4) :: circle
    logical  flag
    if (allocated(this%c11il)) deallocate(this%c11il)
    node = this%get_node_c11()
    do i = 1 , node%length()
      j = node%value(i)
      circle = this%get_c11p(j)
      if(this%judge_circle(circle)) then
        allocate(temp(size(this%c11il)+1), &
        & source = [this%c11il, new_int_list(circle)])
        call move_alloc(temp, this%c11il)
      end if
    end do
  end subroutine circle_list_set_c11il

  !subroutine circle_list_set_cil(this, rs)
  subroutine circle_list_set_cil(this)
    class(circle_list), intent(inout)      :: this
    !character(1),optional                   :: rs
    if (allocated(this%cil)) deallocate(this%cil)
    !if(present(rs))then
    !    allocate(this%cil(this%nc))
    !    !this%cil = [this%c1il, this%c2il, this%c3il, this%c4il, this%c5il, this%c6il]
    !    this%cil = [this%c1il, this%c3il, this%c4il, this%c5il, this%c6il, this%c7il]
    !else
    !    allocate(this%cil(this%nc))
    !    this%cil = [this%c1il, this%c2il]
    !endif
    allocate(this%cil(this%nc))
    this%c_1st_il = [this%c1il]
    this%c_rh_il = [this%c3il, this%c4il, this%c5il, this%c6il, this%c7il]
    this%c_2nd_il = [this%c8il, this%c9il]
    this%c_z_il = [this%c10il, this%c11il]
    this%cil = [this%c_1st_il, this%c_rh_il, this%c_2nd_il, this%c_z_il]
    !this%cil = [this%c_1st_il, this%c_rh_il, this%c_z_il]
  end subroutine circle_list_set_cil


  !end cil init
  function circle_list_get_delta_phase(this, phase) result(r)
    class(circle_list), intent(in)   :: this
    real(8), dimension(this%n), intent(in) :: phase
    real(8), dimension(this%np)            :: r
    integer                                :: i
    type(path)                             :: p
    do i = 1, this%np
      p = this%cpl%value(i)
      r(i) = phase(p%f)-phase(p%i)
    end do
  end function circle_list_get_delta_phase

  ! if cil(l) contain j or -j, return 1 or -1
  integer function circle_list_check_cil(this, l, j) result(r)!
    class(circle_list), intent(in) :: this
    integer, intent(in)          :: l, j
    if ((l < 1 .or. this%nc < l) .or. (j < 1 .or. this%np < j)) stop
    r = sign2(1, this%cil(l)%seek2(j))
  end function circle_list_check_cil

  !function circle_list_get_winding_number(this, phase, hole_num) result(r)!
  !  class(circle_list), intent(in)           :: this
  !  real(8), dimension(this%n), intent(in) :: phase
  !  integer, dimension(:), intent(in)      :: hole_num
  !  real(8), dimension(size(hole_num))     :: r
  !  integer                                :: i, j, a, x, y
  !  integer, dimension(9)                  :: s
  !  r = -8.d0 * pi
  !  do j = 1, size(r)
  !    a = this%hole(hole_num(j),1)
  !    x = this%x(a)
  !    y = this%y(a)
  !    s(1) = this%site(x-1,y-1)
  !    s(2) = this%site(x  ,y-1)
  !    s(3) = this%site(x+1,y-1)
  !    s(4) = this%site(x+1,y  )
  !    s(5) = this%site(x+1,y+1)
  !    s(6) = this%site(x  ,y+1)
  !    s(7) = this%site(x-1,y+1)
  !    s(8) = this%site(x-1,y  )
  !    s(9) = s(1)
  !    do i = 1, 8
  !      r(j) = r(j) + principal3(phase(s(i+1)) - phase(s(i)))
  !    end do
  !    r(j) = r(j) / (2.d0 * pi)
  !  end do
  !end function circle_list_get_winding_number

  function circle_list_judge_circle(this, circle) result(r)
    class(circle_list), intent(in)   :: this
    integer, dimension(:), intent(in) :: circle
    logical r
    integer i,j
    r = .false.
    do j = 1, size(circle)
        if(circle(j) == 0 ) return
    end do
    do i = 1, this%nh
      do j = 1, size(circle)
        if(this%hole(i,1) == this%cpl%value(abs(circle(j)))%i) return
        if(this%hole(i,1) == this%cpl%value(abs(circle(j)))%f) return
      end do
    end do
    r = .true.
    
  end function circle_list_judge_circle
  
  function circle_list_get_node_c_2nd(this) result(r)
    class(circle_list), intent(in)   :: this
    type(int_list), allocatable :: r
    integer i, j, k, l, site
    allocate(r)
    r = new_int_list()
    do k = 1, this%nz
      do j = 1, this%pshape(2,k) - 1
        do i = 1, this%pshape(1,k)-1
          site = xyz_to_site(i+1,j+1,k,this%pshape)
          if(this%ste(site) == 0) then
            cycle
          end if
          site = xyz_to_site(i+1,j,k,this%pshape)
          if(this%ste(site) == 0) then
            cycle
          end if
          site = xyz_to_site(i,j+1,k,this%pshape)
          if(this%ste(site) == 0) then
            cycle
          end if
          site = xyz_to_site(i,j,k,this%pshape)
          if(this%ste(site) /= 0) then
            call r.append(site)
          end if
        end do
      end do
    end do
  end function circle_list_get_node_c_2nd

  function circle_list_get_node_c11(this) result(r)
    class(circle_list), intent(in)   :: this
    type(int_list), allocatable :: r
    integer i, j, k, site, c_site

    allocate(r)
    r = new_int_list()

    do k = 1, this%nz-1
      do j = 1, this%pshape(2,k) - 1
        do i = 1, this%pshape(1,k)
          site = xyz_to_site(i,j,k,this%pshape)
          c_site = xyz_to_site(i-1,j+1,k,this%pshape)
          !print *, i-1,j+1,k,c_site
          !if(c_site < 0 .or. this%ste(c_site) == 0) then
          if(c_site < 0) then
            call r.append(site)
            cycle
          else if(this%ste(c_site) == 0) then
            call r.append(site)
            cycle
          end if
          c_site = xyz_to_site(i-1,j+1,k+1,this%pshape)
          !if(c_site < 0 .or. this%ste(c_site) == 0) then
          if(c_site < 0) then
            call r.append(site)
            cycle
          else if(this%ste(c_site) == 0) then
            call r.append(site)
            cycle
          end if

        end do
      end do
    end do

  end function circle_list_get_node_c11

  function circle_list_get_circle_center(this, i) result(r)
    class(circle_list), intent(in)   :: this
    integer, intent(in) :: i
    real(8), dimension(3) :: r
    integer :: j, k
    r = 0d0
    do j = 1, this%cil(i)%length()
      k = this%cil(i)%index(j)
      r = r + site_to_r(this%cpl%value(abs(k))%i, this%pshape)
      r = r + site_to_r(this%cpl%value(abs(k))%f, this%pshape)
    end do
    r = r * 0.5d0 / this%cil(i)%length()

  end function circle_list_get_circle_center
  
  
  ! constructor
  !function new_circle_list(pshape, hole) result(r)
  function new_circle_list(pshape, hole, path_2nd) result(r)
    integer, dimension(:,:), allocatable, intent(in)      :: pshape
    integer, dimension(:, :), intent(in)   :: hole
    logical, intent(in), optional   :: path_2nd
    type(circle_list)                        :: r
    !type(cuo2_plane), dimension(:), allocatable :: layers
    !character(1),optional                   :: rs
    !call r%circle_list_init(pshape, hole)
    call r%circle_list_init(pshape, hole, path_2nd)
  end function new_circle_list

end module circle_list_mod
