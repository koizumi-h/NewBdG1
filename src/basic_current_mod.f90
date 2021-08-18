module basic_current_mod
  use base_mod,           only : base
  use math_mod,           only : principal
  use path_list_mod,           only : path
  use hop_iterator_mod,   only : hop_iterator, new_hop_iterator
  use parameter_mod
  implicit none
  type, extends(base), public :: basic_current
    real(8), dimension(:, :), allocatable :: bond
  contains
    procedure :: set_bonds         => basic_current_set_bonds
    procedure :: get_bonds         => basic_current_get_bonds
    procedure :: calc              => basic_current_calc
    procedure :: draw_current      => basic_current_draw_current
    procedure :: draw_current2      => basic_current_draw_current2
    procedure :: draw_rs_current      => basic_current_draw_rs_current
    procedure :: draw_mean_current => basic_current_draw_mean_current
    procedure :: get_default_ratio => basic_current_get_default_ratio
  end type basic_current
contains
  ! public type-bound procedures
  subroutine basic_current_calc(this, chi)
    class(basic_current), intent(inout)         :: this
    real(8), dimension(this%n), intent(in) :: chi
    type(hop_iterator), allocatable        :: iter
    type(path)                             :: p
    if (allocated(this%bond)) deallocate(this%bond)
    allocate(this%bond(this%n, this%n), source = 0.d0)
    allocate(iter, source = new_hop_iterator(this%pshape))
    call iter%set_path([1,2,3,4,5,6], this%hole(:, 1))
    do while(iter%has_next())
      call iter%next(p)
      this%bond(p%f, p%i) = sin(principal(chi(p%f) - chi(p%i))*0.5d0)
    end do
    deallocate(iter)
!    this%bond(:, :) = -0.5d0 *this%bond(:, :)
    this%bond(:, :) = 0.5d0 *this%bond(:, :)
  end subroutine basic_current_calc

  pure subroutine basic_current_set_bonds(this, bond)
    class(basic_current), intent(inout)                 :: this
    real(8), dimension(this%n, this%n), intent(in) :: bond
    if (allocated(this%bond)) deallocate(this%bond)
    allocate(this%bond(this%n, this%n), source = bond)
  end subroutine basic_current_set_bonds

  pure function basic_current_get_bonds(this) result(r)
    class(basic_current), intent(in)        :: this
    real(8), dimension(this%n, this%n) :: r
    r(:, :) = 0.d0
    if (allocated(this%bond)) r(:, :) = this%bond(:, :)
  end function basic_current_get_bonds

  subroutine basic_current_draw_current(this, filename, magnitude)
    class(basic_current), intent(in)     :: this
    character(*), intent(in)        :: filename
    real(8), intent(in), optional   :: magnitude
    type(hop_iterator), allocatable :: iter
    type(path)                      :: p
    !real(8)                         :: ratio = 3.0d0, x, y, delta_x, delta_y
    real(8)                         :: ratio
    real(8),dimension(3)  ::  r, v, t
    integer i,j
    if (present(magnitude))then
      ratio = magnitude
    else
      ratio = 0d0
      !write(*,*) this%n
      do i = 1, this%n
        do j = 1, this%n
          !write(*,*) i,j,this%bond(i,j) 
          if( ratio < abs(this%bond(i,j))) then
            ratio = abs(this%bond(i,j))
          end if
        end do 
      end do
      ratio = 1d0 / ratio
      !write(*,*) "ratio", ratio ," = 1/ ",1d0/ratio 
    end if
    open(90, file = filename)
    !allocate(iter, source = new_hop_iterator(this%pshape))
    !call iter%set_path([1], this%hole(:, 1))
    !do while(iter%has_next())
    !  call iter%next(p)
    !  delta_x = this%bond(p%f, p%i) * ratio
    !  delta_y = 0.d0
    !  x = 0.5d0 * (dble(this%x(p%f)+this%x(p%i))-delta_x)
    !  y = dble(this%y(p%i))
    !  write(90, '(4(f19.15, 1x))') x, y, delta_x, delta_y
    !end do
    !call iter%set_path([2], this%hole(:, 1))
    !do while(iter%has_next())
    !  call iter%next(p)
    !  delta_x = 0.d0
    !  delta_y = this%bond(p%f, p%i) * ratio
    !  x = dble(this%x(p%i))
    !  y = 0.5d0 * (dble(this%y(p%f)+this%y(p%i))-delta_y)
    !  write(90, '(4(f19.15, 1x))') x, y, delta_x, delta_y
    !end do
    !deallocate(iter)
  
    allocate(iter, source = new_hop_iterator(this%pshape))
    !call iter%set_path([1,2,3,4,5,6], this%hole(:, 1))
    call iter%set_path([1,2,3,4,7], this%hole(:, 1))
    do while(iter%has_next())
      call iter%next(p)
      !write(*,*) p%i,p%f
      !if(p%i == 1 .and. p%f == 7) then
      t = (this%r(p%f) - this%r(p%i))
      v = ratio * this%bond(p%f, p%i) * t/(sqrt(t(1)**2 + t(2)**2 + t(3)**2))
      !v_x = 0.5d0
      !v_y = 0.5d0
      r = (this%r(p%f) + this%r(p%i))/2d0 - v / 2d0
      write(90, '(6(f19.15, 1x))') r, v
      !print *, sqrt(v(1)**2 + v(2)**2 + v(3)**2)
      !end if
    end do
    close(90)
  end subroutine basic_current_draw_current

  subroutine basic_current_draw_current2(this, filename)
    class(basic_current), intent(in)     :: this
    character(*), intent(in)        :: filename
    type(hop_iterator), allocatable :: iter
    type(path)                      :: p
    real(8)                         :: ratio
    real(8),dimension(3)  ::  r, v, t
    real(8)   ::  cl
    integer i,j
    ratio = 0d0
    do i = 1, this%n
      do j = 1, this%n
        if( ratio < abs(this%bond(i,j))) then
          ratio = abs(this%bond(i,j))
        end if
      end do 
    end do
    open(90, file = filename)
  
    allocate(iter, source = new_hop_iterator(this%pshape))
    call iter%set_path([1,2,3,4,7], this%hole(:, 1))
    do while(iter%has_next())
      call iter%next(p)
      t = (this%r(p%f) - this%r(p%i))
      v = 0.9d0*sign(1d0,this%bond(p%f,p%i))*t/(sqrt(t(1)**2 + t(2)**2 + t(3)**2))
      r = (this%r(p%f) + this%r(p%i))/2d0 - v / 2d0
      cl = abs(this%bond(p%f,p%i))!/ratio
      write(90, '(7(f19.15, 1x))') r, v, cl
    end do
    close(90)
  end subroutine basic_current_draw_current2

  subroutine basic_current_draw_rs_current(this, filename, magnitude)
    class(basic_current), intent(in)     :: this
    character(*), intent(in)        :: filename
    real(8), intent(in), optional   :: magnitude
    type(hop_iterator), allocatable :: iter!,tmp_iter
    type(path)                      :: p
    real(8)                         :: ratio = 3.0d0, x, y, delta_x, delta_y
    if (present(magnitude)) ratio = magnitude
    open(90, file = filename)
    allocate(iter, source = new_hop_iterator(this%pshape))
    !allocate(tmp_iter, source = new_hop_iterator(this%pshape))

    ! found error------
    call iter%set_path([1], this%hole(:, 1))
    ! found error------
    do while(iter%has_next())
      call iter%next(p)
      delta_x = this%bond(p%f, p%i) * ratio
      delta_y = 0.d0
      x = 0.5d0 * (dble(this%x(p%f)+this%x(p%i))-delta_x)
      y = dble(this%y(p%i))
      write(90, '(4(f19.15, 1x))') x, y, delta_x, delta_y
    end do

    call iter%set_path([2], this%hole(:, 1))
    do while(iter%has_next())
      call iter%next(p)
      delta_x = 0.d0
      delta_y = this%bond(p%f, p%i) * ratio
      x = dble(this%x(p%i))
      y = 0.5d0 * (dble(this%y(p%f)+this%y(p%i))-delta_y)
      write(90, '(4(f19.15, 1x))') x, y, delta_x, delta_y
    end do

    call iter%set_path([3], this%hole(:, 1))
    do while(iter%has_next())
      call iter%next(p)
      delta_x = this%bond(P%f,p%i)  * (1.d0/sqrt(2.d0))*ratio
      delta_y = this%bond(p%f, p%i) * (1.d0/sqrt(2.d0))*ratio
      x = dble(this%x(p%i))
      x = 0.5d0 * (dble(this%x(p%f)+this%x(p%i))-delta_x)
      y = 0.5d0 * (dble(this%y(p%f)+this%y(p%i))-delta_y)
      write(90, '(4(f19.15, 1x))') x, y, delta_x, delta_y
    end do

    call iter%set_path([4], this%hole(:, 1))
    do while(iter%has_next())
      call iter%next(p)
      delta_x = this%bond(p%f,p%i) * (-1.d0/sqrt(2.d0))*ratio
      delta_y = this%bond(p%f,p%i) * ( 1.d0/sqrt(2.d0))*ratio
      x = dble(this%x(p%i))
      x = 0.5d0 * (dble(this%x(p%f)+this%x(p%i))-delta_x)
      y = 0.5d0 * (dble(this%y(p%f)+this%y(p%i))-delta_y)
      write(90, '(4(f19.15, 1x))') x, y, delta_x, delta_y
    end do

    !call iter%set_path([5], this%hole(:, 1))
    !call tmp_iter%set_path([3,4], this%hole(:, 1))
    !call iter%delete_same_path(tmp_iter)
    !do while(iter%has_next())
    !  call iter%next(p)
    !  delta_x = this%bond(P%f,p%i)  * (-1.d0/sqrt(2.d0))*ratio
    !  delta_y = this%bond(p%f, p%i) * ( 1.d0/sqrt(2.d0))*ratio
    !  x = dble(this%x(p%i))
    !  x = 0.5d0 * (dble(this%x(p%f)+this%x(p%i))-delta_x)
    !  y = 0.5d0 * (dble(this%y(p%f)+this%y(p%i))-delta_y)
    !  write(90, '(4(f19.15, 1x))') x, y, delta_x, delta_y
    !end do

    deallocate(iter)
    close(90)
  end subroutine basic_current_draw_rs_current

  subroutine basic_current_draw_mean_current(this, filename, magnitude)
    class(basic_current), intent(in)     :: this
    character(*), intent(in)        :: filename
    real(8), intent(in), optional   :: magnitude
    type(hop_iterator), allocatable :: iter
    type(path)                      :: p
!    real(8)                         :: ratio = 2.0d0
    real(8)                         :: ratio = 2.5d0
    real(8), dimension(this%n)      :: mean_current_x, mean_current_y
    integer                         :: j
    if (present(magnitude)) ratio = magnitude
    mean_current_x = 0.0d0
    mean_current_y = 0.0d0
    allocate(iter, source = new_hop_iterator(this%pshape))
    call iter%set_path([1], this%hole(:, 1))
    do while(iter%has_next())
      call iter%next(p)
      mean_current_x(p%i) = mean_current_x(p%i) + this%bond(p%f, p%i)
      mean_current_x(p%f) = mean_current_x(p%f) - this%bond(p%i, p%f)
    end do
    call iter%set_path([2], this%hole(:, 1))
    do while(iter%has_next())
      call iter%next(p)
      mean_current_y(p%i) = mean_current_y(p%i) + this%bond(p%f, p%i)
      mean_current_y(p%f) = mean_current_y(p%f) - this%bond(p%i, p%f)
    end do
    deallocate(iter)
    mean_current_x = mean_current_x * ratio
    mean_current_y = mean_current_y * ratio
    open(82, file = filename)
    do j = 1, this%n
      write(82, '(4(f18.15, 1x))') &
      & this%x(j)-0.5d0*mean_current_x(j), &
      & this%y(j)-0.5d0*mean_current_y(j), &
      & mean_current_x(j), mean_current_y(j)
    end do
    close(82)
  end subroutine basic_current_draw_mean_current

  function basic_current_get_default_ratio(this) result(ratio)
    class(basic_current), intent(in)     :: this
    real(8)   ratio
    integer   i,j
    ratio = 0d0
    !write(*,*) this%n
    do i = 1, this%n
      do j = 1, this%n
        !write(*,*) i,j,this%bond(i,j) 
        if( ratio < abs(this%bond(i,j))) then
          ratio = abs(this%bond(i,j))
        end if
      end do 
    end do
    ratio = 1d0 / ratio
  end function
  
  
  ! constructor
  function new_basic_current(pshape, hole, bond) result(r)
    integer, dimension(:,:), intent(in)                 :: pshape
    integer, dimension(:, :), intent(in)              :: hole
    real(8), dimension(site_max(pshape), site_max(pshape)), &
    & intent(in), optional :: bond
    type(basic_current)                                    :: r
    call r%base_init_lattice(pshape)
    allocate(r%hole(size(hole, 1), size(hole, 2)), source = hole)
    if (present(bond)) call r%set_bonds(bond)
  end function new_basic_current
end module basic_current_mod
