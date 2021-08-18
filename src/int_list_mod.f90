module int_list_mod
  implicit none
  type, public :: int_list
    integer, dimension(:), allocatable :: value
  contains
    procedure :: int_list_init
    procedure :: append        => int_list_append
    procedure :: append2       => int_list_append2
    procedure :: extend        => int_list_extend
    procedure :: del           => int_list_del
    procedure :: index         => int_list_index
    procedure :: length        => int_list_length
    procedure :: pop           => int_list_pop
    procedure :: pop_sub       => int_list_pop_sub
    procedure :: alloc_value   => int_list_alloc_value
!    procedure :: realloc_value => int_list_realloc_value
    procedure :: display       => int_list_display
    procedure :: remove        => int_list_remove
    procedure :: seek          => int_list_seek
    procedure :: seek2         => int_list_seek2
  end type int_list
contains
  ! public type-bound procedures
  pure subroutine int_list_init(this, int_array)
    class(int_list), intent(inout)              :: this
    integer, dimension(:), intent(in), optional :: int_array
    call this%alloc_value(int_array)
  end subroutine int_list_init

  pure subroutine int_list_append(this, i)
    class(int_list), intent(inout)     :: this
    integer, intent(in)                :: i
    integer, dimension(:), allocatable :: temp
    allocate(temp(this%length()+1), source = [this%value, i])
    deallocate(this%value) !
    call move_alloc(temp, this%value)
  end subroutine int_list_append
  
  pure subroutine int_list_append2(this, i)
    class(int_list), intent(inout)     :: this
    integer, intent(in)                :: i
    integer, dimension(:), allocatable :: temp
    if (.not.allocated(this%value)) allocate(this%value(0))
    allocate(temp(this%length()+1), source = [this%value, i])
    deallocate(this%value) !
    call move_alloc(temp, this%value)
  end subroutine int_list_append2
  
  pure subroutine int_list_extend(this, that)
    class(int_list), intent(inout)     :: this
    class(int_list), intent(in)        :: that
    integer, dimension(:), allocatable :: temp
    allocate(temp(this%length()+that%length()), &
    & source = [this%value, that%value])
    deallocate(this%value)
    call move_alloc(temp, this%value)
  end subroutine int_list_extend
  
  pure subroutine int_list_del(this, m)
    class(int_list), intent(inout)     :: this
    integer, intent(in)                :: m
    integer, dimension(:), allocatable :: temp
    if (m < 1 .or. m > this%length()) return
    allocate(temp(this%length()-1), &
    & source = [this%value(1:m-1), this%value(m+1:this%length())])
    deallocate(this%value)
    call move_alloc(temp, this%value)
  end subroutine int_list_del
  
  pure integer function int_list_index(this, m) result(r)
    class(int_list), intent(in) :: this
    integer, intent(in)         :: m
    if (m < 1 .or. m > this%length()) then
      r = -1
    else
      r = this%value(m)
    end if
  end function int_list_index
  
  pure integer function int_list_length(this) result(r)
    class(int_list), intent(in) :: this
    r = size(this%value)
  end function int_list_length
  
  function int_list_pop(this) result(r)
    class(int_list), intent(inout) :: this
    integer                        :: r
    if (this%length() == 0) then
      r = -1
    else
      r = this%value(1)
      call this%del(1)
    end if
  end function int_list_pop
  
  pure subroutine int_list_pop_sub(this, p)
    class(int_list), intent(inout) :: this
    integer, intent(out)           :: p
    if (this%length() == 0) then
      p = -1
    else
      p = this%value(1)
      call this%del(1)
    end if
  end subroutine int_list_pop_sub
  
  pure subroutine int_list_alloc_value(this, int_array)
    class(int_list), intent(inout)              :: this
    integer, dimension(:), intent(in), optional :: int_array
    if (allocated(this%value)) deallocate(this%value)
    if (present(int_array)) then
      allocate(this%value(size(int_array)), source = int_array)
    else
      allocate(this%value(0))
    end if
  end subroutine int_list_alloc_value
  
  pure subroutine int_list_realloc_value(this, int_array)
    class(int_list), intent(inout)              :: this
    integer, dimension(:), intent(in), optional :: int_array
    if (allocated(this%value)) deallocate(this%value)
    call this%alloc_value(int_array)
  end subroutine int_list_realloc_value
  
  subroutine int_list_display(this)
    class(int_list), intent(in) :: this
    write(6, *) this%value(:)
  end subroutine int_list_display
  
  pure subroutine int_list_remove(this, m)
    class(int_list), intent(inout) :: this
    integer, intent(in)            :: m
    integer                        :: i
    i = 1
    do while(i <= this%length())
      if (this%value(i) == m) then
        call this%del(i)
        cycle
      end if
      i = i + 1
    end do
  end subroutine int_list_remove
  
  pure elemental integer function int_list_seek(this, i) result(r)
    class(int_list), intent(in) :: this
    integer, intent(in)         :: i
    integer                     :: j
    r = 0
    do j = 1, this%length()
      if (i == this%value(j)) then
        r = j
        exit
      end if
    end do
  end function int_list_seek
  
  pure elemental integer function int_list_seek2(this, i) result(r)
    class(int_list), intent(in) :: this
    integer, intent(in)         :: i
    integer                     :: j
    r = 0
    do j = 1, this%length()
      if (i == this%value(j)) then
        r = j
        exit
      else if (i == -this%value(j)) then
        r = -j
        exit
      end if
    end do
  end function int_list_seek2

  ! module subprograms
  subroutine copy_int_list(input, output)
    type(int_list), dimension(:), intent(in)    :: input
    type(int_list), dimension(:), allocatable, intent(out) :: output
    integer                                     :: i, j
    allocate(output(size(input)))
    do i = 1, size(output)
      do j = 1, input(i)%length()
        call output(i)%append2(input(i)%value(j))
      end do
    end do
  end subroutine copy_int_list
  
  ! constructor
  pure type(int_list) function new_int_list(int_array) result(r)
    integer, dimension(:), intent(in), optional :: int_array
    call r%int_list_init(int_array)
  end function new_int_list
end module int_list_mod
