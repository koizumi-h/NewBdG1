module path_list_mod
  implicit none
  type, public :: path
    integer  :: i, f
  end type path

  type, public :: path_list
    type(path), dimension(:), allocatable :: value
  contains
    ! public type-bound procedures
    procedure :: path_list_init
    procedure :: flip          => path_list_flip
    procedure :: append        => path_list_append
    procedure :: extend        => path_list_extend
    procedure :: del           => path_list_del
    procedure :: index         => path_list_index
    procedure :: length        => path_list_length
    procedure :: pop           => path_list_pop
    procedure :: pop_sub       => path_list_pop_sub
    procedure :: alloc_value   => path_list_alloc_value
!    procedure :: realloc_value => path_list_realloc_value
    procedure :: display       => path_list_display
    procedure :: inverse       => path_list_inverse
    procedure :: remove        => path_list_remove
    procedure :: seek          => path_list_seek
    procedure :: seek2         => path_list_seek2
    procedure :: seek3         => path_list_seek3
    procedure :: seek4         => path_list_seek4
  end type path_list

contains
  ! public type-bound procedures
  pure subroutine path_list_init(this, path_array)
    class(path_list), intent(inout)                :: this
    type(path), dimension(:), intent(in), optional :: path_array
    call this%alloc_value(path_array)
  end subroutine path_list_init

   pure elemental  function path_list_flip(this, p) result(r)
    class(path_list), intent(in)                :: this
    type(path), intent(in)                     :: p
    type(path)                                 :: r
    r%i = p%f
    r%f = p%i
  end function path_list_flip

  pure subroutine path_list_append(this, p)
    class(path_list), intent(inout)       :: this
    type(path), intent(in)                :: p
    type(path), dimension(:), allocatable :: temp
    allocate(temp(this%length()+1), source = [this%value, p])
    deallocate(this%value)
    call move_alloc(temp, this%value)
  end subroutine path_list_append

  pure subroutine path_list_extend(this, that)
    class(path_list), intent(inout) :: this
    class(path_list), intent(in)    :: that
    type(path), dimension(:), allocatable :: temp
    if(that%length() /= 0) then
      allocate(temp(this%length()+that%length()), &
      & source = [this%value, that%value])
      deallocate(this%value)
      call move_alloc(temp, this%value)
    end if
  end subroutine path_list_extend

  pure subroutine path_list_del(this, m)
    class(path_list), intent(inout)       :: this
    integer, intent(in)                   :: m
    type(path), dimension(:), allocatable :: temp
    if (m < 1 .or. m > this%length()) return
    allocate(temp(this%length()-1), &
    & source = [this%value(1:m-1), this%value(m+1:this%length())])
    deallocate(this%value)
    call move_alloc(temp, this%value)
  end subroutine path_list_del

  pure type(path) function path_list_index(this, m) result(r)
    class(path_list), intent(in) :: this
    integer, intent(in)                :: m
    if (m < 1 .or. m > this%length()) then
      r = path(-1, -1)
    else
      r = this%value(m)
    end if
  end function path_list_index

  pure type(path) function path_list_index2(this, m) result(r)
    class(path_list), intent(in) :: this
    integer, intent(in)          :: m
    if (abs(m) < 1 .or. abs(m) > this%length()) then
      r = path(-1, -1)
    else
      if (m > 0) r = this%value(m)
!      if (m < 0) r = flip_path(this%value(abs(m)))
      if (m < 0) r = this%flip(this%value(abs(m)))
    end if
  end function path_list_index2

  pure integer function path_list_length(this) result(r)
    class(path_list), intent(in) :: this
    r = size(this%value)
  end function path_list_length

  function path_list_pop(this) result(r)
    class(path_list), intent(inout) :: this
    type(path)                      :: r
    if (this%length() == 0) then
      r = path(-1, -1)
    else
      r = this%value(1)
      call this%del(1)
    end if
  end function path_list_pop

  pure subroutine path_list_pop_sub(this, p)
    class(path_list), intent(inout) :: this
    type(path), intent(out)         :: p
    if (this%length() == 0) then
      p = path(-1, -1)
    else
      p = this%value(1)
      call this%del(1)
    end if
  end subroutine path_list_pop_sub

  pure subroutine path_list_alloc_value(this, path_array)
    class(path_list), intent(inout) :: this
    type(path), dimension(:), intent(in), optional :: path_array
    if (allocated(this%value)) deallocate(this%value)
    if (present(path_array)) then
      allocate(this%value(size(path_array)), source = path_array)
    else
      allocate(this%value(0))
    end if
  end subroutine path_list_alloc_value

  pure subroutine path_list_realloc_value(this, path_array)
    class(path_list), intent(inout) :: this
    type(path), dimension(:), intent(in), optional :: path_array
    if (allocated(this%value)) deallocate(this%value)
    call this%alloc_value(path_array)
  end subroutine path_list_realloc_value

  subroutine path_list_display(this)
    class(path_list), intent(in) :: this
    write(6, *) this%value(:)
  end subroutine path_list_display

  pure subroutine path_list_inverse(this)
    class(path_list), intent(inout) :: this
      this%value(:) = this%flip(this%value(:))
  end subroutine path_list_inverse

  pure subroutine path_list_remove(this, m)
    class(path_list), intent(inout) :: this
    integer, intent(in)             :: m
    integer                         :: i
    i = 1
    do while(i <= this%length())
      if (this%value(i)%i == m .or. this%value(i)%f == m) then
        call this%del(i)
        cycle
      end if
      i = i + 1
    end do
  end subroutine path_list_remove

  pure integer function path_list_seek(this, p) result(r)
    class(path_list), intent(in) :: this
    type(path), intent(in)       :: p
    do r = 1, this%length()
      if(this%value(r)%i == p%i .and. this%value(r)%f == p%f) return
    end do
    r = 0
  end function path_list_seek

  pure integer function path_list_seek2(this, p) result(r)
    class(path_list), intent(in) :: this
    type(path), intent(in)       :: p
    r = this%seek(p)
    if (r == 0) r = -this%seek(this%flip(p))
  end function path_list_seek2

  subroutine path_list_seek3(this, k, p)
    class(path_list),intent(inout) :: this
    integer,intent(in)             :: k
    integer    r
    type(path),allocatable,intent(out) :: p(:)
    type(path),allocatable :: tmp(:)
       allocate(p(0))
    do r=1,this%length()
      if(this%value(r)%i==k .or. this%value(r)%f==k)then
        allocate(tmp(size(p)+1),source=[p,this%value(r)])
        call move_alloc(tmp,p)
      endif
    enddo
  end subroutine path_list_seek3

  subroutine path_list_seek4(this, k, p )
    class(path_list),intent(inout) :: this
    integer,intent(in)             :: k
    integer    r
    type(path),allocatable,intent(out) :: p(:)
    type(path),allocatable :: tmp(:)
       allocate(p(0))
    do r=1,this%length()
      if(this%value(r)%i==k)then
        allocate(tmp(size(p)+1),source=[p,this%flip(this%value(r))])
        call move_alloc(tmp,p)
      elseif(this%value(r)%f==k)then
        allocate(tmp(size(p)+1),source=[p,this%value(r)])
        call move_alloc(tmp,p)
      endif
    enddo
  end subroutine path_list_seek4
  ! constructor
  pure type(path_list) function new_path_list(path_array) result(r)
    type(path), dimension(:), intent(in), optional :: path_array
    call r%path_list_init(path_array)
  end function new_path_list
end module path_list_mod
