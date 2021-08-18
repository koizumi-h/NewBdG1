module conjugate_gradient_spin_mod
  USE nrtype
  use math_mod,  only : ui
  use path_list_mod,  only : path 
  use angular_variable_mod, only : angular_variable
  use conjugate_gradient_mod, only : conjugate_gradient
  use hop_iterator_mod
 implicit none
  !******* nr variables ********!
  !*****************************!
    type  :: near_hole
      type(path)    :: s(4)
    endtype

    !!! this is a value for test
      integer cgspin_test_tmp_iter
    !!! end value
  type, public, extends(conjugate_gradient)  :: conjugate_gradient_spin
    integer                             :: n, nh, nz
    integer, allocatable                :: nl(:),hole_xi(:, :),pshape(:,:)
    real(8)                             :: lm, t(3), jd, u
    real(8), dimension(:), allocatable       :: xi
    !type(path), dimension(:), allocatable     :: pl_1st, pl_2nd, pl_h, pl_z
    type(path), dimension(:), allocatable     :: pl_1st, pl_h, pl_z
  contains
    procedure :: conjugate_gradient_spin_init
    procedure :: calc           => conjugate_gradient_spin_calc
    !procedure :: set_path       => conjugate_gradient_spin_set_path
    !procedure :: set_wf       => conjugate_gradient_spin_set_wf
    procedure :: get_func       => conjugate_gradient_spin_get_func
    procedure :: get_dfunc      => conjugate_gradient_spin_get_dfunc
    final :: conjugate_gradient_spin_finalize
  end type conjugate_gradient_spin
   
  !  type  :: near_hole
  !    type(path)    :: s(4)
  !  endtype
  !*****************************!
contains
  !*****************************!
  ! module subprograms
  ! public type-bound procedures
  subroutine conjugate_gradient_spin_init(this, pshape, hole_xi, t, u, jd, lm, xi)
    class(conjugate_gradient_spin),intent(inout)  :: this
    real(8), intent(in)               :: lm, t(3), u, jd
    integer, intent(in)               :: pshape(:,:), hole_xi(:,:)
    !type(angular_variable), intent(in)       :: an
    !complex(8), dimension(2*n, 2*n), intent(in)         :: wf
    real(8), dimension(:), intent(in)    :: xi!, chi
    !type(near_hole), allocatable      :: nrhl(:)
    type(hop_iterator), allocatable :: iter
    call this%conjugate_gradient_init()
    this%nz = size(pshape(1,:))
    if(allocated(this%nl)) deallocate(this%nl)
    allocate(this%nl(this%nz))
    this%nl = pshape(1,:) * pshape(2,:)
    this%n = sum(this%nl)
    this%nh = size(hole_xi,1)
    this%pshape = pshape
    allocate(this%hole_xi(this%nh,2), source=hole_xi)
    this%lm = lm
    this%jd = jd
    this%t = t
    this%u = u
    allocate(this%xi(this%n), source=xi)
    
    allocate(iter, source=new_hop_iterator(this%pshape))
    call iter%set_path([1,2], this%hole_xi(:, 1))
    call move_alloc(iter%pl%value(:), this%pl_1st)
    deallocate(iter)
    allocate(iter, source = new_hop_iterator(this%pshape))
    call iter%set_path_hole(this%hole_xi(:,1))
    call move_alloc(iter%pl%value(:), this%pl_h)
    deallocate(iter)
    !allocate(iter, source=new_hop_iterator(this%pshape))
    !call iter%set_path([5,6], this%hole_xi(:, 1))
    !call move_alloc(iter%pl%value(:), this%pl_2nd)
    !deallocate(iter)
    allocate(iter, source=new_hop_iterator(this%pshape))
    call iter%set_path([7], this%hole_xi(:, 1))
    call move_alloc(iter%pl%value(:), this%pl_z)
    deallocate(iter)
    !if (allocated(this%hfwf))deallocate(this%hfwf)
    !allocate(this%hfwf(2*n,2*n),source=wf)
    cgspin_test_tmp_iter = 1
  end subroutine conjugate_gradient_spin_init

  !subroutine conjugate_gradient_spin_set_path(this, nrhl)
  !  class(conjugate_gradient_spin),intent(inout)  :: this
  !  type(near_hole),intent(out)     :: nrhl(this%an%nh)
  !  integer    i,j

  !    do i=1,this%an%nh
  !      nrhl(i)%s(1)%i=this%an%hole(i,1)-this%an%nx
  !      nrhl(i)%s(1)%f=this%an%hole(i,1)+1
  !      nrhl(i)%s(2)%i=this%an%hole(i,1)-1
  !      nrhl(i)%s(2)%f=this%an%hole(i,1)+this%an%nx
  !      nrhl(i)%s(3)%i=this%an%hole(i,1)-this%an%nx
  !      nrhl(i)%s(3)%f=this%an%hole(i,1)-1
  !      nrhl(i)%s(4)%i=this%an%hole(i,1)+1
  !      nrhl(i)%s(4)%f=this%an%hole(i,1)+this%an%nx
  !    enddo
  !end subroutine conjugate_gradient_spin_set_path

  !subroutine conjugate_gradient_spin_set_wf(this, xi, wf)
  !  class(conjugate_gradient_spin), intent(inout)        :: this
  !  integer                                    :: i, j
  !  real(8), dimension(this%n), intent(in)          :: xi 
  !  complex(8)                                 :: exi, echi
  !  complex(8), dimension(2*this%n, 2*this%n), intent(out)  :: wf
  !  do i = 1, this%n
  !    exi  = exp(-ui*0.5d0* xi(i))
  !    echi = exp(-ui*0.5d0* this%chi(i))
  !    do j = 1, this%n*2
  !      wf(2*i-1,j) = this%hfwf(2*i-1,j)*exi*echi
  !      wf(2*i  ,j) = this%hfwf(2*i  ,j)*conjg(exi)*echi
  !    end do
  !  end do
  !end subroutine conjugate_gradient_spin_set_wf

  subroutine conjugate_gradient_spin_calc(this, xi)
    class(conjugate_gradient_spin), intent(inout)  :: this
    real(8),intent(inout)     :: xi(this%n)
    real(8)                 :: fret
    integer                 :: iter
    !xi=0.6d0
    call this%frprmn(xi, this%ftol, iter, fret)
  end subroutine conjugate_gradient_spin_calc

  function conjugate_gradient_spin_get_func(this, val) result(r)
    use hop_iterator_mod
    use path_list_mod
    use visual_mod
    class(conjugate_gradient_spin), intent(in)  :: this
    real(8), intent(in) :: val(:)
    real(8)            :: xi(this%n)
    !real(8)                                    :: r, E_nearest, E_h, E_2nd, E_z
    real(8)                                    :: r, E_nearest, E_h, E_z
    integer                            :: i, j, k

    !!! this is a value for test
      character(len=64) tmp
    !!! end value

    r = 0.d0
    E_nearest = 0d0
    E_h = 0d0
    !E_2nd = 0d0
    E_z = 0d0
    xi = val

      !!!!! test
      !if (.false.) then
      !  !write(tmp,"('a,i0,a')"), 'cg/out_', tmp_iter, '.gp'
      !  tmp = 'cg/out_' // string(cgspin_test_tmp_iter) // '.gp'
      !  call draw_phase_3D(tmp, 123 ,this%pshape,val,this%hole_xi(:,1))
      !  cgspin_test_tmp_iter = cgspin_test_tmp_iter + 1
      !  print *, tmp
      !end if
      !!! end test 
    
    ! The first nearest path
    ! E_nearest = (4t^2 / u) sum of inner product of S_j and S_k 
    ! over <k, j> 
    do i = 1, size(this%pl_1st)
      k = this%pl_1st(i)%i
      j = this%pl_1st(i)%f
      E_nearest = E_nearest &
        & + (1d0/4d0) * cos(xi(k) - xi(j))
    end do
    E_nearest = (4d0 * this%t(1)**2 / this%u) * E_nearest

    ! interaction berween holes accross holes
    do i = 1, size(this%pl_h)
      !write(*, *) k, j
      k = this%pl_h(i)%i
      j = this%pl_h(i)%f
      E_h = E_h &
        & + (1d0/4d0) * cos(xi(k) - xi(j))
    end do
    E_h = this%jd * E_h

    ! The 2nd nearest path
    !do i = 1, size(this%pl_2nd)
    !  k = this%pl_2nd(i)%i
    !  j = this%pl_2nd(i)%f
    !  E_2nd = E_2nd &
    !    & + (1d0/4d0) * cos(xi(k) - xi(j))
    !end do
    !!deallocate(iter)
    !E_2nd = (4d0 * this%t(2)**2 / this%u) * E_2nd

    ! z direction  path
    do i = 1, size(this%pl_z)
      !write(*, *) k, j
      k = this%pl_z(i)%i
      j = this%pl_z(i)%f
      E_z = E_z &
        & + (1d0/4d0) * cos(xi(k) - xi(j))
    end do
    !deallocate(iter)
    E_z = (4d0 * this%t(3)**2 / this%u) * E_z
    !print *,E_nearest , E_h , E_2nd , E_z
    !print *,this%t(1), this%jd, this%t(2), this%t(3)
    !r = E_nearest + E_h + E_2nd + E_z
    r = E_nearest + E_h + E_z
  end function conjugate_gradient_spin_get_func

  function conjugate_gradient_spin_get_dfunc(this, val) result(r)
    use hop_iterator_mod
    use path_list_mod
    class(conjugate_gradient_spin), intent(in)  :: this
    real(8) , intent(in) :: val(:)
    real(8)      :: r(size(val)), dE_nearest(size(val)), dE_h(size(val)),&
                !& dE_2nd(size(val)), dE_z(size(val))
                & dE_z(size(val))
    real(8), dimension(this%n)              :: xi
    integer                            :: i, j, k
    r = 0.d0
    dE_nearest = 0d0
    dE_h = 0d0
    !dE_2nd = 0d0
    dE_z = 0d0
    xi = val

    ! The first nearest path
    do i = 1, size(this%pl_1st)
      k = this%pl_1st(i)%i
      j = this%pl_1st(i)%f
      !write(*, *) k, j
      dE_nearest(k) = dE_nearest(k) &
        & - (1d0/4d0) * sin(xi(k) - xi(j)) 
      dE_nearest(j) = dE_nearest(j) &
        & - (1d0/4d0) * sin(xi(j) - xi(k)) 
    end do
    dE_nearest = (4d0 * this%t(1)**2 / this%u) * dE_nearest

    ! interaction berween holes accross holes
    do i = 1, size(this%pl_h)
      k = this%pl_h(i)%i
      j = this%pl_h(i)%f
      !write(*, *) i, k, j
      dE_h(k) = dE_h(k) &
        & - (1d0/4d0) * sin(xi(k) - xi(j)) 
      dE_h(j) = dE_h(j) &
        & - (1d0/4d0) * sin(xi(j) - xi(k)) 
    end do
    dE_h = this%jd * dE_h

    ! The 2nd nearest path
    !do i = 1, size(this%pl_2nd)
    !  k = this%pl_2nd(i)%i
    !  j = this%pl_2nd(i)%f
    !  !write(*, *) k, j
    !  dE_2nd(k) = dE_2nd(k) &
    !    & - (1d0/4d0) * sin(xi(k) - xi(j)) 
    !  dE_2nd(j) = dE_2nd(j) &
    !    & - (1d0/4d0) * sin(xi(j) - xi(k)) 
    !end do
    !dE_2nd = (4d0 * this%t(2)**2 / this%u) * dE_2nd


    ! z direction  path
    do i = 1, size(this%pl_z)
      k = this%pl_z(i)%i
      j = this%pl_z(i)%f
      !write(*, *) k, j
      dE_z(k) = dE_z(k) &
        & - (1d0/4d0) * sin(xi(k) - xi(j)) 
      dE_z(j) = dE_z(j) &
        & - (1d0/4d0) * sin(xi(j) - xi(k)) 
    end do
    dE_z = (4d0 * this%t(3)**2 / this%u) * dE_z
    
    !print *,E_nearest , E_h , E_2nd , E_z
    !r = dE_nearest + dE_h + dE_2nd + dE_z
    r = dE_nearest + dE_h + dE_z
  end function conjugate_gradient_spin_get_dfunc

  function new_conjugate_gradient_spin(pshape, hole_xi, t, u, jd, lm, xi) result(r)
    real(8), intent(in)               :: lm, t(3), u, jd
    integer, intent(in)               :: pshape(:,:), hole_xi(:, :)
    real(8), dimension(:), intent(in)    :: xi
    type(conjugate_gradient_spin)       :: r
    call r%conjugate_gradient_spin_init(pshape, hole_xi, t, u, jd, lm, xi)
  end function new_conjugate_gradient_spin

  subroutine conjugate_gradient_spin_finalize(this)
    type(conjugate_gradient_spin)       :: this
    if (allocated(this%pl_1st)) deallocate(this%pl_1st)
    !if (allocated(this%pl_2nd)) deallocate(this%pl_2nd)
    if (allocated(this%pl_h)) deallocate(this%pl_h)
    if (allocated(this%pl_z)) deallocate(this%pl_z)
  end subroutine conjugate_gradient_spin_finalize
end module conjugate_gradient_spin_mod
