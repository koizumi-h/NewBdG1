module conjugate_gradient_rsoc_mod
  USE nrtype
  use math_mod,  only : pi, ui
  use path_list_mod,  only : path 
  use angular_variable_mod, only : angular_variable
  use conjugate_gradient_mod, only : conjugate_gradient
  use parameter_mod
 implicit none
  !******* nr variables ********!
  !*****************************!
    type  :: near_hole
      type(path)    :: s(4)
    end type

  type, public, extends(conjugate_gradient)  :: conjugate_gradient_rsoc
    integer                             :: n, n_acc
    real(8),allocatable                :: lambda
    type(angular_variable),allocatable       :: an
    complex(8), dimension(:,:),allocatable   :: hfwf
    complex(8), dimension(:,:),allocatable   :: vv_ud, vv_du
    !real(8), dimension(:), allocatable       :: xi_n, chi_n
    real(8), dimension(:), allocatable       :: xi, chi
    type(near_hole), allocatable      :: nrhl(:)
  contains
    procedure :: conjugate_gradient_rsoc_init
    procedure :: calc           => conjugate_gradient_rsoc_calc
    procedure :: set_path       => conjugate_gradient_rsoc_set_path
    procedure :: set_wf       => conjugate_gradient_rsoc_set_wf
    procedure :: get_func       => conjugate_gradient_rsoc_get_func
    procedure :: get_dfunc      => conjugate_gradient_rsoc_get_dfunc
    procedure :: set_prod_vv => conjugate_gradient_set_prod_vv
  end type conjugate_gradient_rsoc
   
  !logical, parameter                  :: flag_use_rot_Hros = .false.

  !  type  :: near_hole
  !    type(path)    :: s(4)
  !  endtype
  !*****************************!
contains
  !*****************************!
  ! module subprograms
  ! public type-bound procedures
  subroutine conjugate_gradient_rsoc_init(this, n, n_acc, lambda, an, wf ,xi ,chi)
    class(conjugate_gradient_rsoc),intent(inout)  :: this
    real(8), intent(in)               :: lambda
    integer, intent(in)               :: n, n_acc
    type(angular_variable), intent(in)       :: an
    complex(8), dimension(4*n_acc, 4*n_acc), intent(in)         :: wf
    real(8), dimension(n), intent(in)    :: xi, chi ! index of xi, chi inputed is site number, not number excepted hole
    type(near_hole), allocatable      :: nrhl(:)
    integer i
    call this%conjugate_gradient_init()
    this%ftol = 1.0d-16
    this%n=n
    this%n_acc=n_acc
    if (allocated(this%lambda))deallocate(this%lambda)
    allocate(this%lambda)
    this%lambda=lambda
    if (allocated(this%an))deallocate(this%an)
    allocate(this%an)
    this%an = an
    if (allocated(this%xi))deallocate(this%xi)
    allocate(this%xi(n),source=xi)
    if (allocated(this%chi))deallocate(this%chi)
    allocate(this%chi(n),source=chi)
    !if (allocated(this%xi))deallocate(this%xi)
    !allocate(this%xi(n_acc))
    !if (allocated(this%chi))deallocate(this%chi)
    !allocate(this%chi(n_acc))
    !do i = 1, n
    !  this%xi( this%an%ste(i) ) = xi(i)
    !  this%chi( this%an%ste(i) ) = chi(i)
    !end do
    if (allocated(nrhl))deallocate(nrhl)
    allocate(nrhl(this%an%nh))
    call this%set_path(nrhl)
    if (allocated(this%nrhl))deallocate(this%nrhl)
    allocate(this%nrhl(this%an%nh))
    this%nrhl=nrhl
    if (allocated(this%hfwf))deallocate(this%hfwf)
    allocate(this%hfwf(4*n_acc,4*n_acc),source=wf)
    if (allocated(this%vv_ud))deallocate(this%vv_ud)
    allocate(this%vv_ud(4, this%an%nh))
    if (allocated(this%vv_du))deallocate(this%vv_du)
    allocate(this%vv_du(4, this%an%nh))
    call this%set_prod_vv()
  end subroutine conjugate_gradient_rsoc_init

  subroutine conjugate_gradient_rsoc_set_path(this, nrhl)
    class(conjugate_gradient_rsoc),intent(inout)  :: this
    type(near_hole),intent(out)     :: nrhl(this%an%nh)
    integer    i,j
    integer   z

      do i=1,this%an%nh
        z = this%an%z(this%an%hole(i,1))
        nrhl(i)%s(1)%i = this%an%ste( this%an%hole(i,1)-this%an%nx(z) )
        nrhl(i)%s(1)%f = this%an%ste( this%an%hole(i,1)+1 )
        nrhl(i)%s(2)%i = this%an%ste( this%an%hole(i,1)-1 )
        nrhl(i)%s(2)%f = this%an%ste( this%an%hole(i,1)+this%an%nx(z) )
        nrhl(i)%s(3)%i = this%an%ste( this%an%hole(i,1)-this%an%nx(z) )
        nrhl(i)%s(3)%f = this%an%ste( this%an%hole(i,1)-1 )
        nrhl(i)%s(4)%i = this%an%ste( this%an%hole(i,1)+1 )
        nrhl(i)%s(4)%f = this%an%ste( this%an%hole(i,1)+this%an%nx(z) )
      enddo
  end subroutine conjugate_gradient_rsoc_set_path

  subroutine conjugate_gradient_rsoc_set_wf(this, xi, wf)
    class(conjugate_gradient_rsoc), intent(in)        :: this
    integer                                    :: i, j
    real(8), dimension(this%n), intent(in)          :: xi 
    complex(8)                                 :: exi, echi
    complex(8), dimension(4*this%n_acc, 4*this%n_acc), intent(out)  :: wf
    do i = 1, this%n_acc
      exi  = exp(-ui*0.5d0* xi( this%an%ets(i) ))
      echi = exp(-ui*0.5d0* this%chi( this%an%ets(i) ))
      do j = 1, this%n_acc*4
        wf(4*i-3,j) = this%hfwf(4*i-3,j)*conjg(exi)*echi
        wf(4*i-2,j) = this%hfwf(4*i-2,j)*exi*echi
        wf(4*i-1,j) = this%hfwf(4*i-1,j)*exi*conjg(echi)
        wf(4*i  ,j) = this%hfwf(4*i  ,j)*conjg(exi)*conjg(echi)
      end do
    end do
  end subroutine conjugate_gradient_rsoc_set_wf

  subroutine conjugate_gradient_rsoc_calc(this, p, optimized)
    class(conjugate_gradient_rsoc), intent(inout)  :: this
    real(8),intent(out)     :: p(1)
    real(8)                 :: fret
    integer                 :: iter
    logical                 :: optimized
    call this%frprmn(p, this%ftol, iter, fret, optimized)
  end subroutine conjugate_gradient_rsoc_calc

  function conjugate_gradient_rsoc_get_func(this, val) result(r)
    class(conjugate_gradient_rsoc), intent(in)  :: this
    real(8), intent(in) :: val(:)
    real(8)                                    :: r
    complex(8)                                :: reg
    complex(8)                                 :: exi, echi
    complex(8), dimension(4*this%n_acc, 4*this%n_acc)        :: wf
    real(8), dimension(this%n)              :: xi
    integer                            :: i, j, k, l, a, js, ks
    complex(8)                         :: exichi_du, exichi_ud
    r = 0.d0
    reg=0.d0
    xi=this%an%rebuilt_angval_xi(this%xi, val(1))
  !  do i = 1, this%n
  !    exi  = exp(-ui*0.5d0* xi(i))
  !    echi = exp(-ui*0.5d0* this%chi(i))
  !    do j = 1, this%n*2
  !      wf(2*i-1,j) = this%hfwf(2*i-1,j)*exi*echi
  !      wf(2*i  ,j) = this%hfwf(2*i  ,j)*conjg(exi)*echi
  !    end do
  !  end do
   
   !  call this%set_wf(xi, wf)
   ! do i=1,this%an%nh
   !   do a = 2*this%n_acc + 1, 4*this%n_acc
   !     do l=1,2
   !       !k = this%nrhl(i)%s(l)%i
   !       !j = this%nrhl(i)%s(l)%f

   !      ! !reg = reg + exp( ui*0.25d0*pi)*conjg(wf(2*j,a)  )*wf(2*k-1,a) &
   !       !        & + exp( ui*0.75d0*pi)*conjg(wf(2*j-1,a))*wf(2*k,a)   &
   !       !        & + exp(-ui*0.25d0*pi)*conjg(wf(2*k-1,a))*wf(2*j,a)   &
   !       !        & + exp(-ui*0.75d0*pi)*conjg(wf(2*k,a)  )*wf(2*j-1,a)

   !       !reg=reg+conjg(wf(4*j,a))*wf(4*k-1,a)-conjg(wf(4*j-1,a))*wf(4*k,a)&
   !       !      &+conjg(wf(4*k-1,a))*wf(4*j,a)-conjg(wf(4*k,a))*wf(4*j-1,a)
   !       ! reg=reg+conjg(this%wf(2*j,a))*this%wf(2*k-1,a)-conjg(this%wf(2*j-1,a))*this%wf(2*k,a)&
   !       !     &+conjg(this%wf(2*k-1,a))*this%wf(2*j,a)-conjg(this%wf(2*k,a))*this%wf(2*j-1,a)
   !       
   !       j = this%nrhl(i)%s(l)%i
   !       k = this%nrhl(i)%s(l)%f
   !       reg = reg &
   !         & - wf(4*k  , a)*conjg(wf(4*j-1, a)) &
   !         & + wf(4*k-1, a)*conjg(wf(4*j  , a)) &
   !         & - wf(4*j-1, a)*conjg(wf(4*k  , a)) &
   !         & + wf(4*j  , a)*conjg(wf(4*k-1, a))
   !     enddo
   !     do l=3,4
   !       !k = this%nrhl(i)%s(l)%i
   !       !j = this%nrhl(i)%s(l)%f
   !       
   !       !reg = reg + exp( ui*0.75d0*pi)*conjg(wf(2*j,a)  )*wf(2*k-1,a) &
   !       !        & + exp( ui*0.25d0*pi)*conjg(wf(2*j-1,a))*wf(2*k,a)   &
   !       !        & + exp(-ui*0.75d0*pi)*conjg(wf(2*k-1,a))*wf(2*j,a)   &
   !       !        & + exp(-ui*0.25d0*pi)*conjg(wf(2*k,a)  )*wf(2*j-1,a)
   !       !reg=reg-ui*(conjg(wf(4*j,a))*wf(4*k-1,a)+conjg(wf(4*j-1,a))*wf(4*k,a))&
   !       !      &+ui*(conjg(wf(4*k-1,a))*wf(4*j,a)+conjg(wf(4*k,a))*wf(4*j-1,a))
   !       ! reg=reg+ui*(conjg(this%wf(2*j,a))*this%wf(2*k-1,a)+conjg(this%wf(2*j-1,a))*this%wf(2*k,a))&
   !       !       &-ui*(conjg(this%wf(2*k-1,a))*this%wf(2*j,a)+conjg(this%wf(2*k,a))*this%wf(2*j-1,a))
   !       
   !       j = this%nrhl(i)%s(l)%i
   !       k = this%nrhl(i)%s(l)%f
   !       reg = reg &
   !         & - ui*wf(4*k  , a)*conjg(wf(4*j-1, a)) &
   !         & - ui*wf(4*k-1, a)*conjg(wf(4*j  , a)) &
   !         & + ui*wf(4*j-1, a)*conjg(wf(4*k  , a)) &
   !         & + ui*wf(4*j  , a)*conjg(wf(4*k-1, a))
   !     enddo
   !   enddo
   ! enddo
    do i=1,this%an%nh
      do l=1,2
        j = this%nrhl(i)%s(l)%i
        k = this%nrhl(i)%s(l)%f
        js = this%an%ets(j)
        ks = this%an%ets(k)
        exichi_ud  = exp(-ui*0.5d0* (xi(ks) + xi(js) - this%chi(ks) + this%chi(js)) )
        exichi_du  = exp(-ui*0.5d0* (- xi(ks) - xi(js) - this%chi(ks) + this%chi(js)) )
        reg = reg &
          & - exichi_du*this%vv_du(l, i) &
          & + exichi_ud*this%vv_ud(l, i) &
          & - conjg(exichi_du*this%vv_du(l, i)) &
          & + conjg(exichi_ud*this%vv_ud(l, i))
      end do
      do l=3,4
        j = this%nrhl(i)%s(l)%i
        k = this%nrhl(i)%s(l)%f
        js = this%an%ets(j)
        ks = this%an%ets(k)
        exichi_ud  = exp(-ui*0.5d0* (xi(ks) + xi(js) - this%chi(ks) + this%chi(js)) )
        exichi_du  = exp(-ui*0.5d0* (- xi(ks) - xi(js) - this%chi(ks) + this%chi(js)) )
        reg = reg &
          & - ui*exichi_du*this%vv_du(l, i) &
          & - ui*exichi_ud*this%vv_ud(l, i) &
          & + ui*conjg(exichi_du*this%vv_du(l, i)) &
          & + ui*conjg(exichi_ud*this%vv_ud(l, i))
      end do
    end do
    r=this%lambda*real(reg, 8)
    !print *, "f : ", r
  end function conjugate_gradient_rsoc_get_func

  function conjugate_gradient_rsoc_get_dfunc(this, val) result(r)
    class(conjugate_gradient_rsoc), intent(in)  :: this
    real(8) , intent(in) :: val(:)
    real(8)      :: r(size(val))
    complex(8)                                :: reg
    complex(8)                                 :: exi, echi
    real(8), dimension(this%n)              :: xi
    complex(8), dimension(4*this%n_acc, 4*this%n_acc)        :: wf
    integer                            :: i, j, k, l, a, js, ks
    complex(8)                         :: exichi_du, exichi_ud
    r = 0.d0
    reg = 0.d0
    xi=this%an%rebuilt_angval_xi(this%xi, val(1))
 !   do i = 1, this%n
 !     exi  = exp(-ui*0.5d0* xi(i))
 !     echi = exp(-ui*0.5d0* this%chi(i))
 !     do j = 1, this%n*2
 !       wf(2*i-1,j) = this%hfwf(2*i-1,j)*exi*echi
 !       wf(2*i  ,j) = this%hfwf(2*i  ,j)*conjg(exi)*echi
 !     end do
 !   end do
    
    !call this%set_wf(xi, wf)
    !do i = 1, this%an%nh 
    !  do a = 2*this%n_acc + 1, 4*this%n_acc
    !    do l=1,2
    !    !  k = this%nrhl(i)%s(l)%i
    !    !  j = this%nrhl(i)%s(l)%f
    !      !if(flag_use_rot_Hros) then
    !      !  reg = reg - ui*exp(-ui*0.75d0*pi)*conjg(wf(2*j,a)  )*wf(2*k-1,a)  &
    !      !        &  + ui*exp( ui*0.75d0*pi)*conjg(wf(2*j-1,a))*wf(2*k,a)    &
    !      !        &  + ui*exp( ui*0.75d0*pi)*conjg(wf(2*k-1,a))*wf(2*j,a)    &
    !      !        &  - ui*exp(-ui*0.75d0*pi)*conjg(wf(2*k,a))*wf(2*j-1,a)
    !      !else
    !    !    reg=reg-ui*conjg(wf(4*j,a))*wf(4*k-1,a)-ui*conjg(wf(4*j-1,a))*wf(4*k,a)&
    !    !         &+ui*conjg(wf(4*k-1,a))*wf(4*j,a)-ui*conjg(wf(4*k,a))*wf(4*j-1,a)
    !    ! !  reg=reg-ui*conjg(this%wf(2*j,a))*this%wf(2*k-1,a)-ui*conjg(this%wf(2*j-1,a))*this%wf(2*k,a)&
    !    ! !        &+ui*conjg(this%wf(2*k-1,a))*this%wf(2*j,a)-ui*conjg(this%wf(2*k,a))*this%wf(2*j-1,a)
    !     ! end if
    !      
    !      j = this%nrhl(i)%s(l)%i
    !      k = this%nrhl(i)%s(l)%f
    !      reg = reg &
    !        & -   ui *wf(4*k  , a)*conjg(wf(4*j-1, a)) &
    !        & + (-ui)*wf(4*k-1, a)*conjg(wf(4*j  , a)) &
    !        & - (-ui)*wf(4*j-1, a)*conjg(wf(4*k  , a)) &
    !        & +   ui *wf(4*j  , a)*conjg(wf(4*k-1, a))
    !    enddo
    !    do l=3,4
    !    !  k = this%nrhl(i)%s(l)%i
    !    !  j = this%nrhl(i)%s(l)%f
    !      !if(flag_use_rot_Hros) then
    !      !  reg = reg + ui*exp(-ui*0.25d0*pi)*conjg(wf(2*j,a)  )*wf(2*k-1,a)   &
    !      !        & - ui*exp( ui*0.25d0*pi)*conjg(wf(2*j-1,a))*wf(2*k,a)     &
    !      !        & - ui*exp( ui*0.25d0*pi)*conjg(wf(2*k-1,a))*wf(2*j,a)     &
    !      !        & + ui*exp(-ui*0.25d0*pi)*conjg(wf(2*k,a)  )*wf(2*j-1,a)
    !      !else
    !    !    reg=reg-ui*(-ui*conjg(wf(4*j,a))*wf(4*k-1,a)+ui*conjg(wf(4*j-1,a))*wf(4*k,a))&
    !    !         &+ui*(ui*conjg(wf(4*k-1,a))*wf(4*j,a)-ui*conjg(wf(4*k,a))*wf(4*j-1,a))
    !    ! !  reg=reg+ui*(-ui*conjg(this%wf(2*j,a))*this%wf(2*k-1,a)+ui*conjg(this%wf(2*j-1,a))*this%wf(2*k,a))&
    !    ! !        &-ui*(ui*conjg(this%wf(2*k-1,a))*this%wf(2*j,a)-ui*conjg(this%wf(2*k,a))*this%wf(2*j-1,a))
    !      !end if
    !      
    !      j = this%nrhl(i)%s(l)%i
    !      k = this%nrhl(i)%s(l)%f
    !      reg = reg &
    !        & - ui*   ui*wf(4*k  , a)*conjg(wf(4*j-1, a)) &
    !        & - ui*(-ui)*wf(4*k-1, a)*conjg(wf(4*j  , a)) &
    !        & + ui*(-ui)*wf(4*j-1, a)*conjg(wf(4*k  , a)) &
    !        & + ui*   ui*wf(4*j  , a)*conjg(wf(4*k-1, a))
    !    enddo
    !  enddo
    !end do
    do i=1,this%an%nh
      do l=1,2
        j = this%nrhl(i)%s(l)%i
        k = this%nrhl(i)%s(l)%f
        js = this%an%ets(j)
        ks = this%an%ets(k)
        exichi_ud  = exp(-ui*0.5d0* (  xi(ks) + xi(js) - this%chi(ks) + this%chi(js)) )
        exichi_du  = exp(-ui*0.5d0* (- xi(ks) - xi(js) - this%chi(ks) + this%chi(js)) )
        reg = reg &
          & - ui*exichi_du*this%vv_du(l, i) &
          & + (-ui)*exichi_ud*this%vv_ud(l, i) &
          & - (-ui)*conjg(exichi_du*this%vv_du(l, i)) &
          & + ui*conjg(exichi_ud*this%vv_ud(l, i))
      end do
      do l=3,4
        j = this%nrhl(i)%s(l)%i
        k = this%nrhl(i)%s(l)%f
        js = this%an%ets(j)
        ks = this%an%ets(k)
        exichi_ud  = exp(-ui*0.5d0* (  xi(ks) + xi(js) - this%chi(ks) + this%chi(js)) )
        exichi_du  = exp(-ui*0.5d0* (- xi(ks) - xi(js) - this%chi(ks) + this%chi(js)) )
        reg = reg &
          & - ui*ui*exichi_du*this%vv_du(l, i) &
          & - ui*(-ui)*exichi_ud*this%vv_ud(l, i) &
          & + ui*(-ui)*conjg(exichi_du*this%vv_du(l, i)) &
          & + ui*ui*conjg(exichi_ud*this%vv_ud(l, i))
      end do
    end do
    r=this%lambda*real(reg, 8)
    !print *, "df : ", r
  end function conjugate_gradient_rsoc_get_dfunc
  
  subroutine conjugate_gradient_set_prod_vv(this)
    class(conjugate_gradient_rsoc), intent(inout)  :: this
    integer a, i, j, k, l
    this%vv_ud(:, :) = 0d0
    this%vv_du(:, :) = 0d0
    do i = 1, this%an%nh 
      do l=1,4
        j = this%nrhl(i)%s(l)%i
        k = this%nrhl(i)%s(l)%f
        do a = 2*this%n_acc + 1, 4*this%n_acc
          this%vv_du(l, i) = this%vv_du(l, i) + this%hfwf(4*k  , a)*conjg(this%hfwf(4*j-1, a))
          this%vv_ud(l, i) = this%vv_ud(l, i) + this%hfwf(4*k-1, a)*conjg(this%hfwf(4*j  , a))
        enddo
      enddo
    end do
  end subroutine conjugate_gradient_set_prod_vv

  function new_conjugate_gradient_rsoc(n, n_acc, lambda, an, wf, xi, chi) result(r)
    real(8), intent(in)               :: lambda
    integer, intent(in)               :: n, n_acc
    type(angular_variable), intent(in)       :: an
    complex(8), dimension(4*n_acc, 4*n_acc), intent(in)         :: wf
    real(8), dimension(n), intent(in)    :: xi, chi
    type(conjugate_gradient_rsoc)       :: r
    call r%conjugate_gradient_rsoc_init(n, n_acc, lambda, an, wf, xi, chi)
  end function new_conjugate_gradient_rsoc
end module conjugate_gradient_rsoc_mod
