module wave_function_mod
  use math_mod,          only : ui
  use path_list_mod,     only : path
  use hop_iterator_mod,  only : hop_iterator, new_hop_iterator
  use base_mod,          only : base
  use parameter_mod
  implicit none
  type, extends(base) :: wave_function
    real(8), dimension(:), allocatable       :: eg
    complex(8), dimension(:, :), allocatable :: wf1
    complex(8)                               :: teg
  contains
    procedure :: wave_function_init
    procedure :: wave_function_init_cp
    procedure :: current =>  wave_function_current
    procedure :: rs_current =>  wave_function_rs_current
    procedure :: set_wf =>  wave_function_set_wf
    procedure :: set_eg =>  wave_function_set_eg
    procedure :: total_energy => wave_function_total_energy
    procedure :: get_wf =>  wave_function_get_wf
    procedure :: get_eg =>  wave_function_get_eg
    procedure :: calc_delta =>  wave_function_calc_delta
    procedure :: calc_S =>  wave_function_calc_S
  end type wave_function
contains
!  subroutine wave_function_()
!    class(wave_function), intent(inout)        :: this
!  end subroutine wave_function_
!  function wave_function_() result(r)
!    class(wave_function), intent(in)        :: this
!  end function wave_function_

  subroutine wave_function_init(this, pshape, hole, ne, t, u, jd, lm, xi, chi, wf)
    class(wave_function), intent(inout)        :: this
    integer, dimension(:, :), intent(in)    :: pshape
    integer, dimension(:, :), intent(in) :: hole
    real(8), dimension(:), intent(in)          :: ne
    real(8), dimension(:), intent(in)          :: t
    real(8), intent(in)                        :: u,jd,lm
    real(8), dimension(:), intent(in)    :: xi, chi
    complex(8), dimension(:, :), intent(in),optional    :: wf
    call this%base_init(pshape, hole, ne, t, u, jd, lm, xi, chi, wf)
    this%teg=0.d0
    if(present(wf))then
      call this%set_wf(wf,chi,xi)
      call this%set_eg()
    endif
  end subroutine wave_function_init

  subroutine wave_function_init_cp(this, pshape, hole, ne, t, u, jd, lm, xi, chi, wf)
    class(wave_function), intent(inout)        :: this
    integer, dimension(:, :), intent(in)    :: pshape
    integer, dimension(:, :), intent(in) :: hole
    real(8), dimension(:), intent(in)          :: ne
    real(8), dimension(:), intent(in)          :: t
    real(8), intent(in)                        :: u,jd,lm
    real(8), dimension(:), intent(in)    :: xi, chi
    !real(8), dimension(pshape(1)*pshape(2))     :: temp
    !complex(8), dimension(2*pshape(1)*pshape(2), 2*pshape(1)*pshape(2)), &
    !& intent(in),optional                                  :: wf
    complex(8), dimension(:, :), intent(in), optional :: wf
    !real(8), dimension(:),allocatable     :: temp
    complex(8), dimension(:,:),allocatable     :: t_wf
    integer i, site, k,l
    call this%base_init(pshape, hole, ne, t, u, jd, lm, xi, chi, wf)
    this%teg=0.d0
    if(present(wf))then
      !allocate(temp(this%n))
      !temp=0.d0
      !call this%set_wf(wf,chi,temp) ! if (xi /= 0) then w -> -w
      !call this%set_wf(wf,chi,xi) ! if (xi /= 0) then w -> -w
      !call this%set_eg()
    else
      allocate(t_wf(4*this%n, 4*this%n))
      !t_wf = 0d0
      !do i = 1, sum(int(this%ne))
      !!  site = this%ets(i)
      !  site =i
      !  t_wf(4*site-3,i) = exp(-0.5d0*ui*this%xi(site))/sqrt(2.d0)
      !  t_wf(4*site-2,i) = exp(0.5d0*ui*this%xi(site))/sqrt(2.d0)
      !  t_wf(4*site-1,i) = exp(-0.5d0*ui*this%xi(site))/sqrt(2.d0)
      !  t_wf(4*site,i) = exp(0.5d0*ui*this%xi(site))/sqrt(2.d0)
      !end do
      t_wf = 0d0
      do i=1,4*this%n
        t_wf(i,i)=1d0
      end do
      this%wf = t_wf
      !call this%set_wf(t_wf,chi,xi)
      call this%set_eg()
    !do k=1,sum(int(this%ne))
    !  do l = 1, sum(int(this%ne))
    !    print *, k, l, DOT_PRODUCT(this%wf(:,l),this%wf(:,k))
    !  end do
    !end do
    !stop
    end if
  end subroutine wave_function_init_cp
  
  function wave_function_current(this) result(r)
    class(wave_function), intent(in)          :: this
    real(8), dimension(this%n,this%n)         :: r
    type(hop_iterator), allocatable           :: iter
    type(path)                                :: p
    integer                                   :: i, j, k, a
    r(:,:) = 0.d0
    allocate(iter, source = new_hop_iterator(this%pshape))
    call iter%set_path([1,2], this%hole(:,1))
    i=0
    do while(iter%has_next())
    i=i+1
      call iter%next(p)
      k = p%i
      j = p%f
      do a = 1, sum(int(this%ne))
        r(k,j) = r(k,j) -ui*(conjg(this%wf(4*k-1,a))*this%wf(4*j-1,a)  &
        &                   +conjg(this%wf(4*k  ,a))*this%wf(4*j  ,a)  &
        &                   -conjg(this%wf(4*j-1,a))*this%wf(4*k-1,a)  &
        &                   -conjg(this%wf(4*j  ,a))*this%wf(4*k  ,a)  &
        &                   +conjg(this%wf(4*k-1,a))*this%wf(4*j-1,a)  &
        &                   +conjg(this%wf(4*k  ,a))*this%wf(4*j  ,a)  &
        &                   -conjg(this%wf(4*j-1,a))*this%wf(4*k-1,a)  &
        &                   -conjg(this%wf(4*j  ,a))*this%wf(4*k  ,a))

      end do
      r(j,k) = -r(k,j)
   !   write(*,*)k,'to',j,r(j,k)
!      write(*,*)i,r(j,k)
    end do
    deallocate(iter)
 !   stop
  end function wave_function_current

  function wave_function_rs_current(this) result(r)
    class(wave_function), intent(in)          :: this
    real(8), dimension(this%n,this%n)         :: r
    type(hop_iterator), allocatable           :: iter
    type(path)                                :: p
    integer                                   :: i, j, k, a,np,npt
    r(:,:) = 0.d0
    allocate(iter, source = new_hop_iterator(this%pshape))
    call iter%set_path([1,2,3,4], this%hole(:,1))
    np = iter%pl%length()
    npt= np-4*this%nh
    do i = 1, npt
      k = iter%pl%value(i)%i
      j = iter%pl%value(i)%f
      do a = 1, sum(int(this%ne))
        r(k,j) = r(k,j) -ui*this%t(1)* &
        &                   (conjg(this%wf(4*k-1,a))*this%wf(4*j-1,a)  &
        &                   +conjg(this%wf(4*k  ,a))*this%wf(4*j  ,a)  &
        &                   -conjg(this%wf(4*j-1,a))*this%wf(4*k-1,a)  &
        &                   -conjg(this%wf(4*j  ,a))*this%wf(4*k  ,a)  &
        &                   +conjg(this%wf(4*k-1,a))*this%wf(4*j-1,a)  &
        &                   +conjg(this%wf(4*k  ,a))*this%wf(4*j  ,a)  &
        &                   -conjg(this%wf(4*j-1,a))*this%wf(4*k-1,a)  &
        &                   -conjg(this%wf(4*j  ,a))*this%wf(4*k  ,a))
     end do
      r(j,k) = -r(k,j)
   !   write(*,*)k,'to',j,r(j,k)
!      write(*,*)i,r(j,k)
    end do
    if (this%lm == 0.d0 )return
    do i = npt+1, np-(this%nh*2)
      p = iter%pl%index(i)
      k = p%f
      j = p%i
      do a = 1, sum(int(this%ne))
        r(j,k) = r(j,k) +ui*this%lm* &
              &      ( conjg(this%wf(2*k  ,a))*this%wf(2*j-1,a)  &
              &       -conjg(this%wf(2*k-1,a))*this%wf(2*j  ,a)  &
              &       -conjg(this%wf(2*j-1,a))*this%wf(2*k  ,a)  &
              &       +conjg(this%wf(2*j  ,a))*this%wf(2*k-1,a) )
      end do
      r(k,j) = -r(j,k)
    end do

    do i = np-(this%nh*2)+1,np
      p = iter%pl%index(i)
      k = p%f
      j = p%i
      do a = 1, sum(int(this%ne))
        r(j,k) = r(j,k) - this%lm* &
              &    ( conjg(this%wf(2*k  ,a))*this%wf(2*j-1,a)  &
              &     +conjg(this%wf(2*k-1,a))*this%wf(2*j  ,a)  &
              &     +conjg(this%wf(2*j-1,a))*this%wf(2*k  ,a)  &
              &     +conjg(this%wf(2*j  ,a))*this%wf(2*k-1,a)  )
      end do
      r(k,j) = -r(j,k)
    end do
    deallocate(iter)
 !   stop
  end function wave_function_rs_current

  subroutine wave_function_set_wf(this, wf, chi, xi)
    class(wave_function), intent(inout)        :: this
    complex(8), dimension(:,:), intent(in)     :: wf
    real(8), dimension(:), intent(in)          :: xi, chi
    integer                                    :: i, j
    complex(8)                                 :: exi, echi
    if (allocated(this%wf)) deallocate(this%wf)
    allocate(this%wf(this%n*2,this%n*2))
    allocate(this%wf1(this%n*2,this%n*2))
    if (allocated(this%xi)) deallocate(this%xi)
    allocate(this%xi(this%n))
    if (allocated(this%chi)) deallocate(this%chi)
    allocate(this%chi(this%n))
    this%wf1 = wf
    this%xi = xi
    this%chi = chi
    do i = 1, this%n
      exi  = exp(-ui*0.5d0* xi(i))
      echi = exp(-ui*0.5d0*this%chi(i))
      do j = 1, this%n*2
        this%wf(2*i-1,j) = wf(2*i-1,j)*exi*echi
        this%wf(2*i  ,j) = wf(2*i  ,j)*conjg(exi)*echi
      end do
    end do
  end subroutine wave_function_set_wf

  subroutine wave_function_set_eg(this)
    class(wave_function), intent(inout)        :: this
    integer                                    :: i, j, k, a, l, b
    real(8), dimension(this%n)                 :: density, Sx, Sy, Sz
    complex(8)                                 :: com,Hh,Hrsoc,Hrsoc_tot
    !complex(8), dimension(2*this%n)            :: ctemp
    real(8)                                    :: lambda !,rtemp
    !type(hop_iterator), allocatable           :: iter
    type(path)                               :: p
    type(path), allocatable                  :: p_rashba(:)
    real(8)                                   :: tmp_eg(4*this%n)
    integer z
    call this%calc_S(density, Sx, Sy, Sz)
    if(allocated(this%eg))deallocate(this%eg)
    !allocate(this%eg(sum(int(this%ne))))
    allocate(this%eg(4*this%n))
    Hrsoc_tot = 0.d0
    Hrsoc = 0.d0
    lambda=this%lm    ! coupling constant lambda=0.01
    !write(*,*) 'w_f_mod lm:', lambda
!    allocate(iter, source = new_hop_iterator(this%pshape))
    if(allocated(p_rashba)) deallocate(p_rashba)
    allocate(p_rashba(4*this%nh))
    do k=1,this%nh
      z = site_to_z(this%hole(k,1),this%pshape)
      p_rashba(1+(k-1)*4)%i=this%hole(k,1)-this%nx(z)
      p_rashba(1+(k-1)*4)%f=this%hole(k,1)+1
      p_rashba(2+(k-1)*4)%i=this%hole(k,1)-1
      p_rashba(2+(k-1)*4)%f=this%hole(k,1)+this%nx(z)
      p_rashba(3+(k-1)*4)%i=this%hole(k,1)-this%nx(z)
      p_rashba(3+(k-1)*4)%f=this%hole(k,1)-1
      p_rashba(4+(k-1)*4)%i=this%hole(k,1)+1
      p_rashba(4+(k-1)*4)%f=this%hole(k,1)+this%nx(z)
    enddo

    this%eg = 0.d0
    tmp_eg = 0d0
    !allocate(iter, source = new_hop_iterator(this%pshape))
    !call iter%set_path([1,2], this%hole(:,1))
    !do while(iter%has_next()) !NN hopping energy
    do i = 1, this%pl_1st%length()
      !call iter%next(p)
      p = this%pl_1st%value(i)
      k = p%i
      j = p%f
!      do a = 1, sum(int(this%ne))*2
      do a = 1, 4*this%n
        tmp_eg(a) = tmp_eg(a) + (conjg(this%wf(4*k-3,a))*this%wf(4*j-3,a) &
        &                       + conjg(this%wf(4*k-2 ,a))*this%wf(4*j-2,a) &
        &                       + conjg(this%wf(4*j-3,a))*this%wf(4*k-3,a) &
        &                       + conjg(this%wf(4*j-2,a))*this%wf(4*k-2,a) )&
        &                        -(conjg(this%wf(4*k-1,a))*this%wf(4*j-1,a) &
        &                       + conjg(this%wf(4*k  ,a))*this%wf(4*j  ,a) &
        &                       + conjg(this%wf(4*j-1,a))*this%wf(4*k-1,a) &
        &                       + conjg(this%wf(4*j  ,a))*this%wf(4*k  ,a) )
      end do
    end do
    this%eg(1:4*this%n) = this%eg(1:4*this%n) - tmp_eg * this%t(1) 
    !deallocate(iter)
    
    !tmp_eg = 0d0
    !allocate(iter, source = new_hop_iterator(this%pshape))
    !call iter%set_path([5,6], this%hole(:,1))
    !do while(iter%has_next()) !2nd hopping energy
    !!do i = 1, this%pl_2nd%length()
    !!  !call iter%next(p)
    !!  p = this%pl_2nd%value(i)
    !!  k = p%i
    !!  j = p%f
!   !!   do a = 1, sum(int(this%ne))*2
    !!  do a = 1, sum(int(this%ne))
    !!    tmp_eg(a) = tmp_eg(a) + (conjg(this%wf(2*k-1,a))*this%wf(2*j-1,a) &
    !!    &                       + conjg(this%wf(2*k  ,a))*this%wf(2*j  ,a) &
    !!    &                       + conjg(this%wf(2*j-1,a))*this%wf(2*k-1,a) &
    !!    &                       + conjg(this%wf(2*j  ,a))*this%wf(2*k  ,a) )
    !!  end do
    !!end do
    !!deallocate(iter)
    !this%eg(1:sum(int(this%ne))) = this%eg(1:sum(int(this%ne))) - tmp_eg * this%t(2)
     if(.false.)then
    tmp_eg = 0d0
    do i = 1, this%pl_z%length() ! z hop energy
      p = this%pl_z%value(i)
      k = p%i
      j = p%f
!      do a = 1, sum(int(this%ne))*2
      do a = 1, 4*this%n
        tmp_eg(a) = tmp_eg(a) + (conjg(this%wf(4*k-3,a))*this%wf(4*j-3,a) &
        &                       + conjg(this%wf(4*k-2,a))*this%wf(4*j-2,a) &
        &                       + conjg(this%wf(4*j-3,a))*this%wf(4*k-3,a) &
        &                       + conjg(this%wf(4*j-2,a))*this%wf(4*k-2,a) )&
        &                       -(conjg(this%wf(4*k-1,a))*this%wf(4*j-1,a) &
        &                       + conjg(this%wf(4*k  ,a))*this%wf(4*j  ,a) &
        &                       + conjg(this%wf(4*j-1,a))*this%wf(4*k-1,a) &
        &                       + conjg(this%wf(4*j  ,a))*this%wf(4*k  ,a) )
      end do
    end do
    !deallocate(iter)
    this%eg(1:4*this%n) = this%eg(1:4*this%n) - tmp_eg * this%t(3)


 
 !   write(*,*)sum(this%eg(1:sum(int(this%ne))))
!    do a = 1, sum(int(this%ne))*2 !Coulomb energy
    do a = 1, 4*this%n !Coulomb energy
      com = (0.d0, 0.d0)
      do j = 1, this%n
        com = com + (density(j)*0.5d0 - Sz(j)) * conjg(this%wf(4*j-3,a))*this%wf(4*j-3,a) &
        &         + (density(j)*0.5d0 + Sz(j)) * conjg(this%wf(4*j-2,a))*this%wf(4*j-2,a) &
        &         - (Sx(j) - ui *Sy(j))        * conjg(this%wf(4*j-3,a))*this%wf(4*j-2,a) &
        &         - (Sx(j) + ui *Sy(j))        * conjg(this%wf(4*j-2,a))*this%wf(4*j-3,a) &
        &         - (density(j)*0.5d0 - Sz(j)) * conjg(this%wf(4*j-1,a))*this%wf(4*j-1,a) &
        &         - (density(j)*0.5d0 + Sz(j)) * conjg(this%wf(4*j  ,a))*this%wf(4*j  ,a) &
        &         + (Sx(j) - ui *Sy(j))        * conjg(this%wf(4*j-1,a))*this%wf(4*j  ,a) &
        &         + (Sx(j) + ui *Sy(j))        * conjg(this%wf(4*j  ,a))*this%wf(4*j-1,a)
      end do
      this%eg(a) = this%eg(a) + this%u * com
    end do
    !tmp_eg = 0d0
    !do a = 1, sum(int(this%ne)) !Coulomb energy
    !  com = 0d0
    !  do j = 1, this%n
    !    com = com  &
    !    !tmp_eg(a) = tmp_eg(a)  &
    !         &   + (-(2d0/3d0)*Sz(j) + 0.5d0)   * CONJG(this%wf(2*j-1, a)) * this%wf(2*j-1, a) &
    !         &   + ( (2d0/3d0)*Sz(j) + 0.5d0)   * CONJG(this%wf(2*j  , a)) * this%wf(2*j  , a) &
    !         &   - (2d0/3d0)*(Sx(j) - ui*Sy(j)) * CONJG(this%wf(2*j-1, a)) * this%wf(2*j  , a) &
    !         &   - (2d0/3d0)*(Sx(j) + ui*Sy(j)) * CONJG(this%wf(2*j  , a)) * this%wf(2*j-1, a)
    !  end do
    !  this%eg(a) = this%eg(a) + this%u * com
    !end do
    
    !print *,Sy
    !print *,tmp_eg
    !stop
    !this%eg = this%eg + this%u * tmp_eg

    !allocate(iter, source = new_hop_iterator(this%pshape))
    !call iter%set_path_hole(this%hole(:,1))
!    do a=1 , sum(int(this%ne))*2 !making cp hamiltonian
    do a=1 , 4*this%n !making cp hamiltonian
      Hh = (0.d0, 0.d0) !Jh
      !do k=1,size(iter%pl%value(:))
      do k=1,this%pl_h%length()
        !i=iter%pl%value(k)%i
        !j=iter%pl%value(k)%f
        i=this%pl_h%value(k)%i
        j=this%pl_h%value(k)%f
       ! write(*,*)i,j
        Hh = Hh+ &
        &          (- Sz(j) * conjg(this%wf(4*i-2,a))*this%wf(4*i-2,a) &
        &           + Sz(j) * conjg(this%wf(4*i-3,a))*this%wf(4*i-3,a) &
        &           + (Sx(j)-ui*Sy(j)) * conjg(this%wf(4*i-3,a))*this%wf(4*i-2,a) &
        &           + (Sx(j)+ui*Sy(j)) * conjg(this%wf(4*i-2,a))*this%wf(4*i-3,a) &
        &           - Sz(i) * conjg(this%wf(4*j  ,a))*this%wf(4*j  ,a) &
        &           + Sz(i) * conjg(this%wf(4*j-1,a))*this%wf(4*j-1,a) &
        &           + (Sx(i)-ui*Sy(i)) * conjg(this%wf(4*j-1,a))*this%wf(4*j  ,a) &
        &           + (Sx(i)+ui*Sy(i)) * conjg(this%wf(4*j  ,a))*this%wf(4*j-1,a)  )
      enddo
      this%eg(a) = this%eg(a) + 0.5d0 *this%jd * Hh
    enddo

      !do a=1,sum(int(this%ne))
      !  do b=1,this%nh
      !  do l=1,2
      !    k = p_rashba(l+4*(b-1))%i
      !    j = p_rashba(l+4*(b-1))%f
      !    Hrsoc=Hrsoc+conjg(this%wf(2*j,a))*this%wf(2*k-1,a)-conjg(this%wf(2*j-1,a))*this%wf(2*k,a)&
      !        &+conjg(this%wf(2*k-1,a))*this%wf(2*j,a)-conjg(this%wf(2*k,a))*this%wf(2*j-1,a)
      !  enddo
      !  do l=3,4repare <i,j> for calc Jd
!  !test
!  !    do i=1, size(hop%pl%value(:))
!  !    write(*, *) 'car_parrinello_calc',i, hop%pl%value(i)
!  !    enddo
!  !    stop
!  CALL this%set_eg())+conjg(this%wf(2*k,a))*this%wf(2*j-1,a))
      !  enddo
      !  enddo
      !Hrsoc=this%lm*Hrsoc
      !Hrsoc_tot = Hrsoc_tot + Hrsoc
      !this%eg(a) = this%eg(a) + Hrsoc
!     ! write(*,*) 'a,eg(a)+Hrsoc:', a, this%eg(a)
      !enddo
      !write(*,*) 'Hrsoc_tot', Hrsoc_tot

 !   do i = 1, size(this%eg)-1
 !     do j = i, size(this%eg)
 !       if (this%eg(i) > this%eg(j)) then
 !         rtemp = this%eg(i)
 !         this%eg(i) = this%eg(j)
 !         this%eg(j) = rtemp
 !         ctemp = this%wf(:,i)
 !         this%wf(:,i) = this%wf(:,j)
 !         this%wf(:,j) = ctemp
 !       end if
 !     end do
 !   end do
  end if
  end subroutine wave_function_set_eg

  subroutine wave_function_total_energy(this)
    class(wave_function), intent(inout)        :: this
    integer                                    :: i,j,k,l,a,b
    complex(8)                                :: Hk,Hu,Hh,Hrsoc
    real(8)                                   :: eg(4*this%n),lambda
    type(hop_iterator), allocatable           :: iter
    type(path)                                :: p1
    type(path),allocatable   :: p(:)
    integer z
    eg = 0.d0
    K=0.d0
    Hu=0.d0
    Hh=0.d0
    Hrsoc=0.d0
      do l=1,this%nh
        this%wf(4*this%hole(l,1)-3,:)=0.d0
        this%wf(4*this%hole(l,1)-2,:)=0.d0
        this%wf(4*this%hole(l,1)-1,:)=0.d0
        this%wf(4*this%hole(l,1)  ,:)=0.d0
      enddo
    !expectation value of K
    allocate(iter, source = new_hop_iterator(this%pshape))
    call iter%set_path([1,2], this%hole(:,1))
    do while(iter%has_next())
      call iter%next(p1)
      k = p1%i
      j = p1%f
      do a = 1, this%n*4
        eg(a) = eg(a)  +(conjg(this%wf(4*k-3,a))*this%wf(4*j-3,a) &
        &              + conjg(this%wf(4*k-2,a))*this%wf(4*j-2,a) &
        &              - conjg(this%wf(4*j-1,a))*this%wf(4*k-1,a) &
        &              - conjg(this%wf(4*j  ,a))*this%wf(4*k  ,a) )
      end do
    end do
    deallocate(iter)
    Hk = -this%t(1)* sum(eg(1:sum(int(this%ne))))

   if(.false.)then
    !expectation value of Hu
    do a=1,sum(int(this%ne))
      do b=a+1,sum(int(this%ne))
        do j=1,this%n
          Hu=Hu+&
            &this%u*(-conjg(this%wf(4*j-2,a))  *conjg(this%wf(4*j-3,b))*this%wf(4*j-2,b)  *this%wf(4*j-3,a) &
            &        +conjg(this%wf(4*j-3,a))*conjg(this%wf(4*j-2,b))  *this%wf(4*j-2,b)  *this%wf(4*j-3,a) &
            &        +conjg(this%wf(4*j,a))  *conjg(this%wf(4*j-1,b))*this%wf(4*j-1,b)*this%wf(4*j,a) &
            &        -conjg(this%wf(4*j-1,a))*conjg(this%wf(4*j,b))  *this%wf(4*j-1,b)*this%wf(4*j,a))
        enddo
      enddo
    enddo

    !expectation value of Hh
    allocate(iter, source = new_hop_iterator(this%pshape))
    call iter%set_path_hole(this%hole(:,1))
    do a=1,sum(int(this%ne))
      do b=a+1,sum(int(this%ne))
        do k=1,size(iter%pl%value(:))
          i=iter%pl%value(k)%i
          j=iter%pl%value(k)%f
          Hh=Hh+&
            &0.5*this%jd*(-conjg(this%wf(4*j-3,a))*conjg(this%wf(4*i-2,b))  *this%wf(4*j-2,b)  *this%wf(4*i-3,a) &
            &        +conjg(this%wf(4*i-2,a))  *conjg(this%wf(4*j-3,b))*this%wf(4*j-2,b)  *this%wf(4*i-3,a) &
            &        +conjg(this%wf(4*j-3,a))*conjg(this%wf(4*i-2,b))  *this%wf(4*i-3,b)*this%wf(4*j-2,a) &
            &        -conjg(this%wf(4*i-2,a))  *conjg(this%wf(4*j-3,b))*this%wf(4*i-3,b)*this%wf(4*j-2,a) &
            &        -conjg(this%wf(4*j,a))  *conjg(this%wf(4*i-1,b))*this%wf(4*j-1,b)*this%wf(4*i,a) &
            &        +conjg(this%wf(4*i-1,a))*conjg(this%wf(4*j,b))  *this%wf(4*j-1,b)*this%wf(4*i,a) &
            &        +conjg(this%wf(4*j,a))  *conjg(this%wf(4*i-1,b))*this%wf(4*i,b)  *this%wf(4*j-1,a) &
            &        -conjg(this%wf(4*i-1,a))*conjg(this%wf(4*j,b))  *this%wf(4*i,b)  *this%wf(4*j-1,a))
        enddo
      enddo
    enddo
    deallocate(iter)

    !expectation value of Hrsoc
    lambda=0.01d0    ! coupling constant lambda=0.01
    allocate(iter, source = new_hop_iterator(this%pshape))
    if(allocated(p)) deallocate(p)
    allocate(p(4))
    do i=1,this%nh
      z = site_to_z(this%hole(i,1),this%pshape)
      p(1)%i=this%hole(i,1)-this%nx(z)
      p(1)%f=this%hole(i,1)+1
      p(2)%i=this%hole(i,1)-1
      p(2)%f=this%hole(i,1)+this%nx(z)
      p(3)%i=this%hole(i,1)-this%nx(z)
      p(3)%f=this%hole(i,1)-1
      p(4)%i=this%hole(i,1)+1
      p(4)%f=this%hole(i,1)+this%nx(z)
      do a=1,sum(int(this%ne))
        do l=1,2
          k = p(l)%i
          j = p(l)%f
          Hrsoc=Hrsoc+conjg(this%wf(4*j-2,a))*this%wf(4*k-3,a)-conjg(this%wf(4*j-3,a))*this%wf(4*k-2,a)&
              &+conjg(this%wf(4*k-3,a))*this%wf(4*j-2,a)-conjg(this%wf(4*k-2,a))*this%wf(4*j-3,a)&
              &-conjg(this%wf(4*j,a))*this%wf(4*k-1,a)+conjg(this%wf(4*j-1,a))*this%wf(4*k,a)&
              &-conjg(this%wf(4*k-1,a))*this%wf(4*j,a)+conjg(this%wf(4*k,a))*this%wf(4*j-1,a)
        enddo
        do l=3,4
          k = p(l)%i
          j = p(l)%f
          Hrsoc=Hrsoc+ui*(conjg(this%wf(4*j-2,a))*this%wf(4*k-3,a)+conjg(this%wf(4*j-3,a))*this%wf(4*k-2,a))&
            &-ui*(conjg(this%wf(4*k-3,a))*this%wf(4*j-2,a)+conjg(this%wf(4*k-2,a))*this%wf(4*j-3,a))&
            &-ui*(conjg(this%wf(4*j,a))*this%wf(4*k-1,a)+conjg(this%wf(4*j-1,a))*this%wf(4*k,a))&
            &-ui*(conjg(this%wf(4*k-1,a))*this%wf(4*j,a)+conjg(this%wf(4*k,a))*this%wf(4*j-1,a))
        enddo
      enddo
    enddo
    Hrsoc=this%lm*Hrsoc
    end if

    write(*,*)'K energy',Hk
    write(*,*)'Hu energy',Hu
    write(*,*)'Hh energy',Hh
    write(*,*)'Hrsoc energy',Hrsoc
    this%teg=Hk+Hu+Hh+Hrsoc
    write(*,*)'total energy',this%teg
    deallocate(iter)
  end subroutine wave_function_total_energy

  subroutine wave_function_calc_S(this, density, Sx, Sy, Sz)
    class(wave_function), intent(in)            :: this
    real(8), dimension(this%n), intent(out)     :: density, Sx, Sy, Sz
    integer                                       :: i, a
    complex(8)                                    :: uu,ud,vu, vd, Sp
    density = 0.d0
    Sx = 0.d0
    Sy = 0.d0
    Sz = 0.d0
    do i = 1, this%n
      Sp = (0.d0, 0.d0)
      do a = 2*this%n + 1, 4*this%n 
        uu = this%wf(4*i-3, a)
        ud = this%wf(4*i-2, a)
        vu = this%wf(4*i-1, a)
        vd = this%wf(4*i  , a)
        Sx(i) = Sx(i) + real(vu*conjg(vd))
        Sy(i) = Sy(i) - aimag(vu*conjg(vd))
        Sz(i) = Sz(i) + real(vu*conjg(vu)-vd*conjg(vd))/2
        density(i) = density(i) + conjg(uu)*uu + conjg(ud)*ud + conjg(vu)*vu + conjg(vd)*vd
      end do
    end do
  end subroutine wave_function_calc_S

  subroutine wave_function_calc_delta(this,density_up, density_down, delta)
    class(wave_function), intent(inout)            :: this
    real(8), dimension(this%n), intent(out)     :: density_up, density_down
    complex(8),dimension(size(this%pl_1st%value)),intent(out)            :: delta
    integer                                       :: i, j,a,k
    complex(8)                                    ::  uu, ud, vu, vd, Sp
    type(path)            :: p
    density_up   = 0.d0
    density_down = 0.d0
    delta        = 0.d0
    do i = 1, this%n
      do a = 1, 4*this%n
        uu = this%wf(4*i-3, a)
        ud = this%wf(4*i-2, a)
        vu = this%wf(4*i-1, a)
        vd = this%wf(4*i  , a)
        density_up(i) = density_up(i) + conjg(uu)*uu + conjg(vu)*vu
        density_down(i) = density_down(i) + conjg(ud)*ud + conjg(vd)*vd
        end do
      end do
    CALL this%set_eg()
    do a = 1, 4*this%n
      if(this%eg(a)>0)then
      do k = 1, this%pl_1st%length()
        p = this%pl_1st%value(k)
        i = p%i
        j = p%f
        delta(k) = this%jd*(this%wf(4*i-3, a)*conjg(this%wf(4*j, a))+this%wf(4*j-2, a)*conjg(this%wf(4*i-1, a))&
                &+this%wf(4*i-2, a)*conjg(this%wf(4*j-1, a))+this%wf(4*j-3, a)*conjg(this%wf(4*i  , a)))+delta(k)
      end do
      end if
     end do
  end subroutine wave_function_calc_delta

  function wave_function_get_wf(this) result(r)
    class(wave_function), intent(in)        :: this
    complex(8), dimension(this%n*4, this%n*4) :: r
    r = this%wf
  end function wave_function_get_wf

  function wave_function_get_eg(this) result(r)
    class(wave_function), intent(in)        :: this
    complex(8), dimension(this%n*4) :: r
    r = this%eg
  end function wave_function_get_eg

  ! constructor
  function new_wave_function(pshape, hole, ne, t, u, jd, lm,xi, chi, hfdwf) result(r)
    INTEGER, DIMENSION(:, :), INTENT(in)    :: pshape
    INTEGER, DIMENSION(:, :), INTENT(in) :: hole
    real(8), dimension(:), intent(in)    :: ne
    real(8), dimension(:), intent(in)    :: t
    real(8), intent(in)                  :: u,jd,lm
    real(8), dimension(:), intent(in)    :: xi, chi
    !complex(8), dimension(2*pshape(1)*pshape(2), 2*pshape(1)*pshape(2)), &
    !& intent(in)                                           :: hfdwf
    complex(8), dimension(:, :), intent(in),optional  :: hfdwf
    type(wave_function)               :: r
    call r%wave_function_init(pshape, hole, ne, t, u, jd, lm, xi, chi, hfdwf)
  end function new_wave_function
end module wave_function_mod