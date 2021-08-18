module fchi_optimizer_Aem_mod
  use math_mod, only:pi, ui, principal, sign2, principal2, fermi
  use int_list_mod
  !use io_mod, only:save_txt, load_txt, create_vortex
  use path_list_mod, only:path, path_list
  use circle_list_mod, only:circle_list, new_circle_list
  use visual_mod
  use parameter_mod
  implicit none
  type   ::  near_path
    integer, allocatable  :: p(:)
  end type near_path
  
  !logical, parameter                  :: flag_use_rot_Hros = .false.

  type, extends(circle_list), public :: fchi_optimizer_Aem
    ! c1:small loop, c2:large loop, c:all loop
    real(8), dimension(:), allocatable        ::  delta_chi, jump, u_Aeff, Aem_dphase
    real(8)                            :: tol = 1.d-10, ratio = 0.2d0
    integer                            :: itermax = 20
    integer, dimension(:, :), allocatable     :: lcp
    type(near_path), allocatable               :: npath(:)
    logical                                   :: conv, rspath_flag
    real(8), dimension(:), allocatable              :: wconv
    integer, dimension(:), allocatable        :: w, w1, w2
    real(8), dimension(:), allocatable              :: jex_mat
    complex(8), allocatable                   :: EV(:)
    real(8), dimension(:), allocatable        :: ne_u, ne_d
    real(8), dimension(:), allocatable        :: eg
    !real(8)                                   :: kbt = 0.000001d0
    complex(8), dimension(:), allocatable        :: vv_uu, vv_dd, vv_ud, vv_du
    logical                                    :: pflag
    real(8)                                    :: Bz

  contains
    procedure :: fchi_optimizer_Aem_init
    !procedure :: draw_path    ! I made this subroutine just for debugging
    !procedure :: draw_int     ! for debugging
    procedure :: modify_init => fchi_optimizer_Aem_modify_init
    procedure :: load_cranking => fchi_optimizer_Aem_load_cranking
    procedure :: d_angval       => fchi_optimizer_Aem_d_angval
    procedure :: self_consistent => fchi_optimizer_Aem_self_consistent
    !procedure :: self_consistent2 => fchi_optimizer_Aem_self_consistent2
    procedure :: oneshot => fchi_optimizer_Aem_oneshot
    procedure :: get_Jp => fchi_optimizer_Aem_get_Jp
    !procedure :: get_Jp2 => fchi_optimizer_Aem_get_Jp2
    procedure :: get_dJdf => fchi_optimizer_Aem_get_dJdf
    !procedure :: get_dJdf2 => fchi_optimizer_Aem_get_dJdf2
    procedure :: get_n => fchi_optimizer_Aem_get_n
    procedure :: get_w => fchi_optimizer_Aem_get_w
    procedure :: set_lcp => fchi_optimizer_Aem_set_lcp
    procedure :: adjacent_path => fchi_optimizer_Aem_adjacent_path
    procedure :: calc_fdelta_chi => fchi_optimizer_Aem_calc_fdelta_chi
    procedure :: mix_delta_dchi => fchi_optimizer_Aem_mix_delta_dchi
    procedure :: set_jump => fchi_optimizer_Aem_set_jump
    procedure :: set_w1 => fchi_optimizer_Aem_set_w1
    procedure :: set_w2 => fchi_optimizer_Aem_set_w2
    procedure :: set_w => fchi_optimizer_Aem_set_w
    procedure :: set_param => fchi_optimizer_Aem_set_param
    procedure :: set_jex => fchi_optimizer_Aem_set_jex
    procedure :: set_prod_vv => fchi_optimizer_Aem_set_prod_vv
    procedure :: calc_only_Aem_dphase => fchi_optimizer_Aem_calc_only_Aem_dphase
    procedure :: set_vector_potential => fchi_optimizer_Aem_set_vector_potential
    procedure :: unset_vector_potential => fchi_optimizer_Aem_unset_vector_potential
    procedure :: calc_A_em => fchi_optimizer_Aem_calc_A_em
    procedure :: my_dgetrs
    procedure :: get_Jp_matrix => fchi_optimizer_Aem_get_Jp_matrix
    procedure :: approximated_J => fchi_optimizer_Aem_approximated_J
    procedure :: J0 => fchi_optimizer_Aem_J0
    procedure :: make_jmatrix => fchi_optimizer_Aem_make_jmatrix
  end type fchi_optimizer_Aem
contains

  subroutine fchi_optimizer_Aem_init(this, pshape, hole, ne, t, u, jd, lm, Bz, xi, chi, wf, eg, ne_u, ne_d)
    class(fchi_optimizer_Aem), intent(inout)                :: this
    integer, dimension(:, :), intent(in)                   :: pshape
    integer, dimension(:, :), intent(in)                :: hole
    real(8), dimension(:), intent(in)                   :: ne, t
    real(8), intent(in)                                 :: u, jd, lm, Bz
    real(8), dimension(site_max(pshape)), intent(in) :: xi, chi
    complex(8), dimension(:, :), intent(in)             :: wf
    real(8), dimension(:), intent(in)                   :: eg, ne_u, ne_d
    !type(cuo2_plane), dimension(:), allocatable         :: layers
    !character(1), optional                   :: rs
    !character(6)                             :: head = '2hole_'
    call this%base_init(pshape, hole, ne, t, u, jd, lm, xi, chi, wf)
    !call this%draw_path(this%cpl%value, head//'cpl')
    !call this%circle_list_init(pshape, hole, .false.) ! not 2nd hop
    call this%circle_list_init(pshape, hole) ! contain 2nd hop
    this%rspath_flag = .true.
    call this%set_lcp()
    call this%adjacent_path()
    call this%set_w1(this%nc1)
    call this%set_w2(hole(:, 2))
    call this%set_w(this%w1, this%w2)
    !allocate(this%wconv(this%np))
    allocate (this%wconv(this%n_acc - 1))
    !write(*,*) 'path, diffxi'
    !do i = 1, this%np
      !write(*,*) i, this%xi(this%cpl%value(i)%f) - this%xi(this%cpl%value(i)%i)
    !enddo
    if (allocated(this%delta_chi)) deallocate (this%delta_chi)
    this%delta_chi = this%d_angval(this%chi)
    call this%set_jump()
    if (allocated(this%jex_mat)) deallocate(this%jex_mat)
    allocate(this%jex_mat(this%np),source=0d0)
    if (allocated(this%ne_u)) deallocate(this%ne_u)
    allocate(this%ne_u(this%n),source=ne_u)
    if (allocated(this%ne_d)) deallocate(this%ne_d)
    allocate(this%ne_d(this%n),source=ne_d)
    if (allocated(this%eg)) deallocate(this%eg)
    allocate(this%eg(4*this%n_acc),source=eg)
    if (allocated(this%vv_uu)) deallocate(this%vv_uu)
    if (allocated(this%vv_dd)) deallocate(this%vv_dd)
    if (allocated(this%vv_ud)) deallocate(this%vv_ud)
    if (allocated(this%vv_du)) deallocate(this%vv_du)
    allocate(this%vv_uu(this%np))
    allocate(this%vv_dd(this%np))
    allocate(this%vv_ud(this%np))
    allocate(this%vv_du(this%np))
    call this%set_prod_vv()
    this%pflag = .true.
    this%Bz = Bz
    call this%calc_only_Aem_dphase()

    !write(*,*) 'np_1st', this%np_1st
    !write(*,*) 'np_rh', this%np_rh
    !write(*,*) 'np_2nd', this%np_2nd
    !write(*,*) 'np_z', this%np_z
    !! check whether the circles are constructed correctly
    !write(*,*)  "nc", this%nc
    !write(*,*)  "n_acc", this%n_acc
    !write(*,*)  "nc + n_acc - 1 : ", this%nc + this%n_acc -1
    !write(*,*)  "np : ", this%np
    if (this%nc + this%n_acc - 1 /= this%np) stop'error: nc+ne_acc-1/=np'
    !stop
  end subroutine fchi_optimizer_Aem_init


  function fchi_optimizer_Aem_d_angval(this, val) result(r)
    class(fchi_optimizer_Aem), intent(in)        :: this
    real(8), dimension(this%n), intent(in)     :: val
    real(8), dimension(this%np)                :: r
    integer                                    :: i
    type(path)                                 :: p
    do i = 1, this%np
      !write(*, *) i
      p = this%cpl%value(i)
      r(i) = modulo(val(p%f) - val(p%i), 2.d0*pi)
     if (abs(r(i)) > pi) r(i)=r(i)-sign(1.d0,r(i))*2.d0*pi
    end do
  end function fchi_optimizer_Aem_d_angval


  subroutine fchi_optimizer_Aem_modify_init(this, hole, xi, chi, wf)
    !------------
    ! this subroutine modifies fchi_optimizer_Aem_init but do not executes
    ! circle_list_mod we've allready executed
    !---------
    class(fchi_optimizer_Aem), intent(inout)            :: this
    integer, intent(in)                             :: hole(this%nh, 2)
    real(8), intent(in)                             :: xi(this%n), chi(this%n)
    complex(8), intent(in)                          :: wf(4*this%n_acc, 4*this%n_acc)
    call this%base_modify_init(hole, xi, chi, wf)
    call this%adjacent_path()
    call this%set_w1(this%nc1)
    call this%set_w2(hole(:, 2))
    call this%set_w(this%w1, this%w2)
    call this%set_jump()
    call this%set_prod_vv()
  end subroutine fchi_optimizer_Aem_modify_init


  subroutine fchi_optimizer_Aem_load_cranking(this)
    class(fchi_optimizer_Aem), intent(inout)            :: this
    integer                                         :: access
    if(access("crank_EV.fbin", " ") /= 0) then
      stop 'cannot find crank_EV.fbin, please run crank.exe'
    endif
    allocate(this%EV(21))
    open(31,file="crank_EV.fbin",status='old',form='unformatted')
    read(31) this%EV
    close(31)
  end subroutine fchi_optimizer_Aem_load_cranking

  subroutine fchi_optimizer_Aem_self_consistent(this, optimized)
    ! This subroutine is main program of fchi_optimizer_Aem,
    ! and calculate d_chi = chi_k - chi_j by the self-consistent culculation
    ! of current conservation.
    ! The chi is consist of the mult-valued function chi_0 and the single valued function fchi.
    class(fchi_optimizer_Aem), intent(inout)     :: this
    !real(8), dimension(this%np), &
    !& intent(in)                              :: init
    logical, intent(out), optional            :: optimized
    real(8), dimension(this%np)               :: dchi_0
    real(8), dimension(this%np)               :: df, dchi, dchi1, odf, odchi
    real(8), dimension(this%np)               :: u_Aeff  & ! u_Aeff = tau - intgral{A^{em} dr}
                                              &, o_u_Aeff  ! old of u_Aeff
    real(8), dimension(this%np)               :: Jp, dJdf
    real(8)                                   :: temp, temp1, val = 0d0, norm_odf, norm_df
    logical                                   :: flag
    integer                                   :: iter
    
    !call this%set_jump(init)
    !do iter = 1, this%np
    !write(*,*) 'i, dchi_0(i)', iter, dchi_0(iter)
    !enddo
    !do iter = 1, this%np
    !write(*,*) 'i, dxi(i)', iter, this%xi(this%cpl%value(iter)%f)-this%xi(this%cpl%value(iter)%i)
    !enddo
    !write(*, *) this%np
    !write(*, *) dchi_0
    
    dchi_0 = this%d_angval(this%chi)
    !write(*, *) dchi_02
    call this%set_jump()
    if (allocated(this%delta_chi)) deallocate(this%delta_chi)
    if (allocated(this%u_Aeff)) deallocate(this%u_Aeff)
    if (allocated(this%Aem_dphase)) deallocate(this%Aem_dphase)
    flag = .false.
    !dchi(:) = init(:)
    !odchi(:) = dchi_0(:)
    dchi(:) = dchi_0(:)
    
    !!!! u = tau - integral{ A^em dr }
    u_Aeff = this%set_vector_potential(dchi)
    
    !Jp = this%get_Jp(dchi)
    !dJdf = this%get_dJdf(dchi)
    
    Jp = this%get_Jp(u_Aeff)
    dJdf = this%get_dJdf(u_Aeff)

    ! for debug
    !print *,"---- wf ----"
    !print *, this%wf
    !print *,"---- wf ----"
    !print *,"---- dchi_0 ----"
    !print *,dchi_0
    !print *,"---- dchi_0 ----"
    !print *,"---- jp ----"
    !print *,jp
    !print *,"---- jp ----"
    !print *,"----- djdf ----"
    !print *,djdf
    !print *,"----- djdf ----"
    df = this%calc_fdelta_chi(Jp, dJdf)

    !for debug    
    !print *,"----- df --------"
    !print *,df
    !print *,"----- df --------"
    !stop

    odchi(:) = dchi(:)
    o_u_Aeff = u_Aeff
    odf = df
    norm_odf = dot_product(odf, odf)
    !write(*, *) 'init norm_odf', norm_odf
    dchi(:) = dchi(:) + df(:)
    u_Aeff(:) = u_Aeff(:) + df(:)

    !dchi = this%unset_vector_potential(u_Aeff)
    
    !odf(:) = 0.d0
    !odchi(:) = 0d0
    do iter = 1, this%itermax
      if(this%pflag) then
        write(6, '(a7, i4, a3, i7, a2)') '- iter ', iter, ' / ', this%itermax, ' -'
      end if
      dchi(:) = this%mix_delta_dchi(odchi, dchi, this%ratio)
      !Jp = this%get_Jp(dchi)
      !dJdf = this%get_dJdf(dchi)
      Jp = this%get_Jp(u_Aeff)
      dJdf = this%get_dJdf(u_Aeff)
      df = this%calc_fdelta_chi(Jp, dJdf)
      !odchi(:) = dchi(:)
      norm_df = sqrt(dot_product(df, df))
      norm_odf = sqrt(dot_product(odf, odf))
      if (norm_df < norm_odf) then
        dchi(:) = dchi(:) + df(:)
        u_Aeff(:) = u_Aeff(:) + df(:)
      else
        dchi(:) = dchi(:) + df(:) * (norm_odf/norm_df) * 0.8d0
        u_Aeff(:) = u_Aeff(:) + df(:) * (norm_odf/norm_df) * 0.8d0
      endif
      if(iter == 1) dchi1(:) = dchi(:)
      temp = maxval(abs(this%wconv))
      if(this%pflag) then
        write(6, *) 'flag', temp, norm_df, norm_df*(norm_odf/norm_df)*0.2d0
      end if
      !write(6, *) norm_odf, norm_df
      if (temp < this%tol) then
      !if (temp < this%tol .and. iter > 10) then
        !print *,temp
        flag = .true.
        exit
      end if
      !val = sum(odchi) - sum(dchi(:)) !
      odchi(:) = dchi(:) !
      o_u_Aeff = u_Aeff
      odf=df
      !write(*,*) 'odchi-dchi, df^2', val !
      !write(*,*) ''
      !write(*, *) maxval(abs(this%wconv))
    end do
      !write(*, *) 'd val', val
      !val = sum(odf) - sum(df) !
      !write(*, *) 'f val', val
      !write(*, *) 'flag', flag
    !allocate(this%delta_chi(this%np), source = dchi + df)
    this%conv = .false.
    allocate(this%delta_chi(this%np), source=dchi)
    allocate(this%u_Aeff(this%np), source=u_Aeff)
    
    allocate(this%Aem_dphase(this%np))
    df = 0d0
    this%Aem_dphase = this%set_vector_potential(df)
    
    if (flag) then
      !write (6, *) 'SCO-calculation has been converged!'
      if(this%pflag) then
        write (6, *) 'fchi calculation has been converged!'
      end if
      this%conv = .true.
      if(present(optimized)) optimized = .true.
    else
      ! If the self-consistent calculation is not converged,
      ! we use the dchi with first cycle
      !allocate (this%delta_chi(this%np), source=dchi)
      if(this%pflag) then
        write(6, '(a30)') '-------- ! warning ! ---------'
        !write (6, *) 'SCO-calculation is not converged.'
        write (6, *) 'fchi calculation is not converged.'
        !write(6, *) temp1, temp
      end if
      !if(temp>1) write(*,*) this%hole(:,2)
     if(present(optimized)) optimized = .false.
    end if
  end subroutine fchi_optimizer_Aem_self_consistent

  !subroutine fchi_optimizer_Aem_self_consistent2(this,EV, optimized)
  !  ! This subroutine is main program of fchi_optimizer_Aem,
  !  ! and calculate d_chi = chi_k - chi_j by the self-consistent culculation
  !  ! of current conservation.
  !  ! The chi is consist of the mult-valued function chi_0 and the single valued function fchi.
  !  class(fchi_optimizer_Aem), intent(inout)     :: this
  !  !real(8), dimension(this%np), &
  !  !& intent(in)                              :: init
  !  real(8), dimension(this%np)               :: dchi_0
  !  logical, intent(out)                      :: optimized
  !  real(8), dimension(this%np)               :: df, dchi, dchi1, odf, odchi
  !  real(8), dimension(this%np)               :: Jp, dJdf
  !  real(8)                                   :: temp, temp1, val = 0d0, norm_odf, norm_df
  !  complex(8), intent(in)                    :: EV
  !  logical                                   :: flag
  !  integer                                   :: iter
  !  !call this%set_jump(init)
  !  !do iter = 1, this%np
  !  !write(*,*) 'i, dchi_0(i)', iter, dchi_0(iter)
  !  !enddo
  !  !do iter = 1, this%np
  !  !write(*,*) 'i, dxi(i)', iter, this%xi(this%cpl%value(iter)%f)-this%xi(this%cpl%value(iter)%i)
  !  !enddo
  !  !write(*, *) this%np
  !  !write(*, *) dchi_0
  !  dchi_0 = this%d_angval(this%chi)
  !  !write(*, *) dchi_02
  !  call this%set_jump()
  !  if (allocated(this%delta_chi)) deallocate (this%delta_chi)
  !  flag = .false.
  !  !dchi(:) = init(:)
  !  !odchi(:) = dchi_0(:)
  !  dchi(:) = dchi_0(:)
  !  Jp = this%get_Jp2(dchi,EV)
  !  dJdf = this%get_dJdf2(dchi, EV)
  !  df = this%calc_fdelta_chi(Jp, dJdf)
  !  odchi(:) = dchi(:)
  !  odf = df
  !  norm_odf = dot_product(odf, odf)
  !  !write(*, *) 'init norm_odf', norm_odf
  !  dchi(:) = dchi(:) + df(:)

  !  !odf(:) = 0.d0
  !  !odchi(:) = 0d0
  !  do iter = 1, this%itermax
  !    !write(6, '(a7, i4, a3, i4, a2)') '- iter ', iter, ' / ', this%itermax, ' -'
  !    dchi(:) = this%mix_delta_dchi(odchi, dchi, this%ratio)
  !    Jp = this%get_Jp2(dchi, EV)
  !    dJdf = this%get_dJdf2(dchi, EV)
  !    df = this%calc_fdelta_chi(Jp, dJdf)
  !    !odchi(:) = dchi(:)
  !    norm_df = sqrt(dot_product(df, df))
  !    norm_odf = sqrt(dot_product(odf, odf))
  !    if (norm_df < norm_odf) then
  !      dchi(:) = dchi(:) + df(:)
  !    else
  !      dchi(:) = dchi(:) + df(:) * (norm_odf/norm_df) * 0.8d0
  !    endif
  !    if(iter == 1) dchi1(:) = dchi(:)
  !    temp = maxval(abs(this%wconv))
  !    !write(6, *) 'flag', temp, norm_df, norm_df*(norm_odf/norm_df)*0.2d0
  !    if (temp < this%tol) then
  !      flag = .true.
  !      exit
  !    end if
  !    !val = sum(odchi) - sum(dchi(:)) !
  !    odchi(:) = dchi(:) !
  !    odf=df
  !    !write(*,*) 'odchi-dchi, df^2', val !
  !    !write(*,*) ''
  !    !write(*, *) maxval(abs(this%wconv))
  !  end do
  !    !write(*, *) 'd val', val
  !    !val = sum(odf) - sum(df) !
  !    !write(*, *) 'f val', val
  !    !write(*, *) 'flag', flag
  !  !allocate(this%delta_chi(this%np), source = dchi + df)
  !  this%conv = .false.
  !  if (flag) then
  !    allocate (this%delta_chi(this%np), source=dchi)
  !    !write (6, *) 'SCO-calculation has been converged!'
  !    this%conv = .true.
  !    optimized = .true.
  !  else
  !    ! If the self-consistent calculation is not converged,
  !    ! we use the dchi with first cycle
  !    allocate (this%delta_chi(this%np), source=dchi)
  !    !write(6, '(a30)') '-------- ! warning ! ---------'
  !    !write (6, *) 'SCO-calculation is not converged.'
  !    !write(6, *) temp1, temp
  !    !if(temp>1) write(*,*) this%hole(:,2)
  !    optimized = .false.
  !  end if
  !end subroutine fchi_optimizer_Aem_self_consistent2


  !subroutine fchi_optimizer_Aem_self_consistent_oneshot(this, dchi_0, df_result)
  subroutine fchi_optimizer_Aem_oneshot(this)
    ! This subroutine is main program of fchi_optimizer_Aem,
    ! and calculate d_chi = chi_k - chi_j by the self-consistent culculation
    ! of current conservation.
    ! The chi is consist of the mult-valued function chi_0 and the single valued function fchi.
    class(fchi_optimizer_Aem), intent(inout)     :: this
    !real(8), dimension(this%np), &
    !& intent(in)                              :: init
    !real(8), dimension(this%np), intent(out)  :: df_result
    real(8), dimension(this%np)               :: df, dchi, u_Aeff
    real(8), dimension(this%np)               :: Jp, dJdf
    !call this%set_jump(init)
    !do iter = 1, this%np
    !write(*,*) 'i, dchi_0(i)', iter, dchi_0(iter)
    !enddo
    !do iter = 1, this%np
    !write(*,*) 'i, dxi(i)', iter, this%xi(this%cpl%value(iter)%f)-this%xi(this%cpl%value(iter)%i)
    !enddo
    call this%set_jump()
    if (allocated(this%delta_chi)) deallocate (this%delta_chi)
    if (allocated(this%u_Aeff)) deallocate(this%u_Aeff)
    if (allocated(this%Aem_dphase)) deallocate(this%Aem_dphase)
    !dchi(:) = init(:)
    !odchi(:) = dchi_0(:)
    !dchi(:) = dchi_0(:)
    dchi = this%d_angval(this%chi)
    u_Aeff = this%set_vector_potential(dchi)
    !Jp = this%get_Jp(dchi)
    !dJdf = this%get_dJdf(dchi)
    Jp = this%get_Jp(u_Aeff)
    dJdf = this%get_dJdf(u_Aeff)
    !df_result = this%calc_fdelta_chi(Jp, dJdf)
    df = this%calc_fdelta_chi(Jp, dJdf)
    dchi = dchi + df
    u_Aeff = u_Aeff + df
    !dchi = dchi - this%calc_fdelta_chi(Jp, dJdf)*this%ratio
    allocate(this%delta_chi(this%np), source=dchi)
    allocate(this%u_Aeff(this%np), source=u_Aeff)
    allocate(this%Aem_dphase(this%np))
    this%Aem_dphase = u_Aeff - dchi !this%set_vector_potential([0d0,....,0d0])
  end subroutine fchi_optimizer_Aem_oneshot


  function fchi_optimizer_Aem_get_Jp(this, dchi) result(r)
    class(fchi_optimizer_Aem), intent(in),target         :: this
    real(8), dimension(this%np), intent(in)     :: dchi
    real(8), dimension(this%np)                :: r
    !real(8), dimension(4*this%n_acc)                :: r_temp
    integer                                   :: i, j, k, a, npt, ie, fe, is, fs, iter_start
    complex(8)                                :: exichipm(this%np), exichimp(this%np), exichipp(this%np), exichimm(this%np)
    type path_pointer
      integer, pointer :: i, f
    end type path_pointer
    type(path_pointer)                       :: p(this%np)
    real(8)                :: r_temp1
    !   !DIR$ VECTOR ALWAYS
    do k = 1, this%np
      p(k)%i => this%cpl%value(k)%i
      p(k)%f => this%cpl%value(k)%f
    enddo
    !if (this%rspath_flag) then
    !  npt = this%np - 4*this%nh
    !else
    !  npt = this%np
    !endif
    !   !DIR$ VECTOR ALWAYS
    do k = 1, this%np
      ie = this%ste(p(k)%i)
      fe = this%ste(p(k)%f)
      is = p(k)%i
      fs = p(k)%f
      exichipm(k) = exp(0.5d0*ui*( this%xi(fs) - this%xi(is) + dchi(k) + this%jump(k)))
      exichimp(k) = exp(0.5d0*ui*(-this%xi(fs) + this%xi(is) + dchi(k) + this%jump(k)))
      exichipp(k) = exp(0.5d0*ui*( this%xi(fs) + this%xi(is) + dchi(k) + this%jump(k)))
      exichimm(k) = exp(0.5d0*ui*(-this%xi(fs) - this%xi(is) + dchi(k) + this%jump(k)))
    end do
    r = 0.d0
! path[1,2] the nearest path
    !do i = 1, npt
    do k = 1, this%np_1st
      ie = this%ste(p(k)%i)
      fe = this%ste(p(k)%f)
      is = p(k)%i
      fs = p(k)%f
      !r_temp(:) =  0d0
      if(this%z(is) == 1) then ! surface
        !   !DIR$ VECTOR ALWAYS
        !do a = 2*this%n_acc+1, 4*this%n_acc
        !  r_temp(a) =  - ui* ( &
        !  !&       (       exichimp(k) *conjg(this%wf(4*fe-3, a))*this%wf(4*ie-3, a)*(1d0-this%ne_d(is))*(1d0-this%ne_d(fs))  &
        !  !&       +       exichipm(k) *conjg(this%wf(4*fe-2, a))*this%wf(4*ie-2, a)*(1d0-this%ne_u(is))*(1d0-this%ne_u(fs))  &
        !  !&       - conjg(exichimp(k))*conjg(this%wf(4*ie-3, a))*this%wf(4*fe-3, a)*(1d0-this%ne_d(is))*(1d0-this%ne_d(fs))  &
        !  !&       - conjg(exichipm(k))*conjg(this%wf(4*ie-2, a))*this%wf(4*fe-2, a)*(1d0-this%ne_u(is))*(1d0-this%ne_u(fs))  &
        !  !&       )*fermi(this%eg(a), this%kbt) & 
        !  &     + (       exichimp(k) *this%wf(4*fe-1, a)*conjg(this%wf(4*ie-1, a))*(1d0-this%ne_d(is))*(1d0-this%ne_d(fs))  &
        !  &       +       exichipm(k) *this%wf(4*fe  , a)*conjg(this%wf(4*ie  , a))*(1d0-this%ne_u(is))*(1d0-this%ne_u(fs))  &
        !  &       - conjg(exichimp(k))*this%wf(4*ie-1, a)*conjg(this%wf(4*fe-1, a))*(1d0-this%ne_d(is))*(1d0-this%ne_d(fs))  &
        !  &       - conjg(exichipm(k))*this%wf(4*ie  , a)*conjg(this%wf(4*fe  , a))*(1d0-this%ne_u(is))*(1d0-this%ne_u(fs))  &
        !  !&       )*fermi(-this%eg(a), this%kbt) )
        !  &       ) )
        !enddo
        r_temp1 =  - ui* ( &
          &     + (       exichipm(k) *this%vv_uu(k)*(1d0-this%ne_d(is))*(1d0-this%ne_d(fs))  &
          &       +       exichimp(k) *this%vv_dd(k)*(1d0-this%ne_u(is))*(1d0-this%ne_u(fs))  &
          &       - conjg(exichipm(k))*conjg(this%vv_uu(k))*(1d0-this%ne_d(is))*(1d0-this%ne_d(fs))  &
          &       - conjg(exichimp(k))*conjg(this%vv_dd(k))*(1d0-this%ne_u(is))*(1d0-this%ne_u(fs)) ) )
          
        !r_temp(1) =  - ui* ( &
          !&     + (       exichimp(k) *this%vv_uu(k)*(1d0-this%ne_d(is))*(1d0-this%ne_d(fs))  &
          !&       +       exichipm(k) *this%vv_dd(k)*(1d0-this%ne_u(is))*(1d0-this%ne_u(fs))  &
          !&       - conjg(exichimp(k))*conjg(this%vv_uu(k))*(1d0-this%ne_d(is))*(1d0-this%ne_d(fs))  &
          !&       - conjg(exichipm(k))*conjg(this%vv_dd(k))*(1d0-this%ne_u(is))*(1d0-this%ne_u(fs)) ) )
      else    ! bulk
        !   !DIR$ VECTOR ALWAYS
        !do a = 2*this%n_acc+1, 4*this%n_acc
        !  r_temp(a) =  - ui* ( &
        !  !&       (       exichimp(k) *conjg(this%wf(4*fe-3, a))*this%wf(4*ie-3, a)  &
        !  !&       +       exichipm(k) *conjg(this%wf(4*fe-2, a))*this%wf(4*ie-2, a)  &
        !  !&       - conjg(exichimp(k))*conjg(this%wf(4*ie-3, a))*this%wf(4*fe-3, a)  &
        !  !&       - conjg(exichipm(k))*conjg(this%wf(4*ie-2, a))*this%wf(4*fe-2, a)  &
        !  !&       )*fermi(this%eg(a), this%kbt) & 
        !  &     + (       exichimp(k) *this%wf(4*fe-1, a)*conjg(this%wf(4*ie-1, a))  &
        !  &       +       exichipm(k) *this%wf(4*fe  , a)*conjg(this%wf(4*ie  , a))  &
        !  &       - conjg(exichimp(k))*this%wf(4*ie-1, a)*conjg(this%wf(4*fe-1, a))  &
        !  &       - conjg(exichipm(k))*this%wf(4*ie  , a)*conjg(this%wf(4*fe  , a))  &
        !  !&       )*fermi(-this%eg(a), this%kbt) )
        !  &       ) )
        !enddo
        r_temp1 =  - ui* ( &
          &     + (       exichipm(k) *this%vv_uu(k)  &
          &       +       exichimp(k) *this%vv_dd(k)  &
          &       - conjg(exichipm(k))*conjg(this%vv_uu(k))  &
          &       - conjg(exichimp(k))*conjg(this%vv_dd(k))  ) )
        
        !r_temp(1) =  - ui* ( &
        !  &     + (       exichimp(k) *this%vv_uu(k)  &
        !  &       +       exichipm(k) *this%vv_dd(k)  &
        !  &       - conjg(exichimp(k))*conjg(this%vv_uu(k))  &
        !  &       - conjg(exichipm(k))*conjg(this%vv_dd(k))  ) )
      end if
      !r(k) = this%t(1) * SUM(r_temp)
      r(k) = this%t(1) * r_temp1
    enddo

! path[3] rashba (upper-right direction path around hole)
    !do i = npt + 1, this%np - (this%nh*2)
    do k = this%np_1st + 1, this%np_1st + this%np_rh / 2
      ie = this%ste(p(k)%i)
      fe = this%ste(p(k)%f)
      is = p(k)%i
      fs = p(k)%f
      !print *, "[3]", this%z(is), this%z(fs)
      !r_temp(:) =  0d0
      if(this%z(is) == 1) then ! surface -> there are no small-polarons on the surface
        !   !DIR$ VECTOR ALWAYS
        !do a = 2*this%n_acc+1, 4*this%n_acc
        !enddo
      else  ! bulk
        !   !DIR$ VECTOR ALWAYS
        !do a = 2*this%n_acc+1, 4*this%n_acc
        !  r_temp(a) =  ui*( &
        !  !&       (       exichipp(k) *conjg(this%wf(4*fe-2, a))*this%wf(4*ie-3, a)  &
        !  !&       -       exichimm(k) *conjg(this%wf(4*fe-3, a))*this%wf(4*ie-2, a)  &
        !  !&       - conjg(exichipp(k))*conjg(this%wf(4*ie-3, a))*this%wf(4*fe-2, a)  &
        !  !&       + conjg(exichimm(k))*conjg(this%wf(4*ie-2, a))*this%wf(4*fe-3, a)  &
        !  !&       )*fermi(this%eg(a), this%kbt) & 
        !  &     - (       exichipp(k) *this%wf(4*fe  , a)*conjg(this%wf(4*ie-1, a))  &
        !  &       -       exichimm(k) *this%wf(4*fe-1, a)*conjg(this%wf(4*ie  , a))  &
        !  &       - conjg(exichipp(k))*this%wf(4*ie-1, a)*conjg(this%wf(4*fe  , a))  &
        !  &       + conjg(exichimm(k))*this%wf(4*ie  , a)*conjg(this%wf(4*fe-1, a))  &
        !  !&       )*fermi(-this%eg(a), this%kbt) )
        !  &       ) )
        !end do
        
        r_temp1 =   ui* ( &
          &     - (       exichimm(k) *this%vv_du(k)  &
          &       -       exichipp(k) *this%vv_ud(k)  &
          &       - conjg(exichimm(k))*conjg(this%vv_du(k))  &
          &       + conjg(exichipp(k))*conjg(this%vv_ud(k))  ) )
        
        !r_temp(1) =   ui* ( &
        !  &     - (       exichipp(k) *this%vv_du(k)  &
        !  &       -       exichimm(k) *this%vv_ud(k)  &
        !  &       - conjg(exichipp(k))*conjg(this%vv_du(k))  &
        !  &       + conjg(exichimm(k))*conjg(this%vv_ud(k))  ) )
      end if
      !r(k) = this%lm * SUM(r_temp)
      r(k) = this%lm * r_temp1
    end do

  ! path[4] rashba (upper-left direction path around hole)
    !do i = this%np - (this%nh*2) + 1, this%np
    do k = this%np_1st + this%np_rh/2 + 1, this%np_1st + this%np_rh
      ie = this%ste(p(k)%i)
      fe = this%ste(p(k)%f)
      !r_temp(:) =  0d0
      !print *, "[4]", this%z(is), this%z(fs)
    !   !DIR$ VECTOR ALWAYS
       !do a = 2*this%n_acc+1, 4*this%n_acc
       ! r_temp(a) = - ( &
       !   !&       (       exichipp(k) *conjg(this%wf(4*fe-2, a))*this%wf(4*ie-3, a)  &
       !   !&       +       exichimm(k) *conjg(this%wf(4*fe-3, a))*this%wf(4*ie-2, a)  &
       !   !&       + conjg(exichipp(k))*conjg(this%wf(4*ie-3, a))*this%wf(4*fe-2, a)  &
       !   !&       + conjg(exichimm(k))*conjg(this%wf(4*ie-2, a))*this%wf(4*fe-3, a)  &
       !   !&       )*fermi(this%eg(a), this%kbt) & 
       !   &     - (       exichipp(k) *this%wf(4*fe  , a)*conjg(this%wf(4*ie-1, a))  &
       !   &       +       exichimm(k) *this%wf(4*fe-1, a)*conjg(this%wf(4*ie  , a))  &
       !   &       + conjg(exichipp(k))*this%wf(4*ie-1, a)*conjg(this%wf(4*fe  , a))  &
       !   &       + conjg(exichimm(k))*this%wf(4*ie  , a)*conjg(this%wf(4*fe-1, a))  &
       !   !&       )*fermi(-this%eg(a), this%kbt) )
       !   &       ) )
       !end do
        r_temp1 =   - ( &
          &     - (       exichimm(k) *this%vv_du(k)  &
          &       +       exichipp(k) *this%vv_ud(k)  &
          &       + conjg(exichimm(k))*conjg(this%vv_du(k))  &
          &       + conjg(exichipp(k))*conjg(this%vv_ud(k))  ) )
        !r_temp(1) =   - ( &
        !  &     - (       exichipp(k) *this%vv_du(k)  &
        !  &       +       exichimm(k) *this%vv_ud(k)  &
        !  &       + conjg(exichipp(k))*conjg(this%vv_du(k))  &
        !  &       + conjg(exichimm(k))*conjg(this%vv_ud(k))  ) )
      !r(k) = this%lm * SUM(r_temp)
      r(k) = this%lm * r_temp1
    end do

 ! path[3,4,5,6] p_rh and p_2nd   * p_rh are conteined in p_2nd
    !do i = this%np_1st + 1, this%np_1st + this%np_2nd
    do k = this%np_1st + 1, this%np_1st + this%np_2nd
      ie = this%ste(p(k)%i)
      fe = this%ste(p(k)%f)
      is = p(k)%i
      fs = p(k)%f
      !print *, "[3,4,5,6]", this%z(is), this%z(fs)
      !r_temp(:) =  0d0
      if(this%z(is) == 1) then ! surface
      !do a = 2*this%n_acc+1, 4*this%n_acc
      !  r_temp(a) =  - ui*&
      !  &      (exichipm(i)*conjg(this%wf(2*p(i)%f - 1, a))*this%wf(2*p(i)%i - 1, a)  &
      !  &       + exichimp(i)*conjg(this%wf(2*p(i)%f, a))*this%wf(2*p(i)%i , a)  &
      !  &       - conjg(exichipm(i))*conjg(this%wf(2*p(i)%i - 1, a))*this%wf(2*p(i)%f - 1, a)  &
      !  &       - conjg(exichimp(i))*conjg(this%wf(2*p(i)%i , a))*this%wf(2*p(i)%f, a))
      !enddo
        r_temp1 =  - ui* ( &
          &     + (       exichipm(k) *this%vv_uu(k)*(1d0-this%ne_d(is))*(1d0-this%ne_d(fs))  &
          &       +       exichimp(k) *this%vv_dd(k)*(1d0-this%ne_u(is))*(1d0-this%ne_u(fs))  &
          &       - conjg(exichipm(k))*conjg(this%vv_uu(k))*(1d0-this%ne_d(is))*(1d0-this%ne_d(fs))  &
          &       - conjg(exichimp(k))*conjg(this%vv_dd(k))*(1d0-this%ne_u(is))*(1d0-this%ne_u(fs)) ) )
        !r_temp(1) =  - ui* ( &
        !  &     + (       exichimp(k) *this%vv_uu(k)*(1d0-this%ne_d(is))*(1d0-this%ne_d(fs))  &
        !  &       +       exichipm(k) *this%vv_dd(k)*(1d0-this%ne_u(is))*(1d0-this%ne_u(fs))  &
        !  &       - conjg(exichimp(k))*conjg(this%vv_uu(k))*(1d0-this%ne_d(is))*(1d0-this%ne_d(fs))  &
        !  &       - conjg(exichipm(k))*conjg(this%vv_dd(k))*(1d0-this%ne_u(is))*(1d0-this%ne_u(fs)) ) )
      else ! bulk
        r_temp1 =  - ui* ( &
          &     + (       exichipm(k) *this%vv_uu(k)  &
          &       +       exichimp(k) *this%vv_dd(k)  &
          &       - conjg(exichipm(k))*conjg(this%vv_uu(k))  &
          &       - conjg(exichimp(k))*conjg(this%vv_dd(k))  ) )
       ! r_temp(1) =  - ui* ( &
       !   &     + (       exichimp(k) *this%vv_uu(k)  &
       !   &       +       exichipm(k) *this%vv_dd(k)  &
       !   &       - conjg(exichimp(k))*conjg(this%vv_uu(k))  &
       !   &       - conjg(exichipm(k))*conjg(this%vv_dd(k))  ) )
      end if
      !r(k) = r(k) + this%t(2) * SUM(r_temp)
      r(k) = r(k) + this%t(2) * r_temp1
    enddo

  ! path[7]  inter-layer path
    if(this%np_2nd /= 0) then
      iter_start = this%np_1st + this%np_2nd + 1
    else
      iter_start = this%np_1st + this%np_rh + 1 !! this is ignored p_2nd : when p_2nd = 0
    end if 
    do k = iter_start, this%np
    !do k = this%np_1st + this%np_2nd + 1, this%np
    !do k = this%np_1st + this%np_rh + 1, this%np !! this is ignored p_2nd : when p_2nd = 0
      ie = this%ste(p(k)%i)
      fe = this%ste(p(k)%f)
      is = p(k)%i
      fs = p(k)%f
      !r_temp(:) =  0d0
      !print *, "[7]", this%z(is), this%z(fs)
      if(this%z(is) == 1) then ! path of from surface to bulk
      ! do a = 2*this%n_acc+1, 4*this%n_acc
      !   r_temp(a) =  - ui*(&
      !   !&       (       exichimp(k) *conjg(this%wf(4*fe-3, a))*this%wf(4*ie-3, a)*(1d0-this%ne_d(is))*(1d0-this%ne_d(fs))  &
      !   !&       +       exichipm(k) *conjg(this%wf(4*fe-2, a))*this%wf(4*ie-2, a)*(1d0-this%ne_u(is))*(1d0-this%ne_u(fs))  &
      !   !&       - conjg(exichimp(k))*conjg(this%wf(4*ie-3, a))*this%wf(4*fe-3, a)*(1d0-this%ne_d(is))*(1d0-this%ne_d(fs))  &
      !   !&       - conjg(exichipm(k))*conjg(this%wf(4*ie-2, a))*this%wf(4*fe-2, a)*(1d0-this%ne_u(is))*(1d0-this%ne_u(fs))  &
      !   !&       )*fermi(this%eg(a), this%kbt) & 
      !   &     + (       exichimp(k) *this%wf(4*fe-1, a)*conjg(this%wf(4*ie-1, a))*(1d0-this%ne_d(is))*(1d0-this%ne_d(fs))  &
      !   &       +       exichipm(k) *this%wf(4*fe  , a)*conjg(this%wf(4*ie  , a))*(1d0-this%ne_u(is))*(1d0-this%ne_u(fs))  &
      !   &       - conjg(exichimp(k))*this%wf(4*ie-1, a)*conjg(this%wf(4*fe-1, a))*(1d0-this%ne_d(is))*(1d0-this%ne_d(fs))  &
      !   &       - conjg(exichipm(k))*this%wf(4*ie  , a)*conjg(this%wf(4*fe  , a))*(1d0-this%ne_u(is))*(1d0-this%ne_u(fs))  &
      !   !&       )*fermi(-this%eg(a), this%kbt) )
      !   &       ) )
      ! enddo
        r_temp1 =  - ui* ( &
          &     + (       exichipm(k) *this%vv_uu(k)*(1d0-this%ne_d(is))*(1d0-this%ne_d(fs))  &
          &       +       exichimp(k) *this%vv_dd(k)*(1d0-this%ne_u(is))*(1d0-this%ne_u(fs))  &
          &       - conjg(exichipm(k))*conjg(this%vv_uu(k))*(1d0-this%ne_d(is))*(1d0-this%ne_d(fs))  &
          &       - conjg(exichimp(k))*conjg(this%vv_dd(k))*(1d0-this%ne_u(is))*(1d0-this%ne_u(fs)) ) )
        !r_temp(1) =  - ui* ( &
        !  &     + (       exichimp(k) *this%vv_uu(k)*(1d0-this%ne_d(is))*(1d0-this%ne_d(fs))  &
        !  &       +       exichipm(k) *this%vv_dd(k)*(1d0-this%ne_u(is))*(1d0-this%ne_u(fs))  &
        !  &       - conjg(exichimp(k))*conjg(this%vv_uu(k))*(1d0-this%ne_d(is))*(1d0-this%ne_d(fs))  &
        !  &       - conjg(exichipm(k))*conjg(this%vv_dd(k))*(1d0-this%ne_u(is))*(1d0-this%ne_u(fs)) ) )
      else ! path of from bulk to bulk -> now, not considering
        !do a = 2*this%n_acc+1, 4*this%n_acc
        !    r_temp(a) =  - ui* ( &
        !    &       (       exichimp(k) *conjg(this%wf(4*fe-3, a))*this%wf(4*ie-3, a)  &
        !    &       +       exichipm(k) *conjg(this%wf(4*fe-2, a))*this%wf(4*ie-2, a)  &
        !    &       - conjg(exichimp(k))*conjg(this%wf(4*ie-3, a))*this%wf(4*fe-3, a)  &
        !    &       - conjg(exichipm(k))*conjg(this%wf(4*ie-2, a))*this%wf(4*fe-2, a)  &
        !    &       )*fermi(this%eg(a), this%kbt) & 
        !    &     + (       exichimp(k) *this%wf(4*fe-1, a)*conjg(this%wf(4*ie-1, a))  &
        !    &       +       exichipm(k) *this%wf(4*fe  , a)*conjg(this%wf(4*ie  , a))  &
        !    &       - conjg(exichimp(k))*this%wf(4*ie-1, a)*conjg(this%wf(4*fe-1, a))  &
        !    &       - conjg(exichipm(k))*this%wf(4*ie  , a)*conjg(this%wf(4*fe  , a))  &
        !    &       )*fermi(-this%eg(a), this%kbt) )
        !end do
      end if
      !r(k) = this%t(3) * SUM(r_temp)
      r(k) = this%t(3) * r_temp1
    enddo
  end function fchi_optimizer_Aem_get_Jp

 ! function fchi_optimizer_Aem_get_Jp2(this, dchi, EV) result(r)
 !   class(fchi_optimizer_Aem), intent(in),target         :: this
 !   real(8), dimension(this%np), intent(in)     :: dchi
 !   real(8), dimension(this%np)                :: r
 !   real(8), dimension(this%ne)                :: r_temp
 !   integer                                   :: i, j, k, a, npt
 !   complex(8)                                :: exichip(this%np), exichim(this%np), echi(this%np), val, EV
 !   type path_pointer
 !     integer, pointer :: i, f
 !   end type path_pointer
 !   type(path_pointer)                       :: p(this%np)
 !   !DIR$ VECTOR ALWAYS
 !   do i = 1, this%np
 !     p(i)%i => this%cpl%value(i)%i
 !     p(i)%f => this%cpl%value(i)%f
 !   enddo
 !   if (this%rspath_flag) then
 !     npt = this%np - 4*this%nh
 !   else
 !     npt = this%np
 !   endif
 !   !DIR$ VECTOR ALWAYS
 !   do i = 1, this%np
 !     exichip(i) = exp(0.5d0*ui*(this%xi(p(i)%f) - this%xi(p(i)%i) + dchi(i) + this%jump(i)))
 !     exichim(i) = exp(0.5d0*ui*(-this%xi(p(i)%f) + this%xi(p(i)%i) + dchi(i) + this%jump(i)))
 !     echi(i) = exp(0.5d0*ui*(dchi(i) + this%jump(i)))
 !   enddo
 !   r = 0.d0
 !   ! path[1,2]
 !   do i = 1, npt
 !   !DIR$ VECTOR ALWAYS
 !     do a = 1, this%ne
 !       r_temp(a) =  - ui*&
 !      &      (exichip(i)*conjg(this%wf(2*p(i)%f - 1, a))*this%wf(2*p(i)%i - 1, a)  &
 !      &       + exichim(i)*conjg(this%wf(2*p(i)%f, a))*this%wf(2*p(i)%i , a)  &
 !      &       - conjg(exichip(i))*conjg(this%wf(2*p(i)%i - 1, a))*this%wf(2*p(i)%f - 1, a)  &
 !      &       - conjg(exichim(i))*conjg(this%wf(2*p(i)%i , a))*this%wf(2*p(i)%f, a))
 !     enddo
 !     r(i) = this%t(1) * SUM(r_temp(1:this%ne))
 !   enddo

 !   ! path[3,4]
 !   do i = npt + 1, npt + this%nh*2, 2
 !   !DIR$ VECTOR ALWAYS
 !       !write(*, *) i, 'a'
 !       !write(*, *) 'path', p(i)%i, p(i)%f
 !     do a = 1, this%ne
 !       r_temp(a) =  ui* &
 !             &      (echi(i) &
 !             &       * exp(-0.5d0*ui*(this%xi(p(i)%f)-this%xi(p(i)%f-this%nx)+this%xi(p(i)%i)-this%xi(p(i)%i+1))) &
 !             &       * conjg(this%wf(2*p(i)%f, a))*this%wf(2*p(i)%i - 1, a)  &
 !             &       * conjg(EV) * exp(-ui*this%xi(p(i)%f-this%nx)) &! Dd* Du, h-y -> h+x
 !             &       - echi(i) &
 !             &       * exp(0.5d0*ui*(this%xi(p(i)%f)-this%xi(p(i)%f-this%nx)+this%xi(p(i)%i)-this%xi(p(i)%i+1))) &
 !             &       * conjg(this%wf(2*p(i)%f-1, a))*this%wf(2*p(i)%i, a)  &
 !             &       * EV * exp(ui*this%xi(p(i)%f-this%nx)) &! Du* Dd, h-y -> h+x
 !             &       - conjg(echi(i)) &
 !             &       * exp(0.5d0*ui*(this%xi(p(i)%f)-this%xi(p(i)%f-this%nx)+this%xi(p(i)%i)-this%xi(p(i)%i+1))) &
 !             &       * this%wf(2*p(i)%f, a)*conjg(this%wf(2*p(i)%i - 1, a))  &
 !             &       * EV * exp(ui*this%xi(p(i)%f-this%nx)) &! Dd* Du, h-y -> h+x
 !             &       + conjg(echi(i)) &
 !             &       * exp(-0.5d0*ui*(this%xi(p(i)%f)-this%xi(p(i)%f-this%nx)+this%xi(p(i)%i)-this%xi(p(i)%i+1))) &
 !             &       * this%wf(2*p(i)%f-1, a)*conjg(this%wf(2*p(i)%i, a))  &
 !             &       * conjg(EV) * exp(-ui*this%xi(p(i)%f-this%nx)) &! Du* Dd, h-y -> h+x
 !           & )
 !     end do
 !     val = sum(r_temp(1:this%ne))
 !     r(i) = this%lm * val
 !   end do
 !   do i = npt + 2, npt + 2*this%nh, 2
 !   !DIR$ VECTOR ALWAYS
 !       !write(*, *) i, 'b'
 !       !write(*, *) 'path', p(i)%i, p(i)%f
 !     do a = 1, this%ne
 !       r_temp(a) =  ui* &
 !             &       ( echi(i) &
 !             &       * exp(-0.5d0*ui*(this%xi(p(i)%f)-this%xi(p(i)%f-1)+this%xi(p(i)%i)-this%xi(p(i)%i+this%nx))) &
 !             &       * conjg(this%wf(2*p(i)%f, a))*this%wf(2*p(i)%i - 1, a)  &
 !             &       * conjg(EV) * exp(-ui*this%xi(p(i)%f-1)) &! Dd* Du, h-x -> h+y
 !             &       - echi(i) &
 !             &       * exp(0.5d0*ui*(this%xi(p(i)%f)-this%xi(p(i)%f-1)+this%xi(p(i)%i)-this%xi(p(i)%i+this%nx))) &
 !             &       * conjg(this%wf(2*p(i)%f-1, a))*this%wf(2*p(i)%i, a)  &
 !             &       * EV * exp(ui*this%xi(p(i)%f-1)) &! , h-x,d -> h+y,u
 !             &       - conjg(echi(i)) &
 !             &       * exp(0.5d0*ui*(this%xi(p(i)%f)-this%xi(p(i)%f-1)+this%xi(p(i)%i)-this%xi(p(i)%i+this%nx))) &
 !             &       * this%wf(2*p(i)%f, a)*conjg(this%wf(2*p(i)%i - 1, a))  &
 !             &       * EV * exp(ui*this%xi(p(i)%f-1)) &! Dd* Du, h-x -> h+y
 !             &       + conjg(echi(i)) &
 !             &       * exp(-0.5d0*ui*(this%xi(p(i)%f)-this%xi(p(i)%f-1)+this%xi(p(i)%i)-this%xi(p(i)%i+this%nx))) &
 !             &       * this%wf(2*p(i)%f-1, a)*conjg(this%wf(2*p(i)%i, a))  &
 !             &       * conjg(EV) * exp(-ui*this%xi(p(i)%f-1)) &! , h-x,d -> h+y,u
 !           & )
 !     end do
 !     val = sum(r_temp(1:this%ne))
 !     r(i) = this%lm * val
 !   end do

 !   do i = npt + this%nh*2 + 1, npt + 4*this%nh, 2
 !   !DIR$ VECTOR ALWAYS
 !       !write(*, *) i, 'c'
 !       !write(*, *) 'path', p(i)%i, p(i)%f
 !     do a = 1, this%ne
 !       r_temp(a) = -  &
 !             &     (echi(i) &
 !             &       * exp(-0.5d0*ui*(this%xi(p(i)%f)-this%xi(p(i)%f-this%nx)+this%xi(p(i)%i)-this%xi(p(i)%i-1))) & !!
 !             &       * conjg(this%wf(2*p(i)%f, a))*this%wf(2*p(i)%i - 1, a)  &
 !             &       * conjg(EV) * exp(-ui*this%xi(p(i)%f-this%nx)) &! h-y, u -> h-x,d
 !             &       + echi(i) &
 !             &       * exp(0.5d0*ui*(this%xi(p(i)%f)-this%xi(p(i)%f-this%nx)+this%xi(p(i)%i)-this%xi(p(i)%i-1))) &
 !             &       * conjg(this%wf(2*p(i)%f-1, a))*this%wf(2*p(i)%i, a)  &
 !             &       * EV * exp(ui*this%xi(p(i)%f-this%nx)) &! h-y d -> h-x u
 !             &       + conjg(echi(i)) &
 !             &       * exp(0.5d0*ui*(this%xi(p(i)%f)-this%xi(p(i)%f-this%nx)+this%xi(p(i)%i)-this%xi(p(i)%i-1))) & !!
 !             &       * this%wf(2*p(i)%f, a)*conjg(this%wf(2*p(i)%i - 1, a))  &
 !             &       * EV * exp(ui*this%xi(p(i)%f-this%nx)) &! h-y, u -> h-x,d
 !             &       + conjg(echi(i)) &
 !             &       * exp(-0.5d0*ui*(this%xi(p(i)%f)-this%xi(p(i)%f-this%nx)+this%xi(p(i)%i)-this%xi(p(i)%i-1))) &
 !             &       * this%wf(2*p(i)%f-1, a)*conjg(this%wf(2*p(i)%i, a))  &
 !             &       * conjg(EV) * exp(-ui*this%xi(p(i)%f-this%nx)) &! h-y d -> h-x u
 !             & )
 !     end do
 !     val = sum(r_temp(1:this%ne))
 !     r(i) = this%lm * val
 !   end do
 !   do i = npt + this%nh*2 + 2, npt + 4*this%nh, 2
 !   !DIR$ VECTOR ALWAYS
 !       !write(*, *) i, 'd'
 !       !write(*, *) 'path', p(i)%i, p(i)%f
 !     do a = 1, this%ne
 !       r_temp(a) = -  &
 !             &       ( echi(i) &
 !             &       * exp(-0.5d0*ui*(this%xi(p(i)%f)-this%xi(p(i)%f+1)+this%xi(p(i)%i)-this%xi(p(i)%i+this%nx))) & !!
 !             &       * conjg(this%wf(2*p(i)%f, a))*this%wf(2*p(i)%i - 1, a)  & !!
 !             &       * conjg(EV) * exp(-ui*this%xi(p(i)%f+1)) &! h+x u -> h+y d !!
 !             &       + echi(i) &
 !             &       * exp(0.5d0*ui*(this%xi(p(i)%f)-this%xi(p(i)%f+1)+this%xi(p(i)%i)-this%xi(p(i)%i+this%nx))) & !!
 !             &       * conjg(this%wf(2*p(i)%f-1, a))*this%wf(2*p(i)%i, a)  & !!
 !             &       * EV * exp(ui*this%xi(p(i)%f+1)) &! , h+x d -> h+y u !!
 !             &       + conjg(echi(i)) &
 !             &       * exp(0.5d0*ui*(this%xi(p(i)%f)-this%xi(p(i)%f+1)+this%xi(p(i)%i)-this%xi(p(i)%i+this%nx))) & !!
 !             &       * this%wf(2*p(i)%f, a)*conjg(this%wf(2*p(i)%i - 1, a))  & !!
 !             &       * EV * exp(ui*this%xi(p(i)%f+1)) &! h+x u -> h+y d !!
 !             &       + conjg(echi(i)) &
 !             &       * exp(-0.5d0*ui*(this%xi(p(i)%f)-this%xi(p(i)%f+1)+this%xi(p(i)%i)-this%xi(p(i)%i+this%nx))) & !!
 !             &       * this%wf(2*p(i)%f-1, a)*conjg(this%wf(2*p(i)%i, a))  & !!
 !             &       * conjg(EV) * exp(-ui*this%xi(p(i)%f+1)) &! , h+x d -> h+y u !!
 !             & )
 !     end do
 !     val = sum(r_temp(1:this%ne))
 !     r(i) = this%lm * val
 !   end do

 !   !do i=1,this%np
 !   !p = this%cpl%index(i)
 !   !p(i)%f = p%f
 !   !p(i)%i = p%i
 !   !write(*,*)i,p(i)%i ,p(i)%f,r(i)
 !   !enddo
 ! end function fchi_optimizer_Aem_get_Jp2

  function fchi_optimizer_Aem_get_dJdf(this, dchi) result(r)
    class(fchi_optimizer_Aem), intent(in), target  :: this
    real(8), dimension(this%np), intent(in)     :: dchi
    real(8), dimension(this%np)                :: r
    !real(8), dimension(4*this%n_acc)               :: r_temp
    integer                                   :: i, j, k, a, npt, ie, is,fe ,fs, iter_start
    complex(8)                                :: exichipm(this%np), exichimp(this%np), exichipp(this%np), exichimm(this%np)
    type path_pointer
      integer, pointer :: i, f
    end type path_pointer
    type(path_pointer)                       :: p(this%np)
    real(8)               :: r_temp1
    !   !DIR$ VECTOR ALWAYS
    do k = 1, this%np
      p(k)%i => this%cpl%value(k)%i
      p(k)%f => this%cpl%value(k)%f
    enddo
    !   !DIR$ VECTOR ALWAYS
    do k = 1, this%np
      exichipm(k) = exp(0.5d0*ui*( this%xi(p(k)%f) - this%xi(p(k)%i) + dchi(k) + this%jump(k)))
      exichimp(k) = exp(0.5d0*ui*(-this%xi(p(k)%f) + this%xi(p(k)%i) + dchi(k) + this%jump(k)))
      exichipp(k) = exp(0.5d0*ui*( this%xi(p(k)%f) + this%xi(p(k)%i) + dchi(k) + this%jump(k)))
      exichimm(k) = exp(0.5d0*ui*(-this%xi(p(k)%f) - this%xi(p(k)%i) + dchi(k) + this%jump(k)))
    enddo
    !if (this%rspath_flag) then
    !  npt = this%np - 4*this%nh
    !else
    !  npt = this%np
    !endif
    r = 0.d0
! path[1,2] the nearest path
    !do i = 1, npt
    do k = 1, this%np_1st
      ie = this%ste(p(k)%i)
      fe = this%ste(p(k)%f)
      is = p(k)%i
      fs = p(k)%f
      !r_temp(:) =  0d0
    !   !DIR$ VECTOR ALWAYS
      if(this%z(is) == 1) then ! surface
        !do a = 2*this%n_acc+1, 4*this%n_acc
        !  r_temp(a) =  0.5d0*this%t(1)*( &
        !  !&       (       exichimp(k) *conjg(this%wf(4*fe-3, a))*this%wf(4*ie-3, a)*(1d0-this%ne_d(is))*(1d0-this%ne_d(fs))  &
        !  !&       +       exichipm(k) *conjg(this%wf(4*fe-2, a))*this%wf(4*ie-2, a)*(1d0-this%ne_u(is))*(1d0-this%ne_u(fs))  &
        !  !&       + conjg(exichimp(k))*conjg(this%wf(4*ie-3, a))*this%wf(4*fe-3, a)*(1d0-this%ne_d(is))*(1d0-this%ne_d(fs))  &
        !  !&       + conjg(exichipm(k))*conjg(this%wf(4*ie-2, a))*this%wf(4*fe-2, a)*(1d0-this%ne_u(is))*(1d0-this%ne_u(fs))  &
        !  !&       )*fermi(this%eg(a), this%kbt) & 
        !  &     + (       exichimp(k) *this%wf(4*fe-1, a)*conjg(this%wf(4*ie-1, a))*(1d0-this%ne_d(is))*(1d0-this%ne_d(fs))  &
        !  &       +       exichipm(k) *this%wf(4*fe  , a)*conjg(this%wf(4*ie  , a))*(1d0-this%ne_u(is))*(1d0-this%ne_u(fs))  &
        !  &       + conjg(exichimp(k))*this%wf(4*ie-1, a)*conjg(this%wf(4*fe-1, a))*(1d0-this%ne_d(is))*(1d0-this%ne_d(fs))  &
        !  &       + conjg(exichipm(k))*this%wf(4*ie  , a)*conjg(this%wf(4*fe  , a))*(1d0-this%ne_u(is))*(1d0-this%ne_u(fs))  &
        !  !&       )*fermi(-this%eg(a), this%kbt) )
        !  &       ) )
        !enddo
        r_temp1 =  0.5d0*this%t(1)* ( &
          &     + (       exichipm(k) *this%vv_uu(k)*(1d0-this%ne_d(is))*(1d0-this%ne_d(fs))  &
          &       +       exichimp(k) *this%vv_dd(k)*(1d0-this%ne_u(is))*(1d0-this%ne_u(fs))  &
          &       + conjg(exichipm(k))*conjg(this%vv_uu(k))*(1d0-this%ne_d(is))*(1d0-this%ne_d(fs))  &
          &       + conjg(exichimp(k))*conjg(this%vv_dd(k))*(1d0-this%ne_u(is))*(1d0-this%ne_u(fs)) ) )
       ! r_temp(1) =  0.5d0*this%t(1)* ( &
       !   &     + (       exichimp(k) *this%vv_uu(k)*(1d0-this%ne_d(is))*(1d0-this%ne_d(fs))  &
       !   &       +       exichipm(k) *this%vv_dd(k)*(1d0-this%ne_u(is))*(1d0-this%ne_u(fs))  &
       !   &       + conjg(exichimp(k))*conjg(this%vv_uu(k))*(1d0-this%ne_d(is))*(1d0-this%ne_d(fs))  &
       !   &       + conjg(exichipm(k))*conjg(this%vv_dd(k))*(1d0-this%ne_u(is))*(1d0-this%ne_u(fs)) ) )
      else
        !do a = 2*this%n_acc+1, 4*this%n_acc
        !  r_temp(a) =  0.5d0*this%t(1)*( &
        !  !&       (       exichimp(k) *conjg(this%wf(4*fe-3, a))*this%wf(4*ie-3, a)  &
        !  !&       +       exichipm(k) *conjg(this%wf(4*fe-2, a))*this%wf(4*ie-2, a)  &
        !  !&       + conjg(exichimp(k))*conjg(this%wf(4*ie-3, a))*this%wf(4*fe-3, a)  &
        !  !&       + conjg(exichipm(k))*conjg(this%wf(4*ie-2, a))*this%wf(4*fe-2, a)  &
        !  !&       )*fermi(this%eg(a), this%kbt) & 
        !  &     + (       exichimp(k) *this%wf(4*fe-1, a)*conjg(this%wf(4*ie-1, a))  &
        !  &       +       exichipm(k) *this%wf(4*fe  , a)*conjg(this%wf(4*ie  , a))  &
        !  &       + conjg(exichimp(k))*this%wf(4*ie-1, a)*conjg(this%wf(4*fe-1, a))  &
        !  &       + conjg(exichipm(k))*this%wf(4*ie  , a)*conjg(this%wf(4*fe  , a))  &
        !  !&       )*fermi(-this%eg(a), this%kbt) )
        !  &       ) )
        !enddo
        r_temp1 =  0.5d0*this%t(1)* ( &
          &     + (       exichipm(k) *this%vv_uu(k)  &
          &       +       exichimp(k) *this%vv_dd(k)  &
          &       + conjg(exichipm(k))*conjg(this%vv_uu(k))  &
          &       + conjg(exichimp(k))*conjg(this%vv_dd(k)) ) )
        !r_temp(1) =  0.5d0*this%t(1)* ( &
        !  &     + (       exichimp(k) *this%vv_uu(k)  &
        !  &       +       exichipm(k) *this%vv_dd(k)  &
        !  &       + conjg(exichimp(k))*conjg(this%vv_uu(k))  &
        !  &       + conjg(exichipm(k))*conjg(this%vv_dd(k)) ) )
      end if
      !r(k) = SUM(r_temp)
      r(k) = r_temp1
    enddo

    ! path[3] rashba (upper-right direction path around hole)
    !do i = npt + 1, this%np - (this%nh*2)
    do k = this%np_1st + 1, this%np_1st + this%np_rh / 2
      ie = this%ste(p(k)%i)
      fe = this%ste(p(k)%f)
      is = p(k)%i
      fs = p(k)%f
    !   !DIR$ VECTOR ALWAYS
      !do a = 2*this%n_acc+1, 4*this%n_acc
      !  r_temp(a) =  - 0.5d0*this%lm*( &
      !    !&       (       exichipp(k) *conjg(this%wf(4*fe-2, a))*this%wf(4*ie-3, a)  &
      !    !&       -       exichimm(k) *conjg(this%wf(4*fe-3, a))*this%wf(4*ie-2, a)  &
      !    !&       + conjg(exichipp(k))*conjg(this%wf(4*ie-3, a))*this%wf(4*fe-2, a)  &
      !    !&       - conjg(exichimm(k))*conjg(this%wf(4*ie-2, a))*this%wf(4*fe-3, a)  &
      !    !&       )*fermi(this%eg(a), this%kbt) & 
      !    &     - (       exichipp(k) *this%wf(4*fe  , a)*conjg(this%wf(4*ie-1, a))  &
      !    &       -       exichimm(k) *this%wf(4*fe-1, a)*conjg(this%wf(4*ie  , a))  &
      !    &       + conjg(exichipp(k))*this%wf(4*ie-1, a)*conjg(this%wf(4*fe  , a))  &
      !    &       - conjg(exichimm(k))*this%wf(4*ie  , a)*conjg(this%wf(4*fe-1, a))  &
      !    !&       )*fermi(-this%eg(a), this%kbt) )
      !    &       ) )
      !  end do
        r_temp1 =  - 0.5d0*this%lm*( &
          &     - (       exichimm(k) *this%vv_du(k)  &
          &       -       exichipp(k) *this%vv_ud(k)  &
          &       + conjg(exichimm(k))*conjg(this%vv_du(k))  &
          &       - conjg(exichipp(k))*conjg(this%vv_ud(k))  ) )
       ! r_temp(1) =  - 0.5d0*this%lm*( &
       !   &     - (       exichipp(k) *this%vv_du(k)  &
       !   &       -       exichimm(k) *this%vv_ud(k)  &
       !   &       + conjg(exichipp(k))*conjg(this%vv_du(k))  &
       !   &       - conjg(exichimm(k))*conjg(this%vv_ud(k))  ) )
      !r(k) = SUM(r_temp)
      r(k) = r_temp1
    end do

    ! path[4] rashba (upper-left direction path around hole)
    !do i = this%np - (this%nh*2) + 1, this%np
    do k = this%np_1st + this%np_rh/2 + 1, this%np_1st + this%np_rh
      ie = this%ste(p(k)%i)
      fe = this%ste(p(k)%f)
      is = p(k)%i
      fs = p(k)%f
    !   !DIR$ VECTOR ALWAYS
      !do a = 2*this%n_acc+1, 4*this%n_acc
      !  r_temp(a) =   ui*0.5d0*this%lm*( &
      !    !&       (       exichipp(k) *conjg(this%wf(4*fe-2, a))*this%wf(4*ie-3, a)  &
      !    !&       +       exichimm(k) *conjg(this%wf(4*fe-3, a))*this%wf(4*ie-2, a)  &
      !    !&       - conjg(exichipp(k))*conjg(this%wf(4*ie-3, a))*this%wf(4*fe-2, a)  &
      !    !&       - conjg(exichimm(k))*conjg(this%wf(4*ie-2, a))*this%wf(4*fe-3, a)  &
      !    !&       )*fermi(this%eg(a), this%kbt) & 
      !    &     - (       exichipp(k) *this%wf(4*fe  , a)*conjg(this%wf(4*ie-1, a))  &
      !    &       +       exichimm(k) *this%wf(4*fe-1, a)*conjg(this%wf(4*ie  , a))  &
      !    &       - conjg(exichipp(k))*this%wf(4*ie-1, a)*conjg(this%wf(4*fe  , a))  &
      !    &       - conjg(exichimm(k))*this%wf(4*ie  , a)*conjg(this%wf(4*fe-1, a))  &
      !    !&       )*fermi(-this%eg(a), this%kbt) )
      !    &       ) )
      !  end do
        r_temp1 =  ui*0.5d0*this%lm*( &
          &     - (       exichimm(k) *this%vv_du(k)  &
          &       +       exichipp(k) *this%vv_ud(k)  &
          &       - conjg(exichimm(k))*conjg(this%vv_du(k))  &
          &       - conjg(exichipp(k))*conjg(this%vv_ud(k))  ) )
        !r_temp(1) =  ui*0.5d0*this%lm*( &
        !  &     - (       exichipp(k) *this%vv_du(k)  &
        !  &       +       exichimm(k) *this%vv_ud(k)  &
        !  &       - conjg(exichipp(k))*conjg(this%vv_du(k))  &
        !  &       - conjg(exichimm(k))*conjg(this%vv_ud(k))  ) )
      !r(k) = SUM(r_temp)
      r(k) = r_temp1
    end do

 ! path[3,4,5,6] p_rh and p_2nd   * p_rh are conteined in p_2nd
!    !do i = this%np_1st + 1, this%np_1st + this%np_2nd
!    !  do a = 1, this%ne
!    !    r_temp(a) =  0.5d0*this%t(2)*&
!    !   &      (exichipm(i)*conjg(this%wf(2*p(i)%f - 1, a))*this%wf(2*p(i)%i - 1, a)  &
!    !   &       + exichimp(i)*conjg(this%wf(2*p(i)%f, a))*this%wf(2*p(i)%i , a)  &
!    !   &       + conjg(exichipm(i))*conjg(this%wf(2*p(i)%i - 1, a))*this%wf(2*p(i)%f - 1, a)  &
!    !   &       + conjg(exichimp(i))*conjg(this%wf(2*p(i)%i , a))*this%wf(2*p(i)%f, a))
!    !  enddo
!    !  r(i) = r(i) + SUM(r_temp(1:this%ne))
!    !enddo
    do k = this%np_1st + 1, this%np_1st + this%np_2nd
      ie = this%ste(p(k)%i)
      fe = this%ste(p(k)%f)
      is = p(k)%i
      fs = p(k)%f
      !r_temp(:) =  0d0
      if(this%z(is) == 1) then ! surface
        r_temp1 =  0.5d0*this%t(2)* ( &
          &     + (       exichipm(k) *this%vv_uu(k)*(1d0-this%ne_d(is))*(1d0-this%ne_d(fs))  &
          &       +       exichimp(k) *this%vv_dd(k)*(1d0-this%ne_u(is))*(1d0-this%ne_u(fs))  &
          &       + conjg(exichipm(k))*conjg(this%vv_uu(k))*(1d0-this%ne_d(is))*(1d0-this%ne_d(fs))  &
          &       + conjg(exichimp(k))*conjg(this%vv_dd(k))*(1d0-this%ne_u(is))*(1d0-this%ne_u(fs)) ) )
        !r_temp(1) =  0.5d0*this%t(2)* ( &
        !  &     + (       exichimp(k) *this%vv_uu(k)*(1d0-this%ne_d(is))*(1d0-this%ne_d(fs))  &
        !  &       +       exichipm(k) *this%vv_dd(k)*(1d0-this%ne_u(is))*(1d0-this%ne_u(fs))  &
        !  &       + conjg(exichimp(k))*conjg(this%vv_uu(k))*(1d0-this%ne_d(is))*(1d0-this%ne_d(fs))  &
        !  &       + conjg(exichipm(k))*conjg(this%vv_dd(k))*(1d0-this%ne_u(is))*(1d0-this%ne_u(fs)) ) )
      else ! bulk
        r_temp1 =  0.5d0*this%t(2)* ( &
          &     + (       exichipm(k) *this%vv_uu(k)  &
          &       +       exichimp(k) *this%vv_dd(k)  &
          &       + conjg(exichipm(k))*conjg(this%vv_uu(k))  &
          &       + conjg(exichimp(k))*conjg(this%vv_dd(k))  ) )
        !r_temp(1) =  0.5d0*this%t(2)* ( &
        !  &     + (       exichimp(k) *this%vv_uu(k)  &
        !  &       +       exichipm(k) *this%vv_dd(k)  &
        !  &       + conjg(exichimp(k))*conjg(this%vv_uu(k))  &
        !  &       + conjg(exichipm(k))*conjg(this%vv_dd(k))  ) )
      end if
      !r(k) = r(k) + SUM(r_temp)
      r(k) = r(k) + r_temp1
    enddo

! path[7]  inter-layer path
    if(this%np_2nd /= 0) then
      iter_start = this%np_1st + this%np_2nd + 1
    else
      iter_start = this%np_1st + this%np_rh + 1 !! this is ignored p_2nd : when p_2nd = 0
    end if 
    do k = iter_start, this%np
    !do k = this%np_1st + this%np_2nd + 1, this%np
    !do k = this%np_1st + this%np_rh + 1, this%np !! this is ignored p_2nd : when p_2nd = 0
      ie = this%ste(p(k)%i)
      fe = this%ste(p(k)%f)
      is = p(k)%i
      fs = p(k)%f
      !r_temp(:) =  0d0
      if(this%z(is) == 1) then ! path of from surface to bulk
        !   !DIR$ VECTOR ALWAYS
        !do a = 2*this%n_acc+1, 4*this%n_acc
        !  r_temp(a) =  0.5d0*this%t(3)*(&
        !  !&       (       exichimp(k) *conjg(this%wf(4*fe-3, a))*this%wf(4*ie-3, a)*(1d0-this%ne_d(is))*(1d0-this%ne_d(fs))  &
        !  !&       +       exichipm(k) *conjg(this%wf(4*fe-2, a))*this%wf(4*ie-2, a)*(1d0-this%ne_u(is))*(1d0-this%ne_u(fs))  &
        !  !&       + conjg(exichimp(k))*conjg(this%wf(4*ie-3, a))*this%wf(4*fe-3, a)*(1d0-this%ne_d(is))*(1d0-this%ne_d(fs))  &
        !  !&       + conjg(exichipm(k))*conjg(this%wf(4*ie-2, a))*this%wf(4*fe-2, a)*(1d0-this%ne_u(is))*(1d0-this%ne_u(fs))  &
        !  !&       )*fermi(this%eg(a), this%kbt) & 
        !  &     + (       exichimp(k) *this%wf(4*fe-1, a)*conjg(this%wf(4*ie-1, a))*(1d0-this%ne_d(is))*(1d0-this%ne_d(fs))  &
        !  &       +       exichipm(k) *this%wf(4*fe  , a)*conjg(this%wf(4*ie  , a))*(1d0-this%ne_u(is))*(1d0-this%ne_u(fs))  &
        !  &       + conjg(exichimp(k))*this%wf(4*ie-1, a)*conjg(this%wf(4*fe-1, a))*(1d0-this%ne_d(is))*(1d0-this%ne_d(fs))  &
        !  &       + conjg(exichipm(k))*this%wf(4*ie  , a)*conjg(this%wf(4*fe  , a))*(1d0-this%ne_u(is))*(1d0-this%ne_u(fs))  &
        !  !&       )*fermi(-this%eg(a), this%kbt) )
        !  &       ) )
        !enddo
        r_temp1 =  0.5d0*this%t(3)* ( &
          &     + (       exichipm(k) *this%vv_uu(k)*(1d0-this%ne_d(is))*(1d0-this%ne_d(fs))  &
          &       +       exichimp(k) *this%vv_dd(k)*(1d0-this%ne_u(is))*(1d0-this%ne_u(fs))  &
          &       + conjg(exichipm(k))*conjg(this%vv_uu(k))*(1d0-this%ne_d(is))*(1d0-this%ne_d(fs))  &
          &       + conjg(exichimp(k))*conjg(this%vv_dd(k))*(1d0-this%ne_u(is))*(1d0-this%ne_u(fs)) ) )
        !r_temp(1) =  0.5d0*this%t(3)* ( &
        !  &     + (       exichimp(k) *this%vv_uu(k)*(1d0-this%ne_d(is))*(1d0-this%ne_d(fs))  &
        !  &       +       exichipm(k) *this%vv_dd(k)*(1d0-this%ne_u(is))*(1d0-this%ne_u(fs))  &
        !  &       + conjg(exichimp(k))*conjg(this%vv_uu(k))*(1d0-this%ne_d(is))*(1d0-this%ne_d(fs))  &
        !  &       + conjg(exichipm(k))*conjg(this%vv_dd(k))*(1d0-this%ne_u(is))*(1d0-this%ne_u(fs)) ) )
      else ! path of from bulk to bulk -> now, not considering
        !   !DIR$ VECTOR ALWAYS
        !do a = 2*this%n_acc+1, 4*this%n_acc
        !  r_temp(a) =  0.5d0*this%t(3)*(&
        !  &       (       exichimp(k) *conjg(this%wf(4*fe-3, a))*this%wf(4*ie-3, a)  &
        !  &       +       exichipm(k) *conjg(this%wf(4*fe-2, a))*this%wf(4*ie-2, a)  &
        !  &       + conjg(exichimp(k))*conjg(this%wf(4*ie-3, a))*this%wf(4*fe-3, a)  &
        !  &       + conjg(exichipm(k))*conjg(this%wf(4*ie-2, a))*this%wf(4*fe-2, a)  &
        !  &       )*fermi(this%eg(a), this%kbt) & 
        !  &     + (       exichimp(k) *this%wf(4*fe-1, a)*conjg(this%wf(4*ie-1, a))  &
        !  &       +       exichipm(k) *this%wf(4*fe  , a)*conjg(this%wf(4*ie  , a))  &
        !  &       + conjg(exichimp(k))*this%wf(4*ie-1, a)*conjg(this%wf(4*fe-1, a))  &
        !  &       + conjg(exichipm(k))*this%wf(4*ie  , a)*conjg(this%wf(4*fe  , a))  &
        !  &       )*fermi(-this%eg(a), this%kbt) )
        !enddo
      end if
      !r(k) = SUM(r_temp)
      r(k) = r_temp1
    enddo
  end function fchi_optimizer_Aem_get_dJdf

  !function fchi_optimizer_Aem_get_dJdf2(this, dchi, EV) result(r)
  !  class(fchi_optimizer_Aem), intent(in), target  :: this
  !  real(8), dimension(this%np), intent(in)     :: dchi
  !  real(8), dimension(this%np)                :: r
  !  real(8), dimension(this%ne)                :: r_temp
  !  integer                                   :: i, j, k, a, npt
  !  complex(8)                                :: echi(this%np),exichip(this%np), exichim(this%np), EV, val
  !  type path_pointer
  !    integer, pointer :: i, f
  !  end type path_pointer
  !  type(path_pointer)                       :: p(this%np)
  !  !DIR$ VECTOR ALWAYS
  !  do i = 1, this%np
  !    p(i)%i => this%cpl%value(i)%i
  !    p(i)%f => this%cpl%value(i)%f
  !  enddo
  !  !DIR$ VECTOR ALWAYS
  !  do i = 1, this%np
  !    echi(i) = exp(0.5d0*ui*(dchi(i) + this%jump(i)))
  !  enddo
  !  if (this%rspath_flag) then
  !    npt = this%np - 4*this%nh
  !  else
  !    npt = this%np
  !  endif
  !  r = 0.d0
  !  ! path[1,2]
  !  do i = 1, npt
  !  !DIR$ VECTOR ALWAYS
  !    do a = 1, this%ne
  !      r_temp(a) =  0.5d0*this%t(1)*&
  !      &      (exichip(i)*conjg(this%wf(2*p(i)%f - 1, a))*this%wf(2*p(i)%i- 1, a)  &
  !      &       + exichim(i)*conjg(this%wf(2*p(i)%f, a))*this%wf(2*p(i)%i, a)  &
  !      &       + conjg(exichip(i))*conjg(this%wf(2*p(i)%i- 1, a))*this%wf(2*p(i)%f - 1, a)  &
  !      &       + conjg(exichim(i))*conjg(this%wf(2*p(i)%i, a))*this%wf(2*p(i)%f, a))
  !    enddo
  !    r(i) = SUM(r_temp(1:this%ne))
  !  enddo

! !path[3,4]
  !  do i = npt + 1, npt + this%nh*2, 2
  !  !DIR$ VECTOR ALWAYS
  !    do a = 1, this%ne
  !      r_temp(a) =  - 0.5d0*this%lm* &
  !            &      (echi(i) &
  !            &       * exp(-0.5d0*ui*(this%xi(p(i)%f)-this%xi(p(i)%f-this%nx)+this%xi(p(i)%i)-this%xi(p(i)%i+1))) &
  !            &       * conjg(this%wf(2*p(i)%f, a))*this%wf(2*p(i)%i - 1, a)  &
  !            &       * conjg(EV) * exp(-ui*this%xi(p(i)%f-this%nx)) &! Dd* Du, h-y -> h+x
  !            &       - echi(i) &
  !            &       * exp(0.5d0*ui*(this%xi(p(i)%f)-this%xi(p(i)%f-this%nx)+this%xi(p(i)%i)-this%xi(p(i)%i+1))) &
  !            &       * conjg(this%wf(2*p(i)%f-1, a))*this%wf(2*p(i)%i, a)  &
  !            &       * EV * exp(ui*this%xi(p(i)%f-this%nx)) &! Du* Dd, h-y -> h+x
  !            &       + conjg(echi(i)) &
  !            &       * exp(0.5d0*ui*(this%xi(p(i)%f)-this%xi(p(i)%f-this%nx)+this%xi(p(i)%i)-this%xi(p(i)%i+1))) &
  !            &       * this%wf(2*p(i)%f, a)*conjg(this%wf(2*p(i)%i - 1, a))  &
  !            &       * EV * exp(ui*this%xi(p(i)%f-this%nx)) &! Dd* Du, h-y -> h+x
  !            &       - conjg(echi(i)) &
  !            &       * exp(-0.5d0*ui*(this%xi(p(i)%f)-this%xi(p(i)%f-this%nx)+this%xi(p(i)%i)-this%xi(p(i)%i+1))) &
  !            &       * this%wf(2*p(i)%f-1, a)*conjg(this%wf(2*p(i)%i, a))  &
  !            &       * conjg(EV) * exp(-ui*this%xi(p(i)%f-this%nx)) &! Du* Dd, h-y -> h+x
  !          & )
  !    end do
  !    val= SUM(r_temp(1:this%ne))
  !    r(i) = val
  !  end do
  !  do i = npt + 2, npt + 2*this%nh, 2
  !  !DIR$ VECTOR ALWAYS
  !    do a = 1, this%ne
  !      r_temp(a) =  - 0.5d0*this%lm* &
  !            &       ( echi(i) &
  !            &       * exp(-0.5d0*ui*(this%xi(p(i)%f)-this%xi(p(i)%f-1)+this%xi(p(i)%i)-this%xi(p(i)%i+this%nx))) &
  !            &       * conjg(this%wf(2*p(i)%f, a))*this%wf(2*p(i)%i - 1, a)  &
  !            &       * conjg(EV) * exp(-ui*this%xi(p(i)%f-1)) &! Dd* Du, h-x -> h+y
  !            &       - echi(i) &
  !            &       * exp(0.5d0*ui*(this%xi(p(i)%f)-this%xi(p(i)%f-1)+this%xi(p(i)%i)-this%xi(p(i)%i+this%nx))) &
  !            &       * conjg(this%wf(2*p(i)%f-1, a))*this%wf(2*p(i)%i, a)  &
  !            &       * EV * exp(ui*this%xi(p(i)%f-1)) &! , h-x,d -> h+y,u
  !            &       + conjg(echi(i)) &
  !            &       * exp(0.5d0*ui*(this%xi(p(i)%f)-this%xi(p(i)%f-1)+this%xi(p(i)%i)-this%xi(p(i)%i+this%nx))) &
  !            &       * this%wf(2*p(i)%f, a)*conjg(this%wf(2*p(i)%i - 1, a))  &
  !            &       * EV * exp(ui*this%xi(p(i)%f-1)) &! Dd* Du, h-x -> h+y
  !            &       - conjg(echi(i)) &
  !            &       * exp(-0.5d0*ui*(this%xi(p(i)%f)-this%xi(p(i)%f-1)+this%xi(p(i)%i)-this%xi(p(i)%i+this%nx))) &
  !            &       * this%wf(2*p(i)%f-1, a)*conjg(this%wf(2*p(i)%i, a))  &
  !            &       * conjg(EV) * exp(-ui*this%xi(p(i)%f-1)) &! , h-x,d -> h+y,u
  !          & )
  !    end do
  !    val= SUM(r_temp(1:this%ne))
  !    r(i) = val
  !  end do

  !  do i = npt + this%nh*2 + 1, npt + this%nh*4, 2
  !  !DIR$ VECTOR ALWAYS
  !    do a = 1, this%ne
  !      r_temp(a) =   ui*0.5d0*this%lm* &
  !            &     (echi(i) &
  !            &       * exp(-0.5d0*ui*(this%xi(p(i)%f)-this%xi(p(i)%f-this%nx)+this%xi(p(i)%i)-this%xi(p(i)%i-1))) & !!
  !            &       * conjg(this%wf(2*p(i)%f, a))*this%wf(2*p(i)%i - 1, a)  &
  !            &       * conjg(EV) * exp(-ui*this%xi(p(i)%f-this%nx)) &! h-y, u -> h-x,d
  !            &       + echi(i) &
  !            &       * exp(0.5d0*ui*(this%xi(p(i)%f)-this%xi(p(i)%f-this%nx)+this%xi(p(i)%i)-this%xi(p(i)%i-1))) &
  !            &       * conjg(this%wf(2*p(i)%f-1, a))*this%wf(2*p(i)%i, a)  &
  !            &       * EV * exp(ui*this%xi(p(i)%f-this%nx)) &! h-y d -> h-x u
  !            &       - conjg(echi(i)) &
  !            &       * exp(0.5d0*ui*(this%xi(p(i)%f)-this%xi(p(i)%f-this%nx)+this%xi(p(i)%i)-this%xi(p(i)%i-1))) & !!
  !            &       * this%wf(2*p(i)%f, a)*conjg(this%wf(2*p(i)%i - 1, a))  &
  !            &       * EV * exp(ui*this%xi(p(i)%f-this%nx)) &! h-y, u -> h-x,d
  !            &       - conjg(echi(i)) &
  !            &       * exp(-0.5d0*ui*(this%xi(p(i)%f)-this%xi(p(i)%f-this%nx)+this%xi(p(i)%i)-this%xi(p(i)%i-1))) &
  !            &       * this%wf(2*p(i)%f-1, a)*conjg(this%wf(2*p(i)%i, a))  &
  !            &       * conjg(EV) * exp(-ui*this%xi(p(i)%f-this%nx)) &! h-y d -> h-x u
  !            & )
  !    end do
  !    val= SUM(r_temp(1:this%ne))
  !    r(i) = val
  !  end do
  !  do i = npt + this%nh*2 + 2, npt + this%nh*4, 2
  !  !DIR$ VECTOR ALWAYS
  !    do a = 1, this%ne
  !      r_temp(a) =   ui*0.5d0*this%lm* &
  !            &       ( echi(i) &
  !            &       * exp(-0.5d0*ui*(this%xi(p(i)%f)-this%xi(p(i)%f+1)+this%xi(p(i)%i)-this%xi(p(i)%i+this%nx))) & !!
  !            &       * conjg(this%wf(2*p(i)%f, a))*this%wf(2*p(i)%i - 1, a)  & !!
  !            &       * conjg(EV) * exp(-ui*this%xi(p(i)%f+1)) &! h+x u -> h+y d !!
  !            &       + echi(i) &
  !            &       * exp(0.5d0*ui*(this%xi(p(i)%f)-this%xi(p(i)%f+1)+this%xi(p(i)%i)-this%xi(p(i)%i+this%nx))) & !!
  !            &       * conjg(this%wf(2*p(i)%f-1, a))*this%wf(2*p(i)%i, a)  & !!
  !            &       * EV * exp(ui*this%xi(p(i)%f+1)) &! , h+x d -> h+y u !!
  !            &       - conjg(echi(i)) &
  !            &       * exp(0.5d0*ui*(this%xi(p(i)%f)-this%xi(p(i)%f+1)+this%xi(p(i)%i)-this%xi(p(i)%i+this%nx))) & !!
  !            &       * this%wf(2*p(i)%f, a)*conjg(this%wf(2*p(i)%i - 1, a))  & !!
  !            &       * EV * exp(ui*this%xi(p(i)%f+1)) &! h+x u -> h+y d !!
  !            &       - conjg(echi(i)) &
  !            &       * exp(-0.5d0*ui*(this%xi(p(i)%f)-this%xi(p(i)%f+1)+this%xi(p(i)%i)-this%xi(p(i)%i+this%nx))) & !!
  !            &       * this%wf(2*p(i)%f-1, a)*conjg(this%wf(2*p(i)%i, a))  & !!
  !            &       * conjg(EV) * exp(-ui*this%xi(p(i)%f+1)) &! , h+x d -> h+y u !!
  !            & )
  !    end do
  !    val= SUM(r_temp(1:this%ne))
  !    r(i) = val
  !  end do
  !end function fchi_optimizer_Aem_get_dJdf2

  subroutine fchi_optimizer_Aem_set_lcp(this)
    class(fchi_optimizer_Aem), intent(inout) :: this
    integer                              :: i, j, k, l
    if (allocated(this%lcp)) deallocate (this%lcp)
    allocate (this%lcp(this%nc, this%np))
    this%lcp = 0
    do i = 1, this%nc
      do j = 1, this%cil(i)%length()
        k = this%cil(i)%index(j)
        l = 1
        if (k < 0) l = -1
        this%lcp(i, abs(k)) = l
      end do
    end do
    end subroutine fchi_optimizer_Aem_set_lcp

  function fchi_optimizer_Aem_get_n(this, dJdf) result(r)
    use ieee_arithmetic
    class(fchi_optimizer_Aem), intent(in)       :: this
    real(8), dimension(this%np), intent(in)  :: dJdf
    real(8), dimension(this%np, this%np)     :: r
    integer                                  :: i, j, l
    r(:, :) = 0.d0
    do i = 1, this%n_acc - 1
      do j = 1, size(this%npath(i)%p(:))
        l = this%npath(i)%p(j)
        r(i, abs(l)) = r(i, abs(l)) + sign(1.d0, real(l, 8))*dJdf(abs(l))
      enddo
    enddo
    do i = 1, this%nc
      !do j =1,this%np
      do j = 1, size(this%cil(i)%value, 1)
        l = this%cil(i)%value(j)
        !r(i+this%ne-1,l) = r(i+this%ne-1,j) + this%lcp(i,j)
        r(i + this%n_acc - 1, abs(l)) = r(i + this%n_acc - 1, abs(l)) + sign(1.d0, real(l, 8))
      enddo
    enddo
  end function fchi_optimizer_Aem_get_n

  function fchi_optimizer_Aem_get_w(this, Jp) result(r)
    class(fchi_optimizer_Aem), intent(in)       :: this
    real(8), dimension(this%np), intent(in)  :: Jp
    real(8), dimension(this%np)              :: r
    integer                                  :: i, j, l
    r(:) = 0.d0
    do i = 1, this%n_acc - 1
      do j = 1, size(this%npath(i)%p(:))
        l = this%npath(i)%p(j)
        r(i) = r(i) + sign(1.d0, real(l, 8))*Jp(abs(l))
      enddo
    enddo
    !do i=1,this%nh
    !l = this%nc - this%nc7 + i
    !r(l+this%ne-1) = 2d0*pi*this%hole(i,2)
    !enddo
  end function fchi_optimizer_Aem_get_w

  subroutine fchi_optimizer_Aem_adjacent_path(this)
    class(fchi_optimizer_Aem), intent(inout)        :: this
    integer i, j, k
    type(path), allocatable  ::  p(:)
    if (allocated(this%npath)) deallocate (this%npath)
    allocate(this%npath(this%n_acc))
    do i = 1, this%n_acc
      j = this%ets(i)
      call this%cpl%seek4(j, p)
      allocate(this%npath(i)%p(size(p)))
      do k = 1, size(p)
        this%npath(i)%p(k) = this%cpl%seek2(p(k))
      enddo
      !write(*,*) 'i, ets, p(k)',i, j, this%npath(i)%p
      deallocate (p)
    enddo
  end subroutine fchi_optimizer_Aem_adjacent_path

  function fchi_optimizer_Aem_calc_fdelta_chi(this, Jp, dJdf) result(r)
    class(fchi_optimizer_Aem), intent(inout)       :: this
    real(8), dimension(this%np), intent(in)  :: Jp
    real(8), dimension(this%np), intent(in)  :: dJdf
    real(8), dimension(this%np)              :: r
    real(8), dimension(this%np)              :: mat_A
    real(8), dimension(this%np, this%np)     :: mat_B
    integer                                 :: info
    integer i
    real(8), dimension(this%np)              :: test
    mat_B = this%get_n(dJdf)
    mat_A = -this%get_w(Jp) - this%jex_mat
    this%wconv(1:this%n_acc - 1) = mat_A(1:this%n_acc - 1)
    call this%my_dgetrs(mat_B, r, mat_A, info)
    !write(*, *) 'info', info
  end function fchi_optimizer_Aem_calc_fdelta_chi

  function fchi_optimizer_Aem_mix_delta_dchi(this, odf, df, ratio) &
  & result(r)
    class(fchi_optimizer_Aem), intent(in)      :: this
    real(8), dimension(this%np), intent(in) :: odf, df
    real(8), intent(in)                     :: ratio
    real(8), dimension(this%np)             :: r
    r(:) = (1d0 - ratio)*odf(:) + ratio*df(:)
  end function fchi_optimizer_Aem_mix_delta_dchi

  subroutine fchi_optimizer_Aem_set_jump(this)
    class(fchi_optimizer_Aem), intent(inout) :: this
    !real(8), dimension(this%np)          :: dchi_0
    real(8)                              :: jmp
    integer                              :: i, k, j

    if (allocated(this%jump)) then
      this%jump(:) = 0.d0
    ELSE
      allocate (this%jump(this%np), source=0.d0)
    endif
    !do i = 1, this%nc2
    !this%jump(abs(this%c2il(i)%value(5))) = 2d0*pi*this%hole(i,2)
    !enddo
    !! the jumps for c6il
    !do i = 1, this%nc6
    !this%jump(abs(this%c6il(i)%value(2))) = - this%jump(abs(this%c6il(i)%value(1))) + this%jump(abs(this%c6il(i)%value(3)))
    !write(*,*) 'path',this%c6il(i)%value(1:3)
    !enddo
    !do i = 1, this%nc4
    !this%jump(abs(this%c4il(i)%value(2))) = - this%jump(abs(this%c4il(i)%value(1))) - this%jump(abs(this%c4il(i)%value(3)))
    !!write(*,*) 'path',this%c4il(i)%value(1:3)
    !write(*,*) 'before i,this%jump(i)!,chi_f-chi_i'
    !enddo
    !do i=1, this%np
    !write(*,*) i, this%jump(i), this%chi(this%cpl%value(i)%f) - this%chi(this%cpl%value(i)%i)
    !enddo
    !print *, "----jump of chi ----" 
    do i = 1, this%np ! - 4*this%nh
      k = this%cpl%value(i)%f
      j = this%cpl%value(i)%i
      !jmp = (this%chi(k)-this%chi(j))-dchi_0(i)
      jmp = nint((this%chi(k) - this%chi(j))/(2.d0*pi))
      jmp = 2.d0*pi*jmp
      this%jump(i) = jmp
      !write(*,*) 'after i,this%jump(i)',i, jmp
      !print "(f,f6.2,3i4,'  ->',3i4)",jmp,jmp/2d0/pi,site_to_r(j,this%pshape),site_to_r(k,this%pshape)

    enddo
   ! !print *, "----end jump of chi ----" 
   ! !print *, "----jump of xi ----" 
   ! do i = 1, this%np ! - 4*this%nh
   !   k = this%cpl%value(i)%f
   !   j = this%cpl%value(i)%i
   !   !jmp = nint((this%xi(k) - this%xi(j))/(4.d0*pi))
   !   jmp = this%xi(k) - this%xi(j) - modulo(this%xi(k) - this%xi(j),2*pi)
   !   !jmp = this%xi(k) - this%xi(j)
   !   !jmp = 2.d0*pi*jmp
   !   !this%jump(i) = jmp
   !   !write(*,*) 'after i,this%jump(i)',i, jmp
   !   !print "(f,f6.2,3i4,'  ->',3i4)",jmp,jmp/2d0/pi,site_to_r(j,this%pshape),site_to_r(k,this%pshape)
   ! enddo
   ! !print *, "----end jump of xi ----" 
   ! !print *,this%xi
   ! !print *,this%xi/2/pi

   !stop
  end subroutine fchi_optimizer_Aem_set_jump

  ! protected type-bound procedures
  subroutine fchi_optimizer_Aem_set_w1(this, nc1)
    class(fchi_optimizer_Aem), intent(inout) :: this
    integer, intent(in)                  :: nc1
    if (allocated(this%w1)) deallocate (this%w1)
    allocate (this%w1(nc1), source=0)
  end subroutine fchi_optimizer_Aem_set_w1

  subroutine fchi_optimizer_Aem_set_w2(this, hole_w)
    class(fchi_optimizer_Aem), intent(inout) :: this
    integer, dimension(:), intent(in)    :: hole_w
    if (allocated(this%w2)) deallocate (this%w2)
    allocate (this%w2(size(hole_w)), source=hole_w)
  end subroutine fchi_optimizer_Aem_set_w2

  subroutine fchi_optimizer_Aem_set_w(this, w1, w2)
    class(fchi_optimizer_Aem), intent(inout) :: this
    integer, dimension(:), intent(in)    :: w1, w2
    if (allocated(this%w)) deallocate (this%w)
    allocate (this%w(size(w1) + size(w2)), source=[w1, w2])
  end subroutine fchi_optimizer_Aem_set_w

  subroutine fchi_optimizer_Aem_set_jex(this, jex)
    class(fchi_optimizer_Aem), intent(inout) :: this
    real(8), dimension(:), intent(in)    :: jex
    integer                              :: i
    if (allocated(this%jex_mat)) deallocate(this%jex_mat)
    allocate(this%jex_mat(this%np))
    this%jex_mat = 0d0
    do i=1, this%n-1
      if(this%ste(i) < 1) cycle
      this%jex_mat(this%ste(i)) = jex(i)
    end do

  end subroutine fchi_optimizer_Aem_set_jex
  
  subroutine fchi_optimizer_Aem_set_param(this, xi, chi, wf, itermax, ratio, tol, pflag, Bz)
    class(fchi_optimizer_Aem), intent(inout) :: this
    real(8), dimension(this%n), intent(in),optional :: xi, chi
    complex(8), dimension(2*this%n, 2*this%n), intent(in), optional :: wf
    integer, intent(in), optional :: itermax
    real(8), intent(in), optional :: ratio, tol ,Bz
    logical, intent(in), optional :: pflag
    if(present(xi)) this%xi = xi
    if(present(chi)) this%chi = chi
    if(present(wf)) this%wf = wf
    if(present(itermax)) this%itermax = itermax
    if(present(ratio)) this%ratio = ratio
    if(present(tol)) this%tol = tol
    if(present(pflag)) this%pflag = pflag
    if(present(Bz)) this%Bz = Bz
  end subroutine fchi_optimizer_Aem_set_param

  subroutine fchi_optimizer_Aem_set_prod_vv(this)
    class(fchi_optimizer_Aem), intent(inout) :: this
    integer k, a, i, f
    this%vv_uu = 0d0
    this%vv_dd = 0d0
    this%vv_ud = 0d0
    this%vv_du = 0d0
    do k = 1, this%np
      i = this%ste(this%cpl%value(k)%i)
      f = this%ste(this%cpl%value(k)%f)
      do a = 2*this%n_acc+1, 4*this%n_acc
        this%vv_uu(k) = this%vv_uu(k) + this%wf(4*f-1, a)*conjg(this%wf(4*i-1, a))
        this%vv_dd(k) = this%vv_dd(k) + this%wf(4*f  , a)*conjg(this%wf(4*i  , a))
        this%vv_ud(k) = this%vv_ud(k) + this%wf(4*f-1, a)*conjg(this%wf(4*i  , a))
        this%vv_du(k) = this%vv_du(k) + this%wf(4*f  , a)*conjg(this%wf(4*i-1, a))
      enddo
    enddo
    
!    do k = 1, this%np_1st
!      i = this%ste(this%cpl%value(k)%i)
!      f = this%ste(this%cpl%value(k)%f)
!      do a = 2*this%n_acc+1, 4*this%n_acc
!        this%vv_uu(k) = this%vv_uu(k) + this%wf(4*f-1, a)*conjg(this%wf(4*i-1, a))
!        this%vv_dd(k) = this%vv_dd(k) + this%wf(4*f  , a)*conjg(this%wf(4*i  , a))
!      enddo
!    enddo

!! path[3, 4] rashba
!    !do i = npt + 1, this%np - (this%nh*2)
!    do k = this%np_1st + 1, this%np_1st + this%np_rh / 2
!      i = this%ste(this%cpl%value(k)%i)
!      f = this%ste(this%cpl%value(k)%f)
!      do a = 2*this%n_acc+1, 4*this%n_acc
!          this%vv_ud(k) = this%vv_ud(k) + this%wf(4*f-1, a)*conjg(this%wf(4*i  , a))
!          this%vv_du(k) = this%vv_du(k) + this%wf(4*f  , a)*conjg(this%wf(4*i-1, a))
!      end do
!    end do

  end subroutine fchi_optimizer_Aem_set_prod_vv
  
  subroutine fchi_optimizer_Aem_calc_only_Aem_dphase(this)
    class(fchi_optimizer_Aem), intent(inout)     :: this
    integer                                    :: i
    if(allocated(this%Aem_dphase)) deallocate(this%Aem_dphase)
    allocate(this%Aem_dphase(this%np))
    do i = 1, this%np
      this%Aem_dphase(i) = - this%calc_A_em(this%cpl%value(i)%i, this%cpl%value(i)%f, this%Bz)
    end do
    !stop
  end subroutine fchi_optimizer_Aem_calc_only_Aem_dphase
  
  function fchi_optimizer_Aem_set_vector_potential(this, tau) result(u_Aeff)
    ! set vector potential
    class(fchi_optimizer_Aem), intent(in)     :: this
    real(8), dimension(this%np), intent(in) :: tau
    real(8), dimension(this%np)             :: u_Aeff
    integer                                    :: i
    do i = 1, this%np
      u_Aeff(i) = tau(i) - this%calc_A_em(this%cpl%value(i)%i, this%cpl%value(i)%f, this%Bz)
    end do
    !stop
  end function fchi_optimizer_Aem_set_vector_potential
  
  function fchi_optimizer_Aem_unset_vector_potential(this, u_Aeff) result(tau)
    ! unset vector potential (make tau from u)
    class(fchi_optimizer_Aem), intent(in)     :: this
    real(8), dimension(this%np), intent(in)    :: u_Aeff
    real(8), dimension(this%np)                :: tau
    integer                                    :: i
    do i = 1, this%np
      tau(i) = u_Aeff(i) + this%calc_A_em(this%cpl%value(i)%i, this%cpl%value(i)%f, this%Bz)
    end do
  end function fchi_optimizer_Aem_unset_vector_potential

  function fchi_optimizer_Aem_calc_A_em(this, pos_i, pos_f, Bz) result(r)
    ! calc integral_{from r_i to r_f} [ Aem * dr]
    !      ~= 1/2 (Aem(r_i) + Aem(r_f)) * (r_f  - r_i)
    class(fchi_optimizer_Aem), intent(in)       :: this
    real(8)                                        :: r
    integer, intent(in)                            :: pos_i, pos_f
    integer                                        :: pos_i_acc, pos_f_acc
    real(8), intent(in)                            :: Bz

    ! A = Bz(-y, 0, 0)
    r = (Bz*0.5d0) * ((this%y(pos_i) + this%y(pos_f)) &
                         * (this%x(pos_i) - this%x(pos_f)))

    ! A = Bz(0, x, 0)
    !r = - (Bz*0.5d0) * ((this%x(pos_i) + this%x(pos_f)) &
    !                     * (this%y(pos_i) - this%y(pos_f)))
    
    !print *, pos_i, pos_f, r
    !print *, this%x(pos_i), this%y(pos_i), this%x(pos_f), this%y(pos_f)
  end function fchi_optimizer_Aem_calc_A_em

  subroutine my_dgetrs(this, m, lambda, w, info)
    class(fchi_optimizer_Aem), intent(in) :: this
    real(8), dimension(this%np, this%np), intent(in) :: m
    real(8), dimension(this%np), intent(inout)       :: w
    real(8), dimension(this%np), intent(out)         :: lambda
    !real(8), dimension(this%np)                      :: temp
    real(8), dimension(this%np, this%np)             :: mm
    integer, dimension(this%np)                      :: ipiv
    integer                                          :: info, lwork
    integer                                           :: i
    real(8), dimension(this%np)                      :: t
    lwork = this%np**2
    mm = m
    !call check_matrix(m)
    call dgetrf(this%np, this%np, mm, this%np, ipiv, info)
    call dgetrs('N', this%np, 1, mm, this%np, ipiv, w, &
    & this%np, info)
    lambda = w
  end subroutine my_dgetrs

  subroutine fchi_optimizer_Aem_get_Jp_matrix(this, jp_mat)
    class(fchi_optimizer_Aem), intent(in) :: this
    real(8), dimension(this%n,this%n), intent(out)         :: jp_mat
    real(8), dimension(this%np)               :: Jp
    integer :: k,j,i
    real(8) :: tmp_jp
    !if(allocated(jp_mat)) deallocate(jp_mat)
    !allocate(jp_mat(this%ne,this%ne))
    jp_mat = 0d0
    !Jp = this%get_Jp(this%delta_chi)
    Jp = this%get_Jp(this%u_Aeff)
    !write(*,*) "------------------d_chi--------------"
    !write(*,*) this%delta_chi
    !write(*,*) "--------------end d_chi--------------"
    !write(*,*) "------------------Jp--------------"
    !write(*,*) Jp
    !write(*,*) "--------------end Jp--------------"
    do i = 1, this%np
      k = this%cpl%value(i)%i 
      j = this%cpl%value(i)%f 
      jp_mat(k,j) = - Jp(i)
      jp_mat(j,k) = - jp_mat(k,j)
    end do
    !write(*,*) jp_mat

    !tmp_jp = 0d0
    !do i = this%np_1st+this%np_rh+1,this%np
    !  tmp_jp = tmp_jp + Jp(i)
    !end do
    !print *,tmp_jp
  end subroutine fchi_optimizer_Aem_get_Jp_matrix

  function fchi_optimizer_Aem_approximated_J(this) result(J)
    ! this J is the approximated current in order to fit the Rondon eqn.
    class(fchi_optimizer_Aem), intent(in) :: this
    real(8)                                 :: J(this%np)
    real(8)                                 :: dJdf_0(this%np), dchi_0(this%np)
    dchi_0 = 0d0
    dJdf_0 = this%get_dJdf(dchi_0)
   !J(:) = dJdf_0(:) * this%delta_chi(:)
   J(:) = dJdf_0(:) * this%u_Aeff(:)
  end function fchi_optimizer_Aem_approximated_J

  function fchi_optimizer_Aem_J0(this) result(J)
    ! this J is the approximated current in order to fit the Rondon eqn.
    class(fchi_optimizer_Aem), intent(in) :: this
    real(8)                                 :: J(this%np)
    real(8)                                 :: dJdf_0(this%np), dchi_0(this%np)
    dchi_0 = 0d0
    dJdf_0 = this%get_Jp(dchi_0)
   J(:) = dJdf_0(:)
  end function fchi_optimizer_Aem_J0
  
  subroutine fchi_optimizer_Aem_make_jmatrix(this, jp, j_mat)
    class(fchi_optimizer_Aem), intent(in) :: this
    real(8), dimension(this%n,this%n), intent(out)        :: j_mat
    real(8), dimension(this%np), intent(in)               :: jp
    integer :: k,j,i
   
    do i = 1, this%np
      k = this%cpl%value(i)%i 
      j = this%cpl%value(i)%f 
      j_mat(k,j) = - Jp(i)
      j_mat(j,k) = - j_mat(k,j)
    end do
    
  end subroutine fchi_optimizer_Aem_make_jmatrix

 ! subroutine check_matrix(mat)
 ! USE math_mod, ONLY:my_zheev
 !   real(8), dimension(:,:), intent(in) :: mat
 !   real(8)  :: norm, det, dprod, norm_2    
 !   integer   :: n,i,j
 !   integer, dimension(:),allocatable                 :: ipiv
 !   integer                                          :: info,id
 !   complex(8), dimension(:,:), allocatable :: mat_c
 !   complex(8), dimension(:,:), allocatable :: vec
 !   real(8), dimension(:), allocatable :: value
!
 !   n = size(mat,1)
 !   allocate(ipiv(n))
 !   allocate(vec(n,n))
 !   allocate(value(n))
 !   allocate(mat_c(n,n))
 !   mat_c = mat
 !   mat_c(574,:) = mat_c(1344,:)
 !   print *, n
 !   !do i = 1, n
 !   !  !norm = sqrt(dot_product(mat(:,i),mat(:,i)))
 !   !  norm = sqrt(dot_product(mat(1:731,i),mat(1:731,i)))
 !   !  print *, i, norm
 !   !end do
 !   
 !   !do i = 1, n
 !   !    norm = sqrt(dot_product(mat(i,1:731),mat(i,1:731)))
 !   !  do j = 1, n
 !   !  dprod = dot_product(mat(i,1:731),mat(j,1:731))
 !   !  if (i /= j .and. dprod /= 0d0) then
 !   !    !print *,i,j,dprod
 !   !    norm_2 = sqrt(dot_product(mat(j,1:731),mat(j,1:731)))
 !   !    if(dprod/(norm*norm_2) > 0.5)then
 !   !      print *,i,j,dprod/(norm*norm_2)
 !   !    end if
 !   !  end if
 !   !  end do
 !   !end do
 !   !call f03baf(n, mat, n, ipiv, det, id,info)
!
 !   CALL my_zheev('l', mat_c, value, vec)
 !   do i = 1, n
 !     print *,value(i)
 !   end do
 !   stop "check fchi"
 ! end subroutine check_matrix

  ! constructor
  function new_fchi_optimizer_Aem(pshape, hole, ne, t, u, jd, lm, Bz, xi, chi, wf, eg, ne_u, ne_d) result(r)
    integer, dimension(:, :), intent(in)                   :: pshape
    integer, dimension(:, :), intent(in)                :: hole
    real(8), dimension(:), intent(in)                   :: ne, t
    real(8), intent(in)                                 :: u, jd, lm, Bz
    real(8), dimension(site_max(pshape)), intent(in) :: xi, chi
    complex(8), dimension(:, :), intent(in)      :: wf
    real(8), dimension(:), intent(in)                   :: eg, ne_u, ne_d
    type(fchi_optimizer_Aem)                                :: r
    !type(cuo2_plane), dimension(:), allocatable         :: layers
    !character(1), optional                   :: rs
    !if (present(rs)) then
      !call r%fchi_optimizer_Aem_init(pshape, hole, ne, t, u, jd, lm, xi, chi, wf )
      call r%fchi_optimizer_Aem_init(pshape, hole, ne, t, u, jd, lm, Bz, xi, chi, wf, eg, ne_u, ne_d)
    !else
    !  call r%fchi_optimizer_Aem_init(pshape, hole, t, u, jd, lm, xi, chi, wf)
    !endif
  end function new_fchi_optimizer_Aem
end module fchi_optimizer_Aem_mod
