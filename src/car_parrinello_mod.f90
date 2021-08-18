module car_parrinello_mod
  use io_mod, only:save_txt
  use math_mod, only:ui, pi,fermi,dfermi
  !use wave_function_mod, only:wave_function
  use utility_mod, only:check_sorted_eg, sort_wf, gram_schmidt, modified_gram_schmidt
  use path_list_mod, only:path, path_list
  use hop_iterator_mod, only:get_hopping_path, hop_iterator, new_hop_iterator
  !use site_list_mod
  !use parameter_mod
  use visual_mod, only : draw_phase_3d
  use base_mod, only : base
  !!!! !$ use omp_lib
  use utility_mod, only : gram_schmidt_3, gram_schmidt_3_half
  implicit none

  !logical :: g__debug_out = .false.

  !TYPE   ::  near_site
  !  INTEGER, ALLOCATABLE  :: s(:)
  !END TYPE near_site

  TYPE, extends(base) ::  car_parrinello
    !complex(8), dimension(:, :),  allocatable :: nwf, owf
    complex(8), dimension(:, :),  allocatable :: owf
    !real(8),    dimension(:),     allocatable :: eg, oeg
    real(8),    dimension(:),     allocatable :: eg
    real(8),    dimension(:),     allocatable :: Sx, Sy, Sz, ne_u, ne_d
    complex(8), dimension(:),     allocatable :: delta
    real(8)         :: h = 0.1d0,   & !time step
                  &    eta = 1.d0,  & !friction coefficient
                  &    mass = 1.d0, & !mass
                  &    tol = 1.d-7    !threshold
    real(8)         :: kbt, mu(2), tzU ! BdG 2 leyers additional parameters
    real(8)         :: c, ocoef, h2, h2t2, h2t, h2tz, h2u, h2jd, h2lm, h2tzU ! car-parinnelo parameters
    integer         :: itermax = 50000
    integer         :: delta_size
  contains
    procedure :: car_parrinello_init
    procedure :: set_initial_value => car_parrinello_set_initial_value
    procedure :: set_initial_wf => car_parrinello_set_initial_wf
    procedure :: run => car_parrinello_run
    procedure :: calc_only_wf => car_parrinello_calc_only_wf
    procedure :: calc => car_parrinello_calc
    procedure :: calc_next_wf => car_parrinello_calc_next_wf
    procedure :: calc_eigen_energy => car_parrinello_calc_eigen_energy
    procedure :: calc_spin => car_parrinello_calc_spin
    procedure :: calc_delta => car_parrinello_calc_delta
    procedure :: calc_ne => car_parrinello_calc_ne
    procedure :: calc_xi => car_parrinello_calc_xi
    procedure :: converged => car_parrinello_converged
  END TYPE car_parrinello

CONTAINS

  subroutine car_parrinello_init(this, pshape, hole, ne, t, u, jd, lm, kbt, mu, xi, chi, wf)
    class(car_parrinello), intent(inout)        :: this
    integer, dimension(:, :), intent(in)    :: pshape
    integer, dimension(:, :), intent(in) :: hole
    real(8), dimension(:), intent(in)    :: ne
    real(8), dimension(:), intent(in)    :: t
    real(8), intent(in)                  :: u, jd, lm, kbt, mu(2)
    real(8), dimension(:), intent(in)    :: xi, chi
    complex(8), dimension(:,:), intent(in), optional         :: wf
    integer i
    type(path) p

    CALL this%base_init(pshape, hole, ne, t, u, jd, lm, xi, chi, wf)
    allocate(this%eg(4*this%n_acc), source = 0d0)
    allocate(this%Sx(this%n), source = 0d0)
    allocate(this%Sy(this%n), source = 0d0)
    allocate(this%Sz(this%n), source = 0d0)
    allocate(this%ne_u(this%n))
    allocate(this%ne_d(this%n))
    this%delta_size = this%pl_1st_l(1)%length()
    allocate(this%delta(this%delta_size), source = (0d0,0d0))
    this%kbt = kbt
    this%mu = mu
    this%tzU = 2d0*this%t(3)**2d0/this%U
    this%c = 1.d0/(1d0 + this%h*this%eta*0.5d0)
    this%ocoef = this%h*this%eta*0.5d0 - 1.d0
    this%h2 = this%h**2d0
    this%h2t   = this%h2*this%t(1)
    this%h2t2  = this%h2*this%t(2)
    this%h2tz = this%h2*this%t(3)
    this%h2u = this%h2*this%u
    this%h2jd = this%jd*this%h2
    this%h2lm = this%lm*this%h2
    this%h2tzU = this%tzU*this%h2

    if(present(wf)) then
      this%wf = wf
      call this%calc_eigen_energy()
      call this%calc_spin(this%Sx, this%Sy, this%Sz, this%kbt)
      call this%calc_ne(this%ne_u, this%ne_d, this%kbt)
      call this%calc_delta(this%delta, this%kbt)
    else
      this%wf(:,:) = 0d0
      do i = 1, 4*this%n_acc
        this%wf(i,i) = 1d0
      end do

      ! initial Sx, Sy, Sz, ne_u, ne_d, delta set default values
      ! they can set any values again using subroutine "set_initial_value"
      this%Sx = 0d0
      this%Sy = 0d0
      this%Sz = 0d0
      do i = this%n_start(2), this%n_end(2) ! bulk
        this%Sx(i) = 0.5d0*cos(xi(i))
        this%Sy(i) = 0.5d0*sin(xi(i))
      end do
      this%ne_u(this%n_start(1):this%n_end(1)) = 0.5d0 * this%ne(1) / this%nl(1)
      this%ne_d(this%n_start(1):this%n_end(1)) = 0.5d0 * this%ne(1) / this%nl(1)
      !this%ne_u(this%n_start(1):this%n_end(1)) = 0.5d0
      !this%ne_d(this%n_start(1):this%n_end(1)) = 0.5d0
      this%ne_u(this%n_start(2):this%n_end(2)) = 0.5d0
      this%ne_d(this%n_start(2):this%n_end(2)) = 0.5d0
      do i = 1, this%nh
        this%ne_u(this%hole(i, 1)) = 0d0
        this%ne_d(this%hole(i, 1)) = 0d0
      end do
      do i = 1, this%pl_1st_l(1)%length()
        p = this%pl_1st_l(1)%index(i)
        if(this%y(p%i) == this%y(p%f)) then
          this%delta(i) = 0.2d0 ! x-dir
        else
          this%delta(i) = -0.2d0 ! y-dir
        end if
      end do

      call this%calc_eigen_energy()
    end if

    !allocate(this%oeg(4*this%n_acc), source = this%eg)
    !allocate(this%nwf(4*this%n_acc, 4*this%n_acc), source = (0d0,0d0))
    allocate(this%owf(4*this%n_acc, 4*this%n_acc), source = this%wf)

  end subroutine car_parrinello_init
  
  subroutine car_parrinello_set_initial_value(this, Sx, Sy, Sz, delta, ne_u, ne_d, wf)
    class(car_parrinello), intent(inout)        :: this
    real(8),    dimension(this%n), intent(in), optional :: Sx, Sy, Sz, ne_u, ne_d
    complex(8), dimension(this%delta_size), intent(in), optional :: delta
    complex(8), dimension(4*this%n_acc, 4*this%n_acc), intent(in), optional :: wf
    if(present(Sx))  this%Sx = Sx
    if(present(Sy))  this%Sy = Sy
    if(present(Sz))  this%Sz = Sz
    if(present(delta))  this%delta = delta
    if(present(ne_u))  this%ne_u = ne_u
    if(present(ne_d))  this%ne_d = ne_d
    if(present(wf))  this%wf = wf
    call this%calc_eigen_energy()
  end subroutine car_parrinello_set_initial_value
  
  subroutine car_parrinello_set_initial_wf(this, wf, eg)
    class(car_parrinello), intent(inout)        :: this
    complex(8), dimension(4*this%n_acc, 4*this%n_acc), intent(in) :: wf
    real(8),    dimension(4*this%n_acc), intent(in) :: eg
    this%wf = wf
    this%eg = eg
    call this%calc_ne(this%ne_u, this%ne_d, this%kbt)
    call this%calc_delta(this%delta, this%kbt)
    call this%calc_spin(this%Sx, this%Sy, this%Sz, this%kbt)
  end subroutine car_parrinello_set_initial_wf
  
  subroutine car_parrinello_run(this)
    class(car_parrinello), intent(inout)        :: this
    !call this%calc(0d0, 0d0, 0d0) ! only wf
    !call this%calc(0.1d0, 0d0, 0d0) ! ne
    !call this%calc(0d0, 0.1d0, 0d0) ! delta
    !call this%calc(0d0, 0d0, 0.1d0) ! spin
    call this%calc(1d0, 1d0, 1d0) ! all
  end subroutine car_parrinello_run
  
  subroutine car_parrinello_calc_only_wf(this)
    class(car_parrinello), intent(inout)        :: this
    call this%calc(0d0, 0d0, 0d0)
  end subroutine car_parrinello_calc_only_wf
  
  subroutine car_parrinello_calc(this, ratio_ne, ratio_delta, ratio_spin)
    class(car_parrinello), intent(inout)        :: this
    real(8),    intent(in)  :: ratio_ne, ratio_delta, ratio_spin
    !complex(8), dimension(4*this%n_acc,4*this%n_acc) :: nwf
    complex(8), dimension(:, :), allocatable :: nwf
    real(8),    dimension(4*this%n_acc)     :: oeg
    real(8),    dimension(this%n)           :: oSx, oSy, oSz, old_ne_u, old_ne_d
    real(8),    dimension(this%n)           :: nSx, nSy, nSz, new_ne_u, new_ne_d
    complex(8), dimension(this%delta_size)  :: old_delta, new_delta
    LOGICAL conv_eg, conv_ne_u, conv_ne_d, conv_delta_real, conv_delta_img, conv_Sx, conv_Sy, conv_Sz
    integer :: iter, i, a
    character(1) :: cheg, chne, chdelta, chspin
    allocate(nwf(4*this%n_acc, 4*this%n_acc))

    this%owf = this%wf
    oeg = this%eg
    oSx = this%Sx
    oSy = this%Sy
    oSz = this%Sz
    old_ne_u = this%ne_u
    old_ne_d = this%ne_d
    old_delta = this%delta

    if(ratio_ne == 0d0) chne = '-'
    if(ratio_delta == 0d0) chdelta = '-'
    if(ratio_spin == 0d0) chspin = '-'

    do iter = 1, this%itermax
      WRITE (*, '(a6,i7,a2,i7,a2)') '- iter', iter, ' /', this%itermax, ' -'
      
      call this%calc_next_wf(nwf)
      
      this%owf = this%wf
      this%wf = nwf
      if (.false.) then !!!! othogonization for all 
        !call gram_schmidt(this%wf)
        !call modified_gram_schmidt(this%wf)
        call gram_schmidt_3(this%wf)
        call this%calc_eigen_energy()
      else  !!!! othogonization for half
        !call gram_schmidt(this%wf, 2*this%n_acc)
        call gram_schmidt_3_half(this%wf)
        call this%calc_eigen_energy()
        do a = 1, 2*this%n_acc      ! set E_{-n} & phi_{-n}
          this%eg(4*this%n_acc - a + 1) = - this%eg(a)
          do i = 1, this%n_acc
            this%wf(4*i - 3, 4*this%n_acc - a + 1) = - conjg(this%wf(4*i - 1, a))
            this%wf(4*i - 2, 4*this%n_acc - a + 1) =   conjg(this%wf(4*i    , a))
            this%wf(4*i - 1, 4*this%n_acc - a + 1) = - conjg(this%wf(4*i - 3, a))
            this%wf(4*i    , 4*this%n_acc - a + 1) =   conjg(this%wf(4*i - 2, a))
          end do
        end do
      end if
      
     !if(iter <= 1000 .and. mod(iter, 100) == 0 .and. not(check_sorted_eg(eg))) then
     if(mod(iter, 1000) == 0 .and. (.not. check_sorted_eg(this%eg))) then
     !if(not(check_sorted_eg(eg))) then
        print *, "sorted"
        call sort_wf(this%eg, this%wf, oeg, this%owf)
      end if
      
     ! call sort_wf(this%eg, this%wf, oeg, this%owf)



      conv_eg = this%converged(this%eg, oeg, this%tol)
      oeg = this%eg
      if(conv_eg) then
        cheg = 'T'
      else
        cheg = 'F'
      end if
      if(ratio_ne > 0d0) then
        old_ne_u  =this%ne_u
        old_ne_d  =this%ne_d
        call this%calc_ne(new_ne_u, new_ne_d, this%kbt)
        conv_ne_u = this%converged(new_ne_u, old_ne_u, this%tol)
        conv_ne_d = this%converged(new_ne_d, old_ne_d, this%tol)
        this%ne_u = ratio_ne*new_ne_u  + (1d0-ratio_ne)*old_ne_u
        this%ne_d = ratio_ne*new_ne_d  + (1d0-ratio_ne)*old_ne_d
        if(conv_ne_u .and. conv_ne_d) then
          chne = 'T'
        else
          chne = 'F'
        end if
      else
        conv_ne_u = .true.
        conv_ne_d = .true.
      end if
      if(ratio_delta > 0d0) then
        old_delta = this%delta
        call this%calc_delta(new_delta, this%kbt)
        conv_delta_real = this%converged(real(new_delta, 8), real(old_delta, 8), this%tol)
        conv_delta_img  = this%converged(aimag(new_delta), aimag(old_delta), this%tol)
        !!!conv_delta_abs  = this%converged(abs(delta), abs(old_delta), this%tol)
        this%delta = ratio_delta*new_delta  + (1d0-ratio_delta)*old_delta
        if(conv_delta_real .and. conv_delta_img) then
          chdelta = 'T'
        else
          chdelta = 'F'
        end if
      else
        conv_delta_real = .true.
        conv_delta_img = .true.
      end if
      if(ratio_spin > 0d0) then
        oSx = this%Sx
        oSy = this%Sy
        oSz = this%Sz
        call this%calc_spin(nSx, nSy, nSz, this%kbt)
        conv_Sx = this%converged(nSx, oSx, this%tol)
        conv_Sy = this%converged(nSy, oSy, this%tol)
        conv_Sz = this%converged(nSz, oSz, this%tol)
        this%Sx = ratio_spin*nSx  + (1d0-ratio_spin)*oSx
        this%Sy = ratio_spin*nSy  + (1d0-ratio_spin)*oSy
        this%Sz = ratio_spin*nSz  + (1d0-ratio_spin)*oSz
        if(conv_Sx .and. conv_Sy .and. conv_Sz) then
          chspin = 'T'
        else
          chspin = 'F'
        end if
      else
        conv_Sx = .true.
        conv_Sy = .true.
        conv_Sz = .true.
      end if
      
      WRITE (*, *) 'sum eigen energy', SUM(this%eg(2*this%n_acc+1:4*this%n_acc))
     print *, "bulk ne : ", sum(this%ne_u(this%nl(1)+1:)) + sum(this%ne_d(this%nl(1)+1:)), &
            & "surface ne : ", sum(this%ne_u(1:this%nl(1))) + sum(this%ne_d(1:this%nl(1)))
      print *, "ave abs(delta) : ", sum(abs(this%delta))/this%delta_size
      print "(x,3(a,f))", "Sx:", sum(abs(this%Sx)), "  Sy:", sum(abs(this%Sy)), "  Sz:", sum(abs(this%Sz))
      !WRITE (*, *) 'sum eigen energy', SUM(this%eg(1:4*this%n_acc))
      !WRITE (*, *) 'sum d energy', SUM(ABS((eg(:) - oeg(:))))
      !WRITE (*, *) 'max abs(delta)', maxval(ABS((delta(:))))
      !write(*,"(a,f10.5,a,f10.5,a,f8.5)") "ne; up : " &
      !      &  ,sum(this%ne_u), "  down : ", sum(this%ne_d), "  max ne : ", maxval(this%ne_u + this%ne_d)
      !print *, conv_eg, conv_ne_u, conv_ne_d, conv_delta_real, conv_delta_img, conv_Sx, conv_Sy, conv_Sz
      !print "(4(a,l,2x))", "eg:", conv_eg, "ne:", conv_ne_u .and. conv_ne_d, &
      !                &    "delta:", conv_delta_real .and. conv_delta_img, &
      !                &    "spin:", conv_Sx .and. conv_Sy .and. conv_Sz
      print "(x,4(a,a,2x))", "eg:", cheg, "ne:", chne, "delta:", chdelta, "spin:", chspin
      
      if( iter >= 50 .and. ALL([conv_eg, conv_ne_u, conv_ne_d, &
                              & conv_delta_real, conv_delta_img, conv_Sx, conv_Sy, conv_Sz])) then
        exit    
      end if
    end do

    if(.not.check_sorted_eg(this%eg)) then
      print *, "sorted"
      call sort_wf(this%eg, this%wf, oeg, this%owf)
    end if


  end subroutine car_parrinello_calc
  
  subroutine car_parrinello_calc_next_wf(this, nwf)
    class(car_parrinello), intent(in)        :: this
    complex(8), dimension(4*this%n_acc, 4*this%n_acc), intent(out) :: nwf
    integer :: r, i, j, k, is, js, ks, l
    nwf(:, :) = 0d0

    !$omp parallel do
    do r = 1, 4*this%n_acc
      !!!!!!!!!!!!!!!!!!    common    !!!!!!!!!!!!!!!!!!
      do i = 1, this%n_acc
        nwf(4*i - 3, r) = ( 2.d0 + this%h2*this%eg(r))*this%wf(4*i - 3, r) + this%ocoef*this%owf(4*i - 3, r)
        nwf(4*i - 2, r) = ( 2.d0 + this%h2*this%eg(r))*this%wf(4*i - 2, r) + this%ocoef*this%owf(4*i - 2 ,r)
        nwf(4*i - 1, r) = ( 2.d0 + this%h2*this%eg(r))*this%wf(4*i - 1, r) + this%ocoef*this%owf(4*i - 1, r)
        nwf(4*i    , r) = ( 2.d0 + this%h2*this%eg(r))*this%wf(4*i    , r) + this%ocoef*this%owf(4*i    , r)
      end do
      !!!!!!!!!!!!!!!!!!    end common    !!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!    surface    !!!!!!!!!!!!!!!!!!
      do i = this%n_acc_start(1), this%n_acc_end(1)
        is = this%ets(i)
        !nwf(4*i - 3, r) = ( 2.d0 + this%h2*(this%eg(r) + this%mu(1)) )*this%wf(4*i - 3, r) + this%ocoef*this%owf(4*i - 3, r)
        !nwf(4*i - 2, r) = ( 2.d0 + this%h2*(this%eg(r) + this%mu(1)) )*this%wf(4*i - 2, r) + this%ocoef*this%owf(4*i - 2 ,r)
        !nwf(4*i - 1, r) = ( 2.d0 + this%h2*(this%eg(r) - this%mu(1)) )*this%wf(4*i - 1, r) + this%ocoef*this%owf(4*i - 1, r)
        !nwf(4*i    , r) = ( 2.d0 + this%h2*(this%eg(r) - this%mu(1)) )*this%wf(4*i    , r) + this%ocoef*this%owf(4*i    , r)
        nwf(4*i-3, r) = nwf(4*i-3, r) + this%h2*this%mu(1)*this%wf(4*i-3, r)
        nwf(4*i-2, r) = nwf(4*i-2, r) + this%h2*this%mu(1)*this%wf(4*i-2, r)
        nwf(4*i-1, r) = nwf(4*i-1, r) - this%h2*this%mu(1)*this%wf(4*i-1, r)
        nwf(4*i  , r) = nwf(4*i  , r) - this%h2*this%mu(1)*this%wf(4*i  , r)
        do l = 1, size(this%ns_1st(is)%s)
          js = this%ns_1st(is)%s(l)
          j = this%ste(js)
          nwf(4*i - 3, r) = nwf(4*i - 3, r) + this%h2t*this%wf(4*j - 3, r)*(1.d0 - this%ne_d(is))*(1.d0 - this%ne_d(js))
          nwf(4*i - 2, r) = nwf(4*i - 2, r) + this%h2t*this%wf(4*j - 2, r)*(1.d0 - this%ne_u(is))*(1.d0 - this%ne_u(js))
          nwf(4*i - 1, r) = nwf(4*i - 1, r) - this%h2t*this%wf(4*j - 1, r)*(1.d0 - this%ne_d(is))*(1.d0 - this%ne_d(js))
          nwf(4*i    , r) = nwf(4*i    , r) - this%h2t*this%wf(4*j    , r)*(1.d0 - this%ne_u(is))*(1.d0 - this%ne_u(js))
        end do ! l
        do l = 1, size(this%ns_2nd(is)%s)
          js = this%ns_2nd(is)%s(l)
          j = this%ste(js)
          nwf(4*i - 3, r) = nwf(4*i - 3, r) + this%h2t2*this%wf(4*j - 3, r)*(1.d0 - this%ne_d(is))*(1.d0 - this%ne_d(js))
          nwf(4*i - 2, r) = nwf(4*i - 2, r) + this%h2t2*this%wf(4*j - 2, r)*(1.d0 - this%ne_u(is))*(1.d0 - this%ne_u(js))
          nwf(4*i - 1, r) = nwf(4*i - 1, r) - this%h2t2*this%wf(4*j - 1, r)*(1.d0 - this%ne_d(is))*(1.d0 - this%ne_d(js))
          nwf(4*i    , r) = nwf(4*i    , r) - this%h2t2*this%wf(4*j    , r)*(1.d0 - this%ne_u(is))*(1.d0 - this%ne_u(js))
        end do ! l    
      end do!i
      do i = 1, this%pl_1st_l(1)%length()
        ks = this%pl_1st_l(1)%value(i)%i
        js = this%pl_1st_l(1)%value(i)%f
        j = this%ste(js)
        k = this%ste(ks)
        nwf(4*k - 3, r) = nwf(4*k - 3, r)  -       this%delta(i) * this%h2*this%wf(4*j    , r)
        nwf(4*k - 2, r) = nwf(4*k - 2, r)  -       this%delta(i) * this%h2*this%wf(4*j - 1, r)
        nwf(4*k - 1, r) = nwf(4*k - 1, r)  - conjg(this%delta(i))* this%h2*this%wf(4*j - 2, r)
        nwf(4*k    , r) = nwf(4*k    , r)  - conjg(this%delta(i))* this%h2*this%wf(4*j - 3, r) 
        nwf(4*j - 3, r) = nwf(4*j - 3, r)  -       this%delta(i) * this%h2*this%wf(4*k    , r)
        nwf(4*j - 2, r) = nwf(4*j - 2, r)  -       this%delta(i) * this%h2*this%wf(4*k - 1, r)
        nwf(4*j - 1, r) = nwf(4*j - 1, r)  - conjg(this%delta(i))* this%h2*this%wf(4*k - 2, r)
        nwf(4*j    , r) = nwf(4*j    , r)  - conjg(this%delta(i))* this%h2*this%wf(4*k - 3, r) 
      end do ! i
      !!!!!!!!!!!!!!    end surface    !!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!    bulk    !!!!!!!!!!!!!!!!
      do i = this%n_acc_start(2), this%n_acc_end(2)
        is = this%ets(i)
        nwf(4*i-3, r) = nwf(4*i-3, r) + this%h2*this%mu(2)*this%wf(4*i-3, r)
        nwf(4*i-2, r) = nwf(4*i-2, r) + this%h2*this%mu(2)*this%wf(4*i-2, r)
        nwf(4*i-1, r) = nwf(4*i-1, r) - this%h2*this%mu(2)*this%wf(4*i-1, r)
        nwf(4*i  , r) = nwf(4*i  , r) - this%h2*this%mu(2)*this%wf(4*i  , r)
        nwf(4*i-3, r) = nwf(4*i-3, r) - this%h2U*(-2d0/3d0*this%Sz(is)+0.5d0)*this%wf(4*i-3, r)
        nwf(4*i-2, r) = nwf(4*i-2, r) - this%h2U*( 2d0/3d0*this%Sz(is)+0.5d0)*this%wf(4*i-2, r)
        nwf(4*i-1, r) = nwf(4*i-1, r) + this%h2U*(-2d0/3d0*this%Sz(is)+0.5d0)*this%wf(4*i-1, r)
        nwf(4*i  , r) = nwf(4*i  , r) + this%h2U*( 2d0/3d0*this%Sz(is)+0.5d0)*this%wf(4*i  , r)
        nwf(4*i-3, r) = nwf(4*i-3, r) + this%h2U*2d0/3d0*(this%Sx(is)-ui*this%Sy(is))*this%wf(4*i-2, r)
        nwf(4*i-2, r) = nwf(4*i-2, r) + this%h2U*2d0/3d0*(this%Sx(is)+ui*this%Sy(is))*this%wf(4*i-3, r)
        nwf(4*i-1, r) = nwf(4*i-1, r) + this%h2U*2d0/3d0*(this%Sx(is)+ui*this%Sy(is))*this%wf(4*i  , r)
        nwf(4*i  , r) = nwf(4*i  , r) + this%h2U*2d0/3d0*(this%Sx(is)-ui*this%Sy(is))*this%wf(4*i-1, r)
        do l = 1, size(this%ns_1st(is)%s)
          js = this%ns_1st(is)%s(l)
          j = this%ste(js)
          nwf(4*i - 3, r) = nwf(4*i - 3, r) + this%h2t*this%wf(4*j - 3, r)
          nwf(4*i - 2, r) = nwf(4*i - 2, r) + this%h2t*this%wf(4*j - 2, r)
          nwf(4*i - 1, r) = nwf(4*i - 1, r) - this%h2t*this%wf(4*j - 1, r)
          nwf(4*i    , r) = nwf(4*i    , r) - this%h2t*this%wf(4*j    , r)
        end do ! l
        do l = 1, size(this%ns_2nd(is)%s)
          js = this%ns_2nd(is)%s(l)
          j = this%ste(js)
          nwf(4*i - 3, r) = nwf(4*i - 3, r) + this%h2t2*this%wf(4*j - 3, r)
          nwf(4*i - 2, r) = nwf(4*i - 2, r) + this%h2t2*this%wf(4*j - 2, r)
          nwf(4*i - 1, r) = nwf(4*i - 1, r) - this%h2t2*this%wf(4*j - 1, r)
          nwf(4*i    , r) = nwf(4*i    , r) - this%h2t2*this%wf(4*j    , r)
        end do ! l    
      end do!i
      do i = 1, this%pl_h%length()
        js = this%pl_h%value(i)%i
        ks = this%pl_h%value(i)%f
        j = this%ste(js)
        k = this%ste(ks)
        nwf(4*j-3, r) = nwf(4*j-3, r) - 0.5d0*this%h2jd*this%Sz(ks)*this%wf(4*j-3, r)
        nwf(4*j-2, r) = nwf(4*j-2, r) + 0.5d0*this%h2jd*this%Sz(ks)*this%wf(4*j-2, r)
        nwf(4*j-1, r) = nwf(4*j-1, r) + 0.5d0*this%h2jd*this%Sz(ks)*this%wf(4*j-1, r)
        nwf(4*j  , r) = nwf(4*j  , r) - 0.5d0*this%h2jd*this%Sz(ks)*this%wf(4*j  , r)
        nwf(4*j-3, r) = nwf(4*j-3, r) - 0.5d0*this%h2jd*(this%Sx(ks) - ui*this%Sy(ks))*this%wf(4*j-2, r)
        nwf(4*j-2, r) = nwf(4*j-2, r) - 0.5d0*this%h2jd*(this%Sx(ks) + ui*this%Sy(ks))*this%wf(4*j-3, r)
        nwf(4*j-1, r) = nwf(4*j-1, r) - 0.5d0*this%h2jd*(this%Sx(ks) + ui*this%Sy(ks))*this%wf(4*j  , r)
        nwf(4*j  , r) = nwf(4*j  , r) - 0.5d0*this%h2jd*(this%Sx(ks) - ui*this%Sy(ks))*this%wf(4*j-1, r)
        
        nwf(4*k-3, r) = nwf(4*k-3, r) - 0.5d0*this%h2jd*this%Sz(js)*this%wf(4*k-3, r)
        nwf(4*k-2, r) = nwf(4*k-2, r) + 0.5d0*this%h2jd*this%Sz(js)*this%wf(4*k-2, r)
        nwf(4*k-1, r) = nwf(4*k-1, r) + 0.5d0*this%h2jd*this%Sz(js)*this%wf(4*k-1, r)
        nwf(4*k  , r) = nwf(4*k  , r) - 0.5d0*this%h2jd*this%Sz(js)*this%wf(4*k  , r)
        nwf(4*k-3, r) = nwf(4*k-3, r) - 0.5d0*this%h2jd*(this%Sx(js) - ui*this%Sy(js))*this%wf(4*k-2, r)
        nwf(4*k-2, r) = nwf(4*k-2, r) - 0.5d0*this%h2jd*(this%Sx(js) + ui*this%Sy(js))*this%wf(4*k-3, r)
        nwf(4*k-1, r) = nwf(4*k-1, r) - 0.5d0*this%h2jd*(this%Sx(js) + ui*this%Sy(js))*this%wf(4*k  , r)
        nwf(4*k  , r) = nwf(4*k  , r) - 0.5d0*this%h2jd*(this%Sx(js) - ui*this%Sy(js))*this%wf(4*k-1, r)
      end do ! i
      !!!!!!!!!!!!!!    end bulk    !!!!!!!!!!!!!!!!!!
      
      !!!!!!!!!!!!!!!!!!    inter-layer    !!!!!!!!!!!!!!!!
      ! pl_z_l(i) is z-dir path list between ith layer and (i+1)th layer.
      ! here, pl_z_l(1) is interlayer path ( between 1st(surfaca) and 2nd(bulk) ).
      !do i = this%n_acc_start(1), this%n_acc_end(1)
      !  is = this%ets(i)
      !  do l = 1, size(this%ns_1st(is)%s)
      !    js = this%ns_1st(is)%s(l)
      !    j = this%ste(js)
      
      !    nwf(4*i - 3, r) = nwf(4*i - 3, r) + h2t2*wf(4*j - 3, r)*(1.d0 - ne_d(i))*(1.d0 - ne_d(j))
      !    nwf(4*i - 2, r) = nwf(4*i - 2, r) + h2t2*wf(4*j - 2, r)*(1.d0 - ne_u(i))*(1.d0 - ne_u(j))
      !    nwf(4*i - 1, r) = nwf(4*i - 1, r) - h2t2*wf(4*j - 1, r)*(1.d0 - ne_d(i))*(1.d0 - ne_d(j))
      !    nwf(4*i    , r) = nwf(4*i    , r) - h2t2*wf(4*j    , r)*(1.d0 - ne_u(i))*(1.d0 - ne_u(j))

      !end do
      !do i = this%n_acc_start(2), this%n_acc_end(2)
      !  
      !end do

      do i = 1, this%pl_z_l(1)%length()
        js = this%pl_z_l(1)%value(i)%i
        ks = this%pl_z_l(1)%value(i)%f
        j = this%ste(js)
        k = this%ste(ks)
        nwf(4*j-3, r) = nwf(4*j-3, r) + this%h2tz*(1d0-this%ne_d(js))*(1d0-this%ne_d(ks))*this%wf(4*k-3, r)
        nwf(4*j-2, r) = nwf(4*j-2, r) + this%h2tz*(1d0-this%ne_u(js))*(1d0-this%ne_u(ks))*this%wf(4*k-2, r)
        nwf(4*j-1, r) = nwf(4*j-1, r) - this%h2tz*(1d0-this%ne_d(js))*(1d0-this%ne_d(ks))*this%wf(4*k-1, r)
        nwf(4*j  , r) = nwf(4*j  , r) - this%h2tz*(1d0-this%ne_u(js))*(1d0-this%ne_u(ks))*this%wf(4*k  , r)
        nwf(4*k-3, r) = nwf(4*k-3, r) + this%h2tz*(1d0-this%ne_d(js))*(1d0-this%ne_d(ks))*this%wf(4*j-3, r)
        nwf(4*k-2, r) = nwf(4*k-2, r) + this%h2tz*(1d0-this%ne_u(js))*(1d0-this%ne_u(ks))*this%wf(4*j-2, r)
        nwf(4*k-1, r) = nwf(4*k-1, r) - this%h2tz*(1d0-this%ne_d(js))*(1d0-this%ne_d(ks))*this%wf(4*j-1, r)
        nwf(4*k  , r) = nwf(4*k  , r) - this%h2tz*(1d0-this%ne_u(js))*(1d0-this%ne_u(ks))*this%wf(4*j  , r)
        
        nwf(4*j-3, r) = nwf(4*j-3, r) - this%h2tzU*(this%Sz(ks)-0.5d0*(this%ne_u(ks)+this%ne_d(ks)))*this%wf(4*j-3, r)
        nwf(4*j-2, r) = nwf(4*j-2, r) + this%h2tzU*(this%Sz(ks)+0.5d0*(this%ne_u(ks)+this%ne_d(ks)))*this%wf(4*j-2, r)
        nwf(4*j-1, r) = nwf(4*j-1, r) + this%h2tzU*(this%Sz(ks)-0.5d0*(this%ne_u(ks)+this%ne_d(ks)))*this%wf(4*j-1, r)
        nwf(4*j  , r) = nwf(4*j  , r) - this%h2tzU*(this%Sz(ks)+0.5d0*(this%ne_u(ks)+this%ne_d(ks)))*this%wf(4*j  , r)
        nwf(4*j-3, r) = nwf(4*j-3, r) - this%h2tzU*(this%Sx(ks) - ui*this%Sy(ks))*this%wf(4*j-2, r)
        nwf(4*j-2, r) = nwf(4*j-2, r) - this%h2tzU*(this%Sx(ks) + ui*this%Sy(ks))*this%wf(4*j-3, r)
        nwf(4*j-1, r) = nwf(4*j-1, r) - this%h2tzU*(this%Sx(ks) + ui*this%Sy(ks))*this%wf(4*j  , r)
        nwf(4*j  , r) = nwf(4*j  , r) - this%h2tzU*(this%Sx(ks) - ui*this%Sy(ks))*this%wf(4*j-1, r)
        
        nwf(4*k-3, r) = nwf(4*k-3, r) - this%h2tzU*(this%Sz(js)-0.5d0*(this%ne_u(js)+this%ne_d(js)))*this%wf(4*k-3, r)
        nwf(4*k-2, r) = nwf(4*k-2, r) + this%h2tzU*(this%Sz(js)+0.5d0*(this%ne_u(js)+this%ne_d(js)))*this%wf(4*k-2, r)
        nwf(4*k-1, r) = nwf(4*k-1, r) + this%h2tzU*(this%Sz(js)-0.5d0*(this%ne_u(js)+this%ne_d(js)))*this%wf(4*k-1, r)
        nwf(4*k  , r) = nwf(4*k  , r) - this%h2tzU*(this%Sz(js)+0.5d0*(this%ne_u(js)+this%ne_d(js)))*this%wf(4*k  , r)
        nwf(4*k-3, r) = nwf(4*k-3, r) - this%h2tzU*(this%Sx(js) - ui*this%Sy(js))*this%wf(4*k-2, r)
        nwf(4*k-2, r) = nwf(4*k-2, r) - this%h2tzU*(this%Sx(js) + ui*this%Sy(js))*this%wf(4*k-3, r)
        nwf(4*k-1, r) = nwf(4*k-1, r) - this%h2tzU*(this%Sx(js) + ui*this%Sy(js))*this%wf(4*k  , r)
        nwf(4*k  , r) = nwf(4*k  , r) - this%h2tzU*(this%Sx(js) - ui*this%Sy(js))*this%wf(4*k-1, r)
      end do ! i
      !!!!!!!!!!!!!!!!!!    end inter-layer    !!!!!!!!!!!!!!!!
    end do!r
    !$omp end parallel do
  nwf = this%c*nwf
  end subroutine car_parrinello_calc_next_wf
  
  subroutine car_parrinello_calc_eigen_energy(this)
    class(car_parrinello), intent(inout)        :: this
    integer :: r, i, j, k, is, js, ks
    real(8) :: site_ne
    this%eg(:) = 0d0
    !$omp parallel do
    do r = 1, 4*this%n_acc
    !!!!!!!!!!!!!!!!!!    surface    !!!!!!!!!!!!!!!!!!
      do i = this%n_acc_start(1), this%n_acc_end(1)
        this%eg(r) = this%eg(r) - this%mu(1)*conjg(this%wf(4*i-3, r))*this%wf(4*i-3, r)
        this%eg(r) = this%eg(r) - this%mu(1)*conjg(this%wf(4*i-2, r))*this%wf(4*i-2, r)
        this%eg(r) = this%eg(r) + this%mu(1)*conjg(this%wf(4*i-1, r))*this%wf(4*i-1, r)
        this%eg(r) = this%eg(r) + this%mu(1)*conjg(this%wf(4*i  , r))*this%wf(4*i  , r)
      end do
      do i = 1, this%pl_1st_l(1)%length()
        js = this%pl_1st_l(1)%value(i)%i
        ks = this%pl_1st_l(1)%value(i)%f
        j = this%ste(js) ! in the 1st layer (surface), js = j, ks = k
        k = this%ste(ks)
        this%eg(r) = this%eg(r) - this%t(1)*(1d0-this%ne_d(js))*(1d0-this%ne_d(ks))*conjg(this%wf(4*j-3, r))*this%wf(4*k-3, r)
        this%eg(r) = this%eg(r) - this%t(1)*(1d0-this%ne_u(js))*(1d0-this%ne_u(ks))*conjg(this%wf(4*j-2, r))*this%wf(4*k-2, r)
        this%eg(r) = this%eg(r) + this%t(1)*(1d0-this%ne_d(js))*(1d0-this%ne_d(ks))*conjg(this%wf(4*j-1, r))*this%wf(4*k-1, r)
        this%eg(r) = this%eg(r) + this%t(1)*(1d0-this%ne_u(js))*(1d0-this%ne_u(ks))*conjg(this%wf(4*j  , r))*this%wf(4*k  , r)
        this%eg(r) = this%eg(r) - this%t(1)*(1d0-this%ne_d(js))*(1d0-this%ne_d(ks))*conjg(this%wf(4*k-3, r))*this%wf(4*j-3, r)
        this%eg(r) = this%eg(r) - this%t(1)*(1d0-this%ne_u(js))*(1d0-this%ne_u(ks))*conjg(this%wf(4*k-2, r))*this%wf(4*j-2, r)
        this%eg(r) = this%eg(r) + this%t(1)*(1d0-this%ne_d(js))*(1d0-this%ne_d(ks))*conjg(this%wf(4*k-1, r))*this%wf(4*j-1, r)
        this%eg(r) = this%eg(r) + this%t(1)*(1d0-this%ne_u(js))*(1d0-this%ne_u(ks))*conjg(this%wf(4*k  , r))*this%wf(4*j  , r)
      end do
      do i = 1, this%pl_2nd_l(1)%length()
        js = this%pl_2nd_l(1)%value(i)%i
        ks = this%pl_2nd_l(1)%value(i)%f
        j = this%ste(js) ! in the 1st layer (surface), js = j, ks = k
        k = this%ste(ks)
        this%eg(r) = this%eg(r) - this%t(2)*(1d0-this%ne_d(js))*(1d0-this%ne_d(ks))*conjg(this%wf(4*j-3, r))*this%wf(4*k-3, r)
        this%eg(r) = this%eg(r) - this%t(2)*(1d0-this%ne_u(js))*(1d0-this%ne_u(ks))*conjg(this%wf(4*j-2, r))*this%wf(4*k-2, r)
        this%eg(r) = this%eg(r) + this%t(2)*(1d0-this%ne_d(js))*(1d0-this%ne_d(ks))*conjg(this%wf(4*j-1, r))*this%wf(4*k-1, r)
        this%eg(r) = this%eg(r) + this%t(2)*(1d0-this%ne_u(js))*(1d0-this%ne_u(ks))*conjg(this%wf(4*j  , r))*this%wf(4*k  , r)
        this%eg(r) = this%eg(r) - this%t(2)*(1d0-this%ne_d(js))*(1d0-this%ne_d(ks))*conjg(this%wf(4*k-3, r))*this%wf(4*j-3, r)
        this%eg(r) = this%eg(r) - this%t(2)*(1d0-this%ne_u(js))*(1d0-this%ne_u(ks))*conjg(this%wf(4*k-2, r))*this%wf(4*j-2, r)
        this%eg(r) = this%eg(r) + this%t(2)*(1d0-this%ne_d(js))*(1d0-this%ne_d(ks))*conjg(this%wf(4*k-1, r))*this%wf(4*j-1, r)
        this%eg(r) = this%eg(r) + this%t(2)*(1d0-this%ne_u(js))*(1d0-this%ne_u(ks))*conjg(this%wf(4*k  , r))*this%wf(4*j  , r)
      end do
      do i = 1, this%pl_1st_l(1)%length()
        ! here, delta's index corresponds to pl_1st_l(1)'s index
        js = this%pl_1st_l(1)%value(i)%i
        ks = this%pl_1st_l(1)%value(i)%f
        j = this%ste(js) ! in the 1st layer ( surface ), js = j, ks = k
        k = this%ste(ks)
        this%eg(r) = this%eg(r) +       this%delta(i) *conjg(this%wf(4*j-3, r))*this%wf(4*k  , r)
        this%eg(r) = this%eg(r) +       this%delta(i) *conjg(this%wf(4*j-2, r))*this%wf(4*k-1, r)
        this%eg(r) = this%eg(r) + conjg(this%delta(i))*conjg(this%wf(4*j-1, r))*this%wf(4*k-2, r)
        this%eg(r) = this%eg(r) + conjg(this%delta(i))*conjg(this%wf(4*j  , r))*this%wf(4*k-3, r)
        this%eg(r) = this%eg(r) +       this%delta(i) *conjg(this%wf(4*k-3, r))*this%wf(4*j  , r)
        this%eg(r) = this%eg(r) +       this%delta(i) *conjg(this%wf(4*k-2, r))*this%wf(4*j-1, r)
        this%eg(r) = this%eg(r) + conjg(this%delta(i))*conjg(this%wf(4*k-1, r))*this%wf(4*j-2, r)
        this%eg(r) = this%eg(r) + conjg(this%delta(i))*conjg(this%wf(4*k  , r))*this%wf(4*j-3, r)
      end do
      !!!!!!!!!!!!!!    end surface    !!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!    bulk    !!!!!!!!!!!!!!!!
      do i = this%n_acc_start(2), this%n_acc_end(2)
        is = this%ets(i)
        this%eg(r) = this%eg(r) + this%U*(-2d0/3d0*this%Sz(is)+0.5d0)*conjg(this%wf(4*i-3, r))*this%wf(4*i-3, r)
        this%eg(r) = this%eg(r) + this%U*( 2d0/3d0*this%Sz(is)+0.5d0)*conjg(this%wf(4*i-2, r))*this%wf(4*i-2, r)
        this%eg(r) = this%eg(r) - this%U*(-2d0/3d0*this%Sz(is)+0.5d0)*conjg(this%wf(4*i-1, r))*this%wf(4*i-1, r)
        this%eg(r) = this%eg(r) - this%U*( 2d0/3d0*this%Sz(is)+0.5d0)*conjg(this%wf(4*i  , r))*this%wf(4*i  , r)
        this%eg(r) = this%eg(r) - this%U*2d0/3d0*(this%Sx(is)-ui*this%Sy(is))*conjg(this%wf(4*i-3, r))*this%wf(4*i-2, r)
        this%eg(r) = this%eg(r) - this%U*2d0/3d0*(this%Sx(is)+ui*this%Sy(is))*conjg(this%wf(4*i-2, r))*this%wf(4*i-3, r)
        this%eg(r) = this%eg(r) - this%U*2d0/3d0*(this%Sx(is)+ui*this%Sy(is))*conjg(this%wf(4*i-1, r))*this%wf(4*i  , r)
        this%eg(r) = this%eg(r) - this%U*2d0/3d0*(this%Sx(is)-ui*this%Sy(is))*conjg(this%wf(4*i  , r))*this%wf(4*i-1, r)
      end do
      do i = this%n_acc_start(2), this%n_acc_end(2)
        this%eg(r) = this%eg(r) - this%mu(2)*conjg(this%wf(4*i-3, r))*this%wf(4*i-3, r)
        this%eg(r) = this%eg(r) - this%mu(2)*conjg(this%wf(4*i-2, r))*this%wf(4*i-2, r)
        this%eg(r) = this%eg(r) + this%mu(2)*conjg(this%wf(4*i-1, r))*this%wf(4*i-1, r)
        this%eg(r) = this%eg(r) + this%mu(2)*conjg(this%wf(4*i  , r))*this%wf(4*i  , r)
      end do
      do i = 1, this%pl_1st_l(2)%length()
        js = this%pl_1st_l(2)%value(i)%i
        ks = this%pl_1st_l(2)%value(i)%f
        j = this%ste(js)
        k = this%ste(ks)
        this%eg(r) = this%eg(r) - this%t(1)*conjg(this%wf(4*j-3, r))*this%wf(4*k-3, r)
        this%eg(r) = this%eg(r) - this%t(1)*conjg(this%wf(4*j-2, r))*this%wf(4*k-2, r)
        this%eg(r) = this%eg(r) + this%t(1)*conjg(this%wf(4*j-1, r))*this%wf(4*k-1, r)
        this%eg(r) = this%eg(r) + this%t(1)*conjg(this%wf(4*j  , r))*this%wf(4*k  , r)
        this%eg(r) = this%eg(r) - this%t(1)*conjg(this%wf(4*k-3, r))*this%wf(4*j-3, r)
        this%eg(r) = this%eg(r) - this%t(1)*conjg(this%wf(4*k-2, r))*this%wf(4*j-2, r)
        this%eg(r) = this%eg(r) + this%t(1)*conjg(this%wf(4*k-1, r))*this%wf(4*j-1, r)
        this%eg(r) = this%eg(r) + this%t(1)*conjg(this%wf(4*k  , r))*this%wf(4*j  , r)
      end do
      do i = 1, this%pl_2nd_l(2)%length()
        js = this%pl_2nd_l(2)%value(i)%i
        ks = this%pl_2nd_l(2)%value(i)%f
        j = this%ste(js)
        k = this%ste(ks)
        this%eg(r) = this%eg(r) - this%t(2)*conjg(this%wf(4*j-3, r))*this%wf(4*k-3, r)
        this%eg(r) = this%eg(r) - this%t(2)*conjg(this%wf(4*j-2, r))*this%wf(4*k-2, r)
        this%eg(r) = this%eg(r) + this%t(2)*conjg(this%wf(4*j-1, r))*this%wf(4*k-1, r)
        this%eg(r) = this%eg(r) + this%t(2)*conjg(this%wf(4*j  , r))*this%wf(4*k  , r)
        this%eg(r) = this%eg(r) - this%t(2)*conjg(this%wf(4*k-3, r))*this%wf(4*j-3, r)
        this%eg(r) = this%eg(r) - this%t(2)*conjg(this%wf(4*k-2, r))*this%wf(4*j-2, r)
        this%eg(r) = this%eg(r) + this%t(2)*conjg(this%wf(4*k-1, r))*this%wf(4*j-1, r)
        this%eg(r) = this%eg(r) + this%t(2)*conjg(this%wf(4*k  , r))*this%wf(4*j  , r)
      end do
      do i = 1, this%pl_h%length()
        js = this%pl_h%value(i)%i
        ks = this%pl_h%value(i)%f
        j = this%ste(js)
        k = this%ste(ks)
        this%eg(r) = this%eg(r) + 0.5d0*this%jd*this%Sz(ks)*conjg(this%wf(4*j-3, r))*this%wf(4*j-3, r)
        this%eg(r) = this%eg(r) - 0.5d0*this%jd*this%Sz(ks)*conjg(this%wf(4*j-2, r))*this%wf(4*j-2, r)
        this%eg(r) = this%eg(r) - 0.5d0*this%jd*this%Sz(ks)*conjg(this%wf(4*j-1, r))*this%wf(4*j-1, r)
        this%eg(r) = this%eg(r) + 0.5d0*this%jd*this%Sz(ks)*conjg(this%wf(4*j  , r))*this%wf(4*j  , r)
        this%eg(r) = this%eg(r) + 0.5d0*this%jd*(this%Sx(ks) - ui*this%Sy(ks))*conjg(this%wf(4*j-3, r))*this%wf(4*j-2, r)
        this%eg(r) = this%eg(r) + 0.5d0*this%jd*(this%Sx(ks) + ui*this%Sy(ks))*conjg(this%wf(4*j-2, r))*this%wf(4*j-3, r)
        this%eg(r) = this%eg(r) + 0.5d0*this%jd*(this%Sx(ks) + ui*this%Sy(ks))*conjg(this%wf(4*j-1, r))*this%wf(4*j  , r)
        this%eg(r) = this%eg(r) + 0.5d0*this%jd*(this%Sx(ks) - ui*this%Sy(ks))*conjg(this%wf(4*j  , r))*this%wf(4*j-1, r)
        
        this%eg(r) = this%eg(r) + 0.5d0*this%jd*this%Sz(js)*conjg(this%wf(4*k-3, r))*this%wf(4*k-3, r)
        this%eg(r) = this%eg(r) - 0.5d0*this%jd*this%Sz(js)*conjg(this%wf(4*k-2, r))*this%wf(4*k-2, r)
        this%eg(r) = this%eg(r) - 0.5d0*this%jd*this%Sz(js)*conjg(this%wf(4*k-1, r))*this%wf(4*k-1, r)
        this%eg(r) = this%eg(r) + 0.5d0*this%jd*this%Sz(js)*conjg(this%wf(4*k  , r))*this%wf(4*k  , r)
        this%eg(r) = this%eg(r) + 0.5d0*this%jd*(this%Sx(js) - ui*this%Sy(js))*conjg(this%wf(4*k-3, r))*this%wf(4*k-2, r)
        this%eg(r) = this%eg(r) + 0.5d0*this%jd*(this%Sx(js) + ui*this%Sy(js))*conjg(this%wf(4*k-2, r))*this%wf(4*k-3, r)
        this%eg(r) = this%eg(r) + 0.5d0*this%jd*(this%Sx(js) + ui*this%Sy(js))*conjg(this%wf(4*k-1, r))*this%wf(4*k  , r)
        this%eg(r) = this%eg(r) + 0.5d0*this%jd*(this%Sx(js) - ui*this%Sy(js))*conjg(this%wf(4*k  , r))*this%wf(4*k-1, r)
      end do
      !!!!!!!!!!!!!!    end bulk    !!!!!!!!!!!!!!!!!!
      
      !!!!!!!!!!!!!!!!!!    interlayer    !!!!!!!!!!!!!!!!
      ! pl_z_l(i) is z-dir path list between ith layer and (i+1)th layer.
      ! here, pl_z_l(1) is interlayer path ( between 1st(surfaca) and 2nd(bulk) ).
      do i = 1, this%pl_z_l(1)%length()
        js = this%pl_z_l(1)%value(i)%i
        ks = this%pl_z_l(1)%value(i)%f
        j = this%ste(js)
        k = this%ste(ks)
        this%eg(r) = this%eg(r) - this%t(3)*(1d0-this%ne_d(js))*(1d0-this%ne_d(ks))*conjg(this%wf(4*j-3, r))*this%wf(4*k-3, r)
        this%eg(r) = this%eg(r) - this%t(3)*(1d0-this%ne_u(js))*(1d0-this%ne_u(ks))*conjg(this%wf(4*j-2, r))*this%wf(4*k-2, r)
        this%eg(r) = this%eg(r) + this%t(3)*(1d0-this%ne_d(js))*(1d0-this%ne_d(ks))*conjg(this%wf(4*j-1, r))*this%wf(4*k-1, r)
        this%eg(r) = this%eg(r) + this%t(3)*(1d0-this%ne_u(js))*(1d0-this%ne_u(ks))*conjg(this%wf(4*j  , r))*this%wf(4*k  , r)
        this%eg(r) = this%eg(r) - this%t(3)*(1d0-this%ne_d(js))*(1d0-this%ne_d(ks))*conjg(this%wf(4*k-3, r))*this%wf(4*j-3, r)
        this%eg(r) = this%eg(r) - this%t(3)*(1d0-this%ne_u(js))*(1d0-this%ne_u(ks))*conjg(this%wf(4*k-2, r))*this%wf(4*j-2, r)
        this%eg(r) = this%eg(r) + this%t(3)*(1d0-this%ne_d(js))*(1d0-this%ne_d(ks))*conjg(this%wf(4*k-1, r))*this%wf(4*j-1, r)
        this%eg(r) = this%eg(r) + this%t(3)*(1d0-this%ne_u(js))*(1d0-this%ne_u(ks))*conjg(this%wf(4*k  , r))*this%wf(4*j  , r)

        site_ne = this%ne_u(ks)+this%ne_d(ks)
        this%eg(r) = this%eg(r) + this%tzU*(this%Sz(ks)-0.5d0*site_ne)*conjg(this%wf(4*j-3, r))*this%wf(4*j-3, r)
        this%eg(r) = this%eg(r) - this%tzU*(this%Sz(ks)+0.5d0*site_ne)*conjg(this%wf(4*j-2, r))*this%wf(4*j-2, r)
        this%eg(r) = this%eg(r) - this%tzU*(this%Sz(ks)-0.5d0*site_ne)*conjg(this%wf(4*j-1, r))*this%wf(4*j-1, r)
        this%eg(r) = this%eg(r) + this%tzU*(this%Sz(ks)+0.5d0*site_ne)*conjg(this%wf(4*j  , r))*this%wf(4*j  , r)
        this%eg(r) = this%eg(r) + this%tzU*(this%Sx(ks) - ui*this%Sy(ks))*conjg(this%wf(4*j-3, r))*this%wf(4*j-2, r)
        this%eg(r) = this%eg(r) + this%tzU*(this%Sx(ks) + ui*this%Sy(ks))*conjg(this%wf(4*j-2, r))*this%wf(4*j-3, r)
        this%eg(r) = this%eg(r) + this%tzU*(this%Sx(ks) + ui*this%Sy(ks))*conjg(this%wf(4*j-1, r))*this%wf(4*j  , r)
        this%eg(r) = this%eg(r) + this%tzU*(this%Sx(ks) - ui*this%Sy(ks))*conjg(this%wf(4*j  , r))*this%wf(4*j-1, r)
        
        site_ne = this%ne_u(js)+this%ne_d(js)
        this%eg(r) = this%eg(r) + this%tzU*(this%Sz(js)-0.5d0*site_ne)*conjg(this%wf(4*k-3, r))*this%wf(4*k-3, r)
        this%eg(r) = this%eg(r) - this%tzU*(this%Sz(js)+0.5d0*site_ne)*conjg(this%wf(4*k-2, r))*this%wf(4*k-2, r)
        this%eg(r) = this%eg(r) - this%tzU*(this%Sz(js)-0.5d0*site_ne)*conjg(this%wf(4*k-1, r))*this%wf(4*k-1, r)
        this%eg(r) = this%eg(r) + this%tzU*(this%Sz(js)+0.5d0*site_ne)*conjg(this%wf(4*k  , r))*this%wf(4*k  , r)
        this%eg(r) = this%eg(r) + this%tzU*(this%Sx(js) - ui*this%Sy(js))*conjg(this%wf(4*k-3, r))*this%wf(4*k-2, r)
        this%eg(r) = this%eg(r) + this%tzU*(this%Sx(js) + ui*this%Sy(js))*conjg(this%wf(4*k-2, r))*this%wf(4*k-3, r)
        this%eg(r) = this%eg(r) + this%tzU*(this%Sx(js) + ui*this%Sy(js))*conjg(this%wf(4*k-1, r))*this%wf(4*k  , r)
        this%eg(r) = this%eg(r) + this%tzU*(this%Sx(js) - ui*this%Sy(js))*conjg(this%wf(4*k  , r))*this%wf(4*k-1, r)
      end do
      !!!!!!!!!!!!!!   end  interlayer    !!!!!!!!!!!!!!!!
    end do ! r
    !$omp end parallel do
  end subroutine car_parrinello_calc_eigen_energy

  subroutine car_parrinello_calc_spin(this, Sx, Sy, Sz, kbt)
    class(car_parrinello), intent(inout)        :: this
    real(8), dimension(this%n), intent(out) :: Sx, Sy, Sz
    real(8), intent(in) :: kbt
    integer :: i, is, a
    complex(8) :: uu, ud, vu, vd
    Sx = 0d0
    Sy = 0d0
    Sz = 0d0
    do i = 1, this%n_acc
      is = this%ets(i)
      do a = 2*this%n_acc + 1, 4*this%n_acc 
        uu = this%wf(4*i-3, a)
        ud = this%wf(4*i-2, a)
        vu = this%wf(4*i-1, a)
        vd = this%wf(4*i  , a)

        !!!! T= 0
        !Sx(is) = Sx(is) - 0.5d0*(vu*conjg(vd) + conjg(vu)*vd)
        !Sy(is) = Sy(is) + 0.5d0*ui*(vu*conjg(vd) - conjg(vu)*vd)
        !Sz(is) = Sz(is) + 0.5d0*(vu*conjg(vu)-vd*conjg(vd))
        
        Sx(is) = Sx(is) + 0.5d0*(uu*conjg(ud) + conjg(uu)*ud)*fermi( this%eg(a), kbt) &
               &        - 0.5d0*(vu*conjg(vd) + conjg(vu)*vd)*fermi(-this%eg(a), kbt) 
        Sy(is) = Sy(is) - 0.5d0*ui*(uu*conjg(ud) - conjg(uu)*ud)*fermi( this%eg(a), kbt) &
               &        + 0.5d0*ui*(vu*conjg(vd) - conjg(vu)*vd)*fermi(-this%eg(a), kbt) 
        Sz(is) = Sz(is) + 0.5d0*(uu*conjg(uu) - ud*conjg(ud))*fermi( this%eg(a), kbt) &
               &        + 0.5d0*(vu*conjg(vu) - vd*conjg(vd))*fermi(-this%eg(a), kbt) 
      end do
    end do
  end subroutine car_parrinello_calc_spin
  
  subroutine car_parrinello_calc_delta(this, delta, kbt)
    class(car_parrinello), intent(inout)        :: this
    complex(8), dimension(this%delta_size), intent(out) :: delta
    real(8), intent(in) :: kbt
    integer :: i, j, k, a
    delta  = 0.d0
    do a = 2*this%n_acc+1, 4*this%n_acc
      do k = 1, this%pl_1st_l(1)%length()
        i = this%pl_1st_l(1)%value(k)%i
        j = this%pl_1st_l(1)%value(k)%f
        ! jd : 2*t(1)**2 / U
        delta(k) = delta(k) + (this%jd*(this%wf(4*i-3, a)*conjg(this%wf(4*j  , a)) &
                  &                   + this%wf(4*j-2, a)*conjg(this%wf(4*i-1, a)) &
                  &                   + this%wf(4*i-2, a)*conjg(this%wf(4*j-1, a)) &
                  &                   + this%wf(4*j-3, a)*conjg(this%wf(4*i  , a))))*tanh(this%eg(a)/kbT*0.5d0)
      end do!k
    end do!a
  end subroutine car_parrinello_calc_delta
  
  subroutine car_parrinello_calc_ne(this, ne_u, ne_d, kbt)
    class(car_parrinello), intent(in)        :: this
    real(8), dimension(this%n), intent(out) :: ne_u, ne_d
    real(8), intent(in) :: kbt
    integer :: i, is, a
    ne_u = 0.d0
    ne_d = 0.d0
    do i = 1, this%n_acc
      is = this%ets(i)
      do a = 2*this%n_acc+1, 4*this%n_acc
        ne_u(is) = ne_u(is) + fermi( this%eg(a), kbt)*conjg(this%wf(4*i-3, a))*this%wf(4*i-3, a) &
                  &         + fermi(-this%eg(a), kbt)*conjg(this%wf(4*i-1, a))*this%wf(4*i-1, a)
        ne_d(is) = ne_d(is) + fermi( this%eg(a), kbt)*conjg(this%wf(4*i-2, a))*this%wf(4*i-2, a) &
                  &         + fermi(-this%eg(a), kbt)*conjg(this%wf(4*i  , a))*this%wf(4*i  , a) 
      end do!a
    end do!i
  end subroutine car_parrinello_calc_ne

  subroutine car_parrinello_calc_xi(this)
    class(car_parrinello), intent(inout)        :: this
  end subroutine car_parrinello_calc_xi
  
  function car_parrinello_converged(this, v, ov, tol) result(r)
    class(car_parrinello), intent(in) :: this
    real(8), dimension(:), intent(in) :: v, ov
    real(8), intent(in) :: tol
    logical :: r
    integer i
    do i = 1, size(v)
      if(abs(ov(i))  <  tol) cycle
      if(ABS((v(i) - ov(i))/ov(i)) >= tol) then
        r = .false.
        return
      end if
    end do
    r = .true.
  end function car_parrinello_converged

  FUNCTION new_car_parrinello(pshape, hole, ne, t, u, jd, lm, kbt, mu, xi, chi, wf) RESULT(r)
    TYPE(car_parrinello)               :: r
    INTEGER, DIMENSION(:,:), INTENT(in)    :: pshape
    INTEGER, DIMENSION(:, :), INTENT(in) :: hole
    REAL(8), DIMENSION(:), INTENT(in)    :: ne
    REAL(8), DIMENSION(:), INTENT(in)    :: t
    REAL(8), INTENT(in)                  :: u, jd, lm, kbt, mu(2)
    REAL(8), DIMENSION(:), INTENT(in)    :: xi, chi
    COMPLEX(8), DIMENSION(:,:), INTENT(in), optional   :: wf
    CALL r%car_parrinello_init(pshape, hole, ne, t, u, jd, lm, kbt, mu, xi, chi, wf)
  END FUNCTION new_car_parrinello

END MODULE car_parrinello_mod