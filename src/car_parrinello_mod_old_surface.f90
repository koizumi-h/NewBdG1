MODULE car_parrinello_mod
  USE io_mod, ONLY:save_txt
  USE math_mod, ONLY:ui, pi,fermi,dfermi
  USE wave_function_mod, ONLY:wave_function
  USE utility_mod, ONLY:check_sorted_eg, sort_wf, gram_schmidt
  USE path_list_mod, ONLY:path, path_list
  USE hop_iterator_mod, ONLY:get_hopping_path, hop_iterator, new_hop_iterator
  USE site_list_mod
  use parameter_mod
  use visual_mod ,only :draw_phase_3D
  use base_mod
  !$ use omp_lib
  IMPLICIT NONE

  logical :: g__debug_out = .false.

  !TYPE   ::  near_site
  !  INTEGER, ALLOCATABLE  :: s(:)
  !END TYPE near_site

  TYPE, extends(wave_function) ::  car_parrinello
    COMPLEX(8), DIMENSION(:, :), ALLOCATABLE :: nwf, owf
    REAL(8)                               :: hbar
    REAL(8), DIMENSION(:), ALLOCATABLE     :: density, Sx, Sy, Sz
    !TYPE(path_list)                       :: pl,pl_2nd
    !TYPE(near_site), ALLOCATABLE           :: ns(:),ns_2nd(:)
    REAL(8)                               :: h = 0.1d0,& !time step
         &eta = 1.d0,& !friction coefficient
         &mu = 1.d0,& !mass
         &tol = 1.d-5 !threshold
    INTEGER                               :: itermax = 10000 !1000000
    LOGICAL                               :: temp
  CONTAINS
    PROCEDURE :: car_parrinello_init
    !PROCEDURE :: adjacent_site => car_parrinello_adjacent_site
    !PROCEDURE :: adjacent_2nd_site => car_parrinello_adjacent_2nd_site
    PROCEDURE :: calc => car_parrinello_calc
    ! PROCEDURE :: car_parrinello_new_calc !!!
    PROCEDURE :: random_wf => car_parrinello_random_wf
    PROCEDURE :: gs => car_parrinello_gs
    PROCEDURE :: normal_vec => car_parrinello_normal_vec
    PROCEDURE :: get_wave_functions => car_parrinello_get_wave_functions
    PROCEDURE :: get_energies => car_parrinello_get_energies
    PROCEDURE :: calc_xi => car_parrinello_calc_xi
    PROCEDURE :: is_converged => car_parrinello_is_converged
    PROCEDURE :: calc_surface
    PROCEDURE :: calc_bulk
    PROCEDURE :: calc_surface_general => car_parinello_calc_surface_general
    PROCEDURE :: calc_mu_surface => car_parinello_calc_mu_surface_at_ne
  END TYPE car_parrinello

CONTAINS
  !  subroutine car_parrinello_()
  !    CLASS(car_parrinello), intent(inout)        :: this
  !  end subroutine car_parrinello_q
  !  function car_parrinello_() result(r)
  !    CLASS(car_parrinello), intent(in)        :: this
  !  end function car_parrinello_

  SUBROUTINE car_parrinello_init(this, pshape, hole, ne, t, u, jd, lm, xi, chi, wf)
    CLASS(car_parrinello), INTENT(inout)        :: this
    INTEGER, DIMENSION(:, :), INTENT(in)    :: pshape
    INTEGER, DIMENSION(:, :), INTENT(in) :: hole
    REAL(8), DIMENSION(:), INTENT(in)    :: ne
    REAL(8), DIMENSION(:), INTENT(in)    :: t
    REAL(8), INTENT(in)                  :: u, jd, lm
    REAL(8), DIMENSION(:), INTENT(in)    :: xi, chi
    !COMPLEX(8), DIMENSION(2*site_max(pshape), 2*site_max(pshape)), &
    !     & INTENT(in), OPTIONAL                             :: wf
    COMPLEX(8), DIMENSION(:,:), INTENT(in), OPTIONAL         :: wf
    integer i
    !    write(*,*) 'hole', hole(:,:)
    !    write(*,*) 't', t(:)
    !    write(*,*) 'u, jd, lm:', u, jd, lm
    !    write(*,*) 'xi:', xi(:)
    !IF (PRESENT(wf)) then
    !!  CALL this%wave_function_init_cp(pshape, hole, t, u, jd, lm, xi, chi, wf)
    !  CALL this%wave_function_init_cp(pshape, hole, t, u, jd, lm, xi, chi)
    !  !CALL this%set_eg()
    !else
      CALL this%wave_function_init_cp(pshape, hole, ne, t, u, jd, lm, xi, chi)
    !end if
    !CALL this%pl%alloc_value()
    !CALL this%pl%extend(get_hopping_path(this%pshape, [1, 2], this%hole(:, 1)))
    !CALL this%pl_2nd%alloc_value()
    !CALL this%pl_2nd%extend(get_hopping_path(this%pshape, [5, 6], this%hole(:, 1)))
    CALL this%set_site_to_e()
    !CALL this%adjacent_site()
    !CALL this%adjacent_2nd_site()
    ALLOCATE (this%nwf(4*this%n, 4*this%n))
    ALLOCATE (this%owf(4*this%n, 4*this%n))
    ALLOCATE (this%density(this%n))
    ALLOCATE (this%Sx(this%n))
    ALLOCATE (this%Sy(this%n))
    ALLOCATE (this%Sz(this%n))
    this%temp = .FALSE.
  END SUBROUTINE car_parrinello_init

  !SUBROUTINE car_parrinello_adjacent_site(this)
  !  CLASS(car_parrinello), INTENT(inout)        :: this
  !  INTEGER i, j, k
  !  TYPE(path), ALLOCATABLE  ::  p(:)
  !  IF (ALLOCATED(this%ns)) DEALLOCATE (this%ns)
  !  ALLOCATE(this%ns(sum(int(this%ne))))
  !  DO i = 1, sum(int(this%ne))
  !    j = this%ets(i)
  !    CALL this%pl%seek3(j, p)
  !    ALLOCATE(this%ns(i)%s(SIZE(p)))
  !    DO k = 1, SIZE(p)
  !      IF (p(k)%i == j) THEN
  !        this%ns(i)%s(k) = p(k)%f
  !      ELSEIF (p(k)%f == j) THEN
  !        this%ns(i)%s(k) = p(k)%i
  !      ELSE
  !        STOP'Error'
  !      ENDIF
  !    ENDDO
  !    DEALLOCATE (p)
  !  ENDDO
  !END SUBROUTINE car_parrinello_adjacent_site

  !SUBROUTINE car_parrinello_adjacent_2nd_site(this)
  !  CLASS(car_parrinello), INTENT(inout)        :: this
  !  INTEGER i, j, k
  !  TYPE(path), ALLOCATABLE  ::  p(:)
  !  IF (ALLOCATED(this%ns_2nd)) DEALLOCATE (this%ns_2nd)
  !  ALLOCATE(this%ns_2nd(sum(int(this%ne))))
  !  DO i = 1, sum(int(this%ne))
  !    j = this%ets(i)
  !    CALL this%pl_2nd%seek3(j, p)
  !    ALLOCATE(this%ns_2nd(i)%s(SIZE(p)))
  !    DO k = 1, SIZE(p)
  !      IF (p(k)%i == j) THEN
  !        this%ns_2nd(i)%s(k) = p(k)%f
  !      ELSEIF (p(k)%f == j) THEN
  !        this%ns_2nd(i)%s(k) = p(k)%i
  !      ELSE
  !        STOP'Error'
  !      ENDIF
  !    ENDDO
  !    DEALLOCATE (p)
  !  ENDDO
  !  !!check pathlist_2nd
  !  !    do i=1, this%pl_2nd%length()
  !  !      write(*, *) '****************',i, this%pl_2nd%value(i)
  !  !    enddo
  !  !!!!!stop
  !  !! check ns_2nd
  !  !  DO i = 1, sum(int(this%ne))
  !  !    do k=1, size(this%ns_2nd(i)%s)
  !  !      write(*, *) '****************',i,"%s(",k,")", this%ns_2nd(i)%s(k)
  !  !    enddo
  !  !  enddo
  !  !  stop
  !END SUBROUTINE car_parrinello_adjacent_2nd_site

  SUBROUTINE car_parrinello_calc(this)
  use visual_mod
    CLASS(car_parrinello), INTENT(inout)      :: this
    REAL(8)      :: c, ocoef, h2, h2t, h2t_z, h2u, h2jd, h2lm, mu, kbT, mu_l, mu_m, mu_h, ne_l, ne_m, ne_h, fe_e
    REAL(8), DIMENSION(this%n)                :: density_up, density_down, density, Sx, Sy, Sz, oSx, oSy, oSz
    REAL(8), DIMENSION(this%n)                :: ne_u, ne_d
    REAL(8), DIMENSION(4*this%n)    :: oeg
    real(8), parameter         :: ratio_s = 1.d0!0.1d0
    INTEGER      :: iter, r, i, j, k, l, f , z, a
    TYPE(path), ALLOCATABLE    :: p(:, :)
    LOGICAL                    :: flag, print_flg
    complex(8) delta(size(this%pl_1st%value))
    character    file*2
    REAL(8), DIMENSION(4*this%n)    :: test_eg
    complex(8), DIMENSION(4*this%n, 4*this%n)    :: test_wf
  !call draw_phase_3D('test1_xi.gp', 44, this%pshape, this%xi, this%hole(:, 1))
  !print *, "saved to [test1_xi.gp]"
    !OPEN (20, file='new_cpeg_h1.dat', status='replace')
    flag = .FALSE.
    !ALLOCATE (hop, source=new_hop_iterator(this%pshape))
    !CALL hop%set_path_hole(this%hole(:, 1)) !prepare <i,j> for calc Jd
    !!test
    !    do i=1, size(hop%pl%value(:))
    !    write(*, *) 'car_parrinello_calc',i, hop%pl%value(i)
    !    enddo
    !    stop
    CALL this%set_eg()
    this%owf = this%wf
    this%nwf = this%wf
    c = 1.d0/(1 + this%h*this%eta*0.5d0)
    ocoef = this%h*this%eta*0.5d0 - 1.d0
    h2 = this%h**2
    h2t = h2*this%t(1)
    !h2t_2nd = h2*this%t(2)
    h2t_z = h2*this%t(3)
    h2u = h2*this%u
    h2jd = this%jd*h2
    h2lm = this%lm*h2
    Sx = 0.d0
    Sy = 0.d0
    Sz = 0.d0
    ! set hole's neighbors
    ALLOCATE (p(this%nh, 4))
    DO i = 1, this%nh
      z = site_to_z(this%hole(i, 1),this%pshape)
      p(i, 1)%i = this%hole(i, 1) - this%nx(z)
      p(i, 1)%f = this%hole(i, 1) + 1
      p(i, 2)%i = this%hole(i, 1) - 1
      p(i, 2)%f = this%hole(i, 1) + this%nx(z)
      p(i, 3)%i = this%hole(i, 1) - this%nx(z)
      p(i, 3)%f = this%hole(i, 1) - 1
      p(i, 4)%i = this%hole(i, 1) + 1
      p(i, 4)%f = this%hole(i, 1) + this%nx(z)
    ENDDO
  !  do i = 1,4*this%n
  !  write(*,*)DOT_PRODUCT(this%wf(:, i), this%wf(:, i))
  !  end do
  !  stop









    if(.false.)then
    DO iter = 1, this%itermax
      if(iter < 10000 .or. (iter > 10000 .and. mod(iter,100) == 0)) then
        print_flg = .true.
      else
        print_flg = .false.
      end if
      if(print_flg) then
        WRITE (*, '(a6,i7,a2,i7,a2)') '-iter', iter, ' /', this%itermax, ' -'
      end if
      oSx = Sx
      oSy = Sy
      oSz = Sz
      oeg = this%eg(1:4*this%n)
      CALL this%calc_S(density, Sx, Sy, Sz)
      CALL this%calc_delta(density_up, density_down, delta)
      ! mix S and oS
      Sx = ratio_s*Sx + (1d0-ratio_s)*oSx
      Sy = ratio_s*Sy + (1d0-ratio_s)*oSy
      Sz = ratio_s*Sz + (1d0-ratio_s)*oSz
      !!!check xi
      if(.false.) then
        if(iter <= 50000) then
          call draw_phase_3D("cp/out_"// string(iter) // ".gp", 71 ,this%pshape,this%calc_xi(),this%hole(:,1))
          print *,"cp/out_"// string(iter) // ".gp"
        end if
      end if
      !!end check xi 
      DO r = 1, 4*this%n
        DO i = 1, this%n
          !Coulomb and friction
    !      this%nwf(4*i - 3, r) = (2.d0 + h2*this%eg(r) - h2u*(0.5d0 - (2d0/3d0)*Sz(k)))*this%wf(4*i - 3, r)&
    !           &+ h2u*((2d0/3d0)*(Sx(k) - ui*Sy(k)))*this%wf(4*i, r) + ocoef*this%owf(4*i - 1, r)
    !      this%nwf(4*i - 2, r) = (2.d0 + h2*this%eg(r) - h2u*(0.5d0 + (2d0/3d0)*Sz(k)))*this%wf(4*i - 2, r)&
    !           &+ h2u*((2d0/3d0)*(Sx(k) + ui*Sy(k)))*this%wf(4*i - 3, r) + ocoef*this%owf(4*i, r)
    !      this%nwf(4*i - 1, r) = (2.d0 + h2*this%eg(r) - h2u*(0.5d0 - (2d0/3d0)*Sz(k)))*this%wf(4*i - 1, r)&
    !           &+ h2u*((2d0/3d0)*(Sx(k) - ui*Sy(k)))*this%wf(4*i, r) + ocoef*this%owf(4*i - 1, r)
    !      this%nwf(4*i, r) = (2.d0 + h2*this%eg(r) - h2u*(0.5d0 + (2d0/3d0)*Sz(k)))*this%wf(4*i, r)&
    !          &+ h2u*((2d0/3d0)*(Sx(k) + ui*Sy(k)))*this%wf(4*i - 1, r) + ocoef*this%owf(4*i, r)
          this%nwf(4*i - 3, r) = (2.d0 + h2*this%eg(r) )*this%wf(4*i - 3, r) + ocoef*this%owf(4*i - 3, r)
          this%nwf(4*i - 2, r) = (2.d0 + h2*this%eg(r) )*this%wf(4*i - 2, r) + ocoef*this%owf(4*i - 2 , r)
          this%nwf(4*i - 1, r) = (2.d0 + h2*this%eg(r) )*this%wf(4*i - 1, r) + ocoef*this%owf(4*i - 1, r)
          this%nwf(4*i, r) = (2.d0 + h2*this%eg(r) )*this%wf(4*i, r) + ocoef*this%owf(4*i, r)
          this%nwf(4*i - 3, r) =  this%nwf(4*i - 3, r) + h2*mu*this%wf(4*i - 3, r)
          this%nwf(4*i - 2, r) =  this%nwf(4*i - 2, r) + h2*mu*this%wf(4*i - 2, r)
          this%nwf(4*i - 1, r) =  this%nwf(4*i - 1, r) - h2*mu*this%wf(4*i - 1, r)
          this%nwf(4*i, r) =  this%nwf(4*i , r) - h2*mu*this%wf(4*i , r)
          ! hopping nearest neighbor
          DO j = 1, SIZE(this%ns_1st(i)%s)
          !  write(*,*)i,this%ns_1st(i)%s(j)
          !  write(*,*)this%n,i,this%ns_1st(i)%s
            this%nwf(4*i - 3, r) = this%nwf(4*i - 3, r) + h2t*this%wf(4*this%ns_1st(i)%s(j) - 3, r)
            this%nwf(4*i - 2, r) = this%nwf(4*i - 2, r) + h2t*this%wf(4*this%ns_1st(i)%s(j) - 2, r)
            this%nwf(4*i - 1, r) = this%nwf(4*i - 1, r) + h2t*this%wf(4*this%ns_1st(i)%s(j) - 1, r)
            this%nwf(4*i, r) = this%nwf(4*i, r) + h2t*this%wf(4*this%ns_1st(i)%s(j), r)
          ENDDO ! j
          ! hopping 2nd nearest neighbor
          !DO j = 1, SIZE(this%ns_2nd(i)%s)
          !  this%nwf(2*k - 1, r) = this%nwf(2*k - 1, r) + h2t_2nd*this%wf(2*this%ns_2nd(i)%s(j) - 1, r)
          !  this%nwf(2*k, r) = this%nwf(2*k, r) + h2t_2nd*this%wf(2*this%ns_2nd(i)%s(j), r)
          !ENDDO  ! j
          ! hopping of z direction
          DO j = 1, SIZE(this%ns_z(i)%s)
            this%nwf(4*i - 3, r) = this%nwf(4*i - 3, r) + h2t_z*this%wf(4*this%ns_z(i)%s(j) - 3, r)
            this%nwf(4*i - 2, r) = this%nwf(4*i - 2, r) + h2t_z*this%wf(4*this%ns_z(i)%s(j) - 2, r)
            this%nwf(4*i - 1, r) = this%nwf(4*i - 1, r) + h2t_z*this%wf(4*this%ns_z(i)%s(j) - 1, r)
            this%nwf(4*i, r) = this%nwf(4*i, r) + h2t_z*this%wf(4*this%ns_z(i)%s(j), r)
          ENDDO  ! j
        ENDDO ! i
        ! interaction of sites across holes
       ! DO k = 1, SIZE(hop%pl%value(:))
        DO k = 1, this%pl_h%length()
          !i = hop%pl%value(k)%i
          !f = hop%pl%value(k)%f
          i = this%pl_h%value(k)%i
          f = this%pl_h%value(k)%f
          this%nwf(4*f - 3, r) = this%nwf(4*f - 3, r) - h2jd*0.5d0*(Sx(i) - ui*Sy(i))*this%wf(4*f - 2, r)
          this%nwf(4*f - 2, r) = this%nwf(4*f - 2, r) - h2jd*0.5d0*(Sx(i) + ui*Sy(i))*this%wf(4*f - 3, r)
          this%nwf(4*i - 1, r) = this%nwf(4*i - 1, r) - h2jd*0.5d0*(Sx(f) - ui*Sy(f))*this%wf(4*i, r)
          this%nwf(4*i, r) = this%nwf(4*i, r) - h2jd*0.5d0*(Sx(f) + ui*Sy(f))*this%wf(4*i - 1, r)
        ENDDO ! k
        ! rashba
        ! DO i = 1, SIZE(p,1)
        !    DO j = 1, 4
        !       pi = p(i,j)%i
        !       pf = p(i,j)%f
        !       this%nwf(2*pf-1,r) = this%nwf(2*pf-1,r) - h2lm * (- this%wf(2*pi,r)  + ui*this%wf(2*pi,r))
        !       this%nwf(2*pf,r) = this%nwf(2*pf,r) - h2lm * (- this%wf(2*pi-1,r)  + ui*this%wf(2*pi-1,r))
        !       this%nwf(2*pi-1,r) = this%nwf(2*pi-1,r) - h2lm * (- this%wf(2*pf,r) + ui*this%wf(2*pf,r))
        !       this%nwf(2*pi,r) = this%nwf(2*pi,r) - h2lm * (- this%wf(2*pf-1,r)  + ui*this%wf(2*pf-1,r))
        !    ENDDO
        ! ENDDO
      ENDDO ! r
      this%owf = this%wf
      this%wf = c*this%nwf
      DO l = 1, this%nh
        this%wf(4*this%hole(l, 1) - 3, :) = 0.d0
        this%wf(4*this%hole(l, 1) - 2, :) = 0.d0
        this%wf(4*this%hole(l, 1) - 1, :) = 0.d0
        this%wf(4*this%hole(l, 1)    , :) = 0.d0
      ENDDO
      CALL this%gs()
      CALL this%set_eg()
      if(print_flg) then
      WRITE (*, *) 'sum energy', SUM(this%eg(1:4*this%n))
      end if
      if(print_flg) then
        WRITE (*, *) 'density', SUM(ABS(density(:))), 'Sx', SUM(ABS(Sx(:)))
        WRITE (*, *) 'Sy', SUM(ABS(Sy(:))), 'Sz', SUM(ABS(Sz(:)))
      end if
      IF (iter > 50) THEN
        flag = this%is_converged(oSx, oSy, oSz, oeg)
        if(flag) EXIT
      ENDIF
      ENDDO!iter
      end if



    ! OPEN(10,file ='dmu_ne_'//file//'.dat')
    !  OPEN(20,file ='display_'//file//'.dat')
    !  OPEN(30,file ='del_'//file//'.dat')
    !  OPEN(40,file ='deg_'//file//'.dat')
    !  OPEN(50,file ='del_all_'//file//'.dat')
    !  OPEN(60,file ='spin_'//file//'.dat')
    !  write(*,*)'file number '//file
    !  write(20,'(4(a10,i7,/),a10,f12.8)')'this%nx = ',this%nx,'this%ny = ',this%ny,'this%n  = ',this%n,'this%nh = ',this%nh,'this%ne = ',this%ne
    !  write(30,*)"p '-' w l"
    !  write(40,*)"p '-' w l"
    !  write(50,*)"p '-' "
    !  write(60,*)"sp '-'  w vec"
    !  do i =0,400
    !    mu_h = real(i)/100
    !    call this%calc_surface(mu_h,ne_h)
    !  end do
    !  close(10)
    !  close(20)
    !  close(30)
    !  close(40)
    !  close(50)
    !  close(60)




    !subroutine car_parinello_calc_surface_general(this, eg, wf, mu, delta, ne_u, ne_d, ratio_delta, ratio_ne)
    test_wf = 0d0
    do i = 1, 4*this%n
      test_wf(i,i) = 1d0
    end do
    !delta = 0d0
    ne_u = 0.5d0
    ne_d = 0.5d0
    delta  = 0.2d0
    delta(size(delta)/2+1:size(delta))  = -0.2d0
    !g__debug_out = .true.
    !read(41) test_wf
    !read(42) test_eg
    !read(43) delta
    !g__debug_out = .true.
    
    !mu_l = -4d0
    !mu_h = 0d0
    !open(101, file = "mu_ne.dat")
    !do i = 1, 4000
    !  mu_m = mu_l + (mu_h - mu_l)/4000.d0*i
    !  ne_u = 0.5d0
    !  ne_d = 0.5d0
    !  delta  = 0.2d0
    !  delta(size(delta)/2+1:size(delta))  = -0.2d0
    !  call this%calc_surface_general(test_eg, test_wf, mu_m, delta, ne_u, ne_d, 0d0, 0d0)
    !  call this%calc_surface_general(test_eg, test_wf, mu_m, delta, ne_u, ne_d, 0.01d0, 0d0)
    !  call this%calc_surface_general(test_eg, test_wf, mu_m, delta, ne_u, ne_d, 0d0, 0.01d0)
    !  call this%calc_surface_general(test_eg, test_wf, mu_m, delta, ne_u, ne_d, 1d0, 1d0)
    !  write(101, *) mu_m, sum(ne_u+ne_d)
    !end do
    !stop

    call this%calc_mu_surface(55d0, test_eg, test_wf, mu_m, delta, ne_u, ne_d)
    stop

    mu_m = -0.2d0
    !call diagonalize_surface_hamiltonian(this, test_eg, test_wf, mu_m, delta, ne_u, ne_d)
    !stop
    
    call this%calc_surface_general(test_eg, test_wf, mu_m, delta, ne_u, ne_d, 0d0, 0d0)
    call this%calc_surface_general(test_eg, test_wf, mu_m, delta, ne_u, ne_d, 0.01d0, 0d0)
    call this%calc_surface_general(test_eg, test_wf, mu_m, delta, ne_u, ne_d, 0d0, 0.01d0)
    call this%calc_surface_general(test_eg, test_wf, mu_m, delta, ne_u, ne_d, 1d0, 1d0)
    write(51) test_wf
    write(52) test_eg
    write(53) delta
    stop "for test, in 424"
    stop
    !g__debug_out = .true.
    !print *, "start cp"
    !call this%calc_surface(-0.2d0,ne_h)
    !!call this%calc_surface(-1d0,ne_h)




      file = '46'
      OPEN(10,file ='dmu_ne_'//file//'.dat')
      OPEN(20,file ='display_'//file//'.dat')
      OPEN(40,file ='deg_'//file//'.dat')
      OPEN(60,file ='spin_'//file//'.dat')
      write(*,*)'file number '//file
    !  write(20,'(4(a10,i7,/),a10,f12.8)')'this%nx = ',this%nx,'this%ny = ',this%ny,'this%n  = ',this%n,'this%nh = ',this%nh,'this%ne = ',this%ne
      write(20,*)'this%nx = ',this%nx,'this%ny = ',this%ny,'this%n  = ',this%n,'this%nh = ',this%nh,'this%ne = ',this%ne
    
    

      write(40,*)"p '-' w l"
      write(60,*)"sp '-'  w vec"
    !  do i =0,400
    !    mu_h = real(i)/100
        mu_h = 0
        call this%calc_bulk(mu_h,ne_h)
      !  call this%calc_surface(mu_h,ne_h)
    !  end do
      close(10)
      close(20)
      close(40)
      close(60)
      stop


      mu_h=  4.d0 
      mu_l= - 4.d0
      call this%calc_bulk(mu_h,ne_h)
      call this%calc_bulk(mu_l,ne_l)
      if( (ne_h - sum(int(this%ne)))*(ne_l - sum(int(this%ne))) > 0.d0  )stop 'car_parinello 373'
      do i =1,1000
        mu_m = (mu_h + mu_l)*5.d-1
        call this%calc_bulk(mu_m ,ne_m)
        write(20,'(2(3(a8,f15.8),/))') 'ne_l',ne_l,'  ne_m',ne_m,'  ne_h',ne_h,'mu_l',mu_l,'  mu_m',mu_m,'  mu_h',mu_h
        if(abs(mu_h - mu_l) < 1.0d-3)then
          exit
        else if((ne_h - sum(int(this%ne)))*(ne_m - sum(int(this%ne))) > 0 )then
          mu_h = mu_m
          ne_h = ne_m
        else 
          mu_l = mu_m
          ne_l = ne_m
        end if
      end do !i




    mu_h=  4.d0 
    mu_l= - 4.d0
    call this%calc_surface(mu_h,ne_h)
    call this%calc_surface(mu_l,ne_l)
    write(20,*)'ne_l = ',ne_l,'ne_h = ',ne_h
    if( (ne_h - sum(int(this%ne)))*(ne_l - sum(int(this%ne))) > 0.d0  )stop 'car_parinello 373'
      do i =1,1000
        mu_m = (mu_h + mu_l)*5.d-1
        call this%calc_surface(mu_m ,ne_m)
        write(20,'(3(a8,f15.8,/))') 'ne_l',ne_l,'  ne_m',ne_m,'  ne_h',ne_h,'mu_l',mu_l,'  mu_m',mu_m,'  mu_h',mu_h
        if(abs(mu_h - mu_l) < 1.0d-3)then
          exit
        else if((ne_h - sum(int(this%ne)))*(ne_m - sum(int(this%ne))) > 0 )then
          mu_h = mu_m
          ne_h = ne_m
        else 
          mu_l = mu_m
          ne_l = ne_m
        end if
      end do !i
      stop
      call this%calc_surface(mu_m ,ne_m,fe_e)
      write(20,*)'result mu = ',mu_m  
      stop
    this%xi = this%calc_xi()
  !call draw_phase_3D('test2_xi.gp', 44, this%pshape, this%xi, this%hole(:, 1))
  !print *, "saved to [test2_xi.gp]"
  !stop
    IF (flag) THEN
      WRITE (*, *) 'CP_caluculation has been converged!'
    !  WRITE (*, *) 'density', SUM(ABS(density(:))), 'Sx', SUM(ABS(Sx(:)))
    !  WRITE (*, *) 'Sy', SUM(ABS(Sy(:))), 'Sz', SUM(ABS(Sz(:)))
    !  WRITE (*, *) 'sum energy', SUM(this%eg(1:sum(int(this%ne))))
    !  WRITE (*, *) 'h2jd', h2jd
    ELSE
      WRITE (6, '(a30)') '-------- ! warning ! ---------'
      WRITE (*, *) 'CP_caluculation is not converged.'
    ENDIF
    !  write(*,*)this%eg
    !  allocate(t_ne(this%n))
  END SUBROUTINE car_parrinello_calc

  subroutine calc_surface(this,mu,ne,fe_e) 
    implicit none
    CLASS(car_parrinello), INTENT(in)   :: this
    INTEGER               :: k, l,r ,i ,j ,a ,iter, time_s,time_e,CountPerSec, CountMax, min,iter_one, iter_two
    real(8),intent(in)    :: mu
    real(8),intent(out)   :: ne
    REAL(8),optional      :: fe_e
    REAL(8)               :: c , ocoef, h2, h2t2, h2t, h2t_z, h2u, h2jd, h2lm, kbT, norm, norm2, tol,sec,ratio
    REAL(8), DIMENSION(this%n) ::dos_u,dos_d, ne_u, ne_d, old_u, old_d
    REAL(8), DIMENSION(4*this%n) :: tmp_eg, eg , oeg
    REAL(8)               :: ratio_delta
    complex(8)            ::  delta(size(this%pl_1st%value)), odelta(size(this%pl_1st%value))
    complex(8),DIMENSION(4*this%n,4*this%n)   :: nwf, owf, temp, wf
    !complex(8),DIMENSION(2*this%n,2*this%n)   :: matQ
    LOGICAL  temp2(4*this%n) !,dflag ,l_n
    !LOGICAL  l_n
    type(path)            :: p
    iter_one = 0
    iter_two = 0
    tol = 1.0d-5
    !kbT =1.0d-2
    kbT =0.03d0
    c = 1.d0/(1 + this%h*this%eta*0.5d0)
    ocoef = this%h*this%eta*0.5d0 - 1.d0
    h2 = this%h**2
    h2t   = h2*this%t(1)
    h2t2  = h2*this%t(2)
    h2t_z = h2*this%t(3)
    h2u = h2*this%u
    h2jd = this%jd*h2
    h2lm = this%lm*h2
    eg = 0d0
    oeg = 0d0
    owf = 0.0d0
    wf = 0.d0
    ratio = 1.d-3
    ratio_delta = 1d-2
    !dflag = .false.
    !l_n = .false.
    call system_clock(time_s, CountPerSec, CountMax)
    do i=1,4*this%n
      wf(i,i)=1.d0
      owf(i,i)=1.d0
    end do!i      
    !  ne_u = 0.5d0
    !  ne_d = 0.5d0
    !old_u = ne_u
    old_u = 0.d0
    old_d = 0.d0
    !old_d = ne_d
    delta  = 0.5d0
    delta(size(delta)/2+1:size(delta))  = -0.5d0
    !delta  = 0d0
    !print *,delta
    !stop
    read(41) wf
    read(42) eg
    read(43) delta
    owf = wf
    oeg = eg
    odelta = delta
    DO iter = 1,60000
        WRITE (*, '(a6,i7,a2,i7,a2)') '-iter', iter, ' /', this%itermax, ' -'
      !if(l_n)then
      ne_u = 0.d0
      ne_d = 0.d0
      do i = 1, this%n
        do a = 2*this%n+1, 4*this%n
        !do a = 1, 4*this%n
          !if(eg(a)  > 0.d0)then
            ne_u(i) = ne_u(i) + fermi( eg(a), kbT)*conjg(wf(4*i-3, a))*wf(4*i-3, a) + fermi(-eg(a), kbT)*conjg(wf(4*i-1, a))*wf(4*i-1, a)
            ne_d(i) = ne_d(i) + fermi( eg(a), kbT)*conjg(wf(4*i-2, a))*wf(4*i-2, a) + fermi(-eg(a), kbT)*conjg(wf(4*i  , a))*wf(4*i  , a) 
          !end if
        end do!a
      end do!i
      !ne_u = ratio*ne_u + (1d0-ratio)*old_u
      !ne_d = ratio*ne_d + (1d0-ratio)*old_d
      !old_u = ne_u
      !old_d = ne_d
      !end if
        !  delta  = 0.5d0
      !if(dflag)then
      odelta = delta
      delta  = 0.d0
      do a = 2*this%n+1, 4*this%n
        !if(eg(a) >= 0.d0)then
        do k = 1, this%pl_1st%length()
          i = this%pl_1st%value(k)%i
          j = this%pl_1st%value(k)%f
          delta(k) = delta(k) + (this%jd*(wf(4*i-3, a)*conjg(wf(4*j  , a)) + wf(4*j-2, a)*conjg(wf(4*i-1, a)) &
          &                 + wf(4*i-2, a)*conjg(wf(4*j-1, a)) + wf(4*j-3, a)*conjg(wf(4*i  , a))))*tanh(eg(a)/kbT*0.5d0)
        end do!k
        !end if
      end do!a
      delta = ratio_delta*delta + (1-ratio_delta)*odelta
      !if(iter == 1000) then
      !  do k = 1, this%pl_1st%length()
      !    print "(i,4f)", k, delta(k), odelta(k)
      !  end do
      !  stop
      !end if
      !end if
      nwf = 0.d0
      DO r = 1, 4*this%n
        DO k = 1, this%n - this%nh
          i = this%ets(k)
          nwf(4*i - 3, r) = (2.d0 + h2*eg(r) )*wf(4*i - 3, r) + ocoef*owf(4*i - 3, r) + h2*mu*wf(4*i - 3, r)
          nwf(4*i - 2, r) = (2.d0 + h2*eg(r) )*wf(4*i - 2, r) + ocoef*owf(4*i - 2 ,r) + h2*mu*wf(4*i - 2, r)
          nwf(4*i - 1, r) = (2.d0 + h2*eg(r) )*wf(4*i - 1, r) + ocoef*owf(4*i - 1, r) - h2*mu*wf(4*i - 1, r)
          nwf(4*i    , r) = (2.d0 + h2*eg(r) )*wf(4*i    , r) + ocoef*owf(4*i    , r) - h2*mu*wf(4*i    , r)
          !nwf(4*k - 3, r) = nwf(4*k - 3, r)  -       delta(1) * h2*wf(4*k    , r)
          !nwf(4*k - 2, r) = nwf(4*k - 2, r)  -       delta(1) * h2*wf(4*k - 1, r)
          !nwf(4*k - 1, r) = nwf(4*k - 1, r)  - conjg(delta(1))* h2*wf(4*k - 2, r)
          !nwf(4*k    , r) = nwf(4*k    , r)  - conjg(delta(1))* h2*wf(4*k - 3, r) 
          DO j = 1, SIZE(this%ns_1st(i)%s)       
            nwf(4*i - 3, r) = nwf(4*i - 3, r) + h2t*wf(4*this%ns_1st(i)%s(j) - 3, r)
            nwf(4*i - 2, r) = nwf(4*i - 2, r) + h2t*wf(4*this%ns_1st(i)%s(j) - 2, r)
            nwf(4*i - 1, r) = nwf(4*i - 1, r) - h2t*wf(4*this%ns_1st(i)%s(j) - 1, r)
            nwf(4*i    , r) = nwf(4*i    , r) - h2t*wf(4*this%ns_1st(i)%s(j)    , r)
            !nwf(4*i - 3, r) = nwf(4*i - 3, r) + h2t*wf(4*this%ns_1st(i)%s(j) - 3, r)*(1.d0 - ne_d(i))*(1.d0 - ne_d(this%ns_1st(i)%s(j)))
            !nwf(4*i - 2, r) = nwf(4*i - 2, r) + h2t*wf(4*this%ns_1st(i)%s(j) - 2, r)*(1.d0 - ne_u(i))*(1.d0 - ne_u(this%ns_1st(i)%s(j)))
            !nwf(4*i - 1, r) = nwf(4*i - 1, r) - h2t*wf(4*this%ns_1st(i)%s(j) - 1, r)*(1.d0 - ne_d(i))*(1.d0 - ne_d(this%ns_1st(i)%s(j)))
            !nwf(4*i    , r) = nwf(4*i    , r) - h2t*wf(4*this%ns_1st(i)%s(j)    , r)*(1.d0 - ne_u(i))*(1.d0 - ne_u(this%ns_1st(i)%s(j)))
          ENDDO ! j
          DO j = 1, SIZE(this%ns_2nd(i)%s)
            nwf(4*i - 3, r) = nwf(4*i - 3, r) + h2t2*wf(4*this%ns_2nd(i)%s(j) - 3, r) 
            nwf(4*i - 2, r) = nwf(4*i - 2, r) + h2t2*wf(4*this%ns_2nd(i)%s(j) - 2, r)
            nwf(4*i - 1, r) = nwf(4*i - 1, r) - h2t2*wf(4*this%ns_2nd(i)%s(j) - 1, r)
            nwf(4*i    , r) = nwf(4*i    , r) - h2t2*wf(4*this%ns_2nd(i)%s(j)    , r)
            !nwf(4*i - 3, r) = nwf(4*i - 3, r) + h2t2*wf(4*this%ns_2nd(i)%s(j) - 3, r)*(1.d0 - ne_d(i))*(1.d0 - ne_d(this%ns_2nd(i)%s(j))) 
            !nwf(4*i - 2, r) = nwf(4*i - 2, r) + h2t2*wf(4*this%ns_2nd(i)%s(j) - 2, r)*(1.d0 - ne_u(i))*(1.d0 - ne_u(this%ns_2nd(i)%s(j)))
            !nwf(4*i - 1, r) = nwf(4*i - 1, r) - h2t2*wf(4*this%ns_2nd(i)%s(j) - 1, r)*(1.d0 - ne_d(i))*(1.d0 - ne_d(this%ns_2nd(i)%s(j)))
            !nwf(4*i    , r) = nwf(4*i    , r) - h2t2*wf(4*this%ns_2nd(i)%s(j)    , r)*(1.d0 - ne_u(i))*(1.d0 - ne_u(this%ns_2nd(i)%s(j)))
          ENDDO ! j    
        end do!k
        !if(dflag ) then
        do i = 1, this%pl_1st%length()
          k = this%pl_1st%value(i)%i
          j = this%pl_1st%value(i)%f
          nwf(4*k - 3, r) = nwf(4*k - 3, r)  -       delta(i) * h2*wf(4*j    , r)
          nwf(4*k - 2, r) = nwf(4*k - 2, r)  -       delta(i) * h2*wf(4*j - 1, r)
          nwf(4*k - 1, r) = nwf(4*k - 1, r)  - conjg(delta(i))* h2*wf(4*j - 2, r)
          nwf(4*k    , r) = nwf(4*k    , r)  - conjg(delta(i))* h2*wf(4*j - 3, r) 
          nwf(4*j - 3, r) = nwf(4*j - 3, r)  -       delta(i) * h2*wf(4*k    , r)
          nwf(4*j - 2, r) = nwf(4*j - 2, r)  -       delta(i) * h2*wf(4*k - 1, r)
          nwf(4*j - 1, r) = nwf(4*j - 1, r)  - conjg(delta(i))* h2*wf(4*k - 2, r)
          nwf(4*j    , r) = nwf(4*j    , r)  - conjg(delta(i))* h2*wf(4*k - 3, r) 
        ENDDO ! i
        !end if
      end do!r

      owf = wf
      wf = c*nwf

    !---  temp = c*nwf
    !---  DO k = 1, 2*this%n
    !---  !DO k = 1, 4*this%n
    !---    DO l = 1, k - 1
    !---      temp(:, k) = temp(:, k) - DOT_PRODUCT(temp(:, l),c*nwf(:, k))*temp(:, l) ! Classic Gram-Schmidt
    !---      !temp(:, k) = temp(:, k) - DOT_PRODUCT(temp(:, l),temp(:, k))*temp(:, l)  ! Modfied Gram- Schmidt
    !---    ENDDO!l
    !---    !norm = SQRT(DOT_PRODUCT(temp(:,k), temp(:,k)))
    !---    norm2 = DOT_PRODUCT(temp(:,k), temp(:,k))
    !---    !IF (norm == 0.d0) THEN
    !---    IF (norm2 == 0.d0) THEN
    !---      temp(:,k) = 0.d0
    !---    ELSE
    !---      !temp(:,k) = temp(:,k)/norm
    !---      temp(:,k) = temp(:,k)*norm2**-0.5d0
    !---    ENDIF
    !---  ENDDO!k
    !---  wf = temp
      
      !call gram_schmidt(wf)
      !call gram_schmidt(wf, 2*this%n)
      
      !DO k = 1, 2*this%n
      !  DO l = 1, k - 1
      !    temp(:, k) = temp(:, k) - DOT_PRODUCT(temp(:, l),c*nwf(:, k))*temp(:, l) ! Classic Gram-Schmidt
      !    !temp(:, k) = temp(:, k) - DOT_PRODUCT(temp(:, l),temp(:, k))*temp(:, l)  ! Modfied Gram- Schmidt
      !  ENDDO!l
      !  norm2 = DOT_PRODUCT(temp(:,k), temp(:,k))
      !  IF (norm2 == 0.d0) THEN
      !    temp(:,k) = 0.d0
      !  ELSE
      !    temp(:,k) = temp(:,k)*norm2**-0.5d0
      !  ENDIF
      !ENDDO!k
      !DO k = 2*this%n + 1, 4*this%n
      !  DO l = 2*this%n+1, k - 1
      !    temp(:, k) = temp(:, k) - DOT_PRODUCT(temp(:, l),c*nwf(:, k))*temp(:, l) ! Classic Gram-Schmidt
      !    !temp(:, k) = temp(:, k) - DOT_PRODUCT(temp(:, l),temp(:, k))*temp(:, l)  ! Modfied Gram- Schmidt
      !  ENDDO!l
      !  norm2 = DOT_PRODUCT(temp(:,k), temp(:,k))
      !  IF (norm2 == 0.d0) THEN
      !    temp(:,k) = 0.d0
      !  ELSE
      !    temp(:,k) = temp(:,k)*norm2**-0.5d0
      !  ENDIF
      !ENDDO!k
      !!DO concurrent(k = 1:2*this%n)
      !!  DO concurrent(l = 1:2*this%n)
      !!$OMP PARALLEL DO
      !DO k = 1, 2*this%n
      !  DO l = 1, 2*this%n
      !   matQ(k,l) = DOT_PRODUCT(temp(:,k), temp(:,l+2*this%n))
      !  end do
      !end do
      !!$OMP END PARALLEL DO
      !!DO concurrent(k = 2*this%n + 1:4*this%n)
      !!$OMP PARALLEL DO
      !DO k = 2*this%n + 1, 4*this%n
      !  DO l = 1, 2*this%n
      !    temp(:, k) = temp(:, k) - matQ(l,k-2*this%n)*temp(:, l) ! Classic Gram-Schmidt
      !    !temp(:, k) = temp(:, k) - DOT_PRODUCT(temp(:, l),temp(:, k))*temp(:, l)  ! Modfied Gram- Schmidt
      !  ENDDO!l
      !  norm2 = DOT_PRODUCT(temp(:,k), temp(:,k))
      !  IF (norm2 == 0.d0) THEN
      !    temp(:,k) = 0.d0
      !  ELSE
      !    temp(:,k) = temp(:,k)*norm2**-0.5d0
      !  ENDIF
      !ENDDO!k
      !!$OMP END PARALLEL DO

      !!DO k = 1, 2*this%n
      !DO k = 1, 4*this%n
      !  DO l = 1, k - 1
      !    temp(:, k) = temp(:, k) - DOT_PRODUCT(temp(:, l),c*nwf(:, k))*temp(:, l) ! Classic Gram-Schmidt
      !    !temp(:, k) = temp(:, k) - DOT_PRODUCT(temp(:, l),temp(:, k))*temp(:, l)  ! Modfied Gram- Schmidt
      !  ENDDO!l
      !  !norm = SQRT(DOT_PRODUCT(temp(:,k), temp(:,k)))
      !  norm2 = DOT_PRODUCT(temp(:,k), temp(:,k))
      !  !IF (norm == 0.d0) THEN
      !  IF (norm2 == 0.d0) THEN
      !    temp(:,k) = 0.d0
      !  ELSE
      !    !temp(:,k) = temp(:,k)/norm
      !    temp(:,k) = temp(:,k)*norm2**-0.5d0
      !  ENDIF
      !ENDDO!k
      !wf = temp
      
      !if(not(check_sorted_eg(eg))) then
      !  print *,"sort T"
      !  call sort_wf(eg,wf)
      !else
      !  print *,"sort F"
      !end if
      
    !---  eg = 0.d0
    !---  tmp_eg = 0d0
    !---  do i = 1, this%pl_1st%length()
    !---    k = this%pl_1st%value(i)%i
    !---    j = this%pl_1st%value(i)%f
    !---    do a = 1, 2*this%n
    !---    !do a = 1, 4*this%n
    !---        tmp_eg(a) = tmp_eg(a) + this%t(1)*( (conjg(wf(4*k-3,a))*wf(4*j-3,a) + conjg(wf(4*j-3,a))*wf(4*k-3,a)) &                           
    !---        &                                 + (conjg(wf(4*k-2,a))*wf(4*j-2,a) + conjg(wf(4*j-2,a))*wf(4*k-2,a)) &
    !---        &                                 - (conjg(wf(4*k-1,a))*wf(4*j-1,a) + conjg(wf(4*j-1,a))*wf(4*k-1,a)) &
    !---        &                                 - (conjg(wf(4*k  ,a))*wf(4*j  ,a) + conjg(wf(4*j  ,a))*wf(4*k  ,a))  )
    !---        !tmp_eg(a) = tmp_eg(a) + this%t(1)*( (conjg(wf(4*k-3,a))*wf(4*j-3,a) + conjg(wf(4*j-3,a))*wf(4*k-3,a))*(1.d0 - ne_d(k))*(1.d0 - ne_d(j)) &                             
    !---        !&                                        + (conjg(wf(4*k-2,a))*wf(4*j-2,a) + conjg(wf(4*j-2,a))*wf(4*k-2,a))*(1.d0 - ne_u(k))*(1.d0 - ne_u(j)) &
    !---        !&                                        - (conjg(wf(4*k-1,a))*wf(4*j-1,a) + conjg(wf(4*j-1,a))*wf(4*k-1,a))*(1.d0 - ne_d(k))*(1.d0 - ne_d(j))&
    !---        !&                                        - (conjg(wf(4*k  ,a))*wf(4*j  ,a) + conjg(wf(4*j  ,a))*wf(4*k  ,a))*(1.d0 - ne_u(k))*(1.d0 - ne_u(j)) )
    !---      !if(dflag)  tmp_eg(a) = tmp_eg(a) + delta(i)*( conjg(wf(4*k-3,a))*wf(4*j  ,a) + conjg(wf(4*j-2,a))*wf(4*k-1,a) &
    !---      ! tmp_eg(a) = tmp_eg(a) + delta(i)*( conjg(wf(4*k-3,a))*wf(4*j  ,a) + conjg(wf(4*j-2,a))*wf(4*k-1,a) &
    !---      !&                                           + conjg(wf(4*j-3,a))*wf(4*k  ,a) + conjg(wf(4*k-2,a))*wf(4*j-1,a))&
    !---      !&                         - conjg(delta(i))*( conjg(wf(4*k-1,a))*wf(4*j-2,a) + conjg(wf(4*j  ,a))*wf(4*k-3,a) &
    !---      !&                                          +  conjg(wf(4*j-1,a))*wf(4*k-2,a) + conjg(wf(4*k  ,a))*wf(4*j-3,a))
    !---       tmp_eg(a) = tmp_eg(a) - delta(i)*( conjg(wf(4*k-3,a))*wf(4*j  ,a) + conjg(wf(4*j-3,a))*wf(4*k  ,a) &
    !---      &                                 + conjg(wf(4*k-2,a))*wf(4*j-1,a) + conjg(wf(4*j-2,a))*wf(4*k-1,a))&
    !---      &               - conjg(delta(i))*( conjg(wf(4*k-1,a))*wf(4*j-2,a) + conjg(wf(4*j-1,a))*wf(4*k-2,a) &
    !---      &                                 + conjg(wf(4*k  ,a))*wf(4*j-3,a) + conjg(wf(4*j  ,a))*wf(4*k-3,a))
    !---    end do!a
    !---  end do!i
    !---  eg(1:2*this%n) = eg(1:2*this%n) - tmp_eg(1:2*this%n) 
    !---  !eg(1:4*this%n) = eg(1:4*this%n) - tmp_eg(1:4*this%n) 
    !---  tmp_eg = 0d0
    !---  do i = 1, this%pl_2nd%length()
    !---    k = this%pl_2nd%value(i)%i
    !---    j = this%pl_2nd%value(i)%f
    !---    do a = 1, 2*this%n
    !---    !do a = 1, 4*this%n
    !---      tmp_eg(a) = tmp_eg(a) + (conjg(wf(4*k-3,a))*wf(4*j-3,a) + conjg(wf(4*j-3,a))*wf(4*k-3,a)) &
    !---      &                     + (conjg(wf(4*k-2,a))*wf(4*j-2,a) + conjg(wf(4*j-2,a))*wf(4*k-2,a)) &
    !---      &                     - (conjg(wf(4*k-1,a))*wf(4*j-1,a) + conjg(wf(4*j-1,a))*wf(4*k-1,a)) &
    !---      &                     - (conjg(wf(4*k  ,a))*wf(4*j  ,a) + conjg(wf(4*j  ,a))*wf(4*k  ,a))  
    !---      !tmp_eg(a) = tmp_eg(a) + (conjg(wf(4*k-3,a))*wf(4*j-3,a) + conjg(wf(4*j-3,a))*wf(4*k-3,a))*(1.d0 - ne_d(k))*(1.d0 - ne_d(j)) &
    !---      !&                                  + (conjg(wf(4*k-2,a))*wf(4*j-2,a) + conjg(wf(4*j-2,a))*wf(4*k-2,a))*(1.d0 - ne_u(k))*(1.d0 - ne_u(j)) &
    !---      !&                                  - (conjg(wf(4*k-1,a))*wf(4*j-1,a) + conjg(wf(4*j-1,a))*wf(4*k-1,a))*(1.d0 - ne_d(k))*(1.d0 - ne_d(j)) &
    !---      !&                                  - (conjg(wf(4*k  ,a))*wf(4*j  ,a) + conjg(wf(4*j  ,a))*wf(4*k  ,a))*(1.d0 - ne_u(k))*(1.d0 - ne_u(j))  
    !---    end do!a
    !---  end do!i
    !---  eg(1:2*this%n) = eg(1:2*this%n) - tmp_eg(1:2*this%n) * this%t(2) 
    !---  !eg(1:4*this%n) = eg(1:4*this%n) - tmp_eg(1:4*this%n) * this%t(2) 
    !---  tmp_eg = 0d0
    !---  !do a = 1, 2*this%n 
    !---  do a = 1, 4*this%n 
    !---    do k = 1, this%n -  this%nh
    !---      j = this%ets(k)
    !---      tmp_eg(a) = tmp_eg(a) + conjg(wf(4*j-3,a))*wf(4*j-3,a) + conjg(wf(4*j-2,a))*wf(4*j-2,a) &
    !---      &                     - conjg(wf(4*j-1,a))*wf(4*j-1,a) - conjg(wf(4*j  ,a))*wf(4*j  ,a) 
    !---    end do!j

    !---  end do!a
    !---  eg(1:2*this%n) = eg(1:2*this%n) - tmp_eg(1:2*this%n) * mu
    !---  !eg(1:4*this%n) = eg(1:4*this%n) - tmp_eg(1:4*this%n) * mu
      
      if (.false.) then !!!! othogonization for all 
        
        call gram_schmidt(wf)
        call calc_surface_eg(this, eg, wf, delta, mu, ne_u, ne_d)
      
      else  !!!! othogonization for half
      
        call gram_schmidt(wf, 2*this%n)
        call calc_surface_eg(this, eg, wf, delta, mu, ne_u, ne_d)

        ! set E_{-n}, phi_{-n}
        do a = 1, 2*this%n 
          eg(4*this%n - a + 1) = - eg(a)
          do i = 1, this%n 
            wf(4*i - 3, 4*this%n - a + 1) = - conjg(wf(4*i - 1, a))
            wf(4*i - 2, 4*this%n - a + 1) =   conjg(wf(4*i    , a))
            wf(4*i - 1, 4*this%n - a + 1) = - conjg(wf(4*i - 3, a))
            wf(4*i    , 4*this%n - a + 1) =   conjg(wf(4*i - 2, a))
          end do
        end do
     
      end if
      
      !if(iter == 100) then
      !call calc_surface_eg(this, oeg, wf, delta, mu, ne_u, ne_d)
      !do i=1, this%n*4
      !  print "(i,3f)", i, eg(i), oeg(i), eg(i)- oeg(i)
      !end do
      !stop
      !end if
  
     !if(iter <= 1000 .and. mod(iter, 100) == 0 .and. not(check_sorted_eg(eg))) then
     if(mod(iter, 100) == 0 .and. not(check_sorted_eg(eg))) then
     !if(not(check_sorted_eg(eg))) then
        print *, "sorted"
        !call sort_wf(eg,wf)
        call sort_wf(eg, wf, oeg, owf)
      end if

      WRITE (*, *) 'sum energy', SUM(eg(2*this%n+1:4*this%n))!, l_n
      WRITE (*, *) 'sum d energy', SUM(ABS((eg(:) - oeg(:))))
      WRITE (*, *) 'max abs(delta)', maxval(ABS((delta(:))))
      write(*,"(a,f10.5,a,f10.5,a,f8.5)") "ne; up : " ,sum(ne_u), "  down : ", sum(ne_d), "  max site : ", maxval(ne_u + ne_d)
      temp2(:) = .FALSE.
      do i = 1,4*this%n
        if(abs(oeg(i))  <  tol .or. ABS((eg(i) - oeg(i))/oeg(i)) < tol)  temp2(i) = .true.
        !if(abs(oeg(i)) >= 1.0d-8)write(60,'(2i,2f18.8,l,f18.8)')iter,i ,eg(i) ,oeg(i), temp2(i) ,ABS((eg(i) - oeg(i))/oeg(i))
        !if(abs(oeg(i)) <  1.0d-8)write(60,'(2i,2f18.8,l)')iter,i ,eg(i) ,oeg(i), temp2(i) 
      end do!i
      oeg(:) = eg(:)
      !if(dflag .and. (iter - iter_two) >= 50 .and. ALL(temp2)) EXIT
      !if((l_n .and. iter - iter_one) >= 50 .and. ALL(temp2)) EXIT
      if( iter >= 50 .and. ALL(temp2)) EXIT
      !if(not(l_n) .and. iter >= 50 .and. ALL(temp2) ) then
      !    l_n = .true.
      !    iter_one = iter
      !end if
     ! if(not(dflag) .and. l_n .and. (iter - iter_one)  >= 50 .and. ALL(temp2) )then
    !  if(not(dflag) .and.  (iter - iter_one)  >= 50 .and. ALL(temp2) )then
     !    dflag = .true.
     !    iter_two = iter
     ! end if

    end do !iter
    !print *,eg
    !do i = 1, 2*this%n
    !print *,eg(2*i-1),eg(2*i), eg(2*i-1) - eg(2*i)
    !end do
    
    !do i = 1, 2*this%n-1
    !  print *,eg(2*i-1),eg(2*(i+1)-1), eg(2*i-1) - eg(2*(i+1)-1)
    !end do
    if(not(check_sorted_eg(eg))) then
      call sort_wf(eg,wf)
    end if

  !  if(present(fe_e))then
  !     dos_u = 0.d0
  !     dos_u = 0.d0
  !  end if
       !   if(present(fe_e))then
       !   dos_u(i) = - ( dfermi( eg(a)-fe_e, kbT)*conjg(wf(4*i-3, a))*wf(4*i-3, a) + dfermi( eg(a)+fe_e, kbT)*conjg(wf(4*i-1, a))*wf(4*i-1, a)  )
       !   dos_d(i) = - ( dfermi( eg(a)-fe_e, kbT)*conjg(wf(4*i-2, a))*wf(4*i-2, a) + dfermi( eg(a)+fe_e, kbT)*conjg(wf(4*i  , a))*wf(4*i  , a)  )
       !   end if
    do i = 1, 4*this%n
    write(40,*)i,eg(i)
    end do
    ne   = 0.d0
    do i = 1, this%n
        do a = 1, 4*this%n
        if(eg(a)  > 0.d0)then
        ne = ne  + fermi( eg(a), kbT)*(conjg(wf(4*i-3, a))*wf(4*i-3, a) + conjg(wf(4*i-2, a))*wf(4*i-2, a))&
              &  + fermi(-eg(a), kbT)*(conjg(wf(4*i-1, a))*wf(4*i-1, a) + conjg(wf(4*i  , a))*wf(4*i  , a)) 
        end if
        end do!a
    end do!i
    do i = 1,this%pl_1st%length()
        write(30,*)i,real(delta(i))
        write(50,*)i,real(delta(i))
    end do!i
    write(30,*)'exit'
    write(30,*)"p '-' w l"
    write(40,*)'exit'
    write(40,*)"p '-' w l"
    write(50,*)''
    call system_clock(time_e)
    sec = real(time_e - time_s)/CountPerSec
    min = int(sec/60)
    write(20,'(a7,f13.5,a3,a7,3i7,2(a7,f16.11))') 'time = ',sec ,'(s) ',' iter = ',iter,iter_two,iter_one,'mu = ',mu,'ne = ',ne
    write(10,*)mu,ne,iter
    if(g__debug_out) then
      write(51) wf
      write(52) eg
      write(53) delta
      stop "for test, in 676"
    end if
  end subroutine calc_surface


  subroutine car_parinello_calc_surface_general(this, eg, wf, mu, delta, ne_u, ne_d, ratio_delta, ratio_ne)
    implicit none
    CLASS(car_parrinello), INTENT(in)   :: this
    REAL(8), DIMENSION(4*this%n), intent(out) :: eg
    complex(8),DIMENSION(4*this%n,4*this%n), intent(inout)   :: wf
    real(8),intent(in)    :: mu
    complex(8), dimension(this%pl_1st%length()), intent(inout) :: delta
    complex(8), dimension(this%pl_1st%length()) :: odelta
    real(8), dimension(this%n), intent(inout)   :: ne_u, ne_d
    real(8), dimension(this%n) :: old_ne_u, old_ne_d
    real(8), intent(in)   :: ratio_delta, ratio_ne
    !logical, intent(in), optional :: slfc_delta, slfc_ne
    REAL(8), DIMENSION(4*this%n) :: oeg
    REAL(8)               :: c , ocoef, h2, h2t2, h2t!, h2t_z, h2u, h2jd, h2lm!, norm, norm2, tol,sec,ratio
    REAL(8)               :: kbT
    !REAL(8), DIMENSION(this%n) ::dos_u,dos_d, ne_u, ne_d, old_u, old_d
    complex(8),DIMENSION(4*this%n,4*this%n)   :: nwf, owf
    !complex(8),DIMENSION(2*this%n,2*this%n)   :: matQ
    !LOGICAL  temp2(4*this%n) !,dflag ,l_n
    !LOGICAL  temp2(4*this%n) 
    integer               :: k, l, r ,i ,j ,a ,iter
    LOGICAL conv_eg, conv_ne_u, conv_ne_d, conv_delta_real, conv_delta_img
    type(path)            :: p
    !tol = 1.0d-5
    kbT =0.01d0
    c = 1.d0/(1 + this%h*this%eta*0.5d0)
    ocoef = this%h*this%eta*0.5d0 - 1.d0
    h2 = this%h**2
    h2t   = h2*this%t(1)
    h2t2  = h2*this%t(2)
    !h2t_z = h2*this%t(3)
    !h2u = h2*this%u
    !h2jd = this%jd*h2
    !h2lm = this%lm*h2
    
    call calc_surface_eg(this, eg, wf, delta, mu, ne_u, ne_d)
    oeg = eg
    owf = wf
    old_ne_u = ne_u
    old_ne_d = ne_d
    odelta = delta
    DO iter = 1, this%itermax
      WRITE (*, '(a6,i7,a2,i7,a2)') '-iter', iter, ' /', this%itermax, ' -'
      if(ratio_ne > 0d0) then
        old_ne_u = ne_u
        old_ne_d = ne_d
        call calc_surface_ne(this, wf, eg, ne_u, ne_d, kbT)
        ne_u = ratio_ne * ne_u + (1d0 - ratio_ne)*old_ne_u
        ne_d = ratio_ne * ne_d + (1d0 - ratio_ne)*old_ne_d
      end if
      
      if(ratio_delta > 0d0) then
        odelta = delta
        call calc_surface_delta(this, wf, eg, delta, kbT)
        delta = ratio_delta*delta + (1-ratio_delta)*odelta
      end if
      
      nwf = 0.d0
      DO r = 1, 4*this%n
        DO k = 1, this%n - this%nh
          i = this%ets(k)
          nwf(4*i - 3, r) = (2.d0 + h2*eg(r) )*wf(4*i - 3, r) + ocoef*owf(4*i - 3, r) + h2*mu*wf(4*i - 3, r)
          nwf(4*i - 2, r) = (2.d0 + h2*eg(r) )*wf(4*i - 2, r) + ocoef*owf(4*i - 2 ,r) + h2*mu*wf(4*i - 2, r)
          nwf(4*i - 1, r) = (2.d0 + h2*eg(r) )*wf(4*i - 1, r) + ocoef*owf(4*i - 1, r) - h2*mu*wf(4*i - 1, r)
          nwf(4*i    , r) = (2.d0 + h2*eg(r) )*wf(4*i    , r) + ocoef*owf(4*i    , r) - h2*mu*wf(4*i    , r)
          DO l = 1, SIZE(this%ns_1st(i)%s)
            j = this%ns_1st(i)%s(l)
            !nwf(4*i - 3, r) = nwf(4*i - 3, r) + h2t*wf(4*j - 3, r)
            !nwf(4*i - 2, r) = nwf(4*i - 2, r) + h2t*wf(4*j - 2, r)
            !nwf(4*i - 1, r) = nwf(4*i - 1, r) - h2t*wf(4*j - 1, r)
            !nwf(4*i    , r) = nwf(4*i    , r) - h2t*wf(4*j    , r)
            nwf(4*i - 3, r) = nwf(4*i - 3, r) + h2t*wf(4*j - 3, r)*(1.d0 - ne_d(i))*(1.d0 - ne_d(j))
            nwf(4*i - 2, r) = nwf(4*i - 2, r) + h2t*wf(4*j - 2, r)*(1.d0 - ne_u(i))*(1.d0 - ne_u(j))
            nwf(4*i - 1, r) = nwf(4*i - 1, r) - h2t*wf(4*j - 1, r)*(1.d0 - ne_d(i))*(1.d0 - ne_d(j))
            nwf(4*i    , r) = nwf(4*i    , r) - h2t*wf(4*j    , r)*(1.d0 - ne_u(i))*(1.d0 - ne_u(j))
          ENDDO ! j
          DO l = 1, SIZE(this%ns_2nd(i)%s)
            j = this%ns_2nd(i)%s(l)
            !nwf(4*i - 3, r) = nwf(4*i - 3, r) + h2t2*wf(4*j - 3, r)
            !nwf(4*i - 2, r) = nwf(4*i - 2, r) + h2t2*wf(4*j - 2, r)
            !nwf(4*i - 1, r) = nwf(4*i - 1, r) - h2t2*wf(4*j - 1, r)
            !nwf(4*i    , r) = nwf(4*i    , r) - h2t2*wf(4*j    , r)
            nwf(4*i - 3, r) = nwf(4*i - 3, r) + h2t2*wf(4*j - 3, r)*(1.d0 - ne_d(i))*(1.d0 - ne_d(j))
            nwf(4*i - 2, r) = nwf(4*i - 2, r) + h2t2*wf(4*j - 2, r)*(1.d0 - ne_u(i))*(1.d0 - ne_u(j))
            nwf(4*i - 1, r) = nwf(4*i - 1, r) - h2t2*wf(4*j - 1, r)*(1.d0 - ne_d(i))*(1.d0 - ne_d(j))
            nwf(4*i    , r) = nwf(4*i    , r) - h2t2*wf(4*j    , r)*(1.d0 - ne_u(i))*(1.d0 - ne_u(j))
          ENDDO ! j    
        end do!k
        !if(dflag ) then
        do i = 1, this%pl_1st%length()
          k = this%pl_1st%value(i)%i
          j = this%pl_1st%value(i)%f
          nwf(4*k - 3, r) = nwf(4*k - 3, r)  -       delta(i) * h2*wf(4*j    , r)
          nwf(4*k - 2, r) = nwf(4*k - 2, r)  -       delta(i) * h2*wf(4*j - 1, r)
          nwf(4*k - 1, r) = nwf(4*k - 1, r)  - conjg(delta(i))* h2*wf(4*j - 2, r)
          nwf(4*k    , r) = nwf(4*k    , r)  - conjg(delta(i))* h2*wf(4*j - 3, r) 
          nwf(4*j - 3, r) = nwf(4*j - 3, r)  -       delta(i) * h2*wf(4*k    , r)
          nwf(4*j - 2, r) = nwf(4*j - 2, r)  -       delta(i) * h2*wf(4*k - 1, r)
          nwf(4*j - 1, r) = nwf(4*j - 1, r)  - conjg(delta(i))* h2*wf(4*k - 2, r)
          nwf(4*j    , r) = nwf(4*j    , r)  - conjg(delta(i))* h2*wf(4*k - 3, r) 
        ENDDO ! i
        !end if
      end do!r

      owf = wf
      wf = c*nwf
      
      if (.false.) then !!!! othogonization for all 
        call gram_schmidt(wf)
        call calc_surface_eg(this, eg, wf, delta, mu, ne_u, ne_d)
      else  !!!! othogonization for half
        call gram_schmidt(wf, 2*this%n)
        call calc_surface_eg(this, eg, wf, delta, mu, ne_u, ne_d)
        do a = 1, 2*this%n      ! set E_{-n} & phi_{-n}
          eg(4*this%n - a + 1) = - eg(a)
          do i = 1, this%n 
            wf(4*i - 3, 4*this%n - a + 1) = - conjg(wf(4*i - 1, a))
            wf(4*i - 2, 4*this%n - a + 1) =   conjg(wf(4*i    , a))
            wf(4*i - 1, 4*this%n - a + 1) = - conjg(wf(4*i - 3, a))
            wf(4*i    , 4*this%n - a + 1) =   conjg(wf(4*i - 2, a))
          end do
        end do
      end if
      
     !if(iter <= 1000 .and. mod(iter, 100) == 0 .and. not(check_sorted_eg(eg))) then
     if(mod(iter, 100) == 0 .and. not(check_sorted_eg(eg))) then
     !if(not(check_sorted_eg(eg))) then
        print *, "sorted"
        call sort_wf(eg, wf, oeg, owf)
      end if

      WRITE (*, *) 'sum energy', SUM(eg(2*this%n+1:4*this%n))
      WRITE (*, *) 'sum d energy', SUM(ABS((eg(:) - oeg(:))))
      WRITE (*, *) 'max abs(delta)', maxval(ABS((delta(:))))
      write(*,"(a,f10.5,a,f10.5,a,f8.5)") "ne; up : " ,sum(ne_u), "  down : ", sum(ne_d), "  max site : ", maxval(ne_u + ne_d)
      
      conv_eg = calc_surface_converged(eg, oeg, this%tol)
      if(ratio_ne > 0d0) then
        conv_ne_u = calc_surface_converged(ne_u, old_ne_u, this%tol)
        conv_ne_d = calc_surface_converged(ne_d, old_ne_d, this%tol)
      else
        conv_ne_u = .true.
        conv_ne_d = .true.
      end if
      if(ratio_delta > 0d0) then
        conv_delta_real = calc_surface_converged(real(delta, 8), real(odelta, 8), this%tol)
        conv_delta_img = calc_surface_converged(aimag(delta), aimag(odelta), this%tol)
      else
        conv_delta_real = .true.
        conv_delta_img = .true.
      end if
      print *, conv_eg, conv_ne_u, conv_ne_d, conv_delta_real, conv_delta_img
      if( iter >= 50 .and. ALL([conv_eg, conv_ne_u, conv_ne_d, conv_delta_real, conv_delta_img])) exit
      oeg = eg

    end do !iter
    
    if(not(check_sorted_eg(eg))) then
      call sort_wf(eg,wf)
    end if
    
    if(ratio_ne > 0d0) then
      call calc_surface_ne(this, wf, eg, ne_u, ne_d, kbT)
    end if
      
    if(ratio_delta > 0d0) then
      call calc_surface_delta(this, wf, eg, delta, kbT)
    end if

  end subroutine car_parinello_calc_surface_general

  subroutine car_parinello_calc_mu_surface_at_ne(this, ne, eg, wf, mu, delta, ne_u, ne_d)
    implicit none
    CLASS(car_parrinello), INTENT(in)   :: this
    real(8), intent(in) :: ne
    REAL(8), DIMENSION(4*this%n), intent(out) :: eg
    complex(8),DIMENSION(4*this%n,4*this%n), intent(out)   :: wf
    real(8),intent(out)    :: mu
    complex(8), dimension(this%pl_1st%length()), intent(out) :: delta
    real(8), dimension(this%n), intent(out)   :: ne_u, ne_d
    real(8) :: mu_l, mu_h
    real(8), parameter :: tol_mu = 1d-6
    mu_l = -4d0
    mu_h = 4d0
    do while(abs(mu_h-mu_l) > tol_mu)
      mu = (mu_h + mu_l)*0.5d0
      delta(1:size(delta)/2) = 0.5d0
      delta(size(delta)/2+1:size(delta)) = -0.5d0
      ne_u = 0.5d0
      ne_d = 0.5d0
      !call diagonalize_surface_hamiltonian(this, eg, wf, mu, delta, ne_u, ne_d)
      call this%calc_surface_general(eg, wf, mu, delta, ne_u, ne_d, 0d0, 0d0)
      call this%calc_surface_general(eg, wf, mu, delta, ne_u, ne_d, 0.01d0, 0d0)
      call this%calc_surface_general(eg, wf, mu, delta, ne_u, ne_d, 0d0, 0.01d0)
      call this%calc_surface_general(eg, wf, mu, delta, ne_u, ne_d, 1d0, 1d0)
      !print *, mu, sum(ne_u + ne_d) 
      write(40, *) mu, sum(ne_u + ne_d) 
      if(sum(ne_u+ne_d) < ne) then
        mu_l = mu
      else
        mu_h = mu
      end if
    end do
  end subroutine car_parinello_calc_mu_surface_at_ne


  !  subroutine calc_bulk(this,mu,ne,fe_e) 
    subroutine calc_bulk(this,mu,ne) 
    implicit none
    CLASS(car_parrinello), INTENT(in)   :: this
    INTEGER               :: i ,j, k, l  ,f, r ,iter,iter_one, time_s,time_e,CountPerSec, CountMax, min,z
    integer, ALLOCATABLE  :: p(:, :)
    real(8),intent(in)    :: mu
    real(8),intent(out)   :: ne
  !  REAL(8),optional      :: fe_e
    REAL(8)               :: c , ocoef, h2, h2t2, h2t, h2t_z, h2u, h2jd, h2lm, kbT, norm, tol, sec, ratio, dix, diy, t
    REAL(8), DIMENSION(this%n) ::Sx, Sy, Sz , oSx ,oSy, oSz, xi, chi, density, dxi, dchi!, dos_u, dos_d
    REAL(8), DIMENSION(4*this%n) :: eg, oeg
    complex(8)            :: exi, echi!, temp(4*this%n)
    complex(8),DIMENSION(4*this%n,4*this%n)   :: wf, nwf, owf
    
    iter_one = 0
    tol = 1.0d-4
    kbT =1.0d-2
    c = 1.d0/(1 + this%h*this%eta*0.5d0)
    ocoef = this%h*this%eta*0.5d0 - 1.d0
    h2 = this%h**2
    h2t   = h2*this%t(1)
    h2t2  = h2*this%t(2)
    h2t_z = h2*this%t(3)
    h2u = h2*this%u
    h2jd = this%jd*h2
    h2lm = this%lm*h2
    eg = 0d0
    wf = 0.d0
    ratio = 1.d-1
    ALLOCATE (p(this%nh, 4))
    do i=1,4*this%n
      wf(i,i)=1.d0
    end do!i  
    do j = 1,this%n - this%nh
      k = this%ets(j)
      xi(k) =0d0
    !  chi(k) = 0d0
      do i = 1, this%nh
        xi(k)  =  xi(k) + this%hole(i, 2) * atan2(real(this%y(k)-this%y(this%hole(i, 1))),real(this%x(k)-this%x(this%hole(i, 1))))
      !  chi(k) =  xi(k) 
        z = site_to_z(this%hole(i,1),this%pshape)
        p(i, 1) = this%hole(i, 1) + 1 
        p(i, 2) = this%hole(i, 1) + this%nx(z)
        p(i, 3) = this%hole(i, 1) - 1
        p(i, 4) = this%hole(i, 1) - this%nx(z)
      enddo
      dxi(k)  = xi(k) + (pi-1.e-2)*dble(this%x(k)+this%y(k)) - xi(1) !   add antiferro
    !  dchi(k) = chi(k) - chi(1)
      Sx(k) = 0.5d0*cos(dxi(k)) 
      Sy(k) = 0.5d0*sin(dxi(k))
      Sz(k) = 0.0d0
      density(k) = 0.75
    enddo
    owf = wf

    call system_clock(time_s, CountPerSec, CountMax)
    DO iter = 1,100000
      DO r = 1, 4*this%n
        nwf(:,r) = 0.d0
        DO k = 1, this%n - this%nh
          i = this%ets(k)
          nwf(4*i-3, r) = (2.d0 + h2*eg(r) )*wf(4*i-3, r) + ocoef*owf(4*i-3, r) - h2u*(2d0/3d0)*(-(Sx(i) - ui*Sy(i))*wf(4*i-2, r) + (density(i) - Sz(i))*wf(4*i-3, r) ) + h2*mu*wf(4*i-3, r)
          nwf(4*i-2, r) = (2.d0 + h2*eg(r) )*wf(4*i-2, r) + ocoef*owf(4*i-2 ,r) - h2u*(2d0/3d0)*(-(Sx(i) + ui*Sy(i))*wf(4*i-3, r) + (density(i) + Sz(i))*wf(4*i-2, r) ) + h2*mu*wf(4*i-2, r)
          nwf(4*i-1, r) = (2.d0 + h2*eg(r) )*wf(4*i-1, r) + ocoef*owf(4*i-1, r) + h2u*(2d0/3d0)*(-(Sx(i) - ui*Sy(i))*wf(4*i  , r) + (density(i) - Sz(i))*wf(4*i-1, r) ) - h2*mu*wf(4*i-1, r)
          nwf(4*i  , r) = (2.d0 + h2*eg(r) )*wf(4*i  , r) + ocoef*owf(4*i  , r) + h2u*(2d0/3d0)*(-(Sx(i) + ui*Sy(i))*wf(4*i-1, r) + (density(i) + Sz(i))*wf(4*i  , r) ) - h2*mu*wf(4*i  , r)
          DO j = 1, SIZE(this%ns_1st(i)%s) 
            nwf(4*i-3, r) = nwf(4*i-3, r) + h2t*wf(4*this%ns_1st(i)%s(j)-3, r)
            nwf(4*i-2, r) = nwf(4*i-2, r) + h2t*wf(4*this%ns_1st(i)%s(j)-2, r)
            nwf(4*i-1, r) = nwf(4*i-1, r) - h2t*wf(4*this%ns_1st(i)%s(j)-1, r)
            nwf(4*i  , r) = nwf(4*i  , r) - h2t*wf(4*this%ns_1st(i)%s(j)  , r)
          ENDDO ! j
          DO j = 1, SIZE(this%ns_2nd(i)%s)
            nwf(4*i-3, r) = nwf(4*i-3, r) + h2t2*wf(4*this%ns_2nd(i)%s(j)-3, r)
            nwf(4*i-2, r) = nwf(4*i-2, r) + h2t2*wf(4*this%ns_2nd(i)%s(j)-2, r)
            nwf(4*i-1, r) = nwf(4*i-1, r) - h2t2*wf(4*this%ns_2nd(i)%s(j)-1, r)
            nwf(4*i  , r) = nwf(4*i  , r) - h2t2*wf(4*this%ns_2nd(i)%s(j)  , r)
          ENDDO ! j    
        end do!k
        ! interaction of sites across holes
        DO j = 1, this%pl_h%length()
          i = this%pl_h%value(j)%i
          f = this%pl_h%value(j)%f
          nwf(4*f-3, r) = nwf(4*f-3, r) - h2jd*0.5d0*( (Sx(i) - ui*Sy(i))*wf(4*f-2, r) + Sz(i)*wf(4*f-3, r) )
          nwf(4*f-2, r) = nwf(4*f-2, r) - h2jd*0.5d0*( (Sx(i) + ui*Sy(i))*wf(4*f-3, r) - Sz(i)*wf(4*f-2, r) )
          nwf(4*f-1, r) = nwf(4*f-1, r) + h2jd*0.5d0*( (Sx(i) - ui*Sy(i))*wf(4*f  , r) + Sz(i)*wf(4*f-1, r) )
          nwf(4*f  , r) = nwf(4*f  , r) + h2jd*0.5d0*( (Sx(i) + ui*Sy(i))*wf(4*f-1, r) - Sz(i)*wf(4*f  , r) )
          nwf(4*i-3, r) = nwf(4*i-3, r) - h2jd*0.5d0*( (Sx(f) - ui*Sy(f))*wf(4*i-2, r) + Sz(f)*wf(4*i-3, r) )
          nwf(4*i-2, r) = nwf(4*i-2, r) - h2jd*0.5d0*( (Sx(f) + ui*Sy(f))*wf(4*i-3, r) - Sz(f)*wf(4*i-2, r) )
          nwf(4*i-1, r) = nwf(4*i-1, r) + h2jd*0.5d0*( (Sx(f) - ui*Sy(f))*wf(4*i  , r) + Sz(f)*wf(4*i-1, r) )
          nwf(4*i  , r) = nwf(4*i  , r) + h2jd*0.5d0*( (Sx(f) + ui*Sy(f))*wf(4*i-1, r) - Sz(f)*wf(4*i  , r) )
        ENDDO
   !     DO i = 1, this%nh
   !         nwf(4*p(i,1)-3,r) = nwf(4*p(i,1)-3,r)  - h2lm * (- wf(4*p(i,4)-2,r) + ui*wf(4*p(i,2)-2,r))
   !         nwf(4*p(i,2)-3,r) = nwf(4*p(i,2)-3,r)  - h2lm * (- wf(4*p(i,3)-2,r) + ui*wf(4*p(i,1)-2,r))
   !         nwf(4*p(i,3)-3,r) = nwf(4*p(i,3)-3,r)  - h2lm * (- wf(4*p(i,2)-2,r) + ui*wf(4*p(i,4)-2,r))
   !         nwf(4*p(i,4)-3,r) = nwf(4*p(i,4)-3,r)  - h2lm * (- wf(4*p(i,1)-2,r) + ui*wf(4*p(i,3)-2,r))

    !        nwf(4*p(i,1)-2,r) = nwf(4*p(i,1)-2,r)  - h2lm * (- wf(4*p(i,4)-3,r) + ui*wf(4*p(i,2)-3,r))
    !        nwf(4*p(i,2)-2,r) = nwf(4*p(i,2)-2,r)  - h2lm * (- wf(4*p(i,3)-3,r) - ui*wf(4*p(i,1)-3,r))
    !        nwf(4*p(i,3)-2,r) = nwf(4*p(i,3)-2,r)  - h2lm * (  wf(4*p(i,2)-3,r) - ui*wf(4*p(i,4)-3,r))
    !        nwf(4*p(i,4)-2,r) = nwf(4*p(i,4)-2,r)  - h2lm * (  wf(4*p(i,1)-3,r) + ui*wf(4*p(i,3)-3,r))

    !        nwf(4*p(i,1)-1,r) = nwf(4*p(i,1)-1,r)  - h2lm * (- wf(4*p(i,4)  ,r) + ui*wf(4*p(i,2)  ,r))
    !        nwf(4*p(i,2)-1,r) = nwf(4*p(i,2)-1,r)  - h2lm * (- wf(4*p(i,3)  ,r) + ui*wf(4*p(i,1)  ,r))
    !        nwf(4*p(i,3)-1,r) = nwf(4*p(i,3)-1,r)  - h2lm * (- wf(4*p(i,2)  ,r) + ui*wf(4*p(i,4)  ,r))            
    !        nwf(4*p(i,4)-1,r) = nwf(4*p(i,4)-1,r)  - h2lm * (- wf(4*p(i,1)  ,r) + ui*wf(4*p(i,3)  ,r))
            
    !        nwf(4*p(i,1)  ,r) = nwf(4*p(i,1)  ,r)  - h2lm * (- wf(4*p(i,4)-1,r) + ui*wf(4*p(i,2)-1,r))
    !        nwf(4*p(i,2)  ,r) = nwf(4*p(i,2)  ,r)  - h2lm * (- wf(4*p(i,3)-1,r) - ui*wf(4*p(i,1)-1,r))
    !        nwf(4*p(i,3)  ,r) = nwf(4*p(i,3)  ,r)  - h2lm * (  wf(4*p(i,2)-1,r) - ui*wf(4*p(i,4)-1,r))
    !        nwf(4*p(i,4)  ,r) = nwf(4*p(i,4)  ,r)  - h2lm * (  wf(4*p(i,1)-1,r) + ui*wf(4*p(i,3)-1,r))
    !    ENDDO

        owf(:,r) = wf(:,r)
        wf(:,r) = c*nwf(:,r)
        oeg(r) = eg(r)
        eg(r) = 0.0d0
        DO l = 1, r-1
          wf(:, r) = wf(:,r) - DOT_PRODUCT(wf(:, l),c*nwf(:, r))*wf(:, l)
        ENDDO!l
        norm = SQRT(DOT_PRODUCT(wf(:,r), wf(:,r)))
        if(norm == 0) cycle
        wf(:,r) = wf(:,r)/norm
        do i = 1, this%pl_1st%length()
          k = this%pl_1st%value(i)%i
          j = this%pl_1st%value(i)%f
            eg(r) = eg(r) - this%t(1)*((conjg(wf(4*k-3,r))*wf(4*j-3,r) + conjg(wf(4*j-3,r))*wf(4*k-3,r))&                             
            &                        + (conjg(wf(4*k-2,r))*wf(4*j-2,r) + conjg(wf(4*j-2,r))*wf(4*k-2,r))&
            &                        - (conjg(wf(4*k-1,r))*wf(4*j-1,r) + conjg(wf(4*j-1,r))*wf(4*k-1,r))&
            &                        - (conjg(wf(4*k  ,r))*wf(4*j  ,r) + conjg(wf(4*j  ,r))*wf(4*k  ,r)))
        end do!i    
        do i = 1, this%pl_2nd%length()
          k = this%pl_2nd%value(i)%i
          j = this%pl_2nd%value(i)%f
          eg(r) = eg(r) -this%t(2)*((conjg(wf(4*k-3,r))*wf(4*j-3,r) + conjg(wf(4*j-3,r))*wf(4*k-3,r)) &
          &                       + (conjg(wf(4*k-2,r))*wf(4*j-2,r) + conjg(wf(4*j-2,r))*wf(4*k-2,r)) &
          &                       - (conjg(wf(4*k-1,r))*wf(4*j-1,r) + conjg(wf(4*j-1,r))*wf(4*k-1,r)) &
          &                       - (conjg(wf(4*k  ,r))*wf(4*j  ,r) + conjg(wf(4*j  ,r))*wf(4*k  ,r)))
        end do!i
        do k = 1, this%n -  this%nh
          i = this%ets(k)
          eg(r) = eg(r) - mu*(conjg(wf(4*i-3,r))*wf(4*i-3,r) + conjg(wf(4*i-2,r))*wf(4*i-2,r) &
      &                     - conjg(wf(4*i-1,r))*wf(4*i-1,r) - conjg(wf(4*i  ,r))*wf(4*i  ,r))&
      & + this%u*(2d0/3d0)*(((0.75d0 - Sz(i))*wf(4*i-3, r)-(Sx(i) - ui*Sy(i))*wf(4*i-2, r))*conjg(wf(4*i-3, r)) &
      &                   + ((0.75d0 + Sz(i))*wf(4*i-2, r)-(Sx(i) + ui*Sy(i))*wf(4*i-3, r))*conjg(wf(4*i-2, r)) &
      &                   - ((0.75d0 - Sz(i))*wf(4*i-1, r)-(Sx(i) - ui*Sy(i))*wf(4*i  , r))*conjg(wf(4*i-1, r)) &
      &                   - ((0.75d0 + Sz(i))*wf(4*i  , r)-(Sx(i) + ui*Sy(i))*wf(4*i-1, r))*conjg(wf(4*i  , r))) 
        end do!k
        do j = 1, this%pl_h%length()
          i = this%pl_h%value(j)%i
          f = this%pl_h%value(j)%f
         eg(r) = eg(r) + this%jd*(( (Sx(i) - ui*Sy(i))*wf(4*f-2, r) + Sz(i)*wf(4*f-3, r) )*conjg(wf(4*f-3, r))&
        &                       + ( (Sx(i) + ui*Sy(i))*wf(4*f-3, r) - Sz(i)*wf(4*f-2, r) )*conjg(wf(4*f-2, r))&
        &                       - ( (Sx(i) + ui*Sy(i))*wf(4*f-1, r) - Sz(i)*wf(4*f  , r) )*conjg(wf(4*f  , r))&
        &                       - ( (Sx(i) - ui*Sy(i))*wf(4*f  , r) + Sz(i)*wf(4*f-1, r) )*conjg(wf(4*f-1, r))&

        &                       + ( (Sx(f) - ui*Sy(f))*wf(4*i-2, r) + Sz(f)*wf(4*i-3, r) )*conjg(wf(4*i-3, r))&
        &                       + ( (Sx(f) + ui*Sy(f))*wf(4*i-3, r) - Sz(f)*wf(4*i-2, r) )*conjg(wf(4*i-2, r))&
        &                       - ( (Sx(f) - ui*Sy(f))*wf(4*i  , r) + Sz(f)*wf(4*i-1, r) )*conjg(wf(4*i-1, r))&
        &                       - ( (Sx(f) + ui*Sy(f))*wf(4*i-1, r) - Sz(f)*wf(4*i  , r) )*conjg(wf(4*i  , r)))
        end do!j
  !    do i = 1, this%nh 
  !           eg(r) = eg(r) &
  !          &  +  (- wf(4*p(i,4)-2,r) + ui*wf(4*p(i,2)-2,r))*conjg(wf(4*p(i,1)-3,r)) + (- wf(4*p(i,3)-2,r) + ui*wf(4*p(i,1)-2,r))*conjg(wf(4*p(i,2)-3,r))&
  !          &  +  (- wf(4*p(i,2)-2,r) + ui*wf(4*p(i,4)-2,r))*conjg(wf(4*p(i,3)-3,r)) + (- wf(4*p(i,1)-2,r) + ui*wf(4*p(i,3)-2,r))*conjg(wf(4*p(i,4)-3,r))& 
  !          &  +  (- wf(4*p(i,4)-3,r) + ui*wf(4*p(i,2)-3,r))*conjg(wf(4*p(i,1)-2,r)) + (- wf(4*p(i,3)-3,r) - ui*wf(4*p(i,1)-3,r))*conjg(wf(4*p(i,2)-2,r))&
  !          &  +  (  wf(4*p(i,2)-3,r) - ui*wf(4*p(i,4)-3,r))*conjg(wf(4*p(i,3)-2,r)) + (  wf(4*p(i,1)-3,r) + ui*wf(4*p(i,3)-3,r))*conjg(wf(4*p(i,4)-2,r))&
  !          &  +  (- wf(4*p(i,4)  ,r) + ui*wf(4*p(i,2)  ,r))*conjg(wf(4*p(i,1)-1,r)) + (- wf(4*p(i,3)  ,r) + ui*wf(4*p(i,1)  ,r))*conjg(wf(4*p(i,2)-1,r))&
  !          &  +  (- wf(4*p(i,2)  ,r) + ui*wf(4*p(i,4)  ,r))*conjg(wf(4*p(i,3)-1,r)) + (- wf(4*p(i,1)  ,r) + ui*wf(4*p(i,3)  ,r))*conjg(wf(4*p(i,4)-1,r))&
  !          &  +  (- wf(4*p(i,4)-1,r) + ui*wf(4*p(i,2)-1,r))*conjg(wf(4*p(i,1)  ,r)) + (- wf(4*p(i,3)-1,r) - ui*wf(4*p(i,1)-1,r))*conjg(wf(4*p(i,2)  ,r))&
  !          &  +  (  wf(4*p(i,2)-1,r) - ui*wf(4*p(i,4)-1,r))*conjg(wf(4*p(i,3)  ,r)) + (  wf(4*p(i,1)-1,r) + ui*wf(4*p(i,3)-1,r))*conjg(wf(4*p(i,4)  ,r))
  !    end do!i
      end do!r
      write(*,*)eg
      stop


      if(iter >= 50 .and. iter_one == 0) then
        do r = 1,4*this%n
          if((abs(oeg(r)) >  1.0d-10 .and. ABS((eg(r) - oeg(r))/oeg(r)) > tol)) goto 200
        end do!i
        iter_one = iter
      !  exit
       else if((iter - iter_one) >= 50 )then
        do i = 1, this%n
          oSx(i) = Sx(i)
          oSy(i) = Sy(i)
          oSz(i) = Sz(i)
          Sx(i) = 0.d0
          Sy(i) = 0.d0
          Sz(i) = 0.d0
          density(i) = 0.d0
          do r = 1, 4*this%n
            if (eg(r) >= 0.d0)then
              Sx(i) = Sx(i) -  real(wf(4*i-1, r)*conjg(wf(4*i  , r)))!+ wf(4*i-3, r)*conjg(wf(4*i-2, r))) 
              Sy(i) = Sy(i) - aimag(wf(4*i-1, r)*conjg(wf(4*i  , r)))!+ wf(4*i-3, r)*conjg(wf(4*i-2, r)))
              Sz(i) = Sz(i) +  real(wf(4*i-1, r)*conjg(wf(4*i-1, r)) - wf(4*i  , r)*conjg(wf(4*i  , r)))/2! + wf(4*i-3, r)*conjg(wf(4*i-3, r)) - wf(4*i-2, r)*conjg(wf(4*i-2, r)))/2
              density(i) = density(i) + real(wf(4*i-1, r)*conjg(wf(4*i-1, r)) + wf(4*i  , r)*conjg(wf(4*i  , r)) + wf(4*i-3, r)*conjg(wf(4*i-3, r)) + wf(4*i-2, r)*conjg(wf(4*i-2, r)))/2
            end if
          end do
        end do
        do i = 1,this%n 
           if((abs(oSx(i))  >  1.0d-10 .and. ABS((Sx(i) - oSx(i))/oSx(i)) > tol) &
      &   .or. (abs(oSy(i))  >  1.0d-10 .and. ABS((Sy(i) - oSy(i))/oSy(i)) > tol) &
      &   .or. (abs(oSz(i))  >  1.0d-10 .and. ABS((Sz(i) - oSz(i))/oSz(i)) > tol))  goto 200
        end do!i
        exit
      end if
    200 continue
    end do !iter

    call system_clock(time_e)
    xi = Sx
    do i = 1,this%n
      IF (xi(i) == 0.d0) xi(i) = 1.d-10
      xi(i) = -0.5d0*(SIGN(1.0d0, xi(i)) - 1.0d0)*pi + ATAN(Sy(i)/xi(i))
    end do!i
    do i = 1, size(this%pshape(1,:))
      do j = 1, this%pshape(1,i) * this%pshape(2,i)
        do l = 1, size(this%hole(:,1))
          if (j == this%hole(l,1)) goto 100
        end do!l
        dix = 0.40d0 * cos(xi(j))
        diy = 0.40d0 * sin(xi(j))
        write(60, '(6(f18.15, 1x))') dble(mod((j-1),this%pshape(1,i))+1)-dix, dble((j-1)/this%pshape(1,i)+1)-diy, dble(i),2.0d0*dix,2.0d0*diy, 0d0
        100 continue
      end do!j 
    end do!i
    write(60,'(2(a12,/))')' exit       ',"sp '-' w vec"
    
    ne = 0.d0
    do r = 1, 4*this%n
      write(40,*)r,eg(r)
      if(eg(r)  > 0.d0)then
        do i = 1, this%n
        if(( eg(r)) < 50*kbt) ne = ne  + (conjg(wf(4*i-3, r))*wf(4*i-3, r) + conjg(wf(4*i-2, r))*wf(4*i-2, r))/(exp(( eg(r))/kbT)+1.d0)
        if((-eg(r)) < 50*kbt) ne = ne  + (conjg(wf(4*i-1, r))*wf(4*i-1, r) + conjg(wf(4*i  , r))*wf(4*i  , r))/(exp((-eg(r))/kbT)+1.d0)
        end do!i
      end if
    end do!r
    write(40,'(2(a10,/))')' exit     ',"p '-' w l"
    sec = real(time_e - time_s)/CountPerSec
  !  min = int(sec/60)
    write(20,'(a7,f13.5,a3,2(a10,i7),2(a11,f13.10))') 'time = ',sec ,'(s) ','    iter = ',iter,' iter_one = ',iter_one,'mu = ',mu,'ne = ',ne
    write(10,*)mu,ne,iter
  end subroutine calc_bulk





    subroutine calc_sb(this,mu,ne) 
    implicit none
    CLASS(car_parrinello), INTENT(in)   :: this
    INTEGER               :: i ,j, k, l , m,f, r ,iter,iter_two,iter_one, time_s,time_e,CountPerSec, CountMax, min,z
    integer, ALLOCATABLE  :: p(:, :)
    real(8),intent(in)    :: mu
    real(8),intent(out)   :: ne
  !  REAL(8),optional      :: fe_e
    REAL(8)               :: c , ocoef, h2, h2t2, h2t, h2t_z, h2u, h2jd, h2lm, kbT, norm, tol, sec, ratio, dix, diy, t
    REAL(8), DIMENSION(this%n) ::Sx, Sy, Sz , oSx ,oSy, oSz, xi, chi, density, dxi, dchi!, dos_u, dos_d
    REAL(8), DIMENSION(4*this%n) :: tmp_eg, eg, oeg
    complex(8)            :: exi, echi!, temp(4*this%n)
    complex(8),DIMENSION(4*this%n,4*this%n)   :: wf, nwf, owf, temp
    LOGICAL  temp2(4*this%n) ,dflag ,l_n
   complex(8)            ::  delta(size(this%pl_1st%value))
    REAL(8), DIMENSION(this%nx(1)*this%ny(1)) ::dos_u,dos_d, ne_u, ne_d, old_u, old_d
    
    iter_one = 0
    tol = 1.0d-4
    kbT =1.0d-2
    c = 1.d0/(1 + this%h*this%eta*0.5d0)
    ocoef = this%h*this%eta*0.5d0 - 1.d0
    h2 = this%h**2
    h2t   = h2*this%t(1)
    h2t2  = h2*this%t(2)
    h2t_z = h2*this%t(3)
    h2u = h2*this%u
    h2jd = this%jd*h2
    h2lm = this%lm*h2
    eg = 0d0
    wf = 0.d0
    ratio = 1.d-1
    ALLOCATE (p(this%nh, 4))
    do i=1,4*this%n
      wf(i,i)=1.d0
    end do!i  
    do j = 1,this%n - this%nh
      k = this%ets(j)
      xi(k) =0d0
    !  chi(k) = 0d0
      do i = 1, this%nh
        xi(k)  =  xi(k) + this%hole(i, 2) * atan2(real(this%y(k)-this%y(this%hole(i, 1))),real(this%x(k)-this%x(this%hole(i, 1))))
      !  chi(k) =  xi(k) 
        z = site_to_z(this%hole(i,1),this%pshape)
        p(i, 1) = this%hole(i, 1) + 1 
        p(i, 2) = this%hole(i, 1) + this%nx(z)
        p(i, 3) = this%hole(i, 1) - 1
        p(i, 4) = this%hole(i, 1) - this%nx(z)
      enddo
      dxi(k)  = xi(k) + (pi-1.e-2)*dble(this%x(k)+this%y(k)) - xi(1) !   add antiferro
    !  dchi(k) = chi(k) - chi(1)
      Sx(k) = 0.5d0*cos(dxi(k)) 
      Sy(k) = 0.5d0*sin(dxi(k))
      Sz(k) = 0.0d0
      density(k) = 0.75
    enddo
    owf = wf

    call system_clock(time_s, CountPerSec, CountMax)
    DO iter = 1,1000
      nwf = 0.d0
      DO r = 1, 4*this%n
        DO k = 1, int(this%ne(1))
          i = this%ets(k)
          nwf(4*i - 3, r) = (2.d0 + h2*oeg(r) )*wf(4*i - 3, r) + ocoef*owf(4*i - 3, r) + h2*mu*wf(4*i - 3, r)
          nwf(4*i - 2, r) = (2.d0 + h2*oeg(r) )*wf(4*i - 2, r) + ocoef*owf(4*i - 2 ,r) + h2*mu*wf(4*i - 2, r)
          nwf(4*i - 1, r) = (2.d0 + h2*oeg(r) )*wf(4*i - 1, r) + ocoef*owf(4*i - 1, r) - h2*mu*wf(4*i - 1, r)
          nwf(4*i    , r) = (2.d0 + h2*oeg(r) )*wf(4*i    , r) + ocoef*owf(4*i    , r) - h2*mu*wf(4*i    , r)
          DO j = 1, SIZE(this%ns_1st(i)%s)       
            nwf(4*i - 3, r) = nwf(4*i - 3, r) + h2t*wf(4*this%ns_1st(i)%s(j) - 3, r)*(1.d0 - ne_d(i))*(1.d0 - ne_d(this%ns_1st(i)%s(j)))
            nwf(4*i - 2, r) = nwf(4*i - 2, r) + h2t*wf(4*this%ns_1st(i)%s(j) - 2, r)*(1.d0 - ne_u(i))*(1.d0 - ne_u(this%ns_1st(i)%s(j)))
            nwf(4*i - 1, r) = nwf(4*i - 1, r) - h2t*wf(4*this%ns_1st(i)%s(j) - 1, r)*(1.d0 - ne_d(i))*(1.d0 - ne_d(this%ns_1st(i)%s(j)))
            nwf(4*i    , r) = nwf(4*i    , r) - h2t*wf(4*this%ns_1st(i)%s(j)    , r)*(1.d0 - ne_u(i))*(1.d0 - ne_u(this%ns_1st(i)%s(j)))
          ENDDO ! j
          DO j = 1, SIZE(this%ns_2nd(i)%s)
            nwf(4*i - 3, r) = nwf(4*i - 3, r) + h2t2*wf(4*this%ns_2nd(i)%s(j) - 3, r)*(1.d0 - ne_d(i))*(1.d0 - ne_d(this%ns_2nd(i)%s(j))) 
            nwf(4*i - 2, r) = nwf(4*i - 2, r) + h2t2*wf(4*this%ns_2nd(i)%s(j) - 2, r)*(1.d0 - ne_u(i))*(1.d0 - ne_u(this%ns_2nd(i)%s(j)))
            nwf(4*i - 1, r) = nwf(4*i - 1, r) - h2t2*wf(4*this%ns_2nd(i)%s(j) - 1, r)*(1.d0 - ne_d(i))*(1.d0 - ne_d(this%ns_2nd(i)%s(j)))
            nwf(4*i    , r) = nwf(4*i    , r) - h2t2*wf(4*this%ns_2nd(i)%s(j)    , r)*(1.d0 - ne_u(i))*(1.d0 - ne_u(this%ns_2nd(i)%s(j)))
          ENDDO ! j    
        end do!k
        if(dflag ) then
      !  do i = 1, this%pl_1st%length()
        DO k = 1, int(this%ne(1))
          i = this%ets(k)
          nwf(4*i - 3, r) = nwf(4*i - 3, r)  -       delta(i) * h2*wf(4*j    , r)
          nwf(4*i - 2, r) = nwf(4*i - 2, r)  -       delta(i) * h2*wf(4*j - 1, r)
          nwf(4*i - 1, r) = nwf(4*i - 1, r)  - conjg(delta(i))* h2*wf(4*j - 2, r)
          nwf(4*i    , r) = nwf(4*i    , r)  - conjg(delta(i))* h2*wf(4*j - 3, r) 
        !  k = this%pl_1st%value(i)%i
        !  j = this%pl_1st%value(i)%f
        !  nwf(4*k - 3, r) = nwf(4*k - 3, r)  -       delta(i) * h2*wf(4*j    , r)
        !  nwf(4*k - 2, r) = nwf(4*k - 2, r)  -       delta(i) * h2*wf(4*j - 1, r)
        !  nwf(4*k - 1, r) = nwf(4*k - 1, r)  - conjg(delta(i))* h2*wf(4*j - 2, r)
        !  nwf(4*k    , r) = nwf(4*k    , r)  - conjg(delta(i))* h2*wf(4*j - 3, r) 
        !  nwf(4*j - 2, r) = nwf(4*j - 2, r)  -       delta(i) * h2*wf(4*k - 1, r)
        !  nwf(4*j - 1, r) = nwf(4*j - 1, r)  - conjg(delta(i))* h2*wf(4*k - 2, r)
        !  nwf(4*j    , r) = nwf(4*j    , r)  - conjg(delta(i))* h2*wf(4*k - 3, r) 
        ENDDO ! i
        end if
      end do!r

      eg = 0.d0
      tmp_eg = 0d0
      do i = 1, this%pl_1st%length()
        k = this%pl_1st%value(i)%i
        j = this%pl_1st%value(i)%f
        do r = 1, 4*this%n
            tmp_eg(r) = tmp_eg(r) + this%t(1)*( (conjg(wf(4*k-3,r))*wf(4*j-3,r) + conjg(wf(4*j-3,r))*wf(4*k-3,r))*(1.d0 - ne_d(k))*(1.d0 - ne_d(j)) &                             
            &                                        + (conjg(wf(4*k-2,r))*wf(4*j-2,r) + conjg(wf(4*j-2,r))*wf(4*k-2,r))*(1.d0 - ne_u(k))*(1.d0 - ne_u(j)) &
            &                                        - (conjg(wf(4*k-1,r))*wf(4*j-1,r) + conjg(wf(4*j-1,r))*wf(4*k-1,r))*(1.d0 - ne_d(k))*(1.d0 - ne_d(j))&
            &                                        - (conjg(wf(4*k  ,r))*wf(4*j  ,r) + conjg(wf(4*j  ,r))*wf(4*k  ,r))*(1.d0 - ne_u(k))*(1.d0 - ne_u(j)) )
          if(dflag)  tmp_eg(r) = tmp_eg(r) + delta(i)*( conjg(wf(4*k-3,r))*wf(4*j  ,r) + conjg(wf(4*j-2,r))*wf(4*k-1,r) &
          &                                           + conjg(wf(4*j-3,r))*wf(4*k  ,r) + conjg(wf(4*k-2,r))*wf(4*j-1,r))&
          &                         + conjg(delta(i))*( conjg(wf(4*k-1,r))*wf(4*j-2,r) + conjg(wf(4*j  ,r))*wf(4*k-3,r) &
          &                                          +  conjg(wf(4*j-1,r))*wf(4*k-2,r) + conjg(wf(4*k  ,r))*wf(4*j-3,r))
        end do!r
      end do!i
      eg(1:4*this%n) = eg(1:4*this%n) - tmp_eg(1:4*this%n) 
      tmp_eg = 0d0
      do i = 1, this%pl_2nd%length()
        k = this%pl_2nd%value(i)%i
        j = this%pl_2nd%value(i)%f
        do r = 1, 4*this%n
          tmp_eg(r) = tmp_eg(r) + (conjg(wf(4*k-3,r))*wf(4*j-3,r) + conjg(wf(4*j-3,r))*wf(4*k-3,r))*(1.d0 - ne_d(k))*(1.d0 - ne_d(j)) &
          &                                  + (conjg(wf(4*k-2,r))*wf(4*j-2,r) + conjg(wf(4*j-2,r))*wf(4*k-2,r))*(1.d0 - ne_u(k))*(1.d0 - ne_u(j)) &
          &                                  - (conjg(wf(4*k-1,r))*wf(4*j-1,r) + conjg(wf(4*j-1,r))*wf(4*k-1,r))*(1.d0 - ne_d(k))*(1.d0 - ne_d(j)) &
          &                                  - (conjg(wf(4*k  ,r))*wf(4*j  ,r) + conjg(wf(4*j  ,r))*wf(4*k  ,r))*(1.d0 - ne_u(k))*(1.d0 - ne_u(j))  
        end do!r
      end do!i
      eg(1:4*this%n) = eg(1:4*this%n) - tmp_eg(1:4*this%n) * this%t(2) 
      tmp_eg = 0d0
      do r = 1, 4*this%n 
        do k = 1, this%n -  this%nh
          j = this%ets(k)
          tmp_eg(r) = tmp_eg(r) + conjg(wf(4*j-3,r))*wf(4*j-3,r) + conjg(wf(4*j-2,r))*wf(4*j-2,r) &
          &                     - conjg(wf(4*j-1,r))*wf(4*j-1,r) - conjg(wf(4*j  ,r))*wf(4*j  ,r) 
        end do!j

      end do!r
      eg(1:4*this%n) = eg(1:4*this%n) - tmp_eg(1:4*this%n) * mu

      owf = wf
      temp = c*nwf
      DO k = 1, 4*this%n
        DO l = 1, k - 1
          temp(:, k) = temp(:, k) - DOT_PRODUCT(temp(:, l),c*nwf(:, k))*temp(:, l)
        ENDDO!l
        norm = SQRT(DOT_PRODUCT(temp(:,k), temp(:,k)))
        IF (norm == 0.d0) THEN
          temp(:,k) = 0.d0
        ELSE
          temp(:,k) = temp(:,k)/norm
        ENDIF
      ENDDO!k
      wf = temp

      temp2(:) = .FALSE.
      do i = 1,4*this%n
        if(abs(oeg(i))  <  tol .or. ABS((eg(i) - oeg(i))/oeg(i)) < tol)  temp2(i) = .true.
        if(abs(oeg(i)) >= 1.0d-8)write(60,'(2i,2f18.8,l,f18.8)')iter,i ,eg(i) ,oeg(i), temp2(i) ,ABS((eg(i) - oeg(i))/oeg(i))
        if(abs(oeg(i)) <  1.0d-8)write(60,'(2i,2f18.8,l)')iter,i ,eg(i) ,oeg(i), temp2(i) 
        oeg(i) = eg(i)
      end do!i
      if(dflag .and. (iter - iter_two) >= 50 .and. ALL(temp2)) EXIT
      if(not(l_n) .and. iter >= 50 .and. ALL(temp2) ) then
          l_n = .true.
          iter_one = iter
      end if
      if(not(dflag) .and. l_n .and. (iter - iter_one)  >= 50 .and. ALL(temp2) )then
         dflag = .true.
         iter_two = iter
      end if

      ne_u = 0.d0
      ne_d = 0.d0
      if(l_n)then
        do r = 1, 4*this%n
          if(eg(r)  > 0.d0)then
          do i = 1, this%nx(1)*this%ny(1)
            ne_u(i) = ne_u(i) + fermi( eg(r), kbT)*conjg(wf(4*i-3, r))*wf(4*i-3, r) + fermi(-eg(r), kbT)*conjg(wf(4*i-1, r))*wf(4*i-1, r)
            ne_d(i) = ne_d(i) + fermi( eg(r), kbT)*conjg(wf(4*i-2, r))*wf(4*i-2, r) + fermi(-eg(r), kbT)*conjg(wf(4*i  , r))*wf(4*i  , r) 
          end do!i
          end if
        end do!a
      ne_u = ratio*ne_u + (1d0-ratio)*old_u
      ne_d = ratio*ne_d + (1d0-ratio)*old_d
      old_u = ne_u
      old_d = ne_d
      end if

      if(dflag)then
        delta  = 0.d0
        do k = 1, this%pl_1st%length()
          do r = 1, 4*this%n
            if(eg(r) >= 0.d0)then
            i = this%pl_1st%value(k)%i
            j = this%pl_1st%value(k)%f
            delta(k) = this%jd*(wf(4*i-3, r)*conjg(wf(4*j  , r)) + wf(4*j-2, r)*conjg(wf(4*i-1, r)) &
            &                 + wf(4*i-2, r)*conjg(wf(4*j-1, r)) + wf(4*j-3, r)*conjg(wf(4*i  , r))) + delta(k)
            end if
          end do!a
        end do!k
      end if














      DO r = 1, 4*this%n
        nwf(:,r) = 0.d0
      do m= 1,size(this%ne)-1
        DO k = 1, int(this%ne(m))
          i = this%ets(k)
          nwf(4*i-3, r) = (2.d0 + h2*eg(r) )*wf(4*i-3, r) + ocoef*owf(4*i-3, r) + h2*mu*wf(4*i-3, r) - h2u*(2d0/3d0)*(-(Sx(i) - ui*Sy(i))*wf(4*i-2, r) + (density(i) - Sz(i))*wf(4*i-3, r) )
          nwf(4*i-2, r) = (2.d0 + h2*eg(r) )*wf(4*i-2, r) + ocoef*owf(4*i-2 ,r) + h2*mu*wf(4*i-2, r) - h2u*(2d0/3d0)*(-(Sx(i) + ui*Sy(i))*wf(4*i-3, r) + (density(i) + Sz(i))*wf(4*i-2, r) )
          nwf(4*i-1, r) = (2.d0 + h2*eg(r) )*wf(4*i-1, r) + ocoef*owf(4*i-1, r) - h2*mu*wf(4*i-1, r) + h2u*(2d0/3d0)*(-(Sx(i) - ui*Sy(i))*wf(4*i  , r) + (density(i) - Sz(i))*wf(4*i-1, r) )
          nwf(4*i  , r) = (2.d0 + h2*eg(r) )*wf(4*i  , r) + ocoef*owf(4*i  , r) - h2*mu*wf(4*i  , r) + h2u*(2d0/3d0)*(-(Sx(i) + ui*Sy(i))*wf(4*i-1, r) + (density(i) + Sz(i))*wf(4*i  , r) )
          DO j = 1, SIZE(this%ns_1st(i)%s) 
            nwf(4*i-3, r) = nwf(4*i-3, r) + h2t*wf(4*this%ns_1st(i)%s(j)-3, r)
            nwf(4*i-2, r) = nwf(4*i-2, r) + h2t*wf(4*this%ns_1st(i)%s(j)-2, r)
            nwf(4*i-1, r) = nwf(4*i-1, r) - h2t*wf(4*this%ns_1st(i)%s(j)-1, r)
            nwf(4*i  , r) = nwf(4*i  , r) - h2t*wf(4*this%ns_1st(i)%s(j)  , r)
          ENDDO ! j
          DO j = 1, SIZE(this%ns_2nd(i)%s)
            nwf(4*i-3, r) = nwf(4*i-3, r) + h2t2*wf(4*this%ns_2nd(i)%s(j)-3, r)
            nwf(4*i-2, r) = nwf(4*i-2, r) + h2t2*wf(4*this%ns_2nd(i)%s(j)-2, r)
            nwf(4*i-1, r) = nwf(4*i-1, r) - h2t2*wf(4*this%ns_2nd(i)%s(j)-1, r)
            nwf(4*i  , r) = nwf(4*i  , r) - h2t2*wf(4*this%ns_2nd(i)%s(j)  , r)
          ENDDO ! j    
        end do!k
        ! interaction of sites across holes
        DO j = 1, this%pl_h%length()
          i = this%pl_h%value(j)%i
          f = this%pl_h%value(j)%f
          nwf(4*f-3, r) = nwf(4*f-3, r) - h2jd*0.5d0*( (Sx(i) - ui*Sy(i))*wf(4*f-2, r) + Sz(i)*wf(4*f-3, r) )
          nwf(4*f-2, r) = nwf(4*f-2, r) - h2jd*0.5d0*( (Sx(i) + ui*Sy(i))*wf(4*f-3, r) - Sz(i)*wf(4*f-2, r) )
          nwf(4*f-1, r) = nwf(4*f-1, r) + h2jd*0.5d0*( (Sx(i) - ui*Sy(i))*wf(4*f  , r) + Sz(i)*wf(4*f-1, r) )
          nwf(4*f  , r) = nwf(4*f  , r) + h2jd*0.5d0*( (Sx(i) + ui*Sy(i))*wf(4*f-1, r) - Sz(i)*wf(4*f  , r) )
          nwf(4*i-3, r) = nwf(4*i-3, r) - h2jd*0.5d0*( (Sx(f) - ui*Sy(f))*wf(4*i-2, r) + Sz(f)*wf(4*i-3, r) )
          nwf(4*i-2, r) = nwf(4*i-2, r) - h2jd*0.5d0*( (Sx(f) + ui*Sy(f))*wf(4*i-3, r) - Sz(f)*wf(4*i-2, r) )
          nwf(4*i-1, r) = nwf(4*i-1, r) + h2jd*0.5d0*( (Sx(f) - ui*Sy(f))*wf(4*i  , r) + Sz(f)*wf(4*i-1, r) )
          nwf(4*i  , r) = nwf(4*i  , r) + h2jd*0.5d0*( (Sx(f) + ui*Sy(f))*wf(4*i-1, r) - Sz(f)*wf(4*i  , r) )
        ENDDO

        owf(:,r) = wf(:,r)
        wf(:,r) = c*nwf(:,r)
        oeg(r) = eg(r)
        eg(r) = 0.0d0
        DO l = 1, r-1
          wf(:, r) = wf(:,r) - DOT_PRODUCT(wf(:, l),c*nwf(:, r))*wf(:, l)
        ENDDO!l
        norm = SQRT(DOT_PRODUCT(wf(:,r), wf(:,r)))
        if(norm == 0) cycle
        wf(:,r) = wf(:,r)/norm
        do i = 1, this%pl_1st%length()
          k = this%pl_1st%value(i)%i
          j = this%pl_1st%value(i)%f
            eg(r) = eg(r) - this%t(1)*((conjg(wf(4*k-3,r))*wf(4*j-3,r) + conjg(wf(4*j-3,r))*wf(4*k-3,r))&                             
            &                        + (conjg(wf(4*k-2,r))*wf(4*j-2,r) + conjg(wf(4*j-2,r))*wf(4*k-2,r))&
            &                        - (conjg(wf(4*k-1,r))*wf(4*j-1,r) + conjg(wf(4*j-1,r))*wf(4*k-1,r))&
            &                        - (conjg(wf(4*k  ,r))*wf(4*j  ,r) + conjg(wf(4*j  ,r))*wf(4*k  ,r)))
        end do!i    
        do i = 1, this%pl_2nd%length()
          k = this%pl_2nd%value(i)%i
          j = this%pl_2nd%value(i)%f
          eg(r) = eg(r) -this%t(2)*((conjg(wf(4*k-3,r))*wf(4*j-3,r) + conjg(wf(4*j-3,r))*wf(4*k-3,r)) &
          &                       + (conjg(wf(4*k-2,r))*wf(4*j-2,r) + conjg(wf(4*j-2,r))*wf(4*k-2,r)) &
          &                       - (conjg(wf(4*k-1,r))*wf(4*j-1,r) + conjg(wf(4*j-1,r))*wf(4*k-1,r)) &
          &                       - (conjg(wf(4*k  ,r))*wf(4*j  ,r) + conjg(wf(4*j  ,r))*wf(4*k  ,r)))
        end do!i
        do k = 1, this%n -  this%nh
          i = this%ets(k)
          eg(r) = eg(r) - mu*(conjg(wf(4*i-3,r))*wf(4*i-3,r) + conjg(wf(4*i-2,r))*wf(4*i-2,r) &
      &                     - conjg(wf(4*i-1,r))*wf(4*i-1,r) - conjg(wf(4*i  ,r))*wf(4*i  ,r))&
      & + this%u*(2d0/3d0)*(((0.75d0 - Sz(i))*wf(4*i-3, r)-(Sx(i) - ui*Sy(i))*wf(4*i-2, r))*conjg(wf(4*i-3, r)) &
      &                   + ((0.75d0 + Sz(i))*wf(4*i-2, r)-(Sx(i) + ui*Sy(i))*wf(4*i-3, r))*conjg(wf(4*i-2, r)) &
      &                   - ((0.75d0 - Sz(i))*wf(4*i-1, r)-(Sx(i) - ui*Sy(i))*wf(4*i  , r))*conjg(wf(4*i-1, r)) &
      &                   - ((0.75d0 + Sz(i))*wf(4*i  , r)-(Sx(i) + ui*Sy(i))*wf(4*i-1, r))*conjg(wf(4*i  , r))) 
        end do!k
        do j = 1, this%pl_h%length()
          i = this%pl_h%value(j)%i
          f = this%pl_h%value(j)%f
         eg(r) = eg(r) + this%jd*(( (Sx(i) - ui*Sy(i))*wf(4*f-2, r) + Sz(i)*wf(4*f-3, r) )*conjg(wf(4*f-3, r))&
        &                       + ( (Sx(i) + ui*Sy(i))*wf(4*f-3, r) - Sz(i)*wf(4*f-2, r) )*conjg(wf(4*f-2, r))&
        &                       - ( (Sx(i) + ui*Sy(i))*wf(4*f-1, r) - Sz(i)*wf(4*f  , r) )*conjg(wf(4*f  , r))&
        &                       - ( (Sx(i) - ui*Sy(i))*wf(4*f  , r) + Sz(i)*wf(4*f-1, r) )*conjg(wf(4*f-1, r))&

        &                       + ( (Sx(f) - ui*Sy(f))*wf(4*i-2, r) + Sz(f)*wf(4*i-3, r) )*conjg(wf(4*i-3, r))&
        &                       + ( (Sx(f) + ui*Sy(f))*wf(4*i-3, r) - Sz(f)*wf(4*i-2, r) )*conjg(wf(4*i-2, r))&
        &                       - ( (Sx(f) - ui*Sy(f))*wf(4*i  , r) + Sz(f)*wf(4*i-1, r) )*conjg(wf(4*i-1, r))&
        &                       - ( (Sx(f) + ui*Sy(f))*wf(4*i-1, r) - Sz(f)*wf(4*i  , r) )*conjg(wf(4*i  , r)))
        end do!j
  !    do i = 1, this%nh 
  !           eg(r) = eg(r) &
  !          &  +  (- wf(4*p(i,4)-2,r) + ui*wf(4*p(i,2)-2,r))*conjg(wf(4*p(i,1)-3,r)) + (- wf(4*p(i,3)-2,r) + ui*wf(4*p(i,1)-2,r))*conjg(wf(4*p(i,2)-3,r))&
  !          &  +  (- wf(4*p(i,2)-2,r) + ui*wf(4*p(i,4)-2,r))*conjg(wf(4*p(i,3)-3,r)) + (- wf(4*p(i,1)-2,r) + ui*wf(4*p(i,3)-2,r))*conjg(wf(4*p(i,4)-3,r))& 
  !          &  +  (- wf(4*p(i,4)-3,r) + ui*wf(4*p(i,2)-3,r))*conjg(wf(4*p(i,1)-2,r)) + (- wf(4*p(i,3)-3,r) - ui*wf(4*p(i,1)-3,r))*conjg(wf(4*p(i,2)-2,r))&
  !          &  +  (  wf(4*p(i,2)-3,r) - ui*wf(4*p(i,4)-3,r))*conjg(wf(4*p(i,3)-2,r)) + (  wf(4*p(i,1)-3,r) + ui*wf(4*p(i,3)-3,r))*conjg(wf(4*p(i,4)-2,r))&
  !          &  +  (- wf(4*p(i,4)  ,r) + ui*wf(4*p(i,2)  ,r))*conjg(wf(4*p(i,1)-1,r)) + (- wf(4*p(i,3)  ,r) + ui*wf(4*p(i,1)  ,r))*conjg(wf(4*p(i,2)-1,r))&
  !          &  +  (- wf(4*p(i,2)  ,r) + ui*wf(4*p(i,4)  ,r))*conjg(wf(4*p(i,3)-1,r)) + (- wf(4*p(i,1)  ,r) + ui*wf(4*p(i,3)  ,r))*conjg(wf(4*p(i,4)-1,r))&
  !          &  +  (- wf(4*p(i,4)-1,r) + ui*wf(4*p(i,2)-1,r))*conjg(wf(4*p(i,1)  ,r)) + (- wf(4*p(i,3)-1,r) - ui*wf(4*p(i,1)-1,r))*conjg(wf(4*p(i,2)  ,r))&
  !          &  +  (  wf(4*p(i,2)-1,r) - ui*wf(4*p(i,4)-1,r))*conjg(wf(4*p(i,3)  ,r)) + (  wf(4*p(i,1)-1,r) + ui*wf(4*p(i,3)-1,r))*conjg(wf(4*p(i,4)  ,r))
  !    end do!i
      end do!l
      end do!r

      if(iter >= 50 .and. iter_one == 0) then
        do r = 1,4*this%n
          if((abs(oeg(r)) >  1.0d-10 .and. ABS((eg(r) - oeg(r))/oeg(r)) > tol)) goto 200
        end do!i
        iter_one = iter
      !  exit
       else if((iter - iter_one) >= 50 )then
        do i = 1, this%n
          oSx(i) = Sx(i)
          oSy(i) = Sy(i)
          oSz(i) = Sz(i)
          Sx(i) = 0.d0
          Sy(i) = 0.d0
          Sz(i) = 0.d0
          density(i) = 0.d0
          do r = 1, 4*this%n
            if (eg(r) >= 0.d0)then
              Sx(i) = Sx(i) -  real(wf(4*i-1, r)*conjg(wf(4*i  , r)))!+ wf(4*i-3, r)*conjg(wf(4*i-2, r))) 
              Sy(i) = Sy(i) - aimag(wf(4*i-1, r)*conjg(wf(4*i  , r)))!+ wf(4*i-3, r)*conjg(wf(4*i-2, r)))
              Sz(i) = Sz(i) +  real(wf(4*i-1, r)*conjg(wf(4*i-1, r)) - wf(4*i  , r)*conjg(wf(4*i  , r)))/2! + wf(4*i-3, r)*conjg(wf(4*i-3, r)) - wf(4*i-2, r)*conjg(wf(4*i-2, r)))/2
              density(i) = density(i) + real(wf(4*i-1, r)*conjg(wf(4*i-1, r)) + wf(4*i  , r)*conjg(wf(4*i  , r)) + wf(4*i-3, r)*conjg(wf(4*i-3, r)) + wf(4*i-2, r)*conjg(wf(4*i-2, r)))/2
            end if
          end do
        end do
        do i = 1,this%n 
           if((abs(oSx(i))  >  1.0d-10 .and. ABS((Sx(i) - oSx(i))/oSx(i)) > tol) &
      &   .or. (abs(oSy(i))  >  1.0d-10 .and. ABS((Sy(i) - oSy(i))/oSy(i)) > tol) &
      &   .or. (abs(oSz(i))  >  1.0d-10 .and. ABS((Sz(i) - oSz(i))/oSz(i)) > tol))  goto 200
        end do!i
        exit
      end if
    200 continue
    end do !iter

    call system_clock(time_e)
    xi = Sx
    do i = 1,this%n
      IF (xi(i) == 0.d0) xi(i) = 1.d-10
      xi(i) = -0.5d0*(SIGN(1.0d0, xi(i)) - 1.0d0)*pi + ATAN(Sy(i)/xi(i))
    end do!i
    do i = 1, size(this%pshape(1,:))
      do j = 1, this%pshape(1,i) * this%pshape(2,i)
        do l = 1, size(this%hole(:,1))
          if (j == this%hole(l,1)) goto 100
        end do!l
        dix = 0.40d0 * cos(xi(j))
        diy = 0.40d0 * sin(xi(j))
        write(60, '(6(f18.15, 1x))') dble(mod((j-1),this%pshape(1,i))+1)-dix, dble((j-1)/this%pshape(1,i)+1)-diy, dble(i),2.0d0*dix,2.0d0*diy, 0d0
        100 continue
      end do!j 
    end do!i
    write(60,'(2(a12,/))')' exit       ',"sp '-' w vec"
    
    ne = 0.d0
    do r = 1, 4*this%n
      write(40,*)r,eg(r)
      if(eg(r)  > 0.d0)then
        do i = 1, this%n
        if(( eg(r)) < 50*kbt) ne = ne  + (conjg(wf(4*i-3, r))*wf(4*i-3, r) + conjg(wf(4*i-2, r))*wf(4*i-2, r))/(exp(( eg(r))/kbT)+1.d0)
        if((-eg(r)) < 50*kbt) ne = ne  + (conjg(wf(4*i-1, r))*wf(4*i-1, r) + conjg(wf(4*i  , r))*wf(4*i  , r))/(exp((-eg(r))/kbT)+1.d0)
        end do!i
      end if
    end do!r
    write(40,'(2(a10,/))')' exit     ',"p '-' w l"
    sec = real(time_e - time_s)/CountPerSec
  !  min = int(sec/60)
    write(20,'(a7,f13.5,a3,2(a10,i7),2(a11,f13.10))') 'time = ',sec ,'(s) ','    iter = ',iter,' iter_one = ',iter_one,'mu = ',mu,'ne = ',ne
    write(10,*)mu,ne,iter
  end subroutine calc_sb


! SUBROUTINE car_parrinello_calc(this)
!  CLASS(car_parrinello), INTENT(inout)      :: this
!  REAL(8), DIMENSION(this%n)                :: density, Sx, Sy, Sz, oSx,oSy,oSz
!  REAL(8)      :: c, ocoef, h2, h2t, h2u,h2jd
!  REAL(8),DIMENSION(sum(int(this%ne))*2)    :: oeg
!  INTEGER      :: iter,r,i,j,k,l,f
!  TYPE(hop_iterator),ALLOCATABLE            :: hop
!  LOGICAL         :: flag
!  OPEN(10,file='cpeg.dat',status='replace')
!  flag=.FALSE.
!
!  ALLOCATE(hop, source = new_hop_iterator(this%pshape))
!  CALL hop%set_path_hole(this%hole(:,1)) !prepare <i,j> for calc Jd
!  !test
!  !    do i=1, size(hop%pl%value(:))
!  !    write(*, *) 'car_parrinello_calc',i, hop%pl%value(i)
!  !    enddo
!  !    stop
!  CALL this%set_eg()
!  this%owf = this%wf
!  this%nwf = this%wf
!  c = 1.d0/(1+this%h*this%eta*0.5d0)
!  ocoef = this%h*this%eta*0.5d0-1.d0
!  h2 = this%h**2
!  h2t = h2*this%t(1)
!  h2u = h2*this%u
!  h2jd = this%jd*h2
!  Sx = 0.d0
!  Sy = 0.d0
!  Sz = 0.d0
!  DO iter = 1, this%itermax
!     WRITE(*,'(a6,i5,a2,i6,a2)')'-iter',iter,' /',this%itermax,' -'
!     oSx = Sx
!     oSy = Sy
!     oSz = Sz
!     oeg = this%eg
!     CALL this%calc_S(density, Sx, Sy, Sz)
!     DO r=1,sum(int(this%ne))
!        DO i=1,sum(int(this%ne))
!           k=this%ets(i)
!           this%nwf(2*k-1,r) =(2.d0+h2*this%eg(r)-h2u*(density(k)*0.5d0-Sz(k)))*this%wf(2*k-1,r)&
!                &+h2u*(Sx(k)-ui*Sy(k))         *this%wf(2*k,r)&
!                &+ocoef                        *this%owf(2*k-1,r)
!           this%nwf(2*k,r)   =(2.d0+h2*this%eg(r)-h2u*(density(k)*0.5d0+Sz(k)))*this%wf(2*k,r)&
!                &+h2u*(Sx(k)+ui*Sy(k))         *this%wf(2*k-1,r)&
!                &+ocoef                        *this%owf(2*k,r)
!
!           DO j=1,SIZE(this%ns(i)%s)
!              this%nwf(2*k-1,r) =this%nwf(2*k-1,r)&
!                   &+h2t*this%wf(2*this%ns(i)%s(j)-1,r)
!              this%nwf(2*k,r)   =this%nwf(2*k,r)&
!                   &+h2t*this%wf(2*this%ns(i)%s(j),r)
!           ENDDO
!        ENDDO
!
!        DO k=1,SIZE(hop%pl%value(:))
!           i=hop%pl%value(k)%i
!           f=hop%pl%value(k)%f
!           this%nwf(2*f-1,r) = this%nwf(2*f-1,r)-h2jd*0.5d0*(Sx(i)-ui*Sy(i))*this%wf(2*f,r)
!           this%nwf(2*f  ,r) = this%nwf(2*f,  r)-h2jd*0.5d0*(Sx(i)+ui*Sy(i))*this%wf(2*f-1,r)
!           this%nwf(2*i-1,r) = this%nwf(2*i-1,r)-h2jd*0.5d0*(Sx(f)-ui*Sy(f))*this%wf(2*i,r)
!           this%nwf(2*i  ,r) = this%nwf(2*i,  r)-h2jd*0.5d0*(Sx(f)+ui*Sy(f))*this%wf(2*i-1,r)
!        ENDDO
!
!     ENDDO
!     this%owf=this%wf
!     this%wf=c*this%nwf
!     DO l=1,this%nh
!        this%wf(2*this%hole(l,1)-1,:)=0.d0
!        this%wf(2*this%hole(l,1)  ,:)=0.d0
!     ENDDO
!     CALL this%gs()
!     CALL this%set_eg()
!     WRITE(*,*)'sum energy',SUM(this%eg(1:sum(int(this%ne))))
!     WRITE(10,*)iter,SUM(this%eg(1:sum(int(this%ne))))
!     this%density=density
!     this%Sx=Sx
!     this%Sy=Sy
!     this%Sz=Sz
!     WRITE(*,*)'density',SUM(ABS(density(:))),'Sx',SUM(ABS(Sx(:)))
!     WRITE(*,*)'Sy',SUM(ABS(Sy(:))),'Sz',SUM(ABS(Sz(:)))
!     IF( this%is_converged(oSx, oSy, oSz, oeg) ) THEN
!        !     if(abs(sum(this%eg(1:sum(int(this%ne))))-sum(oeg(1:sum(int(this%ne)))))<1.d-4)then
!        flag=.TRUE.
!        EXIT
!     ENDIF
!  ENDDO
!  this%xi=this%calc_xi()
!  IF(flag)THEN
!     WRITE(*,*)'CP_caluculation has been converged!'
!     WRITE(*,*)'density',SUM(ABS(density(:))),'Sx',SUM(ABS(Sx(:)))
!     WRITE(*,*)'Sy',SUM(ABS(Sy(:))),'Sz',SUM(ABS(Sz(:)))
!     WRITE(*,*)'sum energy',SUM(this%eg(1:sum(int(this%ne))))
!     WRITE(*,*)'h2jd', h2jd
!  ELSE
!     WRITE(6, '(a30)') '-------- ! warning ! ---------'
!     WRITE(*,*)'CP_caluculation is not converged.'
!  ENDIF
!  CLOSE(10)
! END SUBROUTINE car_parrinello_calc

  SUBROUTINE car_parrinello_random_wf(this)
    CLASS(car_parrinello), INTENT(inout)        :: this
    COMPLEX(8)                         ::r(4*this%n, 4*this%n)
    REAL(8)                         ::c(4*this%n, 4*this%n), d(4*this%n, 4*this%n)
    CALL random_seed
    CALL RANDOM_NUMBER(c)
    CALL RANDOM_NUMBER(d)
    r = CMPLX(c, d)
    this%wf = r
  END SUBROUTINE car_parrinello_random_wf

  SUBROUTINE car_parrinello_gs(this) !gram schmitt
    !Gram-Schmidt orthonormalizaton
    CLASS(car_parrinello), INTENT(inout)        :: this
    COMPLEX(8)                     :: temp(4*this%n, 4*this%n)
    INTEGER k, l
    temp(:, 1) = this%normal_vec(this%wf(:, 1))
    DO k = 2, 4*this%n
      temp(:, k) = this%wf(:, k)
      DO l = 1, k - 1
        temp(:, k) = temp(:, k) - DOT_PRODUCT(temp(:, l), this%wf(:, k))*temp(:, l)
      ENDDO
      temp(:, k) = this%normal_vec(temp(:, k))
    ENDDO
    this%wf = temp
    !do k=1,sum(int(this%ne))
    !  do l = 1, sum(int(this%ne))
    !    print *, k, l, DOT_PRODUCT(this%wf(:,l),this%wf(:,k))
    !  end do
    !end do
    !stop
  END SUBROUTINE car_parrinello_gs

  FUNCTION car_parrinello_normal_vec(this, z) RESULT(r)
    CLASS(car_parrinello), INTENT(in)        :: this
    COMPLEX(8)                                  :: z(this%n*4)
    COMPLEX(8)                                  :: r(this%n*4)
    REAL(8)                                     :: norm
    norm = SQRT(DOT_PRODUCT(z(:), z(:)))
    IF (norm == 0.d0) THEN
      r = 0.d0
    ELSE
      r = z/norm
    ENDIF
  END FUNCTION car_parrinello_normal_vec

  FUNCTION car_parrinello_calc_xi(this) RESULT(r)
    CLASS(car_parrinello), INTENT(in)       :: this
    REAL(8), DIMENSION(this%n)       :: r
    INTEGER                          :: j
    r(:) = this%Sx(:)
    DO j = 1, this%n
      IF (r(j) == 0.d0) r(j) = 1.d-10
      r(j) = -0.5d0*(SIGN(1.0d0, r(j)) - 1.0d0)*pi + ATAN(this%Sy(j)/r(j))
    END DO
  END FUNCTION car_parrinello_calc_xi

  FUNCTION car_parrinello_get_wave_functions(this) RESULT(r)
    CLASS(car_parrinello), INTENT(in)           :: this
    COMPLEX(8), DIMENSION(4*this%n, 4*this%n) :: r
    r(:, :) = (0.d0, 0.d0)
    IF (ALLOCATED(this%wf)) r(:, :) = this%wf(:, :)
  END FUNCTION car_parrinello_get_wave_functions

  FUNCTION car_parrinello_get_energies(this) RESULT(r)
    CLASS(car_parrinello), INTENT(in) :: this
    REAL(8), DIMENSION(4*this%n)    :: r
    !REAL(8), DIMENSION(sum(int(this%ne)))    :: r
    integer :: i
    r(:) = 0.d0
    IF (ALLOCATED(this%eg)) then
      do i = 1, 4*this%n
      !do i = 1, sum(int(this%ne))
        r(i) = this%eg(i)
      enddo
    endif
    !print *, this%eg
  END FUNCTION car_parrinello_get_energies

  LOGICAL FUNCTION car_parrinello_is_converged(this, oSx, oSy, oSz, oeg) RESULT(r)
    CLASS(car_parrinello), INTENT(inout)             :: this
    LOGICAL, DIMENSION(this%n)             :: temp1, temp2, temp3
    LOGICAL, DIMENSION(4*this%n)             :: temp4
    REAL(8), DIMENSION(this%n), INTENT(in)  :: oSx, oSy, oSz
    REAL(8), DIMENSION(sum(int(this%ne))*4), INTENT(in) :: oeg
    REAL(8), DIMENSION(sum(int(this%ne))*4) :: seg
    !  real(8), intent(in)                    :: tol
    INTEGER                                :: i, j
    temp1(:) = .FALSE.
    temp2(:) = .FALSE.
    temp3(:) = .FALSE.
    temp4(:) = .FALSE.
    !SELECT CASE (this%nh)
    !CASE (1:)
!    do i = 1, this%n
!      if(oSx(i) == 0 .or. ABS((this%Sx(i) - oSx(i))/oSx(i)) < this%tol)  temp1(i) = .true.
!      if(oSy(i) == 0 .or. ABS((this%Sy(i) - oSy(i))/oSy(i)) < this%tol)  temp2(i) = .true.
!      if(oSz(i) == 0 .or. ABS((this%Sz(i) - oSz(i))/oSz(i)) < this%tol)  temp3(i) = .true.
!    end do
    seg(:)=0
    do i = 1,4*this%n
      if(abs(oeg(i))  <  this%tol .or. ABS((this%eg(i) - oeg(i))/oeg(i)) < this%tol)  temp4(i) = .true.
      if(oeg(i) /= 0) seg(i)=abs(oeg(i))
    end do
    write(*,*) sum(abs(seg))
    r = ALL(temp4)

    !CASE (0)
    !  IF (SUM(this%eg(1:sum(int(this%ne)))) - SUM(oeg(1:sum(int(this%ne)))) > this%tol) this%temp = .TRUE.
    !  !    if(this%temp .and. (sum(this%eg(1:sum(int(this%ne))))-sum(oeg(1:sum(int(this%ne)))))/sum(oeg(1:sum(int(this%ne)))) < 1.d-6) temp1=.true.
    !  WHERE (ABS((this%Sx - oSx)/oSx) < this%tol .OR. ABS(this%Sx) < this%tol) temp1 = .TRUE.
    !  FORALL (j=1:this%nh)
    !  temp1(this%hole(j, 1)) = .TRUE.
    !  END FORALL
    !  WHERE (ABS((this%Sy - oSy)/oSy) < this%tol .OR. ABS(this%Sy) < this%tol) temp2 = .TRUE.
    !  FORALL (j=1:this%nh)
    !  temp2(this%hole(j, 1)) = .TRUE.
    !  END FORALL
    !  WHERE (ABS((this%Sz - oSz)/oSz) < this%tol .OR. ABS(this%Sz) < this%tol) temp3 = .TRUE.
    !  FORALL (j=1:this%nh)
    !  temp3(this%hole(j, 1)) = .TRUE.
    !  END FORALL
    !  r = ALL([this%temp, temp1, temp2, temp3])
    !END SELECT
  END FUNCTION car_parrinello_is_converged

  FUNCTION new_car_parrinello(pshape, hole, ne, t, u, jd, lm, xi, chi, wf) RESULT(r)
    TYPE(car_parrinello)               :: r
    INTEGER, DIMENSION(:,:), INTENT(in)    :: pshape
    INTEGER, DIMENSION(:, :), INTENT(in) :: hole
    REAL(8), DIMENSION(:), INTENT(in)    :: ne
    REAL(8), DIMENSION(:), INTENT(in)    :: t
    REAL(8), INTENT(in)                  :: u, jd, lm
    REAL(8), DIMENSION(:), INTENT(in)    :: xi, chi
    COMPLEX(8), DIMENSION(:,:), INTENT(in), optional   :: wf
    CALL r%car_parrinello_init(pshape, hole, ne, t, u, jd, lm, xi, chi, wf)
  END FUNCTION new_car_parrinello


  ! public subroutines & functions
  !subroutine calc_surface_eg(sys, eg, wf, delta, mu, ne_u, ne_d, calc_n_max)
  subroutine calc_surface_eg(sys, eg, wf, delta, mu, ne_u, ne_d)
    class(base), intent(in):: sys
    real(8), dimension(:), intent(out) :: eg
    real(8), dimension(:), intent(in) :: ne_u, ne_d
    real(8), intent(in) :: mu
    real(8), dimension(size(eg)) :: tmp_eg
    complex(8), dimension(:, :), intent(in) :: wf
    complex(8), dimension(:), intent(in) :: delta
    !integer, intent(in), optional :: calc_n_max
    integer :: i, j, k, a, n_eg

    n_eg = size(wf(1, :))
    eg = 0.d0
    
    tmp_eg = 0d0
    do i = 1, sys%pl_1st%length()
      k = sys%pl_1st%value(i)%i
      j = sys%pl_1st%value(i)%f
      do a = 1, n_eg
          !tmp_eg(a) = tmp_eg(a) + (conjg(wf(4*k-3,a))*wf(4*j-3,a) + conjg(wf(4*j-3,a))*wf(4*k-3,a)) &                           
          !&                     + (conjg(wf(4*k-2,a))*wf(4*j-2,a) + conjg(wf(4*j-2,a))*wf(4*k-2,a)) &
          !&                     - (conjg(wf(4*k-1,a))*wf(4*j-1,a) + conjg(wf(4*j-1,a))*wf(4*k-1,a)) &
          !&                     - (conjg(wf(4*k  ,a))*wf(4*j  ,a) + conjg(wf(4*j  ,a))*wf(4*k  ,a)) 
          tmp_eg(a) = tmp_eg(a) + (conjg(wf(4*k-3,a))*wf(4*j-3,a) + conjg(wf(4*j-3,a))*wf(4*k-3,a))*(1d0-ne_d(k))*(1d0-ne_d(j)) &                           
          &                     + (conjg(wf(4*k-2,a))*wf(4*j-2,a) + conjg(wf(4*j-2,a))*wf(4*k-2,a))*(1d0-ne_u(k))*(1d0-ne_u(j)) &
          &                     - (conjg(wf(4*k-1,a))*wf(4*j-1,a) + conjg(wf(4*j-1,a))*wf(4*k-1,a))*(1d0-ne_d(k))*(1d0-ne_d(j)) &
          &                     - (conjg(wf(4*k  ,a))*wf(4*j  ,a) + conjg(wf(4*j  ,a))*wf(4*k  ,a))*(1d0-ne_u(k))*(1d0-ne_u(j))
          !tmp_eg(a) = tmp_eg(a) + this%t(1)*( (conjg(wf(4*k-3,a))*wf(4*j-3,a) + conjg(wf(4*j-3,a))*wf(4*k-3,a))*(1.d0 - ne_d(k))*(1.d0 - ne_d(j)) &                             
          !&                                        + (conjg(wf(4*k-2,a))*wf(4*j-2,a) + conjg(wf(4*j-2,a))*wf(4*k-2,a))*(1.d0 - ne_u(k))*(1.d0 - ne_u(j)) &
          !&                                        - (conjg(wf(4*k-1,a))*wf(4*j-1,a) + conjg(wf(4*j-1,a))*wf(4*k-1,a))*(1.d0 - ne_d(k))*(1.d0 - ne_d(j))&
          !&                                        - (conjg(wf(4*k  ,a))*wf(4*j  ,a) + conjg(wf(4*j  ,a))*wf(4*k  ,a))*(1.d0 - ne_u(k))*(1.d0 - ne_u(j)) )
        !if(dflag)  tmp_eg(a) = tmp_eg(a) + delta(i)*( conjg(wf(4*k-3,a))*wf(4*j  ,a) + conjg(wf(4*j-2,a))*wf(4*k-1,a) &
        ! tmp_eg(a) = tmp_eg(a) + delta(i)*( conjg(wf(4*k-3,a))*wf(4*j  ,a) + conjg(wf(4*j-2,a))*wf(4*k-1,a) &
        !&                                           + conjg(wf(4*j-3,a))*wf(4*k  ,a) + conjg(wf(4*k-2,a))*wf(4*j-1,a))&
        !&                         - conjg(delta(i))*( conjg(wf(4*k-1,a))*wf(4*j-2,a) + conjg(wf(4*j  ,a))*wf(4*k-3,a) &
        !&                                          +  conjg(wf(4*j-1,a))*wf(4*k-2,a) + conjg(wf(4*k  ,a))*wf(4*j-3,a))
      end do!a
    end do!i
    eg(:) = eg(:) - sys%t(1)*tmp_eg(:) 
    
    tmp_eg = 0d0
    do i = 1, sys%pl_1st%length()
      k = sys%pl_1st%value(i)%i
      j = sys%pl_1st%value(i)%f
      do a = 1, n_eg
          tmp_eg(a) = tmp_eg(a) + delta(i)*(  conjg(wf(4*k-3,a))*wf(4*j  ,a) + conjg(wf(4*j-3,a))*wf(4*k  ,a) &
        &                                   + conjg(wf(4*k-2,a))*wf(4*j-1,a) + conjg(wf(4*j-2,a))*wf(4*k-1,a))&
        &                 + conjg(delta(i))*( conjg(wf(4*k-1,a))*wf(4*j-2,a) + conjg(wf(4*j-1,a))*wf(4*k-2,a) &
        &                                   + conjg(wf(4*k  ,a))*wf(4*j-3,a) + conjg(wf(4*j  ,a))*wf(4*k-3,a))
      end do!a
    end do!i
    eg(:) = eg(:) + tmp_eg(:) 
    
    tmp_eg = 0d0
    do i = 1, sys%pl_2nd%length()
      k = sys%pl_2nd%value(i)%i
      j = sys%pl_2nd%value(i)%f
      do a = 1, n_eg
      !do a = 1, 4*this%n
        tmp_eg(a) = tmp_eg(a) + (conjg(wf(4*k-3,a))*wf(4*j-3,a) + conjg(wf(4*j-3,a))*wf(4*k-3,a))*(1d0-ne_d(k))*(1d0-ne_d(j)) &
        &                     + (conjg(wf(4*k-2,a))*wf(4*j-2,a) + conjg(wf(4*j-2,a))*wf(4*k-2,a))*(1d0-ne_u(k))*(1d0-ne_u(j)) &
        &                     - (conjg(wf(4*k-1,a))*wf(4*j-1,a) + conjg(wf(4*j-1,a))*wf(4*k-1,a))*(1d0-ne_d(k))*(1d0-ne_d(j)) &
        &                     - (conjg(wf(4*k  ,a))*wf(4*j  ,a) + conjg(wf(4*j  ,a))*wf(4*k  ,a))*(1d0-ne_u(k))*(1d0-ne_u(j))  
        !tmp_eg(a) = tmp_eg(a) + (conjg(wf(4*k-3,a))*wf(4*j-3,a) + conjg(wf(4*j-3,a))*wf(4*k-3,a))*(1.d0 - ne_d(k))*(1.d0 - ne_d(j)) &
        !&                                  + (conjg(wf(4*k-2,a))*wf(4*j-2,a) + conjg(wf(4*j-2,a))*wf(4*k-2,a))*(1.d0 - ne_u(k))*(1.d0 - ne_u(j)) &
        !&                                  - (conjg(wf(4*k-1,a))*wf(4*j-1,a) + conjg(wf(4*j-1,a))*wf(4*k-1,a))*(1.d0 - ne_d(k))*(1.d0 - ne_d(j)) &
        !&                                  - (conjg(wf(4*k  ,a))*wf(4*j  ,a) + conjg(wf(4*j  ,a))*wf(4*k  ,a))*(1.d0 - ne_u(k))*(1.d0 - ne_u(j))  
      end do!a
    end do!i
    eg(:) = eg(:) - sys%t(2)*tmp_eg(:) 
    !eg(1:4*this%n) = eg(1:4*this%n) - tmp_eg(1:4*this%n) * this%t(2) 
    
    tmp_eg = 0d0
    !do a = 1, 2*this%n 
    do a = 1, n_eg 
      do j = 1, sys%n ! n_acc
        !j = this%ets(k)
        !tmp_eg(a) = tmp_eg(a) - mu*( conjg(wf(4*j-3,a))*wf(4*j-3,a) + conjg(wf(4*j-2,a))*wf(4*j-2,a) ) &
        !&                     + mu*( conjg(wf(4*j-1,a))*wf(4*j-1,a) + conjg(wf(4*j  ,a))*wf(4*j  ,a) )
        tmp_eg(a) = tmp_eg(a)  + conjg(wf(4*j-3,a))*wf(4*j-3,a) + conjg(wf(4*j-2,a))*wf(4*j-2,a) &
        &                      - conjg(wf(4*j-1,a))*wf(4*j-1,a) - conjg(wf(4*j  ,a))*wf(4*j  ,a)
      end do!j
    end do!a
    !eg(:) = eg(:) + tmp_eg(:) 
    eg(:) = eg(:) - mu* tmp_eg(:) 
    !eg(1:2*this%n) = eg(1:2*this%n) - tmp_eg(1:2*this%n) * mu
    !eg(1:4*this%n) = eg(1:4*this%n) - tmp_eg(1:4*this%n) * mu


  end subroutine calc_surface_eg
  
  subroutine calc_surface_ne(sys, wf, eg, ne_u, ne_d, kbt)
    class(base), intent(in) :: sys
    complex(8), dimension(:, :), intent(in) :: wf
    real(8), dimension(:), intent(in) :: eg
    real(8), dimension(:), intent(out) :: ne_u, ne_d
    real(8), intent(in) :: kbt
    integer :: i, a
    ne_u = 0.d0
    ne_d = 0.d0
    do i = 1, sys%n
      do a = 2*sys%n+1, 4*sys%n
        ne_u(i) = ne_u(i) + fermi( eg(a), kbT)*conjg(wf(4*i-3, a))*wf(4*i-3, a) + fermi(-eg(a), kbT)*conjg(wf(4*i-1, a))*wf(4*i-1, a)
        ne_d(i) = ne_d(i) + fermi( eg(a), kbT)*conjg(wf(4*i-2, a))*wf(4*i-2, a) + fermi(-eg(a), kbT)*conjg(wf(4*i  , a))*wf(4*i  , a) 
      end do!a
    end do!i
  end subroutine calc_surface_ne

  subroutine calc_surface_delta(sys, wf, eg, delta, kbt)
    class(base), intent(in) :: sys
    complex(8), dimension(:, :), intent(in) :: wf
    real(8), dimension(:), intent(in) :: eg
    complex(8), dimension(:), intent(out) :: delta
    real(8), intent(in) :: kbt
    integer :: i, j, k, a
    delta  = 0.d0
    do a = 2*sys%n+1, 4*sys%n
      do k = 1, sys%pl_1st%length()
        i = sys%pl_1st%value(k)%i
        j = sys%pl_1st%value(k)%f
        ! je : 2*t(1)**2 / U
        delta(k) = delta(k) + (sys%jd*(wf(4*i-3, a)*conjg(wf(4*j  , a)) + wf(4*j-2, a)*conjg(wf(4*i-1, a)) &
        &                 + wf(4*i-2, a)*conjg(wf(4*j-1, a)) + wf(4*j-3, a)*conjg(wf(4*i  , a))))*tanh(eg(a)/kbT*0.5d0)
      end do!k
    end do!a
  end subroutine calc_surface_delta

  !subroutine calc_surface_S(sys, wf, eg, Sx, Sy, Sz, kbt)
  subroutine calc_surface_S(sys, wf, eg, Sx, Sy, Sz)
    class(base), intent(in) :: sys
    complex(8), dimension(:, :), intent(in) :: wf
    real(8), dimension(:), intent(in) :: eg
    real(8), dimension(:), intent(out) :: Sx, Sy, Sz
    !real(8), intent(in) :: kbt
    integer :: i, j, k, a
    complex(8) :: uu, ud, vu, vd
    Sx = 0d0
    Sy = 0d0
    Sz = 0d0
    do i = 1, sys%n
      do a = 2*sys%n + 1, 4*sys%n 
        uu = wf(4*i-3, a)
        ud = wf(4*i-2, a)
        vu = wf(4*i-1, a)
        vd = wf(4*i  , a)
        Sx(i) = Sx(i) + real(vu*conjg(vd))
        Sy(i) = Sy(i) - aimag(vu*conjg(vd))
        !Sx(i) = Sx(i) + real(uu*conjg(ud))
        !Sy(i) = Sy(i) - aimag(uu*conjg(ud))
        Sz(i) = Sz(i) + real(vu*conjg(vu)-vd*conjg(vd))/2
      end do
    end do
  end subroutine calc_surface_S

  function calc_surface_converged(v, ov, tol) result(r)
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
  end function calc_surface_converged
  
  subroutine diagonalize_surface_hamiltonian(sys, eg, wf, mu, delta, ne_u, ne_d)
    use math_mod, only : my_zheev
    implicit none
    class(base), intent(in)   :: sys
    REAL(8), DIMENSION(4*sys%n), intent(out) :: eg
    complex(8),DIMENSION(4*sys%n,4*sys%n), intent(out)   :: wf
    real(8),intent(in)    :: mu
    complex(8), dimension(sys%pl_1st%length()), intent(in) :: delta
    real(8), dimension(sys%n), intent(in)   :: ne_u, ne_d
    complex(8),DIMENSION(4*sys%n,4*sys%n)   :: ham
    integer :: i, j, k, r

    ham = 0d0
    do i = 1, sys%n
      ham(4*i-3, 4*i-3) = -mu
      ham(4*i-2, 4*i-2) = -mu
      ham(4*i-1, 4*i-1) =  mu
      ham(4*i  , 4*i  ) =  mu
    end do
    do i = 1, sys%pl_1st%length()
      j = sys%pl_1st%value(i)%i
      k = sys%pl_1st%value(i)%f
      ham(4*j-3, 4*k-3) = -sys%t(1)*(1d0-ne_d(j))*(1d0-ne_d(k))
      ham(4*j-2, 4*k-2) = -sys%t(1)*(1d0-ne_u(j))*(1d0-ne_u(k))
      ham(4*j-1, 4*k-1) =  sys%t(1)*(1d0-ne_d(j))*(1d0-ne_d(k))
      ham(4*j  , 4*k  ) =  sys%t(1)*(1d0-ne_u(j))*(1d0-ne_u(k))
      ham(4*k-3, 4*j-3) = -sys%t(1)*(1d0-ne_d(j))*(1d0-ne_d(k))
      ham(4*k-2, 4*j-2) = -sys%t(1)*(1d0-ne_u(j))*(1d0-ne_u(k))
      ham(4*k-1, 4*j-1) =  sys%t(1)*(1d0-ne_d(j))*(1d0-ne_d(k))
      ham(4*k  , 4*j  ) =  sys%t(1)*(1d0-ne_u(j))*(1d0-ne_u(k))
    end do
    do i = 1, sys%pl_2nd%length()
      j = sys%pl_2nd%value(i)%i
      k = sys%pl_2nd%value(i)%f
      ham(4*j-3, 4*k-3) = -sys%t(2)*(1d0-ne_d(j))*(1d0-ne_d(k))
      ham(4*j-2, 4*k-2) = -sys%t(2)*(1d0-ne_u(j))*(1d0-ne_u(k))
      ham(4*j-1, 4*k-1) =  sys%t(2)*(1d0-ne_d(j))*(1d0-ne_d(k))
      ham(4*j  , 4*k  ) =  sys%t(2)*(1d0-ne_u(j))*(1d0-ne_u(k))
      ham(4*k-3, 4*j-3) = -sys%t(2)*(1d0-ne_d(j))*(1d0-ne_d(k))
      ham(4*k-2, 4*j-2) = -sys%t(2)*(1d0-ne_u(j))*(1d0-ne_u(k))
      ham(4*k-1, 4*j-1) =  sys%t(2)*(1d0-ne_d(j))*(1d0-ne_d(k))
      ham(4*k  , 4*j  ) =  sys%t(2)*(1d0-ne_u(j))*(1d0-ne_u(k))
    end do
    do i = 1, sys%pl_1st%length()
      j = sys%pl_1st%value(i)%i
      k = sys%pl_1st%value(i)%f
      ham(4*j-3, 4*k  ) = delta(i)
      ham(4*j-2, 4*k-1) = delta(i)
      ham(4*j-1, 4*k-2) = conjg(delta(i))
      ham(4*j  , 4*k-3) = conjg(delta(i))
      ham(4*k-3, 4*j  ) = delta(i)
      ham(4*k-2, 4*j-1) = delta(i)
      ham(4*k-1, 4*j-2) = conjg(delta(i))
      ham(4*k  , 4*j-3) = conjg(delta(i))
    end do
  
    call my_zheev("l", ham, eg, wf)

  end subroutine diagonalize_surface_hamiltonian
  
  subroutine calc_bulk_ne(sys, wf, eg, ne_u, ne_d, kbt)
    class(base), intent(in) :: sys
    complex(8), dimension(4*sys%n_acc, 4*sys%n_acc), intent(in) :: wf
    real(8), dimension(4*sys%n_acc), intent(in) :: eg
    real(8), dimension(sys%n), intent(out) :: ne_u, ne_d
    real(8), intent(in) :: kbt
    integer :: i, is, a
    ne_u = 0.d0
    ne_d = 0.d0
    do i = 1, sys%n_acc
      is = sys%ets(i)
      do a = 2*sys%n_acc+1, 4*sys%n_acc
        ne_u(is) = ne_u(is) + fermi( eg(a), kbT)*conjg(wf(4*i-3, a))*wf(4*i-3, a) + fermi(-eg(a), kbT)*conjg(wf(4*i-1, a))*wf(4*i-1, a)
        ne_d(is) = ne_d(is) + fermi( eg(a), kbT)*conjg(wf(4*i-2, a))*wf(4*i-2, a) + fermi(-eg(a), kbT)*conjg(wf(4*i  , a))*wf(4*i  , a) 
      end do!a
    end do!i
  end subroutine calc_bulk_ne

  !subroutine calc_bulk_S(sys, wf, eg, Sx, Sy, Sz, kbt)
  subroutine calc_bulk_S(sys, wf, eg, Sx, Sy, Sz)
    class(base), intent(in) :: sys
    complex(8), dimension(:, :), intent(in) :: wf
    real(8), dimension(:), intent(in) :: eg
    real(8), dimension(:), intent(out) :: Sx, Sy, Sz
    !real(8), intent(in) :: kbt
    integer :: i, is, k, a
    complex(8) :: uu, ud, vu, vd
    Sx = 0d0
    Sy = 0d0
    Sz = 0d0
    do i = 1, sys%n_acc
      is = sys%ets(i)
      do a = 2*sys%n_acc + 1, 4*sys%n_acc 
        uu = wf(4*i-3, a)
        ud = wf(4*i-2, a)
        vu = wf(4*i-1, a)
        vd = wf(4*i  , a)
        !Sx(is) = Sx(is) + real(vu*conjg(vd))
        !Sy(is) = Sy(is) - aimag(vu*conjg(vd))
        !Sz(is) = Sz(is) + real(vu*conjg(vu)-vd*conjg(vd))/2
        Sx(is) = Sx(is) - 0.5d0*(vu*conjg(vd) + conjg(vu)*vd)
        Sy(is) = Sy(is) + 0.5d0*ui*(vu*conjg(vd) - conjg(vu)*vd)
        Sz(is) = Sz(is) + 0.5d0*(vu*conjg(vu)-vd*conjg(vd))
        
        !Sx(is) = Sx(is) + real(uu*conjg(ud))
        !Sy(is) = Sy(is) - aimag(uu*conjg(ud))
        !Sz(is) = Sz(is) + real(uu*conjg(uu)-ud*conjg(ud))/2
        
        !Sx(is) = Sx(is) + real(uu*conjg(ud))
        !!Sy(is) = Sy(is) - aimag(uu*conjg(ud))
        !!Sz(is) = Sz(is) + real(uu*conjg(vu)-ud*conjg(vd))/2
      end do
    end do
  end subroutine calc_bulk_S
  
  subroutine diagonalize_bulk_hamiltonian(sys, eg, wf, mu, xi, chi)
    use math_mod, only : my_zheev
    implicit none
    class(base), intent(in)   :: sys
    REAL(8), DIMENSION(4*sys%n_acc), intent(out) :: eg
    complex(8),DIMENSION(4*sys%n_acc, 4*sys%n_acc), intent(out)   :: wf
    real(8),intent(in)    :: mu
    real(8), dimension(sys%n), intent(in) :: xi, chi
    complex(8),DIMENSION(4*sys%n_acc, 4*sys%n_acc)   :: ham
    real(8), dimension(sys%n) :: Sx, Sy, Sz
    !real(8), dimension(sys%n) :: ne_u, ne_d
    integer :: i, j, k, r, is, js, ks

    do i = 1 ,sys%n
      Sx(i) = 0.5d0*cos(xi(i))
      Sy(i) = 0.5d0*sin(xi(i))
    end do
    Sz = 0d0

    !ne_u = 0.5d0
    !ne_d = 0.5d0

    ham = 0d0
    do i = 1, sys%n_acc
      is = sys%ets(i)
      ham(4*i-3, 4*i-3) = sys%U*(-2d0/3d0*Sz(is)+0.5d0)
      ham(4*i-2, 4*i-2) = sys%U*( 2d0/3d0*Sz(is)+0.5d0)
      ham(4*i-3, 4*i-2) = - sys%U*2d0/3d0*(Sx(is)-ui*Sy(is))
      ham(4*i-2, 4*i-3) = - sys%U*2d0/3d0*(Sx(is)+ui*Sy(is))
      ham(4*i-1, 4*i-1) = - sys%U*(-2d0/3d0*Sz(is)+0.5d0) ! - conjg
      ham(4*i  , 4*i  ) = - sys%U*( 2d0/3d0*Sz(is)+0.5d0) ! - conjg
      ham(4*i-1, 4*i  ) = conjg(- sys%U*2d0/3d0*(Sx(is)-ui*Sy(is))) ! + conjg
      ham(4*i  , 4*i-1) = conjg(- sys%U*2d0/3d0*(Sx(is)+ui*Sy(is))) ! + conjg
    end do
    do i = 1, sys%n_acc
      ham(4*i-3, 4*i-3) = ham(4*i-3, 4*i-3) - mu
      ham(4*i-2, 4*i-2) = ham(4*i-2, 4*i-2) - mu
      ham(4*i-1, 4*i-1) = ham(4*i-1, 4*i-1) + mu
      ham(4*i  , 4*i  ) = ham(4*i  , 4*i  ) + mu
    end do
    do i = 1, sys%pl_1st%length()
      js = sys%pl_1st%value(i)%i
      ks = sys%pl_1st%value(i)%f
      j = sys%ste(js)
      k = sys%ste(ks)
      ham(4*j-3, 4*k-3) = -sys%t(1)
      ham(4*j-2, 4*k-2) = -sys%t(1)
      ham(4*j-1, 4*k-1) =  sys%t(1)
      ham(4*j  , 4*k  ) =  sys%t(1)
      ham(4*k-3, 4*j-3) = -sys%t(1)
      ham(4*k-2, 4*j-2) = -sys%t(1)
      ham(4*k-1, 4*j-1) =  sys%t(1)
      ham(4*k  , 4*j  ) =  sys%t(1)
    end do
    do i = 1, sys%pl_2nd%length()
      js = sys%pl_2nd%value(i)%i
      ks = sys%pl_2nd%value(i)%f
      j = sys%ste(js)
      k = sys%ste(ks)
      ham(4*j-3, 4*k-3) = -sys%t(2)
      ham(4*j-2, 4*k-2) = -sys%t(2)
      ham(4*j-1, 4*k-1) =  sys%t(2)
      ham(4*j  , 4*k  ) =  sys%t(2)
      ham(4*k-3, 4*j-3) = -sys%t(2)
      ham(4*k-2, 4*j-2) = -sys%t(2)
      ham(4*k-1, 4*j-1) =  sys%t(2)
      ham(4*k  , 4*j  ) =  sys%t(2)
    end do
    do i = 1, sys%pl_h%length()
      js = sys%pl_h%value(i)%i
      ks = sys%pl_h%value(i)%f
      j = sys%ste(js)
      k = sys%ste(ks)
      ham(4*j-3, 4*j-3) = ham(4*j-3, 4*j-3) + 0.5d0*sys%jd*Sz(ks)
      ham(4*j-2, 4*j-2) = ham(4*j-2, 4*j-2) - 0.5d0*sys%jd*Sz(ks)
      ham(4*j-3, 4*j-2) = ham(4*j-3, 4*j-2) + 0.5d0*sys%jd*(Sx(ks) - ui*Sy(ks))
      ham(4*j-2, 4*j-3) = ham(4*j-2, 4*j-3) + 0.5d0*sys%jd*(Sx(ks) + ui*Sy(ks))
      ham(4*j-1, 4*j-1) = ham(4*j-1, 4*j-1) - 0.5d0*sys%jd*Sz(ks) ! -conjg
      ham(4*j  , 4*j  ) = ham(4*j  , 4*j  ) + 0.5d0*sys%jd*Sz(ks) ! -conjg
      ham(4*j-1, 4*j  ) = ham(4*j-1, 4*j  ) + conjg( 0.5d0*sys%jd*(Sx(ks) - ui*Sy(ks)) ) ! + conjg
      ham(4*j  , 4*j-1) = ham(4*j  , 4*j-1) + conjg( 0.5d0*sys%jd*(Sx(ks) + ui*Sy(ks)) ) ! + conjg
      
      ham(4*k-3, 4*k-3) = ham(4*k-3, 4*k-3) + 0.5d0*sys%jd*Sz(js)
      ham(4*k-2, 4*k-2) = ham(4*k-2, 4*k-2) - 0.5d0*sys%jd*Sz(js)
      ham(4*k-3, 4*k-2) = ham(4*k-3, 4*k-2) + 0.5d0*sys%jd*(Sx(js) - ui*Sy(js))
      ham(4*k-2, 4*k-3) = ham(4*k-2, 4*k-3) + 0.5d0*sys%jd*(Sx(js) + ui*Sy(js))
      ham(4*k-1, 4*k-1) = ham(4*k-1, 4*k-1) - 0.5d0*sys%jd*Sz(js) ! -conjg
      ham(4*k  , 4*k  ) = ham(4*k  , 4*k  ) + 0.5d0*sys%jd*Sz(js) ! -conjg
      ham(4*k-1, 4*k  ) = ham(4*k-1, 4*k  ) + conjg( 0.5d0*sys%jd*(Sx(js) - ui*Sy(js)) ) ! + conjg
      ham(4*k  , 4*k-1) = ham(4*k  , 4*k-1) + conjg( 0.5d0*sys%jd*(Sx(js) + ui*Sy(js)) ) ! + conjg
    end do
    !do i = 1, sys%pl_1st%length()
    !  j = sys%pl_1st%value(i)%i
    !  k = sys%pl_1st%value(i)%f
    !  ham(4*j-3, 4*k  ) = delta(i)
    !  ham(4*j-2, 4*k-1) = delta(i)
    !  ham(4*j-1, 4*k-2) = conjg(delta(i))
    !  ham(4*j  , 4*k-3) = conjg(delta(i))
    !  ham(4*k-3, 4*j  ) = delta(i)
    !  ham(4*k-2, 4*j-1) = delta(i)
    !  ham(4*k-1, 4*j-2) = conjg(delta(i))
    !  ham(4*k  , 4*j-3) = conjg(delta(i))
    !end do
  
    call my_zheev("l", ham, eg, wf)

  end subroutine diagonalize_bulk_hamiltonian

  subroutine calc_sb_ne(sys, wf, eg, ne_u, ne_d, kbt)
    class(base), intent(in) :: sys
    complex(8), dimension(4*sys%n_acc, 4*sys%n_acc), intent(in) :: wf
    real(8), dimension(4*sys%n_acc), intent(in) :: eg
    real(8), dimension(sys%n), intent(out) :: ne_u, ne_d
    real(8), intent(in) :: kbt
    integer :: i, is, a
    ne_u = 0.d0
    ne_d = 0.d0
    do i = 1, sys%n_acc
      is = sys%ets(i)
      do a = 2*sys%n_acc+1, 4*sys%n_acc
        ne_u(is) = ne_u(is) + fermi( eg(a), kbT)*conjg(wf(4*i-3, a))*wf(4*i-3, a) + fermi(-eg(a), kbT)*conjg(wf(4*i-1, a))*wf(4*i-1, a)
        ne_d(is) = ne_d(is) + fermi( eg(a), kbT)*conjg(wf(4*i-2, a))*wf(4*i-2, a) + fermi(-eg(a), kbT)*conjg(wf(4*i  , a))*wf(4*i  , a) 
      end do!a
    end do!i
  end subroutine calc_sb_ne

  subroutine diagonalize_sb_2lay_hamiltonian(sys, eg, wf, mu, delta, ne_u, ne_d, Sx, Sy, Sz)
    ! z =1 : surface
    ! z =2 : bulk
    use math_mod, only : my_zheev
    implicit none
    class(base), intent(in)   :: sys
    REAL(8), DIMENSION(4*sys%n_acc), intent(out) :: eg
    complex(8),DIMENSION(4*sys%n_acc, 4*sys%n_acc), intent(out)   :: wf
    real(8),intent(in)    :: mu(2)
    complex(8), dimension(sys%pl_1st_l(1)%length()), intent(in) :: delta
    real(8), dimension(sys%n), intent(in)   :: ne_u, ne_d
    !real(8), dimension(sys%n), intent(in) :: xi, chi
    real(8), dimension(sys%n), intent(in) :: Sx, Sy, Sz
    complex(8),DIMENSION(4*sys%n_acc, 4*sys%n_acc)   :: ham
    !real(8), dimension(sys%n) :: ne_u, ne_d
    integer :: i, j, k, r, is, js, ks
    real(8) :: t3U

    if(sys%nz /= 2) then
      print *, "nz : ", sys%nz
      stop "error : wrong data in diagonalize_sb_2lay_hamiltonian"
      ! this subroutine is available for only 2 layer"
    endif

    t3U = 2*sys%t(3)**2/sys%U

    !ne_u = 0.5d0
    !ne_d = 0.5d0

    ham = 0d0
    
    !!!!!!!!!!!!!!!!!!    surface    !!!!!!!!!!!!!!!!!!
    do i = sys%n_acc_start(1), sys%n_acc_end(1)
      ham(4*i-3, 4*i-3) = -mu(1)
      ham(4*i-2, 4*i-2) = -mu(1)
      ham(4*i-1, 4*i-1) =  mu(1)
      ham(4*i  , 4*i  ) =  mu(1)
    end do
    do i = 1, sys%pl_1st_l(1)%length()
      js = sys%pl_1st_l(1)%value(i)%i
      ks = sys%pl_1st_l(1)%value(i)%f
      j = sys%ste(js) ! in the 1st layer ( surface ), js = j, ks = k
      k = sys%ste(ks)
      ham(4*j-3, 4*k-3) = -sys%t(1)*(1d0-ne_d(js))*(1d0-ne_d(ks))
      ham(4*j-2, 4*k-2) = -sys%t(1)*(1d0-ne_u(js))*(1d0-ne_u(ks))
      ham(4*j-1, 4*k-1) =  sys%t(1)*(1d0-ne_d(js))*(1d0-ne_d(ks))
      ham(4*j  , 4*k  ) =  sys%t(1)*(1d0-ne_u(js))*(1d0-ne_u(ks))
      ham(4*k-3, 4*j-3) = -sys%t(1)*(1d0-ne_d(js))*(1d0-ne_d(ks))
      ham(4*k-2, 4*j-2) = -sys%t(1)*(1d0-ne_u(js))*(1d0-ne_u(ks))
      ham(4*k-1, 4*j-1) =  sys%t(1)*(1d0-ne_d(js))*(1d0-ne_d(ks))
      ham(4*k  , 4*j  ) =  sys%t(1)*(1d0-ne_u(js))*(1d0-ne_u(ks))
    end do
    do i = 1, sys%pl_2nd_l(1)%length()
      js = sys%pl_2nd_l(1)%value(i)%i
      ks = sys%pl_2nd_l(1)%value(i)%f
      j = sys%ste(js) ! in the 1st layer ( surface ), js = j, ks = k
      k = sys%ste(ks)
      ham(4*j-3, 4*k-3) = -sys%t(2)*(1d0-ne_d(js))*(1d0-ne_d(ks))
      ham(4*j-2, 4*k-2) = -sys%t(2)*(1d0-ne_u(js))*(1d0-ne_u(ks))
      ham(4*j-1, 4*k-1) =  sys%t(2)*(1d0-ne_d(js))*(1d0-ne_d(ks))
      ham(4*j  , 4*k  ) =  sys%t(2)*(1d0-ne_u(js))*(1d0-ne_u(ks))
      ham(4*k-3, 4*j-3) = -sys%t(2)*(1d0-ne_d(js))*(1d0-ne_d(ks))
      ham(4*k-2, 4*j-2) = -sys%t(2)*(1d0-ne_u(js))*(1d0-ne_u(ks))
      ham(4*k-1, 4*j-1) =  sys%t(2)*(1d0-ne_d(js))*(1d0-ne_d(ks))
      ham(4*k  , 4*j  ) =  sys%t(2)*(1d0-ne_u(js))*(1d0-ne_u(ks))
    end do
    do i = 1, sys%pl_1st_l(1)%length()
      ! here, delta's index corresponds to pl_1st_l(1)'s index
      js = sys%pl_1st_l(1)%value(i)%i
      ks = sys%pl_1st_l(1)%value(i)%f
      j = sys%ste(js) ! in the 1st layer ( surface ), js = j, ks = k
      k = sys%ste(ks)
      ham(4*j-3, 4*k  ) = delta(i)
      ham(4*j-2, 4*k-1) = delta(i)
      ham(4*j-1, 4*k-2) = conjg(delta(i))
      ham(4*j  , 4*k-3) = conjg(delta(i))
      ham(4*k-3, 4*j  ) = delta(i)
      ham(4*k-2, 4*j-1) = delta(i)
      ham(4*k-1, 4*j-2) = conjg(delta(i))
      ham(4*k  , 4*j-3) = conjg(delta(i))
    end do
    !!!!!!!!!!!!!!    end surface    !!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!    bulk    !!!!!!!!!!!!!!!!
    !do i = 1, sys%n_acc
    do i = sys%n_acc_start(2), sys%n_acc_end(2)
      is = sys%ets(i)
      ham(4*i-3, 4*i-3) = sys%U*(-2d0/3d0*Sz(is)+0.5d0)
      ham(4*i-2, 4*i-2) = sys%U*( 2d0/3d0*Sz(is)+0.5d0)
      ham(4*i-3, 4*i-2) = - sys%U*2d0/3d0*(Sx(is)-ui*Sy(is))
      ham(4*i-2, 4*i-3) = - sys%U*2d0/3d0*(Sx(is)+ui*Sy(is))
      ham(4*i-1, 4*i-1) = - sys%U*(-2d0/3d0*Sz(is)+0.5d0) ! - conjg
      ham(4*i  , 4*i  ) = - sys%U*( 2d0/3d0*Sz(is)+0.5d0) ! - conjg
      ham(4*i-1, 4*i  ) = conjg(- sys%U*2d0/3d0*(Sx(is)-ui*Sy(is))) ! + conjg
      ham(4*i  , 4*i-1) = conjg(- sys%U*2d0/3d0*(Sx(is)+ui*Sy(is))) ! + conjg
    end do
    do i = sys%n_acc_start(2), sys%n_acc_end(2)
      ham(4*i-3, 4*i-3) = ham(4*i-3, 4*i-3) - mu(2)
      ham(4*i-2, 4*i-2) = ham(4*i-2, 4*i-2) - mu(2)
      ham(4*i-1, 4*i-1) = ham(4*i-1, 4*i-1) + mu(2)
      ham(4*i  , 4*i  ) = ham(4*i  , 4*i  ) + mu(2)
    end do
    do i = 1, sys%pl_1st_l(2)%length()
      js = sys%pl_1st_l(2)%value(i)%i
      ks = sys%pl_1st_l(2)%value(i)%f
      j = sys%ste(js)
      k = sys%ste(ks)
      ham(4*j-3, 4*k-3) = -sys%t(1)
      ham(4*j-2, 4*k-2) = -sys%t(1)
      ham(4*j-1, 4*k-1) =  sys%t(1)
      ham(4*j  , 4*k  ) =  sys%t(1)
      ham(4*k-3, 4*j-3) = -sys%t(1)
      ham(4*k-2, 4*j-2) = -sys%t(1)
      ham(4*k-1, 4*j-1) =  sys%t(1)
      ham(4*k  , 4*j  ) =  sys%t(1)
    end do
    do i = 1, sys%pl_2nd_l(2)%length()
      js = sys%pl_2nd_l(2)%value(i)%i
      ks = sys%pl_2nd_l(2)%value(i)%f
      j = sys%ste(js)
      k = sys%ste(ks)
      ham(4*j-3, 4*k-3) = -sys%t(2)
      ham(4*j-2, 4*k-2) = -sys%t(2)
      ham(4*j-1, 4*k-1) =  sys%t(2)
      ham(4*j  , 4*k  ) =  sys%t(2)
      ham(4*k-3, 4*j-3) = -sys%t(2)
      ham(4*k-2, 4*j-2) = -sys%t(2)
      ham(4*k-1, 4*j-1) =  sys%t(2)
      ham(4*k  , 4*j  ) =  sys%t(2)
    end do
    do i = 1, sys%pl_h%length()
      js = sys%pl_h%value(i)%i
      ks = sys%pl_h%value(i)%f
      j = sys%ste(js)
      k = sys%ste(ks)
      ham(4*j-3, 4*j-3) = ham(4*j-3, 4*j-3) + 0.5d0*sys%jd*Sz(ks)
      ham(4*j-2, 4*j-2) = ham(4*j-2, 4*j-2) - 0.5d0*sys%jd*Sz(ks)
      ham(4*j-3, 4*j-2) = ham(4*j-3, 4*j-2) + 0.5d0*sys%jd*(Sx(ks) - ui*Sy(ks))
      ham(4*j-2, 4*j-3) = ham(4*j-2, 4*j-3) + 0.5d0*sys%jd*(Sx(ks) + ui*Sy(ks))
      ham(4*j-1, 4*j-1) = ham(4*j-1, 4*j-1) - 0.5d0*sys%jd*Sz(ks) ! -conjg
      ham(4*j  , 4*j  ) = ham(4*j  , 4*j  ) + 0.5d0*sys%jd*Sz(ks) ! -conjg
      ham(4*j-1, 4*j  ) = ham(4*j-1, 4*j  ) + conjg( 0.5d0*sys%jd*(Sx(ks) - ui*Sy(ks)) ) ! + conjg
      ham(4*j  , 4*j-1) = ham(4*j  , 4*j-1) + conjg( 0.5d0*sys%jd*(Sx(ks) + ui*Sy(ks)) ) ! + conjg
      
      ham(4*k-3, 4*k-3) = ham(4*k-3, 4*k-3) + 0.5d0*sys%jd*Sz(js)
      ham(4*k-2, 4*k-2) = ham(4*k-2, 4*k-2) - 0.5d0*sys%jd*Sz(js)
      ham(4*k-3, 4*k-2) = ham(4*k-3, 4*k-2) + 0.5d0*sys%jd*(Sx(js) - ui*Sy(js))
      ham(4*k-2, 4*k-3) = ham(4*k-2, 4*k-3) + 0.5d0*sys%jd*(Sx(js) + ui*Sy(js))
      ham(4*k-1, 4*k-1) = ham(4*k-1, 4*k-1) - 0.5d0*sys%jd*Sz(js) ! -conjg
      ham(4*k  , 4*k  ) = ham(4*k  , 4*k  ) + 0.5d0*sys%jd*Sz(js) ! -conjg
      ham(4*k-1, 4*k  ) = ham(4*k-1, 4*k  ) + conjg( 0.5d0*sys%jd*(Sx(js) - ui*Sy(js)) ) ! + conjg
      ham(4*k  , 4*k-1) = ham(4*k  , 4*k-1) + conjg( 0.5d0*sys%jd*(Sx(js) + ui*Sy(js)) ) ! + conjg
    end do
    
    !!!!!!!!!!!!!!    end bulk    !!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!    interlayer    !!!!!!!!!!!!!!!!
    do i = 1, sys%pl_z_l(1)%length()
      js = sys%pl_z_l(1)%value(i)%i
      ks = sys%pl_z_l(1)%value(i)%f
      j = sys%ste(js)
      k = sys%ste(ks)
      ham(4*j-3, 4*k-3) = -sys%t(3)*(1d0-ne_d(js))*(1d0-ne_d(ks))
      ham(4*j-2, 4*k-2) = -sys%t(3)*(1d0-ne_u(js))*(1d0-ne_u(ks))
      ham(4*j-1, 4*k-1) = -sys%t(3)*(1d0-ne_d(js))*(1d0-ne_d(ks))
      ham(4*j  , 4*k  ) = -sys%t(3)*(1d0-ne_u(js))*(1d0-ne_u(ks))
      ham(4*k-3, 4*j-3) = -sys%t(3)*(1d0-ne_d(js))*(1d0-ne_d(ks))
      ham(4*k-2, 4*j-2) = -sys%t(3)*(1d0-ne_u(js))*(1d0-ne_u(ks))
      ham(4*k-1, 4*j-1) = -sys%t(3)*(1d0-ne_d(js))*(1d0-ne_d(ks))
      ham(4*k  , 4*j  ) = -sys%t(3)*(1d0-ne_u(js))*(1d0-ne_u(ks))
      
      ham(4*j-3, 4*j-3) = ham(4*j-3, 4*j-3) + t3U*Sz(ks)
      ham(4*j-2, 4*j-2) = ham(4*j-2, 4*j-2) - t3U*Sz(ks)
      ham(4*j-3, 4*j-2) = ham(4*j-3, 4*j-2) + t3U*(Sx(ks) - ui*Sy(ks))
      ham(4*j-2, 4*j-3) = ham(4*j-2, 4*j-3) + t3U*(Sx(ks) + ui*Sy(ks))
      ham(4*j-1, 4*j-1) = ham(4*j-1, 4*j-1) - t3U*Sz(ks) ! -conjg
      ham(4*j  , 4*j  ) = ham(4*j  , 4*j  ) + t3U*Sz(ks) ! -conjg
      ham(4*j-1, 4*j  ) = ham(4*j-1, 4*j  ) + conjg( t3U*(Sx(ks) - ui*Sy(ks)) ) ! + conjg
      ham(4*j  , 4*j-1) = ham(4*j  , 4*j-1) + conjg( t3U*(Sx(ks) + ui*Sy(ks)) ) ! + conjg
      
      ham(4*k-3, 4*k-3) = ham(4*k-3, 4*k-3) + t3U*Sz(js)
      ham(4*k-2, 4*k-2) = ham(4*k-2, 4*k-2) - t3U*Sz(js)
      ham(4*k-3, 4*k-2) = ham(4*k-3, 4*k-2) + t3U*(Sx(js) - ui*Sy(js))
      ham(4*k-2, 4*k-3) = ham(4*k-2, 4*k-3) + t3U*(Sx(js) + ui*Sy(js))
      ham(4*k-1, 4*k-1) = ham(4*k-1, 4*k-1) - t3U*Sz(js) ! -conjg
      ham(4*k  , 4*k  ) = ham(4*k  , 4*k  ) + t3U*Sz(js) ! -conjg
      ham(4*k-1, 4*k  ) = ham(4*k-1, 4*k  ) + conjg( t3U*(Sx(js) - ui*Sy(js)) ) ! + conjg
      ham(4*k  , 4*k-1) = ham(4*k  , 4*k-1) + conjg( t3U*(Sx(js) + ui*Sy(js)) ) ! + conjg
    end do
    !!!!!!!!!!!!!!   end  interlayer    !!!!!!!!!!!!!!!!
  
    call my_zheev("l", ham, eg, wf)

  end subroutine diagonalize_sb_2lay_hamiltonian

END MODULE car_parrinello_mod