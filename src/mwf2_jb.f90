PROGRAM make_wave_function2_jb
  USE math_mod
  !USE energy_mod,            only : total_energy_sb_2lay_xichi_T0, calc_eigen_energy_sb_2lay
  USE energy_mod,            only : total_energy_sb_2lay_xichi_T0_Aem, calc_eigen_energy_sb_2lay_Aem
  USE io_mod,                ONLY : save_txt, load_txt, create_vortex, load_conf_3D_2 ,load_conf_3D, load_conf_chi_w, load_jbconf
  USE visual_mod,            ONLY : draw_phase_3D, write_format_of_xi_3D, draw_path_3D, write_format_of_chi_sb, &
       &                              write_format_of_chi_3D, print_phase,draw_circle, draw_current, draw_winding_circle
  USE angular_variable_mod,  ONLY : angular_variable, new_angular_variable, new_angular_variable_3D
  !USE angular_variable_3D_mod,  ONLY : angular_variable, new_angular_variable
  USE car_parrinello_mod,    ONLY : car_parrinello, new_car_parrinello
  !USE car_parrinello_wf_mod, ONLY : car_parrinello_wf, new_car_parrinello_wf
  !  use sdchi_optimizer_mod,    only : sdchi_optimizer, new_sdchi_optimizer
  !USE fchi_optimizer_mod,    ONLY : fchi_optimizer, new_fchi_optimizer
  !USE basic_current_mod,     ONLY : basic_current, new_basic_current
  !USE single_wf_mod,         ONLY : single_wf, new_single_wf
  USE wave_function_mod,      ONLY : wave_function, new_wave_function
  use conjugate_gradient_rsoc_mod, only : conjugate_gradient_rsoc,new_conjugate_gradient_rsoc
  !use conjugate_gradient_spin_mod
  use parameter_mod
  use hop_iterator_mod
  use path_list_mod

  !USE car_parrinello_mod,    ONLY : diagonalize_surface_hamiltonian, calc_surface_ne, calc_surface_S
  !USE car_parrinello_mod,    ONLY : diagonalize_bulk_hamiltonian, calc_bulk_ne, calc_bulk_S
  !USE car_parrinello_mod,    ONLY : diagonalize_sb_2lay_hamiltonian, calc_sb_ne
  !use diagonalize_mod,   only : diag, new_diag, diagonalize_sb_2lay_xichi_hamiltonian
  use diagonalize_mod,   only : diag, new_diag
  use measurement_mod

  use winding_number_mod
  use circle_list_mod
  !use fchi_optimizer_mod
  use fchi_optimizer_Aem_mod
  !use energy_mod , only : calc_eigen_energy_sb_2lay_xichi
  !use energy_mod , only : calc_eigen_energy_sb_2lay
  !USE energy_mod , only : calc_eigen_energy_sb_2lay_wf_T0
  use utility_mod, only : sort_wf

  IMPLICIT NONE

  INTEGER                                  :: nx, ny, nz, n, nh, n_acc, r(3) , time_s,time_e,CountPerSec, CountMax 
  !real(8)                                  :: rand , ph(3) 
  REAL(8), DIMENSION(:),allocatable        :: ne
  INTEGER, DIMENSION(:,:),allocatable      :: pshape
  REAL(8), DIMENSION(3), PARAMETER         :: t = [1d0, -0.12d0, 0.01d0]  ! hopping coeff[1st, 2nd, layer]
  !REAL(8), DIMENSION(3), PARAMETER         :: t = [1d0, -0.12d0, 0d0]  ! hopping coeff[1st, 2nd, layer]
  !REAL(8), DIMENSION(3), PARAMETER         :: t = [1d0, 0d0, 0.01d0]  ! hopping coeff[1st, 2nd, layer]
  REAL(8), PARAMETER                       :: u = 8.d0         ! coulomb coeff
  REAL(8), PARAMETER                       :: jd = 0.5d0*4.0d0*(t(1)**2)/u ! across hole coeff
  !REAL(8), PARAMETER                       :: jd = 0.0d0 ! across hole coeff
  !REAL(8), PARAMETER                       :: lm = 0.02d0 ! rashba coeff
  REAL(8), PARAMETER                       :: lm = 0.01d0 ! rashba coeff
  !REAL(8), PARAMETER                       :: lm = -0.02d0 ! rashba coeff
  !REAL(8), PARAMETER                       :: lm = 0.2d0 ! rashba coeff
  !REAL(8), PARAMETER                       :: lm = 0.0d0 ! rashba coeff
  !REAL(8), PARAMETER                       :: Bz = 0.01d0 ! z-dir magnetic field
  !REAL(8)                                  :: Bz = 0.01d0 ! z-dir magnetic field
  real(8), parameter                       :: kbt = 0.01d0 ! temperature
  LOGICAL,PARAMETER                        :: sflag=.TRUE.!slanting path
  !  logical,parameter                        :: sflag=.false.!no slanting path
  INTEGER, PARAMETER                       :: unit_num = 10
  INTEGER, DIMENSION(:, :), ALLOCATABLE    :: hole_xi, hole_chi !dim is (holes,2)
  INTEGER, DIMENSION(:, :), ALLOCATABLE    :: tmp_hole_xi !dim is (holes,2)
  !type(hole_info), DIMENSION(:), ALLOCATABLE    :: hole_layer
  REAL(8), DIMENSION(:), ALLOCATABLE       :: xi,chi, Sx, Sy, Sz, dummy, eta, chi0, xi0
  TYPE(angular_variable), ALLOCATABLE      :: an_xi,an_chi
  !TYPE(basic_current), ALLOCATABLE         :: bc
  !TYPE(car_parrinello_wf),ALLOCATABLE      :: cpw
  !TYPE(car_parrinello),ALLOCATABLE         :: cp
  !  type(sdchi_optimizer),allocatable         :: opt
  !TYPE(fchi_optimizer),ALLOCATABLE         :: fopt
  TYPE(fchi_optimizer_Aem),ALLOCATABLE         :: fopta
  !TYPE(single_wf), ALLOCATABLE             :: sw
  !TYPE(wave_function), ALLOCATABLE         :: wf
  !type(conjugate_gradient_rsoc),allocatable:: cgr
  COMPLEX(8), DIMENSION(:, :), ALLOCATABLE :: wf, wf1, wf2
  REAL(8), DIMENSION(:), ALLOCATABLE       :: eg
  INTEGER                                  :: i, j, k, js, ks
  !REAL(8), DIMENSION(:), ALLOCATABLE       :: delta_xi
  logical                                 :: dummy_flag
  !type(conjugate_gradient_spin),allocatable :: cgspin
  type(conjugate_gradient_rsoc),allocatable:: cgr
  logical                                  :: optimized_cgr, optimized_E_tot
  REAL(8)                                  :: p(1)
  real(8)                                  :: total_E, total_E2, xi1, old_E_tot, new_E_tot
  !type(cuo2_plane),allocatable, dimension(:) :: layers
  !type(real_value_plane),allocatable, dimension(:) :: xi_l, chi_l
  !integer                                          :: n_l
  !type(hop_iterator), allocatable :: iter
  !type(path_list), allocatable :: pl
  !type(hop_iterator), allocatable :: iter
  !real(8) :: sum_abs_u,sum_abs_v, mu(2)
  real(8) :: mu(2)
  !real(8), dimension(:), allocatable :: ne_u, ne_d
  complex(8), dimension(:), allocatable :: delta
  complex(8), dimension(:), allocatable :: tdelta_for, tdelta_rev
  !integer :: bs, be
  !REAL(8), DIMENSION(:), ALLOCATABLE       :: bSx, bSy, bSz, sSx, sSy, sSz, bxi, bchi
  !real(8), dimension(:), allocatable :: bne_u, bne_d, sne_u, sne_d
  !COMPLEX(8), DIMENSION(:, :), ALLOCATABLE :: bwf, swf
  !REAL(8), DIMENSION(:), ALLOCATABLE       :: beg, seg
  !real(8) :: vec(3), rvec(3)
  !REAL(8), DIMENSION(:), ALLOCATABLE       :: srho, e_value, rho
  !integer :: stm_site
  !real(8),dimension(:),allocatable        :: jex
  type(circle_list),allocatable :: clist
  integer,dimension(:),allocatable :: wlist
  integer :: sr(3)
  real(8), dimension(3) :: pole
  real(8), dimension(:,:), allocatable :: pole_s, current, current_0, current_diff
  integer, dimension(:), allocatable :: w_s
  real(8), dimension(:), allocatable :: ne_u, ne_d
  real(8), dimension(:,:), allocatable :: pole_s_conv
  integer, dimension(:), allocatable :: w_s_conv

  type(diag), allocatable :: dg
  real(8), dimension(:), allocatable :: tne_u, tne_d
  real(8) :: ne_surf, ne_bulk
  integer :: n_max_surf
  real(8) :: mu_l, mu_c, mu_r

  real(8), allocatable, dimension(:) :: Jp, Japp, Jeff, J0, Jdiff
  real(8) :: current_ratio
  
  real(8) :: Et_s, Et2_s, Et_b, Et2_b, Ersoc, E_interlayer, E_total
  real(8) :: Et_s0, Et2_s0, Et_b0, Et2_b0, Ersoc0, E_interlayer0, E_total0, E_total0_other

  real(8), dimension(:), allocatable        :: jex, jex_from, jex_to
  integer N_jex, site_obs

  real(8) :: Bz

  integer :: jex_iter, b_iter

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!         options        !!!!!!!!!!!!
  
  !   surface_chi_zero
  ! true : set all winding numbers of chi on the surface to zero 
  !        this is expected that surface current does not flow 
  ! false : winding numbers of chi on surface are same one on bulk
  !logical, parameter :: surface_chi_zero = .true.
  logical, parameter :: surface_chi_zero = .false.
  
  !   use_self_consistent
  ! true  : using self-consistent in fchi-optimizer
  ! false : using one-shot in fchi-optimizer
  !logical, parameter :: use_self_consistent = .true.
  logical, parameter :: use_self_consistent = .false.

  !   use_xi1_optimization
  !logical, parameter :: use_xi1_optimization = .true.
  logical, parameter :: use_xi1_optimization = .false.
  !!!!!!!!!!                        !!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !real(8) :: e, kbt(10), dos(10), rho(10), eg_range(2)
  !integer :: n_step, site
  
  print *, "make wave function 2 (calclation of car parrinello)"
  
  ! load hole position, xi winding number and chi winding number
  CALL load_conf_3D_2('./config/vconf', unit_num, pshape, hole_xi,hole_chi, ne)
  print *, "load : ", "./config/vconf"

  !CALL load_jexconf_3D('./config/jexconf', unit_num, pshape, jex)
  CALL load_jbconf('./config/jbconf', unit_num, pshape, jex_from, jex_to, N_jex, site_obs, Bz)
  print *, "load : ", "./config/jbconf"

  nz = size(pshape(1,:))
  !nh =  sum(layers%nh)
  nh =  size(hole_chi(:, 1))
  !nx = pshape(1,1)
  !ny = pshape(2,1)
  !n = nx * ny
  !n = layers(size(layers))%site_end
  n = site_max(pshape)
  n_acc = n - nh

  
  WRITE(*,*) "-----------------------------------------------------------"
  WRITE(*,*) "config:"

  !!check
  !WRITE(*,*) 'nx,ny,nz', nx, ny, nz
  WRITE(*,*) 'n, ne', n, ne
  WRITE(*,*) 'nh(1), nh(2)', SIZE(hole_xi,1), SIZE(hole_xi,2)
  WRITE(*,*) 'hole_xi:', hole_xi(:,:)
  WRITE(*,*) 'hole_chi:', hole_chi(:,:)
  
  WRITE(*,'(a,2(f7.3,",",2x),f7.3)') "t(1st, 2nd, inter-layer) : ", t
  WRITE(*,'(a,f7.3)') "U : ", U
  WRITE(*,'(a,f7.3)') "lambda (rashba) : ", lm
  WRITE(*,'(a,f7.3)') "Jd : ", jd
  WRITE(*,'(a,f7.3)') "Bz : ", Bz

  WRITE(*,'(a,l)') "surface_chi_zero : ", surface_chi_zero
  WRITE(*,'(a,l)') "self consistent  : ", use_self_consistent
  WRITE(*,'(a,l)') "xi1 optimization : ", use_xi1_optimization

  WRITE(*,*) "-----------------------------------------------------------"
  
  ALLOCATE(xi(n))
  ALLOCATE(xi0(n))
  ALLOCATE(chi(n),source=0.d0)
  ALLOCATE(chi0(n),source=0.d0)
  ALLOCATE(eta(n))
  ALLOCATE(wf(4*n_acc, 4*n_acc))
  ALLOCATE(eg(4*n_acc))
  ALLOCATE(current(n,n))
  ALLOCATE(current_0(n,n))
  ALLOCATE(current_diff(n,n))
  ALLOCATE(ne_u(n))
  ALLOCATE(ne_d(n))
  
  allocate(dg)
  dg = new_diag(pshape, hole_xi, t, u, jd, lm)
  ALLOCATE(tdelta_for(dg%pl_1st_l(1)%length()))
  ALLOCATE(tdelta_rev(dg%pl_1st_l(1)%length()))
  deallocate(dg)

  ALLOCATE(Sx(n))
  ALLOCATE(Sy(n))
  ALLOCATE(Sz(n))
  ALLOCATE(tne_u(n))
  ALLOCATE(tne_d(n))
  allocate(jex(n), source=0d0)

  read(50) wf
  read(51) eg
  read(52) xi
  !read(53) chi
  read(54) ne_u
  read(55) ne_d
  read(58) mu

  xi0 = xi

  do i = 1, n
  sr = site_to_r(i, pshape)
    !eta(i) = (sr(1)+sr(2))*pi + xi(i)
    eta(i) = (sr(1)+sr(2)+sr(3))*pi + xi(i)
  end do
  call draw_phase_3D('eta.dat', unit_num, pshape, eta, hole_xi(:, 1))

  !clist = new_circle_list(pshape, hole_chi)
  clist = new_circle_list(pshape, hole_chi, .false.)
  call get_winding_number(eta, clist, wlist)
  !call draw_circle("wl_eta.dat", 50, pshape, clist%cil, clist%cpl, wlist)
  call draw_winding_circle("wp_eta.dat", "wm_eta.dat", 50, 51, pshape, clist%cil, clist%cpl, wlist)
  print *, "save : ""wp_eta.dat"", ""wm_eta.dat"" "
  !wlist(1) = -wlist(1)
  call write_format_of_chi_3D('chi.plt', unit_num, pshape, hole_chi)
  !wlist = - wlist

  ! set surface poles
  k = 0
  do i = 1, size(clist%cil)
    pole = clist%get_circle_center(i)
    if(wlist(i) /= 0 .and. nint(pole(3)) == 1) then
      k = k + 1;
    end if
  end do
  allocate(pole_s(2, k))
  allocate(w_s(k))
  pole_s(:, :) = 0d0
  k = 0
  do i = 1, size(clist%cil)
    pole = clist%get_circle_center(i)
    if(wlist(i) /= 0 .and. nint(pole(3)) == 1) then
      k = k + 1
      pole_s(:, k) = pole(1:2)
      w_s(k) = wlist(i)
      !w_s(k) = -wlist(i)
    end if
  end do

  if(surface_chi_zero) then
    print *, "set all winding numbers of chi on the surface to zero"
    w_s(:) = 0
  end if
  
  call write_format_of_chi_sb('chi_sb_init.plt', unit_num, pshape, hole_chi, pole_s, w_s)

  allocate(an_chi)
  an_chi = new_angular_variable_3D(pshape, hole_chi)
  chi0 = an_chi%get_angular_variable() ! bulk
  chi0 = chi0 + an_chi%get_atan_phase_on_xy(pole_s, w_s, 1) ! surface
  chi0 = an_chi%rebuilt_angval_chi(chi0,0.0d0*pi)

  CALL draw_phase_3D('init_chi.dat', unit_num, pshape, chi0, hole_chi(:, 1))
  
  call get_winding_number(chi0, clist, wlist)
  call draw_winding_circle("wp_chi_init.dat", "wm_chi_init.dat", 50, 51, pshape, clist%cil, clist%cpl, wlist)
  
  open(21, file = "jeg.dat")
  open(22, file = "jeg_diff.dat")
  open(23, file = "jeg_diff_2.dat")
  !write(21, '("# Bz = ", f)') Bz

  !!!! fchi test
  !ALLOCATE(fopta)
  !!xi = xi0
  !fopta = new_fchi_optimizer_Aem(pshape, hole_chi, ne,  t, u, jd, lm, Bz, xi, chi0, wf, eg, ne_u, ne_d) !slanting path
  !jex = jex_from
  !call fopta%set_jex(jex)
  !!call fopta%set_param(pflag = .false.)
  !if (use_self_consistent) then
  !  CALL fopta%self_consistent()
  !else
  !  CALL fopta%oneshot()
  !end if
  !chi = an_chi%get_rebuilt_phase(fopta%delta_chi)
  !if (lm /= 0.d0 .and. use_xi1_optimization == .true.) then
  !!if (lm /= 0.d0) then
  !  old_E_tot = total_energy_sb_2lay_xichi_T0_Aem(fopta, wf, xi, chi, fopta%Aem_dphase)
  !  ALLOCATE(an_xi)
  !  an_xi = new_angular_variable_3D (pshape,hole_xi)
  !  cgr = new_conjugate_gradient_rsoc(n, n_acc, lm, an_xi, wf, xi, chi)
  !  p(1) = 0.6d0*pi ! initial xi1
  !  call cgr%calc(p,  optimized_cgr)
  !  p(1) = modulo(p(1), 2.d0*pi)
  !  xi = an_xi%rebuilt_angval_xi(xi, p(1))
  !  new_E_tot = total_energy_sb_2lay_xichi_T0_Aem(fopta, wf, xi, chi, fopta%Aem_dphase)
  !  DEALLOCATE(an_xi)
  !  call fopta%modify_init(hole_chi, xi, chi0, wf)
  !  if (use_self_consistent) then
  !    CALL fopta%self_consistent()
  !  else
  !    CALL fopta%oneshot()
  !  end if
  !  chi = an_chi%get_rebuilt_phase(fopta%delta_chi)
  !endif
  !E_total = total_energy_sb_2lay_xichi_T0_Aem(fopta, wf, xi, chi, fopta%Aem_dphase, &
  !                  & Et_s0, Et2_s0, Et_b0, Et2_b0, Ersoc0, E_interlayer0 )
  !print *, "E_total : ", E_total 
  !call fopta%get_Jp_matrix(current)
  !call draw_current("current.dat", 101, pshape, fopta%cpl, current)
  !print *, "save : ""current.dat"""
  !DEALLOCATE(fopta)
  !stop "stop : mwf2 test"

  
  !!!! xi1 test
  !open(23, file = "xi1.dat")
  !ALLOCATE(an_xi)
  !an_xi = new_angular_variable_3D (pshape,hole_xi)
  !ALLOCATE(fopta)
  !fopta = new_fchi_optimizer_Aem(pshape, hole_chi, ne,  t, u, jd, lm, Bz, xi, chi0, wf, eg, ne_u, ne_d) !slanting path
  !!CALL fopta%oneshot()
  !CALL fopta%self_consistent()
  !chi = an_chi%get_rebuilt_phase(fopta%delta_chi)
  !do i = 0, 100
  !  p(1) = i/100d0*2d0*pi
  !  xi = an_xi%rebuilt_angval_xi(xi, p(1))
  !  new_E_tot = total_energy_sb_2lay_xichi_T0_Aem(fopta, wf, xi, chi, fopta%Aem_dphase, oHrsoc = Ersoc)
  !  write(23, *) p(1), new_E_tot, Ersoc
  !end do
  !DEALLOCATE(fopta)
  !DEALLOCATE(an_xi)
  !close(23)
  !stop "stop : xi1 test"

  !!!!! London eq. test
  !ALLOCATE(fopta)
  !fopta = new_fchi_optimizer_Aem(pshape, hole_chi, ne,  t, u, jd, lm, 0d0, xi, chi0, wf, eg, ne_u, ne_d) !slanting path
  !if (use_self_consistent) then
  !  CALL fopta%self_consistent()
  !else
  !  CALL fopta%oneshot()
  !end if
  !chi = an_chi%get_rebuilt_phase(fopta%delta_chi)
  !call fopta%get_Jp_matrix(current_0)
  !DEALLOCATE(fopta)
  
  !ALLOCATE(fopta)
  !!xi = xi0
  !fopta = new_fchi_optimizer_Aem(pshape, hole_chi, ne,  t, u, jd, lm, Bz, xi, chi0, wf, eg, ne_u, ne_d) !slanting path
  !if (use_self_consistent) then
  !  CALL fopta%self_consistent()
  !else
  !  CALL fopta%oneshot()
  !end if
  !chi = an_chi%get_rebuilt_phase(fopta%delta_chi)
  !allocate(Jp(fopta%np))
  !allocate(Japp(fopta%np))
  !allocate(Jeff(fopta%np))
  !allocate(J0(fopta%np))
  !allocate(Jdiff(fopta%np))
  !Jp = fopta%get_Jp(fopta%u_Aeff)  
  !Jeff = fopta%approximated_J()
  !J0 = fopta%J0()
  !Japp = J0 + Jeff
  !current_ratio = 1d0 / maxval(abs(Jp))
  
  !call fopta%make_jmatrix(Jp, current)
  !call draw_current_ratio("current_Jp.dat", 101, pshape, fopta%cpl, current, current_ratio)
  !current = current - current_0
  !call draw_current_ratio("current_Jp_diff.dat", 101, pshape, fopta%cpl, current, current_ratio)
  !call draw_current("current_Jp_diff_autoscale.dat", 101, pshape, fopta%cpl, current)
  !call fopta%make_jmatrix(Jeff, current)
  !call draw_current_ratio("current_Jeff.dat", 101, pshape, fopta%cpl, current, current_ratio)
  !call fopta%make_jmatrix(J0, current)
  !call draw_current_ratio("current_J0.dat", 101, pshape, fopta%cpl, current, current_ratio)
  !call draw_current("current_J0_autoscale.dat", 101, pshape, fopta%cpl, current)
  !call fopta%make_jmatrix(Japp, current)
  !call draw_current_ratio("current_Japp.dat", 101, pshape, fopta%cpl, current, current_ratio)
  
  !Jdiff = Jp - Japp
  !call fopta%make_jmatrix(Jdiff, current)
  !call draw_current_ratio("current_Jdiff.dat", 101, pshape, fopta%cpl, current, current_ratio)
  !call draw_current("current_Jdiff_autoscale.dat", 101, pshape, fopta%cpl, current)
  !print *, "max error : ", maxval(abs(Jdiff))
  !print *, "max Jp : ", maxval(abs(Jp))
  !print *, "error ratio: ", maxval(abs(Jdiff))/maxval(abs(Jp))

  !stop "stop : London eq. test"

  !!! when jex = 0
  ALLOCATE(fopta)
  !xi = xi0
  fopta = new_fchi_optimizer_Aem(pshape, hole_chi, ne,  t, u, jd, lm, Bz, xi, chi0, wf, eg, ne_u, ne_d) !slanting path
  call fopta%set_param(pflag = .false.)
  if (use_self_consistent) then
    CALL fopta%self_consistent()
  else
    CALL fopta%oneshot()
  end if
  chi = an_chi%get_rebuilt_phase(fopta%delta_chi)
  if (lm /= 0.d0 .and. use_xi1_optimization == .true.) then
  !if (lm /= 0.d0) then
    old_E_tot = total_energy_sb_2lay_xichi_T0_Aem(fopta, wf, xi, chi, fopta%Aem_dphase)
    ALLOCATE(an_xi)
    an_xi = new_angular_variable_3D (pshape,hole_xi)
    cgr = new_conjugate_gradient_rsoc(n, n_acc, lm, an_xi, wf, xi, chi)
    p(1) = 0.6d0*pi ! initial xi1
    call cgr%calc(p,  optimized_cgr)
    p(1) = modulo(p(1), 2.d0*pi)
    xi = an_xi%rebuilt_angval_xi(xi, p(1))
    new_E_tot = total_energy_sb_2lay_xichi_T0_Aem(fopta, wf, xi, chi, fopta%Aem_dphase)
    DEALLOCATE(an_xi)
    call fopta%modify_init(hole_chi, xi, chi0, wf)
    if (use_self_consistent) then
      CALL fopta%self_consistent()
    else
      CALL fopta%oneshot()
    end if
    chi = an_chi%get_rebuilt_phase(fopta%delta_chi)
  endif
  call fopta%get_Jp_matrix(current_0)
  !call draw_current("current_0.dat", 101, pshape, fopta%cpl, current)
  E_total0 = total_energy_sb_2lay_xichi_T0_Aem(fopta, wf, xi, chi, fopta%Aem_dphase, &
                    & Et_s0, Et2_s0, Et_b0, Et2_b0, Ersoc0, E_interlayer0 )
  DEALLOCATE(fopta)
  
  !write(60) E_total0
  !stop "stop : out E_total0"
  read(60) E_total0_other


  do jex_iter = 0, N_jex
    !xi = xi0
    !Bz = real(b_iter, 8)/10000d0 * 1d0
    jex = jex_from + (jex_to - jex_from)*real(jex_iter,8)/real(N_jex,8)
    !ALLOCATE(fopt)
    ALLOCATE(fopta)
    !fopt = new_fchi_optimizer(pshape, hole_chi, ne,  t, u, jd, lm, xi, chi0, wf, eg, ne_u, ne_d) !slanting path
    fopta = new_fchi_optimizer_Aem(pshape, hole_chi, ne,  t, u, jd, lm, Bz, xi, chi0, wf, eg, ne_u, ne_d) !slanting path
    !call fopt%set_jex(jex)
    call fopta%set_jex(jex)
    call fopta%set_param(pflag = .false.)

    !call fopt%get_Jp_matrix(current)
    !call draw_current("init_current.dat", 101, pshape, fopt%cpl, current)
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!! calc fchi optimizer (1/2) !!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (use_self_consistent) then
      !CALL fopt%self_consistent()
      CALL fopta%self_consistent()
    else
      !CALL fopt%oneshot()
      CALL fopta%oneshot()
    end if
    !chi = an_chi%get_rebuilt_phase(fopt%delta_chi)
    chi = an_chi%get_rebuilt_phase(fopta%delta_chi)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!! optimize xi1         !!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
    if (lm /= 0.d0 .and. use_xi1_optimization == .true.) then
      !write(*, *) '++++++Optimize xi1++++++++++'
      !old_E_tot = total_energy_sb_2lay_xichi_T0(fopt, wf, xi, chi)
      old_E_tot = total_energy_sb_2lay_xichi_T0_Aem(fopta, wf, xi, chi, fopta%Aem_dphase)
      !write(*, *) 'init E_tot', old_E_tot
      ALLOCATE(an_xi)
      an_xi = new_angular_variable_3D (pshape,hole_xi)
      cgr = new_conjugate_gradient_rsoc(n, n_acc, lm, an_xi, wf, xi, chi)
      p(1) = 0.6d0*pi ! initial xi1
      call cgr%calc(p,  optimized_cgr)
      p(1) = modulo(p(1), 2.d0*pi)
      xi = an_xi%rebuilt_angval_xi(xi, p(1))
      !write(*, *) 'optimized xi1', xi(1)
      !new_E_tot = total_energy_sb_2lay_xichi_T0(fopt, wf, xi, chi)
      new_E_tot = total_energy_sb_2lay_xichi_T0_Aem(fopta, wf, xi, chi, fopta%Aem_dphase)
      !write(*, *) 'E_tot', new_E_tot
      !write(*, *) '++++++++end xi1 optimization+++++++++++++++'
      DEALLOCATE(an_xi)
      call fopta%modify_init(hole_chi, xi, chi0, wf)
      if (use_self_consistent) then
        CALL fopta%self_consistent()
      else
        CALL fopta%oneshot()
      end if
      chi = an_chi%get_rebuilt_phase(fopta%delta_chi)
    endif

    allocate(dg)
    dg = new_diag(pshape, hole_xi, t, u, jd, lm)
    
    dg%wf = wf
    dg%eg = eg
    call dg%calc_tdelta_xi(tdelta_for, tdelta_rev, xi, 0.01d0)
    call dg%calc_spin_xi(Sx, Sy, Sz, xi, 0.01d0)

    call calc_eigen_energy_sb_2lay_Aem(fopta, eg, wf, xi, chi, mu, fopta%Aem_dphase, tdelta_for, tdelta_rev, ne_u, ne_d, Sx, Sy, Sz)

    !dg%wf = wf
    !dg%eg = eg
    fopta%wf = wf
    fopta%eg = eg
    !call dg%calc_ne(ne_u, ne_d, kbt)
    !new_E_tot = total_energy_sb_2lay_xichi_T0_Aem(fopta, wf, xi, chi, fopta%Aem_dphase)
    E_total = total_energy_sb_2lay_xichi_T0_Aem(fopta, wf, xi, chi, fopta%Aem_dphase, &
                    & Et_s, Et2_s, Et_b, Et2_b, Ersoc, E_interlayer )
    !print *, 'E_tot(T=0)', E_total
    !write(21, *) Bz, new_E_tot
    !write(21, *) jex(site_obs), new_E_tot
    write(21, '(8f)') jex(site_obs), E_total, Et_s, Et2_s, Et_b, Et2_b, Ersoc, E_interlayer
    write(22, '(8f)') jex(site_obs), E_total - E_total0, Et_s - Et_s0, Et2_s - Et2_s0, &
                    & Et_b - Et_b0, Et2_b - Et2_b0, Ersoc - Ersoc0, E_interlayer - E_interlayer0
    write(23, '(2f)') jex(site_obs), E_total - E_total0_other
    
    if(jex_iter == 0) then
      call fopta%get_Jp_matrix(current)
      call draw_current("current.dat", 101, pshape, fopta%cpl, current)
      call draw_current("current_0.dat", 101, pshape, fopta%cpl, current_0)
      print *, "save : ""current.dat"""
      current_diff = current - current_0
      call draw_current("current_diff.dat", 101, pshape, fopta%cpl, current_diff)
      call draw_current_ratio("current_diff_e.dat", 101, pshape, fopta%cpl, current_diff, 1d0/maxval(abs(current)))
      print *, "save : ""current_diff.dat"", ""current_diff_e.dat"""
    end if

    deallocate(dg)
    deallocate(fopta)
  end do ! jex_iter

  print *, "save : ""jeg.dat"", ""jeg_diff.dat"", ""jeg_diff_2.dat"""
  close(21)
  close(22)
  close(23)

  stop "stop : mwf2_jb"

END PROGRAM make_wave_function2_jb