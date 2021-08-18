PROGRAM make_wave_function2_be
  USE math_mod
  !USE energy_mod,            only : total_energy_sb_2lay_xichi_T0, calc_eigen_energy_sb_2lay
  USE energy_mod,            only : total_energy_sb_2lay_xichi_T0_Aem, calc_eigen_energy_sb_2lay_Aem
  USE io_mod,                ONLY : save_txt, load_txt, create_vortex, load_conf_3D_2 ,load_conf_3D, load_conf_chi_w
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
  real(8)                                  :: E_total ,rand , ph(3) 
  REAL(8), DIMENSION(:),allocatable        :: ne
  INTEGER, DIMENSION(:,:),allocatable      :: pshape
  REAL(8), DIMENSION(3), PARAMETER         :: t = [1d0, -0.12d0, 0.01d0]  ! hopping coeff[1st, 2nd, layer]
  !REAL(8), DIMENSION(3), PARAMETER         :: t = [1d0, -0.12d0, 0d0]  ! hopping coeff[1st, 2nd, layer]
  !REAL(8), DIMENSION(3), PARAMETER         :: t = [1d0, 0d0, 0.1d0]  ! hopping coeff[1st, 2nd, layer]
  REAL(8), PARAMETER                       :: u = 8.d0         ! coulomb coeff
  REAL(8), PARAMETER                       :: jd = 0.5d0*4.0d0*(t(1)**2)/u ! across hole coeff
  !REAL(8), PARAMETER                       :: jd = 0.0d0 ! across hole coeff
  REAL(8), PARAMETER                       :: lm = 0.02d0 ! rashba coeff
  !REAL(8), PARAMETER                       :: lm = 0.2d0 ! rashba coeff
  !REAL(8), PARAMETER                       :: lm = 0.0d0 ! rashba coeff
  !REAL(8), PARAMETER                       :: Bz = 0.01d0 ! z-dir magnetic field
  REAL(8)                                  :: Bz = 0.01d0 ! z-dir magnetic field
  real(8), parameter                       :: kbt = 0.01d0 ! temperature
  LOGICAL,PARAMETER                        :: sflag=.TRUE.!slanting path
  !  logical,parameter                        :: sflag=.false.!no slanting path
  INTEGER, PARAMETER                       :: unit_num = 10
  INTEGER, DIMENSION(:, :), ALLOCATABLE    :: hole_xi, hole_chi !dim is (holes,2)
  INTEGER, DIMENSION(:, :), ALLOCATABLE    :: tmp_hole_xi !dim is (holes,2)
  !type(hole_info), DIMENSION(:), ALLOCATABLE    :: hole_layer
  REAL(8), DIMENSION(:), ALLOCATABLE       :: xi,chi, Sx, Sy, Sz, dummy, eta, chi0
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
  real(8), dimension(:,:), allocatable :: pole_s, current
  integer, dimension(:), allocatable :: w_s
  real(8), dimension(:), allocatable :: ne_u, ne_d
  real(8), dimension(:,:), allocatable :: pole_s_conv
  integer, dimension(:), allocatable :: w_s_conv

  type(diag), allocatable :: dg
  real(8), dimension(:), allocatable :: tne_u, tne_d
  real(8) :: ne_surf, ne_bulk
  integer :: n_max_surf
  real(8) :: mu_l, mu_c, mu_r

  integer :: b_iter

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!         options        !!!!!!!!!!!!
  
  !   surface_chi_zero
  ! true : set all winding numbers of chi on the surface to zero 
  !        this is expected that surface current does not flow 
  ! false : winding numbers of chi on surface are same one on bulk
  logical, parameter :: surface_chi_zero = .false.
  
  !   use_self_consistent
  ! true  : using self-consistent in fchi-optimizer
  ! false : using one-shot in fchi-optimizer
  logical, parameter :: use_self_consistent = .false.
  !logical, parameter :: use_self_consistent = .true.

  !!!!!!!!!!                        !!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !real(8) :: e, kbt(10), dos(10), rho(10), eg_range(2)
  !integer :: n_step, site
  
  print *, "make wave function 2 (calclation of car parrinello)"
  ! load hole position, xi winding number and chi winding number
  CALL load_conf_3D_2('./config/vconf', unit_num, pshape, hole_xi,hole_chi, ne)

  !CALL load_jexconf_3D('./config/jexconf', unit_num, pshape, jex)

  nz = size(pshape(1,:))
  !nh =  sum(layers%nh)
  nh =  size(hole_chi(:, 1))
  !nx = pshape(1,1)
  !ny = pshape(2,1)
  !n = nx * ny
  !n = layers(size(layers))%site_end
  n = site_max(pshape)
  n_acc = n - nh

  !!check
  !WRITE(*,*) 'nx,ny,nz', nx, ny, nz
  !WRITE(*,*) 'n, ne', n, ne
  !WRITE(*,*) 'nh(1), nh(2)', SIZE(hole_xi,1), SIZE(hole_xi,2)
  !WRITE(*,*) 'hole_xi:', hole_xi(:,:)
  !WRITE(*,*) 'hole_chi:', hole_chi(:,:)
  ALLOCATE(xi(n))
  ALLOCATE(chi(n),source=0.d0)
  ALLOCATE(chi0(n),source=0.d0)
  ALLOCATE(eta(n))
  ALLOCATE(wf(4*n_acc, 4*n_acc))
  ALLOCATE(eg(4*n_acc))
  ALLOCATE(current(n,n))
  ALLOCATE(ne_u(n))
  ALLOCATE(ne_d(n))

  read(50) wf
  read(51) eg
  read(52) xi
  !read(53) chi
  read(54) ne_u
  read(55) ne_d
  
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
  print *, "saved to ""wp_eta.dat"", ""wm_eta.dat"" "
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
    end if
  end do

  if(surface_chi_zero) then
    print *, "set all winding numbers of chi on the surface to zero"
    w_s(:) = 0
  end if
  
  !w_s = -w_s
  !w_s(1) = -w_s(1)
  !w_s(2) = -w_s(2)
  !w_s(3) = -w_s(3)
  !w_s(4) = -w_s(4)
  !w_s(5) = -w_s(5)
  !w_s(8) = -w_s(8)
  !w_s(2) = 1
  !w_s(3) = -1
  !w_s(7) = -1
  !w_s(8) = 1
  !w_s(5) = -1
  !w_s(8) = 0
  !w_s(5) = 0
  !w_s(10) = 1
  !w_s(11) = -1
  !w_s(3) = -1
  
  !do i = 1, size(pole_s(1, :))-1
  !  do j = i+1, size(pole_s(1, :))
  !   pole(1:2) = pole_s(1:2, i) - pole_s(1:2, j)
  !    if(pole(1)**2 + pole(2)**2 < 1.1**2) then
  !      ! near pole
  !      pole_s(:, i) = pole_s(:, i) + 0.01d0*pole(1:2)
  !      pole_s(:, j) = pole_s(:, j) - 0.01d0*pole(1:2)
  !    end if
  !  end do
  !end do
  
  call write_format_of_chi_sb('chi_sb_init.plt', unit_num, pshape, hole_chi, pole_s, w_s)

  !w_s(1) = -w_s(1)
  !print *, pole_s
  !print *, w_s
  ! set init chi
  allocate(an_chi)
  an_chi = new_angular_variable_3D(pshape, hole_chi)
  chi0 = an_chi%get_angular_variable() ! bulk
  chi0 = chi0 + an_chi%get_atan_phase_on_xy(pole_s, w_s, 1) ! surface
  !ph = [5.5d0,4d0,1.5d0]
  !do i=1,n
  !  r = site_to_r(i,pshape)
  !  if (r(2) == 4) then
  !    !chi0(i) = chi0(i) + atan2(10*(ph(3)-r(3)),(ph(1)-r(1)))
  !  end if
  !end do
  !ph = [6d0,6.5d0,1.5d0]
  !do i=1,n
  !  r = site_to_r(i,pshape)
  !  if (r(1) == 6) then
  !    !chi0(i) = chi0(i) + 1.d0*atan2((ph(3)-r(3)),(ph(2)-r(2)))
  !  end if
  !end do
  !!chi0 = eta
  chi0 = an_chi%rebuilt_angval_chi(chi0,0.0d0*pi)

  CALL draw_phase_3D('init_chi.dat', unit_num, pshape, chi0, hole_chi(:, 1))
  
  call get_winding_number(chi0, clist, wlist)
  call draw_winding_circle("wp_chi_init.dat", "wm_chi_init.dat", 50, 51, pshape, clist%cil, clist%cpl, wlist)
  
  open(21, file = "beg.dat")

  do b_iter = 0, 10000
    Bz = real(b_iter, 8)/10000d0 * 1d0
    !ALLOCATE(fopt)
    ALLOCATE(fopta)
    !fopt = new_fchi_optimizer(pshape, hole_chi, ne,  t, u, jd, lm, xi, chi0, wf, eg, ne_u, ne_d) !slanting path
    fopta = new_fchi_optimizer_Aem(pshape, hole_chi, ne,  t, u, jd, lm, Bz, xi, chi0, wf, eg, ne_u, ne_d) !slanting path
    !!call fopt%set_jex(jex)
    !call fopta%set_jex(jex)

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

  ! ALLOCATE(an_xi)
  ! an_xi = new_angular_variable_3D(pshape,hole_xi)
  ! open(120, file = "xi1_eg.dat")
  ! do i = 1, 100
  !   xi = an_xi%rebuilt_angval_xi(xi, 2d0*pi*(i-1d0)/100d0 )
  !   write(120, *) xi(1) , total_energy_sb_2lay_xichi_T0(fopt, wf, xi, chi)
  ! end do
  ! close(120)
  ! DEALLOCATE(an_xi)
  
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!! optimize xi1         !!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
    write(*, *) '++++++Optimize xi1++++++++++'
    !old_E_tot = total_energy_sb_2lay_xichi_T0(fopt, wf, xi, chi)
    old_E_tot = total_energy_sb_2lay_xichi_T0_Aem(fopta, wf, xi, chi, fopta%Aem_dphase)
    write(*, *) 'init E_tot', old_E_tot
    if (lm /= 0.d0) then
      ALLOCATE(an_xi)
      an_xi = new_angular_variable_3D (pshape,hole_xi)
      cgr = new_conjugate_gradient_rsoc(n, n_acc, lm, an_xi, wf, xi, chi)
      p(1) = 0.6d0*pi ! initial xi1
      call cgr%calc(p,  optimized_cgr)
      p(1) = modulo(p(1), 2.d0*pi)
      xi = an_xi%rebuilt_angval_xi(xi, p(1))
      write(*, *) 'optimized xi1', xi(1)
    !new_E_tot = total_energy_sb_2lay_xichi_T0(fopt, wf, xi, chi)
    new_E_tot = total_energy_sb_2lay_xichi_T0_Aem(fopta, wf, xi, chi, fopta%Aem_dphase)
    write(*, *) 'E_tot', new_E_tot
    write(*, *) '++++++++end xi1 optimization+++++++++++++++'

      DEALLOCATE(an_xi)
      call fopta%modify_init(hole_chi, xi, chi0, wf)
      if (use_self_consistent) then
        CALL fopta%self_consistent()
      else
        CALL fopta%oneshot()
      end if
      chi = an_chi%get_rebuilt_phase(fopta%delta_chi)

    endif

    if(b_iter == 0) then 
      ALLOCATE(tdelta_for(fopta%pl_1st_l(1)%length()))
      ALLOCATE(tdelta_rev(fopta%pl_1st_l(1)%length()))
      ALLOCATE(Sx(n))
      ALLOCATE(Sy(n))
      ALLOCATE(Sz(n))
      ALLOCATE(tne_u(n))
      ALLOCATE(tne_d(n))
      !read(56) delta
      !read(57) Sx, Sy, Sz
      !read(58) Sy
      !read(59) Sz
      !mu(1) = -0.2d0
      !mu(2) = 4d0
      read(58) mu
      !mu(:) = mu(:) - 0.1
      ! print *, (abs(delta))
    end if
    
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
    new_E_tot = total_energy_sb_2lay_xichi_T0_Aem(fopta, wf, xi, chi, fopta%Aem_dphase)
    print *, 'E_tot(T=0)', new_E_tot
    write(21, *) Bz, new_E_tot

    deallocate(dg)
    deallocate(fopta)
  end do ! b_iter
   
  stop

  !CALL draw_phase_3D('fchi.dat', unit_num, pshape, chi, hole_chi(:, 1))
  
  !!call fopt%get_Jp_matrix(current)
  !call fopta%get_Jp_matrix(current)
  !!call draw_current("current.dat", 101, pshape, fopt%cpl, current)
  !call draw_current("current.dat", 101, pshape, fopta%cpl, current)
  
  !call get_winding_number(chi, clist, wlist)
  !call draw_winding_circle("wp_chi.dat", "wm_chi.dat", 50, 51, pshape, clist%cil, clist%cpl, wlist)
  !print *, "saved to ""wp_chi.dat"", ""wm_chi.dat"" "
  

  !!call get_winding_point_surface(chi, clist, pole_s_conv, w_s_conv)
  !call get_chi_winding_number_with_eta_winding_point_surface(eta, chi, clist, pole_s_conv, w_s_conv)
  !call write_format_of_chi_sb('chi_sb_conv.plt', unit_num, pshape, hole_chi, pole_s_conv, w_s_conv)
  
  !!! current on only bulk
  !!current(1:fopt%nl(1),:) = 0d0
  !!current(:,1:fopt%nl(1)) = 0d0
  !!call draw_current("current_b.dat", 101, pshape, fopt%cpl, current)
  !!!print *, current

  !!!! out jump
  !!open(108, file = "jump.dat") 
  !!do i = 1, fopt%np ! - 4*this%nh
  !!  k = fopt%cpl%value(i)%f
  !!  j = fopt%cpl%value(i)%i
  !!  if(fopt%jump(i) /= 0d0) then
  !!    write(108,"(6f)") (site_to_r(j, pshape) + site_to_r(k, pshape) )*0.5d0 &
  !!                    &   - (site_to_r(k, pshape) - site_to_r(j, pshape) )*0.8d0*0.5d0, &
  !!                    & (site_to_r(k, pshape) - site_to_r(j, pshape) )*0.8d0
  !!  end if
  !!end do
  !!close(108)
  !!WRITE(*,*)'end'
  

  !!!!!!!  calc 1-electron-energy  !!!!!!!
  !ALLOCATE(delta(fopt%pl_1st_l(1)%length()))
  !ALLOCATE(tdelta_for(fopt%pl_1st_l(1)%length()))
  !ALLOCATE(tdelta_rev(fopt%pl_1st_l(1)%length()))
  ALLOCATE(tdelta_for(fopta%pl_1st_l(1)%length()))
  ALLOCATE(tdelta_rev(fopta%pl_1st_l(1)%length()))
  ALLOCATE(Sx(n))
  ALLOCATE(Sy(n))
  ALLOCATE(Sz(n))
  ALLOCATE(tne_u(n))
  ALLOCATE(tne_d(n))
  !read(56) delta
  !read(57) Sx, Sy, Sz
  !read(58) Sy
  !read(59) Sz
  !mu(1) = -0.2d0
  !mu(2) = 4d0
  read(58) mu
  !mu(:) = mu(:) - 0.1
  ! print *, (abs(delta))
  allocate(dg)
  dg = new_diag(pshape, hole_xi, t, u, jd, lm)
  
  
  dg%wf = wf
  dg%eg = eg
  call dg%calc_tdelta_xi(tdelta_for, tdelta_rev, xi, 0.01d0)
  !tdelta_for = delta
  !tdelta_rev = delta
  !do i = 1 dg%pl_1st_l(1)%length()
  !
  !end do
  !call dg%calc_ne(tne_u, tne_d, 0.01d0)
  !ALLOCATE(an_xi)
  !an_xi = new_angular_variable_3D(pshape,hole_xi)
  !xi = an_xi%rebuilt_angval_xi(xi, 0d0)
  call dg%calc_spin_xi(Sx, Sy, Sz, xi, 0.01d0)
  !print *, "diag, bulk ne : ", sum(tne_u(dg%nl(1)+1:) + tne_d(dg%nl(1)+1:))
  !print *, "diag, surface ne : ", sum(tne_u(1:dg%nl(1)) + tne_d(1:dg%nl(1)))

  !ne_u = tne_u
  !ne_d = tne_d

  !chi = 0d0
  !delta = 0d0
  !chi(1:dg%nl(1)) = 0d0
  !do i = 1, dg%pl_1st_l(1)%length()
  !  js = dg%pl_1st_l(1)%value(i)%i
  !  ks = dg%pl_1st_l(1)%value(i)%f
  !  j = dg%ste(js)
  !  k = dg%ste(ks)
  !  delta = exp(0.5d0*ui*(chi(js)+chi(ks)))*delta
  !end do
  !dg%lm = 0d0
  !dg%lm = 0d0

  !call calc_eigen_energy_sb_2lay(dg, eg, wf, xi, chi, mu, tdelta_for, tdelta_rev, ne_u, ne_d, Sx, Sy, Sz)
  call calc_eigen_energy_sb_2lay_Aem(dg, eg, wf, xi, chi, mu, fopta%Aem_dphase, tdelta_for, tdelta_rev, ne_u, ne_d, Sx, Sy, Sz)

  !call calc_eigen_energy_sb_2lay_xichi(dg, eg, wf, mu, delta, ne_u, ne_d, Sx, Sy, Sz, xi, chi)
  !!call calc_eigen_energy_sb_2lay_xi(dg, eg, wf, mu, delta, ne_u, ne_d, Sx, Sy, Sz, xi, chi)
  !call dg%calc_eg_sb_2lay_xi(mu, delta, ne_u, ne_d, Sx, Sy, Sz, xi)
  !eg = dg%eg
  !wf = dg%wf
  open(105,file = "fchi_eg0.dat")
  do i = 1, size(eg)
      write(105, *) i, eg(i)
  end do
  close(105)
  call sort_wf(eg, wf)
  open(106,file = "fchi_eg.dat")
  do i = 1, size(eg)
      write(106, *) i, eg(i)
  end do
  close(106)
  print *, "saved to ""fchi_eg.dat"" "
  
  !call calc_eigen_energy_sb_2lay_wf_T0(dg, eg, wf, xi, chi, mu)
  !open(105,file = "fchi_eg2.dat")
  !do i = 1, size(eg)
  !    write(105, *) i, eg(i)
  !end do
  !close(105)
  !stop

  dg%wf = wf
  dg%eg = eg
  call dg%calc_ne(ne_u, ne_d, kbt)
  !call dg%calc_ne(tne_u, tne_d, 0.01d0)
  !call dg%calc_ne(tne_u, tne_d, 0.0000001d0)
  !call dg%calc_spin(Sx, Sy, Sz)
  print *, "fchi, bulk ne : ", sum(ne_u(dg%nl(1)+1:) + ne_d(dg%nl(1)+1:))
  print *, "fchi, surface ne : ", sum(ne_u(1:dg%nl(1)) + ne_d(1:dg%nl(1)))
  print "(a10,f6.3,a10,f10.3)", "mu surf : ", mu(1), "   bulk : ", mu(2)
  !print *, tne_u + tne_d
  !print *,  tne_u(dg%nl(1)+1:) + tne_d(dg%nl(1)+1:)
  !print *, 'E_tot(T=0)', total_energy_sb_2lay_xichi_T0(dg, wf, xi, chi)
  print *, 'E_tot(T=0)', total_energy_sb_2lay_xichi_T0_Aem(dg, wf, xi, chi, fopta%Aem_dphase)
  !stop
  
  
  call dg%calc_spin_xi(Sx, Sy, Sz, xi, kbt)
  call draw_spin("spin_fchi.dat", 114, pshape, Sx, Sy, Sz)
  call draw_spin_normalize("spin_n_fchi.dat", 115, pshape, Sx, Sy, Sz)
  print *, "saved to ""spin_fchi.dat"", ""spin_n_fchi.dat"" "

  
  !! for monte-carlo
  !!write(60) wf
  !!write(61) eg
  !!write(62) xi
  !!write(63) chi
  !!write(64) ne_u
  !!write(65) ne_d
  !!write(66) tdelta_for, tdelta_rev
  !!write(67) Sx, Sy, Sz
  !!write(68) mu

  ! for monte-carlo
  write(70) wf
  write(71) eg
  write(72) xi
  write(73) chi
  write(74) ne_u
  write(75) ne_d
  write(76) pole_s
  write(77) w_s
  !write(77) w_s_conv
  write(78) size(pole_s(1,:))  ! n_pole
  write(79) mu
  
  ALLOCATE(wf1(4*n_acc, 4*n_acc))
  ALLOCATE(wf2(4*n_acc, 4*n_acc))
  do i = 1, dg%n_acc
    j = dg%ets(i)
    wf1(4*i-3, :) = exp(-0.5d0*ui*chi(j))*exp(-0.5d0*ui*xi(j))*wf(4*i-3, :)
    wf1(4*i-2, :) = exp(-0.5d0*ui*chi(j))*exp( 0.5d0*ui*xi(j))*wf(4*i-2, :)
    wf1(4*i-1, :) = exp( 0.5d0*ui*chi(j))*exp( 0.5d0*ui*xi(j))*wf(4*i-1, :)
    wf1(4*i  , :) = exp( 0.5d0*ui*chi(j))*exp(-0.5d0*ui*xi(j))*wf(4*i  , :)
  
    wf2(4*i-3, :) = exp(-0.5d0*ui*xi(j))*wf(4*i-3, :)
    wf2(4*i-2, :) = exp( 0.5d0*ui*xi(j))*wf(4*i-2, :)
    wf2(4*i-1, :) = exp( 0.5d0*ui*xi(j))*wf(4*i-1, :)
    wf2(4*i  , :) = exp(-0.5d0*ui*xi(j))*wf(4*i  , :)
  end do

 ! for stm, arpes, etc.
  write(120) eg
  write(121) wf1 ! chiの因子を含む波動関数
  write(122) wf2 ! chiの因子を含まない波動関数(xiは含む)
  !write(123) wf ! chiとxiの因子を含まない波動関数

  !chi = 0d0
  !ne_u = 0.1d0
  !ne_d = 0.1d0
  !ne_u = 0.5d0
  !ne_d = 0.5d0
  !ne_u(1:dg%nl(1)) = 0.5d0*dfloat(dg%nl(1)-dg%nh)/dfloat(dg%nl(1))
  !ne_d(1:dg%nl(1)) = 0.5d0*dfloat(dg%nl(1)-dg%nh)/dfloat(dg%nl(1))
  !delta(1:size(delta)/2) = 0.2d0 ! x-dir
  !delta(size(delta)/2+1:size(delta)) = -0.2d0 ! y-dir
  
  !
  !call dg%diag_sb_2lay_xichi(mu, delta, ne_u, ne_d, Sx, Sy, Sz, xi, chi)
  !open(106,file = "fchi_eg2.dat")
  !do i = 1, size(dg%eg)
  !    write(106, *) i, dg%eg(i)
  !end do
  !close(106)
  !print *, "saved to ""fchi_eg2.dat"" "
  !call dg%calc_ne(ne_u, ne_d, 0.01d0)
  !!call dg%calc_spin(Sx, Sy, Sz)
  !print *, "re-diag, bulk ne : ", sum(ne_u(dg%nl(1)+1:) + ne_d(dg%nl(1)+1:))
  !print *, "re-diag, surface ne : ", sum(ne_u(1:dg%nl(1)) + ne_d(1:dg%nl(1)))
  
  print *, "Bz = ", Bz

  stop "stop : mwf2"

  !!!!!!!!!!!!!!!!    stm    !!!!!!!!!!!!!!!
  !n_step = 5000
  !eg_range(1) = -4d0
  !eg_range(2) =  4d0
  !wf = dg%wf
  !eg = dg%eg
  !do i = 1, dg%n_acc
  !  j = dg%ets(i)
  !  wf(4*i-3, :) = exp(-0.5d0*ui*chi(j))*exp( 0.5d0*ui*xi(j))*wf(4*i-3, :)
  !  wf(4*i-2, :) = exp(-0.5d0*ui*chi(j))*exp(-0.5d0*ui*xi(j))*wf(4*i-2, :)
  !  wf(4*i-1, :) = exp( 0.5d0*ui*chi(j))*exp(-0.5d0*ui*xi(j))*wf(4*i-1, :)
  !  wf(4*i  , :) = exp( 0.5d0*ui*chi(j))*exp( 0.5d0*ui*xi(j))*wf(4*i  , :)
  !end do
  !do i = 1, 10
  !  kbt(i) = 0.01d0*i
  !end do
  !open(102,file = "rho.dat")
  !do i = 1, n_step
  !  e = eg_range(1) + (eg_range(2)-eg_range(1))/float(n_step)*i
  !  rho = 0d0
  !  do j = 2*dg%n_acc+1,4*dg%n_acc
  !    site = dg%ste(xyz_to_site(4,4,1,pshape))
  !    do k = 1, 10
  !      rho(k) = rho(k) - conjg(wf(4*site-3, j))*wf(4*site-3, j)*dfermi(eg(j) - e, kbt(k))
  !      rho(k) = rho(k) - conjg(wf(4*site-2, j))*wf(4*site-2, j)*dfermi(eg(j) - e, kbt(k))
  !      rho(k) = rho(k) - conjg(wf(4*site-1, j))*wf(4*site-1, j)*dfermi(eg(j) + e, kbt(k))
  !      rho(k) = rho(k) - conjg(wf(4*site  , j))*wf(4*site  , j)*dfermi(eg(j) + e, kbt(k))
  !    end do
  !  end do
  !  write(102,"(11f)") e, rho
  !end do
  !close(102)
  !print *,"saved to ""rho.dat"""

END PROGRAM make_wave_function2_be