PROGRAM main
  USE energy_mod
  USE math_mod,              only : pi
  USE io_mod,                ONLY : save_txt, load_txt, create_vortex, load_conf_3D, load_jexconf_3D
  USE visual_mod
  USE angular_variable_mod,  ONLY : angular_variable, new_angular_variable, new_angular_variable_3D
  !use arpes_mod
  !USE car_parrinello_mod,    ONLY : car_parrinello,new_car_parrinello
  !USE car_parrinello_wf_mod, ONLY : car_parrinello_wf,new_car_parrinello_wf
  USE conjugate_gradient_rsoc_mod
  !USE conjugate_gradient_fchi_mod 
  !  use sdchi_optimizer_mod,    only : sdchi_optimizer, new_sdchi_optimizer 
  USE fchi_optimizer_mod,    ONLY : fchi_optimizer, new_fchi_optimizer
  USE basic_current_mod,     ONLY : basic_current, new_basic_current
  !USE single_wf_mod,         ONLY : single_wf, new_single_wf
  !USE site_list_mod
  !USE monte_carlo_zeta_mod !!
  !USE wave_function_mod,      ONLY : wave_function, new_wave_function
  ! use conjugate_gradient_rsoc_mod, only : conjugate_gradient_rsoc,new_conjugate_gradient_rsoc
  use parameter_mod
  use circle_list_mod
  
  IMPLICIT NONE
  INTEGER                                  :: nx, ny, n, ne, nh, iter, nz
  INTEGER, DIMENSION(:,:) , allocatable                   :: pshape
  REAL(8), DIMENSION(3), PARAMETER         :: t = [1.d0, 0d0, 0.1d0]  ! hopping coeff
  REAL(8), PARAMETER                       :: u = 8.d0         ! coulomb coeff
  !  real(8), parameter                       :: jd = 0.d0 ! across hole coeff
  REAL(8), PARAMETER                       :: jd = 0.5d0*4.0d0*(t(1)**2)/u ! across hole coeff
  REAL(8), PARAMETER                       :: lm = 0.02d0 ! rashba coeff
  LOGICAL,PARAMETER                        :: sflag=.TRUE.!slanting path
  !  logical,parameter                        :: sflag=.false.!no slanting path
  logical                                  :: optimized_cgr, optimized_E_tot, dummy_flag
  INTEGER, PARAMETER                       :: unit_num = 10
  INTEGER, DIMENSION(:, :), ALLOCATABLE    :: hole_xi, hole_chi
  REAL(8), DIMENSION(:), ALLOCATABLE       :: xi, chi, ochi, oxi, zeta, chi0
  TYPE(angular_variable), ALLOCATABLE      :: an_xi,an_chi
  !TYPE(basic_current), ALLOCATABLE         :: bc
  !TYPE(car_parrinello_wf),ALLOCATABLE      :: cpw
  !TYPE(car_parrinello),ALLOCATABLE         :: cp
  !  type(sdchi_optimizer),allocatable         :: opt
  TYPE(fchi_optimizer),ALLOCATABLE         :: fopt, fopt_t
  !TYPE(fchi_optimizer_ex),ALLOCATABLE         :: fopt_ex
  !TYPE(single_wf), ALLOCATABLE             :: sw
  !TYPE(site_list), ALLOCATABLE             :: sl
  !TYPE(wave_function), ALLOCATABLE         :: wf
  !type(type_j_ext), allocatable            :: j_ext(:)
  !TYPE(monte_carlo_zeta), ALLOCATABLE      :: mc
  !type(conjugate_gradient_fchi),allocatable:: cgf
  type(conjugate_gradient_rsoc),allocatable:: cgr
  COMPLEX(8), DIMENSION(:, :), ALLOCATABLE :: cpwf, wf1
  real(8), DIMENSION(:, :), ALLOCATABLE    :: rs_current
  REAL(8), DIMENSION(:), ALLOCATABLE       :: cpeg, ranx
  REAL(8)                                  :: p(1)
  INTEGER                                  :: i, ii, nrand, j
  REAL(8), DIMENSION(:), ALLOCATABLE       :: delta_xi, J1, J2, J3
  REAL(8), DIMENSION(:), ALLOCATABLE       :: mcint
  INTEGER, DIMENSION(:), ALLOCATABLE       :: level, seed
  real(8)                                  :: E_total, E_total0, E_hopping, E_hopping0, &
                                            & E_rashba,E_rashba0, total_E2, xi1, old_E_tot, new_E_tot
  real(8), dimension(:), allocatable       :: E_total_list, E_hopping_list, E_rashba_list
  character(len=1)                         :: cha
  real(8)                                 :: j_min, E_min
  real(8),dimension(:),allocatable        :: jex, jex_0 !jex_i,jex_f
  real(8),dimension(:),allocatable        :: jex_list
  integer       :: jex_site
  !real(8),dimension(:),allocatable        :: jex_start, jex_end
  integer                     :: n_jex, Ei

  !! load xi winding number
  !CALL load_conf2('../config/vconf2', unit_num, pshape, hole_xi, hole_chi, j_ext)
  !deallocate(hole_chi)
  !CALL load_conf_chi_w('../config/vconf2', unit_num, pshape, hole_chi, 7)
  !! load chi winding number
  
  ! load xi and chi winding number
  !CALL load_conf_3D('../config/vconf', unit_num, pshape, hole_xi,hole_chi)
  CALL load_conf_3D('./config/vconf', unit_num, pshape, hole_xi,hole_chi)
  ! load jex
  CALL load_jexconf_3D('./config/jexconf', unit_num, pshape, jex)
  !CALL load_jexconf_3D('../config/jexconf', unit_num, pshape, jex_start, jex_end, n_jex)
  n_jex = 50

  nz = size(pshape(1,:))
  n = site_max(pshape)
  ne = n - size(hole_xi(:,1))

  allocate(E_total_list(n_jex+1))
  allocate(E_hopping_list(n_jex+1))
  allocate(E_rashba_list(n_jex+1))
  allocate(jex_list(n_jex+1))

  allocate(xi(n),source=0.d0)
  allocate(chi(n),source=0.d0)
  allocate(oxi(n),source=0.d0)
  allocate(ochi(n),source=0.d0)
  allocate(jex_0(n),source=jex)
  
  allocate(cpwf(2*n, 2*n))
  allocate(wf1(2*n, 2*n))
  
  read(fchi_xi_trans_unit) oxi
  read(chi_trans_unit) ochi
  read(wavefunction_trans_unit) cpwf
  !read(enegy_trans_unit) cpeg

  !E_total0 = energy_total_energy(pshape, hole_chi, t, u, jd, lm, &
  !& oxi, ochi, cpwf, 'T')
  !wf1 = energy_set_wf(n, hole_chi, cpwf, oxi, ochi)
  !E_rashba0 = energy_spin_orbit_interaction_energy(pshape, n, ne, hole_xi, lm, &
  !& wf1) 
  !E_hopping0 = energy_hopping_energy(pshape, n, ne, hole_xi, t, &
  !& wf1)  

    !ALLOCATE(an_xi) !!!!!!!
    !an_xi = new_angular_variable (pshape,hole_xi)
      !xi = an_xi%rebuilt_angval_xi(xi, 6.d0)
      !deallocate(an_xi)
  ALLOCATE(an_chi)
  an_chi = new_angular_variable_3D(pshape,hole_chi)
  ALLOCATE(an_xi)
  an_xi = new_angular_variable_3D(pshape,hole_xi)
  
  ALLOCATE(fopt)
  ALLOCATE(fopt_t)
  !!   if(sflag)then
  fopt = new_fchi_optimizer(pshape, hole_chi, t, u, jd, lm, oxi, ochi, cpwf) !slanting path
  call fopt%set_param(ratio = 0.1d0, tol = 1d-2)
  !call fopt%set_jex(jex)
  CALL fopt%self_consistent(dummy_flag)

  open(88,file = "jex_eg.dat")
  do i=1, n
    !if(jex_start(i) /= 0d0) then
    if(jex(i) /= 0d0) then
      jex_site = i
      exit
    end if
  end do

  Ei = 0
  do j= -n_jex/2, n_jex/2
    Ei = Ei + 1
    !jex = jex_start + (jex_end - jex_start) * j / real(n_jex)
    !print *, jex_start(jex_site) + (jex_end(jex_site) - jex_start(jex_site)) * j / real(n_jex), j, "/", n_jex
    jex = jex_0 * j / n_jex * 2d0
    !print *, j, "/", n_jex
    print *,"jex = ",  jex(jex_site)
    fopt_t = fopt
    chi = ochi
    xi = oxi
    call fopt_t%set_jex(jex)
    CALL fopt_t%self_consistent(dummy_flag)
    !E_total = energy_total_energy(pshape, hole_chi, t, u, jd, lm, &
    !& xi, chi, cpwf, 'T')
    !write(*, *) 'total_E', E_total
    !stop
    if (lm /= 0.d0 .and. .false. ) then ! calc xi1
      old_E_tot = energy_total_energy(pshape, hole_chi, t, u, jd, lm, &
          & xi, chi, cpwf, 'T', .false.)
      optimized_cgr = .false.
      optimized_E_tot = .false.
      allocate(cgr)
      cgr = new_conjugate_gradient_rsoc(n, lm, an_xi, cpwf, xi, chi)
      p(1) = 0.6d0*pi ! initial xi1
      cgr_loop: do  i = 1, 1000 ! Iterlation
        !write(*, *) 'iter:', i
        cgr = new_conjugate_gradient_rsoc(n, lm, an_xi, cpwf, xi, chi)
        !write(*,*) p(1)
        !write(*, *) 'E_tot before cgr', energy_total_energy(pshape, hole_chi, t, u, jd, lm, &
        !& xi, chi, cpwf, 'T', .false.)
        call cgr%calc(p,  optimized_cgr)
        !p(1) = modulo(p(1),2*pi)
        xi = an_xi%rebuilt_angval_xi(xi, p(1))
        !write(*, *) 'E_tot after cgr', energy_total_energy(pshape, hole_chi, t, u, jd, lm, &
        !& xi, chi, cpwf, 'T',.false.)
        write(*, *) 'xi1', xi(1)
        write(*, *) 'chi1',chi(2), optimized_cgr
        !chi = an_chi%get_angular_variable()
        !chi = an_chi%rebuilt_angval_chi(chi)
        chi = ochi
        call fopt_t%modify_init(hole_chi, xi, chi, cpwf)
        CALL fopt_t%self_consistent(dummy_flag)
        !!!CALL fopt%self_consistent(an_chi%d_angval_chi(chi))
        chi = an_chi%get_rebuilt_phase(fopt_t%delta_chi)
        new_E_tot = energy_total_energy(pshape, hole_chi, t, u, jd, lm, &
        & xi, chi, cpwf, 'T', .false.)
        !write(*, *) 'E_tot', new_E_tot
        !write(*, *) 'flag', abs((new_E_tot - old_E_tot) / old_E_tot)
        if(abs((new_E_tot - old_E_tot) / old_E_tot) < 1.d0*10d0**(-14)) optimized_E_tot = .true.
        if(optimized_E_tot) exit
        old_E_tot = new_E_tot
      end do cgr_loop
      deallocate(cgr)

      if(.not.optimized_E_tot) then
        write(*, *) 'fault cgr optimization'
      else
        !write(*, *) '++++++++xi1 optimization succeeded+++++++++++++++'
      endif

      p(1) = modulo(xi(1), 2.d0*pi)
      xi = an_xi%rebuilt_angval_xi(xi, p(1))
      !DEALLOCATE(an_xi)
      call fopt_t%modify_init(hole_chi, xi, chi, cpwf)
      CALL fopt_t%self_consistent(dummy_flag)
    end if !end calc xi1


    !write(*, *) 'E_tot', energy_total_energy(pshape, hole_chi, t, u, jd, lm, &
    !  & xi, chi, cpwf, 'T')

    chi = an_chi%get_rebuilt_phase(fopt_t%delta_chi)
  
    !do i=1, n
    !  wf1(2*i-1,:) = exp(-ui*0.5d0*chi(i))*exp(-ui*0.5d0*xi(i))*cpwf(2*i-1,:)
    !  wf1(2*i,:) = exp(-ui*0.5d0*chi(i))*exp(ui*0.5d0*xi(i))*cpwf(2*i,:)
    !end do
    wf1 = energy_set_wf(n, hole_chi, cpwf, xi, chi)
    E_total = energy_total_energy(pshape, hole_xi, t, u, jd, lm, xi, chi, cpwf, 'T', .false.)
    E_rashba = energy_spin_orbit_interaction_energy(pshape, n, ne, hole_chi, lm, wf1)
    E_hopping = energy_hopping_energy(pshape, n, ne, hole_xi, t, wf1) 
    !E_hopping = E_hopping + energy_hopping_2nd_energy(pshape, n, ne, hole_xi, t, wf1) 
    E_hopping = E_hopping + energy_hopping_z_energy(pshape, n, ne, hole_xi, t, wf1) 
    !write(88, '(4f)') jex_start(jex_site) + (jex_end(jex_site) - jex_start(jex_site)) * j / real(n_jex), &
    !        &  E_total, E_rashba, E_hopping
    !print '(4f)', jex_start(jex_site) + (jex_end(jex_site) - jex_start(jex_site)) * j / real(n_jex), &
    !        &  E_total, E_rashba, E_hopping
    jex_list(Ei) = jex(jex_site)
    E_total_list(Ei) = E_total
    E_rashba_list(Ei) = E_rashba
    E_hopping_list(Ei) = E_hopping
    if(j == 0) then
      E_total0 = E_total
      E_rashba0 = E_rashba
      E_hopping0 = E_hopping
    end if
    !write(88, '(4f)') jex(jex_site), E_total, E_rashba-E_rashba0, E_hopping-E_hopping0
    !print '(4f)', jex(jex_site), E_total, E_rashba-E_rashba0, E_hopping-E_hopping0
    print '(4f)', jex(jex_site), E_total, E_rashba, E_hopping

  end do
  
  do Ei= 1, n_jex+1
    write(88, '(4f)') jex_list(Ei), E_total_list(Ei), E_rashba_list(Ei)-E_rashba0, E_hopping_list(Ei)-E_hopping0
    !print '(4f)', jex_list(Ei), E_total_list(Ei), E_rashba_list(Ei)-E_rashba0, E_hopping_list(Ei)-E_hopping0
  
  end do
  
  
  print *,"finished mwf2 j and E"
  close(88)
END PROGRAM main
