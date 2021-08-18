MODULE monte_carlo_2lay_mod
  USE mpi
  USE angular_variable_mod
  !use energy_mod
  !USE energy_mod,            only : total_energy_sb_2lay_xichi_T0, calc_eigen_energy_sb_2lay_wf_T0
  USE energy_mod,            only : total_energy_sb_2lay_xichi_T0_Aem, calc_eigen_energy_sb_2lay_wf_T0_Aem
  USE math_mod, ONLY:pi, ui, principal, sign2, principal2, principal, &
       &                              inv, matvecmul, avg
  !USE io_mod, ONLY:save_txt
  !USE path_list_mod, ONLY:path_list, path
  !  use dchi_optimizer_mod, only : dchi_optimizer, my_dsysv
  !USE fchi_optimizer_mod
  USE fchi_optimizer_Aem_mod
  USE hop_iterator_mod, ONLY:hop_iterator, new_hop_iterator, get_hopping_path
  !  use monte_carlo_mod,    only : monte_carlo
  USE int_list_mod, ONLY:int_list
  !USE spinwave_mod
  !USE single_wf_mod
  use utility_mod, only : sort_wf
  use winding_number_mod, only : get_eta, get_chi_winding_number_with_eta_winding_point_surface, &
                          &     get_winding_point_surface, get_winding_number
  use visual_mod, only : write_format_of_chi_sb, draw_current, draw_phase_3D, draw_winding_circle
  IMPLICIT NONE


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!         options        !!!!!!!!!!!!
  
  !   stop_found_min_E
  !       when found a state having energy smaller than initial energy,
  ! true  : write/update "config/min_E_state.txt" and stop this program.
  ! false : continue this program. finally, write "min_E_state_($myid).txt".
  !logical, parameter :: stop_found_min_E = .true.
  logical, parameter :: stop_found_min_E = .false.

  !!!!!!!!!!                        !!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  type int_alloc
    integer, allocatable :: v(:)
  end type

  !  type, extends(monte_carlo),public :: monte_carlo_2lay
  TYPE, PUBLIC :: monte_carlo_2lay
    !    type(dchi_optimizer)                      :: opt
    !TYPE(fchi_optimizer), allocatable         :: fopt !!
    TYPE(fchi_optimizer_Aem), allocatable         :: fopt !!
    !TYPE(spinwave), allocatable                :: swave
    TYPE(int_list), DIMENSION(:), ALLOCATABLE   :: near_circle
    !TYPE(single_wf), ALLOCATABLE              :: single_shot
    INTEGER                                   :: myid
    ! INTEGER                                   :: judg,sjudg, swavejudg
    INTEGER                                   :: change_count(6), not_converged_count, mode_count(6) ! count the number of the changing states
    INTEGER                                   :: rcount, nsh, csw_num !, es=0
    INTEGER     :: energydat = 30, specificheatdat = 31, un_loadconf, nm_errorcheck = 32
    INTEGER                                   :: iter, iter_max
    INTEGER                                   :: Titer_e, Titer_ne
    INTEGER                                   :: change_mode
    INTEGER, DIMENSION(4)                     :: cw
    INTEGER, DIMENSION(:), ALLOCATABLE        :: clock, seed
    !INTEGER, DIMENSION(:), ALLOCATABLE        :: wma
    INTEGER, DIMENSION(:), ALLOCATABLE        :: wma
    INTEGER, DIMENSION(:), ALLOCATABLE        :: win, wn
    ! the infomation about which labels of spin-wave excitations are excited
    !INTEGER, DIMENSION(:), ALLOCATABLE        :: swave_info
    !integer, allocatable                      :: hist_xi1(:)  ! histgram of variable xi(1)
    !real(8), allocatable                      :: hist_classval(:)
    !integer                                   :: hist_intervel_number = 100
    REAL(8)                                   :: T_max, T_min, kb = 1.d0, error = 1.d-8
    !REAL(8)                                   :: trial_xi1
    ! sE: sum of energy
    !REAL(8)                                   :: old_E, old_E_sw, old_E_wf, old_E_xi1, teg
    REAL(8)                                   :: old_E
    !REAL(8)                                   :: sE, sE2, sE_sw, sE2_sw, sE_wf, sE2_wf, sE_xi1, sE2_xi1, &
    !  &                                         sxi1, sxi12
    !REAL(8)                                   :: sE, sE2
    REAL(8)                                   :: sum_E, sum_E2 ! E2 means E**2
    !REAL(8)                                   :: w, ws, aw, awc
    REAL(8)                                   :: sum_diff_w, sum_diff_w2, sum_ave_w, sum_ave_w2
    REAL(8)                                   :: sum_diff_w_b, sum_diff_w_b2, sum_ave_w_b, sum_ave_w_b2
    REAL(8)                                   :: sum_diff_w_s, sum_diff_w_s2, sum_ave_w_s, sum_ave_w_s2
    integer                                   :: add_ensemble_count, not_add_ensemble_count
    !REAL(8)                                   :: Mab, Mc, Mhab, Mhc
    !REAL(8)                                   :: Mab2, Mc2, Mhab2, Mhc2
    !REAL(8)                                   :: sMab, sMc, sMhab, sMhc
    !REAL(8)                                   :: sMab2, sMc2, sMhab2, sMhc2
    !REAL(8)                                   :: scrank_level, scrank_level2
    ! counter for statistical parameter___/
    !integer                                   :: SP_wn, SP_sw, SP_xi1
    !REAL(8)                                   :: cwrate = 0.d0, cwsrate = 1d0, &
    !  &                                           cssrate = 0.d0, xirate = 0.3d0, &
    !  &                                           sw_and_wn_rate = 0.d0
    !REAL(8)                                   :: cwrate = 0.5d0, cwsrate = 0.7d0, &
      !&                                           cssrate = 0.d0, xirate = 0.3d0, &
      !&                                           sw_and_wn_rate = 0.d0
    !REAL(8)                                   :: cwrate = 0.5d0, cwsrate = 0.1d0, &
      !&                                           cssrate = 0.d0, xirate = 0.0d0
    !REAL(8)                                   :: change_rate(5) = [0.d0, 1d0, 1d0, 0d0, 1d0]
    !REAL(8)                                   :: change_rate_partitifalse)
    ! variables to store ones accepted chenges in monte-carlo caluculations##
    integer, allocatable                 :: present_hole_chi(:, :), hole_chi_GS(:, :)
    !real(8), allocatable                     :: present_xi(:), present_chi(:), present_zeta(:)
    real(8), allocatable                     :: present_xi(:), present_chi(:), present_eg(:), &
                                             &  present_ne_u(:), present_ne_d(:), present_mu(:)
    REAL(8), ALLOCATABLE                      :: Sx(:), Sy(:), Sz(:)
    complex(8), allocatable                  :: present_wf(:, :)
    !complex(8), allocatable                 :: crank_EV(:)
    !real(8), allocatable                 :: crank_eg(:)
    !integer                                 :: trial_crank_level, present_crank_level, crank_valid_level
    ! #######################################################################
    REAL(8), DIMENSION(:), ALLOCATABLE        :: cseta, cavs
    REAL(8), DIMENSION(:, :), ALLOCATABLE     :: m
    COMPLEX(8), DIMENSION(:, :), ALLOCATABLE     :: Dp, Dm
    !    integer                                   :: un_judg=32, un_E2=33, un_lw=34
    !    & un_nw = 35, un_nm = 36, un_nw2 = 37, un_nm2 = 38
    LOGICAL                                   :: cwf, csf
    !    integer, dimension(:), allocatable        :: wa, wm

    !!! BdG ---additonal---
    integer                                   :: n_pole
    !integer, dimension(:), allocatable        :: wma_s
    integer, dimension(:), allocatable        :: wma0_s
    real(8), dimension(:,:), allocatable      :: present_pole_s
    integer, dimension(:), allocatable        :: present_w_s, w_s
    integer, dimension(:), allocatable        :: wn_s
    type(int_alloc), dimension(:), allocatable:: near_pole
    real(8)                                   :: near_range = 3d0

    ! snap shot
    real(8), allocatable                     :: chi(:), eta(:)
    !integer, allocatable                     :: present_hole_chi(:, :), hole_chi_GS(:, :)

    ! check groud state
    real(8)                                   :: present_E_tot
    real(8)                                   :: min_E_tot
    integer, dimension(:), allocatable        :: min_w_b, min_w_s
    real(8), dimension(:), allocatable        :: min_chi
    
    real(8)                                   :: Bz

  CONTAINS
    ! public type-bound procedures
    PROCEDURE :: monte_carlo_2lay_init
    !    procedure :: get_energy             => monte_carlo_2lay_get_energy
    !    procedure :: get_total_energy       => monte_carlo_2lay_get_total_energy
    PROCEDURE :: set_mcint => monte_carlo_2lay_set_mcint
    PROCEDURE :: result_mc => monte_carlo_2lay_result_mc
    ! protected type-bound procedures
    PROCEDURE :: calc => monte_carlo_2lay_calc
    PROCEDURE :: Metropolis => monte_carlo_2lay_Metropolis
    PROCEDURE :: Judgment => monte_carlo_2lay_Judgment
    PROCEDURE :: set_ran_seed => monte_carlo_2lay_set_ran_seed
    !PROCEDURE :: set_near_circle => monte_carlo_2lay_set_near_circle
    !PROCEDURE :: calc_winding_number => monte_carlo_2lay_calc_winding_number
    !PROCEDURE :: calc_Ms => monte_carlo_2lay_calc_Ms
    !PROCEDURE :: calc_Ms2 => monte_carlo_2lay_calc_Ms2
    PROCEDURE :: calc_init => monte_carlo_2lay_calc_init
    PROCEDURE :: calc_finalize => monte_carlo_2lay_calc_finalize
    !PROCEDURE :: set_change_partition => monte_carlo_2lay_set_change_partition
    PROCEDURE :: change_state => monte_carlo_2lay_change_state
    PROCEDURE :: Metropolis_init => monte_carlo_2lay_Metropolis_init
    PROCEDURE :: Metropolis_finalize => monte_carlo_2lay_Metropolis_finalize
    !PROCEDURE :: xi1_hist_counter => monte_carlo_2lay_xi1_hist_counter
    PROCEDURE :: ensemble_average => monte_carlo_2lay_ensemble_average
    !PROCEDURE :: get_change_w => monte_carlo_2lay_get_change_w
    PROCEDURE :: get_flip_w_b => monte_carlo_2lay_get_flip_w_b
    PROCEDURE :: get_swap_w_b => monte_carlo_2lay_get_swap_w_b
    PROCEDURE :: get_flip_w_s => monte_carlo_2lay_get_flip_w_s
    PROCEDURE :: get_swap_w_s => monte_carlo_2lay_get_swap_w_s
    PROCEDURE :: get_annihilate_w_s => monte_carlo_2lay_get_annihilate_w_s
    PROCEDURE :: get_create_w_s => monte_carlo_2lay_get_create_w_s
    !PROCEDURE :: get_around_hole_pos => monte_carlo_2lay_get_around_hole_pos
    !PROCEDURE :: get_new_wf => monte_carlo_2lay_get_new_wf
    !PROCEDURE :: get_change_spin => monte_carlo_2lay_get_change_spin
    PROCEDURE :: filename_myid => monte_carlo_2lay_filename_myid
    PROCEDURE :: avg_quantity => monte_carlo_2lay_avg_quantity
    !PROCEDURE :: avg_hist => monte_carlo_2lay_avg_hist
    PROCEDURE :: check_min_E => monte_carlo_2lay_check_min_E
    PROCEDURE :: out_min_E_state => monte_carlo_2lay_out_min_E_state
  END TYPE monte_carlo_2lay
CONTAINS

  !SUBROUTINE monte_carlo_zeta_init(this, pshape, hole, t, u, jd, lm, xi, zeta, chi, wf, rs, crank_EV, crank_eg, myid)
  !subroutine monte_carlo_2lay_init(this, pshape, hole, pole_s, w_s, t, u, jd, lm, mu, xi, chi, wf, eg, ne_u, ne_d, myid)
  subroutine monte_carlo_2lay_init(this, pshape, hole, pole_s, w_s, t, u, jd, lm, mu, Bz, xi, chi, wf, eg, ne_u, ne_d, myid)
    CLASS(monte_carlo_2lay), INTENT(inout)                :: this
    !INTEGER, DIMENSION(2), INTENT(in)                   :: pshape
    !INTEGER, DIMENSION(:, :), INTENT(in)                :: hole
    !REAL(8), DIMENSION(:), INTENT(in)                   :: t
    !REAL(8), INTENT(in)                                 :: u, jd, lm !!
    !REAL(8), DIMENSION(pshape(1)*pshape(2)), INTENT(in) :: xi, chi, zeta
    !COMPLEX(8), DIMENSION(2*pshape(1)*pshape(2), 2*pshape(1)*pshape(2)), &
    !     & INTENT(in)                                        :: wf
    !INTEGER, INTENT(in), OPTIONAL                       :: myid
    integer, dimension(:, :), intent(in)                :: pshape
    integer, dimension(:, :), intent(in)                :: hole
    integer, dimension(:), intent(in)                   :: w_s
    integer, intent(in), optional                       :: myid
    real(8), dimension(:, :), intent(in)                :: pole_s
    !real(8), intent(in)                                 :: t(3), lm
    real(8), intent(in)                                 :: t(3), u, jd, lm, mu(2), Bz
    real(8), dimension(:), intent(in)                   :: xi, chi, ne_u, ne_d, eg
    complex(8), dimension(:, :), intent(in)             :: wf
    INTEGER, DIMENSION(:), ALLOCATABLE               :: seed
    integer                                 :: nh
    !REAL(8), DIMENSION(pshape(1)*pshape(2)) :: trial_chi
    REAL(8), DIMENSION(size(chi)) :: trial_chi
    type(angular_variable), allocatable                          :: an_chi
    logical                                         :: optimized_fchi
    INTEGER                                          :: nrand
    !complex(8)                                      :: crank_EV(:)
    !real(8)                                      :: crank_eg(:)
    CHARACTER(1)                                     :: rs !!
    nh = size(hole)
    IF (PRESENT(myid)) this%myid = myid
    write (*, *) 'myid', myid
    CALL RANDOM_SEED(size=nrand)
    ALLOCATE (this%seed(nrand))
    IF (.NOT. PRESENT(myid)) CALL this%set_ran_seed()
    !    call this%monte_carlo_init()
    allocate(this%fopt)
    allocate(an_chi)
    !an_chi = new_angular_variable(pshape, hole)
    an_chi = new_angular_variable_3D(pshape, hole)
    !allocate(this%present_chi(this%fopt%n), source=chi)
    !this%present_chi = an_chi%get_angular_variable()
    !this%present_chi = an_chi%rebuilt_angval_chi(trial_chi)
    !this%fopt = new_fchi_optimizer(pshape, hole, t, u, jd, lm, xi, chi, wf, rs)
    !this%fopt = new_fchi_optimizer(pshape, hole, ne,  t, u, jd, lm, xi, chi, wf, eg, ne_u, ne_d)
    !this%fopt = new_fchi_optimizer(pshape, hole, [0d0, 0d0],  t, 0d0, 0d0, lm, xi, chi, wf, eg, ne_u, ne_d) ! ne, u and jd are not used in fchi
    !this%fopt = new_fchi_optimizer(pshape, hole, [0d0, 0d0],  t, u, jd, lm, xi, chi, wf, eg, ne_u, ne_d) ! ne is not used in fchi(monte-carlo)
    this%fopt = new_fchi_optimizer_Aem(pshape, hole, [0d0, 0d0],  t, u, jd, lm, Bz, xi, chi, wf, eg, ne_u, ne_d) ! ne is not used in fchi(monte-carlo)
    call this%fopt%set_param(pflag = .false.)
    !call this%fopt%self_consistent2(crank_EV(1), optimized_fchi)
    !this%present_chi = an_chi%get_rebuilt_phase(this%fopt%delta_chi)
    allocate(this%present_hole_chi(nh, 2))
    this%present_hole_chi = hole
   ! allocate(this%hole_chi_GS(nh, 2))
   ! this%hole_chi_GS = hole
    allocate(this%present_xi(this%fopt%n), source=xi)
    allocate(this%present_chi(size(chi)), source=chi)
    !!  write(6,*)myid, 3
    !ALLOCATE (this%swave_info(this%fopt%n), source=0)
    !allocate(this%swave)
    !this%swave = new_spinwave(pshape, hole, u, xi, zeta)
    !allocate(this%crank_EV(size(crank_EV, 1)), source=crank_EV)
    !allocate(this%crank_eg(size(crank_eg, 1)), source=crank_eg)
    this%n_pole = size(w_s)
    allocate(this%present_pole_s(2, this%n_pole))
    this%present_pole_s = pole_s
    allocate(this%present_w_s(this%n_pole))
    allocate(this%w_s(this%n_pole))
    this%present_w_s = w_s
    this%w_s = w_s
    allocate(this%present_wf(size(wf,1), size(wf,2)))
    this%present_wf = wf
    allocate(this%present_eg(size(eg)))
    this%present_eg = eg
    allocate(this%present_ne_u(size(ne_u)))
    this%present_ne_u = ne_u
    allocate(this%present_ne_d(size(ne_d)))
    this%present_ne_d = ne_d
    allocate(this%present_mu(size(mu)))
    this%present_mu = mu
    this%Bz = Bz
    
    deallocate(an_chi)
    
    allocate(this%chi(size(chi)), source=chi)
    allocate(this%eta(size(xi)))
    call get_eta(this%eta, xi, pshape)

    call this%fopt%calc_only_Aem_dphase()
    !this%min_E_tot = total_energy_sb_2lay_xichi_T0(this%fopt, this%present_wf, this%present_xi, this%present_chi)
    this%min_E_tot = total_energy_sb_2lay_xichi_T0_Aem(this%fopt, this%present_wf, this%present_xi, this%present_chi,&
                          & this%fopt%Aem_dphase)
    this%present_E_tot = this%min_E_tot
    allocate(this%min_w_b(this%fopt%nh), source = hole(:,2))
    allocate(this%min_w_s(this%n_pole), source = w_s)
    allocate(this%min_chi(size(chi)), source = chi)
  END SUBROUTINE monte_carlo_2lay_init

  SUBROUTINE monte_carlo_2lay_set_mcint(this, mcint)
    CLASS(monte_carlo_2lay), INTENT(inout)  :: this
    REAL(8), DIMENSION(6), INTENT(in)       :: mcint
    this%T_min = mcint(1)
    this%T_max = mcint(2)
    this%iter_max = INT(mcint(3))
    this%Titer_ne = INT(mcint(4))
    this%Titer_e = INT(mcint(5))
    this%nsh = INT(mcint(6))
  END SUBROUTINE monte_carlo_2lay_set_mcint

  SUBROUTINE monte_carlo_2lay_calc(this, RorD)
    CLASS(monte_carlo_2lay), INTENT(inout)    :: this
    CHARACTER, INTENT(in)                :: RorD ! r is rise, u is descent
    INTEGER                              :: iter
    REAL(8)                              :: T, dT
    CALL this%calc_init(RorD, T, dT)
    DO iter = 0, this%iter_max
      this%iter = iter
      T = dT + T
      WRITE (*, *) ''
      !WRITE (*, '(a6,i5,a2,i6,a2)') '-iter', iter, ' /', this%iter_max, ' -'
      WRITE (*, '(a,i5,a,i6,a,i3,a)') '-iter', iter, ' /', this%iter_max, ' (myid:',this%myid,') -'
      WRITE (*, *) 'T', T
      CALL this%Metropolis_init(iter)
      !stop 'a'
      CALL this%Metropolis(T, .FALSE., 't') ! thermodynamic non-equilibrium
      CALL this%Metropolis(T, .TRUE., 't') ! thermodynamic equilibrium
      CALL this%Metropolis_finalize(T)
      !write(*, *) this%hist_xi1
    END DO
    CALL this%calc_finalize()
  END SUBROUTINE monte_carlo_2lay_calc

  SUBROUTINE monte_carlo_2lay_Metropolis(this, T, thermal_equil, rs)
    CLASS(monte_carlo_2lay), INTENT(inout)    :: this
    REAL(8), INTENT(in)                  :: T
    LOGICAL, INTENT(in)                  :: thermal_equil ! after no-mesurement running is finished, thermal_equil=true
    LOGICAL                              :: change_is_accepted, not_converged
    REAL(8)                              :: beta,  ss_eigen_values(2*this%fopt%n), test2, d_xi1, xi_0(this%fopt%n)
    ! trial variables####################
    REAL(8), DIMENSION(this%fopt%n) :: trial_xi, trial_chi, trial_zeta
    real(8)                         :: new_E, new_E_wf, new_E_sw
    INTEGER, DIMENSION(this%fopt%nh, 2) :: trial_hole_chi, trial_hole_xi
    INTEGER, DIMENSION(this%n_pole) :: trial_w_s
    COMPLEX(8), DIMENSION(2*this%fopt%n, 2*this%fopt%n)           :: trial_wf
    ! ###################################
    !INTEGER, DIMENSION(this%fopt%nc2)     :: new_w
    INTEGER, DIMENSION(this%fopt%nh)     :: new_w
    !TYPE(fchi_optimizer), allocatable     :: new_fopt
    TYPE(fchi_optimizer_Aem), allocatable     :: new_fopt
    ! TYPE(spinwave),ALLOCATABLE           :: nswave ! not needed
    INTEGER                              :: i, iter_max,&
         & new_swave_info(this%fopt%n), m, temp, itemp, j, k, a
    ! CHARACTER(1)                         :: ii
    CHARACTER(1), OPTIONAL                         :: rs
    !TYPE(angular_variable)     ::  an_chi, an_xi
    TYPE(angular_variable), allocatable   :: an_xi, an_chi
    logical                               :: optimized_fchi
    integer                               :: not_conv

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer,dimension(:),allocatable :: wlist
    real(8)   , allocatable          :: current(:,:)
    allocate(current(this%fopt%n, this%fopt%n))
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !this%present_chi = an_chi%get_and_rebuilt_chi(this%fopt%hole, this%present_chi)
    !this%present_chi = an_chi%rebuilt_angval_chi(this%present_chi)
    allocate(new_fopt)

    beta = 1.d0/(this%kb*T)
    IF (thermal_equil) THEN
      iter_max = this%Titer_e
    ELSE
      iter_max = this%Titer_ne
    END IF
    ! initalization
    !trial_hole_chi(:, :) = this%fopt%hole(:, :)
    !trial_hole_xi(:, :) = this%fopt%hole(:, :)
    allocate(an_xi)
    allocate(an_chi)
      !an_chi = new_angular_variable(this%fopt%pshape, this%present_hole_chi)
      an_chi = new_angular_variable_3D(this%fopt%pshape, this%present_hole_chi)
    ! debugging-----------
        !jtrial_chi = an_chi%get_angular_variable()
        !jtrial_chi = an_chi%rebuilt_angval_chi(trial_chi)
    !CALL new_fopt%fchi_optimizer_init(this%fopt%pshape, this%present_hole_chi, this%fopt%t, &
    !     & this%fopt%u, this%fopt%jd, this%fopt%lm, this%present_xi, trial_chi, this%present_wf, 't')
    !call new_fopt%fchi_optimizer_init(this%fopt%pshape, this%fopt%hole, this%fopt%ne, this%fopt%t, &
    !      & this%fopt%u, this%fopt%jd, this%fopt%lm, this%present_xi, trial_chi, this%present_wf,  &
    !      & this%present_eg, this%present_ne_u, this%present_ne_d)
    call new_fopt%fchi_optimizer_Aem_init(this%fopt%pshape, this%fopt%hole, this%fopt%ne, this%fopt%t, &
          & this%fopt%u, this%fopt%jd, this%fopt%lm, this%Bz, this%present_xi, trial_chi, this%present_wf,  &
          & this%present_eg, this%present_ne_u, this%present_ne_d)
    call new_fopt%set_param(pflag = .false.)
    !write(*, *) this%fopt%u, this%fopt%jd, this%fopt%lm
    !deallocate(this%fopt)
    !allocate(this%fopt)
    !CALL this%fopt%fchi_optimizer_init(, this%present_hole_chi, this%fopt%t, &
     !    & this%fopt%u, this%fopt%jd, this%fopt%lm, this%present_xi, trial_chi, this%present_wf, 't')
        !CALL this%fopt%modify_init(this%present_hole_chi, &
        !& this%present_xi, trial_chi, this%present_wf)
        !jCALL new_fopt%self_consistent(an_chi%d_angval_chi(trial_chi),optimized_fchi)
        !trial_chi = an_chi%get_rebuilt_phase(new_fopt%delta_chi)
        !jtrial_chi = an_chi%get_rebuilt_phase(new_fopt%delta_chi)
        !jprint *,  "first ene at metropolis", & 
          !j& energy_total_energy(this%fopt%pshape, this%present_hole_chi, &
          !j& this%fopt%t, this%fopt%u, this%fopt%jd, this%fopt%lm, &
          !j!j& this%present_xi, trial_chi, this%present_wf, 'T')
          !write(99) this%present_xi
          !write(98) trial_chi
          !stop 'metropilis init_xi'
    ! debugging-----
    an_xi = new_angular_variable_3D(this%fopt%pshape, this%fopt%hole)
    xi_0 = an_xi%rebuilt_angval_xi(this%present_xi)
    !write(99) xi_0
    !stop 'xi_0'
    trial_hole_chi(:, 1) = this%fopt%hole(:,1)
    trial_hole_xi(:, 1) = this%fopt%hole(:,1)
    !do i = 1, 10
    !  if(allocated(new_fopt)) deallocate(new_fopt)
    !    allocate(new_fopt)
    !  
    !  enddo

    mc_main_loop: DO i = 1, iter_max
    !mc_main_loop: DO i = 1, iter_max + this%not_converged_count
    
    !i = 0
    !not_conv = 0
    !mc_main_loop: DO while(i < iter_max + not_conv )
    !  i = i + 1
    
    !mc_main_loop: DO while(this%add_ensemble_count < iter_max )
      !      write(*,'(a6,i5,a2,i6,a2)')'-iter',i,' /',iter_max,' -'
      !      trial_xi  =this%opt%xi
      !      this%trial_zeta=this%present_zeta
      !write(*, *) 'energy_metropolis', &
      !    & energy_total_energy(this%fopt%pshape, this%fopt%hole, &
      !    & this%fopt%t, this%fopt%u, this%fopt%jd, this%fopt%lm, &
      !    & this%present_xi, this%present_chi, this%present_wf, 'T')
      CALL this%change_state()
      SELECT CASE (this%change_mode)
      CASE (1)
        new_w = this%get_flip_w_b()
        trial_hole_chi(:, 2) = new_w
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !      trial_hole_chi(:,:) = this%present_hole_chi
  !      trial_hole_chi(6,2) = 1
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        trial_w_s = this%w_s
        an_chi = new_angular_variable_3D(this%fopt%pshape, trial_hole_chi)
        trial_chi = an_chi%get_angular_variable() ! bulk
        trial_chi = trial_chi + an_chi%get_atan_phase_on_xy(this%present_pole_s, trial_w_s, 1) ! surface
        trial_chi = an_chi%rebuilt_angval_chi(trial_chi)
        CALL new_fopt%modify_init(trial_hole_chi, this%present_xi, trial_chi, this%present_wf)
        
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !      print *,trial_hole_chi(:,2)
  !      CALL draw_phase_3D('chi0.dat', 10, this%fopt%pshape, trial_chi, trial_hole_chi(:, 1))
  !      call write_format_of_chi_sb('chi0_mc.plt', 11, this%fopt%pshape, trial_hole_chi, &
  !                            &     this%present_pole_s, this%w_s)
  !      call get_winding_number(trial_chi, this%fopt, wlist)
  !      !call draw_circle("wl_eta.dat", 50, pshape, clist%cil, clist%cpl, wlist)
  !      call draw_winding_circle("wp_chi_mc.dat", "wm_chi_mc.dat", 50, 51, this%fopt%pshape, this%fopt%cil, this%fopt%cpl, wlist)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
        !CALL new_fopt%self_consistent(optimized_fchi)
        !if (.not.optimized_fchi) then
        !  this%not_converged_count = this%not_converged_count + 1
        !  not_conv = not_conv + 1
        !  cycle
        !  ! or fchi single shot
        !endif
        
        CALL new_fopt%oneshot()

        trial_chi = an_chi%get_rebuilt_phase(new_fopt%delta_chi)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !      CALL draw_phase_3D('chi_mc.dat', 10, this%fopt%pshape, trial_chi, trial_hole_chi(:, 1))
  !      call new_fopt%get_Jp_matrix(current)
  !      call draw_current("current_mc.dat", 10, this%fopt%pshape, this%fopt%cpl, current)
  !      stop
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !new_E_wf = & 
        !  & energy_total_energy(this%fopt%pshape, this%fopt%hole, &
        !  & this%fopt%t, this%fopt%u, this%fopt%jd, this%fopt%lm, &
        !  & this%present_xi, trial_chi, this%present_wf, 'T')
        !new_E = total_energy_sb_2lay_xichi_T0(this%fopt, this%present_wf, &
        !        &                             this%present_xi, trial_chi)
        new_E = total_energy_sb_2lay_xichi_T0_Aem(this%fopt, this%present_wf, &
                &                             this%present_xi, trial_chi, new_fopt%Aem_dphase)
        CALL this%Judgment(beta, new_E, this%old_E, change_is_accepted)
        IF (change_is_accepted) THEN
          this%old_E = new_E
          this%fopt%w2 = this%get_flip_w_b()
          this%change_count(1) = this%change_count(1) + 1
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !      CALL draw_phase_3D('chi_mc.dat', 10, this%fopt%pshape, trial_chi, trial_hole_chi(:, 1))
  !      call new_fopt%get_Jp_matrix(current)
  !      call draw_current("current_mc.dat", 10, this%fopt%pshape, this%fopt%cpl, current)
  !      stop
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          this%chi = trial_chi
          m = this%wma(-1)
          IF (this%cw(1) <= this%wma(-1)) THEN ! m -> a
            temp = this%wma(this%cw(1))
            this%wma(this%cw(1)) = this%wma(m)
            this%wma(m) = temp
            this%wma(-1) = m - 1
            this%wn(temp) = (this%win(temp) - (-1))**2
          ELSE ! a -> m
            temp = this%wma(this%cw(1))
            this%wma(this%cw(1)) = this%wma(m + 1)
            this%wma(m + 1) = temp
            this%wma(-1) = m + 1
            this%wn(temp) = (this%win(temp) - (1))**2
          END IF
          call this%check_min_E(new_E, trial_chi, this%fopt%w2, this%w_s)
        ENDIF

      CASE (2)
        new_w = this%get_swap_w_b()
        trial_hole_chi(:, 2) = new_w
        trial_w_s = this%w_s
        an_chi = new_angular_variable_3D(this%fopt%pshape, trial_hole_chi)
        trial_chi = an_chi%get_angular_variable() ! bulk
        trial_chi = trial_chi + an_chi%get_atan_phase_on_xy(this%present_pole_s, trial_w_s, 1) ! surface
        trial_chi = an_chi%rebuilt_angval_chi(trial_chi)
        CALL new_fopt%modify_init(trial_hole_chi, this%present_xi, trial_chi, this%present_wf)
        
        !CALL new_fopt%self_consistent(optimized_fchi)
        !if (.not.optimized_fchi) then
        !  this%not_converged_count = this%not_converged_count + 1
        !  not_conv = not_conv + 1
        !  cycle
        !endif
        
        CALL new_fopt%oneshot()
        
        trial_chi = an_chi%get_rebuilt_phase(new_fopt%delta_chi)
        !new_E = total_energy_sb_2lay_xichi_T0(this%fopt, this%present_wf, &
        !        &                             this%present_xi, trial_chi)
        new_E = total_energy_sb_2lay_xichi_T0_Aem(this%fopt, this%present_wf, &
                &                             this%present_xi, trial_chi, new_fopt%Aem_dphase)
        CALL this%Judgment(beta, new_E, this%old_E, change_is_accepted)
        IF (change_is_accepted) THEN
          this%old_E = new_E
          this%fopt%w2 = this%get_swap_w_b()
          ! this%judg = this%judg + 1
          this%change_count(2) = this%change_count(2) + 1
          this%chi = trial_chi
          temp = this%wma(this%cw(2))
          itemp = this%wma(this%cw(3))
          this%wma(this%cw(2)) = itemp
          this%wma(this%cw(3)) = temp
          this%wn(temp) = (this%win(temp) - (this%fopt%w2(temp)))**2
          this%wn(itemp) = (this%win(itemp) - (this%fopt%w2(itemp)))**2
          call this%check_min_E(new_E, trial_chi, this%fopt%w2, this%w_s)
        ENDIF
      
      CASE (3)
        trial_w_s = this%get_flip_w_s()
        trial_hole_chi(:, 2) = this%fopt%w2
        an_chi = new_angular_variable_3D(this%fopt%pshape, trial_hole_chi)
        trial_chi = an_chi%get_angular_variable() ! bulk
        trial_chi = trial_chi + an_chi%get_atan_phase_on_xy(this%present_pole_s, trial_w_s, 1) ! surface
        trial_chi = an_chi%rebuilt_angval_chi(trial_chi)
        CALL new_fopt%modify_init(trial_hole_chi, this%present_xi, trial_chi, this%present_wf)
        
        !CALL new_fopt%self_consistent(optimized_fchi)
        !if (.not.optimized_fchi) then
        !  this%not_converged_count = this%not_converged_count + 1
        !  not_conv = not_conv + 1
        !  cycle
        !endif
        
        CALL new_fopt%oneshot()
        
        trial_chi = an_chi%get_rebuilt_phase(new_fopt%delta_chi)
        !new_E = total_energy_sb_2lay_xichi_T0(this%fopt, this%present_wf, &
        !        &                             this%present_xi, trial_chi)
        new_E = total_energy_sb_2lay_xichi_T0_Aem(this%fopt, this%present_wf, &
                &                             this%present_xi, trial_chi, new_fopt%Aem_dphase)
        CALL this%Judgment(beta, new_E, this%old_E, change_is_accepted)
        IF (change_is_accepted) THEN
          this%old_E = new_E
          this%w_s = this%get_flip_w_s()
          this%change_count(3) = this%change_count(3) + 1
          this%chi = trial_chi
          !m = this%wma_s(-1)
          !IF (this%cw(1) <= this%wma_s(-1)) THEN ! m -> a
          !  temp = this%wma_s(this%cw(1))
          !  this%wma_s(this%cw(1)) = this%wma_s(m)
          !  this%wma_s(m) = temp
          !  this%wma_s(-1) = m - 1
          !  this%wn_s(temp) = (this%present_w_s(temp) - (-1))**2
          !ELSE ! a -> m
          !  temp = this%wma_s(this%cw(1))
          !  this%wma_s(this%cw(1)) = this%wma_s(m + 1)
          !  this%wma_s(m + 1) = temp
          !  this%wma_s(-1) = m + 1
          !  this%wn_s(temp) = (this%present_w_s(temp) - (1))**2
          !END IF
          m = this%wma0_s(-1)
          IF (this%cw(1) <= this%wma0_s(-1)) THEN ! m -> a
            temp = this%wma0_s(this%cw(1))
            this%wma0_s(this%cw(1)) = this%wma0_s(m)
            this%wma0_s(m) = temp
            this%wma0_s(-1) = this%wma0_s(-1) - 1
            this%wma0_s(-2) = this%wma0_s(-2) + 1
            this%wn_s(temp) = (this%present_w_s(temp) - (-1))**2
          ELSE ! a -> m
            temp = this%wma0_s(this%cw(1))
            this%wma0_s(this%cw(1)) = this%wma0_s(m + 1)
            this%wma0_s(m + 1) = temp
            this%wma0_s(-1) = this%wma0_s(-1) + 1
            this%wma0_s(-2) = this%wma0_s(-2) - 1
            this%wn_s(temp) = (this%present_w_s(temp) - (1))**2
          END IF
          call this%check_min_E(new_E, trial_chi, this%fopt%w2, this%w_s)
        ENDIF
      
      CASE (4)
        trial_w_s = this%get_swap_w_s()
        trial_hole_chi(:, 2) = this%fopt%w2
        an_chi = new_angular_variable_3D(this%fopt%pshape, trial_hole_chi)
        trial_chi = an_chi%get_angular_variable() ! bulk
        trial_chi = trial_chi + an_chi%get_atan_phase_on_xy(this%present_pole_s, trial_w_s, 1) ! surface
        trial_chi = an_chi%rebuilt_angval_chi(trial_chi)
        CALL new_fopt%modify_init(trial_hole_chi, this%present_xi, trial_chi, this%present_wf)
        
        !CALL new_fopt%self_consistent(optimized_fchi)
        !if (.not.optimized_fchi) then
        !  this%not_converged_count = this%not_converged_count + 1
        !  not_conv = not_conv + 1
        !  cycle
        !endif
        
        CALL new_fopt%oneshot()
        
        trial_chi = an_chi%get_rebuilt_phase(new_fopt%delta_chi)
        !new_E = total_energy_sb_2lay_xichi_T0(this%fopt, this%present_wf, &
        !        &                             this%present_xi, trial_chi)
        new_E = total_energy_sb_2lay_xichi_T0_Aem(this%fopt, this%present_wf, &
                &                             this%present_xi, trial_chi, new_fopt%Aem_dphase)
        CALL this%Judgment(beta, new_E, this%old_E, change_is_accepted)
        IF (change_is_accepted) THEN
          this%old_E = new_E
          this%w_s = this%get_swap_w_s()
          ! this%judg = this%judg + 1
          this%change_count(4) = this%change_count(4) + 1
          this%chi = trial_chi
          !temp = this%wma_s(this%cw(2))
          !itemp = this%wma_s(this%cw(3))
          !this%wma_s(this%cw(2)) = itemp
          !this%wma_s(this%cw(3)) = temp
          temp = this%wma0_s(this%cw(2))
          itemp = this%wma0_s(this%cw(3))
          this%wma0_s(this%cw(2)) = itemp
          this%wma0_s(this%cw(3)) = temp
          this%wn_s(temp) = (this%present_w_s(temp) - this%w_s(temp))**2
          this%wn_s(itemp) = (this%present_w_s(itemp) - this%w_s(itemp))**2
          call this%check_min_E(new_E, trial_chi, this%fopt%w2, this%w_s)
        ENDIF

      CASE (5)
        trial_w_s = this%get_annihilate_w_s()
        trial_hole_chi(:, 2) = this%fopt%w2
        an_chi = new_angular_variable_3D(this%fopt%pshape, trial_hole_chi)
        trial_chi = an_chi%get_angular_variable() ! bulk
        trial_chi = trial_chi + an_chi%get_atan_phase_on_xy(this%present_pole_s, trial_w_s, 1) ! surface
        trial_chi = an_chi%rebuilt_angval_chi(trial_chi)
        CALL new_fopt%modify_init(trial_hole_chi, this%present_xi, trial_chi, this%present_wf)

        CALL new_fopt%oneshot()
        
        trial_chi = an_chi%get_rebuilt_phase(new_fopt%delta_chi)
        !new_E = total_energy_sb_2lay_xichi_T0(this%fopt, this%present_wf, &
        !        &                             this%present_xi, trial_chi)
        new_E = total_energy_sb_2lay_xichi_T0_Aem(this%fopt, this%present_wf, &
                &                             this%present_xi, trial_chi, new_fopt%Aem_dphase)
        !print *, "5 : ", new_E, this%old_E
        CALL this%Judgment(beta, new_E, this%old_E, change_is_accepted)
        IF (change_is_accepted) THEN
          !print *, "    -> accepted"
          this%old_E = new_E
          this%w_s = trial_w_s
          this%change_count(5) = this%change_count(5) + 1
          this%chi = trial_chi
          
          m = this%wma0_s(-1)
          a = this%wma0_s(-2)

          itemp = this%wma0_s(this%cw(3))
          this%wma0_s(this%cw(3)) = this%wma0_s(m+a)
          this%wma0_s(m+a) = itemp

          temp = this%wma0_s(this%cw(2))
          this%wma0_s(this%cw(2)) = this%wma0_s(m)
          this%wma0_s(m) = temp

          temp = this%wma0_s(m)
          this%wma0_s(m) = this%wma0_s(m+a-1)
          this%wma0_s(m+a-1) = temp
          
          this%wma0_s(-1) = this%wma0_s(-1) - 1
          this%wma0_s(-2) = this%wma0_s(-2) - 1
          this%wma0_s(-3) = this%wma0_s(-3) + 2

          this%wn_s(temp) = (this%present_w_s(temp) - this%w_s(temp))**2
          this%wn_s(itemp) = (this%present_w_s(itemp) - this%w_s(itemp))**2
          
          call this%check_min_E(new_E, trial_chi, this%fopt%w2, this%w_s)
        ENDIF

      CASE (6)
        trial_w_s = this%get_create_w_s()
        trial_hole_chi(:, 2) = this%fopt%w2
        an_chi = new_angular_variable_3D(this%fopt%pshape, trial_hole_chi)
        trial_chi = an_chi%get_angular_variable() ! bulk
        trial_chi = trial_chi + an_chi%get_atan_phase_on_xy(this%present_pole_s, trial_w_s, 1) ! surface
        trial_chi = an_chi%rebuilt_angval_chi(trial_chi)
        CALL new_fopt%modify_init(trial_hole_chi, this%present_xi, trial_chi, this%present_wf)

        CALL new_fopt%oneshot()
        
        trial_chi = an_chi%get_rebuilt_phase(new_fopt%delta_chi)
        !new_E = total_energy_sb_2lay_xichi_T0(this%fopt, this%present_wf, &
        !        &                             this%present_xi, trial_chi)
        new_E = total_energy_sb_2lay_xichi_T0_Aem(this%fopt, this%present_wf, &
                &                             this%present_xi, trial_chi, new_fopt%Aem_dphase)
        !print *, "6 : ", new_E, this%old_E
        CALL this%Judgment(beta, new_E, this%old_E, change_is_accepted)
        IF (change_is_accepted) THEN
          !print *, "    -> accepted"
          this%old_E = new_E
          this%w_s = trial_w_s
          this%change_count(6) = this%change_count(6) + 1
          this%chi = trial_chi
          
          m = this%wma0_s(-1)
          a = this%wma0_s(-2)

          temp = this%wma0_s(this%cw(2))
          this%wma0_s(this%cw(2)) = this%wma0_s(m+a+1)
          this%wma0_s(m+a+1) = temp

          temp = this%wma0_s(m+a+1)
          this%wma0_s(m+a+1) = this%wma0_s(m+1)
          this%wma0_s(m+1) = temp
          
          itemp = this%wma0_s(this%cw(3))
          this%wma0_s(this%cw(3)) = this%wma0_s(m+a+2)
          this%wma0_s(m+a+2) = itemp

          this%wma0_s(-1) = this%wma0_s(-1) + 1
          this%wma0_s(-2) = this%wma0_s(-2) + 1
          this%wma0_s(-3) = this%wma0_s(-3) - 2

          this%wn_s(temp) = (this%present_w_s(temp) - this%w_s(temp))**2
          this%wn_s(itemp) = (this%present_w_s(itemp) - this%w_s(itemp))**2
          
          call this%check_min_E(new_E, trial_chi, this%fopt%w2, this%w_s)
        ENDIF

      END SELECT

      !   write(*,*)this%opt%hole(:,2)
      ! write(*,*)this%opt%chi
      !   write(*,*)new_w
      ! write(*,*)nopt%chi

      !CALL this%Judgment(beta, new_E, this%old_E, change_is_accepted)
      CALL this%ensemble_average(thermal_equil)
    END DO mc_main_loop
      deallocate(an_chi)
      deallocate(an_xi)
      deallocate(new_fopt)
        !write(99) this%present_xi
        !write(98) this%present_chi
        !print *, 'aaa'
  END SUBROUTINE monte_carlo_2lay_Metropolis

  SUBROUTINE monte_carlo_2lay_Judgment(this, beta, new_E, old_E, change_is_accepted)
    CLASS(monte_carlo_2lay), INTENT(inout)     :: this
    REAL(8), INTENT(in)                   :: beta, new_E
    !REAL(8), INTENT(inout)                :: old_E
    REAL(8), INTENT(in)                   :: old_E
    LOGICAL, INTENT(out)                  :: change_is_accepted
    REAL(8)                               :: x, dE
    !character(1)                     :: ii
    !write(*,*) 'new_E',new_E, old_E, beta
    dE = new_E - old_E
    IF (dE < 0d0) THEN
      change_is_accepted = .true.
      return
    else
      CALL RANDOM_NUMBER(x)
      IF (EXP(-beta*dE) > x) THEN
        !old_E = new_E
        !      write(*,*) ii,' ','new_E',new_E
        change_is_accepted = .TRUE.
      ELSE
        change_is_accepted = .FALSE.
      END IF
    end if
    !print *, EXP(-beta*dE), x, change_is_accepted
  END SUBROUTINE monte_carlo_2lay_Judgment

  ! SUBROUTINE monte_carlo_2lay_check_w(this, beta,new_E,old_E,cwflag)
  !   CLASS(monte_carlo_2lay), INTENT(inout)          :: this
  !   REAL(8), INTENT(in)                   :: beta, new_E
  !   LOGICAL, INTENT(out)                            :: cwflag
  !   REAL(8)                                         :: x,old_E
  !   IF (EXP(-beta*(new_E-old_E)) > x) THEN
  !      old_E = new_E
  !      cwflag = .TRUE.
  !   ELSE
  !      cwflag = .FALSE.
  !   END IF
  ! END SUBROUTINE monte_carlo_2lay_check_w

  !SUBROUTINE monte_carlo_2lay_set_near_circle(this)
  !  CLASS(monte_carlo_2lay), INTENT(inout)    :: this
  !  INTEGER                                          ::i, j
  !  IF (ALLOCATED(this%near_circle)) DEALLOCATE (this%near_circle)
  !  ALLOCATE (this%near_circle(this%fopt%nc2))
  !  DO i = 1, this%fopt%nc2
  !    DO j = 1, this%fopt%nc2
  !      IF (this%fopt%hole(i, 1) - 2*this%fopt%nx - 2 == this%fopt%hole(j, 1)) &
  !           & CALL this%near_circle(i)%append2(this%fopt%hole(j, 1))
  !      IF (this%fopt%hole(i, 1) - 2*this%fopt%nx == this%fopt%hole(j, 1)) &
  !           & CALL this%near_circle(i)%append2(this%fopt%hole(j, 1))
  !      IF (this%fopt%hole(i, 1) - 2*this%fopt%nx + 2 == this%fopt%hole(j, 1)) &
  !           & CALL this%near_circle(i)%append2(this%fopt%hole(j, 1))
  !      IF (this%fopt%hole(i, 1) - 2 == this%fopt%hole(j, 1)) &
  !           & CALL this%near_circle(i)%append2(this%fopt%hole(j, 1))
  !      CALL this%near_circle(i)%append2(this%fopt%hole(i, 1))
  !      IF (this%fopt%hole(i, 1) + 2 == this%fopt%hole(j, 1)) &
  !           & CALL this%near_circle(i)%append2(this%fopt%hole(j, 1))
  !      IF (this%fopt%hole(i, 1) + 2*this%fopt%nx - 2 == this%fopt%hole(j, 1)) &
  !           & CALL this%near_circle(i)%append2(this%fopt%hole(j, 1))
  !      IF (this%fopt%hole(i, 1) + 2*this%fopt%nx == this%fopt%hole(j, 1)) &
  !           & CALL this%near_circle(i)%append2(this%fopt%hole(j, 1))
  !      IF (this%fopt%hole(i, 1) + 2*this%fopt%nx + 2 == this%fopt%hole(j, 1)) &
  !           & CALL this%near_circle(i)%append2(this%fopt%hole(j, 1))
  !    ENDDO
  !  ENDDO
  !END SUBROUTINE monte_carlo_2lay_set_near_circle

  SUBROUTINE monte_carlo_2lay_set_ran_seed(this)
    CLASS(monte_carlo_2lay), INTENT(inout)    :: this
    INTEGER, DIMENSION(:), ALLOCATABLE               :: seed
    REAL(8), DIMENSION(:), ALLOCATABLE    :: ranx
    INTEGER, DIMENSION(:), ALLOCATABLE               :: tclock
    INTEGER                                          :: nrand, t1, i
    CHARACTER                                        :: cfile(20)
    CALL RANDOM_SEED(size=nrand)
    ALLOCATE (seed(nrand))
    ALLOCATE (ranx(nrand))
    CALL SYSTEM_CLOCK(count=t1)
    seed = t1
    CALL RANDOM_SEED(put=seed)
    DO i = 1, 4
      CALL RANDOM_NUMBER(ranx)
    END DO
    seed = ranx*(2**30)
    DO i = 1, nrand
      !    write(6,*)i, seed(i), ranx(i)
    END DO
    CALL RANDOM_SEED(put=seed)
    DEALLOCATE (seed)
    DEALLOCATE (ranx)

  END SUBROUTINE monte_carlo_2lay_set_ran_seed


  !SUBROUTINE monte_carlo_2lay_calc_winding_number(this, xi, chi, hole_xi, hole_chi)
  !  CLASS(monte_carlo_2lay), INTENT(in)     :: this
  !  REAL(8), INTENT(in)                     :: xi(this%fopt%n), chi(this%fopt%n)
  !  INTEGER, INTENT(inout)                  :: hole_xi(this%fopt%nh, this%fopt%nh), &
  !       &                                     hole_chi(this%fopt%nh, this%fopt%nh)
  !  INTEGER                                 :: i, j
  !  REAL(8)                                 :: tau(4), tau_tot, amari, eta(this%fopt%n)
  !  TYPE(path)                  :: p(4)
  !  TYPE(hop_iterator), ALLOCATABLE           :: iter
  !
  !  DO i = 1, this%fopt%n
  !    eta(i) = xi(i) - pi*(i - 1)
  !    amari = MODULO(eta(i), 2.d0*pi)
  !    IF (amari > pi) amari = amari - 2.d0*pi
  !    eta(i) = amari
  !  ENDDO
  !  ALLOCATE (iter, source=new_hop_iterator(this%fopt%pshape))
  !  ! Calc the winding numbers for xi
  !  DO i = 1, SIZE(hole_xi(:, 1))
  !    tau_tot = 0.d0
  !    p(1)%i = hole_xi(i, 1) - this%fopt%nx
  !    p(1)%f = hole_xi(i, 1) + 1
  !    p(2)%i = hole_xi(i, 1) - 1
  !    p(2)%f = hole_xi(i, 1) + this%fopt%nx
  !    p(3)%i = hole_xi(i, 1) - this%fopt%nx
  !    p(3)%f = hole_xi(i, 1) - 1
  !    p(4)%i = hole_xi(i, 1) + 1
  !    p(4)%f = hole_xi(i, 1) + this%fopt%nx
  !    tau(1) = eta(p(4)%f) - eta(p(4)%i)
  !    tau(2) = eta(p(2)%i) - eta(p(2)%f)
  !    tau(3) = eta(p(3)%i) - eta(p(3)%f)
  !    tau(4) = eta(p(1)%f) - eta(p(1)%i)
  !    do j = 1, 4
  !      IF (abs(tau(j)) > pi) then
  !        tau(j) = modulo(tau(j), 2.d0*pi)
  !        if (tau(j) > pi) tau(j) = tau(j) - 2.d0*pi
  !      endif
  !    enddo
  !    !write(*,*) i,hole_xi(i,1)
  !    !WRITE(*,*) 'tau1', tau(1)
  !    !WRITE(*,*) 'tau2',tau(2)
  !    !WRITE(*,*) 'tau3',tau(3)
  !    !WRITE(*,*) 'tau4',tau(4)
  !    tau_tot = sum(tau(1:4))
  !    ! WRITE(*,*) 'pos, hole_xi',hole_xi(i,1), NINT(tau_tot/(2.d0*pi))
  !    hole_xi(i, 2) = NINT(tau_tot/(2d0*pi))
  !  ENDDO
  !  ! Calc the winding numbers for chi
  !  DO i = 1, SIZE(hole_chi(:, 1))
  !    IF (hole_xi(i, 2) == 0) THEN
  !      hole_chi(i, 2) = 0
  !      CYCLE
  !    ENDIF
  !    tau_tot = 0.d0
  !    p(1)%i = hole_chi(i, 1) - this%fopt%nx
  !    p(1)%f = hole_chi(i, 1) + 1
  !    p(2)%i = hole_chi(i, 1) - 1
  !    p(2)%f = hole_chi(i, 1) + this%fopt%nx
  !    p(3)%i = hole_chi(i, 1) - this%fopt%nx
  !    p(3)%f = hole_chi(i, 1) - 1
  !    p(4)%i = hole_chi(i, 1) + 1
  !    p(4)%f = hole_chi(i, 1) + this%fopt%nx
  !    tau(1) = chi(p(4)%f) - chi(p(4)%i)
  !    tau(2) = chi(p(2)%i) - chi(p(2)%f)
  !    tau(3) = chi(p(3)%i) - chi(p(3)%f)
  !    tau(4) = chi(p(1)%f) - chi(p(1)%i)
  !    do j = 1, 4
  !      IF (abs(tau(j)) > pi) then
  !        tau(j) = modulo(tau(j), 2.d0*pi)
  !        if (tau(j) > pi) tau(j) = tau(j) - 2.d0*pi
  !      endif
  !    enddo
  !    !write(*,*) i,hole_chi(i,1)
  !    !WRITE(*,*) 'tau1', tau(1)
  !    !WRITE(*,*) 'tau2',tau(2)
  !    !WRITE(*,*) 'tau3',tau(3)
  !    !WRITE(*,*) 'tau4',tau(4)
  !    tau_tot = sum(tau(1:4))
  !    ! WRITE(*,*) 'pos, hole_chi',hole_chi(i,1), NINT(tau_tot/(2d0*pi))
  !    hole_chi(i, 2) = NINT(tau_tot/(2d0*pi))
  !  ENDDO
  !END SUBROUTINE monte_carlo_2lay_calc_winding_number

  !SUBROUTINE monte_carlo_2lay_calc_Ms(this) !!! -> change into calc by spinwave
  !  CLASS(monte_carlo_2lay), INTENT(inout)    :: this
  !  REAL(8), DIMENSION(this%fopt%n)          :: Sx, Sy, Sz
  !  REAL(8), DIMENSION(this%fopt%nh)         :: Shx, Shy, Shz
  !  REAL(8)                                :: S
  !  INTEGER                    :: i, j
  !  INTEGER                   :: r(8)
  !  ! S=0.5d0
  !  ! Sx=0.d0
  !  ! Sy=0.d0
  !  ! Sz=0.d0
  !  ! DO i=1,this%fopt%n
  !  !    Sx(i) = S*COS(this%present_xi(i))*SIN(this%present_zeta(i))
  !  !    Sy(i) = S*SIN(this%present_xi(i))*SIN(this%present_zeta(i))
  !  !    Sz(i) = S*COS(this%present_zeta(i))
  !  ! ENDDO
  !  ! DO i=1,this%fopt%nh
  !  !    Sx(this%fopt%hole(i,1))=0.d0
  !  !    Sy(this%fopt%hole(i,1))=0.d0
  !  !    Sz(this%fopt%hole(i,1))=0.d0
  !  ! ENDDO
  !  ! write(*,*)'Sx',Sx
  !  ! write(*,*)'Sy',Sy
  !  ! write(*,*)'Sz',Sz
  !  Sx(:) = this%Sx(:)
  !  Sy(:) = this%Sy(:)
  !  Sz(:) = this%Sz(:)
  !  this%Mab = 0.d0
  !  this%Mc = 0.d0
  !  DO i = 1, this%fopt%n
  !    this%Mab = this%Mab + Sx(i)**2 + Sy(i)**2
  !    this%Mc = this%Mc + Sz(i)**2
  !  ENDDO

  !  Shx = 0.d0
  !  Shy = 0.d0
  !  Shz = 0.d0
  !  this%Mhab = 0.d0
  !  this%Mhc = 0.d0
  !  DO i = 1, this%fopt%nh
  !    r = this%get_around_hole_pos(this%fopt%hole(i, 1))
  !    DO j = 1, 8
  !      Shx(i) = Shx(i) + S*COS(this%present_xi(r(j)))*SIN(this%present_zeta(r(j)))
  !      Shy(i) = Shy(i) + S*SIN(this%present_xi(r(j)))*SIN(this%present_zeta(r(j)))
  !      Shz(i) = Shz(i) + S*COS(this%present_zeta(r(j)))
  !    ENDDO
  !    this%Mhab = this%Mhab + Shx(i)**2 + Shy(i)**2
  !    this%Mhc = this%Mhc + Shz(i)**2
  !  ENDDO
  !END SUBROUTINE monte_carlo_2lay_calc_Ms

  !SUBROUTINE monte_carlo_2lay_calc_Ms2(this)
  !  CLASS(monte_carlo_2lay), INTENT(inout)    :: this
  !  REAL(8), DIMENSION(this%fopt%n)          :: Sx, Sy, Sz
  !  REAL(8), DIMENSION(this%fopt%nh)         :: Shx, Shy, Shz
  !  REAL(8)                                :: S
  !  INTEGER                    :: i, j
  !  INTEGER                   :: r(8)
  !  S = 0.5d0
  !  Sx = 0.d0
  !  Sy = 0.d0
  !  Sz = 0.d0
  !  DO i = 1, this%fopt%n
  !    Sx(i) = S*COS(this%present_xi(i))*SIN(this%present_zeta(i))
  !    Sy(i) = S*SIN(this%present_xi(i))*SIN(this%present_zeta(i))
  !    Sz(i) = S*COS(this%present_zeta(i))
  !  ENDDO
  !  DO i = 1, this%fopt%nh
  !    Sx(this%fopt%hole(i, 1)) = 0.d0
  !    Sy(this%fopt%hole(i, 1)) = 0.d0
  !    Sz(this%fopt%hole(i, 1)) = 0.d0
  !  ENDDO
  !  ! write(*,*)'Sx',Sx
  !  ! write(*,*)'Sy',Sy
  !  ! write(*,*)'Sz',Sz
  !  this%Mab2 = 0.d0
  !  this%Mc2 = 0.d0
  !  DO i = 1, this%fopt%n
  !    this%Mab2 = this%Mab + Sx(i)**2 + Sy(i)**2
  !    this%Mc2 = this%Mc + Sz(i)**2
  !  ENDDO

  !  Shx = 0.d0
  !  Shy = 0.d0
  !  Shz = 0.d0
  !  this%Mhab2 = 0.d0
  !  this%Mhc2 = 0.d0
  !  DO i = 1, this%fopt%nh
  !    r = this%get_around_hole_pos(this%fopt%hole(i, 1))
  !    DO j = 1, 8
  !      Shx(i) = Shx(i) + S*COS(this%present_xi(r(j)))*SIN(this%present_zeta(r(j)))
  !      Shy(i) = Shy(i) + S*SIN(this%present_xi(r(j)))*SIN(this%present_zeta(r(j)))
  !      Shz(i) = Shz(i) + S*COS(this%present_zeta(r(j)))
  !    ENDDO
  !    this%Mhab2 = this%Mhab + Shx(i)**2 + Shy(i)**2
  !    this%Mhc2 = this%Mhc + Shz(i)**2
  !  ENDDO
  !END SUBROUTINE monte_carlo_2lay_calc_Ms2

  SUBROUTINE monte_carlo_2lay_calc_init(this, RorD, T, dT)
    CLASS(monte_carlo_2lay), INTENT(inout)    :: this
    CHARACTER, INTENT(in)                :: RorD ! R is rise, D is decent (T)
    REAL(8), INTENT(out)                 :: T, dT
    INTEGER                              :: i, m = 0, a = 0
    !INTEGER, DIMENSION(this%fopt%nc2)     :: wm, wa
    INTEGER, DIMENSION(this%fopt%nh)     :: wm, wa
    REAL(8)                              :: s, temp
    INTEGER                              :: m_s = 0, a_s = 0, x_s = 0
    INTEGER, DIMENSION(this%n_pole)      :: wm_s, wa_s, w0_s, near_l
    integer                              :: num_near, j
    real(8)                              :: dx, dy
    SELECT CASE (RorD)
    CASE ("r", "R")
      s = 1.d0
      T = this%T_min
    CASE ("d", "D")
      s = -1.d0
      T = this%T_max
    CASE default
      WRITE (6, *) "Please do not put into RorD other than R and D. "
      STOP
    END SELECT
    WRITE (6, *) "start T is", T
    dT = s*(this%T_max - this%T_min)/REAL(this%iter_max, KIND(0d0))
    T = T - dT
    
    ! set bulk small polaron
    !DO i = 1, this%fopt%nc2
    DO i = 1, this%fopt%nh
      SELECT CASE (this%fopt%w2(i))
      CASE (1)
        m = m + 1
        wm(m) = i
      CASE (-1)
        a = a + 1
        wa(a) = i
      CASE default
        !WRITE (6, *) "Please do not put into w other than 1 and -1."
        WRITE (6, *) "Please do not put into w other than 1 and -1 on bulk."
        WRITE (6, '(a2,i2,4a,1x,i3)') "w(", i, ") is", this%fopt%w2(i)
        STOP
      END SELECT
    END DO
    ALLOCATE (this%wma(-1:a + m), source=[m, 0, wm(1:m), wa(1:a)])
    ALLOCATE (this%win(this%fopt%nh))
    ALLOCATE (this%wn(this%fopt%nh))
    this%win = this%fopt%hole(:, 2)
    this%wn = this%fopt%hole(:, 2)
    this%wn = (this%win - this%wn)**2
    WRITE (6, *) "wn", SUM(ABS(this%wn))
    
    ! set surface pole
    DO i = 1, this%n_pole
      SELECT CASE (this%w_s(i))
      CASE (1)
        m_s = m_s + 1
        wm_s(m_s) = i
      CASE (-1)
        a_s = a_s + 1
        wa_s(a_s) = i
      CASE (0)
        x_s = x_s + 1
        w0_s(x_s) = i
      CASE default
        !WRITE (6, *) "Please do not put into w other than 1 and -1."
        WRITE (6, *) "Please do not put into w other than 1, -1 and 0 on surface."
      !  WRITE (6, '(a2,i2,4a,1x,i3)') "w(", i, ") is", this%fopt%w2(i)
        STOP
      END SELECT
    END DO
    !allocate(this%wma_s(-1:a_s + m_s), source=[m_s, 0, wm_s(1:m_s), wa_s(1:a_s)])
    allocate(this%wma0_s(-3:this%n_pole), source=[x_s, a_s, m_s, 0, wm_s(1:m_s), wa_s(1:a_s), w0_s(1:x_s)])
    allocate(this%wn_s(this%n_pole))
    this%wn_s = (this%present_w_s - this%w_s)**2

    ! set near pole on surface
    allocate(this%near_pole(this%n_pole))
    do i=1, this%n_pole
      num_near = 0
      do j=1, this%n_pole
        if ( i == j ) cycle
        dx = this%present_pole_s(1,i) - this%present_pole_s(1,j)
        dy = this%present_pole_s(2,i) - this%present_pole_s(2,j)
        if(dx**2+dy**2 < this%near_range**2) then
          num_near = num_near + 1
          near_l(num_near) = j
          !print *, i, j
        end if
      end do
      allocate(this%near_pole(i)%v(num_near), source = near_l(1:num_near))
    end do


   ! OPEN (3000+this%myid, file=this%filename_myid("energy", "dat"))
   ! OPEN (3100+this%myid, file=this%filename_myid("specific_heat", "dat"))
   ! OPEN (3200+this%myid, file=this%filename_myid("energy_sw", "dat"))
   ! OPEN (3300+this%myid, file=this%filename_myid("specific_heat_sw", "dat"))
   ! OPEN (3500+this%myid, file=this%filename_myid("energy_wf", "dat"))
   ! OPEN (3600+this%myid, file=this%filename_myid("specific_heat_wf", "dat"))
   ! OPEN (3700+this%myid, file=this%filename_myid("avg_E_xi1", "dat"))
   ! OPEN (3800+this%myid, file=this%filename_myid("dev_E_xi1", "dat"))
   ! OPEN (3400+this%myid, file=this%filename_myid("windingnumber", "dat"))
   ! OPEN (4100+this%myid, file=this%filename_myid("nw", "dat"))
   ! OPEN (4200+this%myid, file=this%filename_myid("wn", "dat"))
   ! OPEN (4300+this%myid, file=this%filename_myid("nwc", "dat"))
   ! OPEN (4400+this%myid, file=this%filename_myid("wnc", "dat"))
   ! OPEN (4500+this%myid, file=this%filename_myid("Mab", "dat"))
   ! OPEN (4600+this%myid, file=this%filename_myid("Mc", "dat"))
   ! OPEN (4700+this%myid, file=this%filename_myid("Mhab", "dat"))
   ! OPEN (4800+this%myid, file=this%filename_myid("Mhc", "dat"))
   ! OPEN (4900+this%myid, file=this%filename_myid("Mab2", "dat"))
   ! OPEN (5000+this%myid, file=this%filename_myid("Mc2", "dat"))
   ! OPEN (5100+this%myid, file=this%filename_myid("Mhab2", "dat"))
   ! OPEN (5200+this%myid, file=this%filename_myid("Mhc2", "dat"))
   ! open (5400+this%myid, file=this%filename_myid("avg_crank_level", "dat"))
   
    OPEN (3000+this%myid, file=this%filename_myid("energy", "dat"))
    OPEN (3100+this%myid, file=this%filename_myid("specific_heat", "dat"))
    OPEN (4000+this%myid, file=this%filename_myid("dwn", "dat"))
    OPEN (4100+this%myid, file=this%filename_myid("variance_dwn", "dat"))
    OPEN (4200+this%myid, file=this%filename_myid("dwn_b", "dat"))
    OPEN (4300+this%myid, file=this%filename_myid("variance_dwn_b", "dat"))
    OPEN (4400+this%myid, file=this%filename_myid("dwn_s", "dat"))
    OPEN (4500+this%myid, file=this%filename_myid("variance_dwn_s", "dat"))
    
    OPEN (5000+this%myid, file=this%filename_myid("awn", "dat"))
    OPEN (5100+this%myid, file=this%filename_myid("variance_awn", "dat"))
    OPEN (5200+this%myid, file=this%filename_myid("awn_b", "dat"))
    OPEN (5300+this%myid, file=this%filename_myid("variance_awn_b", "dat"))
    OPEN (5400+this%myid, file=this%filename_myid("awn_s", "dat"))
    OPEN (5500+this%myid, file=this%filename_myid("variance_awn_s", "dat"))
    !    open(30, file="energy.dat")
    !    open(31, file="specific_heat.dat")
    !    open(34, file="windignumber.dat")
    !    open(41, file="nw.dat")
    !    open(42, file="wn.dat")
    !    open(43, file="nwc.dat")
    !    open(44, file="wnc.dat")
    WRITE (*, *) this%fopt%w2
    !write(*,*)this%opt%chi
    !this%old_E = energy_total_energy(this%fopt%pshape, this%fopt%hole, & ! will be deleted
      !& this%fopt%t, this%fopt%u, this%fopt%jd, this%fopt%lm, &
      !& this%present_xi, this%present_chi, this%fopt%wf, 'T')
  !this%total_energy(this%fopt%wf, this%present_chi, this%fopt%xi, this%Sx, this%Sy, this%Sz)
    !WRITE (*, *) 'total ene init', this%old_E
    !this%old_E = energy_total_energy(this%fopt%pshape, this%fopt%hole, &
      !& this%fopt%t, this%fopt%u, this%fopt%jd, this%fopt%lm, &
      !& this%present_xi, this%present_chi, this%fopt%wf, 'T')
    !ALLOCATE (this%present_wf(2*this%fopt%n, 2*this%fopt%n))
 !++   ALLOCATE (this%present_wf(4*this%fopt%n_acc, 4*this%fopt%n_acc))
 !++   this%present_wf = this%fopt%wf
 !   ! xi1 setting ____/
 !   this%present_crank_level = 1
 !   do i = 1, size(this%crank_eg)
 !     if (this%crank_eg(i) > 0) then
 !       this%crank_valid_level = i
 !       this%crank_eg(i) = 0d0
 !       exit
 !     endif
 !   enddo
 !   write(*, *) this%crank_valid_level
 !   write(*, *) this%crank_eg(1:this%crank_valid_level)
 !   this%old_E_sw = 0d0
 !   this%old_E_wf = &
 !         & energy_total_energy2(this%fopt%pshape, this%fopt%hole, &
 !         & this%fopt%t, this%fopt%u, this%fopt%jd, this%fopt%lm, &
 !         & this%present_xi, this%present_chi, this%present_wf, this%crank_EV(1), 'T')
 !   this%old_E = this%old_E_wf + this%crank_eg(1)
 !   !this%old_E_xi1 = this%old_E
 !   WRITE (*, *) 'total ene init', this%old_E, this%old_E_sw, this%old_E_wf
    !WRITE (*, *) 'Sz^2', DOT_PRODUCT(this%Sz, this%Sz)
    !this%old_E = total_energy_sb_2lay_xichi_T0(this%fopt, this%present_wf, this%present_xi, this%present_chi)
    this%old_E = total_energy_sb_2lay_xichi_T0_Aem(this%fopt, this%present_wf, &
                        & this%present_xi, this%present_chi, this%fopt%Aem_dphase)
    print *, ""
    print *, 'total ene init', this%old_E
    !print "(i3,a,f)", this%myid, ': total ene init', this%old_E
    
    ! histgram
    !allocate(this%hist_xi1(this%hist_intervel_number), source=0)
    !allocate(this%hist_classval(this%hist_intervel_number))
    !temp = 2d0*pi/this%hist_intervel_number
    !this%hist_classval(1) = temp/2d0
    !do i = 2, this%hist_intervel_number
    !  this%hist_classval(i) = this%hist_classval(i-1) + temp
    !enddo

    !write(99) this%present_xi
    !stop 'calc_init'
  END SUBROUTINE monte_carlo_2lay_calc_init

  SUBROUTINE monte_carlo_2lay_calc_finalize(this)
    CLASS(monte_carlo_2lay), INTENT(inout)    :: this
   ! CLOSE (3000+this%myid)
   ! CLOSE (3100+this%myid)
   ! CLOSE (3200+this%myid)
   ! CLOSE (3300+this%myid)
   ! CLOSE (3500+this%myid)
   ! CLOSE (3600+this%myid)
   ! CLOSE (3700+this%myid)
   ! CLOSE (3800+this%myid)
   ! CLOSE (3400+this%myid)
   ! CLOSE (4100+this%myid)
   ! CLOSE (4200+this%myid)
   ! CLOSE (4300+this%myid)
   ! CLOSE (4400+this%myid)
   ! CLOSE (4500+this%myid)
   ! CLOSE (4600+this%myid)
   ! CLOSE (4700+this%myid)
   ! CLOSE (4800+this%myid)
   ! CLOSE (4900+this%myid)
   ! CLOSE (5000+this%myid)
   ! CLOSE (5100+this%myid)
   ! CLOSE (5200+this%myid)
   ! CLOSE (5400+this%myid)
    CLOSE (3000+this%myid)
    CLOSE (3100+this%myid)
    CLOSE (4000+this%myid)
    CLOSE (4100+this%myid)
    CLOSE (4200+this%myid)
    CLOSE (4300+this%myid)
    CLOSE (4400+this%myid)
    CLOSE (4500+this%myid)
  END SUBROUTINE monte_carlo_2lay_calc_finalize

  !subroutine monte_carlo_2lay_set_change_partition(this) CLASS(monte_carlo_2lay), INTENT(inout)    :: this
  !  REAL(8)                                   :: change_probability(size(this%change_rate,1)), sum_rate, temp
  !  INTEGER                                   :: i
  !  sum_rate = sum(this%change_rate)
  !  change_probability(:) = this%change_rate(:) / sum_rate
  !  temp = 0d0
  !  do i = 1, size(this%change_rate, 1)-1
  !    if (this%change_rate(i) < 10d-10) cycle
  !    this%change_rate_partition(i) = temp + change_probability(i)
  !    temp = temp + change_probability(i)
  !  enddo
  !  write(*, *) 'relative prob', this%change_rate(:)
  !  write(*, *) 'absolute prob', change_probability
  !  write(*, *) 'temp', temp
  !  write(*, *) 'parition', this%change_rate_partition
  !end subroutine monte_carlo_2lay_set_change_partition

  SUBROUTINE monte_carlo_2lay_change_state(this)
    !------------------------------------
    ! Last update: Mamabe
    ! change_mode
    ! 1: Change the spin states S + spin-wave -> S'
    ! 2: Invert the winding number
    ! 3: Exchange the winding numbers
    ! 4: Make the spin-wave excitations
    ! 5: Change the angular variable xi(1)
    ! 6: Spin-wave excitations and Exchanging the winding number
    !------------------------------------
    CLASS(monte_carlo_2lay), INTENT(inout)    :: this
    !REAL(8), DIMENSION(7)                      :: x
    REAL(8), DIMENSION(4)                      :: x
    INTEGER                                    :: nh, i, cand, p, q, cand_l(this%n_pole)
    nh = this%fopt%nh

    !REDO (this way of writing is not good... )
    1200 continue

    CALL RANDOM_NUMBER(x)
    ! write(*,*) 'x'
    !write(*,*) 'x6', x(6)
    this%change_mode = 0
  !  do i = 1, nh
  !    if(this%fopt%w2(i) /= this%hole_chi_GS(i, 2)) goto 1100
  !  enddo
  !  !  do i=1,4
  !  !    write(*,*)i,x(i)
  !  !  enddo

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !x(1) = 1d0
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (0d0 <= x(1) .and. x(1) < 0.10d0 ) then
      ! Invert the winding number
      this%change_mode = 1
      this%cw(1) = INT(nh*x(2)) + 1
    
    else if (0.10d0 <= x(1) .and. x(1) < 0.50d0 ) then
      ! Exchange the winding numbers
      this%change_mode = 2
      ! cw(2) is m, cw(3) is a
      this%cw(2) = INT(this%wma(-1)*x(2)) + 1
      this%cw(3) = INT((nh - this%wma(-1))*x(3)) + 1 + this%wma(-1)
      !this%cwf   = .false. !old description
    
    !else if (0.5d0 <= x(1) .and. x(1) < 0.75d0 ) then
    else if (0.5d0 <= x(1) .and. x(1) < 0.6d0 ) then
      ! Invert the winding number on surface
      this%change_mode = 3
      !this%cw(1) = INT(this%n_pole*x(2)) + 1
      this%cw(1) = INT((this%wma0_s(-1)+this%wma0_s(-2))*x(2)) + 1
      if( this%wma0_s(-1) + this%wma0_s(-2) == 0) then
        !print *, "WARNING : m+a=0 (=>REDO)"
        goto 1200 ! REDO
      end if
    !else !if (0.75d0 <= x(1) .and. x(1) <= 1d0 )
    !else if (0.6d0 <= x(1) .and. x(1) < 0.8d0 ) then
    else if (0.6d0 <= x(1) .and. x(1) <= 1d0 ) then
      ! Exchange the winding numbers on surface
      this%change_mode = 4
      !this%cw(2) = INT(this%wma_s(-1)*x(2)) + 1
      !this%cw(3) = INT((this%n_pole - this%wma_s(-1))*x(3)) + 1 + this%wma_s(-1)
      this%cw(2) = INT(this%wma0_s(-1)*x(2)) + 1
      this%cw(3) = INT(this%wma0_s(-2)*x(3)) + 1 + this%wma0_s(-1)
      if( this%wma0_s(-1)==0 .or. this%wma0_s(-2) == 0) then
        !print *, "WARNING : m=0 or a=0 (=>REDO)"
        goto 1200 ! REDO
      endif
    else !if (0.8d0 <= x(1) .and. x(1) <= 1d0 )
      ! Create or Annihilate pole pair
      this%cw(2) = INT(this%n_pole*x(2)) + 1
      cand = 0
      p = this%wma0_s(this%cw(2))
      if(this%cw(2) <= this%wma0_s(-1) + this%wma0_s(-2)) then ! cw(2) is m or a => annihilate
        this%change_mode = 5
        do i = 1, size(this%near_pole(p)%v)
          if(this%w_s(this%near_pole(p)%v(i)) == -this%w_s(p)) then
            cand = cand + 1
            cand_l(cand) = this%near_pole(p)%v(i)
          end if
        end do
        if(cand == 0) then
          !print *, "WARNING : can't choose a pole pair to annihilate (=>REDO)"
          !print *, "     : ", this%present_pole_s(:,p)
          !stop
          goto 1200 ! REDO
        end if
        q = cand_l( INT(cand*x(3)) + 1 )
        do i = 1, this%n_pole
          if(this%wma0_s(i) == q) then
            this%cw(3) = i
          end if
        end do
        if( this%cw(2) > this%wma0_s(-1) ) then ! take only cw(2) = m, cw(3) = a. if cw(2) = a, swap
          i = this%cw(2)
          this%cw(2) = this%cw(3)
          this%cw(3) = i
        end if

      else ! cw(2) is 0 => create (cw(2)->m, cw(3)->a)
        this%change_mode = 6
        do i = 1, size(this%near_pole(p)%v)
          if(this%w_s(this%near_pole(p)%v(i)) == 0) then
            cand = cand + 1
            cand_l(cand) = this%near_pole(p)%v(i)
          end if
        end do
        if(cand == 0) then
          !print *, "WARNING : can't choose a pole pair to create(=>REDO)"
          !print *, "     : ", this%present_pole_s(:,p)
          !stop
          goto 1200 ! REDO
        end if
        q = cand_l( INT(cand*x(3)) + 1 )
        do i = 1, this%n_pole
          if(this%wma0_s(i) == q) then
            this%cw(3) = i
          end if
        end do
        if(x(4) > 0.5d0) then
          i = this%cw(2)
          this%cw(2) = this%cw(3)
          this%cw(3) = i
        end if
      end if
    end if

    !if (x(7) < this%sw_and_wn_rate) then
    !  this%csw_num = INT(x(1)*this%swave%valid_level) + 1 ! the number of label which spinwave excited
    !  this%cw(2) = INT(this%wma(-1)*x(2)) + 1
    !  this%cw(3) = INT((nh - this%wma(-1))*x(3)) + 1 + this%wma(-1)
    !  this%change_mode = 6
    !  return
    !endif
    !IF (x(5) < this%cssrate) THEN ! change spin state
    !  this%change_mode = 1
    !  RETURN
    !ENDIF
    !IF (x(6) < this%xirate) then
    !  this%trial_crank_level = INT(x(2)*this%crank_valid_level) + 1
    !  this%change_mode = 5
    !  return
    !ENDIF
    !1100 continue
    !IF (x(4) < this%cwsrate) THEN ! Change the winding number
    !  ! Invert the winding number
    !  IF (x(1) < this%cwrate .OR. this%wma(-1) == 0 .OR. this%wma(-1) == nh) THEN
      !  this%cw(1) = INT(nh*x(2)) + 1
      !  !this%cwf   = .true. !old description
      !  this%change_mode = 2
      !  RETURN
    !  ELSE ! Exchange the winding numbers
      !  ! cw(2) is m, cw(3) is a
      !  this%cw(2) = INT(this%wma(-1)*x(2)) + 1
      !  this%cw(3) = INT((nh - this%wma(-1))*x(3)) + 1 + this%wma(-1)
      !  !this%cwf   = .false. !old description
      !  this%change_mode = 3
      !  RETURN
    !  END IF
    !  !this%csf = .false. ! old description
    !ELSE ! make spin-wave excitations
    !  this%csw_num = INT(x(2)*this%swave%valid_level) + 1 ! the number of label which spinwave excited
    !  !this%cw(4) = int(nh*x(2))+1 ! old description
    !  !this%csf = .true. ! old description
    !  this%change_mode = 4
    !  RETURN
    !ENDIF
  END SUBROUTINE monte_carlo_2lay_change_state

  SUBROUTINE monte_carlo_2lay_Metropolis_init(this,counter)
    CLASS(monte_carlo_2lay), INTENT(inout)   :: this
    integer, intent(in)               :: counter
    character(3)                      :: num
    !integer, intent(in)                         :: iter
    !call this%monte_carlo_Metropolis_init(T, iter)
    !  if (abs(this%old_E-E) > this%error) then
    !    write(6,*) "error is", abs(this%old_E-E)
    !    stop
    !  end if
    !    call this%set_ran_seed()
    this%change_count(:) = 0
    this%not_converged_count = 0
    this%add_ensemble_count = 0
    this%not_add_ensemble_count = 0
    !this%SP_sw = 0
    !this%SP_wn = 0
    !this%SP_xi1 = 0
    !! this%judg = 0
    !!this%sjudg = 0
    !! this%swavejudg = 0
    !! histgram
    !!write(num, '(i3.3)') counter
    !!open (5300+this%myid, file=this%filename_myid("hist_xi1_"//num//"_","dat"))
    !this%cw = 0
    !this%w = 0.d0
    !this%ws = 0.d0
    !this%aw = 0.d0
    !this%awc = 0.d0
    !this%sE = 0.d0
    !this%sE_sw = 0.d0
    !this%sE_wf = 0.d0
    !!this%sxi1 = 0.d0
    !this%sE_xi1 = 0.d0
    !this%sE2 = 0.d0
    !this%sE2_sw = 0.d0
    !this%sE2_wf = 0.d0
    !!this%sxi12 = 0.d0
    !this%sE2_xi1 = 0.d0
    !this%sMc = 0.d0
    !this%sMab = 0.d0
    !this%sMhc = 0.d0
    !this%sMhab = 0.d0
    !this%sMc2 = 0.d0
    !this%sMab2 = 0.d0
    !this%sMhc2 = 0.d0
    !this%sMhab2 = 0.d0
    !this%scrank_level = 0d0
    !! histgram
    !!this%hist_xi1 = 0
    !!  call this%set_zeta()
    !!  this%old_E = this%total_energy(this%opt%wf,this%opt%chi,this%opt%xi)
    !!  write(*,*)'total ene init',this%old_E
    !!  this%present_wf=this%opt%wf
    !!  this%opt%w2(:)=this%opt%hole(:,2)
    this%sum_E = 0.d0
    this%sum_E2 = 0.d0
    this%sum_diff_w = 0.d0
    this%sum_diff_w2 = 0.d0
    this%sum_ave_w = 0.d0
    this%sum_ave_w2 = 0.d0
    this%sum_diff_w_b = 0.d0
    this%sum_diff_w_b2 = 0.d0
    this%sum_ave_w_b = 0.d0
    this%sum_ave_w_b2 = 0.d0
    this%sum_diff_w_s = 0.d0
    this%sum_diff_w_s2 = 0.d0
    this%sum_ave_w_s = 0.d0
    this%sum_ave_w_s2 = 0.d0
  END SUBROUTINE monte_carlo_2lay_Metropolis_init

  SUBROUTINE monte_carlo_2lay_Metropolis_finalize(this, T)
    CLASS(monte_carlo_2lay), INTENT(inout)    :: this
    REAL(8), INTENT(in)                  :: T
    INTEGER                                :: i
    REAL(8)                              :: avgn, avgn2, aE, aE2, temp, avgM, avgM2, &
      &                                   aE_sw, aE2_sw, aE_wf, aE2_wf, &
      &                                   aE_xi1, aE2_xi1, & !aE: the avarage of Energy
      &                                   axi1, axi12, &
      &                                   acrank_level
    real(8)                               :: N

    real(8)   , allocatable               :: eg(:), current(:,:)
    complex(8), allocatable               :: wf(:, :), wf1(:, :), wf2(:, :)
    real(8), dimension(:,:), allocatable      :: pole_s
    integer, dimension(:), allocatable        :: w_s
    integer                               :: j, hole_chi(this%fopt%nh,2)

    !write(*, *) 'stastical parameter: total, wn, sw, xi1', this%Titer_e, this%SP_wn, this%SP_sw, this%SP_xi1
    
    !N = REAL(this%Titer_e, 8)
    !N = 1d0/REAL(this%Titer_e - this%not_converged_count, 8)
    N = REAL(this%add_ensemble_count, 8)
    
    !aE = this%sE/REAL(this%Titer_e, 8)
    !aE2 = this%sE2/REAL(this%Titer_e, 8)
    aE = this%sum_E/N
    aE2 = this%sum_E2/N
    WRITE (3000+this%myid, '(2(e25.15,1x))') T, aE
    temp = (aE2 - (aE**2))/(this%kb*(T**2))
    !temp = aE2/(this%kb*T**2) - aE**2/(this%kb*T**2)  ! avoid digit loss
    WRITE (3100+this%myid, '(2(e25.15,1x))') T, temp
   
    avgn = this%sum_diff_w/N
    avgn2 = this%sum_diff_w2/N
    WRITE (4000+this%myid, *) T, avgn, avgn2
    WRITE (4100+this%myid, *) T, (avgn2 - (avgn**2))
    avgn = this%sum_ave_w/N
    avgn2 = this%sum_ave_w2/N
    WRITE (5000+this%myid, *) T, avgn, avgn2
    WRITE (5100+this%myid, *) T, (avgn2 - (avgn**2))
    
    avgn = this%sum_diff_w_b/N
    avgn2 = this%sum_diff_w_b2/N
    WRITE (4200+this%myid, *) T, avgn, avgn2
    WRITE (4300+this%myid, *) T, (avgn2 - (avgn**2))
    avgn = this%sum_ave_w_b/N
    avgn2 = this%sum_ave_w_b2/N
    WRITE (5200+this%myid, *) T, avgn, avgn2
    WRITE (5300+this%myid, *) T, (avgn2 - (avgn**2))
    
    avgn = this%sum_diff_w_s/N
    avgn2 = this%sum_diff_w_s2/N
    WRITE (4400+this%myid, *) T, avgn, avgn2
    WRITE (4500+this%myid, *) T, (avgn2 - (avgn**2))
    avgn = this%sum_ave_w_s/N
    avgn2 = this%sum_ave_w_s2/N
    WRITE (5400+this%myid, *) T, avgn, avgn2
    WRITE (5500+this%myid, *) T, (avgn2 - (avgn**2))
   ! if (this%SP_sw /= 0) then
   !   aE_sw = this%sE_sw/REAL(this%SP_sw, 8)
   !   aE2_sw = this%sE2_sw/REAL(this%SP_sw, 8)
   ! else
   !   aE_sw = 0d0
   !   aE2_sw = 0d0
   ! endif
   ! WRITE (3200+this%myid, '(2(e25.15,1x))') T, aE_sw
   ! temp = (aE2_sw - (aE_sw**2))/(this%kb*(T**2))
   ! WRITE (3300+this%myid, '(2(e25.15,1x))') T, temp
   ! if (this%SP_wn /=0 ) then
   !   aE_wf = this%sE_wf/REAL(this%SP_wn, 8)
   !   aE2_wf = this%sE2_wf/REAL(this%SP_wn, 8)
   ! else
   !   aE_wf = 0
   !   aE2_wf = 0
   ! endif
   ! WRITE (3500+this%myid, '(2(e25.15,1x))') T, aE_wf
   ! temp = (aE2_wf - (aE_wf**2))/(this%kb*(T**2))
   ! WRITE (3600+this%myid, '(2(e25.15,1x))') T, temp
   ! !if (this%SP_xi1 /=0 ) then
   ! !  axi1 = this%sxi1/REAL(this%SP_xi1, 8)
   ! !  axi12 = this%sxi12/REAL(this%SP_xi1, 8)
   ! !else
   ! !  axi1 = 0.0d0
   ! !  axi12 = 0.d0
   ! !endif
   ! !WRITE (3700+this%myid, '(2(e25.15,1x))') T, axi1
   ! !temp = sqrt(axi12 - axi1**2)
   ! !WRITE (3800+this%myid, '(2(e25.15,1x))') T, temp

   ! if (this%SP_xi1 /=0 ) then ! old
   !   aE_xi1 = this%sE_xi1/REAL(this%SP_xi1, 8)
   !   aE2_xi1 = this%sE2_xi1/REAL(this%SP_xi1, 8)
   ! else
   !   aE_xi1 = 0
   !   aE2_xi1 = 0
   ! endif
   ! WRITE (3700+this%myid, '(2(e25.15,1x))') T, aE_xi1
   ! temp = (aE2_xi1 - (aE_xi1**2))/(this%kb*(T**2))
   ! WRITE (3800+this%myid, '(2(e25.15,1x))') T, temp

   ! !   write(6,*) "ene", aE**2, this%sE2/real(this%Titer_e,8)
   ! !   write(6,'(3(a3,1x,f20.10, 2x))')"  T", T, "<E>", aE, "C", temp
   ! !   write(6,*)"judg ", this%judg
   ! !   write(6,*)"judg / Titer1", this%judg/real(this%Titer_e,8)
   ! avgn = this%w/REAL(this%Titer_e, KIND(0d0))
   ! avgn2 = this%ws/REAL(this%Titer_e, KIND(0d0))
   ! WRITE (4100+this%myid, *) T, avgn, avgn2
   ! WRITE (4300+this%myid, *) T, (avgn2 - (avgn**2))
   ! avgn = this%aw/REAL(this%Titer_e, KIND(0d0))
   ! avgn2 = this%awc/REAL(this%Titer_e, KIND(0d0))
   ! WRITE (4200+this%myid, *) T, avgn, avgn2
   ! WRITE (4400+this%myid, *) T, (avgn2 - (avgn**2))
   ! WRITE (3400+this%myid, '(2(a3,1x,f20.10, 2x))') "  T", T, "<E>", aE
   ! !DO i = 1, this%fopt%nc2
   ! DO i = 1, this%fopt%nh
   !   WRITE (3400+this%myid, fmt='(1x,i2)', advance="no") this%fopt%w2(i)
   ! END DO
   ! WRITE (3400+this%myid, *)
   ! avgM = this%sMab/REAL(this%Titer_e, 8)
   ! WRITE (4500+this%myid, *) T, avgM
   ! avgM = this%sMc/REAL(this%Titer_e, 8)
   ! WRITE (4600+this%myid, *) T, avgM
   ! avgM = this%sMhab/REAL(this%Titer_e, 8)
   ! WRITE (4700+this%myid, *) T, avgM
   ! avgM = this%sMhc/REAL(this%Titer_e, 8)
   ! WRITE (4800+this%myid, *) T, avgM
   ! avgM2 = this%sMab2/REAL(this%Titer_e, 8)
   ! WRITE (4900+this%myid, *) T, avgM2
   ! avgM2 = this%sMc2/REAL(this%Titer_e, 8)
   ! WRITE (5000+this%myid, *) T, avgM2
   ! avgM2 = this%sMhab2/REAL(this%Titer_e, 8)
   ! WRITE (5100+this%myid, *) T, avgM2
   ! avgM2 = this%sMhc2/REAL(this%Titer_e, 8)
   ! WRITE (5200+this%myid, *) T, avgM2
   ! acrank_level = this%scrank_level/dble(this%SP_xi1)
   ! write(5400+this%myid, *) T, acrank_level
   ! ! histgram
   ! !do i = 1, this%hist_intervel_number
   ! !  write (5300+this%myid, *) this%hist_classval(i), this%hist_xi1(i)
   ! !enddo
   ! close (5300+this%myid)
   ! !write(*, *) T, acrank_level
    !WRITE (*, '(a,7(a2, i0, 3x))') 'num of times', '1:', this%change_count(1), '2:', this%change_count(2), &
    !     & '3:', this%change_count(3), '4:', this%change_count(4), &
    !     & '5:', this%change_count(5), '6:', this%change_count(6), 'f:', this%not_converged_count

    if(this%myid == 0) then
      allocate(eg(4*this%fopt%n_acc)) 
      allocate(wf(4*this%fopt%n_acc, 4*this%fopt%n_acc)) 
      allocate(wf1(4*this%fopt%n_acc, 4*this%fopt%n_acc)) 
      allocate(wf2(4*this%fopt%n_acc, 4*this%fopt%n_acc)) 
      !call calc_eigen_energy_sb_2lay_wf_T0(this%fopt, eg, this%present_wf, this%present_xi, this%chi, this%present_mu)
      call calc_eigen_energy_sb_2lay_wf_T0_Aem(this%fopt, eg, this%present_wf, this%present_xi, this%chi, &
                      & this%present_mu, this%fopt%Aem_dphase)
      wf = this%present_wf
      call sort_wf(eg, wf)
      do i = 1, this%fopt%n_acc
        j = this%fopt%ets(i)
        wf1(4*i-3, :) = exp(-0.5d0*ui*this%chi(j))*exp(-0.5d0*ui*this%present_xi(j))*wf(4*i-3, :)
        wf1(4*i-2, :) = exp(-0.5d0*ui*this%chi(j))*exp( 0.5d0*ui*this%present_xi(j))*wf(4*i-2, :)
        wf1(4*i-1, :) = exp( 0.5d0*ui*this%chi(j))*exp( 0.5d0*ui*this%present_xi(j))*wf(4*i-1, :)
        wf1(4*i  , :) = exp( 0.5d0*ui*this%chi(j))*exp(-0.5d0*ui*this%present_xi(j))*wf(4*i  , :)
      
        wf2(4*i-3, :) = exp(-0.5d0*ui*this%present_xi(j))*wf(4*i-3, :)
        wf2(4*i-2, :) = exp( 0.5d0*ui*this%present_xi(j))*wf(4*i-2, :)
        wf2(4*i-1, :) = exp( 0.5d0*ui*this%present_xi(j))*wf(4*i-1, :)
        wf2(4*i  , :) = exp(-0.5d0*ui*this%present_xi(j))*wf(4*i  , :)
      end do
      

      open(1000 + this%myid, file = "eg_mc_" // string(this%iter) // ".bin", form="unformatted")
      open(1100 + this%myid, file = "wf1_mc_" // string(this%iter) // ".bin", form="unformatted")
      open(1200 + this%myid, file = "wf2_mc_" // string(this%iter) // ".bin", form="unformatted")
      
      write(1000 + this%myid) eg
      write(1100 + this%myid) wf1
      write(1200 + this%myid) wf2

      close(1000 + this%myid)
      close(1100 + this%myid)
      close(1200 + this%myid)
      
      open(1300 + this%myid,file = "eg_mc_" // string(this%iter) // ".dat")
      do i = 1, size(eg)
          write(1300 + this%myid, *) i, eg(i)
      end do
      close(1300 + this%myid)

      !CALL draw_phase_3D('fchi.dat', unit_num, pshape, chi, hole_chi(:, 1))
      
      allocate(current(this%fopt%n, this%fopt%n)) 
      this%fopt%delta_chi = this%fopt%d_angval(this%chi)
      this%fopt%u_Aeff = this%fopt%set_vector_potential(this%fopt%delta_chi)
      call this%fopt%get_Jp_matrix(current)
      call draw_current("current_mc_" // string(this%iter) // ".dat", 1000+this%myid, this%fopt%pshape, this%fopt%cpl, current)

      !call get_winding_point_surface(chi, clist, pole_s_conv, w_s_conv)
      hole_chi = this%present_hole_chi
      hole_chi(:,2) = this%fopt%w2
      call write_format_of_chi_sb('chi0_mc_' // string(this%iter) // '.plt', 1100+this%myid, this%fopt%pshape, hole_chi, &
                              &     this%present_pole_s, this%w_s)
      call get_chi_winding_number_with_eta_winding_point_surface(this%eta, this%chi, this%fopt, pole_s, w_s)
      call write_format_of_chi_sb('chi_mc_' // string(this%iter) // '.plt', 1200+this%myid, this%fopt%pshape, hole_chi, pole_s, w_s)
    end if

    if(this%iter == this%iter_max) then
      call this%out_min_E_state()
    end if
    
    ! output
    ! 1~4 : number of times the mode was accepted
    ! f : number of times not converged (when use oneshot, f = 0)
    ! t : number of times calculation to thermal equilibrium
    ! e : number of times calculation for ensemble
    
    WRITE(*, *) ""
    !WRITE (*, '(a,7(a2, i0, 3x))') 'num of times ', '1:', this%change_count(1), '2:', this%change_count(2), &
    !     & '3:', this%change_count(3), '4:', this%change_count(4), 'f:', this%not_converged_count, &
    !     & 't:', this%not_add_ensemble_count, 'e:', this%add_ensemble_count
    WRITE (*, '(a,9(a2, i0, 3x))') 'num of times ', '1:', this%change_count(1), '2:', this%change_count(2), &
         & '3:', this%change_count(3), '4:', this%change_count(4), &
         & '5:', this%change_count(5), '6:', this%change_count(6), 'f:', this%not_converged_count, &
         & 't:', this%not_add_ensemble_count, 'e:', this%add_ensemble_count
    WRITE (*, *) 'total ene fin', this%old_E
    
  END SUBROUTINE monte_carlo_2lay_Metropolis_finalize

  !subroutine monte_carlo_2lay_xi1_hist_counter(this, xi1)
  !  class(monte_carlo_2lay), intent(inout)        :: this
  !  real(8), intent(in)                           :: xi1
  !  integer                                       :: hist_class
  !  hist_class = int(xi1 / (2d0*pi/this%hist_intervel_number))
  !  this%hist_xi1(hist_class + 1) = this%hist_xi1(hist_class + 1) + 1
  !end subroutine monte_carlo_2lay_xi1_hist_counter

  SUBROUTINE monte_carlo_2lay_ensemble_average(this, thermal_equil)
    CLASS(monte_carlo_2lay), INTENT(inout)    :: this
    LOGICAL, INTENT(in)                        :: thermal_equil
    !    real(8)                                    :: temp
    INTEGER                                    :: itemp, m, temp
    real(8)     :: sdwn_b, sdwn_s, swn_b, swn_s

    IF (thermal_equil) THEN
      !this%sE = this%sE + this%old_E
      !this%sE2 = this%sE2 + this%old_E**2
      !this%sMc = this%sMc + this%Mc
      !this%sMab = this%sMab + this%Mab
      !this%sMhc = this%sMhc + this%Mhc
      !this%sMhab = this%sMhab + this%Mhab
      !this%sMc2 = this%sMc2 + this%Mc2
      !this%sMab2 = this%sMab2 + this%Mab2
      !this%sMhc2 = this%sMhc2 + this%Mhc2
      !this%sMhab2 = this%sMhab2 + this%Mhab2
      
      this%add_ensemble_count = this%add_ensemble_count + 1
      
      ! summation of Enegy
      this%sum_E = this%sum_E + this%old_E
      !this%sum_E2 = this%sum_E2 + this%old_E**2
      this%sum_E2 = this%sum_E2 + this%old_E*this%old_E
      
      ! summation of winding number on surface and bulk layer
      ! diff_w : diffirence from ground state
      ! ave_w : winding number in current state
      sdwn_b = sum(this%wn)
      sdwn_s = sum(this%wn_s)
      swn_b = sum(this%fopt%w2)
      swn_s = sum(this%w_s)
      this%sum_diff_w = this%sum_diff_w + sdwn_b + sdwn_s
      this%sum_diff_w2 = this%sum_diff_w2 + (sdwn_b + sdwn_s)*(sdwn_b + sdwn_s)
      this%sum_ave_w = this%sum_ave_w + swn_b + swn_s
      this%sum_ave_w2 = this%sum_ave_w2 + (swn_b + swn_s)*(swn_b + swn_s)
      
      ! summation of winding number on bulk layer
      this%sum_diff_w_b = this%sum_diff_w_b + sdwn_b
      this%sum_diff_w_b2 = this%sum_diff_w_b2 + sdwn_b*sdwn_b
      this%sum_ave_w_b = this%sum_ave_w_b + swn_b
      this%sum_ave_w_b2 = this%sum_ave_w_b2 + swn_b*swn_b
      
      ! summation of winding number on surface layer
      this%sum_diff_w_s = this%sum_diff_w_s + sdwn_s
      this%sum_diff_w_s2 = this%sum_diff_w_s2 + sdwn_s*sdwn_s
      this%sum_ave_w_s = this%sum_ave_w_s + swn_s
      this%sum_ave_w_s2 = this%sum_ave_w_s2 + swn_s*swn_s
      
 !     select case (this%change_mode)
 !       case(2,3)
 !         this%SP_wn = this%SP_wn + 1
 !         this%sE_wf = this%sE_wf + this%old_E_wf
 !         this%sE2_wf = this%sE2_wf + this%old_E_wf**2
 !         this%w = this%w + SUM(this%wn) ! diff of the winding numbers of chi
 !         this%ws = this%ws + (SUM(this%wn))**2
 !         this%aw = this%aw + SUM(this%fopt%w2)
 !         this%awc = this%awc + (SUM(this%fopt%w2))**2
 !       case(4)
 !         this%SP_sw = this%SP_sw + 1
 !         this%sE_sw = this%sE_sw + this%old_E_sw
 !         this%sE2_sw = this%sE2_sw + this%old_E_sw**2
 !!---       case(5)
 !!---         this%SP_xi1 = this%SP_xi1 + 1
 !!---         if(this%present_crank_level/=this%crank_valid_level) then
 !!---           this%scrank_level = this%scrank_level + 1
 !!---         endif
 !!---         !this%sxi1 = this%sxi1 + this%present_xi(1)
 !!---         !this%sxi12 = this%sxi12 + this%present_xi(1)**2
 !!---         this%sE_xi1 = this%sE_xi1 + this%crank_eg(this%present_crank_level)
 !!---         this%sE2_xi1 = this%sE2_xi1 +  this%crank_eg(this%present_crank_level)**2
 !       case(6)
 !         this%SP_wn = this%SP_wn + 1
 !         this%sE_wf = this%sE_wf + this%old_E_wf
 !         this%sE2_wf = this%sE2_wf + this%old_E_wf**2
 !         this%w = this%w + SUM(this%wn) ! diff of the winding numbers of chi
 !         this%ws = this%ws + (SUM(this%wn))**2
 !         this%aw = this%aw + SUM(this%fopt%w2)
 !         this%awc = this%awc + (SUM(this%fopt%w2))**2
 !         this%SP_sw = this%SP_sw + 1
 !         this%sE_sw = this%sE_sw + this%old_E_sw
 !         this%sE2_sw = this%sE2_sw + this%old_E_sw**2
 !     end select
    else
      this%not_add_ensemble_count = this%not_add_ensemble_count + 1
    END IF
  END SUBROUTINE monte_carlo_2lay_ensemble_average

 ! FUNCTION monte_carlo_2lay_get_change_w(this) RESULT(r)
 !   CLASS(monte_carlo_2lay), INTENT(in) :: this
 !   INTEGER, DIMENSION(this%fopt%nc2)     :: r
 !   r(:) = this%fopt%w2(:)
 !   IF (this%change_mode == 2) THEN
 !     r(this%wma(this%cw(1))) = -r(this%wma(this%cw(1)))
 !   ELSE
 !     r(this%wma(this%cw(2))) = -r(this%wma(this%cw(2)))
 !     r(this%wma(this%cw(3))) = -r(this%wma(this%cw(3)))
 !   END IF
 ! END FUNCTION monte_carlo_2lay_get_change_w
  
  function monte_carlo_2lay_get_flip_w_b(this) result(r)
    class(monte_carlo_2lay), intent(in) :: this
    integer, dimension(this%fopt%nh)     :: r
    r(:) = this%fopt%w2(:)
    r(this%wma(this%cw(1))) = -r(this%wma(this%cw(1)))
  end function monte_carlo_2lay_get_flip_w_b
  
  function monte_carlo_2lay_get_swap_w_b(this) result(r)
    class(monte_carlo_2lay), intent(in) :: this
    integer, dimension(this%fopt%nh)     :: r
    r(:) = this%fopt%w2(:)
    r(this%wma(this%cw(2))) = -r(this%wma(this%cw(2)))
    r(this%wma(this%cw(3))) = -r(this%wma(this%cw(3)))
  end function monte_carlo_2lay_get_swap_w_b
  
  function monte_carlo_2lay_get_flip_w_s(this) result(r)
    class(monte_carlo_2lay), intent(in) :: this
    integer, dimension(this%n_pole)     :: r
    r(:) = this%w_s(:)
    !r(this%wma_s(this%cw(1))) = -r(this%wma_s(this%cw(1)))
    r(this%wma0_s(this%cw(1))) = -r(this%wma0_s(this%cw(1)))
  end function monte_carlo_2lay_get_flip_w_s
  
  function monte_carlo_2lay_get_swap_w_s(this) result(r)
    class(monte_carlo_2lay), intent(in) :: this
    integer, dimension(this%n_pole)     :: r
    r(:) = this%w_s(:)
    !r(this%wma_s(this%cw(2))) = -r(this%wma_s(this%cw(2)))
    !r(this%wma_s(this%cw(3))) = -r(this%wma_s(this%cw(3)))
    r(this%wma0_s(this%cw(2))) = -r(this%wma0_s(this%cw(2)))
    r(this%wma0_s(this%cw(3))) = -r(this%wma0_s(this%cw(3)))
  end function monte_carlo_2lay_get_swap_w_s

  function monte_carlo_2lay_get_annihilate_w_s(this) result(r)
    class(monte_carlo_2lay), intent(in) :: this
    integer, dimension(this%n_pole)     :: r
    r(:) = this%w_s(:)
    r(this%wma0_s(this%cw(2))) = 0
    r(this%wma0_s(this%cw(3))) = 0
  end function monte_carlo_2lay_get_annihilate_w_s
  
  function monte_carlo_2lay_get_create_w_s(this) result(r)
    class(monte_carlo_2lay), intent(in) :: this
    integer, dimension(this%n_pole)     :: r
    r(:) = this%w_s(:)
    r(this%wma0_s(this%cw(2))) = 1
    r(this%wma0_s(this%cw(3))) = -1
  end function monte_carlo_2lay_get_create_w_s
 
 ! FUNCTION monte_carlo_2lay_get_around_hole_pos(this, hole_pos) RESULT(r)
 !   CLASS(monte_carlo_2lay), INTENT(in)  :: this
 !   INTEGER, INTENT(in)  :: hole_pos
 !   INTEGER, DIMENSION(:), ALLOCATABLE :: r
 !   ALLOCATE (r(8), source=[ &
 !        & this%fopt%shift_site(7, hole_pos), this%fopt%shift_site(4, hole_pos), &
 !        & this%fopt%shift_site(8, hole_pos), this%fopt%shift_site(1, hole_pos), &
 !        & this%fopt%shift_site(5, hole_pos), this%fopt%shift_site(2, hole_pos), &
 !        & this%fopt%shift_site(6, hole_pos), this%fopt%shift_site(3, hole_pos)])
 ! END FUNCTION monte_carlo_2lay_get_around_hole_pos

  !SUBROUTINE monte_carlo_2lay_get_change_spin(this, trial_zeta)
  !  CLASS(monte_carlo_2lay), INTENT(in) :: this
  !  REAL(8), DIMENSION(this%fopt%pshape(1)*this%fopt%pshape(2)), INTENT(inout) :: trial_zeta
  !  REAL(8)                              :: x, dzeta
  !  INTEGER  :: i, site(8)
  !  CALL RANDOM_NUMBER(x)
  !  dzeta = pi*x
  !  site = this%get_around_hole_pos(this%fopt%hole(this%cw(4), 1))
  !  DO i = 1, SIZE(site)
  !    trial_zeta(site(i)) = dzeta
  !  ENDDO
  !END SUBROUTINE monte_carlo_2lay_get_change_spin

  !FUNCTION monte_carlo_2lay_get_new_wf(this, trial_zeta) RESULT(r)
  !  CLASS(monte_carlo_2lay), INTENT(in) :: this
  !  REAL(8), DIMENSION(this%fopt%pshape(1)*this%fopt%pshape(2)), INTENT(in) :: trial_zeta
  !  COMPLEX(8), DIMENSION(2*this%fopt%n, 2*this%fopt%n)     :: r
  !  INTEGER   :: i
  !  r = 0.d0
  !  DO i = 1, this%fopt%n
  !    r(2*i - 1, :) = this%Dm(i, :)*COS(trial_zeta(i)*0.5d0) - this%Dp(i, :)*SIN(trial_zeta(i)*0.5d0)
  !    r(2*i, :) = this%Dm(i, :)*SIN(trial_zeta(i)*0.5d0) + this%Dp(i, :)*COS(trial_zeta(i)*0.5d0)
  !  ENDDO
  !END FUNCTION monte_carlo_2lay_get_new_wf

  FUNCTION monte_carlo_2lay_filename_myid(this, bname, exe, id) RESULT(r)
    CLASS(monte_carlo_2lay), INTENT(in)        :: this
    CHARACTER(*), INTENT(in)                     :: bname, exe
    INTEGER, INTENT(in), OPTIONAL                 :: id
    CHARACTER(LEN(bname) + LEN(exe) + 4)              :: r
    CHARACTER(4)                                  :: myid
    INTEGER                                       :: i
    IF (PRESENT(id)) THEN
      i = id
    ELSE
      i = this%myid
    END IF
    WRITE (myid, '(i3.3, ".")') i
    r = bname//myid//exe
  END FUNCTION monte_carlo_2lay_filename_myid

  SUBROUTINE monte_carlo_2lay_avg_quantity(this, bname, exe)
    CLASS(monte_carlo_2lay), INTENT(inout)     :: this
    CHARACTER(*), INTENT(in)                     :: bname, exe
    REAL(8), DIMENSION(0:this%iter_max)             :: q, t
    REAL(8)                                       :: d
    INTEGER                                       :: i, j
    q = 0.d0
    DO i = 1, this%nsh
      OPEN (1000, file=this%filename_myid(bname, exe, i - 1))
      DO j = 0, this%iter_max
        READ (1000, *) t(j), d
        q(j) = q(j) + d
      END DO
      CLOSE (1000)
    END DO
    d = REAL(this%nsh, 8)
    OPEN (1000, file=bname//"."//exe)
    DO i = 0, this%iter_max
      WRITE (1000, *) t(i), q(i)/d
    END DO
    CLOSE (1000)
  END SUBROUTINE monte_carlo_2lay_avg_quantity

  !SUBROUTINE monte_carlo_2lay_avg_hist(this, bname, exe)
  !  CLASS(monte_carlo_2lay), INTENT(inout)     :: this
  !  CHARACTER(*), INTENT(in)                     :: bname, exe
  !  REAL(8), DIMENSION(this%hist_intervel_number)             :: q, t
  !  REAL(8)                                       :: d
  !  INTEGER                                       :: i, j
  !  q = 0.d0
  !  DO i = 1, this%nsh
  !    OPEN (1000, file=this%filename_myid(bname, exe, i - 1))
  !    DO j = 1, this%hist_intervel_number
  !      READ (1000, *) t(j), d
  !      q(j) = q(j) + d
  !    END DO
  !    CLOSE (1000)
  !  END DO
  !  d = REAL(this%nsh, 8)
  !  OPEN (1000, file=bname//"."//exe)
  !  DO i = 1, this%hist_intervel_number
  !    WRITE (1000, *) t(i), q(i)/d
  !  END DO
  !  CLOSE (1000)
  !END SUBROUTINE monte_carlo_2lay_avg_hist
  
  subroutine monte_carlo_2lay_check_min_E(this, new_E_tot, new_chi, new_w_b, new_w_s)
    class(monte_carlo_2lay), intent(inout)  :: this
    real(8), intent(in)                     :: new_E_tot
    real(8), dimension(:), intent(in)      :: new_chi
    integer, dimension(:), intent(in)       :: new_w_b, new_w_s
    integer :: i, id, ierr

    if(this%min_E_tot  > new_E_tot) then
      this%min_E_tot = new_E_tot
      this%min_w_b = new_w_b
      this%min_w_s = new_w_s
      this%min_chi = new_chi
      !print *, this%min_E_tot
      
      if(stop_found_min_E) then
        !call this%out_min_E_state()
        id = 1300 + this%myid
        open(id, file = "config/min_E_state.txt", status="replace")
        write(id, "(a,i3)") "myid : ", this%myid
        write(id, "(a,f20.10)") "E_min = ", this%min_E_tot
        write(id, "(a)") "w_b : "
        do i = 1, size(this%min_w_b)
          write(id, "(i3)", advance='no') this%min_w_b(i)
        end do
        write(id, "()")
        write(id, "(a)") "w_s : "
        do i = 1, size(this%min_w_s)
          write(id, "(i3)", advance='no') this%min_w_s(i)
        end do
        close(id)
        print *, "we found a state having energy smaller than initial state."
        print "(a,f20.15)", "itnitial : ", this%present_E_tot
        print "(a,f20.15)", "     new : ", new_E_tot
        print "(a,f20.15)", " delta E : ", new_E_tot - this%present_E_tot
        print *, "wrote/updated ""config/min_E_state.txt."
        print *, "please redo this program."
        call MPI_ABORT(MPI_COMM_WORLD, 0, ierr)
        !stop
      end if

    end if
    
  end subroutine monte_carlo_2lay_check_min_E
  
  subroutine monte_carlo_2lay_out_min_E_state(this)
    class(monte_carlo_2lay), intent(inout)  :: this
    real(8)   , allocatable               :: current(:,:)
    integer                               :: hole_chi(this%fopt%nh,2)
    real(8), dimension(:,:), allocatable      :: pole_s
    integer, dimension(:), allocatable        :: w_s
    integer i, id

    allocate(current(this%fopt%n, this%fopt%n)) 
    this%fopt%delta_chi = this%fopt%d_angval(this%min_chi)
    call this%fopt%get_Jp_matrix(current)
    call draw_current("current_min_" // string(this%myid) // ".dat", 1000 + this%myid, this%fopt%pshape, this%fopt%cpl, current)

    !call get_winding_point_surface(chi, clist, pole_s_conv, w_s_conv)
    hole_chi = this%present_hole_chi
    hole_chi(:,2) = this%min_w_b
    call write_format_of_chi_sb('chi0_min_' // string(this%myid) // '.plt',  1100 + this%myid, this%fopt%pshape, hole_chi, &
                            &     this%present_pole_s, this%min_w_s)
    call get_chi_winding_number_with_eta_winding_point_surface(this%eta, this%min_chi, this%fopt, pole_s, w_s)
    call write_format_of_chi_sb('chi_min_' // string(this%myid) // '.plt', 1200 + this%myid, this%fopt%pshape, hole_chi, pole_s, w_s)

    id = 1300 + this%myid
    open(id, file = "min_E_state_" // string(this%myid) // ".txt")
    write(id, "(a,i3)") "myid : ", this%myid
    write(id, "(a,f20.10)") "E_min = ", this%min_E_tot
    write(id, "(a)") "w_b : "
    do i = 1, size(this%min_w_b)
      write(id, "(i3)", advance='no') this%min_w_b(i)
    end do
    write(id, "()")
    write(id, "(a)") "w_s : "
    do i = 1, size(this%min_w_s)
      write(id, "(i3)", advance='no') this%min_w_s(i)
    end do
    close(id)
    
  end subroutine monte_carlo_2lay_out_min_E_state

  SUBROUTINE monte_carlo_2lay_result_mc(this)
    CLASS(monte_carlo_2lay), INTENT(inout)    :: this
    integer                                   :: i
    character(3)                                  :: num
    CALL this%avg_quantity("energy", "dat")
    CALL this%avg_quantity("specific_heat", "dat")
    !CALL this%avg_quantity("energy_sw", "dat")
    !CALL this%avg_quantity("specific_heat_sw", "dat")
    !CALL this%avg_quantity("energy_wf", "dat")
    !CALL this%avg_quantity("specific_heat_wf", "dat")
    !    call this%avg_quantity("windignumber", "dat")
    !CALL this%avg_quantity("avg_E_xi1", "dat")
    !CALL this%avg_quantity("dev_E_xi1", "dat")
    !CALL this%avg_quantity("nw", "dat")
    !CALL this%avg_quantity("wn", "dat")
    !CALL this%avg_quantity("nwc", "dat")
    !CALL this%avg_quantity("wnc", "dat")
    !CALL this%avg_quantity("Mab", "dat")
    !CALL this%avg_quantity("Mc", "dat")
    !CALL this%avg_quantity("Mhab", "dat")
    !CALL this%avg_quantity("Mhc", "dat")
    call this%avg_quantity("dwn", "dat")
    call this%avg_quantity("variance_dwn", "dat")
    call this%avg_quantity("dwn_b", "dat")
    call this%avg_quantity("variance_dwn_b", "dat")
    call this%avg_quantity("dwn_s", "dat")
    call this%avg_quantity("variance_dwn_s", "dat")
    
    call this%avg_quantity("awn", "dat")
    call this%avg_quantity("variance_awn", "dat")
    call this%avg_quantity("awn_b", "dat")
    call this%avg_quantity("variance_awn_b", "dat")
    call this%avg_quantity("awn_s", "dat")
    call this%avg_quantity("variance_awn_s", "dat")
    ! histgram
    !do i = 1, this%iter_max
    !  write(num, '(i3.3)') i
    !  call this%avg_hist("hist_xi1_"//num//"_", "dat")
    !enddo
  END SUBROUTINE monte_carlo_2lay_result_mc

  !! constructor
  ! function new_monte_carlo_2lay(pshape, hole, t, u, xi, chi, wf, myid) result(r)
  !FUNCTION new_monte_carlo_2lay(pshape, hole, t, u, jd, lm, xi, zeta, chi, wf, rs, crank_EV, crank_eg, myid) RESULT(r) !!
  function new_monte_carlo_2lay(pshape, hole, pole_s, w_s, t, u, jd, lm, mu, Bz, xi, chi, wf, eg, ne_u, ne_d, myid) result(r)
    !INTEGER, DIMENSION(2), INTENT(in)                   :: pshape
    !INTEGER, DIMENSION(:, :), INTENT(in)                :: hole
    !REAL(8), DIMENSION(:), INTENT(in)                   :: t
    !REAL(8), INTENT(in)                                 :: u, jd, lm
    !REAL(8), DIMENSION(pshape(1)*pshape(2)), INTENT(in) :: xi, chi, zeta
    !COMPLEX(8), DIMENSION(2*pshape(1)*pshape(2), 2*pshape(1)*pshape(2)), &
    !     & INTENT(in)                                        :: wf
    !INTEGER, INTENT(in), OPTIONAL                       :: myid
    !complex(8)                                               :: crank_EV(:)
    !real(8)                                               :: crank_eg(:)
    !TYPE(monte_carlo_2lay)                                :: r
    !CHARACTER(1)                                        :: rs
    integer, dimension(:, :), intent(in)                :: pshape
    integer, dimension(:, :), intent(in)                :: hole
    integer, dimension(:), intent(in)                   :: w_s
    integer, intent(in), optional                       :: myid
    real(8), dimension(:, :), intent(in)                :: pole_s
    !real(8), intent(in)                                 :: t(3), lm
    real(8), intent(in)                                 :: t(3), u, jd, lm, mu(2), Bz
    real(8), dimension(:), intent(in)                   :: xi, chi, ne_u, ne_d, eg
    complex(8), dimension(:, :), intent(in)             :: wf
    type(monte_carlo_2lay)                                :: r
    !!   call r%monte_carlo_2lay_init(pshape, hole, t, u, xi, chi, wf, myid)
    !CALL r%monte_carlo_2lay_init(pshape, hole, t, u, jd, lm, xi, zeta, chi, wf, 't', crank_EV, crank_eg, myid) !!
    !call r%monte_carlo_2lay_init(pshape, hole, pole_s, w_s, t, lm, xi, chi, wf, eg, ne_u, ne_d, myid)
    !call r%monte_carlo_2lay_init(pshape, hole, pole_s, w_s, t, u, jd, lm, mu, xi, chi, wf, eg, ne_u, ne_d, myid)
    call r%monte_carlo_2lay_init(pshape, hole, pole_s, w_s, t, u, jd, lm, mu, Bz, xi, chi, wf, eg, ne_u, ne_d, myid)
  END FUNCTION new_monte_carlo_2lay

  ! addition
END MODULE monte_carlo_2lay_mod
