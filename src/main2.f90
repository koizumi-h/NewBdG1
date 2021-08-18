PROGRAM main
  USE math_mod
  USE mpi
  !USE io_mod, ONLY:save_txt, load_txt, create_vortex, load_conf, load_conf_chi_w, load_mc_conf
  USE io_mod,                ONLY : save_txt, load_txt, load_conf_3D_2, load_mc_conf, load_min_E_state

  USE visual_mod, ONLY:draw_phase, write_format_of_xi, &
       &                              write_format_of_chi, print_phase
  !USE monte_carlo_zeta_mod2 !!
  USE monte_carlo_2lay_mod
  use diagonalize_mod,   only : diag, new_diag
  !USE energy_mod,            only : total_energy_sb_2lay_xichi_T0
  USE energy_mod,            only : total_energy_sb_2lay_xichi_T0_Aem
  !USE fchi_optimizer_mod
  USE fchi_optimizer_Aem_mod
  USE angular_variable_mod
  ! use conjugate_gradient_rsoc_mod, only : conjugate_gradient_rsoc,new_conjugate_gradient_rsoc
  IMPLICIT NONE
  INTEGER                                  :: nx, ny, n, nh, n_acc, n_pole, nz
  real(8), allocatable                      :: ne(:)
  INTEGER, DIMENSION(:,:), allocatable     :: pshape
  REAL(8), DIMENSION(3), PARAMETER         :: t = [1.d0, -0.12d0, 0.01d0] ! hopping coeff
  REAL(8), PARAMETER                       :: u = 8.d0 ! coulomb coeff
  !  real(8), parameter                       :: jd = 0.d0 ! across hole coeff
  REAL(8), PARAMETER                       :: jd = 0.5d0*4.0d0*(t(1)**2)/u ! across hole coeff
  REAL(8), PARAMETER                       :: lm = 0.02d0 ! rashba coeff
  REAL(8)                                  :: mu(2) ! chemical potential
  REAL(8), PARAMETER                       :: Bz = 0.01d0 ! rashba coeff
  LOGICAL, PARAMETER                        :: sflag = .TRUE. !slanting path
  !  logical,parameter                        :: sflag=.false.!no slanting path
  logical                                  :: optimized_cgr, optimized_E_tot, dummy_flag
  INTEGER, PARAMETER                       :: unit_num = 10
  INTEGER, DIMENSION(:, :), ALLOCATABLE    :: hole_xi, hole_chi
  REAL(8), DIMENSION(:), ALLOCATABLE       :: xi, chi, zeta, d_chi
  !  type(sdchi_optimizer),allocatable         :: opt
  TYPE(monte_carlo_2lay), ALLOCATABLE      :: mc2
  !COMPLEX(8), DIMENSION(:, :), ALLOCATABLE :: cpwf
  COMPLEX(8), DIMENSION(:, :), ALLOCATABLE :: wf
  !REAL(8), DIMENSION(:), ALLOCATABLE       :: cpeg, ranx
  REAL(8), DIMENSION(:), ALLOCATABLE       :: eg, ranx
  REAL(8)                                  :: p(1)
  INTEGER                                  :: i, ii, nrand
  REAL(8), DIMENSION(:), ALLOCATABLE       :: delta_xi
  REAL(8), DIMENSION(:), ALLOCATABLE       :: mcint
  INTEGER, DIMENSION(:), ALLOCATABLE       :: seed
  real(8)                                  :: total_E, total_E2, xi1, old_E_tot, new_E_tot
  INTEGER(kind=4)                          :: ierr, myid, nproc, status(MPI_STATUS_SIZE)
  !complex(8)                               :: EV(21)
  !real(8)                               :: crank_eg(21)
  integer t1, t2, t_rate, t_max, diff
  
  real(8), dimension(:,:), allocatable :: pole_s
  integer, dimension(:), allocatable :: w_s
  real(8), dimension(:), allocatable :: ne_u, ne_d
  type(diag), allocatable :: dg

  integer :: access
  integer :: file_exist
  type(angular_variable), allocatable :: an_chi
  !type(fchi_optimizer), allocatable :: fopt
  type(fchi_optimizer_Aem), allocatable :: fopt

  !-------------mpi start----------------
  CALL MPI_INIT(ierr)
  CALL MPI_COMM_SIZE(MPI_COMM_WORLD, nproc, ierr)
  CALL MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)

  IF (myid == 0) THEN
    call system_clock(t1) ! 開始時を記録
    WRITE (6, *) "nproc", nproc
    !      call load_conf('../conf/vconf', unit_num, pshape, hole1)
    !       call write_format_of_xi('xi.plt', unit_num, pshape, hole1)
    !       call write_format_of_chi('chi.plt', unit_num, pshape, hole2)

    !! load xi winding number
    !CALL load_conf('../config/vconf', unit_num, pshape, hole_xi)
    !! load chi winding number
    !CALL load_conf_chi_w('../config/vconf', unit_num, pshape, hole_chi, 7)
  
    CALL load_conf_3D_2('./config/vconf', unit_num, pshape, hole_xi,hole_chi, ne)
    nx = pshape(1, 1)
    ny = pshape(2, 1)
    nz = size(pshape(1,:))
    n = nx*ny*nz
    nh = SIZE(hole_xi, 1)
    !ne = n - nh
    n_acc = n - nh
    read(78) n_pole
    !open(90, file='../src/crank_EV.fbin', status='old', form='unformatted')
    !read(90) EV
    !close(90)
    !open(91, file='../src/crank_eg.fbin', status='old', form='unformatted')
    !read(91) crank_eg
  ENDIF
  IF (myid /= 0) ALLOCATE(pshape(2, 2))
  
  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) !==================================
  CALL MPI_BCAST(nx, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  CALL MPI_BCAST(ny, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  CALL MPI_BCAST(n, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  CALL MPI_BCAST(nh, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  CALL MPI_BCAST(n_acc, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  CALL MPI_BCAST(n_pole, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  ! broadcast parameters
  !CALL MPI_BCAST(pshape, 2, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  CALL MPI_BCAST(pshape, 4, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) !==================================

  !if (myid == 0) write (*, *) 'electron', nx, ny, n, nh, ne
  if (myid == 0) then 
    !write (*, *) 'electron', nx, ny, n, nh, ne, n_pole
    write (*, *) 'nx, ny, nz : ', nx, ny, nz
    write (*, *) 'n = nx*ny*nz, nh (small polaron) : ', n, nh
    write (*, *) 'ne(surface), ne(bulk) : ', ne
    write (*, *) 'n_pole (surface) : ', n_pole
  endif
  IF (myid /= 0) ALLOCATE (hole_xi(nh, 2), hole_chi(nh, 2))

  !ALLOCATE (xi(n))
  !ALLOCATE (chi(n))
  !ALLOCATE (cpwf(2*n, 2*n))
  !ALLOCATE (cpeg(2*n))
  !ALLOCATE (zeta(n), source=0.5d0*pi)
  allocate(xi(n))
  allocate(chi(n))
  allocate(wf(4*n_acc, 4*n_acc))
  allocate(eg(4*n_acc))
  allocate(ne_u(n))
  allocate(ne_d(n))
  allocate(pole_s(2, n_pole))
  allocate(w_s(n_pole))

 ! mwf2: IF (myid == 0) then
 !   read(70) wf
 !   read(71) eg
 !   read(72) xi
 !   read(73) chi
 !   read(74) ne_u
 !   read(75) ne_d
 !   read(76) pole_s
 !   read(77) w_s
 !   !READ (20) cpwf
 !   !READ (21) cpeg
 !   !READ (32) xi
 !   !READ (33) chi
 !   WRITE (*, *) 'end of read make wave func in wf2'
 !   !write (*, *) 'init E_tot', energy_total_energy(pshape, hole_chi, t, u, jd, lm, &
 !   !& xi, chi, cpwf, 'T')
 !   print *, 'E_tot(T=0)', total_energy_sb_2lay_xichi_T0(dg, wf, xi, chi)
 ! ENDIF mwf2

  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) !==================================
  !! broadcast values
  CALL MPI_BCAST(hole_xi, 2*nh, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  CALL MPI_BCAST(hole_chi, 2*nh, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  !CALL MPI_BCAST(xi, n, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
  !CALL MPI_BCAST(chi, n, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
  !CALL MPI_BCAST(cpeg, 2*n, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) !==================================
  do i = 0, nproc - 1
    if (myid == i) then
      !read (20) cpwf
      read(70) wf
      read(71) eg
      read(72) xi
      read(73) chi
      read(74) ne_u
      read(75) ne_d
      read(76) pole_s
      read(77) w_s
      read(79) mu
      file_exist = access("./config/min_E_state.txt"," ")
      !print *, file_exist
      if(file_exist == 0) then
        call load_min_E_state("./config/min_E_state.txt", 80, old_E_tot, hole_chi, w_s)
        print "(i3,x,a)", myid, ": ***** LOAD  min_E_state.txt *****"
        allocate(an_chi)
        an_chi = new_angular_variable_3D(pshape, hole_chi)
        chi = an_chi%get_angular_variable() ! bulk
        chi = chi + an_chi%get_atan_phase_on_xy(pole_s, w_s, 1) ! surface
        chi = an_chi%rebuilt_angval_chi(chi)
        allocate(fopt)
        !fopt = new_fchi_optimizer(pshape, hole_chi, [0d0, 0d0],  t, u, jd, lm, xi, chi, wf, eg, ne_u, ne_d) ! ne is not used in fchi(monte-carlo)
        fopt = new_fchi_optimizer_Aem(pshape, hole_chi, [0d0, 0d0],  t, u, jd, lm, Bz, xi, chi, wf, eg, ne_u, ne_d) ! ne is not used in fchi(monte-carlo)
        CALL fopt%oneshot()
        chi = an_chi%get_rebuilt_phase(fopt%delta_chi)
      else
        print "(i3,x,a)", myid, "min_E_state.txt does not exist. use fort.77."
        allocate(fopt)
        fopt = new_fchi_optimizer_Aem(pshape, hole_chi, [0d0, 0d0],  t, u, jd, lm, Bz, xi, chi, wf, eg, ne_u, ne_d)
      end if
        call fopt%calc_only_Aem_dphase()
        allocate(dg)
        dg = new_diag(pshape, hole_xi, t, u, jd, lm)
        !old_E_tot = total_energy_sb_2lay_xichi_T0(dg, wf, xi, chi)
        old_E_tot = total_energy_sb_2lay_xichi_T0_Aem(dg, wf, xi, chi, fopt%Aem_dphase)
        deallocate(dg)
      print "(i3,x,a,f)", myid, 'E_tot(T=0) = ', old_E_tot
    endif
    !CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) !==================================
  enddo
  !stop
  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) !==================================
  !CALL MPI_BCAST(cpwf, 4*n*n, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD, ierr)
  !WRITE (*, *) 'end MPI_BCAST(cpwf, 4*n**2, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD,ierr)'
  !CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) !==================================
  !write (6, *) myid, nx, ny, sum(xi), sum(chi), size(cpwf), size(cpeg)
  !write (*, *) myid, size(hole_chi), size(hole_xi)
  !CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) !==================================
  !CALL MPI_BCAST(EV, 21, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD, ierr)
  !write (*, *) myid, "EV", sum(EV)
  !CALL MPI_BCAST(crank_eg, 21, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
  !write (*, *) myid, "EV", sum(EV)
  !CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) !==================================

  !! Monte carlo caluculation---------
  WRITE (*, *) 'Start monte carlo'
  !CALL draw_phase('hfxiv.gp', unit_num, pshape, xi, hole_chi(:, 1))
  ALLOCATE(mc2)
  !WRITE (6, *) 'start mc'
  !mc2 = new_monte_carlo_zeta(pshape, hole_chi, t, u, jd, lm, xi, zeta, chi, cpwf, 'T', EV, crank_eg, myid)
  !mc2 = new_monte_carlo_2lay(pshape, hole_chi, pole_s, w_s, t, lm, xi, chi, wf, eg, ne_u, ne_d, myid)
  !mc2 = new_monte_carlo_2lay(pshape, hole_chi, pole_s, w_s, t, u, jd, lm, mu, xi, chi, wf, eg, ne_u, ne_d, myid)
  mc2 = new_monte_carlo_2lay(pshape, hole_chi, pole_s, w_s, t, u, jd, lm, mu, Bz, xi, chi, wf, eg, ne_u, ne_d, myid)
  WRITE (6, *) "load_mc_conf"
  WRITE (6, *) "load_mc_conf start"
  CALL load_mc_conf('./config/conf_monte_carlo', unit_num, mcint)
  mcint(6) = nproc
  WRITE (6, *) "load_mc_conf end"
  CALL mc2%set_mcint(mcint)

  CALL RANDOM_SEED(size=nrand)
  ALLOCATE (seed(nrand))
  ALLOCATE (ranx(nrand))
  CALL SYSTEM_CLOCK(count=ii)
  !CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) !==================================
  seed = ii + myid*123457
  !print*, "seed : ",seed, ii
  !stop
  CALL RANDOM_SEED(put=seed)
  CALL RANDOM_NUMBER(ranx)
  seed = INT(ranx*(2**30))
  CALL RANDOM_SEED(put=seed)
  !write (*,"(i3,x,a,x,i)"), myid, "seed : ",seed
  write (*,*), myid, "seed : ",seed
  !CALL save_txt("ranx.txt", REAL(INT(ranx*2**30), 8))
  WRITE (6, *) 'mcint', mcint

  !CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) !==================================
  !IF (myid == 0) THEN
  !   DO i = 1, 10
  !      CALL RANDOM_NUMBER(ranx)
  !   END DO
  !   CALL MPI_BCAST(ranx, nrand*nproc, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
  !END IF
!
!  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) !==================================

  !DO i = 1, nrand
  !   seed(i) = ranx((myid)*nrand+i)
  !END DO
  !CALL RANDOM_SEED(put=seed)
  !IF (myid == 0) CALL save_txt("ranx.txt", REAL(INT(ranx*2**30), 8))

  ! call mc%opt%solve
  !  write(6,*) "total ene", mc%get_total_energy()
  !  write(6,*) "kin ene", mc%opt%get_kinetic_energy()
  CALL mc2%calc("r")

  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) !==================================

  IF (myid == 0) then
    CALL mc2%result_mc()
    !  chi = bl%get_rebuilt_phase(mc%tau)
    !  call draw_phase('tchi.gp', unit_num, pshape, chi, hole2(:, 1))
    call system_clock(t2, t_rate, t_max) ! 終了時を記録
    if (t2 < t1) then
      diff = (t_max - t1) + t2 + 1
    else
      diff = t2 - t1
    endif
    print "(A, F10.3)", "time it took was:", diff/dble(t_rate)
  endif
  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) !==================================

 ! DEALLOCATE (mc2)

  !!----------------------------------

  !DEALLOCATE (xi)
  !DEALLOCATE (cpwf)
  !DEALLOCATE (cpeg)
  WRITE (*, *) 'end'
  CALL MPI_FINALIZE(ierr)
END PROGRAM main
