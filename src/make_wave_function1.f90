program make_wave_function1
  use math_mod,              only : ui, pi
  USE io_mod,                ONLY : save_txt, load_txt, create_vortex, load_conf_3D_2 ,load_conf_3D, load_conf_chi_w
  USE visual_mod,            ONLY : draw_phase_3D, write_format_of_xi_3D, draw_path_3D, &
       &                              write_format_of_chi_3D, print_phase, draw_spin, draw_spin_normalize
  USE angular_variable_mod,  ONLY : angular_variable, new_angular_variable, new_angular_variable_3D
  use parameter_mod
  use path_list_mod, only : path_list, path
  use car_parrinello_mod, only: car_parrinello, new_car_parrinello
  use diagonalize_mod,   only : diag, new_diag
 
  use energy_mod
  !use diagonalize_mod,   only : calc_eg_sb_2lay_xi
  use utility_mod
  !$ use omp_lib
  implicit none

  integer                                  :: nx, ny, nz, n, n_acc, nh , r(3) , time_s,time_e,countpersec
  REAL(8), DIMENSION(:),allocatable        :: ne
  INTEGER, DIMENSION(:,:),allocatable      :: pshape
  REAL(8), DIMENSION(3), PARAMETER         :: t = [1d0, -0.12d0, 0.01d0]  ! hopping coeff[1st, 2nd, layer]
  !REAL(8), DIMENSION(3), PARAMETER         :: t = [1d0, 0d0, 0.1d0]  ! hopping coeff[1st, 2nd, layer]
  !REAL(8), DIMENSION(3), PARAMETER         :: t = [1d0, -0.12d0, 0d0]  ! hopping coeff[1st, 2nd, layer]
  !REAL(8), DIMENSION(3), PARAMETER         :: t = [1d0, 0d0, 0d0]  ! hopping coeff[1st, 2nd, layer]
  !REAL(8), DIMENSION(3), PARAMETER         :: t = [0d0, 0d0, 0d0]  ! hopping coeff[1st, 2nd, layer]
  REAL(8), PARAMETER                       :: u = 8.d0         ! coulomb coeff
  !REAL(8), PARAMETER                       :: u = 1000.d0         ! coulomb coeff
  REAL(8), PARAMETER                       :: jd = 0.5d0*4.0d0*(t(1)**2)/u ! across hole coeff
  !REAL(8), PARAMETER                       :: jd = 0.0d0 ! across hole coeff
  REAL(8), PARAMETER                       :: lm = 0.02d0 ! rashba coeff
  real(8), parameter                       :: kbt = 0.01d0 ! temperature
  INTEGER, PARAMETER                       :: unit_num = 10
  INTEGER, DIMENSION(:, :), ALLOCATABLE    :: hole_xi, hole_chi !dim is (holes,2)
  REAL(8), DIMENSION(:), ALLOCATABLE       :: xi,chi, Sx, Sy, Sz, dummy
  TYPE(angular_variable), ALLOCATABLE      :: an_xi,an_chi
  !real(8), dimension(:), allocatable :: ne_u, ne_d
  !complex(8), dimension(:), allocatable :: delta
  real(8), dimension(:), allocatable       :: mu
  !real(8), dimension(:), allocatable       :: eg1, eg2, teg
  !complex(8), dimension(:,:), allocatable    ::twf
  !complex(8), dimension(:), allocatable    ::tdelta
  type(car_parrinello), allocatable        :: cp
  type(diag), allocatable :: dg
  
  INTEGER                                  :: i, j, k, is
  !real(8)                                   :: vec(3), rvec(3)
  
  !type(path)                               :: p
  !integer t1, t2, t_rate, t_max, diff, ndo  
  !complex(8), dimension(:,:), allocatable    ::testwf
  !real(8) :: a, b

  print *, "make wave function 1"
  
  ! load hole position, xi winding number and chi winding number
  CALL load_conf_3D_2('./config/vconf', unit_num, pshape, hole_xi,hole_chi, ne)
  ! output xi 3D plot data
  call write_format_of_xi_3D('xi.plt', unit_num, pshape, hole_xi)

  nz = size(pshape(1,:))
  nh = size(hole_xi(:,1))
  n = site_max(pshape)
  n_acc = n - nh

  ALLOCATE(xi(n))
  ALLOCATE(chi(n),source=0.d0)
  !ALLOCATE(Sx(n))
  !ALLOCATE(Sy(n))
  !ALLOCATE(Sz(n))
  !ALLOCATE(ne_u(n))
  !ALLOCATE(ne_d(n))
  ALLOCATE(mu(nz))

  ! set init xi  
  allocate(an_xi)
  an_xi = new_angular_variable_3D( pshape, hole_xi)
  xi = an_xi%get_af_angular_variable()
  xi = an_xi%rebuilt_angval_xi(xi)
  call draw_phase_3D('init_xi.dat', unit_num, pshape, xi, hole_xi(:, 1))
  
  ! set init mu
  !mu(1) = -0.05d0
  !mu(1) = -0.2d0
  !mu(1) = -0.1d0
  !mu(1) = -0.0775d0
  mu(1) = -0.2d0
  !mu(1) = -0.23d0
  mu(2) = 4d0

  allocate(cp)
  cp = new_car_parrinello(pshape, hole_xi, ne, t, u, jd, lm, kbt, mu, xi, chi) 
  
  allocate(dg)
  dg = new_diag(pshape, hole_xi, t, u, jd, lm)
  
  !! calc & set initial parameters
  !call dg%diag_sb_2lay(cp%mu, cp%delta, cp%ne_u, cp%ne_d, cp%Sx, cp%Sy, cp%Sz)
  !call cp%set_initial_wf(dg%wf, dg%eg)
  
  ! calc initial wf
  call dg%diag_sb_2lay(cp%mu, cp%delta, cp%ne_u, cp%ne_d, cp%Sx, cp%Sy, cp%Sz) ! by LAPACK
  !call cp%calc_only_wf() ! by car-parrinello with constant parameters
  !read(101) cp%wf ! by other program
  
  call cp%set_initial_value(wf = dg%wf)
  !call cp%set_initial_value(wf = wf)
  !call cp%calc_eigen_energy()
  
  ! run self-consistent
  call cp%run()
  

  open(106,file = "cp_eg.dat")
  do i = 1, 4*cp%n_acc
    write(106, *) i, cp%eg(i)
  end do
  close(106)
  print *, "saved to ""cp_eg.dat"" "
  
  call dg%calc_xi(xi, cp%Sx, cp%Sy)
  xi = an_xi%rebuilt_angval_xi(xi)
  call draw_phase_3D('cp_xi.dat', unit_num, pshape, xi, hole_xi(:, 1))
  call draw_spin("cp_spin.dat", 112, pshape, cp%Sx, cp%Sy, cp%Sz)
  call draw_spin_normalize("cp_spin_n.dat", 113, pshape, cp%Sx, cp%Sy, cp%Sz)
  
  !stop "stop before output fort.*"
  
  write(100) cp%eg
  write(101) cp%wf

  print *, "end: mwf1(before single shot)"

  !write(50) dg%wf
  !write(51) dg%eg
  !write(52) xi
  !!write(53) chi
  !write(54) ne_u
  !write(55) ne_d
  !write(56) delta
  !write(57) Sx, Sy, Sz
  !!write(58) Sy
  !!write(59) Sz
  !write(58) mu

end program make_wave_function1