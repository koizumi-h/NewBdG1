program diag_single_shot
  use math_mod,              only : ui, pi
  USE io_mod,                ONLY : save_txt, load_txt, create_vortex, load_conf_3D_2 ,load_conf_3D, load_conf_chi_w
  USE visual_mod,            ONLY : draw_phase_3D, write_format_of_xi_3D, draw_path_3D, &
       &                              write_format_of_chi_3D, print_phase, draw_spin, draw_spin_normalize
  USE angular_variable_mod,  ONLY : angular_variable, new_angular_variable, new_angular_variable_3D
  use parameter_mod
  use path_list_mod, only : path_list, path
  use diagonalize_mod,   only : diag, new_diag
  use utility_mod, only : sort_wf

  !use diagonalize_mod,   only : calc_eg_sb_2lay_xi

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
  REAL(8), PARAMETER                       :: jd = 0.5d0*4.0d0*(t(1)**2)/u ! across hole coeff
  !REAL(8), PARAMETER                       :: jd = 0.0d0 ! across hole coeff
  REAL(8), PARAMETER                       :: lm = 0.02d0 ! rashba coeff
  real(8), parameter                       :: kbt = 0.01d0 ! temperature
  INTEGER, PARAMETER                       :: unit_num = 10
  INTEGER, DIMENSION(:, :), ALLOCATABLE    :: hole_xi, hole_chi !dim is (holes,2)
  REAL(8), DIMENSION(:), ALLOCATABLE       :: xi,chi, Sx, Sy, Sz, dummy
  TYPE(angular_variable), ALLOCATABLE      :: an_xi,an_chi
  type(diag), allocatable :: dg
  real(8), dimension(:), allocatable :: ne_u, ne_d
  complex(8), dimension(:), allocatable :: delta
  real(8), dimension(:), allocatable       :: mu
  !real(8), dimension(:), allocatable       :: eg1, eg2, teg
  complex(8), dimension(:,:), allocatable    ::wf1
  !complex(8), dimension(:), allocatable    ::tdelta
  
  !real(8), dimension(:), allocatable       :: nne_u, nne_d, nSx, nSy, nSz
  !real(8), dimension(:), allocatable       :: oSx, oSy, oSz, oeg
  !logical                                  :: conv
  !complex(8), dimension(:), allocatable    :: ndelta
  !REAL(8), PARAMETER                       :: ratio = 0.2d0 &
  !                                        & , tol = 1d-5
  !type(path)                               :: p
  integer :: i, j

  print *, "start : only single shot"
  
  ! load hole position, xi winding number and chi winding number
  CALL load_conf_3D_2('./config/vconf', unit_num, pshape, hole_xi,hole_chi, ne)
  !! output xi 3D plot data
  !call write_format_of_xi_3D('xi.plt', unit_num, pshape, hole_xi)

  nz = size(pshape(1,:))
  nh = size(hole_xi(:,1))
  n = site_max(pshape)
  n_acc = n - nh

  ALLOCATE(xi(n))
  ALLOCATE(chi(n),source=0.d0)
  ALLOCATE(Sx(n))
  ALLOCATE(Sy(n))
  ALLOCATE(Sz(n))
  ALLOCATE(ne_u(n))
  ALLOCATE(ne_d(n))
  ALLOCATE(mu(nz))
  
  allocate(an_xi)
  an_xi = new_angular_variable_3D( pshape,hole_xi)

  allocate(dg)
  dg = new_diag(pshape, hole_xi, t, u, jd, lm)

  ALLOCATE(delta(dg%pl_1st_l(1)%length()))

  !!!  single shot (xiの取り出し) をする前の wf, eg を入力
  read(100) dg%eg
  read(101) dg%wf

  !!!上記を計算したときと同じ mu, kbt を設定する必要あり。 kbt : line 28
  !!!もしくはsingle shot前後で電子数が合うmuに調整する必要あり?
  !!!(single shot で大きく電子数が変わってしまうため)
  !mu(1) = -0.1d0
  mu(1) = -0.2d0
  !mu(1) = -0.0775d0
  mu(2) = 4d0
  
  call dg%calc_ne(ne_u, ne_d, kbt)
  call dg%calc_spin(Sx, Sy, Sz, kbt)
  call dg%calc_delta(delta, kbt)
  
  !call draw_spin("spin_before.dat", 112, pshape, Sx, Sy, Sz)
  !call draw_spin_normalize("spin_n_before.dat", 113, pshape, Sx, Sy, Sz)
  !print *, ""
  !print *, "save : ""spin_before.dat"", ""spin_n_before.dat"" "

  print *, ""
  print *, "before single shot"
  print *, "bulk ne : ", sum(ne_u(dg%nl(1)+1:) + ne_d(dg%nl(1)+1:))
  print *, "surface ne : ", sum(ne_u(1:dg%nl(1)) + ne_d(1:dg%nl(1)))

  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !                                   !
  !    Below, for fchi_optimizer      !
  !           single shot             !
  !                                   !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  call dg%calc_xi(xi, Sx, Sy)
  xi = an_xi%rebuilt_angval_xi(xi)

  !call draw_phase_3D('sb_init_xi.dat', unit_num, pshape, xi, hole_xi(:, 1))

  !!! single shot 
  call dg%diag_sb_2lay_xi(mu, delta, ne_u, ne_d, Sx, Sy, Sz, xi)

  open(110,file = "ss_eg_xi.dat")
  do i = 1, 4*dg%n_acc
      write(110, *) i, dg%eg(i)
  end do
  close(110)

  !call dg%calc_delta_xi(delta, xi, kbt)
  !print *, "mean(abs(delta)) xi : ", sum(abs(delta))/size(delta)

  call dg%calc_ne(ne_u, ne_d, kbt)
  call dg%calc_spin_xi(Sx, Sy, Sz, xi, kbt)
  print *, ""
  print *, "after single shot"
  print *, "bulk ne : ", sum(ne_u(dg%nl(1)+1:) + ne_d(dg%nl(1)+1:))
  print *, "surface ne : ", sum(ne_u(1:dg%nl(1)) + ne_d(1:dg%nl(1)))
  
  call draw_spin("ss_spin.dat", 114, pshape, Sx, Sy, Sz)
  call draw_spin_normalize("ss_spin_n.dat", 115, pshape, Sx, Sy, Sz)
  print *, ""
  print *, "save : ""ss_spin.dat"", ""ss_spin_n.dat"" "

  ! for mwf2 (fchi optimizer)
  write(50) dg%wf
  write(51) dg%eg
  write(52) xi
  !write(53) chi
  write(54) ne_u
  write(55) ne_d
  write(56) delta
  write(57) Sx, Sy, Sz
  !write(58) Sy
  !write(59) Sz
  write(58) mu

  ALLOCATE(wf1(4*n_acc, 4*n_acc))
  do i = 1, dg%n_acc
    j = dg%ets(i)
    wf1(4*i-3, :) = exp(-0.5d0*ui*xi(j))*dg%wf(4*i-3, :)
    wf1(4*i-2, :) = exp( 0.5d0*ui*xi(j))*dg%wf(4*i-2, :)
    wf1(4*i-1, :) = exp( 0.5d0*ui*xi(j))*dg%wf(4*i-1, :)
    wf1(4*i  , :) = exp(-0.5d0*ui*xi(j))*dg%wf(4*i  , :)
  end do

  ! for arpes, stm, ... etc.
  write(110) dg%eg
  write(111) wf1   !xiの因子を含む波動関数(chiは含まない:chi=0)
  !write(112) dg%wf !xiの因子を含まない波動関数

  print *, "end : only single shot"
end program diag_single_shot