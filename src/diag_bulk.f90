program diag_sb
  USE io_mod,                ONLY : save_txt, load_txt, create_vortex, load_conf_3D_2 ,load_conf_3D, load_conf_chi_w
  USE visual_mod,            ONLY : draw_phase_3D, write_format_of_xi_3D, draw_path_3D, &
       &                              write_format_of_chi_3D, print_phase
  USE angular_variable_mod,  ONLY : angular_variable, new_angular_variable, new_angular_variable_3D
  use parameter_mod
  use path_list_mod, only : path_list
  use diagonalize_mod,   only : diag, new_diag

  !use diagonalize_mod,   only : calc_eg_sb_2lay_xi

  implicit none

  integer                                  :: nx, ny, nz, n, n_acc, nh , r(3) , time_s,time_e,countpersec
  REAL(8), DIMENSION(:),allocatable        :: ne
  INTEGER, DIMENSION(:,:),allocatable      :: pshape
  REAL(8), DIMENSION(3), PARAMETER         :: t = [1d0, -0.d0, 0.01d0]  ! hopping coeff[1st, 2nd, layer]
  !REAL(8), DIMENSION(3), PARAMETER         :: t = [1d0, 0d0, 0.1d0]  ! hopping coeff[1st, 2nd, layer]
  !REAL(8), DIMENSION(3), PARAMETER         :: t = [1d0, -0.12d0, 0d0]  ! hopping coeff[1st, 2nd, layer]
  !REAL(8), DIMENSION(3), PARAMETER         :: t = [1d0, 0d0, 0d0]  ! hopping coeff[1st, 2nd, layer]
  !REAL(8), DIMENSION(3), PARAMETER         :: t = [0d0, 0d0, 0d0]  ! hopping coeff[1st, 2nd, layer]
  REAL(8), PARAMETER                       :: u = 8.d0         ! coulomb coeff
  REAL(8), PARAMETER                       :: jd = 0.5d0*4.0d0*(t(1)**2)/u ! across hole coeff
  !REAL(8), PARAMETER                       :: jd = 0.0d0 ! across hole coeff
  REAL(8), PARAMETER                       :: lm = 0.02d0 ! rashba coeff
  INTEGER, PARAMETER                       :: unit_num = 10
  INTEGER, DIMENSION(:, :), ALLOCATABLE    :: hole_xi, hole_chi !dim is (holes,2)
  REAL(8), DIMENSION(:), ALLOCATABLE       :: xi,chi, Sx, Sy, Sz, dummy
  TYPE(angular_variable), ALLOCATABLE      :: an_xi,an_chi
  type(diag), allocatable :: dg
  real(8), dimension(:), allocatable :: ne_u, ne_d
  !real(8), dimension(:), allocatable       :: mu
  real(8)      :: mu
  real(8), dimension(:), allocatable       :: eg1, eg2, teg
  complex(8), dimension(:,:), allocatable    ::twf
  
  INTEGER                                  :: i, j, k
  real(8)                                   :: vec(3), rvec(3)
  
  print *, "diagnalization for bulk"
  
  ! load hole position, xi winding number and chi winding number
  CALL load_conf_3D_2('./config/vconf_b', unit_num, pshape, hole_xi,hole_chi, ne)
  ! output xi 3D plot data
  call write_format_of_xi_3D('xi.plt', unit_num, pshape, hole_xi)

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
  !ALLOCATE(mu(nz))
  
  allocate(dg)
  dg = new_diag(pshape, hole_xi, t, u, jd, lm)


  ! set init xi  
  allocate(an_xi)
  an_xi = new_angular_variable_3D(pshape, hole_xi)
  xi = an_xi%get_af_angular_variable()
  xi = an_xi%rebuilt_angval_xi(xi)

  print *, "*********************************"
  print *, "load : fort.20"
  read(22) xi

  ! set init S
  Sx = 0d0
  Sy = 0d0
  Sz = 0d0
  do i = 1, dg%n
    Sx(i) = 0.5d0*cos(xi(i))
    Sy(i) = 0.5d0*sin(xi(i))
  end do

  ! set init ne_u, ne_d
  ne_u = 0.5d0
  ne_d = 0.5d0

  !! set init delta
  !delta(1:size(delta)/2) = 0.2d0 ! x-dir
  !delta(size(delta)/2+1:size(delta)) = -0.2d0 ! y-dir

  ! set init mu
  mu = 0d0

  ! diagonalize 
  call dg%diag_bulk(mu, Sx, Sy, Sz)
  ! result >>> dg%eg, dg%wf

  call dg%calc_ne(ne_u, ne_d, 0.01d0)
  call dg%calc_spin(Sx, Sy, Sz, 0.01d0)

  do i = 1, 1
    call dg%diag_bulk(mu, Sx, Sy, Sz)
    call dg%calc_spin(Sx, Sy, Sz, 0.01d0)
    print *, i, sum(Sx), sum(Sy)
  end do
  
  ! output
  print *, "bulk ne : ", sum(ne_u(:) + ne_d(:))
  !print *, "Sx" , Sx
  !print *, "Sy" , Sy
  !print *, "eg" , dg%eg
  !print *, "wf" , dg%wf
  open(104, file = "spin.dat")
  open(105, file = "spin_n.dat") ! normalized vector
  do j = 1, dg%n_acc
    i = dg%ets(j)
    vec(1) = Sx(i)
    vec(2) = Sy(i)
    vec(3) = Sz(i)
    !print *, i, sqrt(vec(1)**2d0+vec(2)**2d0+vec(3)**2d0)
    rvec= site_to_r(i, pshape) - vec * 0.5d0
    write(104, "(6f)") rvec, vec
    if( vec(1) /= 0d0 .and. vec(2) /= 0d0 .and. vec(3)/=0d0) then
      vec = 0.8d0 * vec / sqrt(vec(1)**2d0+vec(2)**2d0+vec(3)**2d0)
    end if
    rvec= site_to_r(i, pshape) - vec * 0.5d0
    write(105, "(6f)") rvec, vec
  end do
  close(104)
  close(105)
  print *, "saved to ""spin.dat"" "
  print *, "saved to ""spin_n.dat"" "

  open(106,file = "bulk_eg.dat")
  k = 1
  do i = 1, 4*dg%n_acc
      write(106, *) i, dg%eg(i)
  end do
  close(106)
  print *, "saved to ""bulk_eg.dat"" "

  !Sx = 0d0
  !Sy = 0d0
  !Sz = 0d0
  !do i = 1, dg%n
  !  Sx(i) = 0.5d0*cos(xi(i))
  !  Sy(i) = 0.5d0*sin(xi(i))
  !end do
  call draw_phase_3D('bulk_init_xi.dat', unit_num, pshape, xi, hole_xi(:, 1))
  print *, 'bulk_init_xi.dat'
 
  call dg%diag_bulk_xi(mu, Sx, Sy, Sz, xi)
  call dg%calc_ne(ne_u, ne_d, 0.01d0)
  call dg%calc_spin_xi(Sx, Sy, Sz, xi, 0.01d0)
  print *, "bulk ne xi: ", sum(ne_u(:) + ne_d(:))
  open(106,file = "bulk_eg_xi.dat")
  k = 1
  do i = 1, 4*dg%n_acc
      write(106, *) i, dg%eg(i)
  end do
  close(106)
  print *, "saved to ""bulk_eg_xi.dat"" "

  open(104, file = "spin_xi.dat")
  open(105, file = "spin_n_xi.dat") ! normalized vector
  do j = 1, dg%n_acc
    i = dg%ets(j)
    vec(1) = Sx(i)
    vec(2) = Sy(i)
    vec(3) = Sz(i)
    !print *, i, sqrt(vec(1)**2d0+vec(2)**2d0+vec(3)**2d0)
    rvec= site_to_r(i, pshape) - vec * 0.5d0
    write(104, "(6f)") rvec, vec
    if( vec(1) /= 0d0 .and. vec(2) /= 0d0 .and. vec(3)/=0d0) then
      vec = 0.8d0 * vec / sqrt(vec(1)**2d0+vec(2)**2d0+vec(3)**2d0)
    end if
    rvec= site_to_r(i, pshape) - vec * 0.5d0
    write(105, "(6f)") rvec, vec
  end do
  close(104)
  close(105)
  print *, "saved to ""spin_xi.dat"" "
  print *, "saved to ""spin_n_xi.dat"" "

  
  !open(110,file = "sb_eg_xi.dat")
  !k = 1
  !do i = 1, 4*dg%n_acc
  !    write(110, *) i, eg1(i), eg2(i)
  !end do
  !close(110)
  !print *, "saved to ""sb_eg_xi.dat"" "
  !call dg%calc_ne(ne_u, ne_d, 0.01d0)
  !call dg%calc_spin_xi(Sx, Sy, Sz, xi)
  !print *, "bulk ne : ", sum(ne_u(dg%nl(1)+1:) + ne_d(dg%nl(1)+1:))
  !print *, "surface ne : ", sum(ne_u(1:dg%nl(1)) + ne_d(1:dg%nl(1)))
  !open(114, file = "spin_xi.dat")
  !open(115, file = "spin_n_xi.dat") ! normalized vector
  !do j = 1, dg%n_acc
  !  i = dg%ets(j)
  !  vec(1) = Sx(i)
  !  vec(2) = Sy(i)
  !  vec(3) = Sz(i)
  !  rvec= site_to_r(i, pshape) - vec * 0.5d0
  !  !print *, i, sqrt(vec(1)**2d0+vec(2)**2d0+vec(3)**2d0)
  !  write(114, "(6f)") rvec, vec
  !  if( vec(1) /= 0d0 .and. vec(2) /= 0d0 .and. vec(3)/=0d0) then
  !    vec = vec / sqrt(vec(1)**2d0+vec(2)**2d0+vec(3)**2d0)
  !  end if
  !  rvec= site_to_r(i, pshape) - vec * 0.5d0
  !  write(115, "(6f)") rvec, vec
  !end do
  !close(114)
  !close(115)
  !print *, "saved to ""spin_xi.dat"" "
  !print *, "saved to ""spin_n_xi.dat"" "

  stop
end program diag_sb