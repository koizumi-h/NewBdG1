program diag_sb_sc
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
  REAL(8), DIMENSION(3), PARAMETER         :: t = [1d0, -0.12d0, 0.1d0]  ! hopping coeff[1st, 2nd, layer]
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
  real(8), dimension(:), allocatable       :: eg1, eg2, teg
  complex(8), dimension(:,:), allocatable    ::twf
  complex(8), dimension(:), allocatable    ::tdelta
  
  INTEGER                                  :: i, j, k, is
  real(8)                                   :: vec(3), rvec(3)
  
  real(8), dimension(:), allocatable       :: nne_u, nne_d, nSx, nSy, nSz
  real(8), dimension(:), allocatable       :: oSx, oSy, oSz, oeg
  logical                                  :: conv
  complex(8), dimension(:), allocatable    :: ndelta
  REAL(8), PARAMETER                       :: ratio = 0.2d0 &
                                          & , tol = 1d-5
  type(path)                               :: p

  print *, "diagnalization for 2 layers of surface and bulk"
  
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
  ALLOCATE(Sx(n))
  ALLOCATE(Sy(n))
  ALLOCATE(Sz(n))
  ALLOCATE(ne_u(n))
  ALLOCATE(ne_d(n))
  ALLOCATE(mu(nz))

  allocate(dg)
  dg = new_diag(pshape, hole_xi, t, u, jd, lm)

  ALLOCATE(delta(dg%pl_1st_l(1)%length()))
  

  ! set init xi  
  allocate(an_xi)
  an_xi = new_angular_variable_3D( pshape,hole_xi)
  xi = an_xi%get_af_angular_variable()
  xi(1:dg%nl(1)) = 0d0 ! surface
  xi = an_xi%rebuilt_angval_xi(xi)

  ! set init S
  Sx = 0d0
  Sy = 0d0
  Sz = 0d0
  do i = dg%n_start(2), dg%n_end(2) ! bulk
    Sx(i) = 0.5d0*cos(xi(i))
    Sy(i) = 0.5d0*sin(xi(i))
  end do

  ! set init ne_u, ne_d
  ne_u = 0.5d0
  ne_d = 0.5d0

  ! set init delta
  !delta(1:size(delta)/2) = 0.2d0 ! x-dir
  !delta(size(delta)/2+1:size(delta)) = -0.2d0 ! y-dir
  do i = 1, dg%pl_1st_l(1)%length()
    p = dg%pl_1st_l(1)%index(i)
    if(dg%y(p%i) == dg%y(p%f)) then
      delta(i) = 0.2d0 ! x-dir
    else
      delta(i) = -0.2d0 ! y-dir
    end if
  end do

  ! set init mu
  !mu(1) = -0.05d0
  mu(1) = -0.1d0
  !mu(1) = -0.23d0
  mu(2) = 4d0

  ! diagonalize 
  call dg%diag_sb_2lay(mu, delta, ne_u, ne_d, Sx, Sy, Sz)
  ! result >>> dg%eg, dg%wf
  
  call dg%calc_ne(ne_u, ne_d, kbt)
  call dg%calc_spin(Sx, Sy, Sz, kbt)
  call dg%calc_delta(delta, kbt)

  print *, "bulk ne : ", sum(ne_u(dg%nl(1)+1:) + ne_d(dg%nl(1)+1:))
  print *, "surface ne : ", sum(ne_u(1:dg%nl(1)) + ne_d(1:dg%nl(1)))

  !! repeat ( self-consistent-like )
  ALLOCATE(nSx(n), nSy(n), nSz(n))
  ALLOCATE(oSx(n), oSy(n), oSz(n))
  ALLOCATE(nne_u(n), nne_d(n))
  ALLOCATE(ndelta(dg%pl_1st_l(1)%length()))
  ALLOCATE(oeg(4*dg%n_acc))
  oSx = Sx
  oSy = Sy
  oSz = Sz
  oeg = dg%eg
  do i = 1, 1000000
    call dg%diag_sb_2lay(mu, delta, ne_u, ne_d, Sx, Sy, Sz)
    call dg%calc_ne(nne_u, nne_d, kbt)
    call dg%calc_spin(nSx, nSy, nSz, kbt)
    call dg%calc_delta(ndelta, kbt)
    Sx = (1d0-ratio)*Sx + ratio*nSx
    Sy = (1d0-ratio)*Sy + ratio*nSy
    Sz = (1d0-ratio)*Sz + ratio*nSz
    ne_u = (1d0-ratio)*ne_u + ratio*nne_u
    ne_d = (1d0-ratio)*ne_d + ratio*nne_d
    delta = (1d0-ratio)*delta + ratio*ndelta
    print *, "repeat : ", i
    print *, "bulk ne : ", sum(ne_u(dg%nl(1)+1:) + ne_d(dg%nl(1)+1:))
    print *, "surface ne : ", sum(ne_u(1:dg%nl(1)) + ne_d(1:dg%nl(1)))
    print *, "max(abs(d eg))", maxval(abs(dg%eg(:)-oeg(:))/(abs(oeg(:))+1d-8))
    print *, "max(abs(d Sx, d Sy))", maxval(abs(Sx(:)-oSx(:))/(abs(oSx(:))+1d-8)) &
                                &, maxval(abs(Sy(:)-oSy(:))/(abs(oSy(:))+1d-8))
    ! judge of converged
    conv = .true.
    do k = 1, n
        if(abs(Sx(k)-oSx(k))/(abs(oSx(k))+1d-8) > tol) conv = .false.
        if(abs(Sy(k)-oSy(k))/(abs(oSy(k))+1d-8) > tol) conv = .false.
    end do
    do k = 1, 4*dg%n_acc
      if(abs(dg%eg(k)-oeg(k))/(abs(oeg(k)) + 1d-8) > tol) conv = .false.
    end do
    if(conv) exit
    oSx = Sx
    oSy = Sy
    oSz = Sz
    oeg = dg%eg
  end do

  ! output
  !print *, "Sx" , Sx
  !print *, "Sy" , Sy
  !print *, "eg" , dg%eg
  !print *, "wf" , dg%wf
  
  call draw_spin("spin.dat", 104, pshape, Sx, Sy, Sz)
  call draw_spin_normalize("spin_n.dat", 105, pshape, Sx, Sy, Sz)
  print *, "saved to ""spin.dat"" "
  print *, "saved to ""spin_n.dat"" "

  open(106,file = "sb_eg.dat")
  do i = 1, 4*dg%n_acc
    write(106, *) i, dg%eg(i)
  end do
  close(106)
  print *, "saved to ""sb_eg.dat"" "

  open(107, file = "delta.txt")
  do i = 1, size(delta)
    write(107, *) delta(i)
  end do
  close(107)
  print *, "saved to ""delta.txt"" "
  
  open(108, file = "ne.dat")
  do i = 1, n
    if(site_to_x(i, pshape) == 1) write(108, "()")
    if(site_to_x(i, pshape) == 1 .and. site_to_y(i, pshape) == 1) write(108, "()")
    write(108, "(3i5,2f20.15)") site_to_r(i, pshape), ne_u(i), ne_d(i)
  end do
  close(108)
  print *, "saved to ""ne.dat"" "
  
  write(100) dg%eg
  write(101) dg%wf

  !stop "stop diag_sb_sc"
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !                                   !
  !    Below, for fchi_optimizer      !
  !                                   !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  call dg%calc_xi(xi, Sx, Sy)
  !xi(1:dg%nl(1)) = 0d0 ! surface
  xi = an_xi%rebuilt_angval_xi(xi)
  !xi = 0d0

  call draw_phase_3D('sb_init_xi.dat', unit_num, pshape, xi, hole_xi(:, 1))
  call dg%diag_sb_2lay_xi(mu, delta, ne_u, ne_d, Sx, Sy, Sz, xi)


  !call calc_eg_sb_2lay_xi(dg, mu, delta, ne_u, ne_d, Sx, Sy, Sz, xi)
  !call sort_wf(dg%eg, dg%wf)
  open(110,file = "sb_eg_xi.dat")
  k = 1
  do i = 1, 4*dg%n_acc
      write(110, *) i, dg%eg(i)
  end do
  close(110)

  !call dg%calc_delta_xi(delta, xi, kbt)
  !print *, "mean(abs(delta)) xi : ", sum(abs(delta))/size(delta)

  call dg%calc_ne(ne_u, ne_d, kbt)
  call dg%calc_spin_xi(Sx, Sy, Sz, xi, kbt)
  print *, "bulk ne : ", sum(ne_u(dg%nl(1)+1:) + ne_d(dg%nl(1)+1:))
  print *, "surface ne : ", sum(ne_u(1:dg%nl(1)) + ne_d(1:dg%nl(1)))
  
  call draw_spin("spin_xi.dat", 114, pshape, Sx, Sy, Sz)
  call draw_spin_normalize("spin_n_xi.dat", 115, pshape, Sx, Sy, Sz)
  print *, "saved to ""spin_xi.dat"" "
  print *, "saved to ""spin_n_xi.dat"" "

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

  stop
end program diag_sb_sc