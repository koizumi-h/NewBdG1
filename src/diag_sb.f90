program diag_sb
  use math_mod,              only : ui, pi
  USE io_mod,                ONLY : save_txt, load_txt, create_vortex, load_conf_3D_2 ,load_conf_3D, load_conf_chi_w
  USE visual_mod,            ONLY : draw_phase_3D, write_format_of_xi_3D, draw_path_3D, &
       &                              write_format_of_chi_3D, print_phase, draw_spin, draw_spin_normalize
  USE angular_variable_mod,  ONLY : angular_variable, new_angular_variable, new_angular_variable_3D
  use parameter_mod
  use path_list_mod, only : path_list
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
  delta(1:size(delta)/2) = 0.2d0 ! x-dir
  delta(size(delta)/2+1:size(delta)) = -0.2d0 ! y-dir

  ! set init mu
  !mu(1) = -0.05d0
  mu(1) = -0.1d0
  mu(1) = -0.2d0
  mu(2) = 4d0

  ! diagonalize 
  call dg%diag_sb_2lay(mu, delta, ne_u, ne_d, Sx, Sy, Sz)
  ! result >>> dg%eg, dg%wf
  
  call dg%calc_ne(ne_u, ne_d, kbt)
  call dg%calc_spin(Sx, Sy, Sz, kbt)

  print *, "bulk ne : ", sum(ne_u(dg%nl(1)+1:) + ne_d(dg%nl(1)+1:))
  print *, "surface ne : ", sum(ne_u(1:dg%nl(1)) + ne_d(1:dg%nl(1)))


  ! output
  !print *, "Sx" , Sx
  !print *, "Sy" , Sy
  !print *, "eg" , dg%eg
  !print *, "wf" , dg%wf
  !open(104, file = "spin.dat")
  !open(105, file = "spin_n.dat") ! normalized vector
  !do j = 1, dg%n_acc
  !  i = dg%ets(j)
  !  vec(1) = Sx(i)
  !  vec(2) = Sy(i)
  !  vec(3) = Sz(i)
  !  !print *, i, sqrt(vec(1)**2d0+vec(2)**2d0+vec(3)**2d0)
  !  rvec= site_to_r(i, pshape) - vec * 0.5d0
  !  write(104, "(6f)") rvec, vec
  !  if( vec(1) /= 0d0 .and. vec(2) /= 0d0 .and. vec(3)/=0d0) then
  !    vec = 0.8d0 * vec / sqrt(vec(1)**2d0+vec(2)**2d0+vec(3)**2d0)
  !  end if
  !  rvec= site_to_r(i, pshape) - vec * 0.5d0
  !  write(105, "(6f)") rvec, vec
  !end do
  !close(104)
  !close(105)
  call draw_spin("spin.dat", 104, pshape, Sx, Sy, Sz)
  call draw_spin_normalize("spin_n.dat", 105, pshape, Sx, Sy, Sz)
  print *, "saved to ""spin.dat"" "
  print *, "saved to ""spin_n.dat"" "

  open(106,file = "sb_eg.dat")
  k = 1
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

  stop "stop diag_sb"
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !                                   !
  !    Below, for fchi_optimizer      !
  !                                   !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !!!! calculation using xi and compare
  allocate(eg1(4*dg%n_acc))
  allocate(eg2(4*dg%n_acc))
  allocate(teg(4*dg%n_acc))
  allocate(twf(4*dg%n_acc, 4*dg%n_acc))
  allocate(tdelta(size(delta)))
  teg = dg%eg
  twf = dg%wf
  tdelta = delta
  
  !ne_u = 0.5d0
  !ne_d = 0.5d0
  call dg%calc_xi(xi, Sx, Sy)
  !xi(1:dg%nl(1)) = 0d0 ! surface
  xi = an_xi%rebuilt_angval_xi(xi)
  !xi = 0d0
  call draw_phase_3D('sb_init_xi.dat', unit_num, pshape, xi, hole_xi(:, 1))
  !Sx = 0d0
  !Sy = 0d0
  !Sz = 0d0
  !do i = dg%n_start(2), dg%n_end(2) ! bulk
  !  Sx(i) = 0.5d0*cos(xi(i))
  !  Sy(i) = 0.5d0*sin(xi(i))
  !end do

  ! previous diagonalize
  !call dg%calc_delta(delta, 0.01d0)
  call dg%diag_sb_2lay(mu, delta, ne_u, ne_d, Sx, Sy, Sz)
  eg2 = dg%eg
  call dg%calc_delta(delta, 0.01d0)
  print *, "mean(abs(delta)) no : ", sum(abs(delta))/size(delta)
  
!  ! calc xi test
!  dg%eg = teg
!  dg%wf = twf
!  !call dg%calc_delta(delta, 0.01d0)
!  !print *, "mean(abs(delta)) tn : ", sum(abs(delta))/size(delta)
!  do i = 1, dg%n_acc
!    is = dg%ets(i)
!    dg%wf(4*i-3, :) = exp(-0.5d0*ui*xi(is))*dg%wf(4*i-3, :) 
!    dg%wf(4*i-2, :) = exp( 0.5d0*ui*xi(is))*dg%wf(4*i-2, :) 
!    dg%wf(4*i-1, :) = exp( 0.5d0*ui*xi(is))*dg%wf(4*i-1, :) 
!    dg%wf(4*i  , :) = exp(-0.5d0*ui*xi(is))*dg%wf(4*i  , :) 
!    
! !   dg%wf(4*i-3, :) = exp( 0.5d0*ui*xi(is))*dg%wf(4*i-3, :) 
! !   dg%wf(4*i-2, :) = exp(-0.5d0*ui*xi(is))*dg%wf(4*i-2, :) 
! !   dg%wf(4*i-1, :) = exp(-0.5d0*ui*xi(is))*dg%wf(4*i-1, :) 
! !   dg%wf(4*i  , :) = exp( 0.5d0*ui*xi(is))*dg%wf(4*i  , :) 
!  end do
!  !call dg%calc_delta_xi(delta, xi, 0.01d0)
!  !print *, "mean(abs(delta)) tx : ", sum(abs(delta))/size(delta)
  
  ! diagonalize using xi
  dg%eg = teg
  dg%wf = twf
  delta = tdelta
  !call dg%calc_delta_xi(delta, xi, 0.01d0)
  !call dg%calc_delta(delta, 0.01d0)
  !delta = delta *2d0
  !delta(1:size(delta)/2) = 0.2d0 ! x-dir
  !delta(size(delta)/2+1:size(delta)) = -0.2d0 ! y-dir
  !call dg%calc_ne(ne_u, ne_d, 0.01d0)
  !ne_u = 0d0
  !ne_d = 0d0
  call dg%diag_sb_2lay_xi(mu, delta, ne_u, ne_d, Sx, Sy, Sz, xi)
  eg1 = dg%eg
  
  
  !! repeat
  !call dg%calc_xi(xi, Sx, Sy)
  !xi = an_xi%rebuilt_angval_xi(xi)
  !do i = 1, 10
  !  call dg%diag_sb_2lay_xi(mu, delta, ne_u, ne_d, Sx, Sy, Sz, xi)
  !  !dg%eg = eg2
  !  call dg%calc_ne(nne_u, nne_d, 0.01d0)
  !  call dg%calc_spin(nSx, nSy, nSz)
  !  call dg%calc_delta(ndelta, 0.01d0)
  !  Sx = (1d0-ratio)*Sx + ratio*nSx
  !  Sy = (1d0-ratio)*Sy + ratio*nSy
  !  Sz = (1d0-ratio)*Sz + ratio*nSz
  !  ne_u = (1d0-ratio)*ne_u + ratio*nne_u
  !  ne_d = (1d0-ratio)*ne_d + ratio*nne_d
  !  delta = (1d0-ratio)*delta + ratio*ndelta
  !  call dg%calc_xi(xi, Sx, Sy)
  !  xi = an_xi%rebuilt_angval_xi(xi)
  !  print *, "repeat : ", i
  !  print *, "bulk ne : ", sum(ne_u(dg%nl(1)+1:) + ne_d(dg%nl(1)+1:))
  !  print *, "surface ne : ", sum(ne_u(1:dg%nl(1)) + ne_d(1:dg%nl(1)))
  !end do
  
  !xi(1:dg%nl(1)) = 0d0 ! surface
  
!  do i = 1, dg%n_acc
!    is = dg%ets(i)
!    dg%wf(4*i-3, :) = exp(-0.5d0*ui*xi(is))*dg%wf(4*i-3, :) 
!    dg%wf(4*i-2, :) = exp( 0.5d0*ui*xi(is))*dg%wf(4*i-2, :) 
!    dg%wf(4*i-1, :) = exp( 0.5d0*ui*xi(is))*dg%wf(4*i-1, :) 
!    dg%wf(4*i  , :) = exp(-0.5d0*ui*xi(is))*dg%wf(4*i  , :) 
!    
! !   dg%wf(4*i-3, :) = exp( 0.5d0*ui*xi(is))*dg%wf(4*i-3, :) 
! !   dg%wf(4*i-2, :) = exp(-0.5d0*ui*xi(is))*dg%wf(4*i-2, :) 
! !   dg%wf(4*i-1, :) = exp(-0.5d0*ui*xi(is))*dg%wf(4*i-1, :) 
! !   dg%wf(4*i  , :) = exp( 0.5d0*ui*xi(is))*dg%wf(4*i  , :) 
!  end do
!  do i = 1, 10
!    do j = 1, 10
!      print *, i, j, dot_product(dg%wf(:,i), twf(:,j))
!      !print *, i, j, dot_product(dg%wf(:,i), dg%wf(:,j))
!      !print *, i, j, dot_product(twf(:,i), twf(:,j))
!    end do
!  end do
!  stop

  open(110,file = "sb_eg_xi.dat")
  k = 1
  do i = 1, 4*dg%n_acc
      write(110, *) i, eg1(i), eg2(i)
  end do
  close(110)
  print *, "saved to ""sb_eg_xi.dat"" "

  !!delta = tdelta
  !call calc_eg_sb_2lay_xi(dg, mu, delta, ne_u, ne_d, Sx, Sy, Sz, xi)
  !!call sort_wf(dg%eg, dg%wf)
  !open(110,file = "sb_eg_xi2.dat")
  !k = 1
  !do i = 1, 4*dg%n_acc
  !    write(110, *) i, dg%eg(i), eg2(i)
  !end do
  !close(110)

  !dg%eg = eg1
  !call dg%calc_delta_xi(tdelta, xi, kbt)
  !print *, "mean(abs(delta)) xi : ", sum(abs(tdelta))/size(tdelta)

  !call dg%calc_ne(ne_u, ne_d, kbt)
  !call dg%calc_spin_xi(Sx, Sy, Sz, xi, kbt)
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
  !write(56) tdelta
  write(57) Sx
  write(58) Sy
  write(59) Sz

  stop
end program diag_sb