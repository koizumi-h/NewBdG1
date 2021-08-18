program test
  use measurement_mod
  use math_mod
  use wave_function_mod
  use visual_mod
  use utility_mod
  use path_list_mod
  use hop_iterator_mod, only : get_hopping_path
  implicit none
  integer, parameter :: n = 256 !64
  integer, dimension(2,1) :: pshape
  complex(8), dimension(4*n,4*n) :: wf
  real(8), dimension(4*n) :: eg
  real(8) :: kbt, eg_range(2)
  integer :: i, j, r
  real(8) :: t, ne, e, de, point(2), radius
  real(8), dimension(:), allocatable :: rho
  integer :: n_sp
  character(len=256) :: form
  type(path_list), allocatable :: pl_1st, pl_2nd
  complex(8), dimension(:), allocatable :: delta
  pshape(1,1) = 16
  pshape(2,1) = 16
  pl_1st = get_hopping_path(pshape,[1,2])
  allocate(delta(pl_1st%length()))
  !read(21) wf ! mu = 0
  !read(22) eg ! mu = 0
  !!read(23) delta
  !read(31) wf ! 32*32
  !read(32) eg ! 32*32
  !read(41) wf !
  !read(42) eg !
  read(51) wf ! selfconsistent
  read(52) eg !
  read(53) delta !
  call sort_wf(eg,wf)
  call draw_delta("delta.dat", 105, pshape, pl_1st, delta)
  print *, "saved to ""delta.dat"""
  call draw_path_delta("path_delta.dat", 105, pshape, pl_1st, delta)
  print *, "saved to ""path_delta.dat"""
  !kbt = 0.15d0
  !kbt = 1d-1
  kbt = 0.02d0
  eg_range(1) = -4d0
  eg_range(2) = 4d0
  point(1) = 8d0
  point(2) = 8d0
  radius = 0.1d0!0.8d0
  n_sp = 10000
  !call save_stm_rho_vs_e_step("rho.dat", 20, n, wf, eg, kbt, eg_range, n_sp)
  call save_stm_point_rho_vs_e_step("rho.dat", 20, point, radius, pshape, n, wf, eg, kbt, eg_range, n_sp)
  print *, "saved to ""rho.dat"""
  
  open(25, file = "dos.dat")
  allocate(rho(n_sp+1))
  de = (eg_range(2) - eg_range(1)) / float(n_sp)
  do i=1, size(rho)
    e = eg_range(1) + de*(i-1)
    !do r = size(eg)/2 + 1, size(eg)
    !    rho(i) = rho(i) - dfermi(eg(r) - e, kbt) - dfermi(eg(r) + e, kbt)
    !end do
    do r = 1, size(eg)
        rho(i) = rho(i) - dfermi(eg(r) - e, kbt)
    end do
    write(25, "(2f30.20)") e, rho(i)
  end do
  deallocate(rho)
  close(25) 
  
  open(22, file = "ne.dat")
  write(form, "(a1,i,a2)") "(", 4*n + 2, "f)"
  !write(*,*) form
  allocate(rho(4*n))
  do i=1, n
    rho = 0d0
    do j=1, 4*n
      rho(j) = conjg(wf(4*i-3,j))*wf(4*i-3,j) + conjg(wf(4*i-2,j))*wf(4*i-2,j)
      rho(j) = rho(j) + conjg(wf(4*i-1,j))*wf(4*i-1,j) + conjg(wf(4*i,j))*wf(4*i,j)
    end do
    write(22,form) float(mod(i-1,pshape(1,1))+1), float((i-1)/pshape(1,1)+1), rho
    if(mod(i,pshape(1,1)) == 0) write(22, *) ""
  end do
  close(22)

  
  open(21, file = "eg.dat")
  do i=1, size(eg)
    write(21,*) i, eg(i)
  end do
  close(21)
  print *, "saved to ""eg.dat"""
  
  open(23, file = "eg2.dat")
  do i=1, size(eg)/2
    write(23,*) eg(i), eg(size(eg/2)-i+1), eg(i) + eg(size(eg/2)-i+1) 
  end do
  close(23)
  
  !n_sp = 100000
  !de = (eg_range(2) - eg_range(1)) / float(n_sp)
  !allocate(rho(n_sp+1))
  !call stm_rho_vs_e_step(n, wf, eg, kbt, eg_range, rho, n_sp)
  !open(22, file = "rhofermi.dat")
  !do i=1, size(rho)
  !  e = eg_range(1) + de*(i-1)
  !  write(22, "(2f30.20)") e, rho(i)*fermi(e, kbt)
  !end do
  !ne = 0d0
  !do i = 1, n_sp
  !  e = eg_range(1) + de*(i-1)
  !  ne = ne + (rho(i)*fermi(e,kbt) + rho(i+1)*fermi(e+de,kbt))*de*0.5d0
  !end do
  !print *,"ne : ",ne

  stop

end program test