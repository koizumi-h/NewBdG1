program diag_test4
  use math_mod, only : ui, pi, my_zheev, fermi, dfermi
  use path_list_mod, only : path_list
  use hop_iterator_mod, only : get_hopping_path
  
  implicit none
  !real(8), parameter :: t(3) = [1d0, 0d0, 0d0]
  real(8), parameter :: t(3) = [1d0, -0.12d0, 0d0]
  integer, parameter :: nx = 48, ny = 48, ne = nx*ny
  integer, parameter :: pshape(2,1) = [[nx,ny]]
  real(8), parameter :: delta = 0.2d0, mu = -1.0d0
  complex(8), dimension(:,:), allocatable :: ham, wf
  integer :: jx, jy, j, k, ix, iy, i, l
  !real(8), dimension(2*nx*ny) :: eg
  real(8), dimension(:), allocatable :: eg
  real(8) :: prod(2*nx*ny)
  character(len=128) :: form
  type(path_list), allocatable :: pl_1st, pl_2nd
  type(path_list), allocatable :: pl_1x, pl_1y
  real(8) :: e, kbt(10), dos(10), rho(10), eg_range(2)
  integer :: n_step, site

  n_step = 5000
  eg_range(1) = -4d0
  eg_range(2) =  4d0
  do i = 1, 10
    kbt(i) = 0.01d0*i
  end do
  allocate(ham(4*nx*ny, 4*nx*ny))
  allocate(wf(4*nx*ny, 4*nx*ny))
  allocate(eg(4*nx*ny))

  pl_1st =  get_hopping_path(pshape, [1,2])
  pl_2nd =  get_hopping_path(pshape, [5,6])
  pl_1x =  get_hopping_path(pshape, [1])
  pl_1y =  get_hopping_path(pshape, [2])
  
  ham = 0d0
  do i = 1, nx*ny
    ham(4*i-3, 4*i-3) = -mu
    ham(4*i-2, 4*i-2) = -mu
    ham(4*i-1, 4*i-1) =  mu
    ham(4*i  , 4*i  ) =  mu
  end do
  do i = 1, pl_1st%length()
    j = pl_1st%value(i)%i
    k = pl_1st%value(i)%f
    ham(4*j-3, 4*k-3) = -t(1)
    ham(4*j-2, 4*k-2) = -t(1)
    ham(4*j-1, 4*k-1) =  t(1)
    ham(4*j  , 4*k  ) =  t(1)
    ham(4*k-3, 4*j-3) = -t(1)
    ham(4*k-2, 4*j-2) = -t(1)
    ham(4*k-1, 4*j-1) =  t(1)
    ham(4*k  , 4*j  ) =  t(1)
  end do
  do i = 1, pl_2nd%length()
    j = pl_2nd%value(i)%i
    k = pl_2nd%value(i)%f
    ham(4*j-3, 4*k-3) = -t(2)
    ham(4*j-2, 4*k-2) = -t(2)
    ham(4*j-1, 4*k-1) =  t(2)
    ham(4*j  , 4*k  ) =  t(2)
    ham(4*k-3, 4*j-3) = -t(2)
    ham(4*k-2, 4*j-2) = -t(2)
    ham(4*k-1, 4*j-1) =  t(2)
    ham(4*k  , 4*j  ) =  t(2)
  end do
  do i = 1, pl_1x%length()
    j = pl_1x%value(i)%i
    k = pl_1x%value(i)%f
    ham(4*j-3, 4*k  ) = delta
    ham(4*j-2, 4*k-1) = delta
    ham(4*j-1, 4*k-2) = delta
    ham(4*j  , 4*k-3) = delta
    ham(4*k-3, 4*j  ) = delta
    ham(4*k-2, 4*j-1) = delta
    ham(4*k-1, 4*j-2) = delta
    ham(4*k  , 4*j-3) = delta
  end do
  do i = 1, pl_1y%length()
    j = pl_1y%value(i)%i
    k = pl_1y%value(i)%f
    ham(4*j-3, 4*k  ) = -delta
    ham(4*j-2, 4*k-1) = -delta
    ham(4*j-1, 4*k-2) = -delta
    ham(4*j  , 4*k-3) = -delta
    ham(4*k-3, 4*j  ) = -delta
    ham(4*k-2, 4*j-1) = -delta
    ham(4*k-1, 4*j-2) = -delta
    ham(4*k  , 4*j-3) = -delta
  end do
  
  call my_zheev("l", ham, eg, wf)

  !print *,eg
  open(100,file = "eg.dat")
  do j = 1,4*nx*ny
    write(100,*) j, eg(j)
  end do
  close(100)
  print *,"saved to ""eg.dat"""
  !open(101,file = "dos.dat")
  !do i = 1, n_step
  !  e = eg_range(1) + (eg_range(2)-eg_range(1))/float(n_step)*i
  !  dos = 0d0
  !  do j = 1,4*nx*ny
  !    do k = 1, 10
  !      dos(k) = dos(k) - dfermi(eg(j) - e, kbt(k))
  !    end do
  !  end do
  !  write(101,"(11f)") e, dos
  !end do
  !close(101)
  !print *,"saved to ""dos.dat"""
  open(102,file = "rho.dat")
  do i = 1, n_step
    e = eg_range(1) + (eg_range(2)-eg_range(1))/float(n_step)*i
    rho = 0d0
    do j = 2*nx*ny+1,4*nx*ny
      site = nx/2 + nx*(ny/2-1)
      !do site = 1, nx*ny
        do k = 1, 10
          rho(k) = rho(k) - conjg(wf(4*site-3, j))*wf(4*site-3, j)*dfermi(eg(j) - e, kbt(k))
          rho(k) = rho(k) - conjg(wf(4*site-2, j))*wf(4*site-2, j)*dfermi(eg(j) - e, kbt(k))
          rho(k) = rho(k) - conjg(wf(4*site-1, j))*wf(4*site-1, j)*dfermi(eg(j) + e, kbt(k))
          rho(k) = rho(k) - conjg(wf(4*site  , j))*wf(4*site  , j)*dfermi(eg(j) + e, kbt(k))
        end do
      !end do
    end do
    write(102,"(11f)") e, rho
  end do
  close(102)
  print *,"saved to ""rho.dat"""

  deallocate(ham)
  deallocate(wf)
  deallocate(eg)
  stop
end program diag_test4