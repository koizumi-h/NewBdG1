program test
  use math_mod
  use path_list_mod
  use hop_iterator_mod
  implicit none
  real(8), parameter :: t(3) = [1d0, 0d0, 0d0]
  !real(8), parameter :: t(3) = [1d0, -0.12d0, 0d0]
  integer, parameter :: nx = 8, ny = 8, ne = nx*ny
  integer, parameter :: pshape(2,1) = [[nx,ny]]
  complex(8), dimension(:,:), allocatable :: ham, wf, wfk 
  integer :: jx, jy, j, k, ix, iy, i, l
  !real(8), dimension(2*nx*ny) :: eg
  real(8), dimension(:), allocatable :: eg
  real(8), dimension(:), allocatable :: kx, ky
  real(8) :: prod(2*nx*ny)
  character(len=128) :: form
  type(path_list), allocatable :: pl_1st, pl_2nd

  allocate(ham(2*nx*ny, 2*nx*ny))
  allocate(wf(2*nx*ny, 2*nx*ny))
  allocate(eg(2*nx*ny))
  allocate(wfk(2*(nx+1)*(ny+1), 2*nx*ny))
  allocate(kx(nx+1))
  allocate(ky(ny+1))

  pl_1st =  get_hopping_path(pshape, [1,2])
  pl_2nd =  get_hopping_path(pshape, [5,6])
  
  ham = 0d0
  do i = 1, pl_1st%length()
    j = pl_1st%value(i)%i
    k = pl_1st%value(i)%f
    !print *, i, j, k
    ham(2*j-1, 2*k-1) = -t(1)
    ham(2*j  , 2*k  ) = -t(1)
    ham(2*k-1, 2*j-1) = -t(1)
    ham(2*k  , 2*j  ) = -t(1)
  end do
  do i = 1, pl_2nd%length()
    j = pl_2nd%value(i)%i
    k = pl_2nd%value(i)%f
    ham(2*j-1, 2*k-1) = -t(2)
    ham(2*j  , 2*k  ) = -t(2)
    ham(2*k-1, 2*j-1) = -t(2)
    ham(2*k  , 2*j  ) = -t(2)
  end do
  CALL my_zheev('l', ham, eg, wf)
  
  !print *,eg
  open(100,file = "eg.dat")
  do j = 1, 2*nx*ny
    write(100,*) j, eg(j)
  end do
  close(100)
  print *,"saved ""eg.dat"""
  stop



  deallocate(ham)
  deallocate(wf)
  deallocate(eg)
  deallocate(wfk)
  deallocate(kx)
  deallocate(ky)
contains 

  function get_site_number(jx, jy, nx, ny) result(j)
    integer, intent(in) :: jx, jy, nx, ny
    integer             :: j
    j = jx + nx*(jy-1)
  end function get_site_number
end program test