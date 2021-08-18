module io_mod
  implicit none
  interface save_txt
    module procedure &
      dsave_txt, zsave_txt, d1asave_txt, z1asave_txt, &
      i2asave_txt, d2asave_txt, z2asave_txt
  end interface
  interface load_txt
    module procedure &
      dload_txt, zload_txt, d1aload_txt, z1aload_txt, &
      i2aload_txt, d2aload_txt, z2aload_txt
  end interface
  !type position
  !  integer :: x,y,z
  !end type position
  type hole_param
    type(integer),dimension(:,:),allocatable :: p
  end type hole_param
contains ! * or '(e25.15)'
  subroutine dsave_txt(filename, v)
    character(*), intent(in)      :: filename
    real(8), intent(in)           :: v
    open(32, file = filename)
    write(32, *) v
    close(32)
  end subroutine dsave_txt
  subroutine dload_txt(filename, v)
    character(*), intent(in) :: filename
    real(8), intent(inout)   :: v
    open(32, file = filename)
    read(32, *) v
    close(32)
  end subroutine dload_txt
  subroutine zsave_txt(filename, v)
    character(*), intent(in)      :: filename
    complex(8), intent(in)        :: v
    open(32, file = filename)
    write(32, *) v
    close(32)
  end subroutine zsave_txt
  subroutine zload_txt(filename, v)
    character(*), intent(in)  :: filename
    complex(8), intent(inout) :: v
    open(32, file = filename)
    read(32, *) v
    close(32)
  end subroutine zload_txt
  subroutine d1asave_txt(filename, a, flag)
    character(*), intent(in)          :: filename
    real(8), dimension(:), intent(in) :: a
    logical, intent(in), optional     :: flag
    integer                           :: i
    open(31, file = filename)
    if (present(flag)) then
      do i = 1, size(a)
        write(31, *) i, a(i)
      end do
    else
      do i = 1, size(a)
        write(31, *) a(i)
      end do
    end if
    close(31)
  end subroutine d1asave_txt
  subroutine d1aload_txt(filename, a, flag)
    character(*), intent(in)             :: filename
    real(8), dimension(:), intent(inout) :: a
    logical, intent(in), optional        :: flag
    integer                              :: i, dummy
    open(31, file = filename)
    if (present(flag)) then
      do i = 1, size(a)
        read(31, *) dummy, a(i)
      end do
    else
      do i = 1, size(a)
        read(31, *) a(i)
      end do
    end if
    close(31)
  end subroutine d1aload_txt
  subroutine z1asave_txt(filename, a, flag)
    character(*), intent(in)             :: filename
    complex(8), dimension(:), intent(in) :: a
    logical, intent(in), optional        :: flag
    integer                              :: i
    open(31, file = filename)
    if (present(flag)) then
      do i = 1, size(a)
        write(31, *) i, a(i)
      end do
    else
      do i = 1, size(a)
        write(31, *) a(i)
      end do
    end if
    close(31)
  end subroutine z1asave_txt
  subroutine z1aload_txt(filename, a, flag)
    character(*), intent(in)                :: filename
    complex(8), dimension(:), intent(inout) :: a
    logical, intent(in), optional           :: flag
    integer                                 :: i, dummy
    open(31, file = filename)
    if (present(flag)) then
      do i = 1, size(a)
        read(31, *) dummy, a(i)
      end do
    else
      do i = 1, size(a)
        read(31, *) a(i)
      end do
    end if
    close(31)
  end subroutine z1aload_txt
  subroutine i2asave_txt(filename, a, flag)
    character(*), intent(in)             :: filename
    integer, dimension(:, :), intent(in) :: a
    logical, intent(in), optional        :: flag
    integer                              :: i, j
    open(31, file = filename)
    if (present(flag)) then
      do j = 1, size(a, 2)
        do i = 1, size(a, 1)
          write(31, *) i, j, a(i, j)
        end do
      end do
    else
      do j = 1, size(a, 2)
        do i = 1, size(a, 1)
          write(31, *) a(i, j)
        end do
      end do
    end if
    close(31)
  end subroutine i2asave_txt
  subroutine i2aload_txt(filename, a, flag)
    character(*), intent(in)                :: filename
    integer, dimension(:, :), intent(inout) :: a
    logical, intent(in), optional           :: flag
    integer                                 :: i, j, dummy
    open(31, file = filename)
    if (present(flag)) then
      do j = 1, size(a, 2)
        do i = 1, size(a, 1)
          read(31, *) dummy, dummy, a(i, j)
        end do
      end do
    else
      do j = 1, size(a, 2)
        do i = 1, size(a, 1)
          read(31, *) a(i, j)
        end do
      end do
    end if
    close(31)
  end subroutine i2aload_txt
  subroutine d2asave_txt(filename, a, flag)
    character(*), intent(in)             :: filename
    real(8), dimension(:, :), intent(in) :: a
    logical, intent(in), optional        :: flag
    integer                              :: i, j
    open(31, file = filename)
    if (present(flag)) then
      do j = 1, size(a, 2)
        do i = 1, size(a, 1)
          write(31, *) i, j, a(i, j)
        end do
      end do
    else
      do j = 1, size(a, 2)
        do i = 1, size(a, 1)
          write(31, *) a(i, j)
        end do
      end do
    end if
    close(31)
  end subroutine d2asave_txt
  subroutine d2aload_txt(filename, a, flag)
    character(*), intent(in)                :: filename
    real(8), dimension(:, :), intent(inout) :: a
    logical, intent(in), optional           :: flag
    integer                                 :: i, j, dummy
    open(31, file = filename)
    if (present(flag)) then
      do j = 1, size(a, 2)
        do i = 1, size(a, 1)
          read(31, *) dummy, dummy, a(i, j)
        end do
      end do
    else
      do j = 1, size(a, 2)
        do i = 1, size(a, 1)
          read(31, *) a(i, j)
        end do
      end do
    end if
    close(31)
  end subroutine d2aload_txt
  subroutine z2asave_txt(filename, a, flag)
    character(*), intent(in)                :: filename
    complex(8), dimension(:, :), intent(in) :: a
    logical, intent(in), optional           :: flag
    integer                                 :: i, j
    open(31, file = filename)
    if (present(flag)) then
      do j = 1, size(a, 2)
        do i = 1, size(a, 1)
          write(31, *) i, j, a(i, j)
        end do
      end do
    else
      do j = 1, size(a, 2)
        do i = 1, size(a, 1)
          write(31, *) a(i, j)
        end do
      end do
    end if
    close(31)
  end subroutine z2asave_txt
  subroutine z2aload_txt(filename, a, flag)
    character(*), intent(in)                   :: filename
    complex(8), dimension(:, :), intent(inout) :: a
    logical, intent(in), optional              :: flag
    integer                                    :: i, j, dummy
    open(31, file = filename)
    if (present(flag)) then
      do j = 1, size(a, 2)
        do i = 1, size(a, 1)
          read(31, *) dummy, dummy, a(i, j)
        end do
      end do
    else
      do j = 1, size(a, 2)
        do i = 1, size(a, 1)
          read(31, *) a(i, j)
        end do
      end do
    end if
    close(31)
  end subroutine z2aload_txt
  subroutine create_vortex(filename, id, confnum, plane, hole)
    character(*), intent(in)                           :: filename
    integer, intent(in)                                :: id, confnum
    integer, dimension(2)                              :: plane
    integer, allocatable, dimension(:, :), intent(out) :: hole
    integer, allocatable, dimension(:)                 :: pos_x, pos_y
    integer                                            :: nh, iost, &
    &                                                     dummy, i, j
    real(8)                                            :: dummyeg
    open(id, file = filename)
    do i = 1, confnum-1
      do j = 1, 4
        read(id, *)
      end do
    end do
    read(id, *) dummy, dummyeg
    nh = 0
    do
      read(id, '(300i3)', advance = 'no', iostat = iost) dummy
      if (iost < 0) exit
      nh = nh + 1
    end do
    rewind(id)
    allocate(hole(nh, 2), pos_x(nh), pos_y(nh))
    do i = 1, confnum-1
      do j = 1, 4
        read(id, *)
      end do
    end do
    read(id, *) dummy, dummyeg
    read(id, *) pos_x(1:nh)
    read(id, *) pos_y(1:nh)
    read(id, *) hole(1:nh, 2)
    forall(i = 1:nh) hole(i, 1) = pos_x(i) + (pos_y(i)-1) * plane(1)
    rewind(id)
    close(id)
    deallocate(pos_x, pos_y)
  end subroutine create_vortex
  subroutine load_conf(filename, id, pshape, hole)
    character(*), intent(in)                           :: filename
    integer, intent(in)                                :: id
    integer, dimension(2), intent(out)                 :: pshape
    integer, allocatable, dimension(:, :), intent(out) :: hole
    integer, allocatable, dimension(:)                 :: pos_x, pos_y
    integer                                            :: nh, iost, &
    &                                                     dummy, i
    open(id, file = filename)
    read(id, '()')
    read(id, *) pshape(1), pshape(2)
    read(id, '()')
    nh = 0
    do
      read(id, '(300i3)', advance = 'no', iostat = iost) dummy
      if (iost < 0) exit
      nh = nh + 1
    end do
    rewind(id)
    if(allocated(hole)) deallocate(hole)
    allocate(hole(nh, 2), pos_x(nh), pos_y(nh))
    read(id, '()')
    read(id, '()')
    read(id, '()')
    read(id, *) pos_x(1:nh)
    read(id, *) pos_y(1:nh)
    read(id, *) hole(1:nh, 2)
    forall(i = 1:nh) hole(i, 1) = pos_x(i) + (pos_y(i)-1) * pshape(1)
    rewind(id)
    close(id)
    deallocate(pos_x, pos_y)
  end subroutine load_conf
  !subroutine load_conf2(filename, id, pshape, hole_xi, hole_chi, j_ext)
  !  use fchi_optimizer_ex_mod
  !  character(*), intent(in)                           :: filename
  !  integer, intent(in)                                :: id
  !  integer, dimension(2), intent(out)                 :: pshape
  !  integer, allocatable, dimension(:, :), intent(out) :: hole_xi, hole_chi
  !  integer, allocatable, dimension(:)                 :: pos_x, pos_y
  !  integer                                            :: nh, nj_ext, iost, &
  !  &                                                     dummy, i
  !  type(type_j_ext), allocatable                      :: j_ext(:)
  !  open(id, file = filename)
  !  read(id, '()') ! skip line 1
  !  read(id, *) pshape(1), pshape(2) ! read line 2
  !  read(id, '()') ! skip line 3
  !  nh = 0
  !  do
  !    read(id, '(300i3)', advance = 'no', iostat = iost) dummy
  !    if (iost < 0) exit
  !    nh = nh + 1
  !  end do
  !  rewind(id)
  !  allocate(hole_xi(nh, 2), pos_x(nh), pos_y(nh))
  !  allocate(hole_chi(nh, 2))
  !  read(id, '()')
  !  read(id, '()')
  !  read(id, '()')
  !  read(id, *) pos_x(1:nh)
  !  read(id, *) pos_y(1:nh)
  !  read(id, *) hole_xi(1:nh, 2)
  !  read(id, *) hole_chi(1:nh, 2)
  !  forall(i = 1:nh) hole_xi(i, 1) = pos_x(i) + (pos_y(i)-1) * pshape(1)
  !  forall(i = 1:nh) hole_chi(i, 1) = pos_x(i) + (pos_y(i)-1) * pshape(1)
  !  read(id, '()')
  !  read(id, *) nj_ext
  !  allocate(j_ext(nj_ext))
  !  deallocate(pos_x, pos_y)
  !  allocate(pos_x(nj_ext), pos_y(nj_ext))
  !  read(id, '()')
  !  read(id, *) pos_x(1:nj_ext)
  !  read(id, *) pos_y(1:nj_ext)
  !  read(id, *) j_ext(1:nj_ext)%current
  !  forall(i = 1:nj_ext) j_ext(i)%pos = pos_x(i) + (pos_y(i)-1) * pshape(1)
  !  rewind(id)
  !  close(id)613         layers(j)%n = pshape(1,j)*pshape(2,j)↲
!614         layers(j)%nh = size(pos_h(j)%p(
  !  deallocate(pos_x, pos_y)
  !end subroutine load_conf2
  subroutine load_conf_chi_w(filename, id, pshape, hole, num)
    character(*), intent(in)                           :: filename
    integer, intent(in)                                :: id
    integer, intent(in)                                :: num ! yomikomitai gyousuu
    integer, dimension(2), intent(out)                 :: pshape
    integer, allocatable, dimension(:, :), intent(out) :: hole
    integer, allocatable, dimension(:)                 :: pos_x, pos_y
    integer                                            :: nh, iost, &
    &                                                     dummy, i
    open(id, file = filename)
    read(id, '()')
    read(id, *) pshape(1), pshape(2)
    read(id, '()')
    nh = 0
    do
      read(id, '(300i3)', advance = 'no', iostat = iost) dummy
      if (iost < 0) exit
      nh = nh + 1
    end do
    rewind(id)
    allocate(hole(nh, 2), pos_x(nh), pos_y(nh))
    read(id, '()')
    read(id, '()')
    read(id, '()')
    read(id, *) pos_x(1:nh)
    read(id, *) pos_y(1:nh)
    do i = 1, num-6
      read(id, '()')
    end do
    read(id, *) hole(1:nh, 2)
    forall(i = 1:nh) hole(i, 1) = pos_x(i) + (pos_y(i)-1) * pshape(1)
    rewind(id)
    close(id)
    deallocate(pos_x, pos_y)
  end subroutine load_conf_chi_w
  subroutine save_vortex(filename, hole)
    character(*), intent(in)             :: filename
    integer, dimension(:, :), intent(in) :: hole
    integer                              :: j
    open(32, file = filename)
    do j = 1, size(hole, 1)
      write(32, *) hole(j, :)
    end do
    close(32)
  end subroutine save_vortex
  subroutine load_vortex(filename, hole)
    character(*), intent(in)                :: filename
    integer, dimension(:, :), intent(inout) :: hole
    integer                                 :: j
    open(32, file = filename)
    do j = 1, size(hole, 1)
      read(32, *) hole(j, :)
    end do
    close(32)
  end subroutine load_vortex
  subroutine read_monte_carlo_dat(filepath, filename, nf, itermax, var, tempra)
    character(*), intent(in)                :: filename, filepath
    integer, intent(in)                     :: nf, itermax
    real(8), dimension(:, :), allocatable, intent(out)   :: var
    real(8), dimension(:), allocatable, intent(out), optional   :: tempra
    real(8), dimension(:), allocatable      :: temp
    character                               :: cfile*60
    integer                                 :: i, j
    allocate(var(itermax, nf))
    allocate(temp(itermax))
    do i = 1, nf
      write(cfile,'(i2.2)')i
      cfile = trim(filepath) // trim(cfile) // "/" // trim(filename)
      write(6,*) cfile
      open(32, file = cfile)
      do j = 1, itermax
        read(32, *) temp(j), var(j, i)
      end do
      close(32)
    end do
    if (present(tempra)) allocate(tempra(itermax), source=temp)
  end subroutine read_monte_carlo_dat
  subroutine load_mc_conf(filename, id, mcint)
  ! mcint  1:T_min, 2:T_max, 3:iter_max, 4:Titer_ne, 5:Titer_e, 6: nsh
    character(*), intent(in)                        :: filename
    integer, intent(in)                             :: id
    real(8), dimension(:), allocatable, intent(out) :: mcint
    if (.not. allocated(mcint)) allocate(mcint(6))
    if (size(mcint) /= 6) then
      write(6,*) "error: size mcint is not 6"
      stop
    end if
    open(id, file=filename)
    read(id, '()')
    read(id, *) mcint(1), mcint(2), mcint(3)
 !   read(id, *) T_min, T_max, iter_max
    read(id, '()')
    read(id, *) mcint(4), mcint(5)
 !   read(id, *) Titer_ne, Titer_e
    read(id, '()')
    read(id, *) mcint(6)
 !   read(id, *) nsh
    close(id)
  end subroutine load_mc_conf
  subroutine load_conf_3D_xyz(filename, id, pshape, hole_xi, hole_chi)
    character(*), intent(in)                           :: filename
    integer, intent(in)                                :: id
    integer, dimension(3), intent(out)                 :: pshape
    integer, allocatable, dimension(:, :), intent(out), optional :: hole_xi, hole_chi
    integer, allocatable, dimension(:)                 :: pos_x, pos_y, pos_z, pos
    integer                                            :: nh, iost, &
    &                                                     dummy, i
    open(id, file = filename)
    read(id, '()')      !skip a line "nx ny nz"
    read(id, *) pshape(1), pshape(2), pshape(3)
    read(id, '()')      !skip a line "4:hole_x, 5:..."
    nh = 0
    do
      read(id, '(300i3)', advance = 'no', iostat = iost) dummy
      if (iost < 0) exit
      nh = nh + 1
    end do
    rewind(id)
      allocate(pos_x(nh), pos_y(nh), pos_z(nh), pos(nh))
    if(present(hole_xi)) then
      if(allocated(hole_xi)) deallocate(hole_xi)
      allocate(hole_xi(nh, 2))
    end if
    if(present(hole_chi)) then
      if(allocated(hole_chi)) deallocate(hole_chi)
      allocate(hole_chi(nh, 2))
    end if
    
    
    read(id, '()')    !skip a line "nx ny nz"
    read(id, '()')    !skip a line nx ny nz
    read(id, '()')    !skip a line "4:hole_x, 5:..."
    read(id, *) pos_x(1:nh)
    read(id, *) pos_y(1:nh)
    read(id, *) pos_z(1:nh)
    if(present(hole_xi))then
      read(id, *) hole_xi(1:nh, 2)
    else
      read(id, '()')    !skip a line hole_xi
    end if
    if(present(hole_chi))then
      read(id, *) hole_chi(1:nh, 2)
    else
      read(id, '()')    !skip a line hole_chi
    end if
    do i = 1, nh
      pos(i) = pos_x(i) + (pos_y(i)-1) * pshape(1) + (pos_z(i)-1)*pshape(1)*pshape(2)
      !hole_xi(i, 1) = pos_x(i) + (pos_y(i)-1) * pshape(1)
      !hole_xi(i, 2) = pos_z(i)
    end do
    if(present(hole_xi))then
      hole_xi(:, 1) = pos
    end if
    if(present(hole_chi))then
      hole_chi(:, 1) = pos
    end if
    !rewind(id)
    close(id)
    deallocate(pos_x, pos_y, pos_z, pos)
  end subroutine load_conf_3D_xyz

  subroutine load_conf_3D(filename, id, pshape, hole_xi, hole_chi,layers)
    use parameter_mod
    character(*), intent(in)                           :: filename
    integer, intent(in)                                :: id
    integer, dimension(:,:), intent(out) ,allocatable, optional  :: pshape
    integer, allocatable, dimension(:, :), intent(out), optional :: hole_xi, hole_chi
    integer, allocatable, dimension(:)  ::  ns
    type(cuo2_plane), allocatable, dimension(:), intent(out), optional :: layers
    type(hole_param), allocatable, dimension(:)         :: pos_h
    integer                                            :: nh, nx, ny, nz,&
    &                                                     i, j, k, nh_all, nh_max
    character(10000)  buff

    open(id, file = filename)
    read(id, '()')      !skip a line "nz"
    read(id, *) nz
    if(allocated(pshape)) deallocate(pshape)
    allocate(pshape(2,nz))
    allocate(ns(nz))
    allocate(pos_h(nz))
    do j = 1, nz
      !write(*,*) "------------",j,"---------------"
      read(id, '()')      !skip a line "nx ny"
      read(id, *) nx, ny
      read(id, '()')      !skip a line "4:hole_x, 5:..."
      if(present(pshape)) then
        pshape(:,j) = [nx,ny]
      end if 
      nh = 0
      read(id,'(a)') buff
      !write(*,*) "buff",trim(buff) 
      if(len_trim(buff) >= 1) then
        if(scan(buff(1:1),"1234567890-") /= 0) then
          nh = nh + 1
        end if
        do i = 2 ,len_trim(buff)
          if(scan(buff(i-1:i-1)," ") /= 0 .and. scan(buff(i:i),"1234567890-") /= 0) then
            nh = nh + 1
          end if
        end do
      end if
      allocate(pos_h(j)%p(nh, 5),source=0)  ! x, y, xi, chi, site_no
    if (nh /= 0)then
      read(buff, *) pos_h(j)%p(:,1)
      write(*,*) pos_h(j)%p(:,1)
      read(id, *) pos_h(j)%p(:,2) 
      write(*,*) pos_h(j)%p(:,2)
      read(id, *) pos_h(j)%p(:,3)
      write(*,*) pos_h(j)%p(:,3)
      read(id, *) pos_h(j)%p(:,4)
      ns(j)=size(pos_h(j)%p(:,1))
    else
      ns(j)=0
    end if
    end do 
    close(id)
    nh_all = 0
    nh_max = 0
    do j = 1, nz
      nh_all = nh_all + ns(j)
      nh_max = max(nh_max,ns(j))
    end do
    do j = 1, nz
      nx = pshape(1,j)
      ny = pshape(2,j)
      do i = 1, ns(j)
       pos_h(j)%p(i,5) = xyz_to_site(pos_h(j)%p(i,1),pos_h(j)%p(i,2),j,pshape)
      end do
    end do
    if(present(hole_xi)) then
      if(allocated(hole_xi)) deallocate(hole_xi)
      allocate(hole_xi(nh_all, 2))
      k = 1
      do j = 1, nz
        do i = 1, ns(j)
          hole_xi(k,1) = pos_h(j)%p(i,5)
          hole_xi(k,2) = pos_h(j)%p(i,3)
          k = k + 1
        end do
      end do
    end if
    if(present(hole_chi)) then
      if(allocated(hole_chi)) deallocate(hole_chi)
      allocate(hole_chi(nh_all, 2))
      k = 1
      do j = 1, nz
        do i = 1, ns(j)
          hole_chi(k,1) = pos_h(j)%p(i,5)
          hole_chi(k,2) = pos_h(j)%p(i,4)
          k = k + 1
        end do
      end do
    end if
    if(present(layers)) then
      if(allocated(layers)) deallocate(layers)
      allocate(layers(nz))
      do j = 1, nz
        layers(j)%nx = pshape(1,j)
        layers(j)%ny = pshape(2,j)
        layers(j)%n = pshape(1,j)*pshape(2,j)
        layers(j)%nh = ns(j)
        layers(j)%ne = layers(j)%n - layers(j)%nh
        layers(j)%site_start = xyz_to_site(1,1,j,pshape)
        layers(j)%site_end = xyz_to_site(layers(j)%nx,layers(j)%ny,j,pshape)
        allocate(layers(j)%hole(layers(j)%nh))
        layers(j)%hole(:)%pos = pos_h(j)%p(:,5)
        layers(j)%hole(:)%xi = pos_h(j)%p(:,3)
    end do
    end if
    do j = 1, nz
      deallocate(pos_h(j)%p)
    end do
    deallocate(pos_h)
  end subroutine load_conf_3D

  subroutine load_conf_3D_2(filename, id, pshape, hole_xi, hole_chi, ne, layers)
    use parameter_mod
    character(*), intent(in)                           :: filename
    integer, intent(in)                                :: id
    integer, dimension(:,:), intent(out) ,allocatable, optional  :: pshape
    integer, allocatable, dimension(:, :), intent(out), optional :: hole_xi, hole_chi
    real(8), allocatable, dimension(:), intent(out), optional  :: ne
    type(cuo2_plane), allocatable, dimension(:), intent(out), optional :: layers
    integer, allocatable, dimension(:)  :: ns
    real(8), allocatable, dimension(:)  :: tne
    type(hole_param), allocatable, dimension(:)         :: pos_h
    integer                                            :: nh, nx, ny, nz, i, j, k, ns_all, ns_max
    character(10000)  buff

    open(id, file = filename)
    read(id, '()')      !skip a line "nz"
    read(id, *) nz
    if(allocated(pshape)) deallocate(pshape)
    allocate(pshape(2,nz),source=0)
    allocate(pos_h(nz))
     allocate(tne(nz),source=0.d0)
   allocate(ns(nz),source=0)
    do j = 1, nz
      read(id, '()') 
      read(id, *) nx, ny
      read(id, '()') 
      read(id, *) tne(j), ns(j)
      read(id, '()') 
      pshape(:,j) = [nx,ny]
    if (ns(j) /= 0)then
      allocate(pos_h(j)%p(ns(j), 5),source=0) 
      read(id, *) pos_h(j)%p(:,1)
      read(id, *) pos_h(j)%p(:,2)
      read(id, *) pos_h(j)%p(:,3)
      read(id, *) pos_h(j)%p(:,4)
      do i = 1, ns(j)
        pos_h(j)%p(i,5) = xyz_to_site(pos_h(j)%p(i,1),pos_h(j)%p(i,2),j,pshape)
      end do
      else 
         read(id, '()')
    end if
   end do
    ns_all = 0
    ns_max = 0
    do j = 1, nz
      ns_all = ns_all + ns(j)
      ns_max = max(ns_max,ns(j))
    end do
    if(allocated(hole_xi)) deallocate(hole_xi)
    allocate(hole_xi(ns_all, 2))
    if(allocated(hole_chi)) deallocate(hole_chi)
    allocate(hole_chi(ns_all, 2))
   close(id)
      k = 1
      do j = 1, nz
        do i = 1, ns(j)
          hole_chi(k,1) = pos_h(j)%p(i,5)
          hole_chi(k,2) = pos_h(j)%p(i,4)
          hole_xi(k,1) = pos_h(j)%p(i,5)
          hole_xi(k,2) = pos_h(j)%p(i,3)
         k = k + 1
        end do
      end do
      if(present(layers)) then
        allocate(layers(nz))
        do j = 1, nz
          layers(j)%nx = pshape(1,j)
          layers(j)%ny = pshape(2,j)
          layers(j)%n = pshape(1,j)*pshape(2,j)
          layers(j)%nh = ns(j)
          layers(j)%ne = tne(j)
          layers(j)%site_start = xyz_to_site(1,1,j,pshape)
          layers(j)%site_end = xyz_to_site(layers(j)%nx,layers(j)%ny,j,pshape)
          if (ns(j) /= 0)then
          allocate(layers(j)%hole(layers(j)%nh))
          layers(j)%hole(:)%pos = pos_h(j)%p(:,5) 
          layers(j)%hole(:)%xi = pos_h(j)%p(:,3)
          end if
        end do
      end if
      if(present(ne)) then
        call move_alloc(tne, ne)
      end if
    deallocate(pos_h)
  end subroutine load_conf_3D_2

  subroutine load_jexconf_3D(filename, id, pshape, jex_initial,jex_final,N)
    use parameter_mod
    character(*), intent(in)                           :: filename
    integer, intent(in)                                :: id
    integer, intent(in), dimension(:,:),allocatable    :: pshape
    real(8), dimension(:), intent(inout) ,allocatable :: jex_initial
    real(8), dimension(:), intent(inout) ,allocatable,optional :: jex_final
    integer , intent(out),optional                             :: N
    real(8), dimension(:),allocatable :: tmp_jex_i,tmp_jex_f
    integer                                            :: i, n_jex, site
    integer, allocatable, dimension(:,:)              :: site_r
    integer, dimension(3)                             :: tmp_r
    character(10000)  buff

    open(id, file = filename)
    read(id, '()')      !skip a line
    
    if(allocated(jex_initial)) deallocate(jex_initial)
    allocate(jex_initial(site_max(pshape)),source = 0d0)
    if(present(jex_final)) then
      if(allocated(jex_final)) deallocate(jex_final)
      allocate(jex_final(site_max(pshape)),source = 0d0)
    end if
    read(id,'(a)') buff
    n_jex = 0
    if(len_trim(buff) >= 1) then
      if(scan(buff(1:1),"1234567890-") /= 0) then
        n_jex = n_jex + 1
      end if
      do i = 2 ,len_trim(buff)
        if(scan(buff(i-1:i-1)," ") /= 0 .and. scan(buff(i:i),"1234567890-") /= 0) then
          n_jex = n_jex + 1
        end if
      end do
    end if
    allocate(site_r(n_jex, 3))
    allocate(tmp_jex_i(n_jex))
    allocate(tmp_jex_f(n_jex))
    if(n_jex == 0) return
    read(buff, *) site_r(:,1)
    read(id, *) site_r(:,2)
    read(id, *) site_r(:,3)
    read(id, *) tmp_jex_i
    if(present(jex_final)) then
      read(id, *) tmp_jex_f
    else
      read(id, "()")
    end if
    if(present(N)) then
      read(id, *) N
    end if
    do i = 1, n_jex
      tmp_r = site_r(i,1:3)
      site = r_to_site(tmp_r,pshape)
      if(site < 1) then
        print *, "error : ", tmp_r , "is not site"
        stop
      end if
      jex_initial(site) = tmp_jex_i(i)
      if(present(jex_final)) then
        jex_final(site) = tmp_jex_f(i)
      end if
    end do 
    close(id)
  end subroutine load_jexconf_3D
  
  subroutine load_min_E_state(filename, id, E, hole_chi, w_s)
    character(*), intent(in)                           :: filename
    integer, intent(in)                                :: id
    real(8), intent(out)                                :: E
    integer , intent(inout)                            :: hole_chi(:,:), w_s(:)
    character(10000)  buff
    
    open(id, file = filename)
    read(id, '()')      !skip a line
    read(id, '(a)' ) buff   ! "E_min =      ###.########...."
    !print *, buff 
    read(buff(8:),*) E
    read(id, '()')      !skip a line
    read(id, * ) hole_chi(:,2)
    read(id, '()')      !skip a line
    read(id, * ) w_s
    close(id)
  end subroutine load_min_E_state
  
  subroutine load_jbconf(filename, id, pshape, jex_from, jex_to, N, site_obs, Bz)
    use parameter_mod
    character(*), intent(in)                           :: filename
    integer, intent(in)                                :: id
    integer, intent(in), dimension(:, :), allocatable  :: pshape
    real(8), dimension(:), intent(inout), allocatable  :: jex_from, jex_to
    integer, intent(out)                               :: N, site_obs
    real(8), intent(out)                               :: Bz
    real(8), dimension(:), allocatable                 :: tmp_jex_a, tmp_jex_b
    integer                                            :: i, n_jex, site
    integer, allocatable, dimension(:,:)               :: site_r
    integer, dimension(3)                              :: tmp_r
    character(10000)  buff

    open(id, file = filename)
    read(id, '()')      !skip a line
    
    if(allocated(jex_from)) deallocate(jex_from)
    allocate(jex_from(site_max(pshape)),source = 0d0)
    if(allocated(jex_to)) deallocate(jex_to)
    allocate(jex_to(site_max(pshape)),source = 0d0)
    read(id,'(a)') buff
    n_jex = 0
    if(len_trim(buff) >= 1) then
      if(scan(buff(1:1),"1234567890-") /= 0) then
        n_jex = n_jex + 1
      end if
      do i = 2 ,len_trim(buff)
        if(scan(buff(i-1:i-1)," ") /= 0 .and. scan(buff(i:i),"1234567890-") /= 0) then
          n_jex = n_jex + 1
        end if
      end do
    end if
    allocate(site_r(n_jex, 3))
    allocate(tmp_jex_a(n_jex))
    allocate(tmp_jex_b(n_jex))
    if(n_jex == 0) return
    read(buff, *) site_r(:,1)
    read(id, *) site_r(:,2)
    read(id, *) site_r(:,3)
    read(id, *) tmp_jex_a
    read(id, *) tmp_jex_b
    read(id, *) N
    !read(id, *) site_obs
    read(id, *) tmp_r
    site_obs = r_to_site(tmp_r, pshape)
    read(id, '()')      !skip a line
    read(id, *) Bz
    do i = 1, n_jex
      tmp_r = site_r(i,1:3)
      site = r_to_site(tmp_r,pshape)
      if(site < 1) then
        print *, "error : ", tmp_r , "is not site"
        stop
      end if
      jex_from(site) = tmp_jex_a(i)
      jex_to(site) = tmp_jex_b(i)
    end do 
    close(id)
  end subroutine load_jbconf
end module io_mod
