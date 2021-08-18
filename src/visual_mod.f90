module visual_mod
  use path_list_mod
  use int_list_mod
  use parameter_mod
  implicit none
contains
  pure integer function digit(n) result(r)
    integer, intent(in) :: n
    integer             :: i
    do i = 1, 100
      if (n/10**i == 0) then
        r = i
        exit
      end if
    end do
    if (n < 0) r = r + 1
  end function digit
  function string(n) result(r)
    integer, intent(in) :: n
    character(digit(n))   :: r
    integer             :: i
    if (n < 0) then
      r(1:1) = '-'
      do i = 2, len(r)
        write(r(i:i), '(i1)') mod(-n, 10**(len(r)-i+1))/10**(len(r)-i)
      end do
    else
      do i = 1, len(r)
        write(r(i:i), '(i1)') mod(n, 10**(len(r)-i+1))/10**(len(r)-i)
      end do
    end if
  end function string
  subroutine print_phase(filename, unit_num, pshape, phase)
    character(*), intent(in)                            :: filename
    integer, intent(in)                                 :: unit_num
    integer, dimension(2), intent(in)                   :: pshape
    real(8), dimension(pshape(1)*pshape(2)), intent(in) :: phase
    integer                                             :: i, j
    open(unit_num, file = filename)
    do i = pshape(2), 1, -1
      do j = pshape(1)*(i-1)+1, pshape(1)*(i-1)+pshape(1)-1
        write(unit_num, '(e25.15)', advance = 'no') phase(j)
      end do
      write(unit_num, '(e25.15)', advance = 'yes') phase(pshape(1)*(i-1)+pshape(1))
    end do
    close(unit_num)
  end subroutine print_phase
  subroutine draw_phase(filename, unit_num, pshape, phase, h)
    ! outputs <p>, xi or chi, in the form of gnuplot
    character(*), intent(in)                            :: filename
    integer, intent(in)                                 :: unit_num
    integer, dimension(2), intent(in)                   :: pshape
    real(8), dimension(pshape(1)*pshape(2)), intent(in) :: phase
    integer, dimension(:), intent(in)                   :: h
    real(8), parameter                                  :: r = 0.40d0
    real(8)                                             :: dix, diy
    integer                                             :: i, j
    open(unit_num, file = filename)
    i_loop: do i = 1, size(phase)
      do j = 1, size(h); if (i == h(j)) cycle i_loop; end do
      dix = r * cos(phase(i))
      diy = r * sin(phase(i))
      write(unit_num, '(4(f18.15, 1x))') &
      & dble(mod((i-1), pshape(1))+1)-dix, dble((i-1)/pshape(1)+1)-diy, &
      & 2.0d0*dix, 2.0d0*diy
    end do i_loop
    close(unit_num)
  end subroutine draw_phase
  subroutine draw_stream(filename, pshape, s)
    character(*), intent(in)                            :: filename
    integer, dimension(2), intent(in)                   :: pshape
    real(8), dimension(pshape(1)*pshape(2)), intent(in) :: s
    integer                                             :: ix, iy
    open(11, file = filename)
    do iy = 1, pshape(2)
      do ix = 1, pshape(1)
        write(11, *) ix, iy, s(ix+(iy-1)*pshape(1))
      end do
      write(11, *)
    end do
    close(11)
  end subroutine draw_stream
  subroutine write_format_of_xi(filename, unit_num, pshape, h)
    !! 'trim' is necessary for bug fixing
    character(*), intent(in)                       :: filename
    integer, intent(in)                            :: unit_num
    integer, dimension(2), intent(in)              :: pshape
    integer, dimension(:, :), intent(in), optional :: h
    character(8)                                   :: str
    integer                                        :: nx, ny, i
    nx = pshape(1); ny = pshape(2)
    open(unit_num, file = filename)
    write(unit_num, *) 'reset'
    write(unit_num, *) 'set xtics 1'
    write(unit_num, *) 'set ytics 1'
    write(unit_num, *) trim('set xrange[0:'//string(nx+1)//']')
    write(unit_num, *) trim('set yrange[0:'//string(ny+1)//']')
    write(unit_num, *) trim('set size ratio '//string(ny+1)//'.0' &
    &                 //'/'//string(nx+1)//'.0')
    write(unit_num, *) 'unset key'
    if (present(h)) then
      do i = 1, size(h, 1)
        select case(h(i, 2))
          case( 1); str = ' "M" at '
          case(-1); str = ' "A" at '
          case( 0); str = ' "X" at '
        end select
        write(unit_num, *) trim('set label '//string(i)//str &
        &     //string(mod((h(i, 1)-1), nx) + 1)//', ' &
        &     //string((h(i, 1)-1)/nx + 1)//' center')
      end do
    end if
    close(unit_num)
  end subroutine write_format_of_xi
  subroutine write_format_of_chi(filename, unit_num, pshape, h)
    !! 'trim' is necessary for bug fixing
    character(*), intent(in)                       :: filename
    integer, intent(in)                            :: unit_num
    integer, dimension(2), intent(in)              :: pshape
    integer, dimension(:, :), intent(in), optional :: h
    character(8)                                   :: str
    integer                                        :: nx, ny, i
    nx = pshape(1); ny = pshape(2)
    open(unit_num, file = filename)
    write(unit_num, *) 'reset'
    write(unit_num, *) 'set xtics 1'
    write(unit_num, *) 'set ytics 1'
    write(unit_num, *) trim('set xrange[0:'//string(nx+1)//']')
    write(unit_num, *) trim('set yrange[0:'//string(ny+1)//']')
    write(unit_num, *) trim('set size ratio '//string(ny+1)//'.0' &
    &                 //'/'//string(nx+1)//'.0')
    write(unit_num, *) 'unset key'
    if (present(h)) then
      do i = 1, size(h, 1)
        select case(h(i, 2))
          case( 1); str = ' "m" at '
          case(-1); str = ' "a" at '
          case( 0); str = ' "x" at '
        end select
        write(unit_num, *) trim('set label '//string(i)//str &
        &     //string(mod((h(i, 1)-1), nx) + 1)//', ' &
        &     //string((h(i, 1)-1)/nx + 1)//' center')
      end do
    end if
    close(unit_num)
  end subroutine write_format_of_chi

  subroutine write_format_of_xi_3D(filename, unit_num, pshape, h)
    !! 'trim' is necessary for bug fixing
    character(*), intent(in)                       :: filename
    integer, intent(in)                            :: unit_num
    integer, dimension(:,:), intent(in)              :: pshape
    integer, dimension(:, :), intent(in), optional :: h
    character(8)                                   :: str
    integer                                        :: nx_max, ny_max, nx, ny, nz, i, j, z
    integer, dimension(:), allocatable :: n_s
    nz = size(pshape(1,:))
    nx_max = maxval(pshape(1,:))
    ny_max = maxval(pshape(2,:))
    allocate(n_s(nz))
    n_s(1) = 0
    do i = 1, nz-1
      n_s(i+1) = n_s(i) + pshape(1,i)*pshape(2,i)
    end do
    open(unit_num, file = filename)
    write(unit_num, *) 'reset'
    write(unit_num, *) 'set xtics 1'
    write(unit_num, *) 'set ytics 1'
    write(unit_num, *) 'set ztics 1'
    write(unit_num, *) trim('set xrange[0:'//string(nx_max+1)//']')
    write(unit_num, *) trim('set yrange[0:'//string(ny_max+1)//']')
    if(nz /= 1) then
      !write(unit_num, *) trim('set zrange[0:'//string(nz+1)//']')
      write(unit_num, *) trim('set zrange[0.9:'//string(nz)//'.1]')
    end if
    write(unit_num, *) trim('set view equal xy')
    !write(unit_num, *) trim('set size ratio '//string(ny+1)//'.0' &
    !&                 //'/'//string(nx+1)//'.0')
    write(unit_num, *) 'unset key'
    if (present(h)) then
      do i = 1, size(h, 1)
        do j = 1, nz
          if(h(i,1) < n_s(j)) then
            exit
          end if
          z = j
        end do
        nx = pshape(1,z)
        ny = pshape(2,z)
        select case(h(i, 2))
          case( 1); str = ' "M" at '
          case(-1); str = ' "A" at '
          case( 0); str = ' "X" at '
        end select
        write(unit_num, *) trim('set label '//string(i)//str &
        &     //string(mod((h(i, 1)-n_s(z)-1), nx) + 1)//', ' &
        &     //string((h(i, 1)-n_s(z)-1)/nx + 1) //','&
        &     //string(z)//' center')
      end do
    end if
    close(unit_num)
  end subroutine write_format_of_xi_3D
  subroutine write_format_of_chi_3D(filename, unit_num, pshape, h)
    !! 'trim' is necessary for bug fixing
    character(*), intent(in)                       :: filename
    integer, intent(in)                            :: unit_num
    integer, dimension(:,:), intent(in)              :: pshape
    integer, dimension(:, :), intent(in), optional :: h
    character(8)                                   :: str
    integer                                        :: nx_max, ny_max, nx, ny, nz, i, j, z
    integer, dimension(:), allocatable :: n_s
    nz = size(pshape(1,:))
    nx_max = maxval(pshape(1,:))
    ny_max = maxval(pshape(2,:))
    allocate(n_s(nz))
    n_s(1) = 0
    do i = 1, nz - 1
      n_s(i+1) = n_s(i) + pshape(1,i)*pshape(2,i)
    end do
    open(unit_num, file = filename)
    write(unit_num, *) 'reset'
    write(unit_num, *) 'set xtics 1'
    write(unit_num, *) 'set ytics 1'
    write(unit_num, *) 'set ztics 1'
    write(unit_num, *) trim('set xrange[0:'//string(nx_max+1)//']')
    write(unit_num, *) trim('set yrange[0:'//string(ny_max+1)//']')
    if(nz /= 1) then
      write(unit_num, *) trim('set zrange[0.9:'//string(nz)//'.1]')
    end if
    write(unit_num, *) trim('set view equal xy')
    !write(unit_num, *) trim('set size ratio '//string(ny+1)//'.0' &
    !&                 //'/'//string(nx+1)//'.0')
    write(unit_num, *) 'unset key'
    if (present(h)) then
      do i = 1, size(h, 1)
        do j = 1, nz
          if(h(i,1) < n_s(j)) then
            exit
          end if
          z = j
        end do
        nx = pshape(1,z)
        ny = pshape(2,z)
        select case(h(i, 2))
          case( 1); str = ' "m" at '
          case(-1); str = ' "a" at '
          case( 0); str = ' "x" at '
        end select
        write(unit_num, *) trim('set label '//string(i)//str &
        &     //string(mod((h(i, 1)-n_s(z)-1), nx) + 1)//', ' &
        &     //string((h(i, 1)-n_s(z)-1)/nx + 1) //','&
        &     //string(z)//' center')
      end do
    end if
    close(unit_num)
  end subroutine write_format_of_chi_3D

  subroutine write_format_of_chi_sb(filename, unit_num, pshape, h, pole_s, wn_s)
    !! 'trim' is necessary for bug fixing
    character(*), intent(in)                       :: filename
    integer, intent(in)                            :: unit_num
    integer, dimension(:,:), intent(in)              :: pshape
    integer, dimension(:, :), intent(in), optional :: h
    real(8), dimension(:, :), intent(in), optional :: pole_s
    integer, dimension(:), intent(in), optional :: wn_s
    character(8)                                   :: str
    integer                                        :: nx_max, ny_max, nx, ny, nz, i, j, z
    integer, dimension(:), allocatable :: n_s
    nz = size(pshape(1,:))
    nx_max = maxval(pshape(1,:))
    ny_max = maxval(pshape(2,:))
    allocate(n_s(nz))
    n_s(1) = 0
    do i = 1, nz - 1
      n_s(i+1) = n_s(i) + pshape(1,i)*pshape(2,i)
    end do
    open(unit_num, file = filename)
    write(unit_num, *) 'reset'
    write(unit_num, *) 'set xtics 1'
    write(unit_num, *) 'set ytics 1'
    write(unit_num, *) 'set ztics 1'
    write(unit_num, *) trim('set xrange[0:'//string(nx_max+1)//']')
    write(unit_num, *) trim('set yrange[0:'//string(ny_max+1)//']')
    if(nz /= 1) then
      write(unit_num, *) trim('set zrange[0.9:'//string(nz)//'.1]')
    end if
    write(unit_num, *) trim('set view equal xy')
    !write(unit_num, *) trim('set size ratio '//string(ny+1)//'.0' &
    !&                 //'/'//string(nx+1)//'.0')
    write(unit_num, *) 'unset key'
    if (present(h)) then
      do i = 1, size(h, 1)
        do j = 1, nz
          if(h(i,1) < n_s(j)) then
            exit
          end if
          z = j
        end do
        nx = pshape(1,z)
        ny = pshape(2,z)
        select case(h(i, 2))
          case( 1); str = ' "m" at '
          case(-1); str = ' "a" at '
          case( 0); str = ' "x" at '
        end select
        write(unit_num, *) trim('set label '//string(i)//str &
        &     //string(mod((h(i, 1)-n_s(z)-1), nx) + 1)//', ' &
        &     //string((h(i, 1)-n_s(z)-1)/nx + 1) //','&
        &     //string(z)//' center')
      end do
    end if
    if (present(pole_s) .and. present(wn_s)) then
      do i = 1, size(pole_s, 2)
        select case(wn_s(i))
          case( 1); str = ' "m" at '
          case(-1); str = ' "a" at '
          case( 0); str = ' "x" at '
        end select
        !write(unit_num, *) trim('set label '//string(i)//str &
        !&     //string( pole_s(1, i) )//', ' &
        !&     //string( pole_s(2, i) )//',' &
        !&     //string(1)//' center')
        write(unit_num, "(a,f6.2,a,f6.2,a)") 'set label '//string(size(h,1)+i)//str &
        &     , pole_s(1, i),  ', ' &
        &     , pole_s(2, i), ', 1 center'
      end do
    end if
    close(unit_num)
  end subroutine write_format_of_chi_sb

  subroutine draw_phase_3D(filename, unit_num, pshape, phase, h)
    ! outputs <p>, xi or chi, in the form of gnuplot
    character(*), intent(in)                            :: filename
    integer, intent(in)                                 :: unit_num
    integer, dimension(:,:), intent(in)                   :: pshape
    real(8), dimension(:), intent(in)                   :: phase
    integer, dimension(:), intent(in)                   :: h
    real(8), parameter                                  :: r = 0.40d0
    real(8)                                             :: dix, diy
    integer                                             :: i, j, k,l, nx, ny, nz, n
    nz = size(pshape(1,:))
    open(unit_num, file = filename)
    k = 0
    do i = 1, nz
      nx = pshape(1,i)
      ny = pshape(2,i)
      n = nx * ny
      j_loop : do j = 1, n
        k = k + 1
        do l = 1, size(h); if (k == h(l)) cycle j_loop; end do ! skip sites of hole
        dix = r * cos(phase(k))
        diy = r * sin(phase(k))
        write(unit_num, '(6(f18.15, 1x))') &
          & dble(mod((j-1), nx)+1)-dix, dble((j-1)/nx+1)-diy, &
          & dble(i), &
          & 2.0d0*dix, 2.0d0*diy, 0d0
      end do j_loop
    end do
    close(unit_num)
  end subroutine draw_phase_3D
  
  ! filename + (z) + extension
  subroutine draw_phase_3D_layerbylayer(filename, extension,unit_num, pshape, phase, h)
    character(*), intent(in)                            :: filename, extension
    integer, intent(in)                                 :: unit_num
    integer, dimension(:,:), intent(in)                   :: pshape
    real(8), dimension(:), intent(in)                   :: phase
    integer, dimension(:), intent(in)                   :: h
    real(8), parameter                                  :: r = 0.40d0
    real(8)                                             :: dix, diy
    integer                                             :: i, j, k,l, nx, ny, nz, n
    nz = size(pshape(1,:))
    k = 0
    do i = 1, nz
      open(unit_num, file = filename//string(i)//extension)
      nx = pshape(1,i)
      ny = pshape(2,i)
      n = nx * ny
      j_loop : do j = 1, n
        k = k + 1
        do l = 1, size(h); if (k == h(l)) cycle j_loop; end do ! skip sites of hole
        dix = r * cos(phase(k))
        diy = r * sin(phase(k))
        write(unit_num, '(6(f20.15, 1x))') &
          & dble(mod((j-1), nx)+1)-dix, dble((j-1)/nx+1)-diy, &
          & dble(i), &
          & 2.0d0*dix, 2.0d0*diy, 0d0
      end do j_loop
      close(unit_num)
    end do
  end subroutine draw_phase_3D_layerbylayer
  !subroutine draw_phase_3D_layers(filename, unit_num, pshape, layers, h)
  !  character(*), intent(in)                            :: filename
  !  integer, intent(in)                                 :: unit_num
  !  integer, dimension(:,:), intent(in)                   :: pshape
  !  type(layers), dimension(:), intent(in)                :: layers
  !  integer, dimension(:), intent(in)                   :: h
  !  real(8), dimension(:), intent(in)                   :: phase
  !  integer i
  !  !do i=1, size(layers)

  !  !end do


  !end subroutine
  subroutine draw_path_3D(filename, unit_num, pshape, pl)
    character(*), intent(in)                            :: filename
    integer, intent(in)                                 :: unit_num
    integer, dimension(:,:), intent(in)                   :: pshape
    type(path_list), intent(in)                         :: pl
    integer                                             :: i
    integer ri(3),rf(3)
    open(unit_num, file = filename)

    do i = 1, pl%length()
      ri = site_to_r(pl%value(i)%i,pshape)
      rf = site_to_r(pl%value(i)%f,pshape)
      write(unit_num,*) ri, rf - ri
    end do

    close(unit_num)
  end subroutine draw_path_3D
  subroutine draw_circle(filename, unit_num, pshape, cl, cpl, wl)
    character(*), intent(in)                            :: filename
    integer, intent(in)                                 :: unit_num
    !integer, dimension(:,:), intent(in)                   :: pshape
    integer, dimension(:,:), intent(in)                   :: pshape
    type(int_list),dimension(:), intent(in)                         :: cl
    type(path_list), intent(in)                         :: cpl
    integer, dimension(:), intent(in), optional   :: wl
    integer                                             :: i, j
    integer si, sf
    real(8) g(3), p(3,5)
    open(unit_num, file = filename)

    do i=1, size(cl)
      g(:) = 0d0
      do j=1, cl(i)%length()
        if(cl(i)%value(j) > 0) then
          si = cpl%value(abs(cl(i)%value(j)))%i
          sf = cpl%value(abs(cl(i)%value(j)))%f
        else
          si = cpl%value(abs(cl(i)%value(j)))%f
          sf = cpl%value(abs(cl(i)%value(j)))%i
        end if
        !p(1,j) = mod((si - 1), pshape(1)) + 1
        !p(2,j) = (si - 1)/ pshape(1) + 1
        p(:,j) = site_to_r(si, pshape)
        g(:) = g(:) + p(:,j)
      end do
      g(:) = g(:) / cl(i)%length()
      do j=1, cl(i)%length()
        p(:,j) = p(:,j) + (g(:) - p(:,j))*0.2d0
      end do
      p(:,cl(i)%length()+1) = p(:,1)
      do j=1, cl(i)%length()
        if(present(wl)) then
        write(unit_num, '(6(f20.15, 1x),i4)') p(:,j), p(:,j+1) - p(:,j), wl(i)
        !write(*,*) p(:,j), p(:,j+1) - p(:,j)
        else
          write(unit_num, '(6(f20.15, 1x))') p(:,j), p(:,j+1) - p(:,j)
        end if
      end do
    end do
    close(unit_num)
  end subroutine draw_circle
  subroutine draw_path_delta(filename, unit_num, pshape, pl, delta)
    character(*), intent(in)                            :: filename
    integer, intent(in)                                 :: unit_num
    integer, dimension(:,:), intent(in)                   :: pshape
    type(path_list), intent(in)                         :: pl
    complex(8), dimension(pl%length()), intent(in)      :: delta
    integer                                             :: i
    integer ri(3),rf(3)
    open(unit_num, file = filename)

    do i = 1, pl%length()
      ri = site_to_r(pl%value(i)%i,pshape)
      rf = site_to_r(pl%value(i)%f,pshape)
      write(unit_num,"(8f20.15)") real(ri), real(rf - ri), real(delta(i)), aimag(delta(i))
    end do

    close(unit_num)
  end subroutine draw_path_delta
  subroutine draw_delta(filename, unit_num, pshape, pl, delta)
    character(*), intent(in)                            :: filename
    integer, intent(in)                                 :: unit_num
    integer, dimension(:,:), intent(in)                   :: pshape
    type(path_list), intent(in)                         :: pl
    complex(8), dimension(pl%length()), intent(in)      :: delta
    integer                                             :: i
    integer ri(3),rf(3)
    open(unit_num, file = filename)

    do i = 1, pl%length()
      ri = site_to_r(pl%value(i)%i,pshape)
      rf = site_to_r(pl%value(i)%f,pshape)
      write(unit_num,"(6f20.15)") real(rf + ri)/2d0, abs(delta(i)), real(delta(i)), aimag(delta(i))
    end do

    close(unit_num)
  end subroutine draw_delta

  subroutine draw_current(filename, unit_num, pshape, pl, current)
    character(*), intent(in)                            :: filename
    integer, intent(in)                                 :: unit_num
    integer, dimension(:,:), intent(in)                 :: pshape
    type(path_list), intent(in)                         :: pl
    real(8), dimension(:,:), intent(in)                 :: current
    type(path)                      :: p
    real(8)                         :: ratio
    real(8),dimension(3)  ::  r, v, t
    integer i,j,n
    n = site_max(pshape)
    !if (present(magnitude))then
    !  ratio = magnitude
    !else
      ratio = 0d0
      !write(*,*) this%n
      do i = 1, n
        do j = 1, n
          !write(*,*) i,j,this%bond(i,j) 
          if( ratio < abs(current(i,j))) then
            ratio = abs(current(i,j))
          end if
        end do 
      end do
      ratio = 1d0 / ratio
      !write(*,*) "ratio", ratio ," = 1/ ",1d0/ratio 
    !end if
    open(unit_num, file = filename)
  
    
    do i = 1, pl%length()
      p = pl%index(i)
      !if(p%i == 1 .and. p%f == 7) then
      t = (site_to_r(p%f, pshape) - site_to_r(p%i, pshape))
      v = ratio * current(p%f, p%i) * t/(sqrt(t(1)**2 + t(2)**2 + t(3)**2))
      !v_x = 0.5d0
      !v_y = 0.5d0
      r = (site_to_r(p%f, pshape) + site_to_r(p%i, pshape))/2d0 - v / 2d0
      write(unit_num, '(6(f19.15, 1x))') r, v
      !print *, sqrt(v(1)**2 + v(2)**2 + v(3)**2)
      !end if
    end do
    close(unit_num)
  end subroutine draw_current
  
  subroutine draw_current_ratio(filename, unit_num, pshape, pl, current, ratio)
    character(*), intent(in)                            :: filename
    integer, intent(in)                                 :: unit_num
    integer, dimension(:,:), intent(in)                 :: pshape
    type(path_list), intent(in)                         :: pl
    real(8), dimension(:,:), intent(in)                 :: current
    type(path)                      :: p
    real(8), intent(in)                         :: ratio
    real(8),dimension(3)  ::  r, v, t
    integer i,j,n
    n = site_max(pshape)
    open(unit_num, file = filename)
    do i = 1, pl%length()
      p = pl%index(i)
      t = (site_to_r(p%f, pshape) - site_to_r(p%i, pshape))
      v = ratio * current(p%f, p%i) * t/(sqrt(t(1)**2 + t(2)**2 + t(3)**2))
      r = (site_to_r(p%f, pshape) + site_to_r(p%i, pshape))/2d0 - v / 2d0
      write(unit_num, '(6(f19.15, 1x))') r, v
    end do
    close(unit_num)
  end subroutine draw_current_ratio
  
  subroutine draw_spin(filename, unit_num, pshape, Sx, Sy, Sz)
    character(*), intent(in)                            :: filename
    integer, intent(in)                                 :: unit_num
    integer, dimension(:,:), intent(in)                 :: pshape
    real(8), dimension(:), intent(in)                   :: Sx, Sy, Sz
    integer                                             :: i, n
    real(8) :: rvec(3),vec(3)
    n = site_max(pshape)
    open(unit_num, file = filename)
    do i = 1, n
      vec(1) = Sx(i)
      vec(2) = Sy(i)
      vec(3) = Sz(i)
      rvec= site_to_r(i, pshape) - vec * 0.5d0
      write(unit_num, "(6f20.15)") rvec, vec
    end do
    close(unit_num)
  end subroutine draw_spin
  
  subroutine draw_spin_normalize(filename, unit_num, pshape, Sx, Sy, Sz)
    character(*), intent(in)                            :: filename
    integer, intent(in)                                 :: unit_num
    integer, dimension(:,:), intent(in)                 :: pshape
    real(8), dimension(:), intent(in)                   :: Sx, Sy, Sz
    integer                                             :: i, n
    real(8) :: rvec(3),vec(3)
    n = site_max(pshape)
    open(unit_num, file = filename)
    do i = 1, n
      vec(1) = Sx(i)
      vec(2) = Sy(i)
      vec(3) = Sz(i)
      !print *, i, sqrt(vec(1)**2d0+vec(2)**2d0+vec(3)**2d0)
      if( vec(1) /= 0d0 .and. vec(2) /= 0d0 .and. vec(3)/=0d0) then
        vec = 0.8d0 * vec / sqrt(vec(1)**2d0+vec(2)**2d0+vec(3)**2d0)
      end if
      rvec= site_to_r(i, pshape) - vec * 0.5d0
      write(unit_num, "(6f20.15)") rvec, vec
    end do
    close(unit_num)
  end subroutine draw_spin_normalize
  
  subroutine draw_winding_circle(filename_plus, filename_minus, unit_num1, unit_num2, pshape, cl, cpl, wl)
    character(*), intent(in)                            :: filename_plus, filename_minus
    integer, intent(in)                                 :: unit_num1, unit_num2
    !integer, dimension(:,:), intent(in)                   :: pshape
    integer, dimension(:,:), intent(in)                   :: pshape
    type(int_list),dimension(:), intent(in)                         :: cl
    type(path_list), intent(in)                         :: cpl
    integer, dimension(:), intent(in)                   :: wl
    integer                                             :: i, j
    integer si, sf
    real(8) g(3), p(3,5)
    open(unit_num1, file = filename_plus)
    open(unit_num2, file = filename_minus)

    do i=1, size(cl)
      if(wl(i) == 0) cycle
      g(:) = 0d0
      do j=1, cl(i)%length()
        if(cl(i)%value(j) > 0) then
          si = cpl%value(abs(cl(i)%value(j)))%i
          sf = cpl%value(abs(cl(i)%value(j)))%f
        else
          si = cpl%value(abs(cl(i)%value(j)))%f
          sf = cpl%value(abs(cl(i)%value(j)))%i
        end if
        !p(1,j) = mod((si - 1), pshape(1)) + 1
        !p(2,j) = (si - 1)/ pshape(1) + 1
        p(:,j) = site_to_r(si, pshape)
        g(:) = g(:) + p(:,j)
      end do
      g(:) = g(:) / cl(i)%length()
      do j=1, cl(i)%length()
        p(:,j) = p(:,j) + (g(:) - p(:,j))*0.2d0
      end do
      p(:,cl(i)%length()+1) = p(:,1)
      do j=1, cl(i)%length()
        if(wl(i)==1) then
          write(unit_num1, '(6(f20.15, 1x))') p(:,j), p(:,j+1) - p(:,j)
        else if(wl(i)==-1) then
          write(unit_num2, '(6(f20.15, 1x))') p(:,j), p(:,j+1) - p(:,j)
        else
          write(*, *) "warning : winding number other than 1 ,0 ,-1"
        end if
      end do
    end do
    close(unit_num1)
    close(unit_num2)
  end subroutine draw_winding_circle
end module visual_mod
