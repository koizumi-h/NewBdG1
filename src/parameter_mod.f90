
module parameter_mod
  !USE angular_variable_mod,  ONLY : angular_variable, new_angular_variable
  implicit none

  integer, parameter  :: xi_trans_unit = 22
  integer, parameter  :: fchi_xi_trans_unit = 32
  integer, parameter  :: chi_trans_unit = 23
  integer, parameter  :: wavefunction_trans_unit = 20
  integer, parameter  :: enegy_trans_unit = 21
  
  type hole_info
    INTEGER     :: pos, xi, chi
  end type hole_info

  type cuo2_plane
    integer                                   :: nx, ny, nh, n
    real(8)                                   :: ne
    type(hole_info),dimension(:),allocatable  :: hole
    !real(8), DIMENSION(:), ALLOCATABLE        :: xi,chi, Sx, Sy, Sz
    !complex(8), dimension(:), allocatable     :: wf
    integer                                   :: site_start, site_end
  end type cuo2_plane

  type real_value_plane
    real(8),dimension(:), allocatable :: value
  end type real_value_plane

  type complex_value_plane
    complex(8),dimension(:), allocatable :: value
  end type complex_value_plane

  !type real_value_site
  !  integer site
  !  real(8) value
  !end type real_value_site
contains
  pure function site_max(pshape) result(site)
    integer, dimension(:,:), intent(in)   :: pshape
    integer  :: site
    integer  :: nz
    nz = size(pshape(1,:))
    if(nz == 0) then
      site = 0
      return
    end if
    site = xyz_to_site(pshape(1,nz),pshape(2,nz),nz,pshape)
  end function site_max
  pure function xyz_to_site(x, y, z, pshape) result(site)
    integer, intent(in) :: x, y, z
    integer, dimension(:,:), intent(in)   :: pshape
    integer  :: site
    integer, dimension(3) :: n
    integer i
    if(x < 1 .or. y<1 .or. z<1 .or. z > size(pshape(1,:))) then
      site = -1
      return
    end if
    site = 0
    do i = 1, z-1
      n(1) = pshape(1,i)
      n(2) = pshape(2,i)
      site = site + n(1)*n(2)
    end do
    n(1) = pshape(1,z)
    n(2) = pshape(2,z)
    site = site + x + (y-1)*n(1)
    if( x > n(1) .or. y > n(2)) site = -1
  end function xyz_to_site
  pure function r_to_site(r, pshape) result(site)
    integer, intent(in) :: r(3)
    integer, dimension(:,:), intent(in)   :: pshape
    integer  :: site
    site = xyz_to_site(r(1),r(2),r(3),pshape)
  end function r_to_site
  pure function site_to_r(site, pshape) result(r)
    integer, intent(in)   :: site
    integer, dimension(:,:), intent(in)   :: pshape
    integer, dimension(3) :: r,n
    integer i , s
    s = site
    n(3) = size(pshape(1,:))
    do i = 1, n(3)
      n(1) = pshape(1,i)
      n(2) = pshape(2,i)
      if(s <= n(1)*n(2)) then
        r(3) = i
        exit
      end if
      s = s - n(1)*n(2)
    end do
    r(1) = mod(s - 1, n(1)) + 1
    r(2) = (s - 1) / n(1) + 1
  end function site_to_r
  pure function site_to_x(site, pshape) result(x)
    integer, intent(in)   :: site
    integer, dimension(:,:), intent(in)   :: pshape
    integer  x
    integer,dimension(3) :: r
    r = site_to_r(site,pshape)
    x = r(1)
  end function site_to_x
  pure function site_to_y(site, pshape) result(y)
    integer, intent(in)   :: site
    integer, dimension(:,:), intent(in)   :: pshape
    integer  y
    integer,dimension(3) :: r
    r = site_to_r(site,pshape)
    y = r(2)
  end function site_to_y
  pure function site_to_z(site, pshape) result(z)
    integer, intent(in)   :: site
    integer, dimension(:,:), intent(in)   :: pshape
    integer  z
    integer,dimension(3) :: r
    r = site_to_r(site,pshape)
    z = r(3)
  end function site_to_z
  pure function make_pshape(layers) result(pshape)
    integer, dimension(:,:) ,allocatable  :: pshape
    type(cuo2_plane), dimension(:),intent(in) :: layers
    integer i
    allocate(pshape(2,size(layers)))
    do i=1, size(layers)
      pshape(1,i) = layers(i)%nx
      pshape(2,i) = layers(i)%ny
    end do
  end function make_pshape
  pure function make_holearray(layers) result(hole)
    integer, dimension(:) ,allocatable  :: hole
    type(cuo2_plane), dimension(:),intent(in) :: layers
    integer i, nh, s, e
    nh = 0
    do i = 1, size(layers)
      nh = nh + size(layers(i)%hole)
    end do
    allocate(hole(nh))
    s = 1
    do i = 1, size(layers)
      e = s - 1 + size(layers(i)%hole)
      hole(s:e) = layers(i)%hole(:)%pos
      s = e + 1
    end do
  end function make_holearray

end module parameter_mod
