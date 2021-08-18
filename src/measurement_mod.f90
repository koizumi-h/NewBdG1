module measurement_mod
use math_mod, only : dfermi
!use io_mod, only : save_txt
use parameter_mod, only : site_to_r
implicit none

contains
  subroutine save_stm_rho_vs_e_step(file, id, n, wf, eg, kbt, e_range, n_split)
    character(*), intent(in) :: file
    integer, intent(in) :: id, n
    complex(8), dimension(4*n,4*n), intent(in) :: wf
    real(8), dimension(4*n), intent(in) :: eg
    real(8), intent(in) :: kbt
    real(8), dimension(2), intent(in) :: e_range   ! range is [from, to]
    integer, intent(in) :: n_split
    real(8), dimension(n_split) :: e_value
    integer i
    do i=0, n_split-1
      e_value(i+1) = e_range(1) + (e_range(2) - e_range(1))/(n_split-1d0) * float(i)
    end do
    call save_stm_rho_vs_e(file, id, n, wf, eg, kbt, e_value)
  end subroutine save_stm_rho_vs_e_step
  
  subroutine save_stm_rho_vs_e(file, id, n, wf, eg, kbt, e_value)
    character(*), intent(in) :: file
    integer, intent(in) :: id, n
    complex(8), dimension(4*n,4*n), intent(in) :: wf
    real(8), dimension(4*n), intent(in) :: eg
    real(8), intent(in) :: kbt
    real(8), dimension(:), intent(in) :: e_value
    real(8), dimension(size(e_value)) :: rho
    integer :: i
    call stm_rho_vs_e(n, wf, eg, kbt, e_value, rho)
    open(id, file = file)
    do i=1, size(rho)
      write(id, "(2f30.20)") e_value(i), rho(i)
    end do
    close(id)
  end subroutine save_stm_rho_vs_e

  subroutine save_stm_point_rho_vs_e_step(file, id, point, radius, pshape, n, wf, eg, kbt, e_range, n_split)
    character(*), intent(in) :: file
    integer, intent(in) :: id, n
    real(8), intent(in) :: point(2), radius
    integer, dimension(:,:), intent(in) :: pshape
    complex(8), dimension(4*n,4*n), intent(in) :: wf
    real(8), dimension(4*n), intent(in) :: eg
    real(8), intent(in) :: kbt
    real(8), dimension(2), intent(in) :: e_range   ! range is [from, to]
    integer, intent(in) :: n_split
    real(8), dimension(n_split) :: e_value
    real(8), dimension(size(e_value)) :: rho
    integer i
    do i=0, n_split-1
      e_value(i+1) = e_range(1) + (e_range(2) - e_range(1))/(n_split-1d0) * float(i)
    end do
    call stm_point_rho_vs_e(point, radius, pshape, n, wf, eg, kbt, e_value, rho)
    open(id, file = file)
    do i=1, size(rho)
      write(id, "(2f30.20)") e_value(i), rho(i)
    end do
    close(id)
  end subroutine save_stm_point_rho_vs_e_step
  
  subroutine stm_rho_vs_e_step(n, wf, eg, kbt, e_range, rho, n_split)
    integer, intent(in) :: n
    complex(8), dimension(4*n,4*n), intent(in) :: wf
    real(8), dimension(4*n), intent(in) :: eg
    real(8), intent(in) :: kbt
    real(8), dimension(2), intent(in) :: e_range   ! range is [from, to]
    real(8), dimension(:), allocatable, intent(out) :: rho
    integer, intent(in) :: n_split
    real(8), dimension(n_split) :: e_value
    integer i
    allocate(rho(n_split+1))
    do i=0, n_split-1
      e_value(i+1) = e_range(1) + (e_range(2) - e_range(1))/(n_split-1d0) * float(i)
    end do
    call stm_rho_vs_e(n, wf, eg, kbt, e_value, rho)
  end subroutine stm_rho_vs_e_step
  
  !subroutine stm_rho_vs_e(n, wf, eg, kbt, e_value, rho)
  !  integer, intent(in) :: n
  !  complex(8), dimension(4*n,4*n), intent(in) :: wf
  !  real(8), dimension(4*n), intent(in) :: eg
  !  real(8), intent(in) :: kbt
  !  real(8), dimension(:), intent(in) :: e_value
  !  real(8), dimension(size(e_value)), intent(out) :: rho
  !  integer :: i, j, r

  !  rho(:) = 0d0
  !  do i = 1, size(e_value)
  !    !do r = 1, size(eg)
  !    !  if(eg(r) <= 0d0) cycle
  !    do r = size(eg)/2 + 1, size(eg)
  !    !do r = 1, size(eg)/2
  !      do j = 1, n
  !        rho(i) = rho(i) - conjg(wf(4*j-3, r))*wf(4*j-3, r)*dfermi(eg(r) - e_value(i), kbt)  & ! u_{j up}
  !                &       - conjg(wf(4*j-2, r))*wf(4*j-2, r)*dfermi(eg(r) - e_value(i), kbt)  & ! u_{j down}
  !                &       - conjg(wf(4*j-1, r))*wf(4*j-1, r)*dfermi(eg(r) + e_value(i), kbt)  & ! v_{j up}
  !                &       - conjg(wf(4*j  , r))*wf(4*j  , r)*dfermi(eg(r) + e_value(i), kbt)    ! v_{j down}
  !      end do
  !    end do
  !  end do

  !end subroutine stm_rho_vs_e
  subroutine stm_rho_vs_e(n, wf, eg, kbt, e_value, rho)
    integer, intent(in) :: n
    complex(8), dimension(4*n,4*n), intent(in) :: wf
    real(8), dimension(4*n), intent(in) :: eg
    real(8), intent(in) :: kbt
    real(8), dimension(:), intent(in) :: e_value
    real(8), dimension(size(e_value)), intent(out) :: rho
    integer, dimension(n) :: sites
    integer :: i
    do i = 1, n
      sites(i) = i
    end do
    call stm_sites_rho_vs_e(sites, n, wf, eg, kbt, e_value, rho)
  end subroutine stm_rho_vs_e

  subroutine stm_point_rho_vs_e(point, radius, pshape, n, wf, eg, kbt, e_value, rho)
    real(8), intent(in) :: point(2), radius
    integer, dimension(:,:), intent(in) :: pshape
    integer, intent(in) :: n
    complex(8), dimension(4*n,4*n), intent(in) :: wf
    real(8), dimension(4*n), intent(in) :: eg
    real(8), intent(in) :: kbt
    real(8), dimension(:), intent(in) :: e_value
    real(8), dimension(size(e_value)), intent(out) :: rho
    integer, dimension(:), allocatable :: sites
    integer :: i, pos(3), n_site
    n_site = 0
    do i=1, pshape(1,1)*pshape(2,1)
      pos = site_to_r(i, pshape)
      if((pos(1)- point(1))**2 + (pos(2)-point(2))**2 <= radius**2 ) then
        n_site = n_site + 1
      end if
    end do
    allocate(sites(n_site))
    print *, n_site
    n_site = 1
    do i=1, pshape(1,1)*pshape(2,1)
      pos = site_to_r(i, pshape)
      if((pos(1)- point(1))**2 + (pos(2)-point(2))**2 <= radius**2 ) then
        sites(n_site) = i
        n_site = n_site + 1
      end if
    end do
    print *, sites
    call stm_sites_rho_vs_e(sites, n, wf, eg, kbt, e_value, rho)
    deallocate(sites)
  end subroutine stm_point_rho_vs_e
  
  subroutine stm_sites_rho_vs_e(sites, n, wf, eg, kbt, e_value, rho)
    integer, dimension(:), intent(in) :: sites
    integer, intent(in) :: n
    complex(8), dimension(4*n,4*n), intent(in) :: wf
    real(8), dimension(4*n), intent(in) :: eg
    real(8), intent(in) :: kbt
    real(8), dimension(:), intent(in) :: e_value
    real(8), dimension(size(e_value)), intent(out) :: rho
    integer :: i, j, k, r
    rho(:) = 0d0
    do i = 1, size(e_value)
      !do r = 1, size(eg)
      !  if(eg(r) <= 0d0) cycle
      do r = size(eg)/2 + 1, size(eg)
        do k = 1, size(sites)
          j = sites(k)
          rho(i) = rho(i) - conjg(wf(4*j-3, r))*wf(4*j-3, r)*dfermi(eg(r) - e_value(i), kbt)  & ! u_{j up}
                  &       - conjg(wf(4*j-2, r))*wf(4*j-2, r)*dfermi(eg(r) - e_value(i), kbt)  & ! u_{j down}
                  &       - conjg(wf(4*j-1, r))*wf(4*j-1, r)*dfermi(eg(r) + e_value(i), kbt)  & ! v_{j up}
                  &       - conjg(wf(4*j  , r))*wf(4*j  , r)*dfermi(eg(r) + e_value(i), kbt)    ! v_{j down}
        end do
      end do
    end do
  end subroutine stm_sites_rho_vs_e


  subroutine stm_alloc_e_value(e_start, e_end, e_step, e_value)
    real(8), intent(in) :: e_start, e_end
    integer, intent(in) :: e_step
    real(8), dimension(:), allocatable, intent(inout) :: e_value
    integer i
    if(allocated(e_value)) deallocate(e_value)
    allocate(e_value(e_step))
    do i=0, e_step-1
      e_value(i+1) = e_start + (e_end - e_start)/(e_step-1d0) * float(i)
    end do
  end subroutine stm_alloc_e_value

  subroutine stm_simple_rho_vs_e(site, wf, eg, kbt, e_value, rho)
    integer, intent(in) :: site
    complex(8), dimension(:,:), intent(in) :: wf
    real(8), dimension(:), intent(in) :: eg
    real(8), intent(in) :: kbt
    real(8), dimension(:), intent(in) :: e_value
    real(8), dimension(size(e_value)), intent(out) :: rho
    integer :: i, j, k, r
    rho(:) = 0d0
    do i = 1, size(e_value)
      do r = size(eg)/2 + 1, size(eg)
        j = site
        rho(i) = rho(i) - conjg(wf(4*j-3, r))*wf(4*j-3, r)*dfermi(eg(r) - e_value(i), kbt)  & ! u_{j up}
                &       - conjg(wf(4*j-2, r))*wf(4*j-2, r)*dfermi(eg(r) - e_value(i), kbt)  & ! u_{j down}
                &       - conjg(wf(4*j-1, r))*wf(4*j-1, r)*dfermi(eg(r) + e_value(i), kbt)  & ! v_{j up}
                &       - conjg(wf(4*j  , r))*wf(4*j  , r)*dfermi(eg(r) + e_value(i), kbt)    ! v_{j down}
      end do
    end do
  end subroutine stm_simple_rho_vs_e
  
  !subroutine measurement_ne_site(sys, wf, eg, kbt, ne_site)
  !  integer, intent(in) :: site
  !  complex(8), dimension(:,:), intent(in) :: wf
  !  real(8), dimension(:), intent(in) :: eg
  !  real(8), intent(in) :: kbt
  !  real(8), dimension(:), intent(in) :: e_value
  !  real(8), dimension(size(e_value)), intent(out) :: rho
  !  integer :: i, j, k, r
  !  rho(:) = 0d0
  !  do i = 1, size(e_value)
  !    do r = size(eg)/2 + 1, size(eg)
  !      j = site
  !      rho(i) = rho(i) - conjg(wf(4*j-3, r))*wf(4*j-3, r)*dfermi(eg(r) - e_value(i), kbt)  & ! u_{j up}
  !              &       - conjg(wf(4*j-2, r))*wf(4*j-2, r)*dfermi(eg(r) - e_value(i), kbt)  & ! u_{j down}
  !              &       - conjg(wf(4*j-1, r))*wf(4*j-1, r)*dfermi(eg(r) + e_value(i), kbt)  & ! v_{j up}
  !              &       - conjg(wf(4*j  , r))*wf(4*j  , r)*dfermi(eg(r) + e_value(i), kbt)    ! v_{j down}
  !    end do
  !  end do
  !end subroutine stm_simple_rho_vs_e
end module measurement_mod