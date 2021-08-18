module winding_number_mod
  use math_mod,        only : pi
  use circle_list_mod,  only : circle_list
  use parameter_mod,  only : site_max, site_to_r
  implicit none
  
  real(8), parameter ::  winding_number_sum_error = 1d-8

  contains

  subroutine get_winding_number(var, clist, wl)
    real(8), dimension(:), intent(in) :: var
    class(circle_list), intent(in) :: clist
    integer, dimension(:), allocatable, intent(inout) :: wl
    integer, dimension(:),allocatable :: loop_site
    integer :: i, j
    if(allocated(wl)) deallocate(wl)
    allocate(wl(size(clist%cil)),source = 0)
    do i=1, size(clist%cil)
      allocate(loop_site(size(clist%cil(i)%value)))
      do j=1, size(clist%cil(i)%value)
        if(clist%cil(i)%value(j) > 0) then
          loop_site(j) = clist%cpl%value( clist%cil(i)%value(j) )%i
        else
          loop_site(j) = clist%cpl%value( -clist%cil(i)%value(j) )%f
        end if
      end do
      call get_winding_number_loop(var, loop_site, wl(i))
      deallocate(loop_site)
    end do
    !print *,wl
  end subroutine get_winding_number
  
  subroutine get_winding_number_loop(var, loop_site, wl)
    real(8), dimension(:), intent(in) :: var
    integer, dimension(:), intent(in) :: loop_site
    integer, intent(out)              :: wl
    integer, dimension(size(loop_site)+1) :: cl
    integer i
    real(8) sum, tmp
    cl(1:size(loop_site)) = loop_site
    cl(size(loop_site)+1) = loop_site(1)
    sum = 0d0
    do i=1, size(loop_site)
      tmp =  modulo(var(cl(i+1)) - var(cl(i)), 2*pi)
      if(tmp >= pi) then !!!!! -pi <= delta var < pi 
        tmp =  tmp - 2*pi
      end if
      sum = sum + tmp
    end do
    if(sum > pi/2d0) then
      sum = sum + winding_number_sum_error
    else if(sum < -pi/2d0) then
      sum = sum - winding_number_sum_error
    end if
    wl = sum / (2*pi)
  end subroutine get_winding_number_loop

  subroutine get_winding_point_surface(var, clist, pole, wn)
    real(8), dimension(:), intent(in) :: var
    class(circle_list), intent(in) :: clist
    real(8), dimension(:,:), allocatable, intent(inout) :: pole
    integer, dimension(:), allocatable, intent(inout) :: wn
    integer, dimension(:), allocatable :: wl
    integer :: i, j, k
    real(8) :: p(3)
  
    if(allocated(pole)) deallocate(pole)
    if(allocated(wn)) deallocate(wn)
    
    allocate(wl(size(clist%cil)),source = 0)
    call get_winding_number(var, clist, wl)
    
    k = 0
    do i = 1, size(clist%cil)
      p = clist%get_circle_center(i)
      !if(wl(i) /= 0 .and. nint(p(3)) == 1) then  ! nint( 1.49999... ) is 1. 0.5 <= p(3) < 1.5
      if(wl(i) /= 0 .and. nint(p(3)*10) == 10) then ! 0.95 <= p(3) < 1.05
        k = k + 1;
      end if
    end do
    allocate(pole(2, k))
    allocate(wn(k))
    pole(:, :) = 0d0
    k = 0
    do i = 1, size(clist%cil)
      p = clist%get_circle_center(i)
      !if(wl(i) /= 0 .and. nint(p(3)) == 1) then
      if(wl(i) /= 0 .and. nint(p(3)*10) == 10) then
        k = k + 1
        pole(:, k) = p(1:2)
        wn(k) = wl(i)
      end if
    end do

  end subroutine get_winding_point_surface

  subroutine get_chi_winding_number_with_eta_winding_point_surface(eta, chi, clist, pole, wn)
    real(8), dimension(:), intent(in) :: eta, chi
    class(circle_list), intent(in) :: clist
    real(8), dimension(:,:), allocatable, intent(inout) :: pole
    integer, dimension(:), allocatable, intent(inout) :: wn
    integer, dimension(:), allocatable :: wl_eta, wl_chi
    integer :: i, j, k
    real(8) :: p(3)
  
    if(allocated(pole)) deallocate(pole)
    if(allocated(wn)) deallocate(wn)
    
    allocate(wl_eta(size(clist%cil)),source = 0)
    allocate(wl_chi(size(clist%cil)),source = 0)
    call get_winding_number(eta, clist, wl_eta)
    call get_winding_number(chi, clist, wl_chi)
    
    k = 0
    do i = 1, size(clist%cil)
      p = clist%get_circle_center(i)
      if(wl_eta(i) /= 0 .and. nint(p(3)*10) == 10) then ! 0.95 <= p(3) < 1.05
        k = k + 1;
      end if
    end do
    allocate(pole(2, k))
    allocate(wn(k))
    pole(:, :) = 0d0
    k = 0
    do i = 1, size(clist%cil)
      p = clist%get_circle_center(i)
      if(wl_eta(i) /= 0 .and. nint(p(3)*10) == 10) then
        k = k + 1
        pole(:, k) = p(1:2)
        wn(k) = wl_chi(i)
      end if
    end do

  end subroutine get_chi_winding_number_with_eta_winding_point_surface

  subroutine get_eta(eta, xi, pshape)
    real(8), dimension(:), allocatable, intent(inout) :: eta
    real(8), dimension(:), intent(in) :: xi
    integer, dimension(:,:), intent(in)  :: pshape
    integer :: i, n, sr(3)
    n = site_max(pshape)
    if(allocated(eta)) deallocate(eta)
    allocate(eta(n), source = 0d0)
    do i = 1, n
      sr = site_to_r(i, pshape)
      !eta(i) = (sr(1)+sr(2))*pi + xi(i)
      eta(i) = (sr(1)+sr(2)+sr(3))*pi + xi(i)
    end do
  end subroutine
end module winding_number_mod
