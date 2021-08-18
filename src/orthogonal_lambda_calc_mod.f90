MODULE orthogonal_lambda_calc_mod
  USE io_mod, ONLY:save_txt
  USE math_mod, ONLY:ui, pi, my_zgetrs2
  use parameter_mod
  use base_mod
  IMPLICIT NONE
    
  integer, parameter, private :: orthogonal_lambda_calc_iter_max = 10000
  real(8), parameter, private :: orthogonal_lambda_calc_tol = 1d-10
  
  private orthogonal_lambda_calc_S
  private orthogonal_lambda_calc_P
  private orthogonal_lambda_calc_lmlm
  !private orthogonal_lambda_calc_prod
  private orthogonal_lambda_set_E
  private orthogonal_lambda_calc_lambda

  public orthogonal_lambda_get_new_wf

contains
  
  subroutine orthogonal_lambda_calc_S(S, bar_wf, wf)
    complex(8), dimension(:, :), intent(out) :: S
    complex(8), dimension(:, :), intent(in) :: bar_wf, wf
    integer i, j, n
    S(:,:) = 0d0
    n = size(bar_wf, 1)
    do i=1, n
      do j=1, n
        S(i,j) = S(i,j) + dot_product(bar_wf(:,i), wf(:,j))
      end do
    end do
  end subroutine orthogonal_lambda_calc_S

  subroutine orthogonal_lambda_calc_P(P, bar_wf)
    complex(8), dimension(:, :), intent(out) :: P
    complex(8), dimension(:, :), intent(in) :: bar_wf
    integer i, j, n
    P(:,:) = 0d0
    n = size(bar_wf, 1)
    do i=1, n
      do j=1, n
        P(i,j) = P(i,j) + dot_product(bar_wf(:,i),bar_wf(:,j))
      end do
    end do
  end subroutine orthogonal_lambda_calc_P
  
  subroutine orthogonal_lambda_calc_lmlm(lmlm, lm)
    complex(8), dimension(:, :), intent(out) :: lmlm
    complex(8), dimension(:, :), intent(in) :: lm
    complex(8), dimension(size(lm), size(lm)) :: trlm
    lmlm = 0d0
    trlm = transpose(lm)
    lmlm = matmul(conjg(lm),trlm)
  end subroutine orthogonal_lambda_calc_lmlm
  
  !function orthogonal_lambda_calc_prod(lD, rD) result(r)
  !  complex(8), dimension(:), intent(in) :: lD, rD
  !  complex(8)                          :: r
  !  integer i
  !  r = 0d0
  !  do i=1, this%n
  !    r = r + conjg(lD(2*i-1))*rD(2*i-1) + conjg(lD(2*i))*rD(2*i)
  !  end do
  !end function orthogonal_lambda_calc_prod

  subroutine orthogonal_lambda_set_E(E)
    complex(8), dimension(:, :), intent(out) :: E
    integer i, n
    E(:, :) = 0d0
    n = size(E, 1)
    do i=1, n
      E(i, i) = 1.0d0
    end do
  end subroutine orthogonal_lambda_set_E
  
  !subroutine orthogonal_lambda_calc_lambda(bar_wf, wf, lambda, init_lm)
  subroutine orthogonal_lambda_calc_lambda(bar_wf, wf, lambda)
    complex(8), dimension(:, :), intent(in) :: bar_wf, wf
    complex(8), dimension(:, :), intent(out) :: lambda
    !complex(8), dimension(:, :), intent(in), optional :: init_lm
    complex(8), dimension(size(wf), size(wf)) :: old_lambda, S, P, lmlm, E, b
    integer i

    lambda = 0d0

    call orthogonal_lambda_calc_S(S, bar_wf, wf)
    call orthogonal_lambda_calc_P(P, bar_wf)
    call orthogonal_lambda_set_E(E)
    !if(present(init_lm)) then
    !  old_lambda = init_lm
    !else
      !!b = 0.5d0*(E - P)  ! initial lm is 0
      !!call my_zgetrs2(S, b, old_lambda)
      !!old_lambda = transpose(old_lambda)
      !old_lambda = 0d0
      old_lambda = E
    !end if
    do i=1, orthogonal_lambda_calc_iter_max
      call orthogonal_lambda_calc_lmlm(lmlm, old_lambda)
      b = 0.5d0*(E - P - lmlm)
      call my_zgetrs2(S, b, lambda)
      lambda = transpose(lambda)
      !print *, "iter is ", i, maxval(abs(old_lambda(:,:) - lambda(:,:)))
      if(maxval(abs(old_lambda(:,:) - lambda(:,:))) < orthogonal_lambda_calc_tol) exit
      old_lambda(:,:) = lambda(:,:)
    end do
  end subroutine orthogonal_lambda_calc_lambda
  
  !subroutine orthogonal_lambda_get_new_wf(new_wf, bar_wf, wf, lambda)
  subroutine orthogonal_lambda_get_new_wf(new_wf, bar_wf, wf)
    complex(8), dimension(:, :), intent(out) :: new_wf
    complex(8), dimension(: ,:), intent(in) :: bar_wf, wf
    !complex(8), dimension(:, :), intent(in) :: lambda
    complex(8), dimension(size(wf), size(wf)) :: lambda
    integer i, j, n
    new_wf(:, :) = bar_wf(:, :)
    call orthogonal_lambda_calc_lambda(bar_wf, wf, lambda)
    n = size(wf, 1)
    do i=1, n
      do j=1, n
        new_wf(:, i) = new_wf(:, i) + lambda(i, j)*wf(:, j)
      end do
    end do
  end subroutine orthogonal_lambda_get_new_wf

END MODULE orthogonal_lambda_calc_mod