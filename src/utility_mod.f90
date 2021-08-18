module utility_mod
  implicit none
  !!!  !$ use omp_lib
  contains
  subroutine gram_schmidt(wf, orthogonal_num)
    complex(8), DIMENSION(:, :), intent(inout)   :: wf
    complex(8), DIMENSION(size(wf(:,1)))   :: twf
    real(8) :: vec2
    integer, intent(in), optional :: orthogonal_num
    integer n_oth, k, l
    n_oth = size(wf(1, :))
    if(present(orthogonal_num)) then
      if(orthogonal_num < n_oth) then
        n_oth = orthogonal_num
      end if
    end if
    DO k = 1, n_oth
      twf = wf(:,k)
      DO l = 1, k - 1
        twf = twf - DOT_PRODUCT(wf(:, l), wf(:, k))*wf(:, l) ! Classic Gram-Schmidt
      ENDDO!l
      vec2 = DOT_PRODUCT(twf, twf)
      IF (vec2 == 0.d0) THEN
        wf(:,k) = 0.d0
      ELSE
        wf(:,k) = twf*vec2**(-0.5d0)
      ENDIF
    ENDDO!k
  end subroutine gram_schmidt

  subroutine modified_gram_schmidt(wf, orthogonal_num)
    complex(8), DIMENSION(:, :), intent(inout)   :: wf
    complex(8), DIMENSION(size(wf(:,1)))   :: twf
    real(8) :: vec2
    integer, intent(in), optional :: orthogonal_num
    integer n_oth, k, l
    n_oth = size(wf(1, :))
    if(present(orthogonal_num)) then
      if(orthogonal_num < n_oth) then
        n_oth = orthogonal_num
      end if
    end if
    DO k = 1, n_oth
      twf = wf(:,k)
      DO l = 1, k - 1
        twf = twf - DOT_PRODUCT(wf(:, l), twf)*wf(:, l)
      ENDDO!l
      vec2 = DOT_PRODUCT(twf, twf)
      IF (vec2 == 0.d0) THEN
        wf(:,k) = 0.d0
      ELSE
        wf(:,k) = twf*vec2**(-0.5d0)
      ENDIF
    ENDDO!k
  end subroutine modified_gram_schmidt
  
  subroutine gram_schmidt_2(wf)
    complex(8), DIMENSION(:, :), intent(inout)   :: wf
    complex(8), DIMENSION(size(wf,1))   :: twf
    real(8) :: vec2
    integer n_oth, k, l
    n_oth = size(wf(1, :))
    DO k = 1, n_oth
      twf = wf(:,k)
      DO l = 1, k - 3, 4
        twf = twf - DOT_PRODUCT(wf(:, l), wf(:, k))*wf(:, l)
        twf = twf - DOT_PRODUCT(wf(:, l+1), wf(:, k))*wf(:, l+1)
        twf = twf - DOT_PRODUCT(wf(:, l+2), wf(:, k))*wf(:, l+2)
        twf = twf - DOT_PRODUCT(wf(:, l+3), wf(:, k))*wf(:, l+3)
      ENDDO!l
      do l = 4*int(k/4) + 1, k
        twf = twf - DOT_PRODUCT(wf(:, l), wf(:, k))*wf(:, l)
      end do
      vec2 = DOT_PRODUCT(twf, twf)
      IF (vec2 == 0.d0) THEN
        wf(:,k) = 0.d0
      ELSE
        wf(:,k) = twf*vec2**(-0.5d0)
      ENDIF
    ENDDO!k
  end subroutine gram_schmidt_2
  
  subroutine gram_schmidt_3(wf)
    complex(8), DIMENSION(:, :), intent(inout)   :: wf
    complex(8), DIMENSION(size(wf,1))   :: twf
    real(8) :: vec2
    integer n_oth, k, l
    n_oth = size(wf(1, :))
    DO k = 1, n_oth
      !twf = wf(:,k)
      twf = 0d0
      !$omp parallel
      !$omp do reduction(+:twf)
      DO l = 1, k-1
        twf = twf + DOT_PRODUCT(wf(:, l), wf(:, k))*wf(:, l)
      ENDDO!l
      !$omp end do
      !$omp end parallel

      twf = wf(:, k) - twf
      !print *,  DOT_PRODUCT(twf, twf)
      
      !twf = wf(:,k)
      !!twf = 0d0
      !DO l = 1, k
      !  twf = twf - DOT_PRODUCT(wf(:, l), wf(:, k))*wf(:, l)
      !ENDDO!l
      !!twf = wf(:, k) - twf
      !print *,  DOT_PRODUCT(twf, twf)
      vec2 = DOT_PRODUCT(twf, twf)
      IF (vec2 == 0.d0) THEN
        wf(:,k) = 0.d0
      ELSE
        wf(:,k) = twf*vec2**(-0.5d0)
      ENDIF
    ENDDO!k

  end subroutine gram_schmidt_3
  
  subroutine gram_schmidt_3_half(wf)
    complex(8), DIMENSION(:, :), intent(inout)   :: wf
    complex(8), DIMENSION(size(wf,1))   :: twf
    real(8) :: vec2
    integer n_oth, k, l
    n_oth = size(wf(1, :))/2
    DO k = 1, n_oth
      twf = 0d0
      !$omp parallel
      !$omp do reduction(+:twf)
      DO l = 1, k-1
        twf = twf + DOT_PRODUCT(wf(:, l), wf(:, k))*wf(:, l)
      ENDDO!l
      !$omp end do
      !$omp end parallel
      twf = wf(:, k) - twf
      vec2 = DOT_PRODUCT(twf, twf)
      IF (vec2 == 0.d0) THEN
        wf(:,k) = 0.d0
      ELSE
        wf(:,k) = twf*vec2**(-0.5d0)
      ENDIF
    ENDDO!k
  end subroutine gram_schmidt_3_half
  
  subroutine gram_schmidt_4(wf)
    complex(8), DIMENSION(:, :), intent(inout)   :: wf
    !complex(8), DIMENSION(:, :), allocatable   :: twf
    real(8) :: vec2
    integer n_oth, k, l
    n_oth = size(wf(1, :))
    !allocate(twf(n_oth, n_oth), source = wf)
    do k = 1, n_oth
      vec2 = dot_product(wf(:, k), wf(:, k))
      !if (vec2 == 0.d0) then
      !  wf(:,k) = 0.d0
      !else
        wf(:,k) = wf(:, k)*vec2**(-0.5d0)
        !wf(:,k) = wf(:, k)/sqrt(vec2)
      !endif
      !$omp parallel
      !$omp do
      do l = k+1, n_oth
        wf(:, l) = wf(:, l) - dot_product(wf(:, l), wf(:, k))*wf(:, l)
      enddo!l
      !$omp end do
      !$omp end parallel
    end do ! k
    !do k = n_oth / 2 + 1, n_oth
    !  !$omp parallel
    !  !$omp do reduction(-:twf)
    !  DO l = n_oth / 2 + 1, k-1
    !    twf = twf + DOT_PRODUCT(wf(:, l), wf(:, k))*wf(:, l)
    !  ENDDO!l
    !  !$omp end do
    !  !$omp end parallel
    !  vec2 = DOT_PRODUCT(twf, twf)
    !  IF (vec2 == 0.d0) THEN
    !    wf(:,k) = 0.d0
    !  ELSE
    !    wf(:,k) = twf*vec2**(-0.5d0)
    !  ENDIF
    !end do
  end subroutine gram_schmidt_4
  
  subroutine gram_schmidt_5(wf)
    complex(8), DIMENSION(:, :), intent(inout)   :: wf
    complex(8), DIMENSION(size(wf,1))   :: twf
    real(8) :: vec2
    integer n_oth, len, k, l, i
    n_oth = size(wf(1, :))
    len = size(wf(:, 1))
    DO k = 1, n_oth
      !$omp parallel
      !$omp do
      do i = 1, len
        twf(i) = wf(i, k)
      end do
      !$omp end do
      !$omp do reduction(-:twf)
      DO l = 1, k-1
        twf = twf + DOT_PRODUCT(wf(:, l), wf(:, k))*wf(:, l)
      ENDDO!l
      !$omp end do
      !!$omp single
      !!vec2 = 0d0
      !!$omp end single
      !!$omp do reduction(+:vec2)
      !do i = 1, len
      !  vec2 = vec2 + conjg(twf(i))*twf(i)
      !end d
      vec2 = dot_product(twf,twf) ** (-0.5d0)
      !!$omp end do
      !!$omp single
      !vec2 = vec2**(-0.5d0)
      !!$omp end single
      !!$omp do
      do i = 1, len
        wf(i,k) = twf(i)*vec2
      end do
      !!$omp end do
      
      !$omp end parallel
      !IF (vec2 == 0.d0) THEN
      !  wf(:,k) = 0.d0
      !ELSE
      !  wf(:,k) = twf*vec2**(-0.5d0)
      !ENDIF
    enddo!k

  end subroutine gram_schmidt_5
  
  pure function check_sorted_eg(eg) result(r)
    real(8), dimension(:), intent(in) :: eg
    logical :: r
    integer i
    do i=1, size(eg)-1
      if(eg(i) > eg(i+1)) then
        r = .false.
        return
      end if
    end do
    r = .true.
  end function check_sorted_eg
  
  subroutine sort_wf(eg, wf, oeg, owf)
    real(8), dimension(:), intent(inout) :: eg
    complex(8), dimension(size(eg),size(eg)), intent(inout) :: wf
    real(8), dimension(size(eg)), intent(inout), optional :: oeg
    complex(8), dimension(size(eg),size(eg)), intent(inout), optional :: owf
    real(8), dimension(size(eg)) :: teg
    complex(8), dimension(size(eg),size(eg)) :: twf
    integer, dimension(size(eg)) :: id
    integer i
    teg = eg
    twf = wf
    !call get_bsort_id_eg(eg,id)
    call get_hsort_id_eg(eg, id)
    do i = 1, size(eg)
      eg(i) = teg(id(i))
      wf(:, i) = twf(:, id(i))
    end do
    if(present(oeg)) then
      teg = oeg
      do i = 1, size(eg)
        oeg(i) = teg(id(i))
      end do
    end if
    if(present(owf)) then
      twf = owf
      do i = 1, size(eg)
        owf(:, i) = twf(:, id(i))
      end do
    end if
  end subroutine sort_wf
  
  subroutine get_bsort_id_eg(eg, id)
    ! bubble sort for eg
    ! id is index of sorted eg
    real(8), dimension(:), intent(in) :: eg
    integer, dimension(size(eg)), intent(out) :: id
    integer :: i, j, n, tmpid
    real(8), dimension(size(eg)) :: teg
    real(8) :: tmpeg
    n = size(eg)
    teg = eg
    do i = 1, n
      id(i) = i
    end do
    do i = n-1, 1, -1
      do j = 1, i
        if(teg(j) > teg(j+1)) then
          tmpeg = teg(j)
          teg(j) = teg(j+1)
          teg(j+1) = tmpeg
          tmpid = id(j)
          id(j) = id(j+1)
          id(j+1) = tmpid
        end if
      end do
    end do
  end subroutine get_bsort_id_eg
  
  subroutine get_hsort_id_eg(eg, id)
    ! heap sort for eg
    ! id is index of sorted eg
    use math_mod, only : swap
    real(8), dimension(:), intent(in) :: eg
    integer, dimension(size(eg)), intent(out) :: id
    integer :: i, n
    real(8), dimension(size(eg)) :: heg
    n = size(eg)
    heg = eg
    do i = 1, n
      id(i) = i
    end do
    do i = n/2, 1, -1
      call pushdown_heap_tree(heg, id, i, n)
    end do
    do i = n, 2, -1
      call swap(heg(1), heg(i))
      call swap(id(1), id(i))
      call pushdown_heap_tree(heg, id, 1, i-1)
    end do
  end subroutine get_hsort_id_eg

  subroutine pushdown_heap_tree(htree, id, first, last)
    use math_mod, only : swap
    real(8), dimension(:), intent(inout) :: htree
    integer, dimension(:), intent(inout) :: id
    integer, intent(in) :: first, last
    integer p, c
    p = first
    c = 2*p
    do while (c <= last)
      !if (c < last .and. htree(c) < htree(c+1) ) then
      !if (c < last .and. htree(c) < htree(c+1) &
      !  & .or. (htree(c) == htree(c+1) .and. id(c) < id(c+1)) ) then
      if (c < last) then
        if(htree(c) < htree(c+1) &
          & .or. (htree(c) == htree(c+1) .and. id(c) < id(c+1)) ) then
          c = c+1;
        end if
      end if
      !if (htree(c) <= htree(p) ) exit
      if (htree(c) < htree(p) &
        & .or. (htree(c) == htree(p) .and. id(c) < id(p)) ) exit
      call swap(htree(c), htree(p))
      call swap(id(c), id(p))
      p = c;
      c = 2*p;
    end do
  end subroutine pushdown_heap_tree
  
end module utility_mod