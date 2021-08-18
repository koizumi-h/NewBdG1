module hop_iterator_mod
  use path_list_mod, only : path_list,path
  use parameter_mod
  implicit none
  type, public :: hop_iterator
    integer, dimension(:,:), allocatable :: pshape
    type(path_list)       :: pl, ss
  contains
    ! public type-bound procedures
    procedure :: hop_iterator_init
    procedure :: set_pshape => hop_iterator_set_pshape
    procedure :: set_path   => hop_iterator_set_path
    procedure :: set_path_hole => hop_iterator_set_path_hole
    procedure :: set_path_site => hop_iterator_set_path_site
    procedure :: next       => hop_iterator_next
    procedure :: has_next   => hop_iterator_has_next
    procedure :: delete_same_path => hop_iterator_delete_same_path
  end type hop_iterator
contains
  ! public type-bound procedures
  pure subroutine hop_iterator_init(this, pshape)
    class(hop_iterator), intent(inout) :: this
    integer, dimension(:,:), intent(in)  :: pshape
    call this%set_pshape(pshape)
  end subroutine hop_iterator_init

  pure subroutine hop_iterator_set_pshape(this, pshape)
    class(hop_iterator), intent(inout) :: this
    integer, dimension(:,:), intent(in)  :: pshape
    if(allocated(this%pshape)) deallocate(this%pshape)
    allocate(this%pshape(2,size(pshape(1,:))))
    this%pshape(:,:) = pshape(:,:)
  end subroutine hop_iterator_set_pshape

  !pure subroutine hop_iterator_set_path(this, hop_dir, hole_pos)
  subroutine hop_iterator_set_path(this, hop_dir, hole_pos)
  !subroutine hop_iterator_set_path(this, hop_dir, hole_pos)
    class(hop_iterator), intent(inout)          :: this
    integer, dimension(:), intent(in)           :: hop_dir
    integer                                     :: i, j
    integer, dimension(:), intent(in), optional :: hole_pos
!    call this%pl%realloc_value()
    call this%pl%alloc_value()
    do i = 1, size(hop_dir)
!      if (hop_dir(i) < 1 .or. hop_dir(i) > 12) cycle
      call this%pl%extend(get_path_list(this%pshape, hop_dir(i), hole_pos))
    end do
    !if (present(hole_pos).and.size(hole_pos,1)>0) then
    if (present(hole_pos)) then
      do i = 1, size(hole_pos)
        call this%pl%remove(hole_pos(i))
      end do
    end if
    ! delete same paths
    call this%delete_same_path()
    !i = 1
    !do
    !  i = i + 1
    !  !print '(I5,I5)', i, this%pl%length()
    !  if(i > this%pl%length() )then 
    !    exit
    !  end if
    !  do j = 1, i - 1
    !    !write(*,*) i,j,this%pl%value(i),this%pl%value(j) 
    !    if(((this%pl%value(i)%i == this%pl%value(j)%i)          &
    !      &   .and. (this%pl%value(i)%f == this%pl%value(j)%f)) &
    !      & .or. ((this%pl%value(i)%f == this%pl%value(j)%i)      &
    !      &   .and. (this%pl%value(i)%i == this%pl%value(j)%f))) then
    !      call this%pl%del(i)
    !      !write(*,*) "ddddddddddddddddddddddddddddd"
    !      i = i - 1
    !      exit
    !    end if
    !  end do
    !end do
    !stop


  end subroutine hop_iterator_set_path

  pure subroutine hop_iterator_set_path_hole(this, hole_pos)
    class(hop_iterator), intent(inout)          :: this
    integer, dimension(:), intent(in)           :: hole_pos
    integer, dimension(4)                       :: hnp !hole_nearest_pos
    integer                                     :: i, j, k, z
    type(path)                                  :: p
!    call this%pl%realloc_value()
    call this%pl%alloc_value()
    do i = 1, size(hole_pos)
      z = site_to_z(hole_pos(i),this%pshape)
      hnp(1) = hole_pos(i) - this%pshape(1,z) !h-y
      hnp(2) = hole_pos(i) - 1              !h-x
      hnp(3) = hole_pos(i) + 1              !h+x
      hnp(4) = hole_pos(i) + this%pshape(1,z) !h+y
      !print "('hole :',I0,'(',3(x1I0),')')",hole_pos(i),site_to_r(hole_pos(i),this%pshape)
      !print *, z
      !print "('hnp = ',4(x1I0))", hnp
      do j = 1, 4
        do k = 1, 4-j
          p%i = hnp(j)
          p%f = hnp(k+j)
!          if (p%f > p%i) call this%pl%append(p)
          call this%pl%append(p)
        end do
      end do
    end do
  end subroutine hop_iterator_set_path_hole

subroutine hop_iterator_set_path_site(this, site_pos, hole_pos)
    class(hop_iterator), intent(inout)          :: this
    integer, dimension(:), intent(in)           :: site_pos, hole_pos
    integer, dimension(4)                       :: snp !site_nearest_pos
    integer                                     :: i, j, k
    type(path)                                  :: p
!    call this%pl%realloc_value()
    call this%pl%alloc_value()
    do i = 1, size(site_pos)
      p%i = site_pos(i)
      snp(1) = site_pos(i) - this%pshape(1,1) !h-y
      if(snp(1) < 1 ) snp(1) = -1
      snp(2) = site_pos(i) - 1              !h-x
      if(i == 1 .or. snp(2) == site_pos(i-1)) snp(2) = -1
      snp(3) = site_pos(i) + 1              !h+x
      if(i == size(site_pos) .or. snp(2) == site_pos(i+1)) snp(3) = -1
      snp(4) = site_pos(i) + this%pshape(1,1) !h+y
      if(snp(4) > size(site_pos)) snp(4) = -1
      do j = 1, 4
        do k = 1, size(hole_pos)
          if(snp(j) == hole_pos(k)) snp(j) = -2
        enddo
          p%f = snp(j)
!          if (p%f > p%i) call this%pl%append(p)
          call this%ss%append(p)
      end do
    end do
    write(*,*)  'path_site', this%ss%value
  end subroutine hop_iterator_set_path_site

  pure subroutine hop_iterator_next(this, p)
    class(hop_iterator), intent(inout) :: this
    type(path), intent(out)            :: p
    if (this%pl%length() == 0) then
      p = path(-1, -1)
    else
      call this%pl%pop_sub(p)
    end if
  end subroutine hop_iterator_next

  pure subroutine hop_iterator_delete_same_path(this, cmp_hop_iter)
    class(hop_iterator), intent(inout) :: this
    class(hop_iterator), intent(in), optional :: cmp_hop_iter
   integer i,j 
    if(present(cmp_hop_iter)) then
      i = 1
      do
        i = i + 1
        if(i > this%pl%length() )then 
          exit
        end if
        do j = 1, cmp_hop_iter%pl%length()
          if(((this%pl%value(i)%i == cmp_hop_iter%pl%value(j)%i)          &
            &   .and. (this%pl%value(i)%f == cmp_hop_iter%pl%value(j)%f)) &
            & .or. ((this%pl%value(i)%f == cmp_hop_iter%pl%value(j)%i)      &
            &   .and. (this%pl%value(i)%i == cmp_hop_iter%pl%value(j)%f))) then
            !write(*,*) "delete : ", this%pl%value(i)
            call this%pl%del(i)
            i = i - 1
            exit
          end if
        end do
      end do
    else
      i = 1
      do
        i = i + 1
        !print '(I5,I5)', i, this%pl%length()
        if(i > this%pl%length() )then 
          exit
        end if
        do j = 1, i - 1
          !write(*,*) i,j,this%pl%value(i),this%pl%value(j) 
          if(((this%pl%value(i)%i == this%pl%value(j)%i)          &
            &   .and. (this%pl%value(i)%f == this%pl%value(j)%f)) &
            & .or. ((this%pl%value(i)%f == this%pl%value(j)%i)      &
            &   .and. (this%pl%value(i)%i == this%pl%value(j)%f))) then
            !write(*,*) "delete : ", this%pl%value(i)
            call this%pl%del(i)
            i = i - 1
            exit
          end if
        end do
      end do
    end if

  end subroutine hop_iterator_delete_same_path

  pure logical function hop_iterator_has_next(this) result(r)
    class(hop_iterator), intent(in) :: this
    if (this%pl%length() == 0) then
      r = .false.
    else
      r = .true.
    end if
  end function hop_iterator_has_next

  ! this will delete
  type(path_list) function get_hopping_path_(pshape, hop_dir, hole_pos) result(r)
    integer, dimension(:), intent(in)           :: pshape
    integer, dimension(:), intent(in)           :: hop_dir
    integer, dimension(:), intent(in), optional :: hole_pos
    type(hop_iterator), allocatable             :: iter
    type(path)                                  :: p
    r = get_hopping_path(reshape(pshape,[2,1]),hop_dir,hole_pos)
  end function get_hopping_path_
  
  ! module subprograms
  type(path_list) function get_hopping_path_old(pshape, hop_dir, hole_pos, except_dir) result(r)
    integer, dimension(:,:), intent(in)           :: pshape
    integer, dimension(:), intent(in)           :: hop_dir
    integer, dimension(:), intent(in), optional :: hole_pos
    integer, dimension(:), intent(in), optional :: except_dir
    type(hop_iterator), allocatable             :: iter
    type(hop_iterator), allocatable             :: except_iter
    type(path)                                  :: p
    allocate(iter, source = new_hop_iterator(pshape))
    call iter%set_path(hop_dir, hole_pos)
    call r%alloc_value()
    do while(iter%has_next())
      call iter%next(p)
      call r%append(p)
    end do
    if(present(except_dir)) then
      allocate(except_iter, source = new_hop_iterator(pshape))
      call except_iter%set_path(except_dir, hole_pos)
      call iter%delete_same_path(except_iter)
    end if
    deallocate(iter)
  end function get_hopping_path_old
  type(path_list) function get_hopping_path(pshape, hop_dir, hole_pos) result(r)
    integer, dimension(:,:), intent(in)           :: pshape
    integer, dimension(:), intent(in)           :: hop_dir
    integer, dimension(:), intent(in), optional :: hole_pos
    !type(hop_iterator), allocatable             :: iter
    !type(hop_iterator), allocatable             :: except_iter
    type(path)                                  :: p
    integer i
    !allocate(iter, source = new_hop_iterator(pshape))
    !call iter%set_path(hop_dir, hole_pos)
    call r%alloc_value()
    !do while(iter%has_next())
    !  call iter%next(p)
    !  call r%append(p)
    !end do
    do i = 1, size(hop_dir)
      call r.extend(get_path_list(pshape, hop_dir(i), hole_pos))
    end do
    if (present(hole_pos)) then
      do i = 1, size(hole_pos)
        call r%remove(hole_pos(i))
      end do
    end if
    ! delete same paths
    call delete_same_path(r)
  end function get_hopping_path
  pure subroutine delete_same_path(pl)
    class(path_list), intent(inout) :: pl
    integer i,j 
    i = 1
    do
      i = i + 1
      !print '(I5,I5)', i, this%pl%length()
      if(i > pl%length() )then 
        exit
      end if
      do j = 1, i - 1
        !write(*,*) i,j,this%pl%value(i),this%pl%value(j) 
        if(((pl%value(i)%i == pl%value(j)%i)          &
          &   .and. (pl%value(i)%f == pl%value(j)%f)) &
          & .or. ((pl%value(i)%f == pl%value(j)%i)      &
          &   .and. (pl%value(i)%i == pl%value(j)%f))) then
          !write(*,*) "delete : ", this%pl%value(i)
          call pl%del(i)
          i = i - 1
          exit
        end if
      end do
    end do
  end subroutine delete_same_path
  !type(path_list) function get_hopping_path_l(layers, hop_dir) result(r)
  !  type(cuo2_plane),intent(in) :: layers(:)
  !  integer, dimension(:), intent(in)           :: hop_dir
  !  integer, dimension(:,:),allocatable         :: pshape
  !  integer, dimension(:),allocatable          :: hole_pos
  !  pshape = make_pshape(layers)
  !  hole_pos = make_holearray(layers)
  !  r = get_hopping_path(pshape, hop_dir, hole_pos)
  !  deallocate(pshape)
  !  deallocate(hole_pos)
  !end function get_hopping_path_l
  !pure type(path_list) function get_path_list(pshape, i, hole_pos) result(r)
  type(path_list) function get_path_list(pshape, i, hole_pos) result(r)
    integer, dimension(:,:), intent(in) :: pshape
    integer, intent(in)               :: i
    integer, dimension(:), intent(in), optional :: hole_pos
    !print *, "size : ",  size(pshape(1,:))
    select case(i)
      case(1);  r = get_path_list1(pshape)
      case(2);  r = get_path_list2(pshape)
      case(3);  r = get_path_list3(pshape,hole_pos)
      case(4);  r = get_path_list4(pshape,hole_pos)
      !case(5);  r = get_path_list_xy(pshape)
      case(5);  r = get_path_list_xy1(pshape)
      !case(6);  r = get_path_list_mxy(pshape)
      case(6);  r = get_path_list_mxy1(pshape)
      case(7);  r = get_path_list_z(pshape)
    end select
  end function get_path_list

  pure type(path_list) function get_path_list1(pshape) result(r)
    integer, dimension(:,:), intent(in) :: pshape
    integer                           :: i, k, n
    call r%alloc_value()
    n  = 0
    do k = 1, size(pshape(1,:))
      do i = 1, pshape(1,k)*pshape(2,k)
        if(mod(i,pshape(1,k)) == 0) cycle
        call r%append(path(i+n, i+n+1))
      end do
      n = n + pshape(1,k)*pshape(2,k)
    end do
  end function get_path_list1

  pure type(path_list) function get_path_list2(pshape) result(r)
    integer, dimension(:,:), intent(in) :: pshape
    integer                           :: i, k, n
    call r%alloc_value()
    n  = 0
    do k = 1, size(pshape(1,:))
      do i = 1, pshape(1,k)*(pshape(2,k)-1)
        call r%append(path(i+n, i+n+pshape(1,k)))
      end do
      n = n + pshape(1,k)*pshape(2,k)
    end do
  end function get_path_list2

  pure type(path_list) function get_path_list3(pshape, hole_pos) result(r)
    integer, dimension(:,:), intent(in) :: pshape
    integer, dimension(:), intent(in) :: hole_pos
    integer                           :: i, z
    call r%alloc_value()
    do i = 1, size(hole_pos)
      z = site_to_z(hole_pos(i), pshape)
      call r%append(path(hole_pos(i)-pshape(1,z),hole_pos(i)+1))
      call r%append(path(hole_pos(i)-1,hole_pos(i)+pshape(1,z)))
    end do
  end function get_path_list3

  pure type(path_list) function get_path_list4(pshape, hole_pos) result(r)
    integer, dimension(:,:), intent(in) :: pshape
    integer, dimension(:), intent(in) :: hole_pos
    integer                           :: i, z
    call r%alloc_value()
    do i = 1, size(hole_pos)
      z = site_to_z(hole_pos(i), pshape)
      call r%append(path(hole_pos(i)-pshape(1,z),hole_pos(i)-1))
      call r%append(path(hole_pos(i)+1,hole_pos(i)+pshape(1,z)))
    end do
  end function get_path_list4

  type(path_list) function get_path_list_xy(pshape) result(r)
    integer, dimension(:,:), intent(in) :: pshape
    integer                           :: i, k, n
    call r%alloc_value()
    n = 0
    do k = 1, size(pshape(1,:))
      do i = 1, pshape(1,k)*(pshape(2,k)-1)
        if(mod(i,pshape(1,k)) == 0) cycle
        call r%append(path(i+n, i+n+1+pshape(1,k)))
      end do
      n = n + pshape(1,k)*pshape(2,k)
    end do
  end function get_path_list_xy

  type(path_list) function get_path_list_mxy(pshape) result(r)
    integer, dimension(:,:), intent(in) :: pshape
    integer                           :: i, k, n
    call r%alloc_value()
    n = 0
    do k = 1, size(pshape(1,:))
      do i = 1, pshape(1,k)*(pshape(2,k)-1)
        if(mod(i-1,pshape(1,k)) == 0) cycle
        call r%append(path(i+n, i+n-1+pshape(1,k)))
      end do
      n = n + pshape(1,k)*pshape(2,k)
    end do
  end function get_path_list_mxy

  type(path_list) function get_path_list_xy1(pshape) result(r)
    integer, dimension(:,:), intent(in) :: pshape
    integer                           :: i, k, n
    call r%alloc_value()
    n = 0
    !do k = 1, size(pshape(1,:))
    k = 1
      do i = 1, pshape(1,k)*(pshape(2,k)-1)
        if(mod(i,pshape(1,k)) == 0) cycle
        call r%append(path(i+n, i+n+1+pshape(1,k)))
      end do
      n = n + pshape(1,k)*pshape(2,k)
    !end do
  end function get_path_list_xy1

  type(path_list) function get_path_list_mxy1(pshape) result(r)
    integer, dimension(:,:), intent(in) :: pshape
    integer                           :: i, k, n
    call r%alloc_value()
    n = 0
    !do k = 1, size(pshape(1,:))
    k = 1
      do i = 1, pshape(1,k)*(pshape(2,k)-1)
        if(mod(i-1,pshape(1,k)) == 0) cycle
        call r%append(path(i+n, i+n-1+pshape(1,k)))
      end do
      n = n + pshape(1,k)*pshape(2,k)
    !!end do
  end function get_path_list_mxy1

  pure type(path_list) function get_path_list_z(pshape) result(r)
  !type(path_list) function get_path_list_z(pshape) result(r)
    integer, dimension(:,:), intent(in) :: pshape
    integer                           :: i, j, n_tmp, n, nx, ny
    integer created
    created = 0
    call r%alloc_value()
    n_tmp = 0
    do i = 1, size(pshape(1,:))-1
      n = pshape(1,i)*pshape(2,i)
      nx = min(pshape(1,i),pshape(1,i+1))
      ny = min(pshape(2,i),pshape(2,i+1))
      do j = 1, nx*ny
        call r%append(path( &
          &    xyz_to_site(mod(j-1,nx)+1,int((j-1)/nx+1),i,pshape), &
          &    xyz_to_site(mod(j-1,nx)+1,int((j-1)/nx+1),i+1,pshape) &
!          &     mod(j,nx) + int((j)/nx)*nx + n_tmp, &
!          &     mod(j,nx) + int((j)/nx)*nx + n_tmp + n &
          &   ))
        created = created + 1
        !print *, j, mod(j-1,nx)+1,int((j-1)/nx)+1,i,  xyz_to_site(mod(j-1,nx)+1,int((j-1)/nx)+1,i,pshape)

        ! test : use only one zpath 
        !return
        !if(created == 2)then
        !  return
        !end if
      end do
      n_tmp = n_tmp + n
    end do
  end function get_path_list_z


  !!  this will delete
  pure type(hop_iterator) function new_hop_iterator_(pshape) result(r)
    integer, dimension(:), intent(in) :: pshape
    call r%hop_iterator_init(reshape(pshape,[2,1]))
  end function new_hop_iterator_

  ! constructor
  pure type(hop_iterator) function new_hop_iterator(pshape) result(r)
    integer, dimension(:,:), intent(in) :: pshape
    call r%hop_iterator_init(pshape)
  end function new_hop_iterator
end module hop_iterator_mod

!-check paths in hop_iterator!!!!!!!!
! open(90,file='path.dat',status = 'replace')
! open(91,file='path.gp',status = 'replace')
! WRITE(91,*) 'reset'
!  WRITE(91,*) 'set nokey'
! ALLOCATE(iter, source = new_hop_iterator(this%pshape))
! WRITE(*,*) 'd1'
! CALL iter%set_path([1,2], this%hole(:, 1))
! count = 1
! DO WHILE(iter%has_next())
!    CALL iter%next(p)
!    x1 = MOD(p%i,this%base%nx)
!    if(x1==0) x1 = this%base%nx
!    y1 = (p%i-x1)/this%base%nx + 1
!    x2 = MOD(p%f,this%base%nx)
!    if(x2==0) x2 = this%base%nx
!    y2 = (p%f-x2)/this%base%nx + 1
!    WRITE(90,'(4I5)') x1,y1,x2-x1,y2-y1
!    WRITE(91,'(A10,I3,A16,f5.2,A1,F5.2,A2,I3,A1)') 'set label ',count, 'center at graph' &
!     &  ,(x2+x1)/(16d0*2d0),',',(y2+y1)/(16d0*2d0) ,'"',count, '"'
!    count = count + 1
! END DO
! WRITE(91,*) 'plot "path.dat" w vec'
! close(90)
! close(91)
! DEALLOCATE(iter)
! !!!!!!!!!!
