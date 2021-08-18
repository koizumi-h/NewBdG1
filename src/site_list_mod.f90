MODULE site_list_mod
  USE path_list_mod, ONLY : path_list,path
  IMPLICIT NONE
  TYPE, PUBLIC :: site_list
     INTEGER               :: n, nx, ny, nh, ne
     !    type(path_list)       :: pl
     INTEGER, DIMENSION(:), ALLOCATABLE  :: hole_pos
     INTEGER, DIMENSION(:,:), ALLOCATABLE :: pos_neighbor, holes_neighbor
   CONTAINS
     ! public type-bound procedures
     PROCEDURE :: site_list_init
     !    procedure :: set_pshape => site_list_set_pshape
     PROCEDURE :: set_neighbor => site_list_set_neighbor
  END TYPE site_list
CONTAINS
  ! public type-bound procedures
  SUBROUTINE site_list_init(this, pshape, hole)
    CLASS(site_list), INTENT(inout)                    :: this
    INTEGER, DIMENSION(2), INTENT(in)                  :: pshape
    INTEGER, DIMENSION(:,:), INTENT(in)                :: hole
    ! call this%set_pshape(pshape)
    this%n = pshape(1) * pshape(2)
    this%nx = pshape(1)
    this%ny = pshape(2)
    this%nh = SIZE(hole, 1)
    this%ne = this%n - this%nh
    IF (this%nh>0) THEN
       ALLOCATE(this%hole_pos(SIZE(hole, 1)))
       this%hole_pos(:) = hole(:,1)
    ELSE
       ALLOCATE(this%hole_pos(1))
       this%hole_pos(1) = 0 !dummy
    ENDIF
    !WRITE(*,*) 'n,nx,ny,nh,ne'
    !WRITE(*,*) this%n,this%nx,this%ny,this%nh,this%ne
    CALL this%set_neighbor()
  END SUBROUTINE site_list_init

  SUBROUTINE site_list_set_neighbor(this)
    CLASS(site_list), INTENT(inout)          :: this
    !    integer, dimension(:), intent(in)           :: hop_dir
    INTEGER                                     :: i, j, k, icount
    !    integer, dimension(:), intent(in), optional     :: hole_pos
    ALLOCATE(this%pos_neighbor(this%ne, 6))
    icount = 1
    DO i = 1, this%n
       k = 3
       DO j = 1, this%nh
          IF(i .EQ. this%hole_pos(j)) GOTO 1000
       ENDDO

       !---current position---
       this%pos_neighbor(icount, 1) = i
       !---ix+1---
       IF(MOD(i, this%nx) .EQ. 0) GOTO 910
       DO j = 1, this%nh
          IF(this%hole_pos(j) .EQ. i+1) GOTO 910
       ENDDO
       this%pos_neighbor(icount, k) = i + 1
       k = k + 1
910    CONTINUE
       !---ix-1---
       IF(MOD(i-1, this%ny) .EQ. 0) GOTO 920
       DO j = 1, this%nh
          IF(this%hole_pos(j) .EQ. i-1) GOTO 920
       ENDDO
       this%pos_neighbor(icount, k) = i - 1 !ix - 1
       k = k + 1
920    CONTINUE
       !---iy+1---
       IF(i+this%nx > this%n) GOTO 930
       DO j = 1, this%nh
          IF(this%hole_pos(j) .EQ. i+this%nx) GOTO 930
       ENDDO
       this%pos_neighbor(icount, k) = i + this%nx !iy + 1
       k = k + 1
930    CONTINUE
       !---iy-1---
       IF(i-this%nx < 1) GOTO 940
       DO j = 1, this%nh
          IF(this%hole_pos(j) .EQ. i-this%nx) GOTO 940
       ENDDO
       this%pos_neighbor(icount, k) = i - this%nx !iy - 1
       k = k + 1
940    CONTINUE
       !---the number of neighbors---
       this%pos_neighbor(icount, 2) = k - 3

       icount = icount +1
1000   CONTINUE
    ENDDO
    !WRITE(*,*) 'the list of site neighbors'
    !DO i = 1, this%ne
    !   WRITE(*,*) this%pos_neighbor(i,1:this%pos_neighbor(i,2)+2)
    !ENDDO
  END SUBROUTINE site_list_set_neighbor


  ! constructor
  FUNCTION new_site_list(pshape, hole) RESULT(r)
    TYPE(site_list)                   :: r
    INTEGER, DIMENSION(2), INTENT(in) :: pshape
    INTEGER, DIMENSION(:,:), INTENT(in)                :: hole
    CALL r%site_list_init(pshape, hole)
  END FUNCTION new_site_list
END MODULE site_list_mod
