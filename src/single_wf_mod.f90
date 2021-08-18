MODULE single_wf_mod
  USE io_mod, ONLY:save_txt, load_txt, create_vortex
  USE base_mod, ONLY:base
  USE math_mod, ONLY:ui, my_zheev, principal, pi
  USE path_list_mod, ONLY:path
  USE hop_iterator_mod, ONLY:hop_iterator, new_hop_iterator
  use parameter_mod
  IMPLICIT NONE
  TYPE, extends(base) :: single_wf
    REAL(8), DIMENSION(:), ALLOCATABLE :: eg
    REAL(8)                            :: epolaron = 100.5d0, tol = 1.d-3
    INTEGER                            :: itermax = 10
    REAL(8), DIMENSION(:), ALLOCATABLE   :: Sx, Sy, Sz, density
  CONTAINS
    ! public type-bound procedures
    PROCEDURE :: single_wf_init
    PROCEDURE :: get_wave_functions => single_wf_get_wave_functions
    PROCEDURE :: get_energies => single_wf_get_energies
    PROCEDURE :: diagonalization => single_wf_diagonalization
    PROCEDURE :: diagonalization_for_newS => single_wf_diagonalization_for_newS
    ! protected type-bound procedures
    PROCEDURE :: calc_S => single_wf_calc_S
    PROCEDURE :: calc_xi => single_wf_calc_xi
    PROCEDURE :: calc_xi_and_zeta => single_wf_calc_xi_and_zeta
    PROCEDURE :: absS => single_wf_absS
    PROCEDURE :: get_density_dependent_ham_d => &
         &              single_wf_get_density_dependent_ham_d
    PROCEDURE :: get_kinetic_ham_d => single_wf_get_kinetic_ham_d
    PROCEDURE :: get_coulomb_ham_d => single_wf_get_coulomb_ham_d
    PROCEDURE :: get_hole_interaction_ham_d => &
         &              single_wf_get_hole_interaction_ham_d
  END TYPE single_wf
CONTAINS
  ! public type-bound procedures
  SUBROUTINE single_wf_init(this, layers, pshape, hole, t, u, jd, lm, xi, chi, wf)
    CLASS(single_wf), INTENT(inout)                            :: this
    INTEGER, DIMENSION(:,:), INTENT(in)                             :: pshape
    INTEGER, DIMENSION(:, :), INTENT(in)                          :: hole
    REAL(8), DIMENSION(:), INTENT(in)                             :: t
    REAL(8), INTENT(in)                                           :: u, jd, lm
    REAL(8), DIMENSION(site_max(pshape)), INTENT(in)           :: xi
    REAL(8), DIMENSION(site_max(pshape)), INTENT(in), OPTIONAL :: chi
    COMPLEX(8), DIMENSION(4*site_max(pshape), 4*site_max(pshape)), &
         & INTENT(in), OPTIONAL                             :: wf
    type(cuo2_plane), dimension(:), allocatable         :: layers
    CALL this%base_init(layers , pshape, hole, t, u, jd, lm, xi, chi, wf)
    IF (ALLOCATED(this%Sx)) DEALLOCATE (this%Sx)
    ALLOCATE (this%Sx(this%n))
    IF (ALLOCATED(this%Sy)) DEALLOCATE (this%Sy)
    ALLOCATE (this%Sy(this%n))
    IF (ALLOCATED(this%Sz)) DEALLOCATE (this%Sz)
    ALLOCATE (this%Sz(this%n))
    IF (ALLOCATED(this%density)) DEALLOCATE (this%density)
    ALLOCATE (this%density(this%n))
  END SUBROUTINE single_wf_init

  FUNCTION single_wf_get_wave_functions(this) RESULT(r)
    CLASS(single_wf), INTENT(in)           :: this
    COMPLEX(8), DIMENSION(4*this%n, 4*this%n) :: r
    r(:, :) = (0.d0, 0.d0)
    IF (ALLOCATED(this%wf)) r(:, :) = this%wf(:, :)
  END FUNCTION single_wf_get_wave_functions

  FUNCTION single_wf_get_energies(this) RESULT(r)
    CLASS(single_wf), INTENT(in) :: this
    REAL(8), DIMENSION(4*this%n)    :: r
    r(:) = 0.d0
    IF (ALLOCATED(this%eg)) r(:) = this%eg(:)
  END FUNCTION single_wf_get_energies

  SUBROUTINE single_wf_diagonalization(this, wf, eg)
    CLASS(single_wf), INTENT(inout)        :: this
    COMPLEX(8), DIMENSION(4*this%n, 4*this%n), &
         & INTENT(inout)                           :: wf
    COMPLEX(8), DIMENSION(4*this%n, 4*this%n) :: ham, kh, dh
    REAL(8), DIMENSION(4*this%n), INTENT(out) :: eg
    REAL(8), DIMENSION(this%n)                :: xi
    CALL this%calc_S(wf)
    xi = this%calc_xi()
    this%xi = xi
    wf(:, :) = (0.d0, 0.d0)
    ham(:, :) = (0.d0, 0.d0)
    kh = this%get_kinetic_ham_d()
    dh = this%get_density_dependent_ham_d()
    ham = kh + dh
    !ham = this%get_kinetic_ham_d() + this%get_density_dependent_ham_d()
    !CALL my_zheev('l', ham, eg, wf)
    !WRITE(6,*)"diagonal etot", SUM(eg(1:sum(int(this%ne))))
  END SUBROUTINE single_wf_diagonalization

  SUBROUTINE single_wf_diagonalization_for_newS(this, wf, eg, Sx, Sy, Sz)
    CLASS(single_wf), INTENT(inout)        :: this
    COMPLEX(8), DIMENSION(4*this%n, 4*this%n), &
         & INTENT(inout)                           :: wf
    COMPLEX(8), DIMENSION(4*this%n, 4*this%n) :: ham
    REAL(8), DIMENSION(4*this%n), INTENT(out) :: eg
    REAL(8), DIMENSION(this%n)                :: xi, Sx, Sy, Sz
    this%Sx(:) = Sx(:)
    this%Sy(:) = Sy(:)
    this%Sz(:) = Sz(:)
    xi = this%calc_xi()
    this%xi = xi
    wf(:, :) = (0.d0, 0.d0)
    ham(:, :) = (0.d0, 0.d0)
    ham = this%get_kinetic_ham_d() + this%get_density_dependent_ham_d()
    !WRITE(*,*) 'next is my_zheev'
    !CALL my_zheev('l', ham, eg, wf)
    !WRITE (6, *) "diagonal etot", SUM(eg(1:sum(int(this%ne))))
  END SUBROUTINE single_wf_diagonalization_for_newS

! protected type-bound procedures
  SUBROUTINE single_wf_calc_S(this, wf)
    CLASS(single_wf), INTENT(inout)                         :: this
    COMPLEX(8), DIMENSION(4*this%n, 4*this%n), INTENT(in) :: wf
    COMPLEX(8)                                            :: Du, Dd, Sp
    REAL(8)                                               :: Duu, Ddd
    INTEGER                                               :: i, j
    DO j = 1, this%n
      this%Sx(j) = 0.d0
      this%Sy(j) = 0.d0
      this%Sz(j) = 0.d0
      this%density(j) = 0.d0
      Sp = 0.d0
      DO i = 1, sum(int(this%ne))
        Du = wf(2*j - 1, i)
        Dd = wf(2*j, i)
        Duu = REAL(CONJG(Du)*Du, KIND(0d0))
        Ddd = REAL(CONJG(Dd)*Dd, KIND(0d0))
        Sp = Sp + CONJG(Du)*Dd
        this%Sz(j) = this%Sz(j) + 0.5d0*(Duu - Ddd)
        this%density(j) = this%density(j) + Duu + Ddd
      END DO
      this%Sx(j) = REAL(Sp, KIND(0d0))
      this%Sy(j) = AIMAG(Sp)
    END DO
  END SUBROUTINE single_wf_calc_S

  ELEMENTAL FUNCTION single_wf_absS(this, i) RESULT(r)
    CLASS(single_wf), INTENT(in)       :: this
    INTEGER, INTENT(in)              :: i
    REAL(8)                          :: r
    r = ABS(this%Sx(i) + ui*this%Sy(i))
  END FUNCTION single_wf_absS

  FUNCTION single_wf_calc_xi(this) RESULT(r)
    CLASS(single_wf), INTENT(in)       :: this
    REAL(8), DIMENSION(this%n)       :: r
    INTEGER                          :: j
    r(:) = this%Sx(:)
    DO j = 1, this%n
      IF (r(j) == 0.d0) r(j) = 1.d-10
      r(j) = -0.5d0*(SIGN(1.0d0, r(j)) - 1.0d0)*pi + ATAN(this%Sy(j)/r(j))
    END DO
  END FUNCTION single_wf_calc_xi

  SUBROUTINE single_wf_calc_xi_and_zeta(this, xi, zeta)
    CLASS(single_wf), INTENT(in)                  :: this
    REAL(8), DIMENSION(this%n), INTENT(out)       :: xi, zeta
    REAL(8)                                       :: norm
    INTEGER                                       :: j
    xi(:) = ATAN2(this%Sy(:), this%Sx(:))
    DO j = 1, this%n
      norm = SQRT(this%Sx(j)**2 + this%Sy(j)**2)
      zeta(j) = ATAN2(norm, this%Sz(j))
    ENDDO
  END SUBROUTINE single_wf_calc_xi_and_zeta

  FUNCTION single_wf_get_kinetic_ham_d(this) RESULT(r)
    CLASS(single_wf), INTENT(in)           :: this
    COMPLEX(8), DIMENSION(4*this%n, 4*this%n) :: r
    TYPE(hop_iterator), ALLOCATABLE           :: iter
    TYPE(path)                                :: p
    COMPLEX(8)                                :: temp
    r(:, :) = (0.d0, 0.d0)
    ALLOCATE (iter, source=new_hop_iterator(this%pshape))
    CALL iter%set_path([1, 2], this%hole(:, 1))
    DO WHILE (iter%has_next())
      CALL iter%next(p)
      temp = EXP(0.5d0*ui*MODULO(this%xi(p%f) - this%xi(p%i), 2.d0*pi))
      r(2*p%f - 1, 2*p%i - 1) = -this%t(1)*temp
      r(2*p%i - 1, 2*p%f - 1) = -this%t(1)*CONJG(temp)
      r(2*p%f, 2*p%i) = -this%t(1)*CONJG(temp)
      r(2*p%i, 2*p%f) = -this%t(1)*temp
    END DO
    DEALLOCATE (iter)
    !ALLOCATE (iter, source=new_hop_iterator(this%pshape))
    !CALL iter%set_path([5, 6], this%hole(:, 1))
    !DO WHILE (iter%has_next())
    !  CALL iter%next(p)
    !  temp = EXP(0.5d0*ui*MODULO(this%xi(p%f) - this%xi(p%i), 2.d0*pi))
    !  r(2*p%f - 1, 2*p%i - 1) = -this%t(2)*temp
    !  r(2*p%i - 1, 2*p%f - 1) = -this%t(2)*CONJG(temp)
    !  r(2*p%f, 2*p%i) = -this%t(2)*CONJG(temp)
    !  r(2*p%i, 2*p%f) = -this%t(2)*temp
    !END DO
    !DEALLOCATE (iter)
    ALLOCATE (iter, source=new_hop_iterator(this%pshape))
    CALL iter%set_path([7], this%hole(:, 1))
    DO WHILE (iter%has_next())
      CALL iter%next(p)
      temp = EXP(0.5d0*ui*MODULO(this%xi(p%f) - this%xi(p%i), 2.d0*pi))
      r(2*p%f - 1, 2*p%i - 1) = -this%t(3)*temp
      r(2*p%i - 1, 2*p%f - 1) = -this%t(3)*CONJG(temp)
      r(2*p%f, 2*p%i) = -this%t(3)*CONJG(temp)
      r(2*p%i, 2*p%f) = -this%t(3)*temp
    END DO
    DEALLOCATE (iter)
  END FUNCTION single_wf_get_kinetic_ham_d

  FUNCTION single_wf_get_density_dependent_ham_d(this) RESULT(r)
    CLASS(single_wf), INTENT(in)           :: this
    COMPLEX(8), DIMENSION(4*this%n, 4*this%n) :: r
    r(:, :) = this%get_coulomb_ham_d()
    r(:, :) = r(:,:) + this%get_hole_interaction_ham_d()
  END FUNCTION single_wf_get_density_dependent_ham_d

  FUNCTION single_wf_get_coulomb_ham_d(this) RESULT(r)
    CLASS(single_wf), INTENT(in)           :: this
    COMPLEX(8), DIMENSION(4*this%n, 4*this%n) :: r
    INTEGER                                   :: j
    r(:, :) = (0.d0, 0.d0)
    FORALL (j=1:this%n)
    r(2*j - 1, 2*j - 1) = this%u*(0.5d0*this%density(j) - this%Sz(j))
    r(2*j, 2*j - 1) = -this%u*this%absS(j)
    r(2*j - 1, 2*j) = -this%u*this%absS(j)
    r(2*j, 2*j) = this%u*(0.5d0*this%density(j) + this%Sz(j))
    END FORALL
    FORALL (j=1:SIZE(this%hole, 1))
    r(2*this%hole(j, 1) - 1, 2*this%hole(j, 1) - 1) &
         & = r(2*this%hole(j, 1) - 1, 2*this%hole(j, 1) - 1) + this%epolaron
    r(2*this%hole(j, 1), 2*this%hole(j, 1)) &
         & = r(2*this%hole(j, 1), 2*this%hole(j, 1)) + this%epolaron
    END FORALL
  END FUNCTION single_wf_get_coulomb_ham_d

  FUNCTION single_wf_get_hole_interaction_ham_d(this) RESULT(r)
    CLASS(single_wf), INTENT(in)           :: this
    TYPE(hop_iterator), ALLOCATABLE           :: iter
    TYPE(path)                                :: p
    COMPLEX(8), DIMENSION(4*this%n, 4*this%n) :: r
    COMPLEX(8)                                :: temp
    r(:, :) = (0.d0, 0.d0)
    ALLOCATE (iter, source=new_hop_iterator(this%pshape))
    CALL iter%set_path_hole(this%hole(:, 1))
    DO WHILE (iter%has_next())
      CALL iter%next(p)
      temp = EXP(0.5d0*ui*MODULO(this%xi(p%f) - this%xi(p%i), 2.d0*pi))
      r(2*p%f - 1, 2*p%f) = this%jd*this%absS(p%i)*temp
      r(2*p%i - 1, 2*p%i) = this%jd*this%absS(p%f)*CONJG(temp)
      r(2*p%f, 2*p%f - 1) = this%jd*this%absS(p%i)*CONJG(temp)
      r(2*p%i, 2*p%i - 1) = this%jd*this%absS(p%f)*temp
      r(2*p%f - 1, 2*p%f - 1) = this%jd*this%Sz(p%i)
      r(2*p%i - 1, 2*p%i - 1) = this%jd*this%Sz(p%f)
      r(2*p%f, 2*p%f) = -this%jd*this%Sz(p%i)
      r(2*p%i, 2*p%i) = -this%jd*this%Sz(p%f)
    END DO
    DEALLOCATE (iter)
  END FUNCTION single_wf_get_hole_interaction_ham_d

! constructor
  FUNCTION new_single_wf(pshape, hole, t, u, jd, lm, xi, chi, wf) RESULT(r)
    INTEGER, DIMENSION(:,:), INTENT(in)                   :: pshape
    INTEGER, DIMENSION(:, :), INTENT(in)                :: hole
    REAL(8), DIMENSION(:), INTENT(in)                   :: t
    REAL(8), INTENT(in)                                 :: u, jd, lm
    REAL(8), DIMENSION(site_max(pshape)), INTENT(in) :: xi
    REAL(8), DIMENSION(site_max(pshape)), INTENT(in), OPTIONAL :: chi
    COMPLEX(8), DIMENSION(4*site_max(pshape), 4*site_max(pshape)), &
         & INTENT(in), OPTIONAL                             :: wf
    TYPE(single_wf)                                  :: r
    type(cuo2_plane), dimension(:), allocatable         :: layers
    CALL r%single_wf_init(layers , pshape, hole, t, u, jd, lm, xi, chi, wf)
  END FUNCTION new_single_wf
END MODULE single_wf_mod
