MODULE energy_mod
  USE math_mod, ONLY:pi, ui, principal, sign2, principal2, principal
  use path_list_mod
  use hop_iterator_mod
  use parameter_mod
  use base_mod, only : base
  IMPLICIT NONE
  !  type, extends(monte_carlo),public :: energy
  !TYPE, PUBLIC :: energy
  !    type(dchi_optimizer)                      :: opt
  ! CONTAINS
  ! public type-bound procedures
  !   PROCEDURE :: set_wf           => energy_set_wf
  !   PROCEDURE :: calc_S           => energy_calc_S
  !   PROCEDURE :: total_energy           => energy_total_energy
  !END TYPE energy
CONTAINS

 ! pure function energy_set_wf(n, hole, cpwf, xi, chi) result(wf)
 !   integer, intent(in)                 :: n, hole(:, :)
 !   REAL(8), INTENT(in)                 :: xi(n), chi(n)
 !   COMPLEX(8), INTENT(in)              :: cpwf(4*n, 4*n)
 !   complex(8)                             :: wf(4*n, 4*n)
 !   INTEGER                             :: i, j,ne
 !   COMPLEX(8)                          :: exi, echi
 !   wf = cpwf
 !   ne = n - size(hole(:,1))
 !   DO i = 1, n
 !     exi = EXP(-ui*0.5d0*xi(i))
 !     echi = EXP(-ui*0.5d0*chi(i))
 !     DO j = 1, 2*ne
 !       wf(4*i - 3, j + 2*ne) = wf(4*i - 3, j + 2*ne)*(exi)*echi
 !       wf(4*i - 2, j + 2*ne) = wf(4*i - 2, j + 2*ne)*CONJG(exi)*echi
 !       wf(4*i - 1, j) = wf(4*i - 1, j)*CONJG(exi)*CONJG(echi)
 !       wf(4*i, j) = wf(4*i, j)*CONJG(echi)*exi
 !     END DO
 !   END DO
 !   DO i = 1, SIZE(hole(:, 1))
 !     wf(4*hole(i, 1) - 1, :) = 0.d0
 !     wf(4*hole(i, 1), :) = 0.d0
 !     wf(4*hole(i, 1) - 1, :) = 0.d0
 !     wf(4*hole(i, 1), :) = 0.d0
 !   ENDDO
 ! END function energy_set_wf

 ! !subroutine energy_calc_S(xi, chi, Sx, Sy, Sz)
 ! pure subroutine energy_calc_S(n, ne, wf, Sx, Sy, Sz)
 !   !real(8), intent(in)         :: xi(:), chi(:)
 !   integer, intent(in)         :: n, ne
 !   complex(8), intent(in)      :: wf(4*n, 4*n) ! both wf, cpwf are okay
 !   real(8), intent(out)        :: Sx(n), Sy(n), Sz(n)
 !   real(8)                     :: r_xi, r_sx, r_sy, r_sz
 !   complex(8)                  :: Sp, uu, ud, vu, vd, density(n)
 !   integer                     :: i, j
 !   Sx(:) = 0.d0
 !   Sy(:) = 0.d0
 !   Sz(:) = 0.d0
 !   density(:) = 0.d0
 !   !do i = 1, n
 !   DO j = 1, n
 !     DO i = 2*ne + 1, 4*ne
 !       uu = wf(4*j - 3, i)
 !       ud = wf(4*j - 2, i)
 !       vu = wf(4*j - 1, i)
 !       vd = wf(4*j, i)
 !       Sx(j) = Sx(j) + real(vu*conjg(vd))
 !       Sy(j) = Sy(j) - aimag(vu*conjg(vd))
 !       Sz(j) = Sz(j) + real(vu*conjg(vu)-vd*conjg(vd))/2
 !       density(j) = density(j)  + real (vu*conjg(vu)+vd*conjg(vd))
 !     END DO
 !   END DO
 !   !write(*, *) 'density', sum(density(1:n))
 !   !do k = 1, n
 !     !write(*, *)  k, density(k)
 !   !enddo
 !   !DO j = 1, n
 !       !r_xi =  ATAN2(Sy(j), Sx(j))
 !     !write(*,*) 'r_xi', j, r_xi
 !   !END DO
 ! end subroutine energy_calc_S

 ! FUNCTION energy_total_energy(pshape, hole, t, u, jd, lm, &
 ! !& xi, chi, cpwf, rs, view) RESULT(total_E)
 ! & xi, chi, cpwf, mu, ne_u, ne_d, delta, rs, view) RESULT(total_E)
 !   integer, intent(in)             :: hole(:, :), pshape(:,:)
 !   REAL(8), INTENT(in)             :: xi(:), chi(size(xi, 1)), &
 !     &                              u, jd, lm, t(3)
 !   COMPLEX(8), INTENT(in)          :: cpwf(:, :)
 !   real(8),intent(in)    :: mu(2)
 !   complex(8), dimension(:), intent(in) :: delta
 !   real(8), dimension(:), intent(in)   :: ne_u, ne_d
 !   !real(8), dimension(:), intent(in) :: Sx, Sy, Sz
 !   real(8)                         :: total_E
 !   COMPLEX(8)                      :: wf1(size(cpwf, 1), size(cpwf, 2)), &
 !     !&                                Hk, Hk_2, Hk_z, Hu, Hh, Hrsoc, Hhz
 !     &                                Hk, Hk_z, Hu, Hh, Hrsoc, Hhz
 !   !REAL(8)                         :: Sx(size(xi, 1)), Sy(size(xi, 1)), Sz(size(xi, 1)), &
 !     !&                                eg(size(xi, 1))
 !   !INTEGER                         :: i, j, k, l, a, b, &
 !   integer                         :: nx, ne, n
 !   character(1), optional          :: rs
 !   logical,intent(in), optional               :: view
 !   logical               :: f_view
 !   !type(path), allocatable         :: p(:)
 !   !type(path)                      :: p1
 !   !type(hop_iterator), allocatable :: iter
 !   f_view = .true.
 !   if(present(view)) then
 !     f_view = view
 !   end if
 !   
 !   n = site_max(pshape)
 !   ne = n - size(hole, 1)
 !   wf1 = energy_set_wf(n, hole, cpwf, xi, chi)
 !   Hk = energy_hopping_energy(pshape, n, ne, hole, t, &
 !   & wf1) 
 !   !Hk_2 = energy_hopping_2nd_energy(pshape, n, ne, hole, t, &
 !   !& wf1) 
 !   Hk_z = energy_hopping_z_energy(pshape, n, ne, hole, t, &
 !   & wf1) 
 !   Hu = energy_coulomb_interaction_energy(pshape, n, ne, hole, u, &
 !   & wf1) 
 !   Hh = energy_across_hole_interaction_energy(pshape, n, ne, hole, jd, &
 !   & wf1) 
 !   if(present(rs)) then
 !     Hrsoc = energy_spin_orbit_interaction_energy(pshape, n, ne, hole, lm, &
 !     & wf1)
 !   endif
 !   !total_E = real(Hk + Hk_2 + Hk_z + Hu + Hh + Hrsoc, 8)
 !   total_E = real(Hk + Hk_z + Hu + Hh + Hrsoc, 8)
 !   if (f_view) then
 !     write(*,*) 'Hu',  Hu
 !     write(*,*) 'Hk',  Hk
 !     !write(*,*) 'Hk_2nd',  Hk_2
 !     write(*,*) 'Hk_z',  Hk_z
 !     write(*,*) 'Hh',  Hh
 !     write(*,*) 'Hrsoc',  Hrsoc
 !   end if
 ! end function energy_total_energy

 ! !FUNCTION energy_total_energy2(pshape, hole, t, u, jd, lm, &
 ! !& xi, chi, cpwf, EV, rs) RESULT(total_E)
 ! !  integer, intent(in)             :: hole(:, :), pshape(2)
 ! !  REAL(8), INTENT(in)             :: xi(:), chi(size(xi, 1)), &
 ! !    &                              u, jd, lm, t(2)
 ! !  COMPLEX(8), INTENT(in)          :: cpwf(:, :)
 ! !  real(8)                         :: total_E
 ! !  COMPLEX(8)                      :: wf1(size(cpwf, 1), size(cpwf, 2)), &
 ! !    &                                Hk, Hk_2, Hu, Hh, Hrsoc, Hhz, EV
 ! !  !REAL(8)                         :: Sx(size(xi, 1)), Sy(size(xi, 1)), Sz(size(xi, 1)), &
 ! !    !&                                eg(size(xi, 1))
 ! !  !INTEGER                         :: i, j, k, l, a, b, &
 ! !  integer                         :: nx, ne, n
 ! !  character(1), optional          :: rs
 ! !  !type(path), allocatable         :: p(:)
 ! !  !type(path)                      :: p1
 ! !  !type(hop_iterator), allocatable :: iter
 ! !  n = pshape(1) * pshape(2)
 ! !  ne = n - size(hole, 1)
 ! !  wf1 = energy_set_wf(n, hole, cpwf, xi, chi)
 ! !  Hk = energy_hopping_energy(pshape, n, ne, hole, t, &
 ! !  & wf1) 
 ! !  !Hk_2 = energy_hopping_2nd_energy(pshape, n, ne, hole, t, &
 ! !  !& wf1) 
 ! !  Hu = energy_coulomb_interaction_energy(pshape, n, ne, hole, u, &
 ! !  & wf1) 
 ! !  Hh = energy_across_hole_interaction_energy(pshape, n, ne, hole, jd, &
 ! !  & wf1) 
 ! !  if(present(rs)) then
 ! !    Hrsoc = energy_spin_orbit_interaction_energy2(pshape, n, ne, hole, lm, &
 ! !    &xi, chi, cpwf, EV)
 ! !  endif
 ! !  total_E = real(Hk + Hk_2+ Hu + Hh + Hrsoc, 8)
 ! !  !write(*,*) 'Hu',  Hu
 ! !  !write(*,*) 'Hk',  Hk
 ! !  !write(*,*) 'Hh',  Hh
 ! !  !write(*,*) 'Hrsoc',  Hrsoc
 ! !end function energy_total_energy2

 ! function energy_hopping_energy(pshape, n, ne, hole, t, &
 !   & wf1) result(Hk)
 !   integer, intent(in)             :: hole(:, :), pshape(:,:), n, ne
 !   REAL(8), INTENT(in)             :: t(2) 
 !   COMPLEX(8), INTENT(in)          :: wf1(:, :) ! wf which set xi and chi are set in
 !   real(8)                         :: total_E
 !   complex(8)                      :: Hk
 !   REAL(8)                         :: eg(ne*4)
 !   INTEGER                         :: j, k, l, a, b, &
 !     &                                nx
 !   type(path)                      :: p1
 !   type(hop_iterator), allocatable :: iter
 !   eg = 0.d0
 !   Hk = 0.d0
 !   !nx = pshape(1)
 !   !expectation value of K
 !   ALLOCATE (iter, source=new_hop_iterator(pshape))
 !   CALL iter%set_path([1, 2], hole(:, 1))
 !   DO WHILE (iter%has_next())
 !     CALL iter%next(p1)
 !     k = p1%i
 !     j = p1%f
 !     !exi  = exp(ui*0.5d0*(chi(k)-chi(j)))
 !     !echi = exp(ui*0.5d0*(xi(k)-xi(j)))
 !     DO a = 2*ne + 1, 4*ne
 !       eg(a) = eg(a) + (- CONJG(wf1(4*k - 3, a)) * wf1(4*j - 3, a) &
 !            &           - CONJG(wf1(4*k - 2, a)) * wf1(4*j - 2, a) &
 !            &           + CONJG(wf1(4*k - 1, a)) * wf1(4*j - 1, a) &
 !            &           + CONJG(wf1(4*k, a))     * wf1(4*j, a) &
 !            &           - CONJG(wf1(4*j - 3, a)) * wf1(4*k - 3, a) &
 !            &           - CONJG(wf1(4*j - 2, a)) * wf1(4*k - 2, a) & 
 !            &           + CONJG(wf1(4*j - 1, a)) * wf1(4*k - 1, a) &
 !            &           + CONJG(wf1(4*j, a))     * wf1(4*k, a)           )
 !     END DO
 !   END DO
 !   DEALLOCATE (iter)
 !   Hk = -t(1)*SUM(eg(2*n+1:4*ne))
 !   !  write(*,*)'K energy',Hk
 ! end function energy_hopping_energy

 ! !pure function energy_hopping_2nd_energy(pshape, n, ne, hole, t, &
 ! !  & wf1) result(Hk_2nd)
 ! !  integer, intent(in)             :: hole(:, :), pshape(:,:), n, ne
 ! !  REAL(8), INTENT(in)             :: t(2) 
 ! !  COMPLEX(8), INTENT(in)          :: wf1(:, :) ! wf which set xi and chi are set in
 ! !  complex(8)                      :: Hk_2nd
 ! !  real(8)                         :: total_E
 ! !  REAL(8)                         :: eg(ne*2)
 ! !  INTEGER                         :: j, k, l, a, b, &
 ! !    &                                nx
 ! !  type(path)                      :: p1
 ! !  type(hop_iterator), allocatable :: iter
 ! !  eg = 0.d0
 ! !  Hk_2nd = 0.d0
 ! !  !nx = pshape(1)
 ! !  !expectation value of K
 ! !  ALLOCATE(iter, source=new_hop_iterator(pshape))
 ! !  CALL iter%set_path([5, 6], hole(:, 1))
 ! !  DO WHILE (iter%has_next())
 ! !    CALL iter%next(p1)
 ! !    k = p1%i
 ! !    j = p1%f
 ! !    !exi  = exp(ui*0.5d0*(chi(k)-chi(j)))
 ! !    !echi = exp(ui*0.5d0*(xi(k)-xi(j)))
 ! !    DO a = 1, ne
 ! !      eg(a) = eg(a) + (  CONJG(wf1(2*k - 1, a)) * wf1(2*j - 1, a) &
 ! !           &           + CONJG(wf1(2*k, a))     * wf1(2*j, a) &
 ! !           &           + CONJG(wf1(2*j - 1, a)) * wf1(2*k - 1, a) &
 ! !           &           + CONJG(wf1(2*j, a))     * wf1(2*k, a)           )
 ! !    END DO
 ! !  END DO
 ! !  DEALLOCATE(iter)
 ! !  Hk_2nd = -t(2)*SUM(eg(1:ne))
 ! !  !  write(*,*)'K energy',Hk
 ! !end function energy_hopping_2nd_energy

 ! function energy_hopping_z_energy(pshape, n, ne, hole, t, &
 !   & wf1) result(Hk_z)
 !   integer, intent(in)             :: hole(:, :), pshape(:,:), n, ne
 !   REAL(8), INTENT(in)             :: t(3) 
 !   COMPLEX(8), INTENT(in)          :: wf1(:, :) ! wf which set xi and chi are set in
 !   real(8)                         :: total_E
 !   complex(8)                      :: Hk_z
 !   REAL(8)                         :: eg(ne*4)
 !   INTEGER                         :: j, k, l, a, b, &
 !     &                                nx
 !   type(path)                      :: p1
 !   type(hop_iterator), allocatable :: iter
 !   eg = 0.d0
 !   Hk_z = 0.d0
 !   !nx = pshape(1)
 !   !expectation value of K
 !   ALLOCATE(iter, source=new_hop_iterator(pshape))
 !   CALL iter%set_path([7], hole(:, 1))
 !   DO WHILE (iter%has_next())
 !     CALL iter%next(p1)
 !     k = p1%i
 !     j = p1%f
 !     !exi  = exp(ui*0.5d0*(chi(k)-chi(j)))
 !     !echi = exp(ui*0.5d0*(xi(k)-xi(j)))
 !     DO a = 2*ne + 1, 4*ne
 !       eg(a) = eg(a) + (  CONJG(wf1(4*k - 3, a)) * wf1(4*j - 3, a) &
 !            &           + CONJG(wf1(4*k - 2, a)) * wf1(4*j - 2, a) &
 !            &           - CONJG(wf1(4*k - 1, a)) * wf1(4*j - 1, a) &
 !            &           - CONJG(wf1(4*k, a))     * wf1(4*j, a) &
 !            &           + CONJG(wf1(4*j - 3, a)) * wf1(4*k - 3, a) &
 !            &           + CONJG(wf1(4*j - 2, a)) * wf1(4*k - 2, a) & 
 !            &           - CONJG(wf1(4*j - 1, a)) * wf1(4*k - 1, a) &
 !            &           - CONJG(wf1(4*j, a))     * wf1(4*k, a)          )
 !     END DO
 !   END DO
 !   DEALLOCATE(iter)
 !   Hk_z = -t(3)*SUM(eg(2*ne+1:4*ne))
 !   !  write(*,*)'K energy',Hk
 ! end function energy_hopping_z_energy

 ! function energy_coulomb_interaction_energy(pshape, n, ne, hole, u, &
 !   & wf1) result(Hu)
 !   integer, intent(in)             :: hole(:, :), pshape(:,:), n, ne
 !   REAL(8), INTENT(in)             :: u 
 !   COMPLEX(8), INTENT(in)          :: wf1(:, :) ! wf which set xi and chi are set in: wf1=cpwf*exi*echi
 !   real(8)                         :: total_E
 !   complex(8)                      :: Hu
 !   REAL(8)                         :: Sx(n), Sy(n), Sz(n)
 !   INTEGER                         :: j, l, a, b
 !   Hu = 0.d0
 !   call energy_calc_S(n, ne, wf1, Sx, Sy, Sz)
 !   !expectation value of Hu
 !   !  teg=0.d0
 !   DO a = 2*ne + 1, 4*ne
 !   !DO a = 1, n
 !     DO j = 1, n
 !       Hu = Hu +&
 !            ! ??miss?? -> !!!&u*(- ((2d0/3d0)*Sz(j) + 0.5d0)     * CONJG(wf1(2*j - 1, a))*wf1(2*j - 1, a) &
 !            &u*(  &
 !            & (-(2d0/3d0)*Sz(j) + 0.5d0)  &
 !            &  * CONJG(wf1(4*j - 3, a))*wf1(4*j - 3, a) &
 !            &   + ((2d0/3d0)*Sz(j) + 0.5d0)   * CONJG(wf1(4*j - 2, a))*wf1(4*j - 2, a) &
 !            &   + ((2d0/3d0)*Sz(j) + 0.5d0)    * CONJG(wf1(4*j - 1, a))*wf1(4*j - 1, a) &
 !            &   - ((2d0/3d0)*Sz(j) + 0.5d0)   * CONJG(wf1(4*j, a))*wf1(4*j, a) &
 !            &   - (2d0/3d0)*(Sx(j) - ui*Sy(j))* CONJG(wf1(4*j - 3, a))*wf1(4*j - 2, a) &
 !            &   - (2d0/3d0)*(Sx(j) + ui*Sy(j))* (CONJG(wf1(4*j -2, a))*wf1(4*j - 3, a)) &
 !            &   + (2d0/3d0)*(Sx(j) - ui*Sy(j))* CONJG(wf1(4*j - 1, a))*wf1(4*j, a) &
 !            &   + (2d0/3d0)*(Sx(j) + ui*Sy(j))* (CONJG(wf1(4*j, a))*wf1(4*j - 1, a)) )
 !     ENDDO
 !   ENDDO
 !   DO j = 1, n
 !     Hu = Hu + (2d0*u/3d0)*(Sx(j)**2 + Sy(j)**2 + Sz(j)**2)
 !   ENDDO
 ! end function energy_coulomb_interaction_energy

 ! pure function energy_across_hole_interaction_energy(pshape, n, ne, hole, jd, &
 !   & wf1) result(Hh)
 !   integer, intent(in)             :: hole(:, :), pshape(:,:), n, ne
 !   REAL(8), INTENT(in)             :: jd 
 !   COMPLEX(8), INTENT(in)          :: wf1(:, :) ! wf which set xi and chi are set in: wf1=cpwf*exi*echi
 !   real(8)                         :: total_E
 !   COMPLEX(8)                      :: Hh, Hhz 
 !   REAL(8)                         :: Sx(n), Sy(n), Sz(n), &
 !     &                                eg(n)
 !   INTEGER                         :: i, j, k, l, a, b
 !   type(hop_iterator), allocatable :: iter
 !   Hh = 0.d0
 !   call energy_calc_S(n, ne, wf1, Sx, Sy, Sz)
 !   !expectation value of Hh
 !   !  teg=0.d0
 !   ! jd= 0.5d0*4.0d0*(t(1)**2)/u   !jd=0.5j
 !   ALLOCATE (iter, source=new_hop_iterator(pshape))
 !   CALL iter%set_path_hole(hole(:, 1))
 !   DO a = 2*ne + 1, 4*ne
 !     DO k = 1, SIZE(iter%pl%value(:))
 !       i = iter%pl%value(k)%i
 !       j = iter%pl%value(k)%f
 !       Hh = Hh +&
 !            &0.5d0*jd*(-Sz(j)*CONJG(wf1(4*i, a))*wf1(4*i, a) &
 !            &                    + Sz(j)*CONJG(wf1(4*i - 1, a))*wf1(4*i - 1, a) &
 !            &                    + (Sx(j) - ui*Sy(j))*CONJG(wf1(4*i - 1, a))*wf1(4*i, a) &
 !            &                    + (Sx(j) + ui*Sy(j))*CONJG(wf1(4*i, a))*wf1(4*i - 1, a) &
 !            &                    + Sz(i)*CONJG(wf1(4*j, a))*wf1(4*j, a) &
 !            &                    - Sz(i)*CONJG(wf1(4*j - 1, a))*wf1(4*j - 1, a) &
 !            &                    - (Sx(i) - ui*Sy(i))*CONJG(wf1(4*j - 1, a))*wf1(4*j, a) &
 !            &                    - (Sx(i) + ui*Sy(i))*CONJG(wf1(4*j, a))*wf1(4*j - 1, a))
 !     ENDDO
 !   ENDDO
 !   DO k = 1, SIZE(iter%pl%value(:))
 !     i = iter%pl%value(k)%i
 !     j = iter%pl%value(k)%f
 !     Hh = Hh - jd*(Sx(i)*Sx(j) + Sy(i)*Sy(j) + Sz(i)*Sz(j)) !constant
 !   ENDDO
 !   DEALLOCATE (iter)
 !   !  write(*,*)'Hz energy',Hhz
 !   !  write(*,*)'Hxy energy',Hh
 !   !  write(*,*)'Hh energy',Hh
 ! end function energy_across_hole_interaction_energy

 ! function energy_spin_orbit_interaction_energy(pshape, n, ne, hole, lm, &
 !   & wf1) result(Hrsoc)
 !   integer, intent(in)             :: hole(:, :), pshape(:,:), n, ne
 !   REAL(8), INTENT(in)             :: lm 
 !   COMPLEX(8), INTENT(in)          :: wf1(:, :) ! wf which set xi and chi are set in: wf1=cpwf*exi*echi
 !   real(8)                         :: total_E
 !   COMPLEX(8)                      :: Hrsoc, temp
 !   REAL(8)                         :: eg(n)
 !   INTEGER                         :: i, j, k, l, a, b, z, nx
 !   type(path), allocatable         :: p(:)
 !   Hrsoc = 0.d0
 !   IF (ALLOCATED(p)) DEALLOCATE (p)
 !   ALLOCATE (p(4))
 !   !expectation value of Hrsoc
 !   do z=1, size(pshape(1,:)) ! for nz
 !     nx = pshape(1,z)
 !     DO i = 1, SIZE(hole(:, 1))
 !       if(site_to_z(hole(i,1), pshape) /= z) cycle
 !       p(1)%i = hole(i, 1) - nx
 !       p(1)%f = hole(i, 1) + 1
 !       p(2)%i = hole(i, 1) - 1
 !       p(2)%f = hole(i, 1) + nx
 !       p(3)%i = hole(i, 1) - nx
 !       p(3)%f = hole(i, 1) - 1
 !       p(4)%i = hole(i, 1) + 1
 !       p(4)%f = hole(i, 1) + nx
 !       DO a = 2*ne + 1, 4*ne
 !         DO l = 1, 2
 !           k = p(l)%i
 !           j = p(l)%f
 !           Hrsoc = Hrsoc + CONJG(wf1(4*j - 2, a))*wf1(4*k - 3, a) - CONJG(wf1(4*j - 3, a))*wf1(4*k - 2, a)&
 !                 &       + CONJG(wf1(4*k - 3, a))*wf1(4*j - 2, a) - CONJG(wf1(4*k - 2, a))*wf1(4*j - 3, a) &
 !                 &       - CONJG(wf1(4*j, a))*wf1(4*k - 1, a) - CONJG(wf1(4*j - 1, a))*wf1(4*k, a)&
 !                 &       - CONJG(wf1(4*k - 1, a))*wf1(4*j, a) - CONJG(wf1(4*k, a))*wf1(4*j - 1, a)
 !         ENDDO
 !         DO l = 3, 4
 !           k = p(l)%i
 !           j = p(l)%f
 !           Hrsoc = Hrsoc + ui*(CONJG(wf1(4*j - 2, a))*wf1(4*k - 3, a) + CONJG(wf1(4*j - 3, a))*wf1(4*k - 2, a))&
 !                 &- ui*(CONJG(wf1(4*k - 3, a))*wf1(4*j - 2, a) + CONJG(wf1(4*k - 2, a))*wf1(4*j - 3, a))&
 !                 &- ui*(CONJG(wf1(4*j, a))*wf1(4*k - 1, a) + CONJG(wf1(4*j - 1, a))*wf1(4*k, a))&
 !                 &+ ui*(CONJG(wf1(4*k - 1, a))*wf1(4*j, a) + CONJG(wf1(4*k, a))*wf1(4*j - 1, a))
 !         ENDDO

 !         !!! Hrsoc of rotation
 !         !DO l = 1, 2
 !         !  k = p(l)%i
 !         !  j = p(l)%f
 !         !  Hrsoc = Hrsoc  &
 !         !       &  + exp( ui*0.25d0*pi)*CONJG(wf1(2*j, a)    )*wf1(2*k - 1, a) &
 !         !       &  + exp( ui*0.75d0*pi)*CONJG(wf1(2*j - 1, a))*wf1(2*k, a)     &
 !         !       &  + exp(-ui*0.25d0*pi)*CONJG(wf1(2*k - 1, a))*wf1(2*j, a)     &
 !         !       &  + exp(-ui*0.75d0*pi)*CONJG(wf1(2*k, a)    )*wf1(2*j - 1, a)
 !         !ENDDO
 !         !DO l = 3, 4
 !         !  k = p(l)%i
 !         !  j = p(l)%f
 !         !  Hrsoc = Hrsoc   &
 !         !       &  + exp( ui*0.75d0*pi)*CONJG(wf1(2*j, a)    )*wf1(2*k - 1, a) &
 !         !       &  + exp( ui*0.25d0*pi)*CONJG(wf1(2*j - 1, a))*wf1(2*k, a)     &
 !         !       &  + exp(-ui*0.75d0*pi)*CONJG(wf1(2*k - 1, a))*wf1(2*j, a)     &
 !         !       &  + exp(-ui*0.25d0*pi)*CONJG(wf1(2*k, a)    )*wf1(2*j - 1, a)
 !         !ENDDO
 !       ENDDO
 !     ENDDO
 !     !temp = 0.d0
 !     !DO i = 1, SIZE(hole(:, 1))
 !     !  p(1)%i = hole(i, 1) - nx
 !     !  p(1)%f = hole(i, 1) + 1
 !     !  p(2)%i = hole(i, 1) - 1
 !     !  p(2)%f = hole(i, 1) + nx
 !     !  p(3)%i = hole(i, 1) - nx
 !     !  p(3)%f = hole(i, 1) - 1
 !     !  p(4)%i = hole(i, 1) + 1
 !     !  p(4)%f = hole(i, 1) + nx
 !     !  DO a = 1, ne
 !     !      k = p(4)%i
 !     !      j = p(4)%f
 !     !      !temp = temp + ui*(CONJG(wf1(2*j, a))*wf1(2*k - 1, a))
 !     !      temp = temp + ui*(CONJG(wf1(2*j - 1, a))*wf1(2*k, a))
 !     !  ENDDO
 !     !ENDDO
 !   end do !z

 !   !write(*, *) 'Hrsoc', lm*Hrsoc
 !   !Hrsoc = lm*(Hrsoc + conjg(Hrsoc))
 !   !write(*, *) '--- c --- Hrsoc', lm*Hrsoc
 !   Hrsoc = lm*Hrsoc
 ! end function energy_spin_orbit_interaction_energy

!  function energy_spin_orbit_interaction_energy2(pshape, n, ne, hole, lm, &
!    & xi, chi, wf, EV) result(Hrsoc)
!    integer, intent(in)             :: hole(:, :), pshape(2), n, ne
!    REAL(8), INTENT(in)             :: lm, xi(n), chi(n)
!    COMPLEX(8), INTENT(in)          :: wf(:, :), EV ! cpwf
!    real(8)                         :: total_E
!    COMPLEX(8)                      :: Hrsoc, temp
!    REAL(8)                         :: eg(n)
!    INTEGER                         :: i, j, k, l, a, b, &
!      &                                nx
!    type(path), allocatable         :: p(:)
!    Hrsoc = 0.d0
!    nx = pshape(1)
!    Hrsoc = 0.d0
!    !expectation value of Hrsoc
!    IF (ALLOCATED(p)) DEALLOCATE (p)
!    ALLOCATE (p(4))
!    DO i = 1, SIZE(hole(:, 1))
!      p(1)%i = hole(i, 1) - nx
!      p(1)%f = hole(i, 1) + 1
!      p(2)%i = hole(i, 1) - 1
!      p(2)%f = hole(i, 1) + nx
!      p(3)%i = hole(i, 1) - nx
!      p(3)%f = hole(i, 1) - 1
!      p(4)%i = hole(i, 1) + 1
!      p(4)%f = hole(i, 1) + nx
!      DO a = 1, ne
!          Hrsoc =  Hrsoc + &
!              &      (exp(ui*0.5d0*(chi(p(1)%f) - chi(p(1)%i))) &
!              &       * exp(-0.5d0*ui*(xi(p(1)%f)-xi(p(1)%f-nx)+xi(p(1)%i)-xi(p(1)%i+1))) &
!              &       * conjg(wf(2*p(1)%f, a))*wf(2*p(1)%i - 1, a)  &
!              &       * conjg(EV) * exp(-ui*xi(p(1)%f-nx)) &! Dd* Du, h-y -> h+x
!              &       - exp(ui*0.5d0*(chi(p(1)%f) - chi(p(1)%i))) &
!              &       * exp(0.5d0*ui*(xi(p(1)%f)-xi(p(1)%f-nx)+xi(p(1)%i)-xi(p(1)%i+1))) &
!              &       * conjg(wf(2*p(1)%f-1, a))*wf(2*p(1)%i, a)  &
!              &       * EV * exp(ui*xi(p(1)%f-nx)) &! Du* Dd, h-y -> h+x
!              &       + exp(ui*0.5d0*(chi(p(2)%f) - chi(p(2)%i))) &
!              &       * exp(-0.5d0*ui*(xi(p(2)%f)-xi(p(2)%f-1)+xi(p(2)%i)-xi(p(2)%i+nx))) &
!              &       * conjg(wf(2*p(2)%f, a))*wf(2*p(2)%i - 1, a)  &
!              &       * conjg(EV) * exp(-ui*xi(p(2)%f-1)) &! Dd* Du, h-x -> h+y
!              &       - exp(ui*0.5d0*(chi(p(2)%f) - chi(p(2)%i))) &
!              &       * exp(0.5d0*ui*(xi(p(2)%f)-xi(p(2)%f-1)+xi(p(2)%i)-xi(p(2)%i+nx))) &
!              &       * conjg(wf(2*p(2)%f-1, a))*wf(2*p(2)%i, a)  &
!              &       * EV * exp(ui*xi(p(2)%f-1)) &! , h-x,d -> h+y,u
!              &       + ui*( exp(ui*0.5d0*(chi(p(3)%f) - chi(p(3)%i)))  &
!              &       * exp(-0.5d0*ui*(xi(p(3)%f)-xi(p(3)%f-nx)+xi(p(3)%i)-xi(p(3)%i-1))) & !!
!              &       * conjg(wf(2*p(3)%f, a))*wf(2*p(3)%i - 1, a)  &
!              &       * conjg(EV) * exp(-ui*xi(p(3)%f-nx)) &! h-y, u -> h-x,d
!              &       + exp(ui*0.5d0*(chi(p(3)%f) - chi(p(3)%i))) &
!              &       * exp(0.5d0*ui*(xi(p(3)%f)-xi(p(3)%f-nx)+xi(p(3)%i)-xi(p(3)%i-1))) &
!              &       * conjg(wf(2*p(3)%f-1, a))*wf(2*p(3)%i, a)  &
!              &       * EV * exp(ui*xi(p(3)%f-nx)) &! h-y d -> h-x u
!              &       + exp(ui*0.5d0*(chi(p(4)%f) - chi(p(4)%i))) &
!              &       * exp(-0.5d0*ui*(xi(p(4)%f)-xi(p(4)%f+1)+xi(p(4)%i)-xi(p(4)%i+nx))) & !!
!              &       * conjg(wf(2*p(4)%f, a))*wf(2*p(4)%i - 1, a)  & !!
!              &       * conjg(EV) * exp(-ui*xi(p(4)%f+1)) &! h+x u -> h+y d !!
!              &       + exp(ui*0.5d0*(chi(p(4)%f) - chi(p(4)%i))) &
!              &       * exp(0.5d0*ui*(xi(p(4)%f)-xi(p(4)%f+1)+xi(p(4)%i)-xi(p(4)%i+nx))) & !!
!              &       * conjg(wf(2*p(4)%f-1, a))*wf(2*p(4)%i, a)  & !!
!              &       * EV * exp(ui*xi(p(4)%f+1))   ) &! , h+x d -> h+y u !!
!            & )
!      ENDDO
!    ENDDO
!    !temp = 0.d0
!    !DO i = 1, SIZE(hole(:, 1))
!    !  p(1)%i = hole(i, 1) - nx
!    !  p(1)%f = hole(i, 1) + 1
!    !  p(2)%i = hole(i, 1) - 1
!    !  p(2)%f = hole(i, 1) + nx
!    !  p(3)%i = hole(i, 1) - nx
!    !  p(3)%f = hole(i, 1) - 1
!    !  p(4)%i = hole(i, 1) + 1
!    !  p(4)%f = hole(i, 1) + nx
!    !  DO a = 1, ne
!    !      temp =  temp  &
!    !          &       + ui*exp(ui*0.5d0*(chi(p(4)%f) - chi(p(4)%i))) &
!    !          &       * exp(0.5d0*ui*(xi(p(4)%f)-xi(p(4)%f+1)+xi(p(4)%i)-xi(p(4)%i+nx))) & !!
!    !          &       * conjg(wf(2*p(4)%f-1, a))*wf(2*p(4)%i, a)  & !!
!    !          &       *  exp(ui*xi(p(4)%f+1))! , h+x d -> h+y u !!
!    !  ENDDO
!    !ENDDO
!    !write(*, *) 'temp2', temp
!    Hrsoc = lm*(Hrsoc + conjg(Hrsoc))
!  end function energy_spin_orbit_interaction_energy2

 ! subroutine calc_eigen_energy_sb_2lay_xichi(sys, eg, wf, mu, delta, ne_u, ne_d, Sx, Sy, Sz, xi, chi)
 !   ! z =1 : surface
 !   ! z =2 : bulk
 !   implicit none
 !   class(base), intent(in)   :: sys
 !   real(8), dimension(:), allocatable, intent(out) :: eg
 !   complex(8),dimension(4*sys%n_acc, 4*sys%n_acc), intent(in)   :: wf
 !   real(8),intent(in)    :: mu(2)
 !   complex(8), dimension(sys%pl_1st_l(1)%length()), intent(in) :: delta
 !   real(8), dimension(sys%n), intent(in)   :: ne_u, ne_d
 !   !real(8), dimension(sys%n), intent(in) :: xi, chi
 !   real(8), dimension(sys%n), intent(in) :: Sx, Sy, Sz, xi, chi
 !   complex(8),DIMENSION(4*sys%n_acc, 4*sys%n_acc)   :: twf
 !   integer :: i, is
 !   twf = wf
 !   do i = 1, sys%n_acc
 !     is = sys%ets(i)
 !     twf(4*i-3, :) = exp(-0.5d0*ui*chi(is))*exp( 0.5d0*ui*xi(is))*twf(4*i-3, :)
 !     twf(4*i-2, :) = exp(-0.5d0*ui*chi(is))*exp(-0.5d0*ui*xi(is))*twf(4*i-2, :)
 !     twf(4*i-1, :) = exp( 0.5d0*ui*chi(is))*exp(-0.5d0*ui*xi(is))*twf(4*i-1, :)
 !     twf(4*i  , :) = exp( 0.5d0*ui*chi(is))*exp( 0.5d0*ui*xi(is))*twf(4*i  , :)
 !   end do
 !   !call calc_eigen_energy_sb_2lay(sys, eg, twf, mu, delta, ne_u, ne_d, Sx, Sy, Sz)
 !   !call calc_eigen_energy_sb_2lay(sys, eg, twf, mu, delta, ne_u, ne_d, Sx, Sy, Sz, wf, xi)
 !   call calc_eigen_energy_sb_2lay(sys, eg, twf, mu, tdelta_for, tdelta_rev, ne_u, ne_d, Sx, Sy, Sz, wf)
 ! end subroutine calc_eigen_energy_sb_2lay_xichi
  
  subroutine calc_eigen_energy_sb_2lay_wf_T0(sys, eg, wf, xi, chi, mu)
    class(base), intent(in)   :: sys
    real(8), dimension(:), allocatable, intent(out) :: eg
    complex(8),dimension(4*sys%n_acc, 4*sys%n_acc), intent(in)   :: wf
    real(8),intent(in)    :: mu(2)
    !complex(8), dimension(sys%pl_1st_l(1)%length()) :: delta
    complex(8), dimension(sys%pl_1st_l(1)%length()) :: tdelta_for, tdelta_rev
    real(8), dimension(sys%n)   :: ne_u, ne_d
    real(8), dimension(sys%n), intent(in) :: xi, chi
    real(8), dimension(sys%n) :: Sx, Sy, Sz
    call energy_mod_calc_spin(sys, wf, Sx, Sy, Sz, xi)
    call energy_mod_calc_ne(sys, wf, ne_u, ne_d)
    call energy_mod_calc_tdelta(sys, wf, tdelta_for, tdelta_rev, xi)
    call calc_eigen_energy_sb_2lay(sys, eg, wf, xi, chi, mu, tdelta_for, tdelta_rev, ne_u, ne_d, Sx, Sy, Sz)
  end subroutine calc_eigen_energy_sb_2lay_wf_T0
  
  subroutine energy_mod_calc_spin(sys, wf, Sx, Sy, Sz, xi)
    class(base), intent(in) :: sys
    complex(8),dimension(4*sys%n_acc, 4*sys%n_acc), intent(in)   :: wf
    real(8), dimension(:), intent(out) :: Sx, Sy, Sz
    real(8), dimension(:), intent(in) :: xi
    !real(8), intent(in) :: kbt
    integer :: i, is, a
    complex(8) :: uu, ud, vu, vd
    Sx = 0d0
    Sy = 0d0
    Sz = 0d0
    do i = 1, sys%n_acc
      is = sys%ets(i)
      do a = 2*sys%n_acc + 1, 4*sys%n_acc 
        !uu = wf(4*i-3, a)*exp(-0.5d0*ui*xi(is))
        !ud = wf(4*i-2, a)*exp( 0.5d0*ui*xi(is))
        vu = wf(4*i-1, a)*exp( 0.5d0*ui*xi(is))
        vd = wf(4*i  , a)*exp(-0.5d0*ui*xi(is))
        Sx(is) = Sx(is) - 0.5d0*(vu*conjg(vd) + conjg(vu)*vd)
        Sy(is) = Sy(is) + 0.5d0*ui*(vu*conjg(vd) - conjg(vu)*vd) 
        Sz(is) = Sz(is) + 0.5d0*(vu*conjg(vu) - vd*conjg(vd)) 
        !Sx(is) = Sx(is) + 0.5d0*(uu*conjg(ud) + conjg(uu)*ud)*fermi( this%eg(a), kbT) &
        !       &        - 0.5d0*(vu*conjg(vd) + conjg(vu)*vd)*fermi(-this%eg(a), kbT) 
        !Sy(is) = Sy(is) - 0.5d0*ui*(uu*conjg(ud) - conjg(uu)*ud)*fermi( this%eg(a), kbT) &
        !       &        + 0.5d0*ui*(vu*conjg(vd) - conjg(vu)*vd)*fermi(-this%eg(a), kbT) 
        !Sz(is) = Sz(is) + 0.5d0*(uu*conjg(uu) - ud*conjg(ud))*fermi( this%eg(a), kbT) &
        !       &        + 0.5d0*(vu*conjg(vu) - vd*conjg(vd))*fermi(-this%eg(a), kbT) 
      end do
    end do
  end subroutine energy_mod_calc_spin
  
  subroutine energy_mod_calc_ne(sys, wf, ne_u, ne_d)
    class(base), intent(in) :: sys
    complex(8),dimension(4*sys%n_acc, 4*sys%n_acc), intent(in)   :: wf
    real(8), dimension(:), intent(out) :: ne_u, ne_d
    !real(8), intent(in) :: kbt
    integer :: i, is, a
    complex(8) :: uu, ud, vu, vd
    ne_u = 0d0
    ne_d = 0d0
    do i = 1, sys%n_acc
      is = sys%ets(i)
      do a = 2*sys%n_acc + 1, 4*sys%n_acc 
        ne_u(is) = ne_u(is) + conjg(wf(4*i-1, a))*wf(4*i-1, a)
        ne_d(is) = ne_d(is) + conjg(wf(4*i  , a))*wf(4*i  , a) 
      end do
    end do
  end subroutine energy_mod_calc_ne
  
  subroutine energy_mod_calc_tdelta(sys, wf, tdelta_for, tdelta_rev, xi)
    class(base), intent(in) :: sys
    complex(8),dimension(4*sys%n_acc, 4*sys%n_acc), intent(in)   :: wf
    !complex(8), dimension(:), intent(out) :: delta
    complex(8), dimension(:), intent(out) :: tdelta_for, tdelta_rev
    real(8), dimension(:), intent(in) :: xi
    !real(8), intent(in) :: kbt
    integer :: i, j, k, a, is, js
    complex(8) :: exi, cexi
    
    tdelta_for  = 0.d0
    tdelta_rev  = 0.d0
    do a = 2*sys%n_acc+1, 4*sys%n_acc
      do k = 1, sys%pl_1st_l(1)%length()
        i = sys%pl_1st_l(1)%value(k)%i
        j = sys%pl_1st_l(1)%value(k)%f
        is = i
        js = j
        exi = exp(-ui*(xi(is)-xi(js)))
        !exi = exp(ui*(xi(is)-xi(js)))
        cexi = conjg(exi)
        ! jd : 2*t(1)**2 / U
        tdelta_for(k) = tdelta_for(k) &
                  & + (sys%jd*(      wf(4*i-3, a)*conjg(wf(4*j  , a)) &
                  &           +      wf(4*j-2, a)*conjg(wf(4*i-1, a)) &
                  &           +  exi*wf(4*i-2, a)*conjg(wf(4*j-1, a)) &
                  &           +  exi*wf(4*j-3, a)*conjg(wf(4*i  , a)))) ! *tanh(this%eg(a)/kbT*0.5d0)
        tdelta_rev(k) = tdelta_rev(k) &
                  & + (sys%jd*( cexi*wf(4*i-3, a)*conjg(wf(4*j  , a)) &
                  &           + cexi*wf(4*j-2, a)*conjg(wf(4*i-1, a)) &
                  &           +      wf(4*i-2, a)*conjg(wf(4*j-1, a)) &
                  &           +      wf(4*j-3, a)*conjg(wf(4*i  , a)))) ! *tanh(this%eg(a)/kbT*0.5d0)
      end do!k
    end do!a
  end subroutine energy_mod_calc_tdelta

  !subroutine calc_eigen_energy_sb_2lay(sys, eg, wf, mu, delta, ne_u, ne_d, Sx, Sy, Sz)
  !subroutine calc_eigen_energy_sb_2lay(sys, eg, wf, mu, delta, ne_u, ne_d, Sx, Sy, Sz, hfwf, xi)
  subroutine calc_eigen_energy_sb_2lay(sys, eg, wf, xi, chi, mu, tdelta_for, tdelta_rev, ne_u, ne_d, Sx, Sy, Sz)
    ! z =1 : surface
    ! z =2 : bulk
    implicit none
    class(base), intent(in)   :: sys
    real(8), dimension(:), allocatable, intent(out) :: eg
    !real(8), dimension(4*sys%n_acc), intent(out) :: eg
    complex(8),dimension(4*sys%n_acc, 4*sys%n_acc), intent(in)   :: wf
    real(8),intent(in)    :: mu(2)
    !complex(8), dimension(sys%pl_1st_l(1)%length()), intent(in) :: delta
    complex(8), dimension(sys%pl_1st_l(1)%length()), intent(in) :: tdelta_for, tdelta_rev
    real(8), dimension(sys%n), intent(in)   :: ne_u, ne_d
    real(8), dimension(sys%n), intent(in) :: xi, chi
    !complex(8),DIMENSION(4*sys%n_acc, 4*sys%n_acc), intent(in)   :: hfwf
    !real(8), dimension(sys%n), intent(in) :: xi
    real(8), dimension(sys%n), intent(in) :: Sx, Sy, Sz
    complex(8),DIMENSION(4*sys%n_acc, 4*sys%n_acc)   :: twf
    !real(8), dimension(sys%n) :: ne_u, ne_d
    integer :: i, j, k, is, js, ks, r
    real(8) :: t3U
    !complex(8) :: exi,echi, cexi

    !! error cheack
    if(sys%nz /= 2 .or. sys%nh_l(1) /= 0) then
      print *, "nz : ", sys%nz
      print *, "nh(1) : ", sys%nh_l(1)
      stop "error : wrong data in calc_eg_sb_2lay_xi_hamiltonian"
      ! sys subroutine is available for only 2 layer"
    endif

    !twf = wf
    do i = 1, sys%n_acc
      is = sys%ets(i)
      twf(4*i-3, :) = exp(-0.5d0*ui*chi(is))*exp(-0.5d0*ui*xi(is))*wf(4*i-3, :)
      twf(4*i-2, :) = exp(-0.5d0*ui*chi(is))*exp( 0.5d0*ui*xi(is))*wf(4*i-2, :)
      twf(4*i-1, :) = exp( 0.5d0*ui*chi(is))*exp( 0.5d0*ui*xi(is))*wf(4*i-1, :)
      twf(4*i  , :) = exp( 0.5d0*ui*chi(is))*exp(-0.5d0*ui*xi(is))*wf(4*i  , :)
    end do

    t3U = 2*sys%t(3)**2/sys%U

    !ne_u = 0.5d0
    !ne_d = 0.5d0
    allocate(eg(4*sys%n_acc))
    eg = 0d0

    do r = 1, 4*sys%n_acc
      !!!!!!!!!!!!!!!!!!    surface    !!!!!!!!!!!!!!!!!!
      do i = sys%n_acc_start(1), sys%n_acc_end(1)
        eg(r) = eg(r) - mu(1)*conjg(twf(4*i-3, r))*twf(4*i-3, r)
        eg(r) = eg(r) - mu(1)*conjg(twf(4*i-2, r))*twf(4*i-2, r)
        eg(r) = eg(r) + mu(1)*conjg(twf(4*i-1, r))*twf(4*i-1, r)
        eg(r) = eg(r) + mu(1)*conjg(twf(4*i  , r))*twf(4*i  , r)
      end do
      do i = 1, sys%pl_1st_l(1)%length()
        js = sys%pl_1st_l(1)%value(i)%i
        ks = sys%pl_1st_l(1)%value(i)%f
        j = sys%ste(js) ! in the 1st layer (surface), js = j, ks = k
        k = sys%ste(ks)
        eg(r) = eg(r) - sys%t(1)*(1d0-ne_d(js))*(1d0-ne_d(ks))*conjg(twf(4*j-3, r))*twf(4*k-3, r)
        eg(r) = eg(r) - sys%t(1)*(1d0-ne_u(js))*(1d0-ne_u(ks))*conjg(twf(4*j-2, r))*twf(4*k-2, r)
        eg(r) = eg(r) + sys%t(1)*(1d0-ne_d(js))*(1d0-ne_d(ks))*conjg(twf(4*j-1, r))*twf(4*k-1, r)
        eg(r) = eg(r) + sys%t(1)*(1d0-ne_u(js))*(1d0-ne_u(ks))*conjg(twf(4*j  , r))*twf(4*k  , r)
        eg(r) = eg(r) - sys%t(1)*(1d0-ne_d(js))*(1d0-ne_d(ks))*conjg(twf(4*k-3, r))*twf(4*j-3, r)
        eg(r) = eg(r) - sys%t(1)*(1d0-ne_u(js))*(1d0-ne_u(ks))*conjg(twf(4*k-2, r))*twf(4*j-2, r)
        eg(r) = eg(r) + sys%t(1)*(1d0-ne_d(js))*(1d0-ne_d(ks))*conjg(twf(4*k-1, r))*twf(4*j-1, r)
        eg(r) = eg(r) + sys%t(1)*(1d0-ne_u(js))*(1d0-ne_u(ks))*conjg(twf(4*k  , r))*twf(4*j  , r)
      
      end do
      do i = 1, sys%pl_2nd_l(1)%length()
        js = sys%pl_2nd_l(1)%value(i)%i
        ks = sys%pl_2nd_l(1)%value(i)%f
        j = sys%ste(js) ! in the 1st layer (surface), js = j, ks = k
        k = sys%ste(ks)
        eg(r) = eg(r) - sys%t(2)*(1d0-ne_d(js))*(1d0-ne_d(ks))*conjg(twf(4*j-3, r))*twf(4*k-3, r)
        eg(r) = eg(r) - sys%t(2)*(1d0-ne_u(js))*(1d0-ne_u(ks))*conjg(twf(4*j-2, r))*twf(4*k-2, r)
        eg(r) = eg(r) + sys%t(2)*(1d0-ne_d(js))*(1d0-ne_d(ks))*conjg(twf(4*j-1, r))*twf(4*k-1, r)
        eg(r) = eg(r) + sys%t(2)*(1d0-ne_u(js))*(1d0-ne_u(ks))*conjg(twf(4*j  , r))*twf(4*k  , r)
        eg(r) = eg(r) - sys%t(2)*(1d0-ne_d(js))*(1d0-ne_d(ks))*conjg(twf(4*k-3, r))*twf(4*j-3, r)
        eg(r) = eg(r) - sys%t(2)*(1d0-ne_u(js))*(1d0-ne_u(ks))*conjg(twf(4*k-2, r))*twf(4*j-2, r)
        eg(r) = eg(r) + sys%t(2)*(1d0-ne_d(js))*(1d0-ne_d(ks))*conjg(twf(4*k-1, r))*twf(4*j-1, r)
        eg(r) = eg(r) + sys%t(2)*(1d0-ne_u(js))*(1d0-ne_u(ks))*conjg(twf(4*k  , r))*twf(4*j  , r)
      end do
      do i = 1, sys%pl_1st_l(1)%length()
        ! here, delta's index corresponds to pl_1st_l(1)'s index
        js = sys%pl_1st_l(1)%value(i)%i
        ks = sys%pl_1st_l(1)%value(i)%f
        j = sys%ste(js) ! in the 1st layer ( surface ), js = j, ks = k
        k = sys%ste(ks)
       ! exi  = exp(-0.5d0*ui*modulo(xi(js)-xi(ks), 2d0*pi))
       ! !exi  = 1d0
       ! cexi = conjg(exi)
      !  eg(r) = eg(r) +        delta(i)*conjg(twf(4*j-3, r))*twf(4*k  , r)
      !  eg(r) = eg(r) +        delta(i)*conjg(twf(4*j-2, r))*twf(4*k-1, r)
      !  eg(r) = eg(r) + conjg(delta(i))*conjg(twf(4*j-1, r))*twf(4*k-2, r)
      !  eg(r) = eg(r) + conjg(delta(i))*conjg(twf(4*j  , r))*twf(4*k-3, r)
      !  eg(r) = eg(r) +        delta(i)*conjg(twf(4*k-3, r))*twf(4*j  , r)
      !  eg(r) = eg(r) +        delta(i)*conjg(twf(4*k-2, r))*twf(4*j-1, r)
      !  eg(r) = eg(r) + conjg(delta(i))*conjg(twf(4*k-1, r))*twf(4*j-2, r)
      !  eg(r) = eg(r) + conjg(delta(i))*conjg(twf(4*k  , r))*twf(4*j-3, r)
      
      !  eg(r) = eg(r) +        delta(i)*conjg(hfwf(4*j-3, r))*hfwf(4*k  , r)*exi
      !  eg(r) = eg(r) +        delta(i)*conjg(hfwf(4*j-2, r))*hfwf(4*k-1, r)*cexi
      !  eg(r) = eg(r) + conjg(delta(i))*conjg(hfwf(4*j-1, r))*hfwf(4*k-2, r)*cexi
      !  eg(r) = eg(r) + conjg(delta(i))*conjg(hfwf(4*j  , r))*hfwf(4*k-3, r)*exi
      !  eg(r) = eg(r) +        delta(i)*conjg(hfwf(4*k-3, r))*hfwf(4*j  , r)*cexi
      !  eg(r) = eg(r) +        delta(i)*conjg(hfwf(4*k-2, r))*hfwf(4*j-1, r)*exi
      !  eg(r) = eg(r) + conjg(delta(i))*conjg(hfwf(4*k-1, r))*hfwf(4*j-2, r)*exi
      !  eg(r) = eg(r) + conjg(delta(i))*conjg(hfwf(4*k  , r))*hfwf(4*j-3, r)*cexi
      
        eg(r) = eg(r) +       tdelta_for(i) *conjg(wf(4*j-3, r))*wf(4*k  , r)
        eg(r) = eg(r) +       tdelta_rev(i) *conjg(wf(4*j-2, r))*wf(4*k-1, r)
        eg(r) = eg(r) + conjg(tdelta_for(i))*conjg(wf(4*j-1, r))*wf(4*k-2, r)
        eg(r) = eg(r) + conjg(tdelta_rev(i))*conjg(wf(4*j  , r))*wf(4*k-3, r)
        eg(r) = eg(r) +       tdelta_rev(i) *conjg(wf(4*k-3, r))*wf(4*j  , r)
        eg(r) = eg(r) +       tdelta_for(i) *conjg(wf(4*k-2, r))*wf(4*j-1, r)
        eg(r) = eg(r) + conjg(tdelta_rev(i))*conjg(wf(4*k-1, r))*wf(4*j-2, r)
        eg(r) = eg(r) + conjg(tdelta_for(i))*conjg(wf(4*k  , r))*wf(4*j-3, r)
        
        !eg(r) = eg(r) +       tdelta_for(i) *conjg(wf(4*j-3, r))*wf(4*k  , r)
        !eg(r) = eg(r) +       tdelta_rev(i) *conjg(wf(4*j-2, r))*wf(4*k-1, r)
        !eg(r) = eg(r) + conjg(tdelta_rev(i))*conjg(wf(4*j-1, r))*wf(4*k-2, r)
        !eg(r) = eg(r) + conjg(tdelta_for(i))*conjg(wf(4*j  , r))*wf(4*k-3, r)
        !eg(r) = eg(r) +       tdelta_rev(i) *conjg(wf(4*k-3, r))*wf(4*j  , r)
        !eg(r) = eg(r) +       tdelta_for(i) *conjg(wf(4*k-2, r))*wf(4*j-1, r)
        !eg(r) = eg(r) + conjg(tdelta_for(i))*conjg(wf(4*k-1, r))*wf(4*j-2, r)
        !eg(r) = eg(r) + conjg(tdelta_rev(i))*conjg(wf(4*k  , r))*wf(4*j-3, r)
      end do
      !!!!!!!!!!!!!!    end surface    !!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!    bulk    !!!!!!!!!!!!!!!!
      do i = sys%n_acc_start(2), sys%n_acc_end(2)
        is = sys%ets(i)
        eg(r) = eg(r) + sys%U*(-2d0/3d0*Sz(is)+0.5d0)*conjg(twf(4*i-3, r))*twf(4*i-3, r)
        eg(r) = eg(r) + sys%U*( 2d0/3d0*Sz(is)+0.5d0)*conjg(twf(4*i-2, r))*twf(4*i-2, r)
        eg(r) = eg(r) - sys%U*(-2d0/3d0*Sz(is)+0.5d0)*conjg(twf(4*i-1, r))*twf(4*i-1, r)
        eg(r) = eg(r) - sys%U*( 2d0/3d0*Sz(is)+0.5d0)*conjg(twf(4*i  , r))*twf(4*i  , r)
        eg(r) = eg(r) - sys%U*2d0/3d0*(Sx(is)-ui*Sy(is))*conjg(twf(4*i-3, r))*twf(4*i-2, r)
        eg(r) = eg(r) - sys%U*2d0/3d0*(Sx(is)+ui*Sy(is))*conjg(twf(4*i-2, r))*twf(4*i-3, r)
        eg(r) = eg(r) - sys%U*2d0/3d0*(Sx(is)+ui*Sy(is))*conjg(twf(4*i-1, r))*twf(4*i  , r)
        eg(r) = eg(r) - sys%U*2d0/3d0*(Sx(is)-ui*Sy(is))*conjg(twf(4*i  , r))*twf(4*i-1, r)
      end do
      do i = sys%n_acc_start(2), sys%n_acc_end(2)
        eg(r) = eg(r) - mu(2)*conjg(twf(4*i-3, r))*twf(4*i-3, r)
        eg(r) = eg(r) - mu(2)*conjg(twf(4*i-2, r))*twf(4*i-2, r)
        eg(r) = eg(r) + mu(2)*conjg(twf(4*i-1, r))*twf(4*i-1, r)
        eg(r) = eg(r) + mu(2)*conjg(twf(4*i  , r))*twf(4*i  , r)
      end do
      do i = 1, sys%pl_1st_l(2)%length()
        js = sys%pl_1st_l(2)%value(i)%i
        ks = sys%pl_1st_l(2)%value(i)%f
        j = sys%ste(js)
        k = sys%ste(ks)
        eg(r) = eg(r) - sys%t(1)*conjg(twf(4*j-3, r))*twf(4*k-3, r)
        eg(r) = eg(r) - sys%t(1)*conjg(twf(4*j-2, r))*twf(4*k-2, r)
        eg(r) = eg(r) + sys%t(1)*conjg(twf(4*j-1, r))*twf(4*k-1, r)
        eg(r) = eg(r) + sys%t(1)*conjg(twf(4*j  , r))*twf(4*k  , r)
        eg(r) = eg(r) - sys%t(1)*conjg(twf(4*k-3, r))*twf(4*j-3, r)
        eg(r) = eg(r) - sys%t(1)*conjg(twf(4*k-2, r))*twf(4*j-2, r)
        eg(r) = eg(r) + sys%t(1)*conjg(twf(4*k-1, r))*twf(4*j-1, r)
        eg(r) = eg(r) + sys%t(1)*conjg(twf(4*k  , r))*twf(4*j  , r)
      end do
      do i = 1, sys%pl_2nd_l(2)%length()
        js = sys%pl_2nd_l(2)%value(i)%i
        ks = sys%pl_2nd_l(2)%value(i)%f
        j = sys%ste(js)
        k = sys%ste(ks)
        eg(r) = eg(r) - sys%t(2)*conjg(twf(4*j-3, r))*twf(4*k-3, r)
        eg(r) = eg(r) - sys%t(2)*conjg(twf(4*j-2, r))*twf(4*k-2, r)
        eg(r) = eg(r) + sys%t(2)*conjg(twf(4*j-1, r))*twf(4*k-1, r)
        eg(r) = eg(r) + sys%t(2)*conjg(twf(4*j  , r))*twf(4*k  , r)
        eg(r) = eg(r) - sys%t(2)*conjg(twf(4*k-3, r))*twf(4*j-3, r)
        eg(r) = eg(r) - sys%t(2)*conjg(twf(4*k-2, r))*twf(4*j-2, r)
        eg(r) = eg(r) + sys%t(2)*conjg(twf(4*k-1, r))*twf(4*j-1, r)
        eg(r) = eg(r) + sys%t(2)*conjg(twf(4*k  , r))*twf(4*j  , r)
      end do
      do i = 1, sys%pl_h%length()
        js = sys%pl_h%value(i)%i
        ks = sys%pl_h%value(i)%f
        j = sys%ste(js)
        k = sys%ste(ks)
        eg(r) = eg(r) + 0.5d0*sys%jd*Sz(ks)*conjg(twf(4*j-3, r))*twf(4*j-3, r)
        eg(r) = eg(r) - 0.5d0*sys%jd*Sz(ks)*conjg(twf(4*j-2, r))*twf(4*j-2, r)
        eg(r) = eg(r) - 0.5d0*sys%jd*Sz(ks)*conjg(twf(4*j-1, r))*twf(4*j-1, r)
        eg(r) = eg(r) + 0.5d0*sys%jd*Sz(ks)*conjg(twf(4*j  , r))*twf(4*j  , r)
        eg(r) = eg(r) + 0.5d0*sys%jd*(Sx(ks) - ui*Sy(ks))*conjg(twf(4*j-3, r))*twf(4*j-2, r)
        eg(r) = eg(r) + 0.5d0*sys%jd*(Sx(ks) + ui*Sy(ks))*conjg(twf(4*j-2, r))*twf(4*j-3, r)
        eg(r) = eg(r) + 0.5d0*sys%jd*(Sx(ks) + ui*Sy(ks))*conjg(twf(4*j-1, r))*twf(4*j  , r)
        eg(r) = eg(r) + 0.5d0*sys%jd*(Sx(ks) - ui*Sy(ks))*conjg(twf(4*j  , r))*twf(4*j-1, r)
      
        eg(r) = eg(r) + 0.5d0*sys%jd*Sz(js)*conjg(twf(4*k-3, r))*twf(4*k-3, r)
        eg(r) = eg(r) - 0.5d0*sys%jd*Sz(js)*conjg(twf(4*k-2, r))*twf(4*k-2, r)
        eg(r) = eg(r) - 0.5d0*sys%jd*Sz(js)*conjg(twf(4*k-1, r))*twf(4*k-1, r)
        eg(r) = eg(r) + 0.5d0*sys%jd*Sz(js)*conjg(twf(4*k  , r))*twf(4*k  , r)
        eg(r) = eg(r) + 0.5d0*sys%jd*(Sx(js) - ui*Sy(js))*conjg(twf(4*k-3, r))*twf(4*k-2, r)
        eg(r) = eg(r) + 0.5d0*sys%jd*(Sx(js) + ui*Sy(js))*conjg(twf(4*k-2, r))*twf(4*k-3, r)
        eg(r) = eg(r) + 0.5d0*sys%jd*(Sx(js) + ui*Sy(js))*conjg(twf(4*k-1, r))*twf(4*k  , r)
        eg(r) = eg(r) + 0.5d0*sys%jd*(Sx(js) - ui*Sy(js))*conjg(twf(4*k  , r))*twf(4*k-1, r)
      end do
    
      do i = 1, sys%pl_rh_1%length()  ! rashba (upper-right path)
        js = sys%pl_rh_1%value(i)%i ! h-y   |   h-x
        ks = sys%pl_rh_1%value(i)%f ! h+x   |   h+y
        j = sys%ste(js)
        k = sys%ste(ks)
        eg(r) = eg(r) + sys%lm*conjg(twf(4*j-2, r))*twf(4*k-3, r)
        eg(r) = eg(r) - sys%lm*conjg(twf(4*j-3, r))*twf(4*k-2, r)
        eg(r) = eg(r) + sys%lm*conjg(twf(4*j  , r))*twf(4*k-1, r)
        eg(r) = eg(r) - sys%lm*conjg(twf(4*j-1, r))*twf(4*k  , r)
        eg(r) = eg(r) + sys%lm*conjg(twf(4*k-2, r))*twf(4*j-3, r)
        eg(r) = eg(r) - sys%lm*conjg(twf(4*k-3, r))*twf(4*j-2, r)
        eg(r) = eg(r) + sys%lm*conjg(twf(4*k  , r))*twf(4*j-1, r)
        eg(r) = eg(r) - sys%lm*conjg(twf(4*k-1, r))*twf(4*j  , r)
      end do
      do i = 1, sys%pl_rh_2%length()  ! rashba (upper-left path)
        js = sys%pl_rh_2%value(i)%i ! h-y   |   h+x
        ks = sys%pl_rh_2%value(i)%f ! h-x   |   h+y
        j = sys%ste(js)
        k = sys%ste(ks)
        eg(r) = eg(r) + ui*sys%lm*conjg(twf(4*j-2, r))*twf(4*k-3, r)
        eg(r) = eg(r) + ui*sys%lm*conjg(twf(4*j-3, r))*twf(4*k-2, r)
        eg(r) = eg(r) - ui*sys%lm*conjg(twf(4*j  , r))*twf(4*k-1, r)
        eg(r) = eg(r) - ui*sys%lm*conjg(twf(4*j-1, r))*twf(4*k  , r)
        eg(r) = eg(r) - ui*sys%lm*conjg(twf(4*k-2, r))*twf(4*j-3, r)
        eg(r) = eg(r) - ui*sys%lm*conjg(twf(4*k-3, r))*twf(4*j-2, r)
        eg(r) = eg(r) + ui*sys%lm*conjg(twf(4*k  , r))*twf(4*j-1, r)
        eg(r) = eg(r) + ui*sys%lm*conjg(twf(4*k-1, r))*twf(4*j  , r)
      end do
      !!!!!!!!!!!!!!    end bulk    !!!!!!!!!!!!!!!!!!
  !
      !!!!!!!!!!!!!!!!!!    interlayer    !!!!!!!!!!!!!!!!
      ! pl_z_l(i) is z-dir path list between ith layer and (i+1)th layer.
      ! here, pl_z_l(1) is interlayer path ( between 1st(surface) and 2nd(bulk) ).
      do i = 1, sys%pl_z_l(1)%length()
        js = sys%pl_z_l(1)%value(i)%i
        ks = sys%pl_z_l(1)%value(i)%f
        j = sys%ste(js)
        k = sys%ste(ks)
        eg(r) = eg(r) - sys%t(3)*(1d0-ne_d(js))*(1d0-ne_d(ks))*conjg(twf(4*j-3, r))*twf(4*k-3, r)
        eg(r) = eg(r) - sys%t(3)*(1d0-ne_u(js))*(1d0-ne_u(ks))*conjg(twf(4*j-2, r))*twf(4*k-2, r)
        eg(r) = eg(r) + sys%t(3)*(1d0-ne_d(js))*(1d0-ne_d(ks))*conjg(twf(4*j-1, r))*twf(4*k-1, r)
        eg(r) = eg(r) + sys%t(3)*(1d0-ne_u(js))*(1d0-ne_u(ks))*conjg(twf(4*j  , r))*twf(4*k  , r)
        eg(r) = eg(r) - sys%t(3)*(1d0-ne_d(js))*(1d0-ne_d(ks))*conjg(twf(4*k-3, r))*twf(4*j-3, r)
        eg(r) = eg(r) - sys%t(3)*(1d0-ne_u(js))*(1d0-ne_u(ks))*conjg(twf(4*k-2, r))*twf(4*j-2, r)
        eg(r) = eg(r) + sys%t(3)*(1d0-ne_d(js))*(1d0-ne_d(ks))*conjg(twf(4*k-1, r))*twf(4*j-1, r)
        eg(r) = eg(r) + sys%t(3)*(1d0-ne_u(js))*(1d0-ne_u(ks))*conjg(twf(4*k  , r))*twf(4*j  , r)
      
        eg(r) = eg(r) + t3U*(Sz(ks)-0.5d0*(ne_u(ks)+ne_d(ks)))*conjg(twf(4*j-3, r))*twf(4*j-3, r)
        eg(r) = eg(r) - t3U*(Sz(ks)+0.5d0*(ne_u(ks)+ne_d(ks)))*conjg(twf(4*j-2, r))*twf(4*j-2, r)
        eg(r) = eg(r) - t3U*(Sz(ks)-0.5d0*(ne_u(ks)+ne_d(ks)))*conjg(twf(4*j-1, r))*twf(4*j-1, r)
        eg(r) = eg(r) + t3U*(Sz(ks)+0.5d0*(ne_u(ks)+ne_d(ks)))*conjg(twf(4*j  , r))*twf(4*j  , r)
        eg(r) = eg(r) + t3U*(Sx(ks) - ui*Sy(ks))*conjg(twf(4*j-3, r))*twf(4*j-2, r)
        eg(r) = eg(r) + t3U*(Sx(ks) + ui*Sy(ks))*conjg(twf(4*j-2, r))*twf(4*j-3, r)
        eg(r) = eg(r) + t3U*(Sx(ks) + ui*Sy(ks))*conjg(twf(4*j-1, r))*twf(4*j  , r)
        eg(r) = eg(r) + t3U*(Sx(ks) - ui*Sy(ks))*conjg(twf(4*j  , r))*twf(4*j-1, r)
      
        eg(r) = eg(r) + t3U*(Sz(js)-0.5d0*(ne_u(js)+ne_d(js)))*conjg(twf(4*k-3, r))*twf(4*k-3, r)
        eg(r) = eg(r) - t3U*(Sz(js)+0.5d0*(ne_u(js)+ne_d(js)))*conjg(twf(4*k-2, r))*twf(4*k-2, r)
        eg(r) = eg(r) - t3U*(Sz(js)-0.5d0*(ne_u(js)+ne_d(js)))*conjg(twf(4*k-1, r))*twf(4*k-1, r)
        eg(r) = eg(r) + t3U*(Sz(js)+0.5d0*(ne_u(js)+ne_d(js)))*conjg(twf(4*k  , r))*twf(4*k  , r)
        eg(r) = eg(r) + t3U*(Sx(js) - ui*Sy(js))*conjg(twf(4*k-3, r))*twf(4*k-2, r)
        eg(r) = eg(r) + t3U*(Sx(js) + ui*Sy(js))*conjg(twf(4*k-2, r))*twf(4*k-3, r)
        eg(r) = eg(r) + t3U*(Sx(js) + ui*Sy(js))*conjg(twf(4*k-1, r))*twf(4*k  , r)
        eg(r) = eg(r) + t3U*(Sx(js) - ui*Sy(js))*conjg(twf(4*k  , r))*twf(4*k-1, r)
      end do
      !!!!!!!!!!!!!!   end  interlayer    !!!!!!!!!!!!!!!!
  !  
  !    call my_zheev("l", ham, sys%eg, sys%wf)

    end do
  end subroutine calc_eigen_energy_sb_2lay
  
  function total_energy_sb_2lay_xichi_T0(sys, wf, xi, chi) result(r)
    implicit none
    class(base), intent(in)   :: sys
    complex(8),dimension(4*sys%n_acc, 4*sys%n_acc), intent(in)   :: wf
    !real(8),intent(in)    :: mu(2)
    !complex(8), dimension(sys%pl_1st_l(1)%length()), intent(in) :: delta
    !real(8), dimension(sys%n), intent(in)   :: ne_u, ne_d
    real(8), dimension(sys%n)             :: ne_u, ne_d
    real(8), dimension(sys%n), intent(in) :: xi, chi
    !real(8), dimension(sys%n), intent(in) :: Sx, Sy, Sz, xi, chi
    complex(8),DIMENSION(4*sys%n_acc, 4*sys%n_acc)   :: twf
    !real(8), dimension(:), allocatable :: eigen_eg
    real(8) :: r
    real(8) :: Ht_s, Ht2_s, Ht_b, Ht2_b, Hrsoc, H_interlayer
    integer :: i, is, j, js, k, ks, a
    r = 0d0
    twf = wf
    do i = 1, sys%n_acc
      is = sys%ets(i)
      twf(4*i-3, :) = exp(-0.5d0*ui*chi(is))*exp(-0.5d0*ui*xi(is))*twf(4*i-3, :)
      twf(4*i-2, :) = exp(-0.5d0*ui*chi(is))*exp( 0.5d0*ui*xi(is))*twf(4*i-2, :)
      twf(4*i-1, :) = exp( 0.5d0*ui*chi(is))*exp( 0.5d0*ui*xi(is))*twf(4*i-1, :)
      twf(4*i  , :) = exp( 0.5d0*ui*chi(is))*exp(-0.5d0*ui*xi(is))*twf(4*i  , :)
    end do    
    
    ne_u = 0d0
    ne_d = 0d0
    do i = 1, sys%n_acc
      is = sys%ets(i)
      do j = 2*sys%n_acc+1, 4*sys%n_acc
        ne_u(is) = ne_u(is) + conjg(twf(4*i-1, j))*twf(4*i-1, j)
        ne_d(is) = ne_d(is) + conjg(twf(4*i  , j))*twf(4*i  , j)
      end do
    end do
    !print *, ne_u(:)
    !stop    
   
    ! surface, the nearest hopping energy
    Ht_s = 0d0
    do i = 1, sys%pl_1st_l(1)%length()
      js = sys%pl_1st_l(1)%value(i)%i
      ks = sys%pl_1st_l(1)%value(i)%f
      j = sys%ste(js) ! in the 1st layer (surface), js = j, ks = k
      k = sys%ste(ks)
      do a = 2*sys%n_acc+1, 4*sys%n_acc
        Ht_s = Ht_s - sys%t(1)*(1d0-ne_d(js))*(1d0-ne_d(ks))*twf(4*j-1, a)*conjg(twf(4*k-1, a))
        Ht_s = Ht_s - sys%t(1)*(1d0-ne_u(js))*(1d0-ne_u(ks))*twf(4*j  , a)*conjg(twf(4*k  , a))
        Ht_s = Ht_s - sys%t(1)*(1d0-ne_d(js))*(1d0-ne_d(ks))*twf(4*k-1, a)*conjg(twf(4*j-1, a))
        Ht_s = Ht_s - sys%t(1)*(1d0-ne_u(js))*(1d0-ne_u(ks))*twf(4*k  , a)*conjg(twf(4*j  , a))
      end do
    end do

    ! surface, the 2nd nearest hopping energy
    Ht2_s = 0d0
    do i = 1, sys%pl_2nd_l(1)%length()
      js = sys%pl_2nd_l(1)%value(i)%i
      ks = sys%pl_2nd_l(1)%value(i)%f
      j = sys%ste(js) ! in the 1st layer (surface), js = j, ks = k
      k = sys%ste(ks)
      do a = 2*sys%n_acc+1, 4*sys%n_acc
        Ht2_s = Ht2_s - sys%t(2)*(1d0-ne_d(js))*(1d0-ne_d(ks))*twf(4*j-1, a)*conjg(twf(4*k-1, a))
        Ht2_s = Ht2_s - sys%t(2)*(1d0-ne_u(js))*(1d0-ne_u(ks))*twf(4*j  , a)*conjg(twf(4*k  , a))
        Ht2_s = Ht2_s - sys%t(2)*(1d0-ne_d(js))*(1d0-ne_d(ks))*twf(4*k-1, a)*conjg(twf(4*j-1, a))
        Ht2_s = Ht2_s - sys%t(2)*(1d0-ne_u(js))*(1d0-ne_u(ks))*twf(4*k  , a)*conjg(twf(4*j  , a))
      end do
    end do
    
    ! bulk, the nearest hopping energy
    Ht_b = 0d0
    do i = 1, sys%pl_1st_l(2)%length()
      js = sys%pl_1st_l(2)%value(i)%i
      ks = sys%pl_1st_l(2)%value(i)%f
      j = sys%ste(js) ! in the 1st layer (surface), js = j, ks = k
      k = sys%ste(ks)
      do a = 2*sys%n_acc+1, 4*sys%n_acc
        Ht_b = Ht_b - sys%t(1)*twf(4*j-1, a)*conjg(twf(4*k-1, a))
        Ht_b = Ht_b - sys%t(1)*twf(4*j  , a)*conjg(twf(4*k  , a))
        Ht_b = Ht_b - sys%t(1)*twf(4*k-1, a)*conjg(twf(4*j-1, a))
        Ht_b = Ht_b - sys%t(1)*twf(4*k  , a)*conjg(twf(4*j  , a))
      end do
    end do
    
    ! bulk, the 2nd nearest hopping energy
    Ht2_b = 0d0
    do i = 1, sys%pl_2nd_l(2)%length()
      js = sys%pl_2nd_l(2)%value(i)%i
      ks = sys%pl_2nd_l(2)%value(i)%f
      j = sys%ste(js) ! in the 1st layer (surface), js = j, ks = k
      k = sys%ste(ks)
      do a = 2*sys%n_acc+1, 4*sys%n_acc
        Ht2_b = Ht2_b - sys%t(2)*twf(4*j-1, a)*conjg(twf(4*k-1, a))
        Ht2_b = Ht2_b - sys%t(2)*twf(4*j  , a)*conjg(twf(4*k  , a))
        Ht2_b = Ht2_b - sys%t(2)*twf(4*k-1, a)*conjg(twf(4*j-1, a))
        Ht2_b = Ht2_b - sys%t(2)*twf(4*k  , a)*conjg(twf(4*j  , a))
      end do
    end do

    ! bulk, rashba energy
    Hrsoc = 0d0
    do i = 1, sys%pl_rh_1%length()
      js = sys%pl_rh_1%value(i)%i
      ks = sys%pl_rh_1%value(i)%f
      j = sys%ste(js)
      k = sys%ste(ks)
      do a = 2*sys%n_acc+1, 4*sys%n_acc
        Hrsoc = Hrsoc - sys%lm*twf(4*k  , a)*conjg(twf(4*j-1, a))
        Hrsoc = Hrsoc + sys%lm*twf(4*k-1, a)*conjg(twf(4*j  , a))
        Hrsoc = Hrsoc - sys%lm*twf(4*j-1, a)*conjg(twf(4*k  , a))
        Hrsoc = Hrsoc + sys%lm*twf(4*j  , a)*conjg(twf(4*k-1, a))
      end do
    end do
    do i = 1, sys%pl_rh_2%length() 
      js = sys%pl_rh_2%value(i)%i
      ks = sys%pl_rh_2%value(i)%f
      j = sys%ste(js)
      k = sys%ste(ks)
      do a = 2*sys%n_acc+1, 4*sys%n_acc
        Hrsoc = Hrsoc - ui*sys%lm*twf(4*k  , a)*conjg(twf(4*j-1, a))
        Hrsoc = Hrsoc - ui*sys%lm*twf(4*k-1, a)*conjg(twf(4*j  , a))
        Hrsoc = Hrsoc + ui*sys%lm*twf(4*j-1, a)*conjg(twf(4*k  , a))
        Hrsoc = Hrsoc + ui*sys%lm*twf(4*j  , a)*conjg(twf(4*k-1, a))
      end do
    end do
    
    ! inter-layer, hopping energy
    H_interlayer = 0d0
    do i = 1, sys%pl_z_l(1)%length()
      js = sys%pl_z_l(1)%value(i)%i
      ks = sys%pl_z_l(1)%value(i)%f
      j = sys%ste(js)
      k = sys%ste(ks)
      do a = 2*sys%n_acc+1, 4*sys%n_acc
        H_interlayer = H_interlayer - sys%t(3)*(1d0-ne_d(js))*(1d0-ne_d(ks))*conjg(twf(4*j-1, a))*twf(4*k-1, a)
        H_interlayer = H_interlayer - sys%t(3)*(1d0-ne_u(js))*(1d0-ne_u(ks))*conjg(twf(4*j  , a))*twf(4*k  , a)
        H_interlayer = H_interlayer - sys%t(3)*(1d0-ne_d(js))*(1d0-ne_d(ks))*conjg(twf(4*k-1, a))*twf(4*j-1, a)
        H_interlayer = H_interlayer - sys%t(3)*(1d0-ne_u(js))*(1d0-ne_u(ks))*conjg(twf(4*k  , a))*twf(4*j  , a)
      end do
    end do

    r = Ht_s + Ht2_s + Ht_b + Ht2_b + Hrsoc + H_interlayer

  end function total_energy_sb_2lay_xichi_T0
  
  subroutine calc_eigen_energy_sb_2lay_wf_T0_Aem(sys, eg, wf, xi, chi, mu, Aemdr)
    class(base), intent(in)   :: sys
    real(8), dimension(:), allocatable, intent(out) :: eg
    complex(8),dimension(4*sys%n_acc, 4*sys%n_acc), intent(in)   :: wf
    real(8),intent(in)    :: mu(2)
    real(8), dimension(sys%np), intent(in)    :: Aemdr    ! u_Aeff = d(chi) + Aemdr
                                                          ! Aemdr = - intgral{A^{em} dr}
    !complex(8), dimension(sys%pl_1st_l(1)%length()) :: delta
    complex(8), dimension(sys%pl_1st_l(1)%length()) :: tdelta_for, tdelta_rev
    real(8), dimension(sys%n)   :: ne_u, ne_d
    real(8), dimension(sys%n), intent(in) :: xi, chi
    real(8), dimension(sys%n) :: Sx, Sy, Sz
    call energy_mod_calc_spin(sys, wf, Sx, Sy, Sz, xi)
    call energy_mod_calc_ne(sys, wf, ne_u, ne_d)
    call energy_mod_calc_tdelta(sys, wf, tdelta_for, tdelta_rev, xi)
    call calc_eigen_energy_sb_2lay_Aem(sys, eg, wf, xi, chi, mu, Aemdr, tdelta_for, tdelta_rev, ne_u, ne_d, Sx, Sy, Sz)
  end subroutine calc_eigen_energy_sb_2lay_wf_T0_Aem

  subroutine calc_eigen_energy_sb_2lay_Aem(sys, eg, wf, xi, chi, mu, Aemdr, tdelta_for, tdelta_rev, ne_u, ne_d, Sx, Sy, Sz)
    ! z =1 : surface
    ! z =2 : bulk
    implicit none
    class(base), intent(in)   :: sys
    real(8), dimension(:), allocatable, intent(out) :: eg
    !real(8), dimension(4*sys%n_acc), intent(out) :: eg
    complex(8),dimension(4*sys%n_acc, 4*sys%n_acc), intent(in)   :: wf
    real(8),intent(in)                        :: mu(2)
    real(8), dimension(sys%np), intent(in)    :: Aemdr    ! u_Aeff = d(chi) + Aemdr
                                                          !Aemdr = - intgral{A^{em} dr}
    !complex(8), dimension(sys%pl_1st_l(1)%length()), intent(in) :: delta
    complex(8), dimension(sys%pl_1st_l(1)%length()), intent(in) :: tdelta_for, tdelta_rev
    real(8), dimension(sys%n), intent(in)   :: ne_u, ne_d
    real(8), dimension(sys%n), intent(in) :: xi, chi
    !complex(8),DIMENSION(4*sys%n_acc, 4*sys%n_acc), intent(in)   :: hfwf
    !real(8), dimension(sys%n), intent(in) :: xi
    real(8), dimension(sys%n), intent(in) :: Sx, Sy, Sz
    complex(8),DIMENSION(4*sys%n_acc, 4*sys%n_acc)   :: twf
    !real(8), dimension(sys%n) :: ne_u, ne_d
    integer :: i, j, k, is, js, ks, r, m
    real(8) :: t3U
    !complex(8) :: exi,echi, cexi
    complex(8) :: expAem, cexpAem

    !! error cheack
    if(sys%nz /= 2 .or. sys%nh_l(1) /= 0) then
      print *, "nz : ", sys%nz
      print *, "nh(1) : ", sys%nh_l(1)
      stop "error : wrong data in calc_eg_sb_2lay_xi_hamiltonian"
      ! sys subroutine is available for only 2 layer"
    endif

    !twf = wf
    do i = 1, sys%n_acc
      is = sys%ets(i)
      twf(4*i-3, :) = exp(-0.5d0*ui*chi(is))*exp(-0.5d0*ui*xi(is))*wf(4*i-3, :)
      twf(4*i-2, :) = exp(-0.5d0*ui*chi(is))*exp( 0.5d0*ui*xi(is))*wf(4*i-2, :)
      twf(4*i-1, :) = exp( 0.5d0*ui*chi(is))*exp( 0.5d0*ui*xi(is))*wf(4*i-1, :)
      twf(4*i  , :) = exp( 0.5d0*ui*chi(is))*exp(-0.5d0*ui*xi(is))*wf(4*i  , :)
    end do

    t3U = 2*sys%t(3)**2/sys%U

    !ne_u = 0.5d0
    !ne_d = 0.5d0
    allocate(eg(4*sys%n_acc))
    eg = 0d0

    do r = 1, 4*sys%n_acc
      !!!!!!!!!!!!!!!!!!    surface    !!!!!!!!!!!!!!!!!!
      do i = sys%n_acc_start(1), sys%n_acc_end(1)
        eg(r) = eg(r) - mu(1)*conjg(twf(4*i-3, r))*twf(4*i-3, r)
        eg(r) = eg(r) - mu(1)*conjg(twf(4*i-2, r))*twf(4*i-2, r)
        eg(r) = eg(r) + mu(1)*conjg(twf(4*i-1, r))*twf(4*i-1, r)
        eg(r) = eg(r) + mu(1)*conjg(twf(4*i  , r))*twf(4*i  , r)
      end do
      do i = 1, sys%pl_1st_l(1)%length()
        js = sys%pl_1st_l(1)%value(i)%i
        ks = sys%pl_1st_l(1)%value(i)%f
        j = sys%ste(js) ! in the 1st layer (surface), js = j, ks = k
        k = sys%ste(ks)
        m = sys%cpl%seek( sys%pl_1st_l(1)%value(i) )
        expAem = exp(0.5d0*ui*Aemdr(m))
        cexpAem = exp(-0.5d0*ui*Aemdr(m))
        eg(r) = eg(r) - sys%t(1)*(1d0-ne_d(js))*(1d0-ne_d(ks))*cexpAem*conjg(twf(4*j-3, r))*twf(4*k-3, r)
        eg(r) = eg(r) - sys%t(1)*(1d0-ne_u(js))*(1d0-ne_u(ks))*cexpAem*conjg(twf(4*j-2, r))*twf(4*k-2, r)
        eg(r) = eg(r) + sys%t(1)*(1d0-ne_d(js))*(1d0-ne_d(ks))*cexpAem*conjg(twf(4*j-1, r))*twf(4*k-1, r)
        eg(r) = eg(r) + sys%t(1)*(1d0-ne_u(js))*(1d0-ne_u(ks))*cexpAem*conjg(twf(4*j  , r))*twf(4*k  , r)
        eg(r) = eg(r) - sys%t(1)*(1d0-ne_d(js))*(1d0-ne_d(ks))* expAem*conjg(twf(4*k-3, r))*twf(4*j-3, r)
        eg(r) = eg(r) - sys%t(1)*(1d0-ne_u(js))*(1d0-ne_u(ks))* expAem*conjg(twf(4*k-2, r))*twf(4*j-2, r)
        eg(r) = eg(r) + sys%t(1)*(1d0-ne_d(js))*(1d0-ne_d(ks))* expAem*conjg(twf(4*k-1, r))*twf(4*j-1, r)
        eg(r) = eg(r) + sys%t(1)*(1d0-ne_u(js))*(1d0-ne_u(ks))* expAem*conjg(twf(4*k  , r))*twf(4*j  , r)
      
      end do
      do i = 1, sys%pl_2nd_l(1)%length()
        js = sys%pl_2nd_l(1)%value(i)%i
        ks = sys%pl_2nd_l(1)%value(i)%f
        j = sys%ste(js) ! in the 1st layer (surface), js = j, ks = k
        k = sys%ste(ks)
        !m = sys%cpl%seek( sys%pl_2nd_l(1)%value(i) )
        !expAem = exp(0.5d0*ui*Aemdr(m))
        !cexpAem = exp(-0.5d0*ui*Aemdr(m))
        expAem = 1d0
        cexpAem = 1d0
        eg(r) = eg(r) - sys%t(2)*(1d0-ne_d(js))*(1d0-ne_d(ks))*cexpAem*conjg(twf(4*j-3, r))*twf(4*k-3, r)
        eg(r) = eg(r) - sys%t(2)*(1d0-ne_u(js))*(1d0-ne_u(ks))*cexpAem*conjg(twf(4*j-2, r))*twf(4*k-2, r)
        eg(r) = eg(r) + sys%t(2)*(1d0-ne_d(js))*(1d0-ne_d(ks))*cexpAem*conjg(twf(4*j-1, r))*twf(4*k-1, r)
        eg(r) = eg(r) + sys%t(2)*(1d0-ne_u(js))*(1d0-ne_u(ks))*cexpAem*conjg(twf(4*j  , r))*twf(4*k  , r)
        eg(r) = eg(r) - sys%t(2)*(1d0-ne_d(js))*(1d0-ne_d(ks))* expAem*conjg(twf(4*k-3, r))*twf(4*j-3, r)
        eg(r) = eg(r) - sys%t(2)*(1d0-ne_u(js))*(1d0-ne_u(ks))* expAem*conjg(twf(4*k-2, r))*twf(4*j-2, r)
        eg(r) = eg(r) + sys%t(2)*(1d0-ne_d(js))*(1d0-ne_d(ks))* expAem*conjg(twf(4*k-1, r))*twf(4*j-1, r)
        eg(r) = eg(r) + sys%t(2)*(1d0-ne_u(js))*(1d0-ne_u(ks))* expAem*conjg(twf(4*k  , r))*twf(4*j  , r)
      end do
      do i = 1, sys%pl_1st_l(1)%length()
        ! here, delta's index corresponds to pl_1st_l(1)'s index
        js = sys%pl_1st_l(1)%value(i)%i
        ks = sys%pl_1st_l(1)%value(i)%f
        j = sys%ste(js) ! in the 1st layer ( surface ), js = j, ks = k
        k = sys%ste(ks)
        
        !m = sys%cpl%seek( sys%pl_1st_l(1)%value(i) )
        !expAem = exp(0.5d0*ui*Aemdr(m))
        !cexpAem = exp(-0.5d0*ui*Aemdr(m))
      
        !eg(r) = eg(r) +       tdelta_for(i) *cexpAem*conjg(wf(4*j-3, r))*wf(4*k  , r)
        !eg(r) = eg(r) +       tdelta_rev(i) *cexpAem*conjg(wf(4*j-2, r))*wf(4*k-1, r)
        !eg(r) = eg(r) + conjg(tdelta_for(i))*cexpAem*conjg(wf(4*j-1, r))*wf(4*k-2, r)
        !eg(r) = eg(r) + conjg(tdelta_rev(i))*cexpAem*conjg(wf(4*j  , r))*wf(4*k-3, r)
        !eg(r) = eg(r) +       tdelta_rev(i) * expAem*conjg(wf(4*k-3, r))*wf(4*j  , r)
        !eg(r) = eg(r) +       tdelta_for(i) * expAem*conjg(wf(4*k-2, r))*wf(4*j-1, r)
        !eg(r) = eg(r) + conjg(tdelta_rev(i))* expAem*conjg(wf(4*k-1, r))*wf(4*j-2, r)
        !eg(r) = eg(r) + conjg(tdelta_for(i))* expAem*conjg(wf(4*k  , r))*wf(4*j-3, r)
        
        eg(r) = eg(r) +       tdelta_for(i) *conjg(wf(4*j-3, r))*wf(4*k  , r)
        eg(r) = eg(r) +       tdelta_rev(i) *conjg(wf(4*j-2, r))*wf(4*k-1, r)
        eg(r) = eg(r) + conjg(tdelta_for(i))*conjg(wf(4*j-1, r))*wf(4*k-2, r)
        eg(r) = eg(r) + conjg(tdelta_rev(i))*conjg(wf(4*j  , r))*wf(4*k-3, r)
        eg(r) = eg(r) +       tdelta_rev(i) *conjg(wf(4*k-3, r))*wf(4*j  , r)
        eg(r) = eg(r) +       tdelta_for(i) *conjg(wf(4*k-2, r))*wf(4*j-1, r)
        eg(r) = eg(r) + conjg(tdelta_rev(i))*conjg(wf(4*k-1, r))*wf(4*j-2, r)
        eg(r) = eg(r) + conjg(tdelta_for(i))*conjg(wf(4*k  , r))*wf(4*j-3, r)
      end do
      !!!!!!!!!!!!!!    end surface    !!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!    bulk    !!!!!!!!!!!!!!!!
      do i = sys%n_acc_start(2), sys%n_acc_end(2)
        is = sys%ets(i)
        eg(r) = eg(r) + sys%U*(-2d0/3d0*Sz(is)+0.5d0)*conjg(twf(4*i-3, r))*twf(4*i-3, r)
        eg(r) = eg(r) + sys%U*( 2d0/3d0*Sz(is)+0.5d0)*conjg(twf(4*i-2, r))*twf(4*i-2, r)
        eg(r) = eg(r) - sys%U*(-2d0/3d0*Sz(is)+0.5d0)*conjg(twf(4*i-1, r))*twf(4*i-1, r)
        eg(r) = eg(r) - sys%U*( 2d0/3d0*Sz(is)+0.5d0)*conjg(twf(4*i  , r))*twf(4*i  , r)
        eg(r) = eg(r) - sys%U*2d0/3d0*(Sx(is)-ui*Sy(is))*conjg(twf(4*i-3, r))*twf(4*i-2, r)
        eg(r) = eg(r) - sys%U*2d0/3d0*(Sx(is)+ui*Sy(is))*conjg(twf(4*i-2, r))*twf(4*i-3, r)
        eg(r) = eg(r) - sys%U*2d0/3d0*(Sx(is)+ui*Sy(is))*conjg(twf(4*i-1, r))*twf(4*i  , r)
        eg(r) = eg(r) - sys%U*2d0/3d0*(Sx(is)-ui*Sy(is))*conjg(twf(4*i  , r))*twf(4*i-1, r)
      end do
      do i = sys%n_acc_start(2), sys%n_acc_end(2)
        eg(r) = eg(r) - mu(2)*conjg(twf(4*i-3, r))*twf(4*i-3, r)
        eg(r) = eg(r) - mu(2)*conjg(twf(4*i-2, r))*twf(4*i-2, r)
        eg(r) = eg(r) + mu(2)*conjg(twf(4*i-1, r))*twf(4*i-1, r)
        eg(r) = eg(r) + mu(2)*conjg(twf(4*i  , r))*twf(4*i  , r)
      end do
      do i = 1, sys%pl_1st_l(2)%length()
        js = sys%pl_1st_l(2)%value(i)%i
        ks = sys%pl_1st_l(2)%value(i)%f
        j = sys%ste(js)
        k = sys%ste(ks)
        m = sys%cpl%seek( sys%pl_1st_l(2)%value(i) )
        expAem = exp(0.5d0*ui*Aemdr(m))
        cexpAem = exp(-0.5d0*ui*Aemdr(m))
        eg(r) = eg(r) - sys%t(1)*cexpAem*conjg(twf(4*j-3, r))*twf(4*k-3, r)
        eg(r) = eg(r) - sys%t(1)*cexpAem*conjg(twf(4*j-2, r))*twf(4*k-2, r)
        eg(r) = eg(r) + sys%t(1)*cexpAem*conjg(twf(4*j-1, r))*twf(4*k-1, r)
        eg(r) = eg(r) + sys%t(1)*cexpAem*conjg(twf(4*j  , r))*twf(4*k  , r)
        eg(r) = eg(r) - sys%t(1)* expAem*conjg(twf(4*k-3, r))*twf(4*j-3, r)
        eg(r) = eg(r) - sys%t(1)* expAem*conjg(twf(4*k-2, r))*twf(4*j-2, r)
        eg(r) = eg(r) + sys%t(1)* expAem*conjg(twf(4*k-1, r))*twf(4*j-1, r)
        eg(r) = eg(r) + sys%t(1)* expAem*conjg(twf(4*k  , r))*twf(4*j  , r)
      end do
      do i = 1, sys%pl_2nd_l(2)%length()
        js = sys%pl_2nd_l(2)%value(i)%i
        ks = sys%pl_2nd_l(2)%value(i)%f
        j = sys%ste(js)
        k = sys%ste(ks)
        !m = sys%cpl%seek( sys%pl_2nd_l(2)%value(i) )
        !expAem = exp(0.5d0*ui*Aemdr(m))
        !cexpAem = exp(-0.5d0*ui*Aemdr(m))
        expAem = 1d0
        cexpAem = 1d0
        eg(r) = eg(r) - sys%t(2)*cexpAem*conjg(twf(4*j-3, r))*twf(4*k-3, r)
        eg(r) = eg(r) - sys%t(2)*cexpAem*conjg(twf(4*j-2, r))*twf(4*k-2, r)
        eg(r) = eg(r) + sys%t(2)*cexpAem*conjg(twf(4*j-1, r))*twf(4*k-1, r)
        eg(r) = eg(r) + sys%t(2)*cexpAem*conjg(twf(4*j  , r))*twf(4*k  , r)
        eg(r) = eg(r) - sys%t(2)* expAem*conjg(twf(4*k-3, r))*twf(4*j-3, r)
        eg(r) = eg(r) - sys%t(2)* expAem*conjg(twf(4*k-2, r))*twf(4*j-2, r)
        eg(r) = eg(r) + sys%t(2)* expAem*conjg(twf(4*k-1, r))*twf(4*j-1, r)
        eg(r) = eg(r) + sys%t(2)* expAem*conjg(twf(4*k  , r))*twf(4*j  , r)
      end do
      do i = 1, sys%pl_h%length()
        js = sys%pl_h%value(i)%i
        ks = sys%pl_h%value(i)%f
        j = sys%ste(js)
        k = sys%ste(ks)
        eg(r) = eg(r) + 0.5d0*sys%jd*Sz(ks)*conjg(twf(4*j-3, r))*twf(4*j-3, r)
        eg(r) = eg(r) - 0.5d0*sys%jd*Sz(ks)*conjg(twf(4*j-2, r))*twf(4*j-2, r)
        eg(r) = eg(r) - 0.5d0*sys%jd*Sz(ks)*conjg(twf(4*j-1, r))*twf(4*j-1, r)
        eg(r) = eg(r) + 0.5d0*sys%jd*Sz(ks)*conjg(twf(4*j  , r))*twf(4*j  , r)
        eg(r) = eg(r) + 0.5d0*sys%jd*(Sx(ks) - ui*Sy(ks))*conjg(twf(4*j-3, r))*twf(4*j-2, r)
        eg(r) = eg(r) + 0.5d0*sys%jd*(Sx(ks) + ui*Sy(ks))*conjg(twf(4*j-2, r))*twf(4*j-3, r)
        eg(r) = eg(r) + 0.5d0*sys%jd*(Sx(ks) + ui*Sy(ks))*conjg(twf(4*j-1, r))*twf(4*j  , r)
        eg(r) = eg(r) + 0.5d0*sys%jd*(Sx(ks) - ui*Sy(ks))*conjg(twf(4*j  , r))*twf(4*j-1, r)
      
        eg(r) = eg(r) + 0.5d0*sys%jd*Sz(js)*conjg(twf(4*k-3, r))*twf(4*k-3, r)
        eg(r) = eg(r) - 0.5d0*sys%jd*Sz(js)*conjg(twf(4*k-2, r))*twf(4*k-2, r)
        eg(r) = eg(r) - 0.5d0*sys%jd*Sz(js)*conjg(twf(4*k-1, r))*twf(4*k-1, r)
        eg(r) = eg(r) + 0.5d0*sys%jd*Sz(js)*conjg(twf(4*k  , r))*twf(4*k  , r)
        eg(r) = eg(r) + 0.5d0*sys%jd*(Sx(js) - ui*Sy(js))*conjg(twf(4*k-3, r))*twf(4*k-2, r)
        eg(r) = eg(r) + 0.5d0*sys%jd*(Sx(js) + ui*Sy(js))*conjg(twf(4*k-2, r))*twf(4*k-3, r)
        eg(r) = eg(r) + 0.5d0*sys%jd*(Sx(js) + ui*Sy(js))*conjg(twf(4*k-1, r))*twf(4*k  , r)
        eg(r) = eg(r) + 0.5d0*sys%jd*(Sx(js) - ui*Sy(js))*conjg(twf(4*k  , r))*twf(4*k-1, r)
      end do
    
      do i = 1, sys%pl_rh_1%length()  ! rashba (upper-right path)
        js = sys%pl_rh_1%value(i)%i ! h-y   |   h-x
        ks = sys%pl_rh_1%value(i)%f ! h+x   |   h+y
        j = sys%ste(js)
        k = sys%ste(ks)
        m = sys%cpl%seek( sys%pl_rh_1%value(i) )
        expAem = exp(0.5d0*ui*Aemdr(m))
        cexpAem = exp(-0.5d0*ui*Aemdr(m))
        eg(r) = eg(r) + sys%lm*cexpAem*conjg(twf(4*j-2, r))*twf(4*k-3, r)
        eg(r) = eg(r) - sys%lm*cexpAem*conjg(twf(4*j-3, r))*twf(4*k-2, r)
        eg(r) = eg(r) + sys%lm*cexpAem*conjg(twf(4*j  , r))*twf(4*k-1, r)
        eg(r) = eg(r) - sys%lm*cexpAem*conjg(twf(4*j-1, r))*twf(4*k  , r)
        eg(r) = eg(r) + sys%lm* expAem*conjg(twf(4*k-2, r))*twf(4*j-3, r)
        eg(r) = eg(r) - sys%lm* expAem*conjg(twf(4*k-3, r))*twf(4*j-2, r)
        eg(r) = eg(r) + sys%lm* expAem*conjg(twf(4*k  , r))*twf(4*j-1, r)
        eg(r) = eg(r) - sys%lm* expAem*conjg(twf(4*k-1, r))*twf(4*j  , r)
      end do
      do i = 1, sys%pl_rh_2%length()  ! rashba (upper-left path)
        js = sys%pl_rh_2%value(i)%i ! h-y   |   h+x
        ks = sys%pl_rh_2%value(i)%f ! h-x   |   h+y
        j = sys%ste(js)
        k = sys%ste(ks)
        m = sys%cpl%seek( sys%pl_rh_2%value(i) )
        expAem = exp(0.5d0*ui*Aemdr(m))
        cexpAem = exp(-0.5d0*ui*Aemdr(m))
        eg(r) = eg(r) + ui*sys%lm*cexpAem*conjg(twf(4*j-2, r))*twf(4*k-3, r)
        eg(r) = eg(r) + ui*sys%lm*cexpAem*conjg(twf(4*j-3, r))*twf(4*k-2, r)
        eg(r) = eg(r) - ui*sys%lm*cexpAem*conjg(twf(4*j  , r))*twf(4*k-1, r)
        eg(r) = eg(r) - ui*sys%lm*cexpAem*conjg(twf(4*j-1, r))*twf(4*k  , r)
        eg(r) = eg(r) - ui*sys%lm* expAem*conjg(twf(4*k-2, r))*twf(4*j-3, r)
        eg(r) = eg(r) - ui*sys%lm* expAem*conjg(twf(4*k-3, r))*twf(4*j-2, r)
        eg(r) = eg(r) + ui*sys%lm* expAem*conjg(twf(4*k  , r))*twf(4*j-1, r)
        eg(r) = eg(r) + ui*sys%lm* expAem*conjg(twf(4*k-1, r))*twf(4*j  , r)
      end do
      !!!!!!!!!!!!!!    end bulk    !!!!!!!!!!!!!!!!!!
  !
      !!!!!!!!!!!!!!!!!!    interlayer    !!!!!!!!!!!!!!!!
      ! pl_z_l(i) is z-dir path list between ith layer and (i+1)th layer.
      ! here, pl_z_l(1) is interlayer path ( between 1st(surface) and 2nd(bulk) ).
      do i = 1, sys%pl_z_l(1)%length()
        js = sys%pl_z_l(1)%value(i)%i
        ks = sys%pl_z_l(1)%value(i)%f
        j = sys%ste(js)
        k = sys%ste(ks)
        m = sys%cpl%seek( sys%pl_z_l(1)%value(i) )
        expAem = exp(0.5d0*ui*Aemdr(m))
        cexpAem = exp(-0.5d0*ui*Aemdr(m))
        eg(r) = eg(r) - sys%t(3)*(1d0-ne_d(js))*(1d0-ne_d(ks))*cexpAem*conjg(twf(4*j-3, r))*twf(4*k-3, r)
        eg(r) = eg(r) - sys%t(3)*(1d0-ne_u(js))*(1d0-ne_u(ks))*cexpAem*conjg(twf(4*j-2, r))*twf(4*k-2, r)
        eg(r) = eg(r) + sys%t(3)*(1d0-ne_d(js))*(1d0-ne_d(ks))*cexpAem*conjg(twf(4*j-1, r))*twf(4*k-1, r)
        eg(r) = eg(r) + sys%t(3)*(1d0-ne_u(js))*(1d0-ne_u(ks))*cexpAem*conjg(twf(4*j  , r))*twf(4*k  , r)
        eg(r) = eg(r) - sys%t(3)*(1d0-ne_d(js))*(1d0-ne_d(ks))* expAem*conjg(twf(4*k-3, r))*twf(4*j-3, r)
        eg(r) = eg(r) - sys%t(3)*(1d0-ne_u(js))*(1d0-ne_u(ks))* expAem*conjg(twf(4*k-2, r))*twf(4*j-2, r)
        eg(r) = eg(r) + sys%t(3)*(1d0-ne_d(js))*(1d0-ne_d(ks))* expAem*conjg(twf(4*k-1, r))*twf(4*j-1, r)
        eg(r) = eg(r) + sys%t(3)*(1d0-ne_u(js))*(1d0-ne_u(ks))* expAem*conjg(twf(4*k  , r))*twf(4*j  , r)
      
        eg(r) = eg(r) + t3U*(Sz(ks)-0.5d0*(ne_u(ks)+ne_d(ks)))*conjg(twf(4*j-3, r))*twf(4*j-3, r)
        eg(r) = eg(r) - t3U*(Sz(ks)+0.5d0*(ne_u(ks)+ne_d(ks)))*conjg(twf(4*j-2, r))*twf(4*j-2, r)
        eg(r) = eg(r) - t3U*(Sz(ks)-0.5d0*(ne_u(ks)+ne_d(ks)))*conjg(twf(4*j-1, r))*twf(4*j-1, r)
        eg(r) = eg(r) + t3U*(Sz(ks)+0.5d0*(ne_u(ks)+ne_d(ks)))*conjg(twf(4*j  , r))*twf(4*j  , r)
        eg(r) = eg(r) + t3U*(Sx(ks) - ui*Sy(ks))*conjg(twf(4*j-3, r))*twf(4*j-2, r)
        eg(r) = eg(r) + t3U*(Sx(ks) + ui*Sy(ks))*conjg(twf(4*j-2, r))*twf(4*j-3, r)
        eg(r) = eg(r) + t3U*(Sx(ks) + ui*Sy(ks))*conjg(twf(4*j-1, r))*twf(4*j  , r)
        eg(r) = eg(r) + t3U*(Sx(ks) - ui*Sy(ks))*conjg(twf(4*j  , r))*twf(4*j-1, r)
      
        eg(r) = eg(r) + t3U*(Sz(js)-0.5d0*(ne_u(js)+ne_d(js)))*conjg(twf(4*k-3, r))*twf(4*k-3, r)
        eg(r) = eg(r) - t3U*(Sz(js)+0.5d0*(ne_u(js)+ne_d(js)))*conjg(twf(4*k-2, r))*twf(4*k-2, r)
        eg(r) = eg(r) - t3U*(Sz(js)-0.5d0*(ne_u(js)+ne_d(js)))*conjg(twf(4*k-1, r))*twf(4*k-1, r)
        eg(r) = eg(r) + t3U*(Sz(js)+0.5d0*(ne_u(js)+ne_d(js)))*conjg(twf(4*k  , r))*twf(4*k  , r)
        eg(r) = eg(r) + t3U*(Sx(js) - ui*Sy(js))*conjg(twf(4*k-3, r))*twf(4*k-2, r)
        eg(r) = eg(r) + t3U*(Sx(js) + ui*Sy(js))*conjg(twf(4*k-2, r))*twf(4*k-3, r)
        eg(r) = eg(r) + t3U*(Sx(js) + ui*Sy(js))*conjg(twf(4*k-1, r))*twf(4*k  , r)
        eg(r) = eg(r) + t3U*(Sx(js) - ui*Sy(js))*conjg(twf(4*k  , r))*twf(4*k-1, r)
      end do
      !!!!!!!!!!!!!!   end  interlayer    !!!!!!!!!!!!!!!!
  !  
  !    call my_zheev("l", ham, sys%eg, sys%wf)

    end do
  end subroutine calc_eigen_energy_sb_2lay_Aem
  
  !function total_energy_sb_2lay_xichi_T0_Aem(sys, wf, xi, chi, Aemdr) result(r)
  function total_energy_sb_2lay_xichi_T0_Aem(sys, wf, xi, chi, Aemdr, &
      & oHt_s, oHt2_s, oHt_b, oHt2_b, oHrsoc, oH_interlayer ) result(r)
    implicit none
    class(base), intent(in)   :: sys
    complex(8),dimension(4*sys%n_acc, 4*sys%n_acc), intent(in)   :: wf
    !real(8),intent(in)    :: mu(2)
    !complex(8), dimension(sys%pl_1st_l(1)%length()), intent(in) :: delta
    !real(8), dimension(sys%n), intent(in)   :: ne_u, ne_d
    real(8), dimension(sys%n)             :: ne_u, ne_d
    real(8), dimension(sys%n), intent(in) :: xi, chi
    real(8), dimension(sys%np), intent(in)    :: Aemdr    ! u_Aeff = d(chi) + Aemdr
                                                          !Aemdr = - intgral{A^{em} dr}
    !real(8), dimension(sys%n), intent(in) :: Sx, Sy, Sz, xi, chi
    complex(8),DIMENSION(4*sys%n_acc, 4*sys%n_acc)   :: twf
    !real(8), dimension(:), allocatable :: eigen_eg
    real(8) :: r
    real(8) :: Ht_s, Ht2_s, Ht_b, Ht2_b, Hrsoc, H_interlayer
    real(8), intent(out), optional :: oHt_s, oHt2_s, oHt_b, oHt2_b, oHrsoc, oH_interlayer
    integer :: i, is, j, js, k, ks, a, m
    complex(8) :: expAem, cexpAem
    r = 0d0
    twf = wf
    do i = 1, sys%n_acc
      is = sys%ets(i)
      twf(4*i-3, :) = exp(-0.5d0*ui*chi(is))*exp(-0.5d0*ui*xi(is))*twf(4*i-3, :)
      twf(4*i-2, :) = exp(-0.5d0*ui*chi(is))*exp( 0.5d0*ui*xi(is))*twf(4*i-2, :)
      twf(4*i-1, :) = exp( 0.5d0*ui*chi(is))*exp( 0.5d0*ui*xi(is))*twf(4*i-1, :)
      twf(4*i  , :) = exp( 0.5d0*ui*chi(is))*exp(-0.5d0*ui*xi(is))*twf(4*i  , :)
    end do    
    
    ne_u = 0d0
    ne_d = 0d0
    do i = 1, sys%n_acc
      is = sys%ets(i)
      do j = 2*sys%n_acc+1, 4*sys%n_acc
        ne_u(is) = ne_u(is) + conjg(twf(4*i-1, j))*twf(4*i-1, j)
        ne_d(is) = ne_d(is) + conjg(twf(4*i  , j))*twf(4*i  , j)
      end do
    end do
    !print *, ne_u(:)
    !stop    
   
    ! surface, the nearest hopping energy
    Ht_s = 0d0
    do i = 1, sys%pl_1st_l(1)%length()
      js = sys%pl_1st_l(1)%value(i)%i
      ks = sys%pl_1st_l(1)%value(i)%f
      j = sys%ste(js) ! in the 1st layer (surface), js = j, ks = k
      k = sys%ste(ks)
      m = sys%cpl%seek( sys%pl_1st_l(1)%value(i) )
      expAem = exp(0.5d0*ui*Aemdr(m))
      cexpAem = exp(-0.5d0*ui*Aemdr(m))
      do a = 2*sys%n_acc+1, 4*sys%n_acc
        Ht_s = Ht_s - sys%t(1)*(1d0-ne_d(js))*(1d0-ne_d(ks))* expAem*twf(4*j-1, a)*conjg(twf(4*k-1, a))
        Ht_s = Ht_s - sys%t(1)*(1d0-ne_u(js))*(1d0-ne_u(ks))* expAem*twf(4*j  , a)*conjg(twf(4*k  , a))
        Ht_s = Ht_s - sys%t(1)*(1d0-ne_d(js))*(1d0-ne_d(ks))*cexpAem*twf(4*k-1, a)*conjg(twf(4*j-1, a))
        Ht_s = Ht_s - sys%t(1)*(1d0-ne_u(js))*(1d0-ne_u(ks))*cexpAem*twf(4*k  , a)*conjg(twf(4*j  , a))
      end do
    end do

    ! surface, the 2nd nearest hopping energy
    Ht2_s = 0d0
    do i = 1, sys%pl_2nd_l(1)%length()
      js = sys%pl_2nd_l(1)%value(i)%i
      ks = sys%pl_2nd_l(1)%value(i)%f
      j = sys%ste(js) ! in the 1st layer (surface), js = j, ks = k
      k = sys%ste(ks)
      ! it need to fix this. Aemdr (made with cpl) does not contain 2nd hop
      !m = sys%cpl%seek( sys%pl_2nd_l(1)%value(i) )
      !expAem = exp(0.5d0*ui*Aemdr(m))
      !cexpAem = exp(-0.5d0*ui*Aemdr(m))
      expAem = 1d0
      cexpAem = 1d0
      do a = 2*sys%n_acc+1, 4*sys%n_acc
        Ht2_s = Ht2_s - sys%t(2)*(1d0-ne_d(js))*(1d0-ne_d(ks))* expAem*twf(4*j-1, a)*conjg(twf(4*k-1, a))
        Ht2_s = Ht2_s - sys%t(2)*(1d0-ne_u(js))*(1d0-ne_u(ks))* expAem*twf(4*j  , a)*conjg(twf(4*k  , a))
        Ht2_s = Ht2_s - sys%t(2)*(1d0-ne_d(js))*(1d0-ne_d(ks))*cexpAem*twf(4*k-1, a)*conjg(twf(4*j-1, a))
        Ht2_s = Ht2_s - sys%t(2)*(1d0-ne_u(js))*(1d0-ne_u(ks))*cexpAem*twf(4*k  , a)*conjg(twf(4*j  , a))
      end do
    end do
    
    ! bulk, the nearest hopping energy
    Ht_b = 0d0
    do i = 1, sys%pl_1st_l(2)%length()
      js = sys%pl_1st_l(2)%value(i)%i
      ks = sys%pl_1st_l(2)%value(i)%f
      j = sys%ste(js) ! in the 1st layer (surface), js = j, ks = k
      k = sys%ste(ks)
      m = sys%cpl%seek( sys%pl_1st_l(2)%value(i) )
      !print *, m
      expAem = exp(0.5d0*ui*Aemdr(m))
      cexpAem = exp(-0.5d0*ui*Aemdr(m))
      do a = 2*sys%n_acc+1, 4*sys%n_acc
        Ht_b = Ht_b - sys%t(1)* expAem*twf(4*j-1, a)*conjg(twf(4*k-1, a))
        Ht_b = Ht_b - sys%t(1)* expAem*twf(4*j  , a)*conjg(twf(4*k  , a))
        Ht_b = Ht_b - sys%t(1)*cexpAem*twf(4*k-1, a)*conjg(twf(4*j-1, a))
        Ht_b = Ht_b - sys%t(1)*cexpAem*twf(4*k  , a)*conjg(twf(4*j  , a))
      end do
    end do
    
    ! bulk, the 2nd nearest hopping energy
    Ht2_b = 0d0
    do i = 1, sys%pl_2nd_l(2)%length()
      js = sys%pl_2nd_l(2)%value(i)%i
      ks = sys%pl_2nd_l(2)%value(i)%f
      j = sys%ste(js) ! in the 1st layer (surface), js = j, ks = k
      k = sys%ste(ks)
      
      ! it need to fix this. Aemdr (made with cpl) does not contain 2nd hop
      !m = sys%cpl%seek( sys%pl_2nd_l(2)%value(i) )
      !expAem = exp(0.5d0*ui*Aemdr(m))
      !cexpAem = exp(-0.5d0*ui*Aemdr(m))
      expAem = 1d0
      cexpAem = 1d0
      do a = 2*sys%n_acc+1, 4*sys%n_acc
        Ht2_b = Ht2_b - sys%t(2)* expAem*twf(4*j-1, a)*conjg(twf(4*k-1, a))
        Ht2_b = Ht2_b - sys%t(2)* expAem*twf(4*j  , a)*conjg(twf(4*k  , a))
        Ht2_b = Ht2_b - sys%t(2)*cexpAem*twf(4*k-1, a)*conjg(twf(4*j-1, a))
        Ht2_b = Ht2_b - sys%t(2)*cexpAem*twf(4*k  , a)*conjg(twf(4*j  , a))
      end do
    end do

    ! bulk, rashba energy
    Hrsoc = 0d0
    do i = 1, sys%pl_rh_1%length()
      js = sys%pl_rh_1%value(i)%i
      ks = sys%pl_rh_1%value(i)%f
      j = sys%ste(js)
      k = sys%ste(ks)
      m = sys%cpl%seek( sys%pl_rh_1%value(i) )
      expAem = exp(0.5d0*ui*Aemdr(m))
      cexpAem = exp(-0.5d0*ui*Aemdr(m))
      do a = 2*sys%n_acc+1, 4*sys%n_acc
        Hrsoc = Hrsoc - sys%lm*cexpAem*twf(4*k  , a)*conjg(twf(4*j-1, a))
        Hrsoc = Hrsoc + sys%lm*cexpAem*twf(4*k-1, a)*conjg(twf(4*j  , a))
        Hrsoc = Hrsoc - sys%lm* expAem*twf(4*j-1, a)*conjg(twf(4*k  , a))
        Hrsoc = Hrsoc + sys%lm* expAem*twf(4*j  , a)*conjg(twf(4*k-1, a))
      end do
    end do
    do i = 1, sys%pl_rh_2%length() 
      js = sys%pl_rh_2%value(i)%i
      ks = sys%pl_rh_2%value(i)%f
      j = sys%ste(js)
      k = sys%ste(ks)
      m = sys%cpl%seek( sys%pl_rh_2%value(i) )
      expAem = exp(0.5d0*ui*Aemdr(m))
      cexpAem = exp(-0.5d0*ui*Aemdr(m))
      do a = 2*sys%n_acc+1, 4*sys%n_acc
        Hrsoc = Hrsoc - ui*sys%lm*cexpAem*twf(4*k  , a)*conjg(twf(4*j-1, a))
        Hrsoc = Hrsoc - ui*sys%lm*cexpAem*twf(4*k-1, a)*conjg(twf(4*j  , a))
        Hrsoc = Hrsoc + ui*sys%lm* expAem*twf(4*j-1, a)*conjg(twf(4*k  , a))
        Hrsoc = Hrsoc + ui*sys%lm* expAem*twf(4*j  , a)*conjg(twf(4*k-1, a))
      end do
    end do
    
    ! inter-layer, hopping energy
    H_interlayer = 0d0
    do i = 1, sys%pl_z_l(1)%length()
      js = sys%pl_z_l(1)%value(i)%i
      ks = sys%pl_z_l(1)%value(i)%f
      j = sys%ste(js)
      k = sys%ste(ks)
      m = sys%cpl%seek( sys%pl_z_l(1)%value(i) )
      expAem = exp(0.5d0*ui*Aemdr(m))
      cexpAem = exp(-0.5d0*ui*Aemdr(m))
      do a = 2*sys%n_acc+1, 4*sys%n_acc
        H_interlayer = H_interlayer - sys%t(3)*(1d0-ne_d(js))*(1d0-ne_d(ks))*cexpAem*conjg(twf(4*j-1, a))*twf(4*k-1, a)
        H_interlayer = H_interlayer - sys%t(3)*(1d0-ne_u(js))*(1d0-ne_u(ks))*cexpAem*conjg(twf(4*j  , a))*twf(4*k  , a)
        H_interlayer = H_interlayer - sys%t(3)*(1d0-ne_d(js))*(1d0-ne_d(ks))*expAem*conjg(twf(4*k-1, a))*twf(4*j-1, a)
        H_interlayer = H_interlayer - sys%t(3)*(1d0-ne_u(js))*(1d0-ne_u(ks))*expAem*conjg(twf(4*k  , a))*twf(4*j  , a)
      end do
    end do

    r = Ht_s + Ht2_s + Ht_b + Ht2_b + Hrsoc + H_interlayer
    
    if(present(oHt_s)) oHt_s = Ht_s
    if(present(oHt2_s)) oHt2_s = Ht2_s
    if(present(oHt_b)) oHt_b = Ht_b
    if(present(oHt2_b)) oHt2_b = Ht2_b
    if(present(oHrsoc)) oHrsoc = Hrsoc
    if(present(oH_interlayer)) oH_interlayer = H_interlayer

  end function total_energy_sb_2lay_xichi_T0_Aem

END MODULE energy_mod
