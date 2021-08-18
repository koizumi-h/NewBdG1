program stm_point


  !!!!!! modules


    use math_mod,         only : ui, pi, my_zheev, fermi, dfermi
    !use path_list_mod,    only : path_list
    !use hop_iterator_mod, only : get_hopping_path
    use io_mod,           only : save_txt, load_txt, create_vortex, load_conf_3D_2 ,load_conf_3D, load_conf_chi_w
    use diagonalize_mod,  only : diag, new_diag
    use parameter_mod,    only : site_max, site_to_r, xyz_to_site


  !!!!!! end modules


  !!!!!! definitions


    implicit none

    type(diag), allocatable                 :: dg

    complex(8), dimension(:,:), allocatable :: wf

    real(8), dimension(:), allocatable      :: eg
    real(8), dimension(:), allocatable      :: ne
    real(8), dimension(3), parameter        :: t = [1d0, -0.12d0, 0.01d0]  ! hopping coeff[1st, 2nd, layer]
    !real(8), dimension(3), parameter        :: t = [1d0, 0d0, 0.1d0]  ! hopping coeff[1st, 2nd, layer]
    real(8), parameter                      :: u = 8.d0 ! coulomb coeff
    real(8), parameter                      :: lm = 0.02d0 ! rashba coeff
    real(8), parameter                      :: jd = 0.5d0*4.0d0*(t(1)**2)/u  ! across hole coeff
    !real(8), parameter                      :: jd = 0.0d0  ! across hole coeff
    !real(8), parameter                      :: kbt = 0.01d0 ! temperature
    real(8)                                 :: kbt = 0.01d0 ! temperature
    !real(8), parameter                      :: set_mu1 = -0.19d0
    real(8)                                 :: e, kbt_rho(10), dos(10), rho, eg_range(2)

    integer, dimension(:,:), allocatable    :: pshape
    integer, dimension(:,:), allocatable    :: hole_xi, hole_chi  ! dim is (holes,2)
    integer, parameter                      :: unit_num = 10
    !integer, parameter                      :: z_calc = 1 ! 1 => first layer, 2 => second layer
    !integer, parameter                      :: fchi_select = 2 ! 0 => sc, 1 => fchi, 2 => mwf2
    !integer, parameter                      :: i_repeat = 1000
    integer                                 :: num, snum, fnum
    integer                                 :: nx, ny, nz, n, n_acc, nh, r(3), time_s, time_e, countpersec
    integer                                 :: n_step, site
    integer                                 :: i, ix, iy, iz, j, jx, jy, k
    integer                                 :: z_calc, data_select

    character(128)                          :: filename_stm
    character(128)                          :: cat_num
    character(128)                          :: cat_z_calc
    character(128)                          :: cat_gnuplot
    character(128)                          :: cat_fort

    !!! input !!!

      !write(6,*) "input snum"
      !read (5,*) snum
      !write(6,*) "input fchi_select : 0 => sc, 1 => fchi, 2 => mwf2"
      !read (5,*) fchi_select
      write(6,*) "input data_select : 1 => mwf1, 2 => single-shot, 3 => mwf2, 4 => mwf2(no chi)"
      read (5,*) data_select
     ! write(6,*) "input z_calc : 1 => first layer, 2 => second layer"
     ! read (5,*) z_calc
     ! write(6,*) "input e"
     ! read (5,*) e
      write(6,*) "input ix iy, iz"
      read (5,*) ix, iy, iz
      write(6,*) "input kbt"
      read (5,*) kbt

    !!! end input !!!

    !!! vconf !!!

      call load_conf_3D_2('./config/vconf',unit_num,pshape,hole_xi,hole_chi,ne)
      nx    = pshape(1,1)
      ny    = pshape(2,1)
      nz    = size(pshape(1,:))
      nh    = size(hole_xi(:,1))
      n     = site_max(pshape)
      n_acc = n - nh
      print *, nx, ny, nz, nh, n, n_acc
      !stop

    !!! end vconf !!!

    !!! allocate !!!

      allocate(dg)
      dg = new_diag(pshape,hole_xi,t,u,jd,lm)
      allocate(wf(4*n_acc,4*n_acc))
      allocate(eg(4*n_acc))

    !!! end allocate !!!

      !if(fchi_select == 1) then
      !  read(100) eg
      !  read(101) wf
      !else
      !  read(110) eg
      !  read(111) wf
      !end if
      select case(data_select)
        case(1)
          read(100) eg
          read(101) wf
        case(2)
          read(110) eg
          read(111) wf
        case(3)
          read(120) eg
          read(121) wf
        case(4)
          read(120) eg
          read(122) wf
        case default
          print *, "ERROR : data_select ", data_select
          stop
      end select
      
  n_step = 5000
  eg_range(1) = -4d0
  eg_range(2) =  4d0
  !wf = dg%wf
  !eg = dg%eg
  !do i = 1, dg%n_acc
  !  j = dg%ets(i)
  !  wf(4*i-3, :) = exp(-0.5d0*ui*chi(j))*exp( 0.5d0*ui*xi(j))*wf(4*i-3, :)
  !  wf(4*i-2, :) = exp(-0.5d0*ui*chi(j))*exp(-0.5d0*ui*xi(j))*wf(4*i-2, :)
  !  wf(4*i-1, :) = exp( 0.5d0*ui*chi(j))*exp(-0.5d0*ui*xi(j))*wf(4*i-1, :)
  !  wf(4*i  , :) = exp( 0.5d0*ui*chi(j))*exp( 0.5d0*ui*xi(j))*wf(4*i  , :)
  !end do
  open(102,file = "rho.dat")
  do i = 1, n_step
    e = eg_range(1) + (eg_range(2)-eg_range(1))/float(n_step)*i
    rho = 0d0
    do j = 2*dg%n_acc+1,4*dg%n_acc
      site = dg%ste(xyz_to_site(ix,iy,iz,pshape))
      do k = 1, 10
        rho = rho - conjg(wf(4*site-3, j))*wf(4*site-3, j)*dfermi(eg(j) - e, kbt)
        rho = rho - conjg(wf(4*site-2, j))*wf(4*site-2, j)*dfermi(eg(j) - e, kbt)
        rho = rho - conjg(wf(4*site-1, j))*wf(4*site-1, j)*dfermi(eg(j) + e, kbt)
        rho = rho - conjg(wf(4*site  , j))*wf(4*site  , j)*dfermi(eg(j) + e, kbt)
      end do
    end do
    write(102,"(2f20.15)") e, rho
  end do
  close(102)
  print *,"saved to ""rho.dat"""
     

  !!!!!! end calc !!!!!!


end program stm_point