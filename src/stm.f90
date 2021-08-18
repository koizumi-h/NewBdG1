program arhooooooooooooo


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
    real(8), parameter                      :: kbt = 0.01d0 ! temperature
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
    integer                                 :: i, ix, iy, j, jx, jy, k
    integer                                 :: z_calc, data_select

    character(128)                          :: filename_stm
    character(128)                          :: cat_num
    character(128)                          :: cat_z_calc
    character(128)                          :: cat_gnuplot
    character(128)                          :: cat_fort


  !!!!!! end definitions


  !!!!!! preferences


    !!! input !!!

      !write(6,*) "input snum"
      !read (5,*) snum
      !write(6,*) "input fchi_select : 0 => sc, 1 => fchi, 2 => mwf2"
      !read (5,*) fchi_select
      write(6,*) "input data_select : 1 => mwf1, 2 => single-shot, 3 => mwf2, 4 => mwf2(no chi)"
      read (5,*) data_select
      write(6,*) "input z_calc : 1 => first layer, 2 => second layer"
      read (5,*) z_calc
      write(6,*) "input e"
      read (5,*) e

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

    !!! num !!!

      !num = 100*snum
      !print *, "num = ", num
      !!stop
      !if    (fchi_select==0) then
      !  fnum = 0
      !elseif(fchi_select==1) then
      !  fnum = 50
      !elseif(fchi_select==2) then
      !  fnum = 60
      !end if
      !print *, "fnum = ", fnum
      !stop

    !!! emd num !!!

    !!! eg&wf !!!

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
      !read(100) eg
      !read(101) wf
      !write(cat_num,*) num+fnum
      !cat_fort = 'fort/'//"fort."//trim(adjustl(cat_num))
      !open(num+fnum,file=cat_fort,form='unformatted')
      !  read(num+fnum) eg
      !close(num+fnum)
      !write(cat_num,*) num+fnum+1
      !cat_fort = 'fort/'//"fort."//trim(adjustl(cat_num))
      !open(num+fnum+1,file=cat_fort,form='unformatted')
      !  read(num+fnum+1) wf
      !close(num+fnum+1)
      !do i = 1,4*n_acc
      !  print *, i, eg(i)
      !end do
      !stop

    !!! end eg&wf !!!

    !!! filename !!!

      !write(cat_num,   *) num+fnum
      !write(cat_z_calc,*) z_calc
      !filename_stm = "stm/stm_"                // &
      !&              trim(adjustl(cat_num))    // &
      !&              "_("                      // &
      !&              trim(adjustl(cat_z_calc)) // &
      !&              ")"                       // &
      !&              ".txt"
      !print *, "filename = ", trim(adjustl(filename_stm))
      !stop

    !!! end filename !!!

    !!! gnuplot !!!
      !cat_gnuplot = "sp "//""""//trim(adjustl(filename_stm))//""""
      !open(20,file="stm.gpi")
      !  write(20,*) "reset"
      !  write(20,*) "set pm3d map"
      !  write(20,*) "set size ratio 1"
      !  write(20,*) trim(adjustl(cat_gnuplot))
      !close(20)
    !!! end gnuplot !!!


  !!!!!! end preferences


  !!!!!! calc !!!!!!


    open(25,file="stm.txt")

      do iy = 1, ny
        do ix = 1, nx

          if(z_calc==1) then
            !site = real(nx)/2+nx*(real(ny)/2-1)
            site = dg%ste(xyz_to_site(ix,iy,1, pshape))
          else if(z_calc==2) then
            !site = dg%ste(real(nx)/2+nx*(real(ny)/2-1)+nx*ny)
            site = dg%ste(xyz_to_site(ix,iy,2,pshape))
          end if
          !print *, site
          !stop
          rho = 0d0

          if(site /= 0) then
            do j = 2*n_acc+1, 4*n_acc

              rho = rho-conjg(wf(4*site-3,j))*wf(4*site-3,j)*dfermi(eg(j)-e,kbt)
              rho = rho-conjg(wf(4*site-2,j))*wf(4*site-2,j)*dfermi(eg(j)-e,kbt)
              rho = rho-conjg(wf(4*site-1,j))*wf(4*site-1,j)*dfermi(eg(j)+e,kbt)
              rho = rho-conjg(wf(4*site  ,j))*wf(4*site  ,j)*dfermi(eg(j)+e,kbt)

            end do
          end if

          write(25,*) ix, iy, rho

        end do

        write(25,*) ""

      end do

    close(25)

    print *, "saved to ""stm.txt"""


  !!!!!! end calc !!!!!!


end program arhooooooooooooo