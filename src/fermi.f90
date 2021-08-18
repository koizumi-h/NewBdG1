program fermiiiiiiiiiiii


  !!!!!! modules !!!!!!


    use math_mod,         only : ui, pi, my_zheev, fermi, dfermi
    !use path_list_mod,    only : path_list
    !use hop_iterator_mod, only : get_hopping_path
    use io_mod,           only : save_txt, load_txt, create_vortex, load_conf_3D_2 ,load_conf_3D, load_conf_chi_w
    use diagonalize_mod,  only : diag, new_diag
    use parameter_mod,    only : site_max, site_to_r


  !!!!!! end modules !!!!!!


  !!!!!! definitions !!!!!!


    implicit none

    type(diag), allocatable                 :: dg

    complex(8), dimension(:,:), allocatable :: phi1, phi2, phi3, phi4
    complex(8), dimension(:,:), allocatable :: wf

    real(8), dimension(:), allocatable      :: eg
    real(8), dimension(:), allocatable      :: kx
    real(8), dimension(:), allocatable      :: ky
    real(8), dimension(:), allocatable      :: ne
    real(8), dimension(3), parameter        :: t = [1d0, -0.12d0, 0.01d0] ! hopping coeff[1st, 2nd, layer]
    !real(8), dimension(3), parameter        :: t = [1d0, 0d0, 0.1d0] ! hopping coeff[1st, 2nd, layer]
    real(8), parameter                      :: u = 8.d0 ! coulomb coeff
    real(8), parameter                      :: lm = 0.02d0 ! rashba coeff
    real(8), parameter                      :: jd = 0.5d0*4.0d0*(t(1)**2)/u ! across hole coeff
    !real(8), parameter                      :: jd = 0.0d0 ! across hole coeff
    !real(8), parameter                      :: kbt = 0.01d0 ! temperature
    !real(8), parameter                      :: set_mu1 = -0.19d0

    integer, dimension(:,:), allocatable    :: pshape
    integer, dimension(:,:), allocatable    :: hole_xi, hole_chi ! dim is (holes,2)
    integer, parameter                      :: unit_num = 10
    integer, parameter                      :: mesh = 100
    !integer, parameter                      :: z_calc = 1 ! 1 => first layer, 2 => second layer
    !integer, parameter                      :: k_range = 1 ! 0 => pi0, 1 => pipi
    !integer, parameter                      :: fchi_select = 2 ! 0 => sc, 1 => fchi, 2 => mwf2
    !integer, parameter                      :: i_repeat = 1000
    integer                                 :: num, snum, fnum
    integer                                 :: nx, ny, nz, n, n_acc, nh, r(3), time_s, time_e, countpersec
    integer                                 :: vec(3), rvec(3)
    integer                                 :: i, j, jx, jy, p
    integer                                 :: j_rangem, j_rangep
    integer                                 :: leig
    integer                                 :: z_calc, data_select, k_range

    character(128)                          :: filename_fermi_3p4
    character(128)                          :: cat_num
    character(128)                          :: cat_z_calc
    character(128)                          :: cat_k_range
    character(128)                          :: cat_gnuplot
    character(128)                          :: cat_fort


  !!!!!! end definitions !!!!!!


  !!!!!! preferences !!!!!!


    !!! input !!!

      !write(6,*) "input snum"
      !read (5,*) snum
      !write(6,*) "input fchi_select : 0 => sc, 1 => fchi, 2 => mwf2"
      write(6,*) "input data_select : 1 => mwf1, 2 => single-shot, 3 => mwf2, 4 => mwf2(no chi)"
      read (5,*) data_select
      write(6,*) "input z_calc : 1 => first layer, 2 => second layer"
      read (5,*) z_calc
      write(6,*) "input k_range : 0 => pi0, 1 => pipi"
      read (5,*) k_range

    !!! end input !!!

    !!! vconf !!!

      call load_conf_3D_2('./config/vconf',unit_num,pshape,hole_xi,hole_chi,ne)
      nx    = pshape(1,1)
      ny    = pshape(2,1)
      nz    = size(pshape(1,:))
      nh    = size(hole_xi(:,1))
      n     = site_max(pshape)
      n_acc = n-nh
      print *, nx, ny, nz, nh, n, n_acc
      !stop

    !!! end vconf !!!

    !!! allocate !!!

      allocate(dg)
      dg = new_diag(pshape,hole_xi,t,u,jd,lm)
      allocate(wf(4*n_acc,4*n_acc))
      allocate(eg(4*n_acc))
      allocate(kx(mesh))
      allocate(ky(mesh))
      allocate(phi1(mesh,4*n_acc),phi2(mesh,4*n_acc),phi3(mesh,4*n_acc),phi4(mesh,4*n_acc))

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

    !!! end num !!!

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

      !write(cat_num,    *) num+fnum
      !write(cat_z_calc, *) z_calc
      !write(cat_k_range,*) k_range
      !if (k_range==1) then
      !  cat_k_range = "pi"
      !end if
      !filename_fermi_3p4 = "fermi/fermi_3&4_"        // &
      !&                    trim(adjustl(cat_num))    // &
      !&                    "_("                      // &
      !&                    trim(adjustl(cat_z_calc)) // &
      !&                    ")_pi"                    // &
      !&                    trim(adjustl(cat_k_range))// &
      !&                    ".txt"
      !print *, "filename = ", trim(adjustl(filename_fermi_3p4))
      !stop

    !!! end filename !!!

    !!! gnuplot !!!

      !cat_gnuplot = "sp "//""""//trim(adjustl(filename_fermi_3p4))//""""
      !open(20,file="fermi.gpi")
      !  write(20,*) "reset"
      !  write(20,*) "set pm3d map"
      !  write(20,*) "set palette @MATLAB"
      !  write(20,*) "set size ratio 1"
      !  write(20,*) "set yrange [-5:1]"
      !  !write(20,*) "set palette defined (0 ""white"", 1 ""red"")"
      !  write(20,*) trim(adjustl(cat_gnuplot))
      !close(20)

    !!! end gnuplot !!!

    !!! z_calc !!!

      if(z_calc==1) then
        j_rangem = 1
        j_rangep = nx*ny
      else if(z_calc==2) then
        j_rangem = nx*ny+1
        j_rangep = n_acc
      end if

    !!! end z_calc !!!

    !!! mesh !!!

      do p = 1, mesh
        kx(p) = (2*pi*p)/mesh-pi
        if(k_range==0) then
          ky(p) = 0
        else if(k_range==1) then
          ky(p) = (2*pi*p)/mesh-pi
        end if
      end do

    !!! end mesh !!!


  !!!!!! end preferences !!!!!!


  !!!!!! fourier !!!!!!


    !phi1 = 0d0
    !phi2 = 0d0
    phi3 = 0d0
    phi4 = 0d0

    do leig = 2*n_acc+1, 4*n_acc
      do p = 1, mesh
        do j = j_rangem, j_rangep

          i = dg%ets(j)
          rvec = site_to_r(i,pshape)
          jx = rvec(1)
          jy = rvec(2)

          !phi1(p,leig) = phi1(p,leig)+(wf(4*j-3,leig)*exp(ui*kx(p)*jx)*exp(ui*ky(p)*jy))
          !phi2(p,leig) = phi2(p,leig)+(wf(4*j-2,leig)*exp(ui*kx(p)*jx)*exp(ui*ky(p)*jy))
          phi3(p,leig) = phi3(p,leig)+(wf(4*j-1,leig)*exp(ui*kx(p)*jx)*exp(ui*ky(p)*jy))
          phi4(p,leig) = phi4(p,leig)+(wf(4*j  ,leig)*exp(ui*kx(p)*jx)*exp(ui*ky(p)*jy))

        end do
      end do

      !print *, "leig count : ", leig-2*n_acc, "/", 4*n_acc-2*n_acc

    end do

    !phi1(:,:) = (1/(nx*ny)**0.5d0)*phi1(:,:)
    !phi2(:,:) = (1/(nx*ny)**0.5d0)*phi2(:,:)
    phi3(:,:) = (1/(nx*ny)**0.5d0)*phi3(:,:)
    phi4(:,:) = (1/(nx*ny)**0.5d0)*phi4(:,:)

    !phi1(:,:) = phi1(:,:)*conjg(phi1(:,:))
    !phi2(:,:) = phi2(:,:)*conjg(phi2(:,:))
    phi3(:,:) = phi3(:,:)*conjg(phi3(:,:))
    phi4(:,:) = phi4(:,:)*conjg(phi4(:,:))

    print *, "end fourier"


  !!!!!! end fourier !!!!!!


  !!!!!! outputs !!!!!!


    !open(28,file=filename_diag_test6_1)
    !  do leig = 2*n_acc, 2*n_acc+400
    !    do p = mesh/2, mesh
    !      write(28,*) p, eg(leig), abs(phi1(p,leig))
    !    enddo
    !    write(28,*) ""
    !  enddo
    !close(28)

    !open(28,file=filename_diag_test6_2)
    !  do leig = 2*n_acc, 2*n_acc+400
    !    do p = mesh/2, mesh
    !      write(28,*) p, eg(leig), abs(phi2(p,leig))
    !    enddo
    !    write(28,*) ""
    !  enddo
    !close(28)

    !open(28,file=filename_diag_test6_3)
    !  do leig = 2*n_acc+1, 4*n_acc
    !    do p = 1, mesh
    !      write(28,*) kx(p), -eg(leig), abs(phi3(p,leig))
    !    enddo
    !    write(28,*) ""
    !  enddo
    !close(28)

    !open(28,file=filename_diag_test6_4)
    !  do leig = 2*n_acc, 2*n_acc+400
    !    do p = mesh/2, mesh
    !      write(28,*) p, eg(leig), abs(phi4(p,leig))
    !    enddo
    !    write(28,*) ""
    !  enddo
    !close(28)

    open(28,file="fermi.txt")
      do leig = 2*n_acc+1, 4*n_acc
        do p = 1, mesh
          write(28,*) kx(p), -eg(leig), abs(phi3(p,leig))+abs(phi4(p,leig))
        enddo
        write(28,*) ""
      enddo
    close(28)

    print *, "save to fermi.txt"


  !!!!!! end outputs !!!!!!


end program fermiiiiiiiiiiii