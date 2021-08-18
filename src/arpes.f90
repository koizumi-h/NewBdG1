program arpes


  !!!!!! modules !!!!!!


    use math_mod,         only : ui, pi!, my_zheev, fermi, dfermi
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

    real(8), dimension(:,:), allocatable    :: sum_phi1, sum_phi2, sum_phi3, sum_phi4
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
    !real(8), parameter                      :: egm = 0d0, egp = 0.20d0
    real(8)                                 :: mineg, maxeg, egm, egp

    integer, dimension(:,:), allocatable    :: pshape
    integer, dimension(:,:), allocatable    :: hole_xi, hole_chi ! dim is (holes,2)
    integer, parameter                      :: unit_num = 10
    !integer, parameter                      :: z_calc = 2 ! 1 => first layer, 2 => second layer
    !integer, parameter                      :: fchi_select = 2 ! 0 => sc, 1 => fchi, 2 => mwf2
    integer, parameter                      :: symmetry = 0 ! 0 => not symmetry, 1 => 4 fold symmetry
    !integer, parameter                      :: i_repeat = 1000
    integer, parameter                      :: mesh = 100
    integer                                 :: num, snum, fnum
    integer                                 :: nx, ny, nz, n, n_acc, nh!, r(3), time_s, time_e, countpersec
    integer                                 :: vec(3), rvec(3)
    integer                                 :: i, j, jx, jy, px, py
    integer                                 :: j_rangem, j_rangep
    integer                                 :: leig
    integer                                 :: leigm, leigp, leigm2, leigp2
    integer                                 :: fchi_select, z_calc

    character(128)                          :: filename_arpes_3p4
    character(128)                          :: filename_fort
    character(128)                          :: cat_num
    character(128)                          :: cat_egm, cat_egp
    character(128)                          :: cat_z_calc
    character(128)                          :: cat_mineg, cat_maxeg
    character(128)                          :: cat_leigm, cat_leigp
    character(128)                          :: cat_leigm2, cat_leigp2
    character(128)                          :: cat_leig_range
    character(128)                          :: cat_gpi, cat_eps


  !!!!!! end definitions !!!!!!


  !!!!!! preferences !!!!!!


    !!! input !!!

      !write(6,*) "input snum"
      !read (5,*) snum
      !write(6,*) "input fchi_select : 0 => sc, 1 => fchi, 2 => mwf2"
      !read (5,*) fchi_select
      write(6,*) "input fchi_select : 1 => mwf1, 2 => mwf2"
      read (5,*) fchi_select
      write(6,*) "input z_calc : 1 => first layer, 2 => second layer"
      read (5,*) z_calc
      write(6,*) "input egm, egp"
      read (5,*) egm, egp

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
      allocate(phi1(mesh,mesh),phi2(mesh,mesh),phi3(mesh,mesh),phi4(mesh,mesh))
      allocate(sum_phi1(mesh,mesh),sum_phi2(mesh,mesh),sum_phi3(mesh,mesh),sum_phi4(mesh,mesh))

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

      if(fchi_select == 1) then
        read(100) eg
        read(101) wf
      else
        read(110) eg
        read(111) wf
      end if
      !read(100) eg
      !read(101) wf
      !write(cat_num,*) num+fnum
      !filename_fort = 'fort/'//"fort."//trim(adjustl(cat_num))
      !open(num+fnum,file=filename_fort,form='unformatted')
      !  read(num+fnum) eg
      !close(num+fnum)
      !write(cat_num,*) num+fnum+1
      !filename_fort = 'fort/'//"fort."//trim(adjustl(cat_num))
      !open(num+fnum+1,file=filename_fort,form='unformatted')
      !  read(num+fnum+1) wf
      !close(num+fnum+1)
      !do i = 1,4*n_acc
      !  print *, i, eg(i)
      !end do
      !stop

    !!! end eg&wf !!!

    !!! filename !!!

      !write(cat_num,   *) num+fnum
      !!write(cat_leigm, *) leigm
      !write(cat_egm,'(f10.2)') egm
      !!write(cat_leigp, *) leigp
      !write(cat_egp,'(f10.2)') egp
      !write(cat_z_calc,*) z_calc
      !filename_arpes_3p4 = "arpes/arpes_3p4_"       // &
      !&                    trim(adjustl(cat_num))   // &
      !&                    "_("                     // &
      !&                    trim(adjustl(cat_z_calc))// &
      !&                    ")_"                     // &
      !!&                    trim(adjustl(cat_leigm)) // &
      !&                    trim(adjustl(cat_egm))   // &
      !&                    "~"                      // &
      !!&                    trim(adjustl(cat_leigp)) // &
      !&                    trim(adjustl(cat_egp))   // &
      !&                    ".txt"
      !print *, "filename = ", trim(adjustl(filename_arpes_3p4))
      !stop

    !!! end filename !!!

    !!! gnuplot !!!

      !cat_gpi = "sp "//""""//trim(adjustl(filename_arpes_3p4))//""""
      !cat_eps = "eps/"                   // &
      !&         "arpes_3"                // &
      !&         "_"                      // &
      !&         trim(adjustl(cat_num))   // &
      !&         "_"                      // &
      !&         trim(adjustl(cat_egm))   // &
      !&         "~"                      // &
      !&         trim(adjustl(cat_egp))   // &
      !&         "("                      // &
      !&         trim(adjustl(cat_z_calc))// &
      !&         ")"                      // &
      !&         ".eps"
      !cat_eps = "set output "//""""//trim(adjustl(cat_eps))//""""
      !open(20,file="arpes.gpi")
      !  write(20,*) "reset"
      !  write(20,*) "set pm3d map"
      !  write(20,*) "set palette @MATLAB"
      !  write(20,*) "set size ratio 1"
      !  write(20,*) trim(adjustl(cat_gpi))
      !  write(20,*) "reset"
      !  !write(20,*) "set terminal postscript eps color enhanced"
      !  !write(20,*) ""
      !  !write(20,*) trim(adjustl(cat_eps))
      !  !write(20,*) trim(adjustl(cat_gpi))
      !  !write(20,*) "set output"
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

    !!! min_eg_range !!!

      do leig = 1, 4*n_acc
        if(eg(leig)>egm) then
          mineg = eg(leig)
          leigm = leig
          goto 1000
        end if
      end do
      1000 continue
      write(cat_mineg,*) mineg
      write(cat_leigm,*) leigm
      cat_mineg = "mineg = "//trim(adjustl(cat_mineg))
      cat_leigm = "leigm = "//trim(adjustl(cat_leigm))
      print *, trim(adjustl(cat_mineg))
      print *, trim(adjustl(cat_leigm))
      !stop

    !!! end min_eg_range !!!

    !!! max_eg_range !!!

      do leig = 1, 4*n_acc
        if(eg(leig)>egp) then
          maxeg = eg(leig-1)
          leigp = leig-1
          goto 1200
        end if
      end do
      1200 continue
      write(cat_maxeg,*) maxeg
      write(cat_leigp,*) leigp
      cat_maxeg = "maxeg = "//trim(adjustl(cat_maxeg))
      cat_leigp = "leigp = "//trim(adjustl(cat_leigp))
      print *, trim(adjustl(cat_maxeg))
      print *, trim(adjustl(cat_leigp))
      !stop

    !!! end max_eg_range !!!

    !!! error !!!

      if(leigm>leigp) then
        print *, "Error!!!!!!!! There is no eg!!!!!"
        stop
      end if

    !!! end error !!!

    !!! leig_range !!!

      leigm2 = leigm-2*n_acc
      leigp2 = leigp-2*n_acc
      write(cat_leigm2,*) leigm2
      write(cat_leigp2,*) leigp2
      if(leigm2>=0) then
        cat_leigm2 = "+"//trim(adjustl(cat_leigm2))
      end if
      if(leigp2>=0) then
        cat_leigp2 = "+"//trim(adjustl(cat_leigp2))
      end if
      cat_leig_range = "2*n_acc"//trim(adjustl(cat_leigm2))//" ~ "//"2*n_acc"//trim(adjustl(cat_leigp2))
      print *, trim(adjustl(cat_leig_range))
      !stop

    !!! end leig_range !!!

    !!! mesh !!!

      do px = 1, mesh
        kx(px) = (2*pi*px)/mesh-pi
      enddo
      do py = 1, mesh
        ky(py) = (2*pi*py)/mesh-pi
      enddo

    !!! end mesh !!!


  !!!!!! end preferences !!!!!!


  !!!!!! fourier !!!!!!


    !phi1 = 0d0
    !phi2 = 0d0
    phi3 = 0d0
    phi4 = 0d0

    !sum_phi1 = 0d0
    !sum_phi2 = 0d0
    sum_phi3 = 0d0
    sum_phi4 = 0d0

    !do leig = 2*n_acc+leigm2, 2*n_acc+leigp2 ! same as next line
    do leig = leigm, leigp
      do py = 1, mesh
        do px = 1, mesh
          do j = j_rangem, j_rangep

            i = dg%ets(j)
            rvec = site_to_r(i,pshape)
            jx = rvec(1)
            jy = rvec(2)

            !phi1(px,py) = phi1(px,py)+(wf(4*j-3,leig)*exp(ui*kx(px)*jx)*exp(ui*ky(py)*jy))
            !phi2(px,py) = phi2(px,py)+(wf(4*j-2,leig)*exp(ui*kx(px)*jx)*exp(ui*ky(py)*jy))
            phi3(px,py) = phi3(px,py)+(wf(4*j-1,leig)*exp(ui*kx(px)*jx)*exp(ui*ky(py)*jy))
            phi4(px,py) = phi4(px,py)+(wf(4*j  ,leig)*exp(ui*kx(px)*jx)*exp(ui*ky(py)*jy))

          end do
        end do
      end do

      !phi1(:,:) = (1/(nx*ny)**0.5d0)*phi1(:,:)
      !phi2(:,:) = (1/(nx*ny)**0.5d0)*phi2(:,:)
      phi3(:,:) = (1/(nx*ny)**0.5d0)*phi3(:,:)
      phi4(:,:) = (1/(nx*ny)**0.5d0)*phi4(:,:)

      !phi1(:,:) = phi3(:,:)*conjg(phi1(:,:))
      !phi2(:,:) = phi3(:,:)*conjg(phi2(:,:))
      phi3(:,:) = phi3(:,:)*conjg(phi3(:,:))
      phi4(:,:) = phi3(:,:)*conjg(phi4(:,:))

      !sum_phi1(:,:) = sum_phi3(:,:)+phi1(:,:)
      !sum_phi2(:,:) = sum_phi3(:,:)+phi2(:,:)
      sum_phi3(:,:) = sum_phi3(:,:)+phi3(:,:)
      sum_phi4(:,:) = sum_phi3(:,:)+phi4(:,:)

      print *, "leig_count : ", leig-2*n_acc-leigm2+1, "/", abs(leigp2-leigm2)+1

    end do

    !print *, "end fourier"


  !!!!!! end fourier !!!!!!


  !!!!!! outputs !!!!!!


    !open(26,file="diag_sb2_5_1.txt")
    !  do py = , mesh
    !    do px = 1, mesh
    !      write(26,*) kx(px), ky(py), sum_phi1(px,py)
    !    enddo
    !    write(26,*) ""
    !  enddo
    !close(26)

    !open(26,file="diag_sb2_5_2.txt")
    !  do py = , mesh
    !    do px = 1, mesh
    !      write(26,*) kx(px), ky(py), sum_phi2(px,py)
    !    enddo
    !    write(26,*) ""
    !  enddo
    !close(26)

    !open(26,file="diag_sb2_5_3.txt")
    !  do py = 1, mesh
    !    do px = 1, mesh
    !      write(26,*) kx(px), ky(py), sum_phi3(px,py)
    !    enddo
    !    write(26,*) ""
    !  enddo
    !close(26)

    !open(26,file="diag_sb2_5_4.txt")
    !  do py = 1, mesh
    !    do px = 1, mesh
    !      write(26,*) kx(px), ky(py), sum_phi4(px,py)
    !    enddo
    !    write(26,*) ""
    !  enddo
    !close(26)

    open(26,file="arpes.txt")
      do py = 1, mesh
        do px = 1, mesh
          if    (symmetry==0) then
            write(26,*) kx(px), ky(py), sum_phi3(px,py)+sum_phi4(px,py)
          elseif(symmetry==1) then
            write(26,*) kx(px), ky(py), sum_phi3(px,py)               &
            &                          +sum_phi3(mesh-px+1,py)        &
            &                          +sum_phi3(px,mesh-py+1)        &
            &                          +sum_phi3(mesh-px+1,mesh-py+1) &
            &                          +sum_phi3(py,px)               &
            &                          +sum_phi3(mesh-py+1,px)        &
            &                          +sum_phi3(py,mesh-px+1)        &
            &                          +sum_phi3(mesh-py+1,mesh-px+1)
          end if
        enddo
        write(26,*) ""
      enddo
    close(26)

    print *, "save to ""arpes.txt"""


  !!!!!! end outputs !!!!!!


end program arpes