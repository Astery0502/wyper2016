!> Magnetic bipolar field
module mod_usr
  use mod_mhd
  implicit none
  integer, parameter :: jmax=5000
  double precision, allocatable :: pbc(:),rbc(:)
  ! 1D solar atmosphere table for pressure, density, and height
  double precision :: pa(jmax),ra(jmax),ya(jmax)
  double precision :: usr_grav,SRadius,rhob,Tiso,dr,gzone,bQ0,heatunit

  double precision :: dh, Bh, dv, Bv, xv, dp ! fan-spine parameters
  double precision :: Bl, Br, cv, ch, htra
  double precision :: trelax, tstop
  double precision :: rm, zfac

contains

  !==============================================================================
  ! Purpose: to include global parameters, set user methods, set coordinate 
  !          system and activate physics module.
  !==============================================================================
  subroutine usr_init()

    unit_length        = 1.d9 ! cm
    unit_temperature   = 1.d6 ! K
    unit_numberdensity = 1.d9 ! cm^-3

    usr_set_parameters  => initglobaldata_usr
    usr_init_one_grid   => initonegrid_usr
    usr_special_bc      => specialbound_usr
    ! usr_source          => specialsource
    usr_gravity         => gravity
    usr_refine_grid     => special_refine_grid
    usr_init_vector_potential=>initvecpot_usr
    usr_aux_output      => specialvar_output
    usr_set_B0          => specialset_B0
    usr_set_J0          => specialset_J0

    ! usr_process_grid    => my_process_grid 
    ! usr_write_analysis  => my_analysis

    call params_read(par_files)
    call set_coordinate_system("Cartesian_2D")
    call mhd_activate()

  end subroutine usr_init

  !> Read parameters from a file
  subroutine params_read(files)
    use mod_global_parameters, only: unitpar
    character(len=*), intent(in) :: files(:)
    integer                      :: n
 
    namelist /my_list/ dh, Bh, dv, Bv, xv, dp, &
            trelax, tstop
 
    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, my_list, end=111)
111    close(unitpar)
    end do
 
  end subroutine params_read

  !==============================================================================
  ! Purpose: to initialize user public parameters and reset global parameters.
  !          Input data are also read here.
  !==============================================================================
  subroutine initglobaldata_usr()

    !unit_density       = 1.4d0*mass_H*unit_numberdensity               ! 2.341668000000000E-015 g*cm^-3
    !unit_pressure      = 2.3d0*unit_numberdensity*k_B*unit_temperature ! 0.317538000000000 erg*cm^-3
    !unit_magneticfield = dsqrt(miu0*unit_pressure)                     ! 1.99757357615242 Gauss
    !unit_velocity      = unit_magneticfield/dsqrt(miu0*unit_density)   ! 1.16448846777562E007 cm/s = 116.45 km/s
    !unit_time          = unit_length/unit_velocity                     ! 85.8746159942810 s 
    heatunit=unit_pressure/unit_time !< 3.697693390805347E-003 erg*cm^-3/s
    if(.not.mhd_energy) then
      ! bottom density
      rhob=2.d0
      ! isothermal uniform temperature
      Tiso= mhd_adiab
    end if

    bQ0=1.d-4/unit_pressure*unit_time          ! 3.697693390805347E-003 erg*cm^-3/s

    gzone=0.2d0+xprobmin2
    ! cell size in 1D solar atmosphere table
    dr=(2.d0*gzone+xprobmax2-xprobmin2)/dble(jmax)
    usr_grav=-2.74d4*unit_length/unit_velocity**2  ! solar gravity
    SRadius=6.955d10/unit_length                   ! Solar radius

    dh = dh + dp
    dv = dv + dp
    ! for fan-spine field use in wyper2016a_field
    cv = (Bv*(dv-dp)**3.d0)*half
    ch = (Bh*(dh-dp)**3.d0)*half

    if(mhd_energy) call inithdstatic

  end subroutine initglobaldata_usr

  !> initialize solar atmosphere table in a vertical line through the global domain
  subroutine inithdstatic
    use mod_global_parameters

    double precision :: Ta(jmax),gg(jmax)
    double precision :: rpho,Tpho,Ttop,wtra,ftra,Ttr,Fc,k_para
    double precision :: res,pb,rhob,invT
    integer :: j,na,ibc,btlevel

    rpho=1.151d15/unit_numberdensity ! number density at the bottom relaxla
    Tpho=8.d3/unit_temperature ! temperature of chromosphere
    Ttop=1.5d6/unit_temperature ! estimated temperature in the top
    htra=2.d8/unit_length ! height of initial transition region
    wtra=2.d7/unit_length ! width of initial transition region 
    Ttr=1.6d5/unit_temperature ! lowest temperature of upper profile
    Fc=2.d5/unit_pressure/unit_velocity  ! constant thermal conduction flux
    ftra=wtra*atanh(2.d0*(Ttr-Tpho)/(Ttop-Tpho)-1.d0)
    ! Spitzer thermal conductivity with cgs units
    k_para=8.d-7*unit_temperature**3.5d0/unit_length/unit_density/unit_velocity**3 

    !! set T distribution with height
    do j=1,jmax
       ya(j)=(dble(j)-0.5d0)*dr+xprobmin2-gzone
       if(ya(j)>htra) then
         Ta(j)=(3.5d0*Fc/k_para*(ya(j)-htra)+Ttr**3.5d0)**(2.d0/7.d0)
       else
         Ta(j)=Tpho+0.5d0*(Ttop-Tpho)*(tanh((ya(j)-htra+ftra)/wtra)+1.d0)
       endif
       gg(j)=usr_grav*(SRadius/(SRadius+ya(j)))**2
    enddo
    !! solution of hydrostatic equation 
    ra(1)=rpho
    pa(1)=rpho*Tpho

    do j=2,jmax
       pa(j)=(pa(j-1)+dr*(gg(j)+gg(j-1))*ra(j-1)/4.d0)/(one-dr*(gg(j)+gg(j-1))/&
              Ta(j)/4.d0)
       ra(j)=pa(j)/Ta(j)
    end do

    !! initialized rho and p in the fixed bottom boundary
    na=floor(gzone/dr+0.5d0)
    res=gzone-(dble(na)-0.5d0)*dr
    rhob=ra(na)+res/dr*(ra(na+1)-ra(na))
    pb=pa(na)+res/dr*(pa(na+1)-pa(na))
    allocate(rbc(nghostcells))
    allocate(pbc(nghostcells))
    btlevel=refine_max_level
    do ibc=nghostcells,1,-1
      na=floor((gzone-dx(ndim,btlevel)*(dble(nghostcells-ibc+1)-0.5d0))/dr+0.5d0)
      res=gzone-dx(ndim,btlevel)*(dble(nghostcells-ibc+1)-0.5d0)-(dble(na)-0.5d0)*dr
      rbc(ibc)=ra(na)+res/dr*(ra(na+1)-ra(na))
      pbc(ibc)=pa(na)+res/dr*(pa(na+1)-pa(na))
    end do
    
    if (mype==0) then
     print*,'minra',minval(ra)
     print*,'maxTa',Ta(jmax)
     print*,'rhob',rhob
     print*,'pb',pb
    endif
  end subroutine inithdstatic

  subroutine initonegrid_usr(ixI^L,ixO^L,w,x)

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    double precision :: A(ixI^S,1:ndim)
    double precision :: Bfr(ixI^S,1:ndir)
    double precision :: res
    integer :: ix^D,na
    logical, save :: first=.true.

    if(first)then
      if(mype==0) then
      write(*,*)'Fan-spine field from Wyper2016a'
      write(*,*)'B0field', B0field
      write(*,*)'nghostcells', nghostcells
    endif
      first=.false.
    endif

    if(B0field) then
      w(ixO^S,mag(:))=zero
    else
      call wyper2016a_field(ixI^L,ixO^L,x,Bfr)
      w(ixO^S,mag(:))=Bfr(ixO^S,:)
    end if

    if(stagger_grid) then
      call b_from_vector_potential(block%ixGs^L,ixI^L,ixO^L,block%ws,x)
      call mhd_face_to_center(ixO^L,block)
    end if

    w(ixO^S,mom(:))=0.d0

    if(mhd_energy) then
      {do ix^DB=ixOmin^DB,ixOmax^DB\}
         na=floor((x(ix^D,ndim)-xprobmin2+gzone)/dr+0.5d0)
         res=x(ix^D,ndim)-xprobmin2+gzone-(dble(na)-0.5d0)*dr
         w(ix^D,rho_)=ra(na)+(one-cos(dpi*res/dr))/two*(ra(na+1)-ra(na))
         w(ix^D,p_)  =pa(na)+(one-cos(dpi*res/dr))/two*(pa(na+1)-pa(na))
      {end do\}
    else if(mhd_adiab/=0) then
      ! isothermal
      w(ixO^S,rho_)=rhob*dexp(usr_grav*SRadius**2/Tiso*&
                    (1.d0/SRadius-1.d0/(x(ixO^S,ndim)+SRadius)))
    else
      ! zero beta
      w(ixO^S,rho_)=sum(w(ixO^S,mag(:))**2,dim=ndim+1)
    end if

    call mhd_to_conserved(ixI^L,ixO^L,w,x)

  end subroutine initonegrid_usr

  subroutine initvecpot_usr(ixI^L, ixC^L, xC, A, idir)
    ! initialize the vectorpotential on the edges
    ! used by b_from_vectorpotential()
    use mod_global_parameters
    integer, intent(in)                :: ixI^L, ixC^L,idir
    double precision, intent(in)       :: xC(ixI^S,1:ndim)
    double precision, intent(out)      :: A(ixI^S)

    ! vector potential
    double precision :: Avec1(ixI^S,1:ndim)

    ! call bipolar_field(ixI^L,ixC^L,xC,Avec1)
    if (B0field) then
      Avec1(ixC^S,:)=0.d0
    else
      call wyper2016a_field(ixI^L,ixC^L,xC,Avec1)
    end if

    if(idir==2) then 
      A(ixC^S)=Avec1(ixC^S,2)
    else
      A(ixC^S)=Avec1(ixC^S,1)
    end if

  end subroutine initvecpot_usr

  subroutine wyper2016a_field(ixI^L,ixO^L,x,Bvec)

    ! from wyper to here: x --> z; y --> x; z --> y
    ! from here to wyper: x --> y; y --> z; z --> x

    integer, intent(in)  :: ixI^L,ixO^L
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision, intent(out) :: Bvec(ixI^S,1:ndir)

    double precision :: xh, zh, zv

    xh = 0.d0
    zh = -1 * dh

    ! xv already defined in amrvac.par
    zv = -1 * dv

    Bvec(ixO^S,1) = -3.d0*ch*((x(ixO^S,2)-zh)**2) / ((x(ixO^S,1)-xh)**2+(x(ixO^S,2)-zh)**2)**(5./2.) + &
                    3.d0*cv*(x(ixO^S,1)-xv)*(x(ixO^S,2)-zv) / ((x(ixO^S,1)-xv)**2+(x(ixO^S,2)-zv)**2)**(5./2.) + &
                    2*ch / ((x(ixO^S,1)-xh)**2+(x(ixO^S,2)-zh)**2)**(3./2.)

    Bvec(ixO^S,2) = 3.d0*ch*(x(ixO^S,1)-xh)*(x(ixO^S,2)-zh) / ((x(ixO^S,1)-xh)**2+(x(ixO^S,2)-zh)**2)**(5./2.) + & 
                    (-3.d0)*cv*((x(ixO^S,1)-xv)**2) / ((x(ixO^S,1)-xv)**2+(x(ixO^S,2)-zv)**2)**(5./2.) + & 
                    2.d0*cv/((x(ixO^S,1)-xv)**2+(x(ixO^S,2)-zv)**2)**(3./2.)

  end subroutine wyper2016a_field

  subroutine specialbound_usr(qt,ixI^L,ixO^L,iB,w,x)
        ! special boundary types, user defined
    integer, intent(in) :: ixO^L, iB, ixI^L
    double precision, intent(in) :: qt, x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    double precision :: pth(ixI^S),tmp(ixI^S),ggrid(ixI^S),invT(ixI^S)
    double precision :: Q(ixI^S),Qp(ixI^S)
    integer :: ix^D,ixOs^L,ixC^L,hxC^L,jxO^L,idir

    select case(iB)
    case(3)
      if(iprob<3) then
        !! fixed zero velocity
        do idir=1,ndir
          w(ixO^S,mom(idir))=-w(ixOmin1:ixOmax1,ixOmax2+nghostcells:ixOmax2+1:-1,mom(idir))&
                     /w(ixOmin1:ixOmax1,ixOmax2+nghostcells:ixOmax2+1:-1,rho_)
        end do
      end if
      !! fixed b1 b2 b3
      if(iprob==0 .or. B0field) then
        w(ixO^S,mag(:))=0.d0
      else if(stagger_grid) then
        do idir=1,nws
          if(idir==2) cycle
          ixOsmax^D=ixOmax^D;
          ixOsmin^D=ixOmin^D-kr(^D,idir);
          ! 2nd order one-sided zero-gradient extrapolation
          !do ix2=ixOsmax2,ixOsmin2,-1
          !   block%ws(ix2^%2ixOs^S,idir)=1.d0/3.d0*&
          !         (-block%ws(ix2+2^%2ixOs^S,idir)&
          !     +4.d0*block%ws(ix2+1^%2ixOs^S,idir))
          !end do
          ! 4th order one-sided equal-gradient extrapolation
          do ix2=ixOsmax2,ixOsmin2,-1
            block%ws(ix2^%2ixOs^S,idir)= &
              0.12d0*block%ws(ix2+5^%2ixOs^S,idir) &
             -0.76d0*block%ws(ix2+4^%2ixOs^S,idir) &
             +2.08d0*block%ws(ix2+3^%2ixOs^S,idir) &
             -3.36d0*block%ws(ix2+2^%2ixOs^S,idir) &
             +2.92d0*block%ws(ix2+1^%2ixOs^S,idir)
          end do
        end do
        ixOs^L=ixO^L-kr(2,^D);
        jxO^L=ixO^L+nghostcells*kr(2,^D);
        block%ws(ixOs^S,2)=zero
        do ix2=ixOsmax2,ixOsmin2,-1
          call get_divb(w,ixI^L,ixO^L,Qp)
          block%ws(ix2^%2ixOs^S,2)=Qp(ix2+1^%2ixO^S)*block%dvolume(ix2+1^%2ixO^S)&
            /block%surfaceC(ix2^%2ixOs^S,2)
        end do
        call mhd_face_to_center(ixO^L,block)
      endif
      !! fixed gravity stratification of density and pressure pre-determined in initial condition
      do ix2=ixOmin2,ixOmax2
        w(ixOmin1:ixOmax1,ix2,rho_)=rbc(ix2)
        w(ixOmin1:ixOmax1,ix2,p_)=pbc(ix2)
      enddo
      if(mhd_hyperbolic_thermal_conduction) then
        do ix2=ixOmin2,ixOmax2
          do ix1=ixOmin1,ixOmax1
            w(ix1,ix2,q_)=w(ix1,ixOmax2+1,q_)
          end do
        end do
      end if
      call mhd_to_conserved(ixI^L,ixO^L,w,x)
    case(4)
      ixOs^L=ixO^L;
      ixOsmin2=ixOmin2-1;ixOsmax2=ixOmin2-1;
      call mhd_get_pthermal(w,x,ixI^L,ixOs^L,pth)
      ixOsmin2=ixOmin2-1;ixOsmax2=ixOmax2;
      call getggrav(ggrid,ixI^L,ixOs^L,x)
      !> fill pth, rho ghost layers according to gravity stratification
      invT(ixOmin2-1^%2ixO^S)=w(ixOmin2-1^%2ixO^S,rho_)/pth(ixOmin2-1^%2ixO^S)
      tmp=0.d0
      do ix2=ixOmin2,ixOmax2
        tmp(ixOmin2-1^%2ixO^S)=tmp(ixOmin2-1^%2ixO^S)+0.5d0*&
            (ggrid(ix2^%2ixO^S)+ggrid(ix2-1^%2ixO^S))*invT(ixOmin2-1^%2ixO^S)
        w(ix2^%2ixO^S,p_)=pth(ixOmin2-1^%2ixO^S)*dexp(tmp(ixOmin2-1^%2ixO^S)*dxlevel(2))
        w(ix2^%2ixO^S,rho_)=w(ix2^%2ixO^S,p_)*invT(ixOmin2-1^%2ixO^S)
      enddo
      !> fixed zero velocity
      do idir=1,ndir
        w(ixO^S,mom(idir)) =-w(ixOmin1:ixOmax1,ixOmin2-1:ixOmin2-nghostcells:-1,mom(idir))&
                     /w(ixOmin1:ixOmax1,ixOmin2-1:ixOmin2-nghostcells:-1,rho_)
      end do
      !> zero normal gradient extrapolation
      if(stagger_grid) then
        do idir=1,nws
          if(idir==2) cycle
          ixOsmax^D=ixOmax^D;
          ixOsmin^D=ixOmin^D-kr(^D,idir);
          do ix2=ixOsmin2,ixOsmax2
             block%ws(ix2^%2ixOs^S,idir)=1.d0/3.d0*&
                   (-block%ws(ix2-2^%2ixOs^S,idir)&
               +4.d0*block%ws(ix2-1^%2ixOs^S,idir))
          end do
        end do
        ixOs^L=ixO^L;
        jxO^L=ixO^L-nghostcells*kr(2,^D);
        block%ws(ixOs^S,2)=zero
        call get_divb(w,ixI^L,jxO^L,Q)
        do ix2=ixOsmin2,ixOsmax2
          call get_divb(w,ixI^L,ixO^L,Qp)
          block%ws(ix2^%2ixOs^S,2)=&
            (Q(jxOmax2^%2jxO^S)*block%dvolume(jxOmax2^%2jxO^S)&
           -Qp(ix2^%2ixO^S)*block%dvolume(ix2^%2ixO^S))&
            /block%surfaceC(ix2^%2ixOs^S,2)
        end do
        call mhd_face_to_center(ixO^L,block)
        do ix2=ixOmin2,ixOmax2
          w(ixOmin1:ixOmax1,ix2,mag(3))=(1.0d0/3.0d0)* &
                      (-w(ixOmin1:ixOmax1,ix2-2,mag(3))&
                 +4.0d0*w(ixOmin1:ixOmax1,ix2-1,mag(3)))
        enddo
      else
        do ix2=ixOmin2,ixOmax2
          w(ixOmin1:ixOmax1,ix2,mag(:))=(1.0d0/3.0d0)* &
                      (-w(ixOmin1:ixOmax1,ix2-2,mag(:))&
                 +4.0d0*w(ixOmin1:ixOmax1,ix2-1,mag(:)))
        enddo
      end if
      if(mhd_hyperbolic_thermal_conduction) then
        do ix2=ixOmin2,ixOmax2
          do ix1=ixOmin1,ixOmax1
            w(ix1,ix2,q_)=w(ix1,ixOmin2-1,q_)
          end do
        end do
      end if
      call mhd_to_conserved(ixI^L,ixO^L,w,x)
    case default
       call mpistop("Special boundary is not defined for this region")
    end select

  end subroutine specialbound_usr

  subroutine getggrav(ggrid,ixI^L,ixO^L,x)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(out)   :: ggrid(ixI^S)

    ggrid(ixO^S)=usr_grav*(SRadius/(SRadius+x(ixO^S,ndim)))**2
  end subroutine getggrav

  !==============================================================================
  ! Purpose: get gravity field
  !==============================================================================
  subroutine gravity(ixI^L,ixO^L,wCT,x,gravity_field)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(in)    :: wCT(ixI^S,1:nw)
    double precision, intent(out)   :: gravity_field(ixI^S,ndim)
    double precision                :: ggrid(ixI^S)

    gravity_field=0.d0
    call getggrav(ggrid,ixI^L,ixO^L,x)
    gravity_field(ixO^S,ndim)=ggrid(ixO^S)
  end subroutine gravity

  subroutine specialsource(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x)
    use mod_global_parameters

    integer, intent(in) :: ixI^L, ixO^L, iw^LIM
    double precision, intent(in) :: qdt, qtC, qt, x(ixI^S,1:ndim), wCT(ixI^S,1:nw)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    double precision :: lQgrid(ixI^S),bQgrid(ixI^S)

    ! add global background heating bQ
    call getbQ(bQgrid,ixI^L,ixO^L,qtC,wCT,x)
    w(ixO^S,e_)=w(ixO^S,e_)+qdt*bQgrid(ixO^S)
    !! add localized heating lQ
    !if(iprob==2) then
    !  call getlQ(lQgrid,ixI^L,ixO^L,qtC,wCT,x)
    !  w(ixO^S,e_)=w(ixO^S,e_)+qdt*lQgrid(ixO^S)
    !endif

  end subroutine specialsource

  subroutine getbQ(bQgrid,ixI^L,ixO^L,qt,w,x)
    !!calculate background heating bQ, Mok 2016 ApJ
    use mod_global_parameters

    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: qt, x(ixI^S,1:ndim), w(ixI^S,1:nw)
    double precision, intent(out) :: bQgrid(ixI^S)

    double precision :: Bmag(ixI^S),bvec(ixI^S,1:ndir),cvec(ixI^S,1:ndir),tmp(ixI^S)
    double precision :: curlo,curhi
    integer :: idims,idir

    bQgrid(ixO^S)=bQ0*dexp(-x(ixO^S,ndim)/6.d0)
    !if(B0field) then
    !  Bmag(ixI^S)=dsqrt(sum((w(ixI^S,mag(:))+block%B0(ixI^S,:,0))**2,dim=ndim+1))
    !  do idir=1,ndir
    !    bvec(ixI^S,idir)=(w(ixI^S,mag(idir))+block%B0(ixI^S,idir,0))/Bmag(ixI^S)
    !  end do
    !else
    !  Bmag(ixI^S)=dsqrt(sum(w(ixI^S,mag(:))**2,dim=ndim+1))
    !  do idir=1,ndir
    !    bvec(ixI^S,idir)=w(ixI^S,mag(idir))/Bmag(ixI^S)
    !  end do
    !endif
    !cvec=0.d0
    !! calculate local curvature of magnetic field
    !do idims=1,ndim
    !  call gradient(bvec(ixI^S,1),ixI^L,ixO^L,idims,tmp) 
    !  cvec(ixO^S,1)=cvec(ixO^S,1)+bvec(ixO^S,idims)*tmp(ixO^S)
    !  call gradient(bvec(ixI^S,2),ixI^L,ixO^L,idims,tmp) 
    !  cvec(ixO^S,2)=cvec(ixO^S,2)+bvec(ixO^S,idims)*tmp(ixO^S)
    !  call gradient(bvec(ixI^S,3),ixI^L,ixO^L,idims,tmp) 
    !  cvec(ixO^S,3)=cvec(ixO^S,3)+bvec(ixO^S,idims)*tmp(ixO^S)
    !end do 
    !tmp(ixO^S)=dsqrt(sum(cvec(ixO^S,:)**2,dim=ndim+1))
    !! set lower and upper limit for curvature
    !curlo=0.1d0 ! 1/10 
    !curhi=2.d0 ! 1/0.5
    !where(tmp(ixO^S)<curlo)
    !  tmp(ixO^S)=curlo
    !elsewhere(tmp(ixO^S)>curhi)
    !  tmp(ixO^S)=curhi
    !end where
    !
    !bQgrid(ixO^S)=bQ0*Bmag(ixO^S)**1.75d0*w(ixO^S,rho_)**0.125d0*tmp(ixO^S)**0.75d0
  
  end subroutine getbQ

  !==============================================================================
  ! Purpose: Enforce additional refinement or coarsening. One can use the
  !          coordinate info in x and/or time qt=t_n and w(t_n) values w.
  !==============================================================================
  subroutine special_refine_grid(igrid,level,ixI^L,ixO^L,qt,w,x,refine,coarsen)
    use mod_global_parameters

    integer, intent(in) :: igrid, level, ixI^L, ixO^L
    double precision, intent(in) :: qt, w(ixI^S,1:nw), x(ixI^S,1:ndim)
    integer, intent(inout) :: refine, coarsen

    refine=-1
    coarsen=-1

  end subroutine special_refine_grid

  !==============================================================================
  ! Purpose: 
  !   this subroutine can be used in convert, to add auxiliary variables to the
  !   converted output file, for further analysis using tecplot, paraview, ....
  !   these auxiliary values need to be stored in the nw+1:nw+nwauxio slots
  !
  !   the array normconv can be filled in the (nw+1:nw+nwauxio) range with
  !   corresponding normalization values (default value 1)
  !==============================================================================
  subroutine specialvar_output(ixI^L,ixO^L,w,x,normconv)
    use mod_global_parameters
    integer, intent(in)                :: ixI^L,ixO^L
    double precision, intent(in)       :: x(ixI^S,1:ndim)
    double precision                   :: w(ixI^S,nw+nwauxio)
    double precision                   :: normconv(0:nw+nwauxio)
    double precision                   :: tmp(ixI^S),dip(ixI^S),divb(ixI^S),B2(ixI^S)
    double precision, dimension(ixI^S,1:ndir) :: Btotal,qvec,curlvec
    integer                            :: ix^D,idirmin,idims,idir,jdir,kdir, nwspecial

    nwspecial = 1

    ! Btotal & B^2
    if(B0field) then
      Btotal(ixI^S,1:ndir)=w(ixI^S,mag(1:ndir))+block%B0(ixI^S,1:ndir,b0i)
    else
      Btotal(ixI^S,1:ndir)=w(ixI^S,mag(1:ndir))
    end if
    B2(ixO^S)=sum((Btotal(ixO^S,:))**2,dim=ndim+1)
    ! output Alfven wave speed B/sqrt(rho)
    ! w(ixO^S,nw+nwspecial)=dsqrt(B2(ixO^S)/w(ixO^S,rho_))
    ! nwspecial = nwspecial + 1
    ! output divB1
    ! call get_divb(w,ixI^L,ixO^L,divb)
    ! w(ixO^S,nw+nwspecial)=divb(ixO^S)
    ! nwspecial = nwspecial + 1
    ! output the plasma beta p*2/B**2
    ! call mhd_get_pthermal(w,x,ixI^L,ixO^L,tmp)
    ! where(B2(ixO^S)/=0.d0)
    !   w(ixO^S,nw+3)=2.d0*tmp(ixO^S)/B2(ixO^S)
    ! else where
    !   w(ixO^S,nw+3)=0.d0
    ! end where
    ! store current
    call curlvector(Btotal,ixI^L,ixO^L,curlvec,idirmin,1,ndir)
    do idir=1,ndir
      w(ixO^S,nw+nwspecial)=curlvec(ixO^S,idir)
      nwspecial = nwspecial + 1
    end do
    ! calculate Lorentz force
    ! qvec(ixO^S,1:ndir)=zero
    ! do idir=1,ndir; do jdir=1,ndir; do kdir=idirmin,3
    !   if(lvc(idir,jdir,kdir)/=0)then
    !     tmp(ixO^S)=curlvec(ixO^S,jdir)*w(ixO^S,mag(kdir))
    !     if(lvc(idir,jdir,kdir)==1)then
    !       qvec(ixO^S,idir)=qvec(ixO^S,idir)+tmp(ixO^S)
    !     else
    !       qvec(ixO^S,idir)=qvec(ixO^S,idir)-tmp(ixO^S)
    !     endif
    !   endif
    ! enddo; enddo; enddo
    ! do idir=1,ndir
    !   w(ixO^S,nw+nwspecial)=qvec(ixO^S,idir)
    !   nwspecial = nwspecial + 1
    ! end do
    ! find magnetic dips
    !dip=0.d0
    !do idir=1,ndir
    !  call gradient(w(ixI^S,mag(3)),ixI^L,ixO^L,idir,tmp)
    !  dip(ixO^S)=dip(ixO^S)+w(ixO^S,b0_+idir)*tmp(ixO^S)
    !end do
    !where(dabs(w(ixO^S,mag(3)))<0.08d0 .and. dip(ixO^S)>=0.d0)
    !  w(ixO^S,nw+8)=1.d0
    !elsewhere
    !  w(ixO^S,nw+8)=0.d0
    !end where
  end subroutine specialvar_output

  subroutine specialset_B0(ixI^L,ixO^L,x,wB0)
  ! Here add a steady (time-independent) potential or 
  ! linear force-free background field
    integer, intent(in)           :: ixI^L,ixO^L
    double precision, intent(in)  :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: wB0(ixI^S,1:ndir)

    double precision :: A(ixI^S,1:ndim)

    call wyper2016a_field(ixI^L,ixO^L,x,wB0)

  end subroutine specialset_B0

  subroutine specialset_J0(ixI^L,ixO^L,x,wJ0)
  ! Here add a time-independent background current density 
    integer, intent(in)           :: ixI^L,ixO^L
    double precision, intent(in)  :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: wJ0(ixI^S,7-2*ndir:ndir)

    wJ0(ixO^S,:)=0.d0

  end subroutine specialset_J0

end module mod_usr
