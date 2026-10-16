program glowbasic_Heion

! This software is part of the GLOW model.  Use is governed by the Open Source
! Academic Research License Agreement contained in the file glowlicense.txt.
! For more information see the file glow.txt.

! Version 0.981, 6/2017

! Adapted from glowdriver by Stan Solomon, 2/2016

! Basic single-processor driver for the GLOW model.
! Uses MSIS/IRI for input.
! Runs GLOW for designated inputs once, or multiple times.
! MPI and netCDF libraries not required.

! For definitions of use-associated variables, see subroutine GLOW and module CGLOW.
! Consider recombination of He+ and e-

! Other definitions:
! f107p   Solar 10.7 cm flux for previous day
! ap      Ap index of geomagnetic activity
! z       altitude array, km

! Array dimensions:
! jmax    number of altitude levels
! nbins   number of energetic electron energy bins
! lmax    number of wavelength intervals for solar flux
! nmaj    number of major species
! nst     number of states produced by photoionization/dissociation
! nei     number of states produced by electron impact
! nex     number of ionized/excited species
! nw      number of airglow emission wavelengths
! nc      number of component production terms for each emission

  use cglow,only: jmax,nbins,lmax,nmaj,nei,nex,nw,nc,nst
  use cglow,only: idate,ut,glat,glong,f107a,f107,f107p,ap,ef,ec
  use cglow,only: iscale,jlocal,kchem,xuvfac
  use cglow,only: sza,dip,efrac,ierr
  use cglow,only: zz,zo,zn2,zo2,zns,znd,zno,ztn,ze,zti,zte
  use cglow,only: ener,del,phitop,wave1,wave2,sflux,pespec,sespec,uflx,dflx,sion
  use cglow,only: photoi,photod,phono,aglw,tei,tpi,tir,ecalc,zxden,zeta,zceta,zlbh
  use cglow,only: cglow_init
  use cglow,only: data_dir
  use write_path

  implicit none

  character(len=1024) :: iri90_dir

  real,allocatable :: z(:)                    ! glow height coordinate in km (jmax)
  real,allocatable :: zun(:), zvn(:)          ! neutral wind components (not in use)
  real,allocatable :: pedcond(:), hallcond(:) ! Pederson and Hall conductivities in S/m (mho)
  real,allocatable :: outf(:,:)               ! iri output (11,jmax)
  real :: rz12,stl,fmono,emono, sza0
  real :: d(8), t(2), sw(25), oarr(30)
  integer :: l,j,jj,ijf,jmag,iday,mmdd,i,ii,n,k,ix,itail
  integer :: instance,iostatus
  !logical :: jf(12)
  data sw/25*1./

      real, allocatable :: z_file(:), uflx_file(:,:), dflx_file(:,:)  ! 用于文件中的临时数据
      integer :: file_j, prog_j, jmax_file, start_index

!  additional varity
  character(len=5) :: inst, doyc
  character(len=4) :: DawnDusk
  real :: cglat(2), cglong(2), ltstart, ltend, ppphitop
  integer :: flag, iii, doyi, year, doy, api, insti
  real,allocatable :: pphitop(:)
  character(len=256) :: filepath
!
! Initialize standard switches:
!
  iscale=1
  xuvfac=3.
  kchem=4
  jlocal=0
  itail=0
  fmono=0.
  emono=0.
!
! Set data directories:
!
  data_dir    = 'data/'
  iri90_dir   = 'data/iri90/'
!
! Set number of altitude levels:
!
  jmax = 138
    year=2024
!
! Allocate local arrays:
!
  allocate(z(jmax))
  allocate(zun(jmax))
  allocate(zvn(jmax))
  allocate(pedcond(jmax))
  allocate(hallcond(jmax))
  allocate(outf(11,jmax))

  allocate(pphitop(jmax))
!
! Call CGLOW_INIT (module CGLOW) to set array dimensions and allocate use-associated variables:
! (This was formerly done using common blocks, including common block /cglow/.)
!
  call cglow_init
!
! Call EGRID to set up electron energy grid:
!
  call egrid (ener, del, nbins)
!
! Loop to call GLOW for designated inputs until end-of-file or any character on standard input:
!

  
  cglat(1)=19.5
  cglong(1)=109.
  !cglat(2)=-5.66
  !cglong(2)=109.24
  !open(10,file = 'C://Users//31217//Desktop//氦激光雷达//IGRF//conjN.txt')
  !  read(10, *) cglat(1)
  !  read(10, *) cglong(1)
  !close(10)
  !open(10,file = 'C://Users//31217//Desktop//氦激光雷达//IGRF//conjS.txt')
  !  read(10, *) cglat(2)
  !  read(10, *) cglong(2)
  !close(10)

  ! 时间
  insti=50
  doyi = 6
  doy=doyi
  ltstart=3.
  ltend=6.0457     ! 5.7-102.3  6.04-98.07 6.3-94.41 6.0457-98.001
  
!地磁数据
ap=3.          !        0106:3.   0102:8.   3-50
f107=159.       !   0106:159.    0102:142   75-250
f107a=149.3         !    81天平均  149.3     (142.8+153.9+151.3)/3=149.3
f107p=153.           !   0106:153.  0102:136

    idate=year*1000+doyi
    glat=cglat(1)
    glong=cglong(1)
    ut=ltstart*3600.-cglong(1)/15.*3600.+(1-1.)*(ltend-ltstart)*3600./(insti-1.)
    call solzen(idate,ut,glat,glong,sza)
    sza0=sza

do instance=1,insti
  !   pphitop(:)=1.
    ppphitop=0.
  !   uflx(:,:)=0.
    ut=ltstart*3600.-cglong(1)/15.*3600.+(instance-1.)*(ltend-ltstart)*3600./(insti-1.)
    if (ut<0.) ut=ut+24.*3600.
    if (ut>=24.*3600.) ut=ut-24.*3600.
    call make_output_path(year, doyi, glat, glong, ut, instance, filepath)
    write(inst,"(i3)")instance
    open(11,file=filepath,status='replace')



  !do instance=1,insti
  !   pphitop(:)=1.
  !   ppphitop=0.
  !   uflx(:,:)=0.
  !   
  !  write(inst,"(i3)")instance
  !  open(unit=11,&      !输出insti个文件
  !        file='HaiNan20240106/'//trim(adjustl(doyc))//'Module'&
  !                //'/output_vs'//trim(adjustl(inst))//'.txt', status='replace')

    do iii=1,100
       flag=2-mod(iii,2)
!
! Get input values:
!
      write(inst,"(i3)")instance+(flag-1)*insti

        !open(unit=11,&      !输出2insti个文件
        !  file='HaiNan20240106/'//trim(adjustl(doyc))//'Dawn'&
        !          //'/output_vs'//trim(adjustl(inst))//'.txt', status='replace')
      
      idate=year*1000+doy
      glat=cglat(flag)
      glong=cglong(flag)
      ut=ltstart*3600.-cglong(1)/15.*3600.+(instance-1.)*(ltend-ltstart)*3600./(insti-1.)
      if (ut<0.) ut=ut+24.*3600.
      if (ut>=24.*3600.) ut=ut-24.*3600.
      ef=0.
      ec=0.
      call solzen(idate,ut,glat,glong,sza)
      sza=sza*3.1415926536/180.

!
! Calculate local solar time:
!
      stl = ut/3600. + glong/15.
      if (stl < 0.) stl = stl + 24.
      if (stl >= 24.) stl = stl - 24.
!
! Call MZGRID to use MSIS/NOEM/IRI inputs on default altitude grid:
!
      call mzgrid (jmax,nex,idate,ut,glat,glong,stl,f107a,f107,f107p,ap,iri90_dir, &
                   z,zo,zo2,zn2,zns,znd,zno,ztn,zun,zvn,ze,zti,zte,zxden)
!
! Call MAXT to put auroral electron flux specified by namelist input into phitop array:
!
      phitop(:) = 0.
      if (ef>.001 .and. ec>1.) call maxt (ef,ec,ener,del,nbins,itail,fmono,emono,phitop)
      phitop=phitop+uflx(:,jmax)

!
! Fill altitude array, converting to cm:
!
      zz(:) = z(:) * 1.e5     ! km to cm at all jmax levels
!
! Call GLOW to calculate ionized and excited species, airglow emission rates,
! and vertical column brightnesses:
!
      call glow
 
!
! Call CONDUCT to calculate Pederson and Hall conductivities:
!
      do j=1,jmax
        call conduct (glat, glong, z(j), zo(j), zo2(j), zn2(j), &
                      zxden(3,j), zxden(6,j), zxden(7,j), ztn(j), zti(j), zte(j), &
                      pedcond(j), hallcond(j))
      enddo
!
! Output section:
!
      
      
    !  write(11,"(1x,2i7,f9.5,i10,f8.1,8f8.2,i6)")jmax,nbins,sza, &      !iii内
    !        idate,ut,glat,glong,f107a,f107,f107p,ap,ef,ec,iii
    !        !open(15,file = 'C:/Users/31217/Desktop/氦激光雷达/地磁场光电子模型/u500to900.txt')
    !    do j = 1, jmax
    !        !if (z(j) >= 400. .AND. z(j) <= 500.) then
    !        !    write(15, "(1x,f6.1,190e12.5)") z(j),(uflx(ii,j),ii=1,nbins)
    !        !end if
    !        write(11,"(1x,f6.1,190e12.5)") z(j), (uflx(ii,j), ii=1,nbins)
    !    end do
    !        !close(15)
    !        
    !        !open(20,file = 'C:/Users/31217/Desktop/氦激光雷达/地磁场光电子模型/d500to900.txt')
    !    do j = 1, jmax
    !        !if (z(j) >= 400. .AND. z(j) <= 500.) then
    !        !    write(20, "(1x,f6.1,190e12.5)") z(j),(dflx(ii,j),ii=1,nbins)
    !        !end if
    !        write(11,"(1x,f6.1,190e12.5)") z(j), (dflx(ii,j), ii=1,nbins)
    !    end do
    !        !close(20)
    !  write(11,"(1x,190f10.3)") (ener(j),j=1,nbins)
    !  write(11,"(1x,190f10.3)") (del(j),j=1,nbins)
    !close(unit=11)

     if(sum(abs(pphitop-uflx(:,jmax)))<1000. .and. flag==1 .and. &
          ppphitop-sum(abs(pphitop-uflx(:,jmax)))<0.001) exit
     if(flag==1) ppphitop=sum(abs(pphitop-uflx(:,jmax)))
     if(flag==1) pphitop=uflx(:,jmax)
     
     if(sum(abs(pphitop-uflx(:,jmax)))<0.01 .and. flag==1) exit
     if(flag==1) pphitop=uflx(:,jmax)
     
     end do
     
     !if (flag==2) write(*,"(1x,A3,2x,A3)")doyc,inst
        
        write(11,"(1x,2i7,f9.5,i10,f8.1,8f8.2,i6)")jmax,nbins,sza, &
            idate,ut,glat,glong,f107a,f107,f107p,ap,ef,ec,iii
            !open(15,file = 'C:/Users/31217/Desktop/氦激光雷达/地磁场光电子模型/u500to900.txt')
        do j = 1, jmax
            !if (z(j) >= 400. .AND. z(j) <= 500.) then
            !    write(15, "(1x,f6.1,190e12.5)") z(j),(uflx(ii,j),ii=1,nbins)
            !end if
            write(11,"(1x,f6.1,190e12.5)") z(j), (uflx(ii,j), ii=1,nbins)
        end do
            !close(15)
            
            !open(20,file = 'C:/Users/31217/Desktop/氦激光雷达/地磁场光电子模型/d500to900.txt')
        do j = 1, jmax
            !if (z(j) >= 400. .AND. z(j) <= 500.) then
            !    write(20, "(1x,f6.1,190e12.5)") z(j),(dflx(ii,j),ii=1,nbins)
            !end if
            write(11,"(1x,f6.1,190e12.5)") z(j), (dflx(ii,j), ii=1,nbins)
        end do
            !close(20)
      write(11,"(1x,190f10.3)") (ener(j),j=1,nbins)
      write(11,"(1x,190f10.3)") (del(j),j=1,nbins)
    
    write(*,*)instance
    
    call mzgrid_elec (jmax,nex,idate,ut,glat,glong,stl,f107a,f107,f107p,ap,iri90_dir, &
                z,zo,zo2,zn2,zns,znd,zno,ztn,zun,zvn,ze,zti,zte,zxden)
    close(unit=11)
  end do  !txt结束

    !open(21,file = 'C://Users//31217//Desktop//氦激光雷达//地磁场光电子模型//sza.txt')
    !write(21, "(f7.3)") sza*180./3.1415926536
    !close(21)
    ! write(*, "(f7.3)") sza0, sza*180./3.1415926536

stop
end program
