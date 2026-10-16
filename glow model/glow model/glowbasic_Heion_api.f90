module glowbasic_Heion_api
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  logical, save :: glow_runtime_initialized = .false.

  type :: glowbasic_Heion_config
    integer :: year
    integer :: doyi
    real :: cglat(2)
    real :: cglong(2)
    real :: utstart
    real :: utend
    integer :: insti
    real :: f107a
    real :: f107
    real :: f107p
    real :: ap
    real :: ap3h(8)
    logical :: has_ap3h
    real :: ef
    real :: ec
    real :: mlat
    real, allocatable :: ef_series(:)
    real, allocatable :: ec_series(:)
    real, allocatable :: mlat_series(:)
    logical :: has_efec_series
    logical :: has_mlat_series
  end type glowbasic_Heion_config

contains

  subroutine run_glowbasic_Heion(cfg)
    use cglow,only: jmax,nbins,nex,idate,ut,glat,glong,f107a,f107,f107p,ap,ef,ec,dip
    use cglow,only: iscale,jlocal,kchem,xuvfac,sza
    use cglow,only: zz,zo,zn2,zo2,zns,znd,zno,ztn,ze,zti,zte
    use cglow,only: ener,del,phitop,uflx,dflx,zxden
    use cglow,only: cglow_init,data_dir
    use write_path

    type(glowbasic_Heion_config), intent(in) :: cfg

    character(len=1024) :: iri90_dir
    character(len=256) :: filepath, meta_filepath, phitop_filepath
    real,allocatable :: z(:), zun(:), zvn(:), pedcond(:), hallcond(:), pphitop(:), feedback_top(:), feedback_prev(:,:)
    real,allocatable :: uflx_prev_valid(:,:), dflx_prev_valid(:,:)
    real,allocatable :: phitop_direct_last(:), phitop_feedback_last(:), phitop_total_last(:)
    real :: stl, fmono, emono, ppphitop, ut_file, mlat_current, coupling
    real :: xf, yf, zf, ff, dec, sdip, raw_dip_deg
    real :: feedback_diff, feedback_norm
    real,parameter :: feedback_relax=0.50
    real,parameter :: output_flux_scale=1.00
    integer :: j, ii, iii, flag, instance, step_denom, itail, loop, next_flag
    logical :: feedback_seen(2), flux_valid_seen

    iscale = 1
    xuvfac = 3.
    kchem = 4
    jlocal = 0
    itail = 0
    fmono = 0.
    emono = 0.

    data_dir = 'data/'
    iri90_dir = 'data/iri90/'
    jmax = 138

    allocate(z(jmax), zun(jmax), zvn(jmax))
    allocate(pedcond(jmax), hallcond(jmax))
    allocate(pphitop(nbins))
    allocate(feedback_top(nbins))
    allocate(feedback_prev(2, nbins))
    allocate(uflx_prev_valid(nbins, jmax))
    allocate(dflx_prev_valid(nbins, jmax))
    allocate(phitop_direct_last(nbins))
    allocate(phitop_feedback_last(nbins))
    allocate(phitop_total_last(nbins))

    if (.not. glow_runtime_initialized) then
      call cglow_init
      call egrid(ener, del, nbins)
      glow_runtime_initialized = .true.
    end if

    idate = cfg%year*1000 + cfg%doyi
    f107a = cfg%f107a
    f107 = cfg%f107
    f107p = cfg%f107p
    ef = cfg%ef
    ec = cfg%ec
    mlat_current = cfg%mlat

    step_denom = max(1, cfg%insti - 1)

    do instance = 1, cfg%insti
      ut_file = cfg%utstart*3600. + &
          (instance-1.)*(cfg%utend-cfg%utstart)*3600./real(step_denom)

      if (ut_file < 0.) ut_file = ut_file + 24.*3600.
      if (ut_file >= 24.*3600.) ut_file = ut_file - 24.*3600.
      if (cfg%has_ap3h) then
        ap = cfg%ap3h(ap_slot_from_ut_seconds(ut_file))
      else
        ap = cfg%ap
      end if
      if (cfg%has_efec_series) then
        ef = cfg%ef_series(instance)
        ec = cfg%ec_series(instance)
      else
        ef = cfg%ef
        ec = cfg%ec
      end if
      if (cfg%has_mlat_series) then
        mlat_current = cfg%mlat_series(instance)
      else
        mlat_current = cfg%mlat
      end if
      coupling = 0.
      loop = 100

      call make_output_path(cfg%year, cfg%doyi, cfg%cglat(1), cfg%cglong(1), ut_file, instance, filepath)
      open(11, file=filepath, status='replace')

      pphitop = 0.
      ppphitop = 0.
      feedback_diff = 0.
      feedback_norm = 1.
      raw_dip_deg = 0.
      pphitop(:)=1.
      feedback_top(:)=0.
      feedback_prev(:,:)=0.
      phitop_direct_last(:)=0.
      phitop_feedback_last(:)=0.
      phitop_total_last(:)=0.
      feedback_seen(:)=.false.
      flux_valid_seen = .false.
      uflx_prev_valid(:,:) = 0.
      dflx_prev_valid(:,:) = 0.
      uflx(:,:)=0.
      dflx(:,:)=0.

      do iii = 1, loop
        flag = 2 - mod(iii, 2)

        glat = cfg%cglat(flag)
        glong = cfg%cglong(flag)

        ut = ut_file
        call solzen(idate, ut, glat, glong, sza)
        sza = sza*3.1415926536/180.

        stl = ut/3600. + glong/15.
        if (stl < 0.) stl = stl + 24.
        if (stl >= 24.) stl = stl - 24.
        if (mlat_current < -999.0) stl = stl

        call mzgrid(jmax, nex, idate, ut, glat, glong, stl, f107a, f107, f107p, ap, iri90_dir, &
                    z, zo, zo2, zn2, zns, znd, zno, ztn, zun, zvn, ze, zti, zte, zxden)
        call fieldm(glat, glong, 300., xf, yf, zf, ff, raw_dip_deg, dec, sdip)
        coupling = conjugate_coupling(dip)

        phitop(:) = 0.
        if (ef > .001 .and. ec > 0.) call maxt(ef, ec, ener, del, nbins, itail, fmono, emono, phitop)
        if (flag == 1) phitop_direct_last = phitop
        feedback_top = safe_feedback_spectrum(coupling * uflx(:, jmax), feedback_prev(flag, :), feedback_seen(flag))
        if (feedback_seen(flag)) then
          feedback_top = feedback_relax * feedback_top + (1.0 - feedback_relax) * feedback_prev(flag, :)
        end if
        if (flag == 1) phitop_feedback_last = feedback_top
        phitop = phitop + feedback_top
        if (flag == 1) phitop_total_last = phitop

        zz(:) = z(:)*1.e5
        call glow
        if (all(ieee_is_finite(uflx)) .and. all(ieee_is_finite(dflx))) then
          uflx_prev_valid(:,:) = uflx(:,:)
          dflx_prev_valid(:,:) = dflx(:,:)
          flux_valid_seen = .true.
        else if (flux_valid_seen) then
          uflx(:,:) = uflx_prev_valid(:,:)
          dflx(:,:) = dflx_prev_valid(:,:)
        else
          uflx(:,:) = 0.
          dflx(:,:) = 0.
        end if

        do j = 1, jmax
          call conduct(glat, glong, z(j), zo(j), zo2(j), zn2(j), &
                       zxden(3,j), zxden(6,j), zxden(7,j), ztn(j), zti(j), zte(j), &
                       pedcond(j), hallcond(j))
        end do

        next_flag = 3 - flag
        feedback_top = safe_feedback_spectrum(coupling * uflx(:, jmax), feedback_prev(next_flag, :), feedback_seen(next_flag))
        if (feedback_seen(next_flag)) then
          feedback_top = feedback_relax * feedback_top + (1.0 - feedback_relax) * feedback_prev(next_flag, :)
        end if
        feedback_diff = sum(abs(feedback_prev(next_flag, :) - feedback_top))
        feedback_norm = max(sum(abs(feedback_prev(next_flag, :))), 1.0)

        if (flag == 1 .and. feedback_seen(next_flag)) then
          if (feedback_diff < 1000. .and. abs(ppphitop - feedback_diff) < 0.001) exit
          if ((feedback_diff / feedback_norm) < 1.0e-3) exit
        end if
        if (flag == 1) ppphitop = feedback_diff
      !  if (flag == 1) pphitop = uflx(:,jmax)
        
        feedback_prev(next_flag, :) = feedback_top
        feedback_seen(next_flag) = .true.
        if (flag == 1) pphitop = feedback_top
      end do

      write(11,"(1x,2i7,f9.5,i10,f8.1,8f8.2,i6)") jmax, nbins, sza, &
          idate, ut, glat, glong, f107a, f107, f107p, ap, ef, ec, iii
      uflx(:,:) = output_flux_scale * uflx(:,:)
      dflx(:,:) = output_flux_scale * dflx(:,:)
      do j = 1, jmax
        write(11,"(1x,f6.1,190e12.5)") z(j), (uflx(ii,j), ii=1,nbins)
      end do
      do j = 1, jmax
        write(11,"(1x,f6.1,190e12.5)") z(j), (dflx(ii,j), ii=1,nbins)
      end do
      write(11,"(1x,190f10.3)") (ener(j), j=1, nbins)
      write(11,"(1x,190f10.3)") (del(j), j=1, nbins)

      call mzgrid_elec(jmax, nex, idate, ut, glat, glong, stl, f107a, f107, f107p, ap, iri90_dir, &
                       z, zo, zo2, zn2, zns, znd, zno, ztn, zun, zvn, ze, zti, zte, zxden)

      close(11)

      meta_filepath = trim(filepath)
      if (len_trim(meta_filepath) > 4) then
        meta_filepath = trim(meta_filepath(:len_trim(meta_filepath)-4))//'_meta.txt'
      else
        meta_filepath = trim(meta_filepath)//'_meta.txt'
      end if

      open(12, file=meta_filepath, status='replace')
      write(12,"(A)") '# GLOW sidecar metadata'
      write(12,"(A,F12.6)") 'dip_rad=', dip
      write(12,"(A,F12.6)") 'dip_deg=', dip*180.0/3.1415926536
      write(12,"(A,F12.6)") 'raw_dip_deg=', raw_dip_deg
      write(12,"(A,F12.6)") 'coupling=', coupling
      write(12,"(A,F12.6)") 'output_flux_scale=', output_flux_scale
      write(12,"(A,F12.6)") 'feedback_relax=', feedback_relax
      write(12,"(A,E15.6)") 'feedback_diff=', feedback_diff
      write(12,"(A,E15.6)") 'feedback_norm=', feedback_norm
      write(12,"(A,E15.6)") 'phitop_direct_sum=', sum(phitop_direct_last)
      write(12,"(A,E15.6)") 'phitop_feedback_sum=', sum(phitop_feedback_last)
      write(12,"(A,E15.6)") 'phitop_total_sum=', sum(phitop_total_last)
      write(12,"(A,I0)") 'iii=', iii
      write(12,"(A,F12.6)") 'src_lat=', cfg%cglat(1)
      write(12,"(A,F12.6)") 'src_lon=', cfg%cglong(1)
      write(12,"(A,F12.6)") 'mlat_current=', mlat_current
      write(12,"(A,F12.6)") 'conj_lat=', cfg%cglat(2)
      write(12,"(A,F12.6)") 'conj_lon=', cfg%cglong(2)
      write(12,"(A,F12.6)") 'glat_written=', glat
      write(12,"(A,F12.6)") 'glong_written=', glong
      write(12,"(A,F12.6)") 'ut_seconds=', ut
      close(12)

      phitop_filepath = trim(filepath)
      if (len_trim(phitop_filepath) > 4) then
        phitop_filepath = trim(phitop_filepath(:len_trim(phitop_filepath)-4))//'_phitop.csv'
      else
        phitop_filepath = trim(phitop_filepath)//'_phitop.csv'
      end if
      open(13, file=phitop_filepath, status='replace')
      write(13,"(A)") 'energy_ev,delta_ev,phitop_direct,phitop_feedback,phitop_total'
      do j = 1, nbins
        write(13,"(E15.6,',',E15.6,',',E15.6,',',E15.6,',',E15.6)") &
            ener(j), del(j), phitop_direct_last(j), phitop_feedback_last(j), phitop_total_last(j)
      end do
      close(13)
    end do

    if (allocated(pphitop)) deallocate(pphitop)
    if (allocated(feedback_top)) deallocate(feedback_top)
    if (allocated(feedback_prev)) deallocate(feedback_prev)
    if (allocated(uflx_prev_valid)) deallocate(uflx_prev_valid)
    if (allocated(dflx_prev_valid)) deallocate(dflx_prev_valid)
    if (allocated(phitop_total_last)) deallocate(phitop_total_last)
    if (allocated(phitop_feedback_last)) deallocate(phitop_feedback_last)
    if (allocated(phitop_direct_last)) deallocate(phitop_direct_last)
    if (allocated(hallcond)) deallocate(hallcond)
    if (allocated(pedcond)) deallocate(pedcond)
    if (allocated(zvn)) deallocate(zvn)
    if (allocated(zun)) deallocate(zun)
    if (allocated(z)) deallocate(z)
  end subroutine run_glowbasic_Heion

  subroutine run_glowbasic_Heion_args(year, doyi, cglat, cglong, utstart, utend, insti, &
                                      f107a, f107, f107p, ap, ef, ec, mlat, ap3h, ef_series, ec_series, mlat_series)
    integer, intent(in) :: year, doyi, insti
    real, intent(in) :: cglat(2), cglong(2)
    real, intent(in) :: utstart, utend
    real, intent(in) :: f107a, f107, f107p, ap, ef, ec, mlat
    real, intent(in), optional :: ap3h(8)
    real, intent(in), optional :: ef_series(:), ec_series(:)
    real, intent(in), optional :: mlat_series(:)
    type(glowbasic_Heion_config) :: cfg

    cfg%year = year
    cfg%doyi = doyi
    cfg%cglat = cglat
    cfg%cglong = cglong
    cfg%utstart = utstart
    cfg%utend = utend
    cfg%insti = insti
    cfg%f107a = f107a
    cfg%f107 = f107
    cfg%f107p = f107p
    cfg%ap = ap
    cfg%ap3h = -1.0
    cfg%has_ap3h = .false.
    if (present(ap3h)) then
      cfg%ap3h = ap3h
      cfg%has_ap3h = .true.
    end if
    cfg%ef = ef
    cfg%ec = ec
    cfg%mlat = mlat
    cfg%has_efec_series = .false.
    cfg%has_mlat_series = .false.
    if (present(ef_series) .and. present(ec_series)) then
      if (size(ef_series) == insti .and. size(ec_series) == insti) then
        allocate(cfg%ef_series(insti), cfg%ec_series(insti))
        cfg%ef_series = ef_series
        cfg%ec_series = ec_series
        cfg%has_efec_series = .true.
      end if
    end if
    if (present(mlat_series)) then
      if (size(mlat_series) == insti) then
        allocate(cfg%mlat_series(insti))
        cfg%mlat_series = mlat_series
        cfg%has_mlat_series = .true.
      end if
    end if
    call run_glowbasic_Heion(cfg)
  end subroutine run_glowbasic_Heion_args

  integer function ap_slot_from_ut_seconds(ut_sec)
    real, intent(in) :: ut_sec
    real :: ut_hour

    ut_hour = modulo(ut_sec / 3600.0, 24.0)
    ap_slot_from_ut_seconds = int(ut_hour / 3.0) + 1
    if (ap_slot_from_ut_seconds < 1) ap_slot_from_ut_seconds = 1
    if (ap_slot_from_ut_seconds > 8) ap_slot_from_ut_seconds = 8
  end function ap_slot_from_ut_seconds

  real function conjugate_coupling(dip_rad)
    real, intent(in) :: dip_rad
    real, parameter :: fixed_coupling = 0.5

    conjugate_coupling = fixed_coupling + 0.0 * dip_rad
  end function conjugate_coupling

  function safe_feedback_spectrum(raw_feedback, prev_feedback, has_prev) result(feedback_out)
    real, intent(in) :: raw_feedback(:)
    real, intent(in) :: prev_feedback(:)
    logical, intent(in) :: has_prev
    real :: feedback_out(size(raw_feedback))
    logical :: bad_mask(size(raw_feedback))

    feedback_out = raw_feedback
    bad_mask = .not. ieee_is_finite(feedback_out)
    if (any(bad_mask)) then
      if (has_prev) then
        feedback_out = merge(prev_feedback, feedback_out, bad_mask)
      else
        feedback_out = merge(0.0, feedback_out, bad_mask)
      end if
    end if
  end function safe_feedback_spectrum

end module glowbasic_Heion_api
