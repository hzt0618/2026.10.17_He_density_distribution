module glow_run_helpers
  implicit none

contains

  subroutine date_to_doy(y, m, d, doy_out, ok)
    integer, intent(in) :: y, m, d
    integer, intent(out) :: doy_out
    logical, intent(out) :: ok

    integer :: month_days(12)
    integer :: i
    logical :: leap

    month_days = (/31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/)
    leap = (mod(y, 4) == 0 .and. mod(y, 100) /= 0) .or. mod(y, 400) == 0
    if (leap) month_days(2) = 29

    ok = .false.
    doy_out = -1

    if (m < 1 .or. m > 12) return
    if (d < 1 .or. d > month_days(m)) return

    doy_out = d
    do i = 1, m - 1
      doy_out = doy_out + month_days(i)
    end do

    ok = .true.
  end subroutine date_to_doy

  subroutine parse_conjugate_line(row, src_lat, src_lon, conj_lat, conj_lon, ok)
    character(len=*), intent(in) :: row
    real, intent(out) :: src_lat, src_lon, conj_lat, conj_lon
    logical, intent(out) :: ok

    integer :: ios_parse, idx_i
    real :: input_lat_i, input_lon_i, input_x_i, input_y_i, input_z_i
    real :: conj_lat_i, conj_lon_i, conj_alt_km_i, b_ratio_i

    ok = .false.
    src_lat = 0.0
    src_lon = 0.0
    conj_lat = 0.0
    conj_lon = 0.0

    ! Preferred CSV format:
    ! idx,input_lat,input_lon,input_x,input_y,input_z,conj_lat,conj_lon,conj_alt_km,B_ratio
    read(row, *, iostat=ios_parse) idx_i, input_lat_i, input_lon_i, &
        input_x_i, input_y_i, input_z_i, conj_lat_i, conj_lon_i, conj_alt_km_i, b_ratio_i
    if (ios_parse == 0) then
      src_lat = input_lat_i
      src_lon = input_lon_i
      conj_lat = conj_lat_i
      conj_lon = conj_lon_i
      ok = .true.
      return
    end if

    ! Legacy CSV format:
    ! idx,input_lat,input_lon,conj_lat,conj_lon,conj_alt_km,B_ratio
    read(row, *, iostat=ios_parse) idx_i, input_lat_i, input_lon_i, &
        conj_lat_i, conj_lon_i, conj_alt_km_i, b_ratio_i
    if (ios_parse == 0) then
      src_lat = input_lat_i
      src_lon = input_lon_i
      conj_lat = conj_lat_i
      conj_lon = conj_lon_i
      ok = .true.
      return
    end if

    ! Backward-compatible format:
    ! src_lat,src_lon,conj_lat,conj_lon
    read(row, *, iostat=ios_parse) src_lat, src_lon, conj_lat, conj_lon
    if (ios_parse == 0) ok = .true.
  end subroutine parse_conjugate_line

  subroutine get_efec_for_lat_utseries(file_path, target_lat, utstart, utend, insti, ef_series, ec_series, found, &
                                       mlat_series, mlon_series, mlt_series, target_lon)
    character(len=*), intent(in) :: file_path
    real, intent(in) :: target_lat, utstart, utend
    integer, intent(in) :: insti
    real, intent(out) :: ef_series(insti), ec_series(insti)
    logical, intent(out) :: found
    real, intent(out), optional :: mlat_series(insti), mlon_series(insti), mlt_series(insti)
    real, intent(in), optional :: target_lon

    integer :: unit_no, ios_local, i, idx_i
    real :: kp_i, lat_i, lon_i, alt_i, mlat_i, mlon_i, mlt_i, ef_i, ec_i
    real :: ut_row, ut_target, dlat, dlon, dut
    character(len=256) :: row
    character(len=32) :: time_utc_i
    character(len=64) :: status_i
    logical :: time_ok
    logical :: found_i(insti)
    real :: best_dut(insti), best_dlat(insti), best_dlon(insti)
    real, parameter :: eps_tie = 1.0e-6
    real, parameter :: lon_tol = 1.0

    found = .false.
    ef_series = 0.0
    ec_series = 0.0
    if (present(mlat_series)) mlat_series = 0.0
    if (present(mlon_series)) mlon_series = 0.0
    if (present(mlt_series)) mlt_series = 0.0
    found_i = .false.
    best_dut = 1.0e30
    best_dlat = 1.0e30
    best_dlon = 1.0e30

    open(newunit=unit_no, file=trim(file_path), status='old', action='read', iostat=ios_local)
    if (ios_local /= 0) return

    do
      read(unit_no, '(A)', iostat=ios_local) row
      if (ios_local /= 0) exit
      if (len_trim(row) == 0) cycle
      if (row(1:1) == '#') cycle
      if (index(row, 'time_utc,kp,lat,lon,alt_km') == 1) cycle
      if (index(row, 'time_utc,kp,idx,lat,lon,alt_km') == 1) cycle

      read(row, *, iostat=ios_local) time_utc_i, kp_i, idx_i, lat_i, lon_i, alt_i, mlat_i, mlon_i, mlt_i, ef_i, ec_i, status_i
      if (ios_local /= 0) then
        read(row, *, iostat=ios_local) time_utc_i, kp_i, lat_i, lon_i, alt_i, mlat_i, mlon_i, mlt_i, ef_i, ec_i, status_i
        if (ios_local /= 0) cycle
      end if
      if (index(adjustl(status_i), 'ok') /= 1) cycle

      call parse_iso_utc_hour(time_utc_i, ut_row, time_ok)
      if (.not. time_ok) cycle

      dlat = abs(lat_i - target_lat)
      if (dlat > 0.26) cycle
      if (present(target_lon)) then
        dlon = lon_cyclic_diff(lon_i, target_lon)
        if (dlon > lon_tol) cycle
      else
        dlon = 0.0
      end if

      do i = 1, insti
        ut_target = ut_for_instance(utstart, utend, insti, i)
        dut = hour_cyclic_diff(ut_row, ut_target)
        if ((dut < best_dut(i)) .or. &
            (abs(dut - best_dut(i)) <= eps_tie .and. dlat < best_dlat(i)) .or. &
            (abs(dut - best_dut(i)) <= eps_tie .and. abs(dlat - best_dlat(i)) <= eps_tie .and. dlon < best_dlon(i))) then
          best_dut(i) = dut
          best_dlat(i) = dlat
          best_dlon(i) = dlon
          ef_series(i) = ef_i
          ec_series(i) = ec_i
          if (present(mlat_series)) mlat_series(i) = mlat_i
          if (present(mlon_series)) mlon_series(i) = mlon_i
          if (present(mlt_series)) mlt_series(i) = mlt_i
          found_i(i) = .true.
        end if
      end do
    end do

    close(unit_no)

    found = all(found_i)
    if (.not. found) then
      ef_series = 0.0
      ec_series = 0.0
      if (present(mlat_series)) mlat_series = 0.0
      if (present(mlon_series)) mlon_series = 0.0
      if (present(mlt_series)) mlt_series = 0.0
    end if
  end subroutine get_efec_for_lat_utseries

  subroutine parse_iso_utc_hour(time_utc, ut_hour, ok)
    character(len=*), intent(in) :: time_utc
    real, intent(out) :: ut_hour
    logical, intent(out) :: ok

    integer :: hh, mm, ss, ios_t

    ut_hour = 0.0
    ok = .false.
    if (len_trim(time_utc) < 19) return

    read(time_utc(12:13), '(i2)', iostat=ios_t) hh
    if (ios_t /= 0) return
    read(time_utc(15:16), '(i2)', iostat=ios_t) mm
    if (ios_t /= 0) return
    read(time_utc(18:19), '(i2)', iostat=ios_t) ss
    if (ios_t /= 0) return

    ut_hour = real(hh) + real(mm)/60.0 + real(ss)/3600.0
    ok = .true.
  end subroutine parse_iso_utc_hour

  real function ut_for_instance(utstart, utend, insti, instance) result(ut_hour)
    real, intent(in) :: utstart, utend
    integer, intent(in) :: insti, instance
    integer :: step_denom

    step_denom = max(1, insti - 1)
    ut_hour = utstart + (real(instance - 1) * (utend - utstart) / real(step_denom))
    ut_hour = modulo(ut_hour, 24.0)
    if (ut_hour < 0.0) ut_hour = ut_hour + 24.0
  end function ut_for_instance

  real function hour_cyclic_diff(h1, h2) result(dh)
    real, intent(in) :: h1, h2
    dh = abs(modulo((h1 - h2) + 12.0, 24.0) - 12.0)
  end function hour_cyclic_diff

  real function lon_cyclic_diff(lon1, lon2) result(dlon)
    real, intent(in) :: lon1, lon2
    dlon = abs(modulo((lon1 - lon2) + 180.0, 360.0) - 180.0)
  end function lon_cyclic_diff

end module glow_run_helpers
