program main
  use glowbasic_Heion_api, only: run_glowbasic_Heion_args
  use glow_run_helpers, only: parse_conjugate_line, get_efec_for_lat_utseries
  use sw_daily_lookup, only: get_spaceweather_by_date_ap3h
  use date_utils, only: doy_to_date
  use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
  implicit none

  real :: cglat(2), cglong(2)
  real :: utstart, utend
  real :: src_lat, src_lon, conj_lat, conj_lon
  real :: ap3h(8), f107, f107a, f107p
  real, allocatable :: ef_series(:), ec_series(:)
  real, allocatable :: mlat_series(:), mlon_series(:), mlt_series(:)
  integer :: ios, u, year, doyi, month, day, insti
  character(len=256) :: line
  character(len=8) :: cdate
  logical :: found_sw, found_efec, ok_conj
  character(len=*), parameter :: conj_file = '../../IGRF/output/1/conjugate_morning.txt'
  character(len=*), parameter :: efec_file = '../ZhangPaxtonAuroralModel-main/output/20251009_24h.txt'

  year = 2025
  doyi = 282

  ! Keep original main.f90 window setup.
  utstart = 11.0
  utend = 12.0
  insti = 10

  allocate(ef_series(insti), ec_series(insti))
  allocate(mlat_series(insti), mlon_series(insti), mlt_series(insti))

  call doy_to_date(year, doyi, cdate)
  read(cdate(5:6), '(i2)') month
  read(cdate(7:8), '(i2)') day

  open(newunit=u, file=conj_file, status='old', action='read', iostat=ios)
  if (ios /= 0) then
    write(*,*) 'Failed to open file: ', trim(conj_file)
    stop 1
  end if

  do
    read(u, '(A)', iostat=ios) line
    if (ios /= 0) exit

    if (len_trim(line) == 0) cycle
    if (line(1:1) == '#') cycle
    if (index(adjustl(line), 'section=') == 1) cycle
    if (index(adjustl(line), 'columns:') == 1) cycle

    call parse_conjugate_line(line, src_lat, src_lon, conj_lat, conj_lon, ok_conj)
    if (.not. ok_conj) cycle
    if (ieee_is_nan(src_lat) .or. ieee_is_nan(src_lon) .or. ieee_is_nan(conj_lat) .or. ieee_is_nan(conj_lon)) cycle

    cglat = [src_lat, conj_lat]
    cglong = [src_lon, conj_lon]

    call get_spaceweather_by_date_ap3h(year, month, day, ap3h, f107, f107a, f107p, found_sw, 'sw_daily_latest.csv')
    if (.not. found_sw) then
      write(*,*) 'No space weather record for date:', year, month, day
      cycle
    end if

    call get_efec_for_lat_utseries(efec_file, src_lat, utstart, utend, insti, ef_series, ec_series, found_efec, &
                                   mlat_series, mlon_series, mlt_series, src_lon)
    if (.not. found_efec) then
      write(*,*) 'No ef/ec time-lat match for src_lat=', src_lat, ' in ', trim(efec_file)
      cycle
    end if

    call run_glowbasic_Heion_args( &
        year, doyi, cglat, cglong, utstart, utend, insti, &
        f107a, f107, f107p, ap3h(1), ef_series(1), ec_series(1), mlat_series(1), &
        ap3h, ef_series, ec_series, mlat_series)
  end do

  close(u)

  if (allocated(mlt_series)) deallocate(mlt_series)
  if (allocated(mlon_series)) deallocate(mlon_series)
  if (allocated(mlat_series)) deallocate(mlat_series)
  if (allocated(ec_series)) deallocate(ec_series)
  if (allocated(ef_series)) deallocate(ef_series)

end program main
