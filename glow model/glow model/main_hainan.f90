program main
  use glowbasic_Heion_api, only: run_glowbasic_Heion_args
  use glow_run_helpers, only: date_to_doy, get_efec_for_lat_utseries
  use sw_daily_lookup, only: get_spaceweather_by_date_ap3h
  implicit none

  real :: cglat(2), cglong(2)
  real :: utstart, utend
  real :: ap3h(8), f107, f107a, f107p
  real, allocatable :: ef_series(:), ec_series(:)
  real, allocatable :: mlat_series(:), mlon_series(:), mlt_series(:)
  integer :: year, month, day, doyi, insti
  logical :: found_sw, found_efec, ok_date
  character(len=8) :: date_stamp
  character(len=256) :: efec_file

  ! Fixed one-day run configuration (input by calendar date).
  year = 2025
  month = 10
  day = 10
  utstart = 0.0
  utend = 23.9166667
  insti = 288
  cglat(1) = 19.5
  cglong(1) = 109.0
  cglat(2) = -5.66
  cglong(2) = 109.24
  
   !cglat(1) = 78.1
   !cglong(1) = 16
  ! cglat(2) = -68.7
  ! cglong(2) = 82

  call date_to_doy(year, month, day, doyi, ok_date)
  if (.not. ok_date) then
    write(*,*) 'Invalid date:', year, month, day
    stop 1
  end if

  write(date_stamp, '(I4.4,I2.2,I2.2)') year, month, day
  efec_file = '../ZhangPaxtonAuroralModel-main/output/north/' // trim(date_stamp) // '_24h.txt'

  allocate(ef_series(insti), ec_series(insti))
  allocate(mlat_series(insti), mlon_series(insti), mlt_series(insti))

  call get_spaceweather_by_date_ap3h(year, month, day, ap3h, f107, f107a, f107p, found_sw, 'sw_daily_latest.csv')
  if (.not. found_sw) then
    write(*,*) 'No space weather record for date:', year, month, day
    stop 1
  end if

  call get_efec_for_lat_utseries(efec_file, cglat(1), utstart, utend, insti, ef_series, ec_series, found_efec, &
                                 mlat_series, mlon_series, mlt_series, cglong(1))
  if (.not. found_efec) then
    write(*,*) 'No ef/ec time-lat match for lat=', cglat(1), ' in ', trim(efec_file)
   ! stop 1
  end if
  
ef_series = 0
ec_series = 0
!ap3h=300

  call run_glowbasic_Heion_args( &
      year, doyi, cglat, cglong, utstart, utend, insti, &
      f107a, f107, f107p, ap3h(1), ef_series(1), ec_series(1), mlat_series(1), &
      ap3h, ef_series, ec_series, mlat_series)

  if (allocated(mlt_series)) deallocate(mlt_series)
  if (allocated(mlon_series)) deallocate(mlon_series)
  if (allocated(mlat_series)) deallocate(mlat_series)
  if (allocated(ec_series)) deallocate(ec_series)
  if (allocated(ef_series)) deallocate(ef_series)

end program main
