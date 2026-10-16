program main
  use glowbasic_Heion_api, only: run_glowbasic_Heion_args
  use glow_run_helpers, only: date_to_doy
  use sw_daily_lookup, only: get_spaceweather_by_date_ap3h
  implicit none

  real, parameter :: lt_target = 6.0
  real, parameter :: ef_fixed = 0.0
  real, parameter :: ec_fixed = 0.0
  real, parameter :: mlat_fixed = -999.0
  real :: cglat(2), cglong(2)
  real :: ut_target
  real :: ap3h(8), f107, f107a, f107p
  integer :: year, month, day, doyi, insti
  integer :: days_this_month
  logical :: found_sw, ok_date

  year = 2026
  insti = 1

  cglat(1) = 19.5
  cglong(1) = 109.0
  cglat(2) = -7.6
  cglong(2) = 108.73

  do month = 1, 5
    days_this_month = days_in_month(year, month)
    do day = 1, days_this_month
      call date_to_doy(year, month, day, doyi, ok_date)
      if (.not. ok_date) then
        write(*,*) 'Invalid date:', year, month, day
        cycle
      end if

      call get_spaceweather_by_date_ap3h(year, month, day, ap3h, f107, f107a, f107p, found_sw, 'sw_daily_latest.csv')
      if (.not. found_sw) then
        write(*,*) 'No space weather record for date:', year, month, day
        cycle
      end if

      ut_target = modulo(lt_target - cglong(1) / 15.0, 24.0)
      if (ut_target < 0.0) ut_target = ut_target + 24.0

      call run_glowbasic_Heion_args( &
          year, doyi, cglat, cglong, ut_target, ut_target, insti, &
          f107a, f107, f107p, ap3h(1), ef_fixed, ec_fixed, mlat_fixed, ap3h)

      write(*,*) 'Finished date:', year, month, day, ' LT=', lt_target, ' UT=', ut_target
    end do
  end do

contains

  integer function days_in_month(y, m)
    integer, intent(in) :: y, m
    logical :: leap

    leap = (mod(y, 4) == 0 .and. mod(y, 100) /= 0) .or. mod(y, 400) == 0

    select case (m)
    case (1, 3, 5, 7, 8, 10, 12)
      days_in_month = 31
    case (4, 6, 9, 11)
      days_in_month = 30
    case (2)
      if (leap) then
        days_in_month = 29
      else
        days_in_month = 28
      end if
    case default
      days_in_month = 0
    end select
  end function days_in_month

end program main
