program test_ap
  use sw_daily_lookup, only: get_spaceweather_by_date
  implicit none

  character(len=256) :: csv_path

  csv_path = 'sw_daily_latest.csv'

  print *, 'CSV path = ', trim(csv_path)
  print *, 'Case A: 1957-10-01, different UT (check Ap slot mapping)'
  call run_case(1957, 10, 1, 1.0, csv_path)    ! Ap0 expected from sample row
  call run_case(1957, 10, 1, 4.0, csv_path)    ! Ap3 expected from sample row
  call run_case(1957, 10, 1, 22.0, csv_path)   ! Ap21 expected from sample row

  print *, 'Case B: date not in file'
  call run_case(1900, 1, 1, 12.0, csv_path)

contains

  subroutine run_case(year, month, day, ut, csv_path)
    integer, intent(in) :: year, month, day
    real, intent(in) :: ut
    character(len=*), intent(in) :: csv_path
    real :: ap, f107, f107a, f107p
    logical :: found

    call get_spaceweather_by_date(year, month, day, ut, ap, f107, f107a, f107p, found, csv_path)

    write(*,'(a,i4.4,a,i2.2,a,i2.2,a,f6.2,a,l1)') &
      'date=', year, '-', month, '-', day, ', ut=', ut, ', found=', found
    write(*,'(a,f10.3,2x,a,f10.3,2x,a,f10.3,2x,a,f10.3)') &
      'ap=', ap, 'f107=', f107, 'f107a=', f107a, 'f107p=', f107p
    print *, '------------------------------------------------------------'
  end subroutine run_case

end program test_ap
