module sw_daily_lookup
  implicit none
  private
  public :: get_spaceweather_by_date, get_spaceweather_by_date_ap3h

contains

  subroutine get_spaceweather_by_date(year, month, day, ut, ap, f107, f107a, f107p, found, csv_path)
    integer, intent(in) :: year, month, day
    real, intent(in) :: ut
    real, intent(out) :: ap, f107, f107a, f107p
    logical, intent(out) :: found
    character(len=*), intent(in), optional :: csv_path

    real :: ap3h(8)
    integer :: ap_slot

    ap = -1.0
    call get_spaceweather_by_date_ap3h(year, month, day, ap3h, f107, f107a, f107p, found, csv_path)
    if (.not. found) return

    ap_slot = ap_slot_from_ut(ut)
    ap = ap3h(ap_slot)
  end subroutine get_spaceweather_by_date

  subroutine get_spaceweather_by_date_ap3h(year, month, day, ap3h, f107, f107a, f107p, found, csv_path)
    integer, intent(in) :: year, month, day
    real, intent(out) :: ap3h(8), f107, f107a, f107p
    logical, intent(out) :: found
    character(len=*), intent(in), optional :: csv_path

    integer :: unit_no, ios, i
    integer :: y, m, d
    character(len=1024) :: line
    character(len=256) :: field
    character(len=512) :: path_to_use

    found = .false.
    ap3h = -1.0
    f107 = -1.0
    f107a = -1.0
    f107p = -1.0

    if (present(csv_path)) then
      path_to_use = trim(csv_path)
    else
      path_to_use = 'sw_daily_latest.csv'
    end if

    open(newunit=unit_no, file=trim(path_to_use), status='old', action='read', iostat=ios)
    if (ios /= 0) return

    ! Skip header line.
    read(unit_no, '(A)', iostat=ios) line
    if (ios /= 0) then
      close(unit_no)
      return
    end if

    do
      read(unit_no, '(A)', iostat=ios) line
      if (ios /= 0) exit

      call get_csv_field(line, 2, field)
      if (.not. parse_int(field, y)) cycle
      call get_csv_field(line, 3, field)
      if (.not. parse_int(field, m)) cycle
      call get_csv_field(line, 4, field)
      if (.not. parse_int(field, d)) cycle

      if (y == year .and. m == month .and. d == day) then
        do i = 1, 8
          call get_csv_field(line, 15 + i, field) ! Ap0..Ap21 -> columns 16..23
          if (.not. parse_real(field, ap3h(i))) ap3h(i) = -1.0
        end do

        call get_csv_field(line, 32, field) ! f107_obs
        if (.not. parse_real(field, f107)) f107 = -1.0

        call get_csv_field(line, 33, field) ! f107_81ctr_obs
        if (.not. parse_real(field, f107a)) f107a = -1.0

        call get_csv_field(line, 34, field) ! f107_81lst_obs
        if (.not. parse_real(field, f107p)) f107p = -1.0

        found = .true.
        exit
      end if
    end do

    close(unit_no)
  end subroutine get_spaceweather_by_date_ap3h

  integer function ap_slot_from_ut(ut)
    real, intent(in) :: ut
    real :: ut_norm

    ut_norm = modulo(ut, 24.0)
    ap_slot_from_ut = int(ut_norm / 3.0) + 1
    if (ap_slot_from_ut < 1) ap_slot_from_ut = 1
    if (ap_slot_from_ut > 8) ap_slot_from_ut = 8
  end function ap_slot_from_ut

  subroutine get_csv_field(line, field_index, field_value)
    character(len=*), intent(in) :: line
    integer, intent(in) :: field_index
    character(len=*), intent(out) :: field_value

    integer :: i, n, start_pos, end_pos, current_field

    n = len_trim(line)
    field_value = ''
    current_field = 1
    start_pos = 1

    if (field_index <= 0) return

    do i = 1, n + 1
      if (i == n + 1 .or. line(i:i) == ',') then
        end_pos = i - 1
        if (current_field == field_index) then
          if (end_pos >= start_pos) then
            field_value = adjustl(line(start_pos:end_pos))
          else
            field_value = ''
          end if
          return
        end if
        current_field = current_field + 1
        start_pos = i + 1
      end if
    end do
  end subroutine get_csv_field

  logical function parse_int(text, value)
    character(len=*), intent(in) :: text
    integer, intent(out) :: value
    integer :: ios

    value = 0
    read(text, *, iostat=ios) value
    parse_int = (ios == 0)
  end function parse_int

  logical function parse_real(text, value)
    character(len=*), intent(in) :: text
    real, intent(out) :: value
    integer :: ios

    value = -1.0
    read(text, *, iostat=ios) value
    parse_real = (ios == 0)
  end function parse_real

end module sw_daily_lookup
