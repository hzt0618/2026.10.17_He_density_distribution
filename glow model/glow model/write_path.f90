module write_path

implicit none

contains

!==============================
! DOY -> YYYYMMDD
!==============================
subroutine doy_to_date(year, doy, cdate)

integer, intent(in) :: year, doy
character(len=8), intent(out) :: cdate

integer :: month_days(12)
integer :: month, day, i
logical :: leap

leap = (mod(year,4)==0 .and. mod(year,100)/=0) .or. mod(year,400)==0

month_days = (/31,28,31,30,31,30,31,31,30,31,30,31/)

if (leap) month_days(2)=29

day = doy
month = 1

do i=1,12
    if(day > month_days(i)) then
        day = day - month_days(i)
        month = month + 1
    else
        exit
    endif
enddo

write(cdate,"(i4.4,i2.2,i2.2)") year,month,day

end subroutine doy_to_date


!==============================
! Build full output path
!==============================
subroutine make_output_path(year, doy, glat, glon, ut_sec, instance, filepath)

integer, intent(in) :: year, doy, instance
real, intent(in) :: glat, glon, ut_sec
character(len=*), intent(out) :: filepath

  character(len=8)   :: cdate
  character(len=3)   :: cdoy
  character(len=20)  :: clat, clon
  character(len=256) :: output_root, base_dir
  character(len=3)   :: cinst
  character(len=256) :: full_dir


! Date conversion
call doy_to_date(year, doy, cdate)

write(cdoy,"(i3.3)") doy

! Latitude/longitude conversion
write(clat,"(f6.2)") glat
write(clon,"(f7.2)") glon

clat = adjustl(clat)
clon = adjustl(clon)

! instance to vsNNN
write(cinst,"(i3.3)") instance

  ! Build directory path (fixed output root)
  output_root = 'D:/1CH2/1study/DXK/zhanshi/spaceborne/code/glow model/glow model/output'
  base_dir = trim(output_root)//'/Lat'//trim(clat)//'_Lon'//trim(clon)//'_'//trim(cdate)
  full_dir = trim(base_dir)//'/'//trim(cdoy)

  ! Create directories
   call system('mkdir "'//trim(output_root)//'"')
   call system('mkdir "'//trim(base_dir)//'"')
   call system('mkdir "'//trim(full_dir)//'"')

  ! Final filename: Latxx.xx_Lonyy.yy_YYYYMMDD_vsNNN.txt
  filepath = trim(full_dir)//'/Lat'//trim(clat)//'_Lon'//trim(clon)//'_'//trim(cdate)//'_vs'//trim(cinst)//'.txt'

  ! Keep ut_sec referenced to avoid unused-argument warnings in strict builds.
  if (ut_sec < 0.0) filepath = trim(filepath)

end subroutine make_output_path
end module write_path

