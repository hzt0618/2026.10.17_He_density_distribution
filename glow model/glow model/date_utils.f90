module date_utils

contains

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

end subroutine

end module date_utils