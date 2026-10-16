module glow_config_module
  implicit none
  type :: glow_config_type
      ! 时间和日期
      integer :: year
      integer :: doy
      real(8) :: ltstart   ! 起始本地时（小时）
      real(8) :: ltend     ! 结束本地时（小时）
      integer :: insti     ! 时间步数

      ! 观测点（已经计算好的共轭点）
      real(8), allocatable :: lat_list(:)
      real(8), allocatable :: lon_list(:)

      ! 极光电子通量参数
      real(8) :: ef
      real(8) :: ec

      ! 太阳/地磁参数
      real(8) :: f107
      real(8) :: f107a
      real(8) :: f107p
      real(8) :: ap
  end type glow_config_type
end module glow_config_module