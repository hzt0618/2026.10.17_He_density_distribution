program glowbasic_data   !数据整合 ! 错误，不使用
    implicit none
    
    ! 定义常量以便修改
    real, parameter :: ALT_START = 850.
    real, parameter :: ALT_END = 950.
    real, parameter :: ALT_LOW = 400.
    real, parameter :: ALT_HIGH = 500.
    integer :: zn = 10, nn = (ALT_END-ALT_START)/10. + 1
    integer :: jmax, nbins
    real :: sza
    integer :: idate
    real :: ut
    real :: glat, glong
    real :: f107a, f107, f107p
    integer :: ap
    real :: ef, ec
    integer :: iii
    integer :: i, j, k
    real, allocatable :: z(:), uflx(:,:), dflx(:,:), ener(:), del(:)
    real, allocatable :: uflx_new(:,:), dflx_new(:,:)
    real :: z_val
    integer :: index_low, index_high, index_start, index_end
    character(len=256) :: input_file, output_file, file_num

    ! 打开固定数据文件
    open(10, file='C:/Users/31217/Desktop/氦激光雷达/地磁场光电子模型/u500to900.txt')
    open(12, file='C:/Users/31217/Desktop/氦激光雷达/地磁场光电子模型/d500to900.txt')

    ! 读取替换数据
    allocate(uflx_new(190, nn))
    allocate(dflx_new(190, nn))
    
    ! 读取 u500to900.txt 数据
    do j = 1, nn
        read(10, "(1x,f6.1,190e12.5)") z_val, (uflx_new(i,j), i=1,190)
    end do

    ! 读取 d500to900.txt 数据
    do j = 1, nn
        read(12, "(1x,f6.1,190e12.5)") z_val, (dflx_new(i,j), i=1,190)
    end do

    ! 关闭固定数据文件
    close(10)
    close(12)

    ! 处理50个文件
    do k = 1, 50
        write(file_num, '(I0)') k
        
        ! 设置输入输出文件名
        input_file = 'HaiNan2024/6Dawn/output_vs' // trim(file_num) // '.txt'
        output_file = 'HaiNan2024/6Module/output_vs' // trim(file_num) // '.txt'
        
        ! 打开输入输出文件
        open(11, file=trim(input_file))
        open(13, file=trim(output_file))

        ! 读取第一行
        read(11, "(1x,2i7,f9.5,i10,f10.1,8f8.2,i6)") jmax, nbins, sza, idate, ut, glat, glong, f107a, f107, f107p, ap, ef, ec, iii
        write(13, "(1x,2i7,f9.5,i10,f10.1,8f8.2,i6)") jmax, nbins, sza, idate, ut, glat, glong, f107a, f107, f107p, ap, ef, ec, iii

        ! 动态分配内存
        allocate(z(jmax))
        allocate(uflx(nbins, jmax))
        allocate(dflx(nbins, jmax))
        allocate(ener(nbins))
        allocate(del(nbins))

        ! 读取 uflx 数据
        do j = 1, jmax
            read(11, "(1x,f6.1,190e12.5)") z(j), (uflx(i,j), i=1,nbins)
        end do

        ! 读取 dflx 数据
        do j = 1, jmax
            read(11, "(1x,f6.1,190e12.5)") z(j), (dflx(i,j), i=1,nbins)
        end do

        ! 读取 ener 数据
        read(11, "(1x,190f10.3)") (ener(j), j=1,nbins)

        ! 读取 del 数据
        read(11, "(1x,190f10.3)") (del(j), j=1,nbins)

        ! 找到 z(j) 为 ALT_LOW, ALT_HIGH, ALT_START, ALT_END 的索引
        index_low = 0
        index_high = 0
        index_start = 0
        index_end = 0
        do j = 1, jmax
            if (abs(z(j) - ALT_LOW) < 0.1) index_low = j
            if (abs(z(j) - ALT_HIGH) < 0.1) index_high = j
            if (abs(z(j) - ALT_START) < 0.1) index_start = j
            if (abs(z(j) - ALT_END) < 0.1) index_end = j
        end do

        ! 检查是否找到所有需要的索引
        if (index_start == 0 .or. index_end == 0) then
            print *, "Error: Could not find altitude indices in file ", trim(input_file)
            cycle
        end if

        ! 写入 ALT_START 之前的 uflx 数据
        do j = 1, index_start - 1
            write(13, "(1x,f6.1,190e12.5)") z(j), (uflx(i,j), i=1,nbins)
        end do

        ! 写入替换后的 ALT_START 到 ALT_END 的 uflx 数据
        z_val = real(ALT_START)
        do j = 1, nn
            write(13, "(1x,f6.1,190e12.5)") z_val, (uflx_new(i,j), i=1,nbins)
            z_val = z_val + zn
        end do

        ! 写入 ALT_END 之后的 uflx 数据
        do j = index_end + 1, jmax
            write(13, "(1x,f6.1,190e12.5)") z(j), (uflx(i,j), i=1,nbins)
        end do

        ! 写入 ALT_START 之前的 dflx 数据
        do j = 1, index_start - 1
            write(13, "(1x,f6.1,190e12.5)") z(j), (dflx(i,j), i=1,nbins)
        end do

        ! 写入替换后的 ALT_START 到 ALT_END 的 dflx 数据
        z_val = real(ALT_START)
        do j = 1, nn
            write(13, "(1x,f6.1,190e12.5)") z_val, (dflx_new(i,j), i=1,nbins)
            z_val = z_val + zn
        end do

        ! 写入 ALT_END 之后的 dflx 数据
        do j = index_end + 1, jmax
            write(13, "(1x,f6.1,190e12.5)") z(j), (dflx(i,j), i=1,nbins)
        end do

        ! 写入 ener 数据
        write(13, "(1x,190f10.3)") (ener(j), j=1,nbins)

        ! 写入 del 数据
        write(13, "(1x,190f10.3)") (del(j), j=1,nbins)

        ! 关闭文件
        close(11)
        close(13)

        ! 释放内存
        deallocate(z, uflx, dflx, ener, del)
        
    end do

    ! 释放替换数据内存
    deallocate(uflx_new, dflx_new)

end program glowbasic_data