module IO_tdmtx
    use Prec
    use IO
    implicit none

contains
        
    subroutine read_tdmtx3(a, fd_a)
        implicit none
        real(mpc), intent(out), dimension(:, :) :: a 
        integer, optional :: fd_a
        integer :: n, i, fd
        n = size(a, 1)
        fd = stdin
        if (present(fd_a)) fd = fd_a
        
        read (fd, *) a(1, 2:3)
        do i = 2, n-1
            read (fd, *) a(i, :)
        end do
        read (fd, *) a(n, 1:2)
    end subroutine read_tdmtx3

    subroutine print_tdmtx3(a, fd_a)
        implicit none
        real(mpc), intent(in), dimension(:, :) :: a
        integer, intent(in), optional :: fd_a
        integer :: n, i, fd
        n = size(a, 1)
        fd = stdout
        if (present(fd_a)) fd = fd_a

        write (*, *) "#", n
        write (fd, *) a(1, 2:3)
        do i = 2, n - 1
            write (fd, *) a(i,:)
        end do
        write (fd, *) a(n, 1:2)
    end subroutine print_tdmtx3
    
    subroutine print_tdmtx5(a, fd_a)
        implicit none
        real(mpc), intent(in), dimension(:, :) :: a
        integer, intent(in), optional :: fd_a
        integer :: n, i, fd
        n = size(a, 1)
        fd = stdout
        if (present(fd_a)) fd = fd_a

        write (*, *) "#", n
        write (fd, *) a(1, 3:5)
        write (fd, *) a(2, 2:5)
        do i = 3, n - 2
            write (fd, *) a(i,:)
        end do
        write (fd, *) a(n-1, 1:4)
        write (fd, *) a(n, 1:3)
    end subroutine print_tdmtx5
end module

