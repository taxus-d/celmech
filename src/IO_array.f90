module IO_array
    use Prec
    use IO
    implicit none
    
    interface prarr
        module procedure print_matrix 
        module procedure print_array
    end interface prarr
        
contains
    subroutine read_mtx(a, a_fd)
        real(mpc), intent(inout), dimension(:, :) :: a
        integer, optional :: a_fd
        integer :: fd
        if (.not. present(a_fd)) fd = stdin
        
        read (fd, *) a
        a = transpose(a)
    end subroutine read_mtx

    subroutine read_exmtx(a,n, a_fd)
        integer, intent(in) :: n
        real(mpc), intent(inout), dimension(n, n+1) :: a
        integer, optional :: a_fd
        integer :: i, fd
        
        if (.not. present(a_fd)) fd = stdin
        
        read (fd, *) a(1:n, 1:n)
        a(1:n, 1:n) = transpose(a(1:n, 1:n))
        do i = 1, n
            read (fd, *) a(i, n+1)
        end do
    end subroutine read_exmtx

    subroutine print_matrix(a,ttl,fd)
        intent(in) :: a,ttl,fd
        optional   :: ttl, fd
        real(mpc)  :: a (:,:)
        integer    :: fd
        character(len=*) :: ttl

        integer :: i, j, n, m, fdc
        m = size(a, 1); n = size(a, 2)

        fdc = stdout
        if (present(fd)) fdc = fd
        
        write(*,*) 
        if (present(ttl)) &
            write(fdc,'(t4,a1,a,a1)') '[',ttl,']' 
        write (fdc, '(a1,i8,8x,i8)') '#', m, n
        if (present(ttl)) &
            write (fdc, *) repeat('-',31)
        do i = 1, m
            do j = 1, n
                write (fdc, '(e16.7)', advance="no") a(i, j)
            end do
            write (fdc, *) 
        end do
    end subroutine print_matrix

    subroutine print_array(a, fd, ttl)
        real(mpc)       , dimension(:) :: a 
        integer         , optional     :: fd
        character(len=*), optional     :: ttl
        intent(in) :: a, fd, ttl
        integer :: i, m, fdc
        m = size(a, 1)
        
        fdc = stdout
        if (present(fd)) fdc = fd

        write(*,*) 
        if (present(ttl)) &
            write(fdc,'(t4,a1,a,a1)') '[',ttl,']' 
        write (fdc, '(a1,i8)') '#', m
        if (present(ttl)) &
            write (fdc, *) repeat('-',13)
        do i = 1, m
            write (fdc, '(e14.7)') a(i)
        end do
    end subroutine print_array
    
end module IO_array

