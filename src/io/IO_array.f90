module IO_array
    use Const
    use FuncIfaces
    use IO
    implicit none
    
    interface prarr
        module procedure print_matrix 
        module procedure print_array
    end interface prarr
       
    interface plot_RnR1_section
        module procedure plot_RnR1_section_startdiredge
        module procedure plot_RnR1_section_startend
    end interface plot_RnR1_section
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
        logical :: titlep, sizep
        m = size(a, 1); n = size(a, 2)

        fdc = stdout
        if (present(fd)) fdc = fd
        
        sizep = .TRUE.; titlep = .FALSE.
        if (present(ttl)) then
            if (ttl == '') then
                sizep = .FALSE.; titlep = .FALSE.
            else
                titlep = .FALSE.
            end if
        end if

        if (titlep) write(fdc,'(t4,a1,a,a1)') '[',ttl,']' 
        if (sizep)  write (fdc, '(a1,i8,8x,i8)') '#', m, n
        if (titlep) write (fdc, *) repeat('-',13)
        
        do i = 1, m
            do j = 1, n
                write (fdc, '(e16.7)', advance="no") a(i, j)
            end do
            write (fdc, *) 
        end do
        write(*,*) ! blank line
    end subroutine print_matrix

    subroutine print_array(a, ttl, fd)
        real(mpc)       , dimension(:) :: a 
        integer         , optional     :: fd
        character(len=*), optional     :: ttl
        intent(in) :: a, fd, ttl
        integer :: i, m, fdc
        logical :: titlep, sizep
        m = size(a, 1)
        
        fdc = stdout
        if (present(fd)) fdc = fd
        sizep = .TRUE.; titlep = .FALSE.
        if (present(ttl)) then
            if (ttl == '') then
                sizep = .FALSE.; titlep = .FALSE.
            else
                titlep = .FALSE.
            end if
        end if

        if (titlep) write(fdc,'(t4,a1,a,a1)') '[',ttl,']' 
        if (sizep)  write (fdc, '(a1,i8)') '#', m
        if (titlep) write (fdc, *) repeat('-',13)
        do i = 1, m
            write (fdc, '(e14.7)') a(i)
        end do
        write(*,*) ! blank line
    end subroutine print_array
    
    subroutine plot_R2R1_func_on_array(f, a, ttl, fd)
        procedure(fRnR1) :: f
        real(mpc)       , dimension(:,:) :: a 
        integer         , optional       :: fd
        character(len=*), optional       :: ttl
        intent(in) ::     a, fd, ttl
        
        logical :: titlep, sizep
        integer :: i, m, fdc, n, j
        m = size(a, 1); n = size(a, 2)
        
        fdc = stdout
        if (present(fd)) fdc = fd
            
        fdc = stdout
        if (present(fd)) fdc = fd
        sizep = .TRUE.; titlep = .FALSE.
        if (present(ttl)) then
            if (ttl == '') then
                sizep = .FALSE.; titlep = .FALSE.
            else
                titlep = .FALSE.
            end if
        end if

        if (titlep) write(fdc,'(t4,a1,a,a1)') '[',ttl,']' 
        if (sizep)  write (fdc, '(a1,i8,8x,i8)') '#', m, n
        if (titlep) write (fdc, *) repeat('-',13)
        do i = 1, m
            do j = 1, n
                write (fdc, '(e16.7)', advance="no") a(i, j)
            end do
            write(fdc ,'(e16.7)') f(a(i,:))
        end do
        write(*,*) 
    end subroutine plot_R2R1_func_on_array
    subroutine plot_RnR1_section_startdiredge(f, x0, dir, edsize, N, ttl, fd)
        procedure(fRnR1) :: f
        real(mpc)       , dimension(:)          :: x0 
        real(mpc)       , dimension(size(x0))   :: dir 
        real(mpc)                               :: edsize
        integer         , optional       :: fd, N
        character(len=*), optional       :: ttl
        intent(in) ::  dir, x0, edsize, fd, ttl
        
        logical :: titlep
        integer :: i, m, fdc, j, N_
        real(mpc), dimension(size(x0)) :: x, edge
        
        fdc = stdout
        if (present(fd)) fdc = fd
            
        titlep = .FALSE.
        if (present(ttl)) then
            titlep = .TRUE.
        end if

        if (titlep) write(fdc,'(t4,a1,a,a1)') '[',ttl,']' 
        if (titlep) write (fdc, *) repeat('-',13)
        
        N_ = int(edsize/norm2(dir))
        if (present(N)) N_ = N
        edge = dir/norm2(dir)*edsize/N_
        x = x0
        do j = 1, N_
            write (fdc, *) (j-1), f(x)
            x = x + edge
        end do
    end subroutine plot_RnR1_section_startdiredge

    subroutine plot_RnR1_section_startend(f, x1, x2, N, ttl, fd)
        procedure(fRnR1) :: f
        real(mpc)       , dimension(:)          :: x1 
        real(mpc)       , dimension(size(x1))   :: x2 
        integer         , optional       :: fd, N
        character(len=*), optional       :: ttl
        intent(in) ::  x1, x2, fd, ttl
        
        logical :: titlep, sizep
        integer :: i, m, fdc, j, N_

        real(mpc) :: dir(size(x1)), edgesize

        edgesize = sqrt(sum((x2-x1)**2))
        dir = (x2-x1)/edgesize

        call plot_RnR1_section_startdiredge(f, x1, (x2-x1)/norm2(x2-x1), norm2(x2-x1), N, ttl, fd)
    end subroutine plot_RnR1_section_startend

    subroutine plot_coordinate_sections(f,x0,offset,Npoints,fd)
        procedure(fRnR1) :: f
        integer, optional :: fd
        integer :: Npoints
        real(mpc), dimension (:) :: x0
        real(mpc) :: offset
        intent(in) :: Npoints, fd, offset, x0
        
        real(mpc) :: x(size(x0)), delta
        integer :: i, k, fdc
        
        delta = offset / Npoints
        
        fdc = stdout
        if (present(fd)) fdc = fd
        
        do i = -Npoints/2, Npoints/2,1
            write(fdc,'(f11.8)', advance='no') delta*i
            do k = 1, size(x0)
                x = x0
                x(k) = x0(k) + delta*i
                write(fdc,'(f14.8)', advance='no') f(x)
            end do
            write(fdc,*)
        end do
    end subroutine plot_coordinate_sections

    subroutine plot_selected_2D(f, x0, indices, offset, resolution, fd)
        procedure(fRnR1) :: f
        integer, optional  :: fd, indices(2)
        real(mpc), dimension (:) :: x0
        real(mpc) :: offset, resolution
        intent(in) ::  fd, offset, x0, resolution, indices
        
        real(mpc) :: x(size(x0)), delta, d1, d2
        integer :: i,j, fdc, Npoints, ind(2)
        
        delta = resolution
        Npoints = int(offset / delta)
        
        fdc = stdout
        if (present(fd)) fdc = fd

        ind = (/1,2/)
        if (present(indices)) ind = indices
        
        x = x0
        do i = -Npoints/2, Npoints/2, 1
            d1 = delta*i
            x(ind(1)) = x0(ind(1)) + d1
            do j = -Npoints/2, Npoints/2, 1
                d2 = delta*j
                x(ind(2)) = x0(ind(2)) + d2
                write(fdc,*) x(ind(1)), x(ind(2)), f(x)
            end do
            write(fdc,*)
        end do
    end subroutine plot_selected_2D
end module IO_array

