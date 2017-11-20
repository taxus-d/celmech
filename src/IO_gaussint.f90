module IO_gaussint
    use Prec
    use LegendrePoly
    use IO
    
contains
    
    subroutine print_gauss_coef_roots(n, fd_out)
        integer, intent(in) :: n, fd_out
        optional :: fd_out
        integer :: fd, i
        real(mpc), allocatable, dimension (:) :: roots, coefs

        if (.not. present(fd_out)) then 
            fd = stdout
        else
            fd = fd_out
        end if
        
        allocate(roots(n), coefs(n))
        write(*,*) 'n = ', n
        roots = legendre_roots(n)
        coefs = int_coefs(roots)

        do i = 1, n
            write(fd, *)  coefs(i), roots(i)
        end do
        
        deallocate(roots, coefs)
        
    end subroutine print_gauss_coef_roots

end module IO_gaussint
