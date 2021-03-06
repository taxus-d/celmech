module NonGradMin
    use Const
    use RandomFill
    use Deriv
    use IO_array
    use Inival
    implicit none
contains
    ! takes eternity to converge
    function randomfind(f,x0,N_max,report_fd) result(root)
        procedure (fRnR1) :: f
        real(mpc), dimension(:) :: x0
        integer :: N_max,i
        intent(in) :: x0, N_max
        real(mpc), dimension(size(x0)) :: root, x, xn, xp, shift
        integer, optional :: report_fd
        real(mpc) :: stepsize, mul

        x = x0; xp = huge(1.0_mpc)
        stepsize = 0.001; mul = 0.9
        i = 0
        do while (norm2(x - xp) > sqrt(eps) .and. i < N_max)
            call randomfill_arr(shift)
            shift = (shift - 0.5_mpc)*2
            shift = shift/norm2(shift)*stepsize
            xn = x + shift
            i = i + 1
            if ( f(xn) > f(x) ) cycle
            xp = x;  x = xn
            stepsize = stepsize * mul
            if (present(report_fd)) write(report_fd, *) x
            write(*,*) '>>' ,f(x)
        end do
        root = x
    end function randomfind

    ! converges slowly
    function mcsearch(f,x0,N_max,report_fd) result(root)
        procedure (fRnR1) :: f
        real(mpc), dimension(:) :: x0
        integer :: N_max,i, Npoints = 100, k
        intent(in) :: x0, N_max
        real(mpc), dimension(:,:), allocatable :: grid
        real(mpc), dimension(size(x0)) :: root, x, Rp, grad, Rc
        integer, optional :: report_fd
        real(mpc) :: mul, edge, Mc, dm
        
        allocate(grid(Npoints, size(x0)))
        edge = 1
        mul = 1
        Rc   = x0; Mc = 0; Rp = huge(1.0_mpc)
        k = 0
        do while (k < N_max) !norm2(Rc-Rp) > 0)
            if (present(report_fd)) write(report_fd, *) Rc
            call randomfill_arr(grid)
            grid = (grid - 0.5_mpc)*2*edge
!             write(*,*)
!             call plot_2d_array(f, grid, '')
            Rp = Rc
            Mc = 0; Rc = 0
            do i =1, Npoints
                dm = f(grid(i,:))
                Mc = Mc + dm
                Rc = Rc + dm*grid(i,:)
            end do
            Rc = Rc / Mc
            edge = newedge(edge,k)
            k = k + 1
        end do
        root = Rc
        deallocate(grid)
        contains
            pure function newedge(e,i) result(e1)
                real(mpc) :: e,e1
                integer :: i
                intent(in) :: e, i
                e1 = e*mul 
            end function newedge
    end function mcsearch

    function neldermead(f,x0,N_max,report_fd) result(root)
        
    end function neldermead

end module NonGradMin
