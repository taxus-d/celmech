module NewtSolve
    use Inival
    use Deriv
    use FuncIfaces
    use LseSolvers
    use IO
    implicit none
    

contains
    
    !-----------------------------------------------------------------------
    ! newtitsolve (f) -> x
    ! f: \R^n \to \R^n
    ! x0 : initial value
    ! N_max  : max iterations number
    ! 	solves nonlinear system iteratively using Newton's method
    !-----------------------------------------------------------------------
    function newtitsolve(f, x_0, N_max) result(root)
        ! NOTE: по поводу интерфейсов : мне показалось, что будет лучше, если сигнатуры (FuncIfaces)
        !+ будут существовать отдельно от образцов функций (Samples).
        ! По крайней мере, больше разных проверок, что всё хорошо.
        ! (может я конечно все страшно усложнил и можно как-то просто и красиво)
        procedure (fRnRn) :: f
        real(mpc), dimension(:) :: x_0
        integer :: N_max
        intent(in) :: x_0, N_max
        real(mpc), dimension(size(x_0)) :: root

        integer :: n, i
        real(mpc), dimension(:), allocatable :: x, x_prev
        real(mpc), dimension(:, :), allocatable :: exmtx ! если вдруг в стек не влезет..
        n = size(x_0)
        
        allocate(x(n), x_prev(n), exmtx(n, n+1)) 
        x_prev = huge(1.0_mpc)
        x = x_0
        i = 1
        do while (i < N_max .and. norm2(x - x_prev) > 2*eps)
            x_prev = x
            exmtx(1:n, 1:n) = jacobimtx(f, x)
            exmtx(:, n+1) = -f(x) 
            x = gauss_solve(exmtx) ! решает относительно разности
            x = x + x_prev
            i = i+ 1
        end do

        root = x
        deallocate(x, x_prev, exmtx)

    end function newtitsolve
    


end module NewtSolve
