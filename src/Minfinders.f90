module Minfinders
    use Const
    use RandomFill
    use Deriv
    use IO_array
    use Inival
    implicit none
contains
    ! breaks
    function fricgraddesc(f,x0,N_max,report_fd) result(root)
        procedure (fRnR1) :: f
        real(mpc), dimension(:) :: x0
        integer :: N_max,i, ordgrad(2)
        intent(in) :: x0, N_max
        real(mpc), dimension(size(x0)) :: root, x, xp, grad
        integer, optional :: report_fd
        real(mpc) :: mul, lambda, lambda0
        real(mpc), dimension(size(x0)) :: shift, vel, accel 
        real(mpc) :: frict,delta
       
        vel = 0; shift = 0
        frict = 0.1_mpc
        ordgrad = 0
        xp = huge(1.0_mpc)
        i = 1
        x = x0
        lambda0 = 100*sqrt(sqrt(eps))/2; lambda = lambda0
        delta = 0.1_mpc
        mul = 2
        do while (norm2(x - xp) > eps .and. i < N_max)
            if (present(report_fd)) write(report_fd, *) x
            write(*,*) '>>' ,f(x)
            grad = ngrad(f,x)
            if (norm2(grad) > 100) write(stderr,*) "Trouble: large grad =", norm2(grad)
            xp = x
!             lambda = lambda0
            accel = -grad - frict*vel
            vel   = vel   + delta*accel
            shift = delta*vel
            x = x + lambda*shift
            i = i + 1
        end do
        root = x
    end function fricgraddesc
    
    function golden_linesearch(f,a0,b0,N_max) result(root)
        procedure (fRnR1) :: f
        real(mpc), dimension(:)       :: a0
        real(mpc), dimension(size(a0)) :: b0
        integer :: N_max
        intent(in) :: a0, b0, N_max

        integer :: i
        real(mpc) :: x(2, size(a0)), delta, vals(2)
        real(mpc), dimension(size(a0)) :: root, a, b
        a = a0; b = b0; i  = 0
        x(1,:) = b - (b-a)/phi
        x(2,:) = a + (b-a)/phi
        do while (i < N_max .and. norm2(a-b) > sqrt(eps)/10.0)
            vals = (/ f(x(1,:)), f(x(2,:)) /)
            if (vals(1) > vals(2)) then
                a    = x(1,:)
                x(1,:) = x(2,:)
                x(2,:) = b - (x(1,:) - a   )
            else
                b    = x(2,:)
                x(2,:) = x(1,:)
                x(1,:) = a + (b    - x(2,:))
            end if
            i = i + 1
        end do
        root = (a + b)/2
    end function golden_linesearch

    function bin_linesearch(f,a0,b0,N_max) result(root)
        procedure (fRnR1) :: f
        real(mpc), dimension(:)       :: a0
        real(mpc), dimension(size(a0)) :: b0
        integer :: N_max
        intent(in) :: a0, b0, N_max

        integer :: i
        real(mpc) :: x(size(a0)), delta, vals(3)
        real(mpc), dimension(size(a0)) :: root, a, b
        a = a0; b = b0; i  = 0
        x = (a + b)/2
        do while (i < N_max .and. norm2(a-b) > sqrt(eps)/10.0)
            x = (a+b)/2
            vals = (/f(a), f(x), f(b)/)
            if (vals(1) < vals(3)) then
                b = x
            else
                a = x
            end if
        end do
        root = (a + b)/2

    end function bin_linesearch
    function conjgraddesc(f,x0,N_max,report_fd, retstat, advancedp) result(root)
        procedure (fRnR1) :: f
        real(mpc), dimension(:), intent(in) :: x0
        integer, intent(in) :: N_max
        integer, optional :: report_fd
        integer, optional, intent(inout) :: retstat
        logical, optional :: advancedp
        real(mpc), dimension(size(x0)) :: root
        
        real(mpc), dimension(size(x0)) :: x, xp, xn
        real(mpc) :: dir(size(x0)), grad(2, size(x0)), x_bnd(2,size(x0))
        real(mpc) :: mul, lambda, lambda0, searchsize
        logical :: updatep, advancedp_
        integer :: i, advmul, penalty, penalty_limit = 5

        updatep = .true.
        advancedp_ = .true.; advmul = 1
        if (present(advancedp)) advancedp_ = advancedp
        if (.not. advancedp_) advmul = 0
        penalty = 0
        xp = huge(1.0_mpc)
        i = 0
        x = x0
        lambda0 = 1; lambda = lambda0
        mul = 1; 
        searchsize = 0.005_mpc
        
        do while (norm2(x - xp) > eps .and. i < N_max)
            grad(1,:) = grad(2,:)
            grad(2,:) = ngrad(f,x)
            if (norm2(grad(2,:)) > 100.0_mpc) write(stderr,*) "Trouble: large grad =",&
                &norm2(grad(2,:))
            if (i > 0) then 
                dir = -grad(2,:) + advmul*beta(-grad(1,:), -grad(2,:),dir)*dir
            else
                dir = -grad(2,:)
            endif 
            if (norm2(dir) > 1.0_mpc) dir = dir/norm2(dir)
            x_bnd(1,:) = x - 1.0_mpc*dir*searchsize*lambda
            x_bnd(2,:) = x + 1.5_mpc*dir*searchsize*lambda
!             write(*,*) dot_product(dir, (x - x0_ideal(1:size(x0_ideal)-1)))&
!             &/(norm2(dir)*norm2(x-x0_ideal(1:size(x0_ideal)-1)))
            xn = golden_linesearch(f, x_bnd(1,:), x_bnd(2,:), 100)

            if (f(xn) > f(x)-eps) penalty = penalty + 1
            xp = x
            x = xn
            if (present(report_fd)) write(report_fd, *) x
            write(stderr,*) '>>' ,f(x), norm2(grad(2,:))

            i = i + 1
            if (penalty > penalty_limit) exit
        end do
        root = x
        if (present(retstat)) then
            if (penalty > penalty_limit ) then
                retstat = EXIT_FAILURE 
            else
                retstat = EXIT_SUCCESS
            end if
        end if 
    contains
        pure function beta(dx_p,dx, s)
            real(mpc), dimension(:) :: dx
            real(mpc), dimension(size(dx)) :: dx_p, s
            real(mpc) :: beta
            intent(in) :: dx_p, dx, s
            if (f()) then
                beta = betaFR(dx_p, dx, s)
                write(*,*) 'little'
            else
                beta = betaDY(dx_p, dx, s)
            end if
        end function beta
    end function conjgraddesc
    pure function betaPR(dx_p, dx, s) result(beta)
        real(mpc), dimension(:) :: dx
        real(mpc), dimension(size(dx)) :: dx_p, s
        real(mpc) :: beta
        intent(in) :: dx_p, dx, s
        beta = dot_product(dx, dx-dx_p)/dot_product(dx_p,dx_p)
    end function
    pure function betaFR(dx_p, dx, s) result(beta)
        real(mpc), dimension(:) :: dx
        real(mpc), dimension(size(dx)) :: dx_p, s
        real(mpc) :: beta
        intent(in) :: dx_p, dx, s
        beta = dot_product(dx, dx)/dot_product(dx_p,dx_p)
    end function
    pure function betaHS(dx_p, dx, s) result(beta)
        real(mpc), dimension(:) :: dx
        real(mpc), dimension(size(dx)) :: dx_p, s
        real(mpc) :: beta
        intent(in) :: dx_p, dx, s
        beta = -dot_product(dx, dx - dx_p)/dot_product(s,dx - dx_p)
    end function
    pure function betaDY(dx_p, dx, s) result(beta)
        real(mpc), dimension(:) :: dx
        real(mpc), dimension(size(dx)) :: dx_p, s
        real(mpc) :: beta
        intent(in) :: dx_p, dx, s
        beta = -dot_product(dx, dx)/dot_product(s,dx - dx_p)
    end function


    
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
end module Minfinders
