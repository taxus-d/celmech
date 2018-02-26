module GradMin
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
        real(mpc) :: mul, lambda, lambda0, searchsize, rms, av_factor
        logical :: updatep, advancedp_
        integer :: i, advmul, penalty, penalty_limit = 100

        updatep = .true.
        advancedp_ = .true.; advmul = 1
        if (present(advancedp)) advancedp_ = advancedp
        if (.not. advancedp_) advmul = 0
        penalty = 0
        xp = huge(1.0_mpc)
        i = 0
        x = x0
        lambda0 = 1_mpc; lambda = lambda0
        mul = 1
        av_factor = 0.99_mpc! as adviced
        searchsize = 0.005_mpc
        
        do while (norm2(x - xp) > eps .and. i < N_max)
            grad(1,:) = grad(2,:)
            grad(2,:) = ngrad(f,x)
            if (norm2(grad(2,:)) > 100.0_mpc) write(stderr,*) "Trouble: large grad =",&
                &norm2(grad(2,:))
            if (i > 0) then 
                dir = -grad(2,:) + advmul*beta(x,-grad(1,:), -grad(2,:),dir)*dir
                rms = rms * av_factor + (1.0_mpc - av_factor) * norm2(x-xp)**2
            else
                dir = -grad(2,:)
                rms = norm2(dir)**2
            endif 
            if (norm2(dir) > 1.0_mpc) dir = dir/norm2(dir)
            lambda = lambda0 / sqrt(rms+eps)
            x_bnd(1,:) = x - 1.0_mpc*dir*searchsize*lambda
            x_bnd(2,:) = x + 1.5_mpc*dir*searchsize*lambda
            
            xn = golden_linesearch(f, x_bnd(1,:), x_bnd(2,:), 100)

            if (f(xn) > f(x)-eps) penalty = penalty + 1
            xp = x
            x = xn
            if (present(report_fd)) write(report_fd, *) x
            write(stderr,*) '>>' ,f(x), norm2(grad(2,:)), rms

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
        function beta(x,dx_p,dx, s)
            real(mpc), dimension(:) :: dx
            real(mpc), dimension(size(dx)) :: x, dx_p, s
            real(mpc) :: beta
            real(mpc), parameter :: somewhat_small_const = sqrt(eps)/1000.0_mpc 
            intent(in) :: dx_p, dx, s
!             if (f(x) < somewhat_small_const) then
!                 beta = betaFR(dx_p, dx, s)
!             else
                beta = betaDY(dx_p, dx, s)
!             end if
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


    
end module GradMin
