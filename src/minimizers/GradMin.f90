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
        real(mpc) :: mul, lambda, lambda0, searchsize, searchsize0, rms, av_factor
        logical :: updatep, advancedp_
        integer :: i, advmul, penalty, penalty_limit = 5
        
        real(mpc) :: ls_ratio

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
        searchsize0 = 0.005_mpc; searchsize=searchsize0
        
        do while (norm2(x - xp) > eps .and. i < N_max)
            grad(1,:) = grad(2,:)
            grad(2,:) = ngrad(f,x)
            if (norm2(grad(2,:)) > 100.0_mpc) write(stderr,*) "Trouble: large grad =",&
                &norm2(grad(2,:))
            if (i > 0) then 
                dir = -grad(2,:) + advmul*beta(x,-grad(1,:), -grad(2,:),dir)*dir
            else
                dir = -grad(2,:)
                rms = norm2(dir)**2
            endif 
            if (norm2(dir) > 1.0_mpc) dir = dir/norm2(dir)
            lambda = lambda0
            x_bnd(1,:) = x - 1.0_mpc*dir*searchsize*lambda
            x_bnd(2,:) = x + 1.5_mpc*dir*searchsize*lambda
            
            xn = golden_linesearch(f, x_bnd(1,:), x_bnd(2,:), 100)
            ls_ratio = norm2(xn - x_bnd(1,:))/norm2(x_bnd(2,:) - x_bnd(1,:))
            if (abs(f(xn)) > abs(f(x))-eps) penalty = penalty + 1
            xp = x
            x = xn
!             searchsize = searchsize0
            if (present(report_fd)) write(report_fd, *) x
            write(stderr,*) '>>' ,f(x), norm2(grad(2,:)), ls_ratio, searchsize, penalty

            i = i + 1
            if (penalty >= penalty_limit) exit
        end do
        root = x
        if (present(retstat)) then
            if (penalty >= penalty_limit ) then
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

    ! x^TAx - b^Tx + c
    ! Axmult -- function (x) -> A*x
    function quad_conjgraddesc(Axmult, b, x0, N_max) result (root) 
        procedure (fRnRn) :: Axmult
        real(mpc), dimension(:), intent(in) :: b
        real(mpc), dimension(:), intent(in) :: x0
        integer, intent(in), optional :: N_max
        real(mpc), dimension(size(x0)) :: root

        integer :: N_max_, i
        real(mpc) , dimension(size(x0)) :: x, dir, mdir, residue, oldresidue
        real(mpc) :: alpha, beta

        ! N (dimension of space) steps should be sufficient
        N_max_ = size(x0)
        if (present(N_max)) N_max_ = N_max
        
        i = 0
        x = x0
        residue = b - Axmult(x0)
        dir = residue
        
        do while (i < N_max_ .and. norm2(residue) > eps)
            mdir = Axmult(dir) ! cached
!             call prarr(joinarr_sbs(dir, mdir))
            ! computed line search
            alpha      = -dot_product(residue, residue) / dot_product(dir, mdir)
            x          =  x + alpha * dir
            oldresidue =  residue
            residue    =  residue - alpha * mdir
            beta       =  dot_product(residue, residue) / dot_product(oldresidue, oldresidue)
            dir        =  residue + beta * dir
            i = i + 1
        end do
        root = x
    end function quad_conjgraddesc 
    
    function newtonhessianfree(f,x0,N_max,report_fd, retstat) result(root)
        procedure (fRnR1) :: f
        real(mpc), dimension(:), intent(in) :: x0
        integer, intent(in) :: N_max
        integer, optional :: report_fd
        integer, optional, intent(inout) :: retstat
        real(mpc), dimension(size(x0)) :: root, nzero
     
        real(mpc) :: point_value, point_grad(size(x0)), point_hess(size(x0), size(x0))
     
        real(mpc), dimension(size(x0)) :: x, xp, xn
        integer i

        xp = huge(1.0_mpc)
        i = 0
        x = x0
        nzero = 0.01!0
     
        do while (norm2(x - xp) > eps .and. i < 1)
            point_value = f(x)
            point_grad  = ngrad(f,x)
            point_hess  = hess_mtx(f,x)
            call prarr(point_hess)
!             xn = x+quad_conjgraddesc(hess_mmult, point_grad, nzero)
            xn = x + conjgraddesc(local_quad_approx, nzero, 100)
            xp = x
            x = xn
            if (present(report_fd)) write(report_fd, *) x
            write(stderr,*) '>>' ,f(x)
            i = i + 1
        end do
        root = x
        contains
            function ngrad_f(x) result (grad)
                real(mpc), intent(in), dimension (:) :: x
                real(mpc), dimension(size(x0)) :: grad
                grad =  ngrad(f,x)
            end function

            function hess_mmult(dx) result(Hdx)
                real(mpc), dimension(:), intent(in) :: dx
                real(mpc), dimension(size(dx)) :: Hdx
                Hdx = dir_deriv(ngrad_f, x, dx)
!                 write(*,*) ":dx:", dx
!                 write(*,*) ":∇f(x):", point_grad
!                 write(*,*) ":∇f(x+dx):", ngrad(f,x+dx*opt_shiftsize(x))
!                 write(*,*) ":res:", Hdx
            end function hess_mmult

            function local_quad_approx(dx) result(q)
                real(mpc), dimension(:), intent(in) :: dx
                real(mpc) :: q
!                 q = point_value + dot_product(point_grad,dx) +
                q = dot_product(dx, matmul(point_hess, dx))/2.0_mpc 
            end function local_quad_approx
    end function newtonhessianfree
end module GradMin
