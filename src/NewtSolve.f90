module NewtSolve
    use Deriv
    use FuncIfaces
    use LseSolvers
    use IO
    implicit none
    
    type :: Ncube
        integer :: dimen
        integer, allocatable :: bdir(:)
        real(mpc),allocatable :: vertex(:)
        real(mpc) :: edge
    end type Ncube

contains
    
    !-----------------------------------------------------------------------
    ! newtitsolve (f) -> x
    ! f: \R^n \to \R^n
    ! x0 : initial value
    ! N_max  : max iterations number
    ! 	solves nonlinear system iteratively using Newton's method
    !-----------------------------------------------------------------------
    function newtitsolve(f, x0, N_max, report_fd) result(root)
        ! NOTE: по поводу интерфейсов :
        !+ мне показалось, что будет лучше, если сигнатуры (FuncIfaces)
        !+ будут существовать отдельно от образцов функций (Samples).
        ! По крайней мере, больше разных проверок, что всё хорошо.
        procedure (fRnRn) :: f
        real(mpc), dimension(:) :: x0
        integer :: N_max
        intent(in) :: x0, N_max
        real(mpc), dimension(size(x0)) :: root
        integer, optional :: report_fd

        integer :: n, i

        real(mpc), dimension(:), allocatable :: x, x_prev
        real(mpc), dimension(:, :), allocatable :: exmtx ! если вдруг в стек не влезет..
        n = size(x0)
        
        allocate(x(n), x_prev(n), exmtx(n, n+1)) 
        x_prev = huge(1.0_mpc)
        x = x0
        i = 1
        do while (i < N_max .and. norm2(x - x_prev) > sqrt(eps))
            if (present (report_fd)) write(report_fd,*) x, norm2(x-x_prev)
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

    function halfdivsolve(f,x0,ini_offset,N_max,report_fd) result(root)
        procedure (fRnRn) :: f
        real(mpc), dimension(:) :: x0
        real(mpc) :: ini_offset
        integer :: N_max
        integer, optional :: report_fd
        intent(in) :: x0, N_max, ini_offset
        
        real(mpc), dimension(size(x0)) :: root, newvert, rootdir
        type(Ncube) :: solvcube
        integer :: direction(size(x0)), D,i, j
        
        D = size(x0)
        forall (i = 1:D)
            direction(i) = 1
        end forall

        allocate(solvcube%bdir(D))
        allocate(solvcube%vertex(D))
        solvcube = Ncube (&
            &dimen  = D,&
            &bdir   = direction,&
            &vertex = x0 - ini_offset,&
            &edge   = 2*ini_offset)

!         solvcube%dimen  = D
!         solvcube%bdir   = direction
!         solvcube%vertex = x0 - ini_offset
!         solvcube%edge   = 2*ini_offset
        
        rootdir = -f(solvcube%vertex)
        
        associate (x=>solvcube%vertex)
        
        i = 0
        do while (i < N_max .and. norm2(rootdir) > sqrt(eps))

        if(present(report_fd)) write(report_fd,*) x 

            solvcube%edge = solvcube%edge / 2.0_mpc
            x = x + solvcube%edge * solvcube%bdir
            rootdir = -f(x)
            do j = 1, D
                if (abs(rootdir(j)) < eps) then
                    solvcube%bdir(j) = 0
                else
                    solvcube%bdir(j) = int(rootdir(j)/abs(rootdir(j)))
                end if
            end do
            i = i + 1

        end do

        end associate
        root = solvcube%vertex
        
        deallocate(solvcube%bdir)
        deallocate(solvcube%vertex)
    end function halfdivsolve
    
    function bisectsolve(f, x0, ini_offset, N_max, report_fd) result(root)
        procedure (fRnR1) :: f
        real(mpc), dimension(:) :: x0
        real(mpc) :: ini_offset
        integer :: N_max, i 
        integer, optional :: report_fd
        intent(in) :: x0, N_max, ini_offset
        real(mpc), dimension(size(x0)):: left, right, mid, root
        left = x0 - ini_offset; right = x0 + ini_offset

        if (f(left)*f(right) > 0) stop 'unable to solve'


        i=0
        do while (i < N_max .and. norm2(left - right) > sqrt(eps))
            mid = 0.5*(left + right)
            if (f(mid)*f(left) > 0) then 
                left = mid
            else
                right = mid
            endif
            i = i+1
        end do
        root = 0.5*(left+right)
    end function bisectsolve

!     function bisectsolve(f, x0, ini_offset, N_max, report_fd) result()
!         procedure (fRnR1) :: f
!         real(mpc), dimension(:) :: x0
!         real(mpc) :: ini_offset
!         integer :: N_max, n, k
!         integer, optional :: report_fd
!         intent(in) :: x0, N_max, ini_offset
!         real(mpc), dimension(size(x0)):: left, right, mid, root
!
!         n = size(x0)
!         do k = 1, n
!
!         end do
!
!
!     end function bisectsolve
    function broydenitsolve(f, x0, N_max, report_fd) result(root)
        ! NOTE: по поводу интерфейсов :
        !+ мне показалось, что будет лучше, если сигнатуры (FuncIfaces)
        !+ будут существовать отдельно от образцов функций (Samples).
        ! По крайней мере, больше разных проверок, что всё хорошо.
        procedure (fRnRn) :: f
        real(mpc), dimension(:) :: x0
        integer :: N_max
        intent(in) :: x0, N_max
        real(mpc), dimension(size(x0)) :: root
        integer, optional :: report_fd

        integer :: n, i

        real(mpc), dimension(:), allocatable :: xc, xp, fc, fp
        real(mpc), dimension(:, :), allocatable :: exmtx ! если вдруг в стек не влезет..
        n = size(x0)
        allocate(xc(n), xp(n), exmtx(n, n+1), fc(n), fp(n)) 
        associate (J=>exmtx(1:n,1:n), B=>exmtx(:,n+1))
        xp = x0
!         if (present(report_fd)) write(report_fd,*) xp
        fp = f(x0)
        J  = jacobimtx(f, xp)
        call prarr(J)
        B  = -fp
        xc = gauss_solve(exmtx) + xp
        i  = 1
        do while (i < N_max .and. norm2(xc - xp) > eps)
            if (present (report_fd)) write(report_fd,*) xc, norm2(f(xc))
            fc = f(xc) 
            J = J + dyad_product(((fc-fp) - matmul(J,(xc-xp)))/norm2(xc-xp)**2,xc-xp) 
            B = -fc 
            xp = xc
            fp = fc
            xc = gauss_solve(exmtx) + xp! решает относительно разности
            i = i + 1
        end do
        end associate
        root = xc
        deallocate(xc, xp, fc,fp, exmtx)
    end function broydenitsolve
    
end module NewtSolve
