module Integrators
    use Inival
    use IO
    use IO_array
    use Poly
    use NewtSolve
    implicit none
   
    logical :: printp

contains

    function euler_step(t, x, h) result(x1)
        real(mpc), intent(in) :: t, x(:), h
        real(mpc) :: x1(size(x))

        x1 = x + h * f(t,x)
    end function euler_step
    
    function runge_step(t, x, h) result(x1)
        real(mpc), intent(in) :: t, x(:), h
        real(mpc) :: x1(size(x)), k(4, size(x))

        k(1,:) = h*f(t, x) 
        k(2,:) = h*f(t + h*0.5_mpc, x + 0.5_mpc*k(1,:)) 
        k(3,:) = h*f(t + h*0.5_mpc, x + 0.5_mpc*k(2,:)) 
        k(4,:) = h*f(t + h, x + k(3,:)) 
        x1 = x + 1.0_mpc/6.0_mpc * & 
            (k(1,:) + 2.0_mpc*k(2,:) + 2.0_mpc*k(3,:) + k(4,:))
    end function runge_step
    function runge_ode(t1, fd) result(x)
        real(mpc), intent(in) :: t1
        real(mpc) :: x(size(x0)), t, step
        integer :: i, N, fd
        optional :: fd
        real(mpc) :: S, S_prev
        N = int((t1-t0)/h)
        t = t0
        x = x0
        S_prev = 0; S = p_S(x)
        do i = 0, N-1
            if (present(fd) .and. printp) write(fd,*) t, x
            x = runge_step(t,x, step)
            t = t + step
            S_prev = S; S = p_S(x)
            if (S_prev*S < 0.0_mpc) then
                step = -S; weirdstep = .TRUE.
                printp = .TRUE.
            else
                step = h; weirdstep = .FALSE.
                printp = .FALSE.
            end if

        end do
    end function runge_ode


    function adams_ex_ode(t1,fd) result(x)
        real(mpc), intent(in) :: t1
        real(mpc) :: x(size(x0)), t, Ac(0:ad_ord-1), &
            Xc(0:ad_ord-1,size(x0)), fc(0:ad_ord-1,size(x0)) ! c -- cached
        integer :: i, Nsteps, n, j, fd
        optional :: fd
        Nsteps = int((t1-t0)/h)
        t = t0
        x = x0
        n = ad_ord ! чуть короче
        do j=0,n-1
            Ac(j) = A(n, j) 
            Xc(j,:) = runge_ode(t0 + h*j)
        end do

        do j=0, n-1
            fc(j,:) = f(t0+h*(n-1-j), Xc(n-1-j, :))
        end do
        t = t0 + h*(n-1)
        do i = n-1, Nsteps
            if (present(fd)) write(fd,*) t, Xc(n-1,:)
            Xc(0:n-2,:) = Xc(1:n-1,:)
            Xc(n-1,:) = Xc(n-1,:) + h * matmul(Ac, fc) 
            t = t + h
            fc(1:n-1,:) = fc(0:n-2,:)
            fc(0,:)     = f(t, Xc(n-1,:))
        end do
        j = min(Nsteps, n-1)
        x= Xc(j,:)

        contains
            function A(n,j)
                integer, intent(in) :: n,j
                real(mpc) :: A, rts(0:n-2), p(0:n-1)
                integer :: i
                ! \int_0^1 ...
                !! корни
                forall(i=0:j-1) rts(i) = i
                forall(i=j+1:n-1) rts(i-1) = i !пропуск jго корня
                !! a0 = 1
                p = poly_by_roots(1.0_mpc, -rts)
                !! интегрируем, уже случайно (почти) получились нужные пределы
                forall (i=0:n-1) p(i) = p(i)/(n-i)
                A = (1-2*mod(j,2)) * sum (p)
                do i=1,j     ; A = A/i; end do
                do i=1,n-1-j ; A = A/i; end do
            end function A
    end function adams_ex_ode

    function adams_in_ode(t1, fd) result(x)
        real(mpc), intent(in) :: t1
        real(mpc) :: x(size(x0)), t, Bc(-1:ad_ord-2), &
            Xc(0:ad_ord-1,size(x0)), fc(-1:ad_ord-2,size(x0)) ! c -- cached
        integer :: i, Nsteps, n, j, fd
        optional :: fd
        Nsteps = int((t1-t0)/h)
        t = t0
        x = x0
        n = ad_ord ! чуть короче
        do j=-1,n-2
            Bc(j) = B(n, j)
        end do
        do j=0, n-1
            Xc(j,:) = runge_ode(t0 + h*j)
        end do
        do j=-1, n-2
            fc(j,:) = f(t0+h*(n-2-j), Xc(n-2-j, :))
        end do
        t = t0 + h*(n-2)
        do i = n-2, Nsteps ! а это место больше всего похоже на дженгу
            if (present(fd)) write(fd,*) t, Xc(n-2,:)

            Xc(n-1,:) = newtitsolve(impl, Xc(n-2, :),100)
            Xc(0:n-2, :) = Xc(1:n-1, :)
            t = t + h
            fc(-1,:)     = f(t, Xc(n-2,:))
            fc(0:n-2,:) = fc(-1:n-3,:)
        end do
        j = min(Nsteps, n-2)
        x= Xc(j,:)
        contains
            ! это не работает в старом фортране и некрасиво
            ! но иначе нужно переписывать Ньютона
            function impl(x) result(y)
                real(mpc) :: x(:) , y(size(x))
                intent(in) :: x
                y = x  - h*Bc(-1)*f(t+h,x) - Xc(n-2,:) - h*matmul(Bc(0:n-2),fc(0:n-2,:))
            end function impl
            function B(n,j)
                integer, intent(in) :: n,j
                real(mpc) :: B, rts(-1:n-3), p(0:n-1)
                integer :: i
                ! \int_0^1 ...
                !! корни
                forall(i=-1:j-1) rts(i) = i
                forall(i=j+1:n-2) rts(i-1) = i !пропуск jго корня
                !! a0 = 1
                p = poly_by_roots(1.0_mpc, -rts)
                !! интегрируем, уже случайно (почти) получились нужные пределы
                forall (i=0:n-1) p(i) = p(i)/(n-i)
                B = (1 - 2*mod(j+1,2)) * sum (p)
                do i=1,j+1     ; B = B/i; end do
                do i=1,n-2-j   ; B = B/i; end do
            end function B
    end function adams_in_ode

end module Integrators

