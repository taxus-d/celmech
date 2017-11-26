module Poincare
    use Integrators
    implicit none

contains
    function run_runge_poincare(t0,t1,x0,fd) result(x)
        real(mpc), intent(in) :: t1, t0, x0(:)
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

end module Integrators

