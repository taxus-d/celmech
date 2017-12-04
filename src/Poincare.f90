module Poincare
    use Celmech
    use Integrators
    implicit none

contains

    subroutine print_poincare_section(integ,t1,fd) result(x)
        class(Integrator) :: integ
        real(mpc) :: t1, S, S_prev
        integer   :: i, N, fd, fd_
        intent(in):: t1, fd
        optional  :: fd
        logical   :: printp
        fd_ = stdout
        if (present(fd)) fd_ = fd
        N=int((t1-integ%t0)/integ%h)
        do i = 1, N
            write(fd,*) integ%time(), integ%val() 
            call integ%step()
        end do

        do i = 0, N-1
            S_prev = S; S = p_S(x)
            if (S_prev*S < 0.0_mpc) then
                write (fd,*) integ%time(), integ%val()
                call integ%step(-S)
                weirdstep = .TRUE.
            else
                call integ%step()
                weirdstep = .FALSE.
            end if

        end do
    end function runge_ode

end module Integrators

