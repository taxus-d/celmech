module Poincare
    use Celmech
    use Integrators
    implicit none

contains

    subroutine print_poincare_section(integ,t1,fd)
        class(Integrator) :: integ
        real(mpc) :: t1, S, S_prev
        integer   :: i, N, fd, fd_
        intent(in):: t1, fd
        optional  :: fd
        logical   :: printp
        fd_ = stdout
        if (present(fd)) fd_ = fd
        N=int((t1-integ%t0)/integ%h)
        
        S_prev = 0; S = p_S(integ%val()) 
        do i = 0, N-1
            if (S_prev*S < 0.0_mpc) then
                weirdstep = .TRUE.
                call integ%step(-S)
                write (fd,*) integ%time(), integ%val()
            else
                weirdstep = .FALSE.
                call integ%step()
            end if
            S_prev = S; S = p_S(integ%val())
        end do
    end subroutine print_poincare_section

end module Poincare

