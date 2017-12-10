module Poincare
    use Celmech
    use NewtSolve
    use Integrators
    implicit none

contains

    subroutine print_poincare_section(integ,t1,fd)
        class(Integrator) :: integ
        real(mpc) :: t1, S, S_prev
        integer   :: i, N, fd, fd_, rets
        intent(in):: t1, fd
        optional  :: fd
        fd_ = stdout
        if (present(fd)) fd_ = fd
        N=int((t1-integ%t0)/integ%h)
        
!         S_prev = 0; S = p_S(integ%val())
!         do i = 0, N-1
!             if (S_prev*S < 0.0_mpc) then
!                 weirdstep = .TRUE.
!                 call integ%step(-S)
!                 write (fd,*) integ%time(), integ%val()
!             else
!                 weirdstep = .FALSE.
!                 call integ%step()
!             end if
!             S_prev = S; S = p_S(integ%val())
!         end do
       
        rets = EXIT_SUCCESS
        do while (rets == EXIT_SUCCESS)
            call shift_to_intersect(integ, t1, rets)
            write(fd,*) integ%time(), integ%val()
        end do

    end subroutine print_poincare_section

    subroutine shift_to_intersect(integ, t1, retstat)
        class(Integrator), intent(inout) :: integ
        integer, intent(out), optional :: retstat
        real(mpc) :: t1, S, S_prev
        integer   :: i, N
        logical :: not_intersectp

        N=int((t1-integ%t0)/integ%h)
        S = p_S(integ%val()) 
        S_prev = S;
        i=0
        weirdstep = .FALSE.
        not_intersectp = .TRUE.
        
        do while (not_intersectp .and. i < N)
            if (S_prev * S < 0) exit
            call integ%step()
            S_prev = S
            S      = p_S(integ%val())
            i      = i + 1
        end do
        
        weirdstep = .TRUE.
        call integ%step(-S)
        weirdstep = .FALSE.
        
        if (present(retstat)) then 
            if (i >= N) then 
                !> call warn("welcome to the end of the time")
                retstat = EXIT_FAILURE
            else
                retstat = EXIT_SUCCESS
            end if
        end if
    end subroutine shift_to_intersect
    subroutine imporove_inipos(integ,x0,t1,fd)
        class(Integrator) :: integ
        real(mpc)  :: t1, x0(:), S, S_prev
        integer    :: i, N, fd, rets
        intent(in) :: t1, fd
        optional   :: fd

    end subroutine imporove_inipos

end module Poincare

