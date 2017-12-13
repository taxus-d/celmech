module Poincare
    use Celmech
    use NewtSolve
    use Integrators
    implicit none

contains

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
        
        do while (integ%t < t1)
            if (S_prev * S < 0) then
                not_intersectp = .FALSE.
                exit
            end if
            call integ%step()
            S_prev = S
            S      = p_S(integ%val())
            i      = i + 1
        end do

        if (.not. not_intersectp) then
            weirdstep = .TRUE.
            call integ%step(-S)
            weirdstep = .FALSE.
        end if
        
        if (present(retstat)) then 
            if (integ%t < t1) then 
                retstat = EXIT_SUCCESS
            else
                !> call warn("welcome to the end of the time")
                retstat = EXIT_FAILURE
            end if
        end if
    end subroutine shift_to_intersect
    
    
    subroutine print_poincare_section(integ,t1,fd)
        class(Integrator) :: integ
        real(mpc) :: t1, S, S_prev, x(integ%dimen)
        integer   :: i, N, fd, fd_, rets
        intent(in):: t1, fd
        optional  :: fd
        fd_ = stdout
        if (present(fd)) [fd_ = fd
       
        call shift_to_intersect(integ, t1, rets)
        do while (rets == EXIT_SUCCESS)
            x = integ%val()
!             if (x (7) > 0) write(fd,*) integ%time(), integ%val()
            write(fd,*) integ%time(), integ%val()
            call integ%step
            call shift_to_intersect(integ, t1, rets)
        end do

    end subroutine print_poincare_section


    subroutine imporove_inipos(integ,t0,x0,t1,fd)
        class(Integrator) :: integ
        real(mpc)  :: t1, t0,x0(:), x(size(x0)-2)
        integer    :: i, N, fd, rets
        intent(in) :: t1, t0, fd
        optional   :: fd
!         write(*,*) newtitsolve(intersection_diff, x0(1:size(x0)-1), 100, report_fd=stdout)
!         call prarr(intersection_diff(x0))

!         write(*,*) broydenitsolve(intersection_diff, x0(1:size(x0)-2), 1000,fd)
!          x = x0(1:size(x0)-1)
!         do i = -100, 100
!             x = x0(1:size(x0)-1) + i*sqrt(sqrt(eps))
!             write(*,*) i, intersection_diff(x)
!         end do
!         write(*,*) intersection_diff(x - 11*sqrt(sqrt(eps)))
    
    contains 
        function intersection_diff(xp) result(x)
            real(mpc), intent(in), dimension(:) :: xp
            real(mpc), dimension(size(xp)) :: x, x1, x2
            real(mpc) :: temp(size(xp)+2)
            integer :: i, cd, retstat
            cd = size(x0)
            call integ%set_inicond((/x0(1),xp,t0/), t0)
            call integ%crewind()
            call shift_to_intersect(integ, t1, retstat)
            temp = integ%val(); x1 = temp(2:cd-1)
            call integ%step()
            call shift_to_intersect(integ, t1, retstat)
            temp = integ%val(); x2 = temp(2:cd-1)
            x = x2 - x1
        end function 
        function intersection_diff_norm(xp) result(x)
            real(mpc), intent(in), dimension(:) :: xp
            real(mpc), dimension(size(xp)) :: x1, x2
            real(mpc) :: x
            integer :: i, cd, retstat
            cd = size(x0)
            call integ%set_inicond(xp, t0)
            call integ%crewind()
            call shift_to_intersect(integ, t1, retstat)
            x1 = integ%val()
            call integ%step()
            call shift_to_intersect(integ, t1, retstat)
            x2 = integ%val()
            x = norm2(x2 - x1)
        end function 
    end subroutine imporove_inipos

end module Poincare
