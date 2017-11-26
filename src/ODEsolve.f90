module ODEsolve
    use IO
    use Integrators
    implicit none
    contains
        subroutine print_solution(integ,t1,fd)
            class(Integrator) :: integ
            real(mpc) :: t1
            integer   :: i, N, fd, fd_
            intent(in):: t1, fd
            optional  :: fd
            fd_ = stdout
            if (present(fd)) fd_ = fd
            N=int((t1-integ%t0)/integ%h)
            do i = 1, N
                write(fd,*) integ%time(), integ%val() 
                call integ%step()
            end do
        end subroutine
end module ODEsolve
