! vile hack to share macros
#include "macros.f90"

module Poincare
    use Utils
    use Celmech
    use Inival
    use NewtSolve
    use Minfinders
    use Integrators
    implicit none
    logical :: weirdstep
    interface
        function eq_fun(t,X) result(f)
            import :: mpc
            real(mpc), intent(in) :: t, X(:)
            real(mpc) :: f(size(X))
        end function eq_fun
        pure function p_sec(X) result(S)
            import :: mpc
            real(mpc), dimension(:),intent(in) :: X
            real(mpc) :: S
        end function p_sec
        pure function p_sec_deriv(X) result(dS)
            import :: mpc
            real(mpc), dimension(:),intent(in) :: X
            real(mpc), dimension(size(X)) :: dS
        end function p_sec_deriv
    end interface
    type :: Psection
        procedure(p_sec),nopass, pointer :: S
        procedure(p_sec_deriv),nopass, pointer :: dS
    end type Psection
    type(Psection) :: current
    type(Psection), dimension(Nbodies) :: sections
contains
    pure function sect_x_1_body(X) result (S)
        real(mpc), dimension(:),intent(in) :: X
        real(mpc) :: S
        S = X(1) - 0.3
    end function
    
    
    pure function sect_x_1_body_deriv(X) result (dS)
        real(mpc), dimension(:) :: X
        real(mpc), dimension(size(X)) :: dS
        intent(in) :: X
        dS = 0
        dS(1) = 1.0_mpc
    end function

    
    pure function sect_x_2_body(X) result (S)
        real(mpc), dimension(:),intent(in) :: X
        real(mpc) :: S
        S = X(3) - 0.3
    end function
    
    
    pure function sect_x_2_body_deriv(X) result (dS)
        real(mpc), dimension(:) :: X
        real(mpc), dimension(size(X)) :: dS
        intent(in) :: X
        dS = 0
        dS(3) = 1.0_mpc
    end function
    
    
    pure function sect_x_3_body(X) result (S)
        real(mpc), dimension(:),intent(in) :: X
        real(mpc) :: S
        S = X(5) - 0.3
    end function
    
    
    pure function sect_x_3_body_deriv(X) result (dS)
        real(mpc), dimension(:) :: X
        real(mpc), dimension(size(X)) :: dS
        intent(in) :: X
        dS = 0
        dS(5) = 1.0_mpc
    end function

    ! a workaround to avoid explicit section pointer passing
    function poincare_section_eq(tt,X) result(f)
        real(mpc), intent(in) :: tt
        real(mpc), dimension(:), intent(in) :: X
        real(mpc), dimension(size(X)) :: f
        real(mpc) :: H
        
        f(1:tDim) = fast_motion_eq(tt, X(1:tDim))
        f(tDim+1) = 1
        H = 1
        if (weirdstep) H=dot_product(current%dS(X(1:tDim)),f(1:tDim))
        f = f/H
    end function
    

    subroutine shift_to_intersect(integ, t1, retstat, fd)
        class(Integrator), intent(inout) :: integ
        integer, intent(out), optional :: retstat
        integer, optional :: fd
        real(mpc) :: t1, S, S_prev
        integer   :: i, N
        logical :: not_intersectp
        
        N=int((t1-integ%t0)/integ%h)
        S = current%S(integ%val()) 
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
            if (present(fd)) write(*,*) integ%val()
            S_prev = S
            S      = current%S(integ%val())
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
        if (present(fd)) fd_ = fd
        if (.not. associated(current%S) .or. .not. associated(current%dS)) &
            &error stop "Specify `currect` section!"
        call shift_to_intersect(integ, t1, rets)
        do while (rets == EXIT_SUCCESS)
            x = integ%val()
            write(fd,*) integ%time(), integ%val(), norm2(integ%val())
            call integ%step
            call shift_to_intersect(integ, t1, rets)
            call integ%step
            call shift_to_intersect(integ, t1, rets)
        end do

    end subroutine print_poincare_section

    ! smart macros, where are you..(
#define std_sect_assign(i,j) \
    sections(i)%S => CONCAT(sect_x_,j)_body; \
    sections(i)%dS => CONCAT(sect_x_,j)_body_deriv


    subroutine imporove_inipos(integ,t0,x0,t1,fd)
        class(Integrator) :: integ
        real(mpc)  :: t1, t0,x0(:), x(size(x0)-1), delta, d1,d2,cell
        integer    :: i, N, fd, fd_, rets, k, j
        intent(in) :: t1, t0, fd
        optional   :: fd
        integer    :: mixp, skipn
        logical    :: advanced_descent_p
        mixp = EXIT_SUCCESS
        
        fd_ = stdout
        if(present (fd)) fd_ = fd

        std_sect_assign(1,2)
        std_sect_assign(2,3)
        std_sect_assign(3,1)
        
        ! extra check
        do i = 1, Nbodies
            if (.not. associated(sections(i)%S) &
                &.or. .not. associated(sections(i)%dS)) &
                & error stop "Define section list before improving!"
        end do
        
        associate(wsp=>x0(1:size(x0)-1), wsp_id=>x0_ideal(1:size(x0)-1))
        
!         x0(2:size(x0)-1) = broydenitsolve(intersection_diff, x, 100,fd_ )
        skipn = 0; advanced_descent_p = .true.
        do i=1,5
            write(*,*) 'conjugate gradient descent ::'
            wsp = conjgraddesc(intersection_diff_scalar,wsp,100,retstat=mixp, advancedp=advanced_descent_p)
            if (mixp == EXIT_FAILURE) then 
                skipn = skipn + 2
                advanced_descent_p = .false.
!                 write(*,*) 'inertion&friction descent simulation ::'
!                 wsp = fricgraddesc(intersection_diff_scalar, wsp, 200)
            end if 
        end do
!         call plot_RnR1_section(intersection_diff_scalar, &
!             &wsp_id, (wsp - wsp_id),&
!             &0.0001_mpc, N=30, fd=stdout)
!         k = 1; delta = 0
!         do i = -100, 100,1
!             delta = sqrt(eps)*i
!             write(*,'(f11.8)', advance='no') delta
!             do k = 1, size(x)
!                 x = x0(1:size(x))
!                 x(k) = x0(k) + delta
!                 write(*,'(f14.8)', advance='no') f(x)
!             end do
!             write(*,*)
!         end do
        
!         x = x0(1:size(x))
!         cell = sqrt(sqrt(eps))/10.0_mpc
!         do i = -100, 100,1
!             d1 = cell*i
!             x(1) = x0(1) + d1
!             do j = -100,100,1
!                 d2 = cell*j
!                 x(2) = x0(2) + d2
!                 write(*,*) x(1), x(2), f(x)
!             end do
!             write(*,*)
!         end do
        end associate
        
    contains 
        function intersection_diff(xp) result(x)
            real(mpc), intent(in), dimension(:) :: xp
            real(mpc), dimension(size(xp)) :: x
            real(mpc) :: temp(3, size(x0)), norms(2)
            integer :: i, cd, retstat, inds(1), b = 1
            cd = size(x0)
            x = 0
            call integ%set_inicond((/xp,t0/), t0)
            call integ%crewind()
                
            current%S => sections(2)%S
            current%dS => sections(2)%dS
            call shift_to_intersect(integ, t1, retstat)
            temp(1, :) = integ%val()
            
            current%S => sections(3)%S
            current%dS => sections(3)%dS
            call shift_to_intersect(integ, t1, retstat)
            temp(2, :) = integ%val()
            temp(2, 1:cd-1) = cyclic_shift_bodies(temp(2, 1:cd-1), 1)
            
            current%S => sections(1)%S
            current%dS => sections(1)%dS
            call shift_to_intersect(integ, t1, retstat)
            temp(3, :) = integ%val()
            temp(3, 1:cd-1) = cyclic_shift_bodies(temp(3, 1:cd-1), -1)
            
            if (retstat == EXIT_FAILURE) then 
                call loud_warn("welcome to the end of the time")
                stop '!'
            end if
            x = temp(3, b:cd-1) + temp(2, b:cd-1) - temp(1, b:cd-1)
        end function 

        function intersection_diff_scalar(xp) result(dev)
            real(mpc), intent(in), dimension(:) :: xp
            real(mpc) :: dev
            real(mpc) :: temp(6, size(x0))
            integer :: i, cd, retstat, b = 1, dir, sectn
            cd = size(x0)
            call integ%set_inicond((/xp,t0/), t0)
            call integ%crewind()
            
            do i = 0, skipn
                current%S => sections(1)%S
                current%dS => sections(1)%dS
                call shift_to_intersect(integ, t1, retstat)
                call integ%step()
                call shift_to_intersect(integ, t1, retstat)
                if (retstat == EXIT_FAILURE) then 
                    call loud_warn("welcome to the end of the time")
                end if
            end do
            do i = 1, 2*Nbodies
                sectn = mod(i-1, Nbodies) + 1
                current%S => sections(sectn)%S
                current%dS => sections(sectn)%dS
                
                call shift_to_intersect(integ, t1, retstat)
                temp(i, :) = integ%val()
               
                dir = mod(i, Nbodies) - 1
                temp(i, 1:cd-1) = cyclic_shift_bodies(temp(i, 1:cd-1), dir)
                
                if (retstat == EXIT_FAILURE) then 
                    call loud_warn("welcome to the end of the time")
                    stop '!'
                end if
            end do

            dev =  &
            &log(1.0_mpc &
                & + norm2(temp(6, 1:cd-1) - temp(3, 1:cd-1))**(2)&
                & + norm2(temp(5, 1:cd-1) - temp(2, 1:cd-1))**(2)&
                & + norm2(temp(4, 1:cd-1) - temp(1, 1:cd-1))**(2)&
                & + norm2(temp(2, 1:cd-1) - temp(1, 1:cd-1))**(2)&
                & + norm2(temp(3, 1:cd-1) - temp(2, 1:cd-1))**(2))
        end function 
        
        function test(x) result(f)
            real(mpc) , intent(in) :: x(:)
            real(mpc) :: f
            f = x(1)**2 + x(2)**2  
        end function test
        function f_2_sq(x) result(y)
        real(mpc),intent(in), dimension(:) :: x
        real(mpc), dimension(size(x)) :: y
        y = (/ (x(1) - 2.0_mpc)**2 + (x(2) - 3.0_mpc)**2 - 13.0_mpc ,&
              &(x(1) - 4.0_mpc)**2 + (x(2) - 1.0_mpc)**2 -17.0_mpc /)
    end function f_2_sq

        function test2(x) result(f)
            real(mpc) , intent(in) :: x(:)
            real(mpc) :: f
            f = 1.0_mpc/(x(1)**2 + x(2)**2)
        end function test2
    end subroutine imporove_inipos

end module Poincare
