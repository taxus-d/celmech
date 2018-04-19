! vile hack to share macros
#include "macros.f90"

module Poincare
    use omp_lib
    use Utils
    use Celmech
    use Inival
    use NewtSolve
    use GradMin
    use Integrators
    implicit none
    logical   :: weirdstep
    real(mpc) :: error_dev = huge(1.0_mpc)
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
        S = X(1) - 0.1
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
        S = X(3) - 0.1
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
        S = X(5) - 0.1
    end function
    
    
    pure function sect_x_3_body_deriv(X) result (dS)
        real(mpc), dimension(:) :: X
        real(mpc), dimension(size(X)) :: dS
        intent(in) :: X
        dS = 0
        dS(5) = 1.0_mpc
    end function


    pure function sect_w_1(X) result (S)
        real(mpc), dimension(:),intent(in) :: X
        real(mpc) :: S
        S = X(1) - 0.0
    end function
    pure function sect_w_1_deriv(X) result (dS)
        real(mpc), dimension(:) :: X
        real(mpc), dimension(size(X)) :: dS
        intent(in) :: X
        dS = 0
        dS(1) = 1.0_mpc
    end function


    pure function sect_w_2(X) result (S)
        real(mpc), dimension(:),intent(in) :: X
        real(mpc) :: S
        S = X(2) - 0.0
    end function
    pure function sect_w_2_deriv(X) result (dS)
        real(mpc), dimension(:) :: X
        real(mpc), dimension(size(X)) :: dS
        intent(in) :: X
        dS = 0
        dS(2) = 1.0_mpc
    end function
    
    
    pure function sect_w_3(X) result (S)
        real(mpc), dimension(:),intent(in) :: X
        real(mpc) :: S
        S = X(3) - 0.0
    end function
    pure function sect_w_3_deriv(X) result (dS)
        real(mpc), dimension(:) :: X
        real(mpc), dimension(size(X)) :: dS
        intent(in) :: X
        dS = 0
        dS(3) = 1.0_mpc
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

    
    subroutine init_workspace(f, x0, edge, wsp)
        procedure (fRnR1) :: f
        real(mpc), dimension(:), intent(in) :: x0
        real(mpc), intent(in) :: edge
        real(mpc), dimension(:,:), intent(out) :: wsp
        
        integer :: i, procp
        i = 1
        do while (i <= size(wsp,1))
            call randomfill_arr(wsp(i,:))
            wsp(i,:) = (wsp(i,:) - 0.5_mpc)*2*edge + x0
            call check_workspace(f, wsp(i:i,:), procp)
            if (procp == EXIT_SUCCESS) i = i + 1
        end do
    end subroutine init_workspace


    subroutine check_workspace(f, wsp, rets, printp)
        procedure (fRnR1) :: f
        real(mpc), dimension(:,:), intent(in) :: wsp
        integer, intent(out) :: rets
        real(mpc) :: res
        integer :: i
        logical, optional,intent(in) :: printp
        
        logical :: printp_
        printp_ = .FALSE.
        if (present(printp)) printp_ = printp
        rets = EXIT_SUCCESS
        do i = 1, size(wsp, 1)
            res = f(wsp(i,:))
            if (abs(res - error_dev) < eps) then
                rets = EXIT_FAILURE
                exit
            end if
            if (printp_) write(*,*) res, wsp(i,:)
        end do
    end subroutine check_workspace

    function reduce_workspace(f, wsp) result(res)
        procedure (fRnR1) :: f
        real(mpc), dimension(:,:), intent(in) :: wsp
        real(mpc), dimension(size(wsp,2)) :: res
        real(mpc) :: f_res
        integer :: i
        res = wsp(1, :)
        f_res = f(res)
        do i = 2, size(wsp, 1)
            if (f_res > f (wsp(i,:))) then
                res = wsp(i,:)
                f_res = f(res)
            end if
        end do
    end function
    

    ! smart macros, where are you..(
#define std_sect_assign(kind,i,j,afterword) \
    sections(i)%S => CONCAT(CONCAT(CONCAT(sect_,kind)_,j),afterword); \
    sections(i)%dS => CONCAT(CONCAT(CONCAT(sect_,kind)_,j),afterword)_deriv


    subroutine imporove_inipos(integ,t0,x0,t1,fd)
        class(Integrator) :: integ
        real(mpc)  :: t1, t0,x0(:), x(size(x0)-1), delta, d1,d2,cell
        integer    :: i, N, fd, fd_, rets, k, j
        intent(in) :: t1, t0, fd
        optional   :: fd
        integer    :: validp, mixp, skipn, wspsize
        real(mpc), allocatable, dimension(:,:)  :: wsp
        logical    :: advanced_descent_p, cycleshiftp
        
        fd_ = stdout
        if(present (fd)) fd_ = fd
        
        mixp = EXIT_SUCCESS; validp = EXIT_SUCCESS
        
        skipn = 0
        wspsize = 10
        cycleshiftp = .FALSE.

        std_sect_assign(w,1,1,)
        std_sect_assign(w,2,2,)
        std_sect_assign(w,3,3,)
        
        allocate(wsp(wspsize,size(x0)-1))

        ! extra check
        do i = 1, Nbodies
            if (.not. associated(sections(i)%S) &
                &.or. .not. associated(sections(i)%dS)) &
                & error stop "Define section list before improving!"
        end do
        
        associate(wsp_id=>x0(1:size(x0)-1))!, f=>intersection_diff_scalar)
        
        call init_workspace(intersection_diff_scalar,x0(1:size(x0)-1), 0.001_mpc, wsp)
        write(*,*) "-- created random workspace"
        call check_workspace(intersection_diff_scalar, wsp, validp)
        if (validp == EXIT_FAILURE) then
            error stop " !- invalid workspace"
        end if
        write(*,*) "-- check successful"
        
        write(*,*) "-- processing workspace"
        do i = 1, size(wsp,1)
            skipn = 0
            do j = 1, 1
                wsp(i,:) = conjgraddesc(intersection_diff_scalar,wsp(i,:),500,retstat=mixp)
                write(*,*) '?>', intersection_diff_scalar(wsp(i,:))
!                 skipn = skipn + 2
!             wsp(i,:) = newtonhessianfree(intersection_diff_scalar, wsp(i,:), 30, report_fd=stdout)
            end do
        end do
        call check_workspace(intersection_diff_scalar, wsp, validp, .TRUE.)
        write(*,*) "-- reducing workspace"
        x0(1:size(x0)-1) = reduce_workspace(intersection_diff_scalar, wsp)
        write(*,*) intersection_diff_scalar(x0(1:size(x0)-1))
!         call plot_RnR1_section(intersection_diff_scalar, &
!             &wsp_id, (wsp - wsp_id),&
!             &0.0001_mpc, N=30, fd=stdout)

!         x0(1:2) = newtonhessianfree(test, zeron(2), 30, report_fd=stdout)
        end associate
        
    contains 
#if 0
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
#endif
        ! as a special value, return -1
        function intersection_diff_scalar(xp) result(dev)
            real(mpc), intent(in), dimension(:) :: xp
            real(mpc) :: dev
            real(mpc) :: temp(6, size(x0)), norms(6)
            integer :: i, cd, retstat, b = 1, dir, sectn
            cd = size(x0)
            call integ%set_inicond((/xp,t0/), t0)
            call integ%crewind()
            retstat = EXIT_SUCCESS 

            ! пропуск нескольких оборотов
            shift: do i = 1, skipn
                current%S => sections(1)%S
                current%dS => sections(1)%dS
                call shift_to_intersect(integ, t1, retstat)
                call integ%step()
                call shift_to_intersect(integ, t1, retstat)
                call integ%step()
                call shift_to_intersect(integ, t1, retstat)
                call integ%step()
                call shift_to_intersect(integ, t1, retstat)
                call integ%step()
                if (retstat == EXIT_FAILURE) then 
                    go to 505! to avoid useless computations
                end if
            end do shift
            
            comp: do i = 1, 2*Nbodies
                sectn = mod(i-1, Nbodies) + 1
                current%S => sections(sectn)%S
                current%dS => sections(sectn)%dS
                
!                 workaround to get right intersections
                if (i == Nbodies+1) then
                    call shift_to_intersect(integ,t1,retstat)
                    call integ%step()
                    call shift_to_intersect(integ,t1,retstat)
                    call integ%step()
                    call shift_to_intersect(integ,t1,retstat)
                    call integ%step()
                end if

                call shift_to_intersect(integ, t1, retstat)
                temp(i, :) = integ%val()
                norms(i) = norm2(temp(i, 1:shapespDim))

                if (cycleshiftp) then
                    dir = mod(i, Nbodies) - 1
                    temp(i, 1:cd-1) = cyclic_shift_bodies(temp(i, 1:cd-1), dir)
                end if
                if (retstat == EXIT_FAILURE) then 
                    go to 505
                end if
            end do comp


            dev =  &
                &( &
                & + norm2(temp(6, 1:cd-1) - temp(3, 1:cd-1))**(2)&
                & + norm2(temp(5, 1:cd-1) - temp(2, 1:cd-1))**(2)&
                & + norm2(temp(4, 1:cd-1) - temp(1, 1:cd-1))**(2)&
                & + sum(norms - 1)**2 &
!                 & + norm2(temp(2, 1:cd-1) - temp(1, 1:cd-1))**(2)&
!                 & + norm2(temp(3, 1:cd-1) - temp(2, 1:cd-1))**(2)&
                &)
            return

505         dev = error_dev
            stop 'EOT'
            return
        end function 
    
#if 0
        function shapediff(xp) result(dev)
            real(mpc), intent(in), dimension(:) :: xp
            real(mpc) :: dev
            integer :: i, cd, retstat, b = 1, dir, sectn
            real(mpc) :: temp(3, size(x0))
            cd = size(x0)
            
            call integ%set_inicond((/xp,t0/), t0)
            call integ%crewind()
            retstat = EXIT_SUCCESS 
            
            do i = 1, 2*Nbodies
                sectn = mod(i-1, Nbodies) + 1
                current%S => sections(sectn)%S
                current%dS => sections(sectn)%dS
                
                call shift_to_intersect(integ, t1, retstat)
                temp(i, :) = integ%val()
               
!                 dir = mod(i, Nbodies) - 1
!                 temp(i, 1:cd-1) = cyclic_shift_bodies(temp(i, 1:cd-1), dir)
                
                if (retstat == EXIT_FAILURE) then 
                    exit comp
                end if
            end do comp
            if (retstat == EXIT_SUCCESS) then 
                dev =  &
                &( &
                & + norm2(temp(3, 1:cd-1) - temp(3, 1:cd-1))**(2)&
                & + norm2(temp(5, 1:cd-1) - temp(2, 1:cd-1))**(2)&
                & + norm2(temp(4, 1:cd-1) - temp(1, 1:cd-1))**(2)&
                & + norm2(temp(2, 1:cd-1) - temp(1, 1:cd-1))**(2)&
                & + norm2(temp(3, 1:cd-1) - temp(2, 1:cd-1))**(2))
            else
                dev = error_dev
                stop 'EOT'
            endif
        end function shapediff
#endif

        function test(x) result(f)
            real(mpc) , intent(in) :: x(:)
            real(mpc) :: f
            f = 10000.0_mpc*(x(1) - 1.0_mpc)**2  + 10000*(x(2)-2.0_mpc)**2  
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
