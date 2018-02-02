module Integrators
    use IO
    use IO_array
    use Debug
    use Poly
    use NewtSolve

!     logical :: weirdstep ! TODO: get rid of
    implicit none
   
    logical :: printp

    type, abstract :: Integrator
        integer :: dimen, ord, cache_size
        real(mpc), allocatable :: x0(:)
        real(mpc) :: t0, t, h
        procedure(int_f), pointer, nopass :: f
        contains
            procedure(int_set_inicond) , deferred :: set_inicond
            procedure(int_crewind)     , deferred :: crewind
            procedure(int_get_cache)   , deferred :: get_cache
            procedure(int_step)        , deferred :: step
            procedure(int_time)        , deferred :: time
            procedure(int_val)         , deferred :: val
    end type
    abstract interface 
        function int_f(t,x) result(f)
            import :: mpc
            real(mpc),intent(in) :: x(:), t
            real(mpc) :: f(size(x))
        end function int_f

        subroutine int_set_inicond(self, x0, t0)
            import Integrator,mpc
            class(Integrator),intent(inout) :: self
            real(mpc) :: x0(:), t0
        end subroutine int_set_inicond
        subroutine int_crewind(self)
            import Integrator, mpc
            class(Integrator), intent(inout) :: self
        end subroutine int_crewind
        
        function int_get_cache(self) result(c)
            import Integrator, mpc
            class(Integrator), intent(in) :: self
            real(mpc) :: c(self%cache_size, self%dimen)
        end function int_get_cache
        
        subroutine int_step(self,h)
            import Integrator,mpc
            class(Integrator),intent(inout) :: self
            real(mpc), optional, intent(in) :: h
        end subroutine int_step
        
        function int_time(self) result(t)
            import Integrator,mpc
            class(Integrator), intent(in) :: self
            real(mpc) :: t
        end function int_time
        function int_val(self) result(v)
            import Integrator,mpc
            class(Integrator), intent(in) :: self
            real(mpc) :: v(self%dimen)
        end function int_val
    end interface  


    type, extends(Integrator) :: RungeKuttaInt
        real(mpc), allocatable :: x(:)
    contains
        procedure :: set_inicond => runge_set_inicond
        procedure :: crewind     => runge_crewind
        procedure :: get_cache   => runge_get_cache
        procedure :: step        => runge_step
        procedure :: time        => runge_time
        procedure :: val         => runge_val
        final :: clear_runge
    end type
    interface RungeKuttaInt
        module procedure :: init_runge
    end interface RungeKuttaInt
    
    
    type, extends(Integrator) :: ExAdamsInt
        real(mpc), allocatable :: Mc(:), Xc(:,:), fc(:,:)
    contains
        procedure :: set_inicond => ex_adams_set_inicond
        procedure :: crewind     => ex_adams_crewind
        procedure :: get_cache   => ex_adams_get_cache
        procedure :: step        => ex_adams_step
        procedure :: time        => ex_adams_time
        procedure :: val         => ex_adams_val
        final :: clear_ex_adams
    end type
    interface ExAdamsInt
        module procedure :: init_ex_adams
    end interface ExAdamsInt
    
    
    type, extends(Integrator) :: ImAdamsInt
        real(mpc), allocatable :: Mc(:), Xc(:,:), fc(:,:)
    contains
        procedure :: set_inicond => im_adams_set_inicond
        procedure :: crewind     => im_adams_crewind
        procedure :: get_cache   => im_adams_get_cache
        procedure :: step        => im_adams_step
        procedure :: time        => im_adams_time
        procedure :: val         => im_adams_val
        final :: clear_im_adams
    end type
    interface ImAdamsInt
        module procedure :: init_im_adams
    end interface ImAdamsInt

contains
! **************************
! Runge-Kutta integrator
! **************************
    subroutine runge_set_inicond(self,x0,t0)
        class(RungeKuttaInt), intent(inout) :: self
        real(mpc) :: x0(:), t0
        self%x0 = x0
        self%t0 = t0
    end subroutine runge_set_inicond

    
    subroutine runge_crewind(self)
        class(RungeKuttaInt), intent(inout) :: self
        self%x = self%x0
        self%t = self%t0
    end subroutine runge_crewind
    
    
    function init_runge(f,x0,t0,step) result(self)
        real(mpc), intent(in) :: x0(:), t0, step
        type(RungeKuttaInt) :: self
        procedure(int_f) :: f 

        self%dimen = size(x0)
        self%ord = 1
        self%cache_size = 1 ! cannot have zero extent
        self%h = step
        self%f => f
        !call step_still_alive()
        allocate(self%x0(self%dimen))
        allocate(self%x(self%dimen))
        
        !call step_still_alive()
        call self%set_inicond(x0,t0)
        call self%crewind()
        !call step_still_alive()
    end function init_runge 
    
    
    function runge_get_cache(self) result(c)
        class(RungeKuttaInt), intent(in) :: self
        real(mpc) :: c(self%cache_size, self%dimen)
        
        c(1,:) = self%x0
    end function runge_get_cache

    function runge_time(self) result(t)
        class(RungeKuttaInt), intent(in) :: self
        real(mpc) :: t
        t = self%t 
    end function runge_time
    
    function runge_val(self) result(v)
        class(RungeKuttaInt), intent(in) :: self
        real(mpc) :: v(self%dimen)
        v = self%x
    end function runge_val


    function internal_euler_step(f,t, x, h) result(x1)
        procedure(int_f) :: f
        real(mpc), intent(in) :: t, x(:), h
        real(mpc) :: x1(size(x))
        x1 = x + h * f(t,x)
    end function internal_euler_step
    

    function internal_runge_step(f,t, x, h) result(x1)
        procedure(int_f) :: f
        real(mpc), intent(in) :: t, x(:), h
        real(mpc) :: x1(size(x)), k(4, size(x))

        k(1,:) = h*f(t, x) 
        k(2,:) = h*f(t + h*0.5_mpc, x + 0.5_mpc*k(1,:)) 
        k(3,:) = h*f(t + h*0.5_mpc, x + 0.5_mpc*k(2,:)) 
        k(4,:) = h*f(t + h, x + k(3,:)) 
        x1 = x + 1.0_mpc/6.0_mpc * & 
            (k(1,:) + 2.0_mpc*k(2,:) + 2.0_mpc*k(3,:) + k(4,:))
    end function internal_runge_step
    
    
    subroutine runge_step (self, h)
        class(RungeKuttaInt), intent(inout) :: self
        real(mpc), optional, intent(in) :: h
        real(mpc) :: h_
        h_ = self%h
        if (present(h)) h_ = h
        
        self%x = internal_runge_step(self%f,self%t, self%x, h_)
        self%t = self%t + h_
    end subroutine runge_step
    
    
    subroutine clear_runge(self)
        type(RungeKuttaInt), intent(inout) :: self
        deallocate(self%x0)
        deallocate(self%x)
    end subroutine clear_runge


! **************************
! Explicit Adams integrator
! **************************
    subroutine ex_adams_set_inicond(self, x0, t0)
        class(ExAdamsInt), intent(inout) :: self
        real(mpc) :: x0(:), t0
        self%Xc(0,:) = x0
        self%t0 = t0
    end subroutine ex_adams_set_inicond
    
    subroutine ex_adams_crewind(self)
        class(ExAdamsInt), intent(inout) :: self
        integer :: j
        self%t = self%t0

        associate(n=>self%ord)
        do j=1,n-1
            self%t = self%t + self%h
            self%Xc(j,:) = internal_runge_step(self%f, self%t, self%Xc(j-1,:), self%h)
        end do
        do j=0, n-1
            self%fc(j,:) = self%f(self%t, self%Xc(n-1-j, :))
            self%t = self%t - self%h
        end do
        end associate
    end subroutine ex_adams_crewind

    function init_ex_adams(f,x0,t0,ord,step) result(self)
        type(ExAdamsInt) :: self
        real(mpc) :: t0, step,x0(:)
        integer :: ord, n, j
        intent(in) :: ord,x0,t0,step
        procedure(int_f) :: f

        self%ord = ord
        self%cache_size = ord
        self%h =  step
        self%dimen = size(x0)
        self%f => f

        allocate (self%Mc(0:self%ord-1))
        allocate (self%Xc(0:self%ord-1, size(x0)))
        allocate (self%fc(0:self%ord-1, size(x0)))

        associate(n => self%ord)
        do j = 0, n-1
            self%Mc(j) = A(n, j) 
        end do
        
        call self%set_inicond(x0, t0)
        call self%crewind()
        
        self%t = self%t0 + self%h*(n-1)
        end associate
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
    end function init_ex_adams

    
    subroutine clear_ex_adams(self)
        type(ExAdamsInt) :: self
        if (allocated(self%Mc)) deallocate (self%Mc)
        if (allocated(self%Xc)) deallocate (self%Xc)
        if (allocated(self%fc)) deallocate (self%fc)
    end subroutine clear_ex_adams
    
    
    function ex_adams_get_cache(self) result(c)
        class(ExAdamsInt), intent(in) :: self
        real(mpc) :: c(self%cache_size, self%dimen)
        
        c(:,:) = self%Xc
    end function ex_adams_get_cache

    
    
    subroutine ex_adams_step(self,h)
        class(ExAdamsInt), intent(inout) :: self
        real(mpc), optional, intent(in) :: h
        real(mpc):: h_
        integer :: n
        h_ = self%h
        if (present(h)) h_ = h
        
        n = self%ord
        self%Xc(0:n-2,:) = self%Xc(1:n-1,:)
        self%Xc(n-1,:) = self%Xc(n-1,:) + h_ * matmul(self%Mc, self%fc) 
        self%t = self%t + h_
        self%fc(1:n-1,:) = self%fc(0:n-2,:)
        self%fc(0,:)     = self%f(self%t, self%Xc(n-1,:))
    end subroutine ex_adams_step

    
    function ex_adams_time(self) result(t)
        class(ExAdamsInt), intent(in) :: self
        real(mpc) :: t
        t = self%t 
    end function ex_adams_time
    
    
    function ex_adams_val(self) result(x)
        class(ExAdamsInt), intent(in) :: self
        real(mpc) :: x(self%dimen)
        x = self%Xc(self%ord-1,:) 
    end function ex_adams_val

! **************************
! Implicit Adams integrator
! **************************
    subroutine im_adams_set_inicond(self, x0, t0)
        class(ImAdamsInt), intent(inout) :: self
        real(mpc) :: x0(:), t0
        self%Xc(0,:) = x0
        self%t0 = t0
    end subroutine im_adams_set_inicond
    
    subroutine im_adams_crewind(self)
        class(ImAdamsInt), intent(inout) :: self
        integer :: j
        self%t = self%t0

        associate(n=>self%ord)
        do j=1, n-1
            self%t = self%t + self%h
            self%Xc(j,:) = internal_runge_step(self%f, self%t, self%Xc(j-1,:), self%h)
        end do
        do j=-1, n-2
            self%fc(j,:) = self%f(self%t0+self%h*(n-2-j), self%Xc(n-2-j, :))
        end do
        end associate
    end subroutine im_adams_crewind
    function init_im_adams(f, x0, t0, ord, step) result(self)
        type(ImAdamsInt) :: self
        real(mpc) :: t0, step,x0(:)
        integer :: ord, n, j
        intent(in) :: ord,x0,t0,step
        procedure(int_f) :: f

        self%ord = ord
        self%cache_size = ord-1
        self%h =  step
        self%dimen = size(x0)
        self%f => f

        allocate (self%Mc(-1:self%ord-2))
        allocate (self%Xc(0:self%ord-1, size(x0)))
        allocate (self%fc(-1:self%ord-2, size(x0)))
        
        associate (n=>self%ord)
        do j=-1,n-2
            self%Mc(j) = B(n, j)
        end do
        end associate
        
        call self%set_inicond(x0, t0)
        call self%crewind()
    contains
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
    end function init_im_adams

    subroutine im_adams_step(self,h)
        class(ImAdamsInt), intent(inout) :: self
        real(mpc), optional, intent(in) :: h
        real(mpc):: h_
        h_ = self%h
        if (present(h)) h_ = h

        associate(n=>self%ord)
        self%Xc(n-1  , :) = newtitsolve(impl, self%Xc(n-2, :),100)
        self%Xc(0:n-2, :) = self%Xc(1:n-1, :)
        self%t = self%t + h_
        self%fc(-1   , :) = self%f(self%t, self%Xc(n-2,:))
        self%fc(0:n-2, :) = self%fc(-1:n-3,:)
        end associate
    contains
        ! но иначе нужно переписывать Ньютона
        function impl(x) result(y)
            real(mpc) :: x(:) , y(size(x))
            intent(in) :: x
            y = x  - h_*self%Mc(-1)*self%f(self%t+h_,x) - self%Xc(self%ord-2,:) - &
                h_*matmul(self%Mc(0:self%ord-2), self%fc(0:self%ord-2,:))
        end function impl
    end subroutine im_adams_step
    
    
    function im_adams_get_cache(self) result(c)
        class(ImAdamsInt), intent(in) :: self
        real(mpc) :: c(self%cache_size, self%dimen)
        
        c(:,:) = self%Xc(0:self%cache_size-1, :)
    end function im_adams_get_cache
    
    
    function im_adams_time(self) result(t)
        class(ImAdamsInt), intent(in) :: self
        real(mpc) :: t
        t = self%t
    end function im_adams_time
    
    
    function im_adams_val(self) result(x)
        class(ImAdamsInt), intent(in) :: self
        real(mpc) :: x(self%dimen)
        x = self%Xc(self%ord-2,:)
    end function im_adams_val
    
    
    subroutine clear_im_adams(self)
        type(ImAdamsInt) :: self
        deallocate (self%Mc)
        deallocate (self%Xc)
        deallocate (self%fc)
    end subroutine clear_im_adams


end module Integrators

