module Inival
    use Const
    use Celmech
    use ShapeTransform
    
    private

    public :: assign_inicond, cleanup_inicond
    public :: Period, h, finalDim

    real(kind=mpc),parameter :: t0=0.00_mpc   
    real(mpc), parameter :: Period = 6.326_mpc;
    real(kind=mpc),parameter :: t1=19.5*Period          ! Конец интервала интегрирования (начало=0)
    integer,parameter :: ad_ord=6                 ! Порядок для методов Адамса
    integer,parameter :: D = tDim+1
    real(kind=mpc),parameter :: h=0.01_mpc           ! Шаг интегрирования
!     integer, parameter :: D = 2*spcdim*Nbodies                 ! Размерность системы
    real(kind=mpc),dimension(D),parameter :: x0_ideal =(/&
        0.97000436_mpc, -0.24308753_mpc,&
        -0.97000436_mpc, 0.24308753_mpc,&
        0.0_mpc, 0.0_mpc,&
        0.5*0.93240737_mpc, 0.5*0.86473146_mpc,&
        0.5*0.93240737_mpc, 0.5*0.86473146_mpc,&
        -0.93240737_mpc, -0.86473146_mpc,&
        t0&
        /) ! Начальные условия задачи Коши
    real(kind=mpc),dimension(D),parameter :: x0 =(/&
        0.97_mpc, -0.24_mpc,&
        -0.97_mpc, 0.24_mpc,&
        0.0_mpc, 0.0_mpc,&
        0.5*0.93_mpc, 0.5*0.86_mpc,&
        0.5*0.93_mpc, 0.5*0.86_mpc,&
        -0.93_mpc, -0.86_mpc,&
        t0&
        /) ! Начальные условия задачи Коши
    
    integer, parameter :: finalDim = shapespDim*2

    procedure(assign_ordinary_inicond), pointer :: assign_inicond => assign_shape_inicond

contains
    subroutine assign_ordinary_inicond(x_id, x_sh, t_b, t_e)
        real(mpc), intent(out), allocatable, dimension(:) :: x_id, x_sh
        real(mpc), intent(out) :: t_b, t_e
        allocate(x_id(tDim+1), x_sh(tDim+1))
        x_id = x0_ideal
        x_sh = x0
        t_b  = t0
        t_e  = t1
    end subroutine assign_ordinary_inicond


    subroutine assign_shape_inicond(x_id, x_sh, t_b, t_e)
        real(mpc), intent(out), allocatable, dimension(:) :: x_id, x_sh
        real(mpc), intent(out) :: t_b, t_e
        real(mpc), dimension(tDim) :: xtemp
        allocate(x_id(2*shapespDim+1), x_sh(2*shapespDim+1))
        
        xtemp = shapecoords(x0_ideal(1:tDim))
        x_id(1:shapespDim*2) = xtemp(1:shapespDim*2); x_id(size(x_id)) = t0
        
        xtemp = shapecoords(x0(1:tDim))
        x_sh(1:shapespDim*2) = xtemp(1:shapespDim*2); x_sh(size(x_id)) = t0
        
        t_b    = t0
        t_e    = t1
    end subroutine assign_shape_inicond
    
    
    subroutine cleanup_inicond(x_id, x_sh)
        real(mpc), allocatable, dimension(:) :: x_id, x_sh
        deallocate(x_id)
        deallocate(x_sh)
    end subroutine
end module Inival
