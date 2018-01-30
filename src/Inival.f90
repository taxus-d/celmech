module Inival
    use Const
    use Celmech
    real(kind=mpc),parameter :: t0=0.00_mpc   
    real(mpc), parameter :: Period = 6.326_mpc;
    real(kind=mpc),parameter :: t1=20*Period          ! Конец интервала интегрирования (начало=0)
    integer,parameter :: ad_ord=6                 ! Порядок для методов Адамса
    integer,parameter :: D = tDim+1
    real(kind=mpc),parameter :: h=0.01_mpc           ! Шаг интегрирования
!     integer, parameter :: D = 2*spcdim*Nbodies                 ! Размерность системы
    real(kind=mpc),dimension(D),parameter :: X0_ideal =(/&
        0.97000436_mpc, -0.24308753_mpc,&
        -0.97000436_mpc, 0.24308753_mpc,&
        0.0_mpc, 0.0_mpc,&
        0.5*0.93240737_mpc, 0.5*0.86473146_mpc,&
        0.5*0.93240737_mpc, 0.5*0.86473146_mpc,&
        -0.93240737_mpc, -0.86473146_mpc,&
        t0&
        /) ! Начальные условия задачи Коши
    real(kind=mpc),dimension(D),parameter :: X0 =(/&
        0.97_mpc, -0.24_mpc,&
        -0.97_mpc, 0.24_mpc,&
        0.0_mpc, 0.0_mpc,&
        0.5*0.93_mpc, 0.5*0.86_mpc,&
        0.5*0.93_mpc, 0.5*0.86_mpc,&
        -0.93_mpc, -0.86_mpc,&
        t0&
        /) ! Начальные условия задачи Коши
end module Inival
