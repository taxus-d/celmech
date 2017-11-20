module Inival
    use Prec
    use Celmech
    real(kind=mpc),parameter :: t0=0.0_mpc   
    real(mpc), parameter :: Period = 6.326_mpc;
    real(kind=mpc),parameter :: t1=20*Period          ! Конец интервала интегрирования (начало=0)
    integer,parameter :: ad_ord=6                 ! Порядок для методов Адамса
    integer,parameter :: D = tDim+1
    real(kind=mpc),parameter :: h=0.1_mpc           ! Шаг интегрирования
!     integer, parameter :: D = 2*spcdim*Nbodies                 ! Размерность системы
!     real(kind=mpc),dimension(D),parameter :: X0 =(/&
!         0.97000436, -0.24308753,&
!         -0.97000436, 0.24308753,&
!         0.0, 0.0,&
!         0.5*0.93240737, 0.5*0.86473146,&
!         0.5*0.93240737, 0.5*0.86473146,&
!         -0.93240737, -0.86473146&
!         /) ! Начальные условия задачи Коши
    real(kind=mpc),dimension(D),parameter :: X0 =(/&
        0.97, -0.24,&
        -0.97, 0.24,&
        0.0, 0.0,&
        0.5*0.93, 0.5*0.86,&
        0.5*0.93, 0.5*0.86,&
        -0.93, -0.86,&
        t0&
        /) ! Начальные условия задачи Коши
    interface
        function eq_fun(tt,X) result(f)
            use Prec; use Celmech
            integer,parameter :: D = tDim +1
            real(mpc):: tt, X(D)
            intent(in) :: tt
            real(mpc) :: f(D)
        end function eq_fun
    end interface
    procedure (eq_fun), pointer :: f => poincare_section_eq
end module Inival
