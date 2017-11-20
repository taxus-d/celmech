module funct
    implicit none
    integer,parameter :: kin=8  ! Параметр разновидности типа для real
    real(kind=kin),parameter :: T=100.0            ! Конец интервала интегрирования (начало=0)
    integer,parameter :: n=6                 ! Порядок для методов Адамса
    integer,parameter :: D=3                 ! Размерность системы
    real(kind=kin),dimension(D),parameter :: X0=(/10.0,10.0,10.0/) ! Начальные условия задачи Коши
    real(kind=kin),parameter :: h=0.0001           ! Шаг интегрирования

contains

    function f(tt,X)
        implicit none
        real(kind=kin),intent(in) :: tt
        real(kind=kin),dimension(D),intent(in) :: X
        real(kind=kin),dimension(D) :: f


        f(1)=10.0*(X(2)-X(1))
        f(2)=X(1)*(28.0-X(3))-X(2)
        f(3)=X(1)*X(2)-8.0*X(3)/3.0
    end function f

end module funct
