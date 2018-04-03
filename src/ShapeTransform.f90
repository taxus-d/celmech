module ShapeTransform
    use Celmech
    use Integrators
    implicit none
    
contains
    function jacobicoords(X) result(Z)
        real(mpc), dimension(Nbodies*spcdim) :: X ! \neq 3 NotImplemented
        real(mpc), dimension(size(X)) :: Z
        integer :: i

        real(mpc) :: mu_1, mu_2

        mu_1 = sqrt(m(1)*m(2)/(m(1)+m(2)))
        mu_2 = sqrt( m(3)*(m(1) + m(2)) / (m(1) + m(2) + m(3)) )

        Z(1:2) = mu_1 * (q(1) - q(2))
        Z(3:4) = mu_2 * (q(3) - (m(1)*q(1) + m(2)*q(2)) / (m(1)+m(2)) )
        Z(5:6) = (m(1)*q(1) + m(2)*q(2) + m(3)*q(3))/sum(m)
        contains
            pure function q(i) result(slice)
                integer,value :: i
                real(mpc), dimension(spcdim) :: slice
                slice = X(2*i-1:2*i)
            end function q
        end function jacobicoords

    function hopfcoords(Z) result(w)
        real(mpc), dimension(spcdim*Nbodies) :: Z
        real(mpc), dimension(size(Z)) :: w
        
        ! > 3 NotImplemented 
        w(1) = 0.5_mpc * (abs(cZ(1))**2 - abs(cZ(2))**2)
        w(2) = realpart(cZ(1) * conjg(cZ(2)))
        w(3) = imagpart(cZ(1) * conjg(cZ(2)))
        contains
            function cZ(i) result(complZ)
                integer, value :: i
                complex(mpc) :: complZ
                complZ = cmplx(Z(2*i-1), Z(2*i), kind=mpc)
            end function cZ
    end function hopfcoords

    function shapecoords(X) result(w)
        real(mpc),intent(in), dimension(:) :: X
        ! -spcdim -- jacobi coordinates
        ! -1      -- hopf map
        real(mpc), dimension(size(X)) :: w
        w = hopfcoords(jacobicoords(X))
    end function shapecoords

end module ShapeTransform
