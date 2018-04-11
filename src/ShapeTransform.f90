module ShapeTransform
    use Celmech
    implicit none
    
    real(mpc) :: mu_1 = sqrt( m(1)*m(2)/(m(1)+m(2)) ), &
        &mu_2 = sqrt( m(3)*(m(1) + m(2)) / (m(1) + m(2) + m(3)) )
    
contains
    function jtransf(X) result(Z)
        real(mpc), dimension(Nbodies*spcdim), intent(in) :: X
        real(mpc), dimension(Nbodies*spcdim) :: Z
        
        Z(1:2) = mu_1 * (q(1) - q(2))! \neq 3 NotImplemente! \neq 3 NotImplementedd
        Z(3:4) = mu_2 * (q(3) - (m(1)*q(1) + m(2)*q(2)) / (m(1)+m(2)) )
        Z(5:6) = (m(1)*q(1) + m(2)*q(2) + m(3)*q(3))/sum(m)
    contains
        pure function q(i) result(slice)
            integer,value :: i
            real(mpc), dimension(spcdim) :: slice
            slice = X(2*i-1:2*i)
        end function q
    end function jtransf

    function jacobicoords(X) result(Z)
        real(mpc), dimension(:), intent(in) :: X 
        real(mpc), dimension(size(X)) :: Z
        Z(1:tDim/2) = jtransf(X(1:tDim/2))
        Z(tDim/2+1:tDim) = jtransf(X(tDim/2+1:tDim))
    end function jacobicoords


    function hopfcoords(Z) result(w)
        real(mpc), dimension(tDim) :: Z
        real(mpc), dimension(size(Z)) :: w
        
        real(mpc), dimension(spcdim) :: z1, z2, dz1, dz2
        
        w = 0
        ! > 3 NotImplemented 
        
        z1  = Z(1:2)
        z2  = Z(3:4)
        dz1 = Z(7:8)
        dz2 = Z(9:10)

        ! coords
        w(1) = 0.5_mpc * (norm2(z1)**2 - norm2(z2)**2)
        w(2) = dot_product(z1,z2)
        w(3) = simpl_product(z1,z2)
        
        ! velocities
        w(4) =   dot_product(dz1,z1) -   dot_product(dz2,z2)
        w(5) =   dot_product(z1,dz2) +   dot_product(dz1,z2)
        w(6) = simpl_product(z1,dz2) + simpl_product(dz1,z2)
!         w(2) = realpart(cZ(1) * conjg(cZ(2)))
!         w(3) = imagpart(cZ(1) * conjg(cZ(2)))
        contains
            pure function simpl_product(a,b) result(p)
                real(mpc), dimension(2), intent(in) :: a, b ! а всё что > 2 -- не скаляр
                real(mpc) :: p
                p = dot_product( (/a(2), -a(1)/), b)
            end function simpl_product
    end function hopfcoords

    function shapecoords(X) result(w)
        real(mpc),intent(in), dimension(:) :: X
        real(mpc), dimension(size(X)) :: w
        w = hopfcoords(jacobicoords(X))
    end function shapecoords

    function shapesphcoords(X) result(w)
        real(mpc),intent(in), dimension(:) :: X
        real(mpc), dimension(size(X)) :: w
        w = hopfcoords(jacobicoords(X))
        w = w / norm2(w)
    end function shapesphcoords

end module ShapeTransform
