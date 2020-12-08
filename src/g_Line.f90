!-------------- cemfDEM code: a dem simulation code ----------------------------
!      D        C enter of
!     D M       E ngineering and
!    D   M      M ultiscale modeling of
!   D     M     F luid flow    
!  EEEEEEEEM    .ir
!------------------------------------------------------------------------------
!  Copyright (C): cemf
!  website: www.cemf.ir
!------------------------------------------------------------------------------  
!  This file is part of cemfDEM code. It is a free software for simulating 
!  granular flow. You can redistribute it and/or modify it under the terms of
!  GNU General Public License version 3 or any other later versions. 
! 
!  cemfDEM code is distributed to help others in their research in the field  
!  of granular flow, but WITHOUT ANY WARRANTY; without even the implied 
!  warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!------------------------------------------------------------------------------
!  file name  : g_Line.f90
!  module name: g_Line
!                     
!  Purpose:
!   A class to define a parametric line with two points p1 and p2. The line has
!    the direction p1 ---> p2
!
!------------------------------------------------------------------------------
    
module g_Line
    
    use g_TypeDef
    implicit none
    
    !// p_line
    type p_line
        type(real3) p1, p2
        type(real3) v
        real(RK) dist
    contains
    
    procedure:: setLine         => pl_setLIne
    procedure:: getPoint        => pl_getPoint
    procedure:: get_Project     => pl_get_Project
    procedure:: get_vect_n      => pl_get_vect_n
    procedure:: get_vect        => pl_get_vect
    procedure:: get_length      => pl_get_length
    procedure:: find_intersect  => pl_find_intersect
    procedure:: rotate_point    => pl_rotate_point
    
    end type
    !// p_line
    
    
    ! constructor for p_line
    interface p_line
        procedure p_line_const
    end interface
    
    !// transfer to z-axis 
    
    type,extends(p_line):: zAxisLine
    
        type(p_line) orgn_line
        logical :: in_aAxis = .false.
        
        real(RK),dimension(4,4):: Trans_P1, Trans_xz, Trans_z, Trans_z_xz_P1
        real(RK),dimension(4,4):: ITrans_P1, ITrans_xz, ITrans_z, ITrans_P1_xz_z
    contains
        procedure:: Rotate_toZ1
        procedure:: Rotate_toZ2
        generic  :: Rotate_toZ  =>  Rotate_toZ2
        
        procedure:: rBack_fromZ2
        generic  :: rBack_fromZ => rBack_fromZ2
        
        procedure,nopass:: zRotate     => zAL_zRotate
        
        
        procedure,private:: TransformationMatrices
    end type
    
contains

!constructor for p_line
function p_line_const( p1, p2 ) result (this)
    implicit none
    type(real3),intent(in):: p1, p2
    type(p_line) this
    
    call this%setLine(p1,p2)
end function


!**********************************************************************
!setting two points of line and calculates other parameters
!**********************************************************************
subroutine pl_setLine( this, p1, p2)
    
    implicit none
    class(p_line) this
    type(real3) p1, p2
    
    !// locals
    type(real3) v
    
    !// body
    
    ! p1 ----> p2
    v = p2 - p1 
    
    this%v = v
    this%p1 = p1
    this%p2 = p2
    this%dist = p1.dist.p2
    
end subroutine 

!**********************************************************************
! returning a point on the line based on the parameter t
!  0 <= t <= 1 : between points p1 and p2
!  t <= 0      : along the line before p1
!  t >= 1      : along the line after p2
!**********************************************************************
type(real3) function pl_getPoint( this, t )
    implicit none
    class(p_line)     this
    real(RK),intent(in):: t
    
    pl_getPoint = this%v * t + this%p1

end function


real(RK) function pl_get_Project(this, p )
    implicit none
    class(p_line)   this
    type(real3),intent(in):: p
    
    !// locals
    type(real3) w
    
    ! locals
    w = p - this%p1 
    pl_get_Project = w .dot. this%v / (this%v .dot. this%v)
    
end function

!**********************************************************************
! returning the normalized direction vector of the line 
!**********************************************************************
type(real3) function pl_get_vect_n(this)
    implicit none
    class(p_line) this
    
    pl_get_vect_n = this%v / this%dist
    
end function

!**********************************************************************
! returning the direction vector of the line 
!**********************************************************************
type(real3) function pl_get_vect(this)
    implicit none
    class(p_line) this
    
    pl_get_vect = this%v 
    
end function    

!**********************************************************************
! returning the length of the line
!**********************************************************************
real(RK) function pl_get_length(this)
    implicit none
    class(p_line) this
    
    pl_get_length = this%dist
    
end function

!**********************************************************************
! finding the intersect of the line with Line L2
!**********************************************************************
type(real3) function pl_find_intersect(this, L2 , res)
    implicit none
    class(p_line) this
    type(p_line), intent(in) :: L2
    integer(IK),      intent(out):: res
    
    !// locals
    real(RK) h, k
    type(real3) vec, crs1, crs2
    
    !// body
    vec = L2%P1 - this%P1 ! the vector that connects two points of two lines
    
    crs1 = L2%v .cross. vec
    crs2 = L2%v .cross. this%v
    
    h = norm(crs1)
    k = norm(crs2)
    
    if ( h .lt.  0.0000000001_RK .or. k .lt. 0.0000000001_RK ) then
        ! the two lines do not have intersect with each other 
        res = -1
        return
    end if
    
    if( (crs1 .dot. crs2) < 0.0 ) then
        
        ! pointing to different directions
        pl_find_intersect = this%P1 - (h/k)*this%v
        
    else
        
        ! pointing to the same direction
        pl_find_intersect = this%P1 + (h/k)*this%v
        
    end if
    
    res = 0
    
end function

!**********************************************************************
! Rotating p0 around the line by angle theta and returning the rotated point 
!**********************************************************************
type(real3) function pl_rotate_point( this, p0 , teta )
    implicit none
    class(p_line) this
    type(real3),intent(in) :: p0
    real(RK),   intent(in) :: teta
    
    !// locals
    type(real3) nv, res, p1
    real cos_tet , sin_tet, u2, v2, w2
    
    !// body
    nv = this%get_vect_n()
    cos_tet = cos(teta)
    sin_tet = sin(teta)
    u2 = nv%x**2.0_RK
    v2 = nv%y**2.0_RK
    w2 = nv%z**2.0_RK
    p1 = this%p1
    
    !          (a(v2+w2) - u( bv + cw - ux - vy - wz)) (1-cos_tet) + x cos_tet + (- cv + bw - wy + vz)sin_tet 
    res%x = ( p1%x*(v2+w2) - (nv%x*(p1%y*nv%y + p1%z*nv%z - nv%x*p0%x - nv%y*p0%y - nv%z*p0%z)) )*(1-cos_tet)  + &
                p0%x * cos_tet + &
            ( -p1%z*nv%y + p1%y*nv%z - nv%z*p0%y + nv%y*p0%z ) * sin_tet
    
    !          ( b(u2+w2) - v( au + cw - ux - vy - wz))(1-cos_tet) + y cos_tet + ( cu - aw + wx - uz ) sin_tet
    res%y = ( p1%y*(u2+w2) - (nv%y*(p1%x*nv%x + p1%z*nv%z - nv%x*p0%x - nv%y*p0%y - nv%z*p0%z)) )*(1-cos_tet)  + &
                p0%y * cos_tet + &
            ( p1%z*nv%x - p1%x*nv%z + nv%z*p0%x - nv%x*p0%z ) * sin_tet
    
    
    !         (c(u2+v2) - w( au + bv - ux - vy - wz ))(1-cos_tet) + z cos_tet + (-bu + av - vx + uy) sin_tet
    res%z = ( p1%z*(u2+v2) - (nv%z*(p1%x*nv%x + p1%y*nv%y - nv%x*p0%x - nv%y*p0%y - nv%z*p0%z)) )*(1-cos_tet)  + &
                p0%z * cos_tet + &
            ( -p1%y*nv%x + p1%x*nv%y - nv%y*p0%x + nv%x*p0%y ) * sin_tet
    
    
    pl_rotate_point = res
    
    end function


!////////////////////////////////////////////////////////////////////////////////

subroutine Rotate_toZ1( this, line )
    implicit none
    class( zAxisLine )              this
    type( p_line), intent(in):: line
    
    !// locals
    type(real3) lp1, lp2
    
    this%orgn_line = line
    call this%TransformationMatrices()
    
    lp1 = this%Rotate_toZ(this%orgn_line%p1)
    lp2 = this%Rotate_toZ(this%orgn_line%p2)
    
    call this%setLine( lp1, lp2 )
    
end subroutine 


function Rotate_toZ2( this , p ) result (res)
    implicit none
    class( zAxisLine )        this
    type(real3), intent(in):: p
    type(real3) res
    
    real(RK),dimension(4):: point
    
    point = (/p%x, p%y, p%z, 1.0_RK/)
    point = matmul(this%Trans_z_xz_P1, point )
    !
    
    res = real3( point(1), point(2), point(3) )
    
end function

function rBack_fromZ2(this, p) result (res)
    implicit none
    class( zAxisLine )        this
    type(real3), intent(in):: p
    type(real3) res
    
    real(RK),dimension(4):: point
    
    point = (/p%x, p%y, p%z, 1.0_RK/)
    point = matmul( this%ITrans_P1_xz_z, point)
    
    res = real3( point(1), point(2), point(3) )
    
end function

function zAL_zRotate( p , teta ) result(res)
    implicit none
    type(real3),intent(in):: p
    real(RK),   intent(in):: teta
    type(real3) res
    
    res = rotz( p, teta )

end function

subroutine TransformationMatrices(this)
    implicit none
    class(zAxisLine) this
    
    !//locals
    real(RK) u, v, w, u2v2
    type(real3) p1, nv
    
    
    p1 = this%orgn_line%p1
    
    ! they are inserted row-major, while FORTRAN needs column major.
    ! so, transpose is used here
    this%Trans_P1 = transpose( reshape(         &
               (/1.0_RK, 0.0_RK, 0.0_RK, -p1%x, &
                 0.0_RK, 1.0_RK, 0.0_RK, -p1%y, &
                 0.0_RK, 0.0_RK, 1.0_RK, -p1%z, &
                 0.0_RK, 0.0_RK, 0.0_RK, 1.0_RK /) , (/4,4/) ) )
    
    nv = this%orgn_line%get_vect_n()
    
    u = nv%x
    v = nv%y
    w = nv%z
    
    u2v2 = sqrt(u*u + v*v )
    u2v2 = max(0.000000000000001_RK , u2v2 ) 

    this%Trans_xz = transpose( reshape(                  &
                    (/ u/u2v2 , v/u2v2, 0.0_RK, 0.0_RK, &
                       -v/u2v2, u/u2v2, 0.0_RK, 0.0_RK, &
                       0.0_RK , 0.0_RK, 1.0_RK, 0.0_RK, &
                       0.0_RK , 0.0_RK, 0.0_RK, 1.0_RK /), (/4,4/) ) )
    
    ! u2 + v2+w2 = 1
    this%Trans_z = transpose( reshape(          &
                    (/ w     , 0.0_RK, -u2v2 , 0.0_RK , &
                      0.0_RK , 1.0_RK, 0.0_RK, 0.0_RK , &
                      u2v2   , 0.0_RK, w     , 0.0_RK , &
                      0.0_RK , 0.0_RK, 0.0_RK, 1.0_Rk/), (/4,4/) ) )
    
    this%Trans_z_xz_P1 = matmul( this%Trans_z, matmul( this%Trans_xz , this%Trans_P1 ) )
    
    ! 
    this%ITrans_P1 = transpose( reshape(         &
               (/1.0_RK, 0.0_RK, 0.0_RK, +p1%x, &
                 0.0_RK, 1.0_RK, 0.0_RK, +p1%y, &
                 0.0_RK, 0.0_RK, 1.0_RK, +p1%z, &
                 0.0_RK, 0.0_RK, 0.0_RK, 1.0_RK /) , (/4,4/) ) )
    
    this%ITrans_xz = transpose( reshape(                  &
                    (/ u/u2v2 , -v/u2v2, 0.0_RK, 0.0_RK, &
                       +v/u2v2, u/u2v2, 0.0_RK, 0.0_RK, &
                       0.0_RK , 0.0_RK, 1.0_RK, 0.0_RK, &
                       0.0_RK , 0.0_RK, 0.0_RK, 1.0_RK /), (/4,4/) ) )
    
    ! u2 + v2+w2 = 1
    this%ITrans_z = transpose( reshape(          &
                    (/ w     , 0.0_RK, +u2v2 , 0.0_RK , &
                      0.0_RK , 1.0_RK, 0.0_RK, 0.0_RK , &
                      -u2v2   , 0.0_RK, w     , 0.0_RK , &
                      0.0_RK , 0.0_RK, 0.0_RK, 1.0_Rk/), (/4,4/) ) )
    
    this%ITrans_P1_xz_z = matmul( this%ITrans_P1, matmul( this%ITrans_xz , this%ITrans_z ) )
    
    ! correcting the transformation matrix in the case of coincidence with z-axis
    if( abs(w-1.0_RK)<0.000000001 ) then
        
        this%Trans_z_xz_P1  = this%Trans_P1
        this%ITrans_P1_xz_z = this%ITrans_P1
        
    end if
    
end subroutine 

end module
