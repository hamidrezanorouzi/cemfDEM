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
!  File name:  g_TypeDef.f90 
!  Module name: g_TypeDef
! 
!  Purpose:
!    1) General definition of intrinsic and derived types that are frequently
!       used in the program.   
! 
!    2) Various operators are overloads for derived types such as sum, product,
!       cross and dot product.
! 
! 
!  External literature used: none              
! 
!  Referenced modules: none                    
! 
!------------------------------------------------------------------------------

module g_TypeDef

    implicit none
    
    ! determining whether 32 or 64 bit integer variable is used
    ! RK = 4 for 32 bit and RK = 8 for 64 bit
    integer,parameter:: IK = 4 ! better to use 4 only
       
    ! Determining whether 32 or 64 bit floating point variable is used.
    ! IK = 4 for 32 bit (single precision) and IK = 8 for 64 bit (double precision)
    integer,parameter:: RK = 8
    
    
    real(RK),parameter:: Pi = 3.14159265359_RK    ! Pi constant
    real(RK),parameter:: epsln = 0.000001_RK      ! a very small value compared to 1.0
    
    
    type real3
        real(RK)    x
        real(RK)    y
        real(RK)    z
    end type

    type integer3
        integer(IK) x
        integer(IK) y
        integer(IK) z
    end type

    type real4
        real(RK)    x
        real(RK)    y
        real(RK)    z
        real(RK)    w
    end type
    
    type irange
        integer(IK) imin
        integer(IK) imax
    end type 
    
    ! a zero vector of type real3
    type(real3),parameter:: zero_r3 = real3(0.0_RK, 0.0_RK, 0.0_RK)
    type(real4),parameter:: zero_r4 = real4(0.0_RK, 0.0_RK, 0.0_RK, 0.0_RK)
    
        
!********************************************************************!
!*********** operator overloading for real3 derived type ************!
    
    interface operator(+)
        module procedure real3_add1 , real3_add2, real3_add2_2
    end interface

    interface operator(-)
        module procedure real3_sub
    end interface

    interface operator(*)      
        module procedure real3_prod1 , real3_prod2, real3_prod2_2
    end interface

    interface operator(/)
        module procedure real3_division1, real3_division2, &
                         real3_division1_2, real3_division2_2
    end interface
    
    interface assignment(=)
        module procedure real3_real4, integer3_real3 
    end interface
    
    ! overloading of dot product operator for two vectors
    interface operator(-)
        module procedure real4_sub
    end interface
    
    interface operator( .dot. )
        module procedure real3_dotProd, real3_dotProd2
    end interface

    ! overloading of cross product operator for two vectors
    interface operator( .cross. )
        module procedure real3_CrossProd
    end interface
    
    ! Whether two real3 vectors are equal or not 
    interface operator( .equal. )
        module procedure real3_equal
    end interface
    
    interface min
        module procedure min_v
    end interface
    
    interface max
        module procedure max_v
    end interface
!*********** operator overloading for real3 derived type ************!    
!********************************************************************!    


!********************************************************************!    
!*********** operator overloading for real4 derived type ************!

    interface operator( .ovlp. )
        module procedure real4_ovlp
    end interface
    
    interface operator (.dist.)
        module procedure real4_dist, real3_dist
    end interface
    
    interface operator(.nv.)
        module procedure real4_nv, real3_nv
    end interface
    
    interface norm
        module procedure real4_norm, real3_norm
    end interface
    
!*********** operator overloading for real4 derived type ************!
!********************************************************************!    

    interface rotx
        module procedure real3_rotx
    end interface
    
    interface roty
        module procedure real3_roty
    end interface
    
    interface rotz
        module procedure real3_rotz
    end interface
       
!********************************************************************!
!********* operator overloading for integer3 derived type ***********!    
    
    interface operator(+)
        module procedure integer3_add
    end interface
    
!********* operator overloading for integer3 derived type ***********!     
!********************************************************************!    
     
    
contains ! module procedures

!************** all methods for real3 derived type ******************!
    !**** sum of two vectors ****
    function real3_add1( p1, p2 ) result( p3)
    
        implicit none
        type(real3),intent(in)::  p1, p2
        type(real3) p3
    
        p3 = real3( p1%x+p2%x , p1%y + p2%y , p1%z + p2%z ) 

    end function

    function real3_add2( p1, v ) result( p3)
    
        implicit none
        type(real3),intent(in)::  p1
        real(RK),intent(in)::     v
        type(real3) p3
    
        p3 = real3( p1%x+ v , p1%y + v , p1%z + v ) 

    end function

    function real3_add2_2( v , p1 ) result( p3)
        implicit none
        type(real3),intent(in)::  p1
        real(RK),intent(in)::     v
        type(real3) p3
    
        p3 = real3( p1%x+ v , p1%y + v , p1%z + v ) 

    end function


    !**** sub of two vectors ****
    function real3_sub( p1, p2 ) result( p3)
        implicit none
        type(real3),intent(in)::  p1, p2
        type(real3) p3
    
        p3 = real3( p1%x-p2%x , p1%y - p2%y , p1%z - p2%z ) 

    end function


    !**** product of two vectors ****
    function real3_prod1( p1, p2 ) result(p3)
        implicit none
        type(real3),intent(in):: p1, p2
        type(real3):: p3

        p3 = real3( p1%x * p2%x, p1%y * p2%y, p1%z * p2%z)

    end function

    function real3_prod2( p1, v ) result( p3)
        implicit none
        type(real3),intent(in)::  p1
        real(RK),intent(in)::     v
        type(real3) p3
    
        p3 = real3( p1%x * v , p1%y * v , p1%z * v ) 

    end function

    function real3_prod2_2( v , p1 ) result( p3)
        implicit none
        type(real3),intent(in)::  p1
        real(RK),intent(in)::     v
        type(real3) p3
    
        p3 = real3( p1%x * v , p1%y * v , p1%z * v ) 

    end function

    !**** division operations ****
    function real3_division1(p1,p2) result (p3)
        implicit none
        type(real3),intent(in)::p1, p2
        type(real3) p3

        p3 = real3( p1%x/p2%x , p1%y/p2%y, p1%z/p2%z )

    end function

    function real3_division2( p1, v) result(p3)
        implicit none
        type(real3),intent(in)  :: p1
        real(RK),intent(in)::      v
        type(real3) p3

        p3 = real3( p1%x/v , p1%y/v, p1%z/v )

    end function

    function real3_division1_2(p1,p2) result (p3)
        implicit none
        type(real3),intent(in),dimension(:)         :: p1
        type(real3),intent(in)                      :: p2
        type(real3),dimension(size(p1,1))           :: p3

        p3(:)%x = p1(:)%x/p2%x
        p3(:)%y = p1(:)%y/p2%y
        p3(:)%z = p1(:)%z/p2%z

    end function

    function real3_division2_2( p1, v) result(p3)
        implicit none
        type(real3),intent(in),dimension(:) :: p1
        real(RK),intent(in)::                  v
        type(real3),dimension(size(p1,1))   :: p3

        p3(:)%x = p1(:)%x/v
        p3(:)%y = p1(:)%y/v
        p3(:)%z = p1(:)%z/v

    end function

    ! assignment operator for x = y
    subroutine real3_real4(x , y)
        implicit none
        type(real3),intent(out) :: x
        type(real4),intent(in)  :: y
        x = real3(y%x, y%y, y%z)
    end subroutine
    
    
    subroutine integer3_real3( x, y)
        implicit none
        type(integer3),intent(out)  :: x
        type(real3),intent(in)      :: y
        x = integer3(y%x , y%y, y%z)
    end subroutine
    
    !**** dot product of two vectors ****
    function real3_dotProd( p1, p2) result( v )
        implicit none
        type(real3),intent(in):: p1, p2
        real(RK)::  v

        v = p1%x*p2%x + p1%y*p2%y + p1%z*p2%z

    end function

    function real3_dotProd2( p1, p2) result( v )
        implicit none
        type(real3),intent(in),dimension(:)         :: p1
        type(real3),intent(in),dimension(size(p1,1)):: p2
        real(RK),dimension(size(p1,1))::    v

        v(:) = p1(:)%x*p2(:)%x + p1(:)%y*p2(:)%y + p1(:)%z*p2(:)%z

    end function

    !**** cross product of two vectors **** 
    function real3_CrossProd( u, v ) result( p )
        implicit none
        type(real3),intent(in):: u, v
        type(real3) p

        p = real3( u%y*v%z - u%z*v%y, &
                    u%z*v%x - u%x*v%z, &
                    u%x*v%y - u%y*v%x )

    end function
    
    function real3_dist(p1, p2 ) result(dist)
        
        implicit none
        type(real3),intent(in):: p1, p2
        real(RK) dist
        
        dist = norm(p1-p2)
            
    end function
    
    ! **** vector operations **** 
    function real3_nv(p1, p2) result(nv)
        
        implicit none
        type(real3),intent(in):: p1, p2
        type(real3) nv
        type(real3) s
        s = p1-p2
        nv = s / norm(s)
        
    end function
    
    function real3_norm(p1) result (norm)
        implicit none
        type(real3),intent(in):: p1
        real(RK) norm
        
        norm = sqrt(p1%x*p1%x + p1%y*p1%y + p1%z*p1%z)
        
    end function
    
    function real3_equal(p1, p2 ) result(eq)
        
        implicit none
        type(real3),intent(in):: p1, p2
        logical eq
        
        if(abs(p1%x-p2%x) > epsln .or. abs(p1%y-p2%y) > epsln .or. abs(p1%z-p2%z) > epsln ) then
            eq = .false.    
            return
        end if
        
        eq =.true.
        return
        
        
    end function
    
    !**** Rotations in right-handed xyz coordinates **** 
    function real3_rotx(p1 , tet) result( rp )
        implicit none
        type(real3),intent(in)  :: p1
        real(RK), intent(in)    :: tet
        
        type(real3) rp
               
        rp = real3( p1%x , p1%y*cos(tet) - p1%z*sin(tet) , p1%y*sin(tet) + p1%z*cos(tet) )
              
    end function
    
    function real3_roty(p1 , tet) result( rp )
    
        implicit none
        type(real3),intent(in)  :: p1
        real(RK), intent(in)    :: tet
        
        type(real3) rp
               
        rp = real3( p1%x* cos(tet) + p1%z*sin(tet) , p1%y , -p1%x* sin(tet) + p1%z*cos(tet) )
              
    end function
    
    function real3_rotz(p1 , tet) result( rp )
    
        implicit none
        type(real3),intent(in)  :: p1
        real(RK), intent(in)    :: tet
        
        type(real3) rp
               
        rp = real3( p1%x* cos(tet) - p1%y*sin(tet) , p1%x* sin(tet) + p1%y*cos(tet) , p1%z )
              
    end function
        
!************ end of all methods for real3 derived type ************!


!************ all methods for real4 derived type *******************!
    function real4_sub( p1, p2 ) result( p3)
        implicit none
        type(real4),intent(in)::  p1, p2
        type(real4) p3
    
        p3 = real4( p1%x-p2%x , p1%y - p2%y , p1%z - p2%z, p1%w - p2%w ) 

    end function    

    function real4_ovlp(p1, p2 ) result(ovlp)
        
        implicit none
        type(real4),intent(in):: p1, p2
        real(RK) ovlp
                
         ovlp = 0.5_RK * (p1%w+p2%w) - norm(p1-p2)
               
    
    end function

    function real4_dist(p1, p2 ) result(dist)
        
        implicit none
        type(real4),intent(in):: p1, p2
        real(RK) dist
        
        dist = norm(p1-p2)
    
    end function
    
    function real4_nv(p1, p2) result(nv)
        
        implicit none
        type(real4),intent(in):: p1, p2
        type(real3) nv
        
        type(real3) s
        
        s = p1-p2
        
        nv = s / norm(s)
        
    end function
    
    function real4_norm(p1) result (norm)
        implicit none
        type(real4),intent(in):: p1
        real(RK) norm
        
        norm = sqrt(p1%x*p1%x + p1%y*p1%y + p1%z*p1%z)
        
    end function
    
!************End of all methods for real3 derived type *************!
    

!************ all methods for integer3 derived type ****************!
    !**** sum of two vectors ****
    function integer3_add( p1, p2 ) result( p3)
        implicit none
        type(integer3),intent(in)::  p1, p2
        type(integer3) p3
    
        p3 = integer3( p1%x+p2%x , p1%y + p2%y , p1%z + p2%z ) 

    end function
!********** end of all methods for integer3 derived type ***********!    


    function det22(mat) result(res)
        implicit none
        real(RK),dimension(2,2),intent(in):: mat
        real(RK) res
        
        res = mat(1,1)*mat(2,2) - (mat(1,2)*mat(2,1))
        
    end function
    
    
    function det33(mat) result(res)
        implicit none
        real(RK),dimension(3,3),intent(in):: mat
        real(RK) res
        
        res = mat(1,1)*( mat(2,2)*mat(3,3)- mat(3,2)*mat(2,3) ) - & 
              mat(1,2)*( mat(2,1)*mat(3,3)- mat(3,1)*mat(2,3) ) + &
              mat(1,3)*( mat(2,1)*mat(3,2)- mat(3,1)*mat(2,2) )
        
    end function

    function Solve33(A,B) result(res)
        implicit none
        real(RK),dimension(3,3),intent(in):: A
        real(RK),dimension(3)  ,intent(in):: B
        real(RK),dimension(3)             :: res
        
        real(RK) detA
        real(RK),dimension(3,3) :: temp
        
        detA = det33(A)
        
        ! first equaiton, x
        temp = A
        temp(:,1) = B
        res(1) = det33(temp)/detA
        
        ! second equation, y
        temp = A
        temp(:,2) = B
        res(2) = det33(temp)/detA
        
        ! second equation, y
        temp = A
        temp(:,3) = B
        res(3) = det33(temp)/detA
        
    end function
    
    
    function min_v( v1, v2) result(v3)
        
        implicit none
        type(real3),intent(in):: v1, v2
        type(real3) v3
        
        v3 = real3( min(v1%x,v2%x) , min(v1%y,v2%y), min(v1%z,v2%z) )
        
    end function
    
    function max_v( v1, v2) result(v3)
        
        implicit none
        type(real3),intent(in):: v1, v2
        type(real3) v3
        
        v3 = real3( max(v1%x,v2%x) , max(v1%y,v2%y), max(v1%z,v2%z) )
        
    end function
    
    
    
    function RPMtoRAD_S( rpm ) result( res )
        implicit none
        real(RK), intent(in) :: rpm
        real(RK) res
        
        res = 2*Pi*rpm/60
        
    end function
    
    subroutine vtk_header( nUnit , mes )
        implicit none
        integer(IK),intent(in):: nUnit
        character(*),intent(in):: mes
    
        write( nUnit, "(A)")"# vtk DataFile Version 2.0"
        write( nUnit, "(A)") trim( mes )
        write( nUnit, "(A)") "ASCII"
    
    end subroutine 
    
    
end module
