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

subroutine ProgramDefinedGeometry( Geom )
    use g_Geometry
    use g_Line
    implicit none
    class(Geometry),intent(in):: Geom
    
    !// locals 
    type(real3) p1, p2, p3, p4
    type(PlaneWall) plane
    type(CylinderWall) cyl,cyl1,cyl2
    logical res
    
    ! cylinder shell, the radius is 0.1 m    
    res =  cyl%CreateCylinder( 0.1_RK, 0.1_RK, p_line( real3(0.0, 0.0, 0.0), real3(0.0,0.0, 0.03) ),24, 1, 1 )
    call Geom%add_Cylinder( cyl ) 

    ! rear wall 
    res =  cyl1%CreateCylinder( 0.001_RK, 0.1_RK, p_line( real3(0.0, 0.0, -0.00001), real3(0.0,0.0, 0.0) ),24, 1, 1 )
    call Geom%add_Cylinder( cyl1 )

    ! front wall 
    res =  cyl2%CreateCylinder( 0.1_RK, 0.001_RK, p_line( real3(0.0, 0.0, 0.03), real3(0.0,0.0, 0.03001) ),24, 1, 1 )
    call Geom%add_Cylinder( cyl2 )
    
end subroutine 
