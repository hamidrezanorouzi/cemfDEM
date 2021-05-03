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
    type(CylinderWall) cyl
    logical res
    
    ! Creates a cylindrical shell with radius 7.5 cm and length 20 cm whose main axis is y-axis. 
    ! the property type is 1
    ! the user_id of this shell wall is 1
    res =  cyl%CreateCylinder( 0.075_RK, 0.075_RK, p_line( real3(0.0, 0.0, 0.0), real3(0.0,0.2,0.0) ),24, 1, 1 )
    call Geom%add_Cylinder( cyl ) 
    
    
    ! a plane with width of 15 cm with normal vector of (0,1,0).
    ! this plane is placed at the bottom of the cylinder shell. 
    ! the property type is 1
    ! the user_id of this plane is 1
    p1 = real3( -0.075,   0, -0.075)
    p2 = real3( -0.075,   0,  0.075)
    p3 = real3(  0.075,   0,  0.075)
    p4 = real3(  0.075,   0, -0.075)   
    call Geom%add_PlaneWall( p1, p2, p3, p4, 1, 1 )
        
    
end subroutine 
