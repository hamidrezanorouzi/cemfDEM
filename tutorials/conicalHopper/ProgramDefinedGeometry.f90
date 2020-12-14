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
    
     !// shell of the hopper, radius = 0.1 m    
    res =  cyl%CreateCylinder( 0.1_RK, 0.1_RK, p_line( real3(0.0, 0.0, 0.0), real3(0.0,0.0, 0.2) ),24, 1, 1 )
    call Geom%add_Cylinder( cyl ) 
    
    !// Conical section of the hopper with a 8-cm exit orifice
    res =  cyl%CreateCylinder( 0.04_RK, 0.1_RK, p_line( real3(0.0, 0.0, -0.059), real3(0.0,0.0,0.0) ),24, 1, 1 )
    call Geom%add_Cylinder( cyl ) 
    
    !// the gate to block the exit orifice, this will be removed after packing stage
    !// the user_id of this wall is 2. This id will be used to remove the wall later in the simulation
    p1 = real3(-0.05, -0.05, -0.059);
    p2 = real3( 0.05, -0.05, -0.059);
    p3 = real3( 0.05,  0.05, -0.059);
    p4 = real3(-0.05,  0.05, -0.059);
    call Geom%add_PlaneWall( p1, p2, p3, p4, 2, 1, .true. )
    
    !// flat plate under the hopper to collect particles for pile formation
    p1 = real3(-0.25, -0.25, -0.2);
    p2 = real3( 0.25, -0.25, -0.2);
    p3 = real3( 0.25,  0.25, -0.2);
    p4 = real3(-0.25,  0.25, -0.2); 
    call Geom%add_PlaneWall( p1, p2, p3, p4, 3, 1, .true. )
        
end subroutine 
