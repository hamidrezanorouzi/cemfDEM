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

function User_Mark( id, flag, ptype, dpos ) result (mark)
    use g_TypeDef
    implicit none 
    
    integer(IK),intent(in):: id, flag, ptype
    type(real4),intent(in):: dpos 
    integer(IK) mark
        
    !// locals
    real minx, maxx, dx
    integer(IK) numx
    real(RK) dist
    real(RK) d1, d2, d3
    type(real3) pos, point1, point2, point3
    
    ! position and diameter of regions which are marked as tracer
    pos = dpos
    d1 = 0.025
    d2 = 0.025
    d3 = 0.025
    point1 = real3(0.0, -0.018, 0.0)
    point2 = real3(0.02, -0.042, 0.0)
    point3 = real3(0.043, -0.066, 0.0)
    
    ! finds particles in each region
    if( sqrt((pos%x- point1%x)**2 + (pos%y- point1%y)**2 ) .le. d1/2) then
        mark = 1
    elseif( sqrt((pos%x- point2%x)**2 + (pos%y- point2%y)**2 ) .le. d2/2) then
        mark = 2
    elseif( sqrt((pos%x- point3%x)**2 + (pos%y- point3%y)**2 ) .le. d3/2) then
        mark = 3
    else
        mark = 0
    end if
        
end function
