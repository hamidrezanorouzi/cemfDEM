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
!  File name:  g_Prtcl_SimDomain.f90 
!  Module name: g_Prtcl_SimDomain 
! 
!  Purpose: 
!    1) Providing definition for simulation domain (The region in which the
!       entire simulation should be done)
! 
!  External literature used: none  
! 
!  Referenced modules:
!           - g_TypeDef
!                      
!------------------------------------------------------------------------------
    
module g_Prtcl_SimDomain
    
    use g_TypeDef
    implicit none
    
    
    type SimDomain
    private    
        type(real3) min_domain ! minimum limit of simulation domain
        type(real3) max_domain ! maximum limit of simulation domain
        
    contains
        procedure:: setDomain     => SimDomain_set
        procedure:: getMinDomain  => SimDomain_getMinDomain
        procedure:: getMaxDomain  => SimDomain_getMaxDomain

    end type
    

contains
    

!**************************************************************!
!******************* SimDomain methods ************************!

    subroutine simDomain_set( this, minl, maxl )
        implicit none
        class(SimDomain)::         this
        type( real3),intent(in) :: minl, maxl

        this%min_domain = minl
        this%max_domain = maxl

    end subroutine

    type(real3) function SimDomain_getMinDomain(this)
        implicit none
        class(SimDomain):: this

        SimDomain_getMinDomain = this%min_Domain

    end function
    
    type(real3) function SimDomain_getMaxDomain(this)
        implicit none
        class(SimDomain):: this

        SimDomain_getMaxDomain = this%max_Domain

    end function


!*************** end of SimDomain methods *********************!
!**************************************************************!
    
end module g_Prtcl_SimDomain
