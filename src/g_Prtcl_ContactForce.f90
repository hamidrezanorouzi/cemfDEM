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
!  file name  : g_Prtcl_ContactForce.f90
!  module name: g_PrtclContactForce
!                                  
!  Purpose:
!   1) Providing the base class for calculating contact forces and torques
!      between particles
!                                                    
!  To develop new contact force model, you need to derive a new class from the
!    base_ContactForce. We provided two contact force models in this code. The
!    same procedure should be followed to develop new contact force models.    
! 
!------------------------------------------------------------------------------
    
module g_Prtcl_ContactForce
           
    use g_Prtcl_Properties
    use g_Prtcl_ContactList
    use g_Geometry
    
    implicit none
  
    
        
    ! the base class for any contact force model developed in this code 

!// base_ContactForce 
    type, abstract:: base_ContactForce
        
        ! contact force type
        integer(IK) CF_Type
        
        ! Torque model
        integer(IK) CT_Model
        
        ! rolling friction type
        integer(IK) RF_Type
        
        ! current number of particles
        integer(IK) numPrtcl
        
        ! maximum number of particles in the system
        integer(IK) max_nPrtcl
        
        !the integration time step for particles
        real(RK)::  dt = d_prtcl_dt
        
        integer,allocatable:: dummy1
        ! pointers to particles property type, bounding box, translational and rotational velocities
        integer(IK),dimension(:),pointer     :: Prtcl_type
        type(real4),dimension(:),pointer     :: bndg_box
        type(real3),dimension(:),pointer     :: trans_vel
        type(real3),dimension(:),pointer     :: rot_vel
        
        !pointers to force and torques  
        type(real3),dimension(:),pointer     :: force
        type(real3),dimension(:),pointer     :: torque
        
        
        ! contact lists for particle-particle and particle-wall contacts 
        class(base_ContactList),pointer      :: b_PP_ContList
        class(base_ContactList),pointer      :: b_PW_ContList
        
        ! pointer to geometry object 
        class(Geometry),pointer              :: m_Geometry
        
    contains
        
    !****** all methods for the base class ******!
        
        ! The method for initializing contact force 
        procedure:: InitContactForce=> BCF_InitContactForce
        
        ! methods for handling particles 
        procedure:: get_Pos         => BCF_get_pos  
        procedure:: get_vel         => BCF_get_vel  
        procedure:: get_rVel        => BCF_get_rVel  
        procedure:: get_pType       => BCF_get_pType
        
        procedure:: rotate_shearPlane   => BCF_rotate_shearPlane
        
        ! methods for handling particle-particle contact forces
        procedure:: Add_Force       => BCF_Add_Force
        procedure:: Add_Torque      => BCF_Add_Torque
        
        ! methods for handling particle-particle contact forces
        procedure:: Add_Force_PW       => BCF_Add_Force_PW
        procedure:: Add_Torque_PW      => BCF_Add_Torque_PW
        
        ! methods for handling particle-particle contact forces
        procedure:: AllContactForce_PP  => BCF_AllContactForce_PP
        procedure:: AllContactForce_PW  => BCF_AllContactForce_PW
        
        ! setting geometry 
        procedure:: set_Geometry        => BCF_set_Geometry

        
        ! setting new number of partiles (current number of particles)  
        procedure:: set_numPrtcl        => BCS_set_numPrtcl
        
    !***** deferred boundings which should be provided in every class derived from this class *****! 
        procedure(Prtcl_Prtcl_Force),deferred:: ContactForce_PP
        procedure(Prtcl_Wall_Force),deferred:: ContactForce_PW
        
        procedure(base_set_ContactLists),deferred    :: set_ContactLists
        procedure(base_set_PysicalProperty),deferred :: set_PysicalProperty
        
    end type
!// base_ContactForce 


    ! abstract interfaces for deferred bindings 
    abstract interface
    
        subroutine Prtcl_Prtcl_Force(this, ind , lnew )
            import base_ContactForce, IK
            class(base_ContactForce) this
            integer(IK),intent(in):: ind
            logical,intent(in)    :: lnew
        end subroutine
        
        subroutine Prtcl_Wall_Force(this, ind , lnew )
            import base_ContactForce, IK
            class(base_ContactForce) this
            integer(IK),intent(in):: ind
            logical,intent(in)    :: lnew
        end subroutine
        
        subroutine base_set_ContactLists(this , PP_CntctList, PW_CntctList)
            import base_ContactForce, base_ContactList
            class(base_ContactForce) this
            class(base_ContactList), pointer, intent(in)    :: PP_CntctList, PW_CntctList
            
        end subroutine 
        
        subroutine base_set_PysicalProperty(this, PhProps )
            import base_ContactForce, base_PhysicalProperty
            class(base_ContactForce) this
            class(base_PhysicalProperty),pointer,intent(in) :: PhProps
        end subroutine
        
    end interface
    

!*****************************************************************************!    
!***************** definitions of methods ************************************!
!*****************************************************************************!    
  
contains

 subroutine BCF_set_Geometry(this, geom )
        implicit none
        class(base_ContactForce):: this
        class(Geometry),pointer:: geom
        
        this%m_Geometry => geom
        
 end subroutine
 
! setting current number of particles for base_ContactForce class 
 subroutine BCS_set_numPrtcl(this , new_numPrtcl )
    implicit none
    class(base_ContactForce)            this
    integer(IK),intent(in)          ::  new_numPrtcl
    
    this%numPrtcl = new_numPrtcl
    
 end subroutine 
 
 ! Initializing the contact force 
 subroutine BCF_InitContactForce(this, CF_Type , CT_Model, &
                                 max_nPrtcl, numPrtcl, dt, &
                                 Prtcl_Type, bndg_box, trans_vel, &
                                 rot_vel,  force, torque )
    implicit none
    class(base_ContactForce)                           this
    integer(IK),intent(in)                          :: CF_Type , CT_Model, max_nPrtcl, numPrtcl
    real(RK),intent(in)                             :: dt
    integer(IK),dimension(:),pointer,   intent(in)  :: Prtcl_Type
    type(real4),dimension(:),pointer,   intent(in)  :: bndg_box
    type(real3),dimension(:),pointer,   intent(in)  :: trans_vel, rot_vel, force, torque
    
    this%CF_Type    = CF_Type
    this%CT_Model   = CT_Model
    this%max_nPrtcl = max_nPrtcl
    this%numPrtcl   = numPrtcl
    this%dt         = dt
    this%Prtcl_Type => Prtcl_Type
    this%bndg_box   => bndg_box
    this%trans_vel  => trans_vel
    this%rot_vel    => rot_vel
 
    this%force      => force
    this%torque     => torque
 
 end subroutine

 
! returning the position and diameter of particle i  
type(real4) function BCF_get_pos  (this, i)
    implicit none
    class(base_ContactForce)   this
    integer(IK),intent(in)  :: i
    
    BCF_get_pos   = this%bndg_box(i)
end function


! returning translation velocity of particle i
type(real3) function BCF_get_vel  (this, i)
    implicit none
    class(base_ContactForce)   this
    integer(IK),intent(in)  :: i
    
    BCF_get_vel   = this%trans_vel(i)
end function


! returning rotational velocity of particle i
type(real3) function BCF_get_rvel(this, i)
    implicit none
    class(base_ContactForce)   this
    integer(IK),intent(in)  :: i
    
    BCF_get_rvel   = this%rot_vel(i)

end function


    
! returning particle property type     
integer(IK) function BCF_get_pType(this, i)
    implicit none
    class(base_ContactForce)   this
    integer(IK),intent(in)  :: i
    
    BCF_get_pType   = this%Prtcl_type(i)
end function
    

!    looping over all contact pairs (particle-particle) in the particle-particle contact list
! to calculate the contact force and torques. 
subroutine BCF_AllContactForce_PP( this  )
    implicit none
    class(base_ContactForce) this
    
    ! locals
    integer(IK) ind
    logical     lnew
     
    
    !   starting the iteration over all contact pairs. lnew stores the state of contact 
    ! which and be new (true) or old (false). 
    ind = this%b_PP_ContList%start_iter(lnew)    
    
    ! continuing the iteration till the end of list
    do while (ind> 0 )
        
        ! calling a proper function to calculate force 
        call this%ContactForce_PP( ind , lnew )
        
        ! next pair in the list
        ind = this%b_PP_ContList%nextItem(lnew)
    end do
    
end subroutine


!    looping over all contact pairs (particle-wall) in the particle-wall contact list
! to calculate the contact force and torques. 
subroutine BCF_AllContactForce_PW( this )
    implicit none
    class(base_ContactForce) this
       
    ! locals
    integer(IK) ind
    logical     lnew
        
    !   starting the iteration over all contact pairs. lnew stores the state of contact 
    ! which and be new (true) or old (false). 
    ind = this%b_PW_ContList%start_iter(lnew)    
    
    ! continuing the iteration till the end of list
    do while (ind> 0 )
        
        ! calling a proper function to calculate force 
        call this%ContactForce_PW( ind , lnew )
        
        ! next pair in the list
        ind = this%b_PW_ContList%nextItem(lnew)

    end do
    
end subroutine



function BCF_rotate_shearPlane( this, vec, nv ) result (res)
    implicit none
    class(base_ContactForce)    this
    type(real3),intent(in):: vec, nv
    type(real3) res
    
    !//local
    real(RK) norm_vec, norm_a
    type(real3) a
    norm_vec = norm(vec)
    
    res = vec - (( vec .dot. nv ) * nv)
    
    return 
    if( norm_vec .gt. 0.0_RK ) then
        
        a = vec - (( vec .dot. nv ) * nv)
        
        norm_a = norm(a)
        if( norm_a .gt. 0.0_RK ) then
            res = (norm_vec/norm_a)*vec
        else
            ! this is expected to never happens
            res = zero_r3
        end if
        
    else
        ! a zero vector on the shear plane is always zero
        res = zero_r3
                        
    end if
    
    
end function



! adding the calculated contact force to particle pair (i,j)
subroutine BCF_Add_Force( this, pari , parj , Fc )
    implicit none
    class(base_ContactForce)   this
    integer(IK),intent(in)  :: pari, parj
    type(real3),intent(in)  :: Fc
    
    this%force(pari) = this%force(pari) + Fc
    this%force(parj) = this%force(parj) - Fc

end subroutine

! adding the calculated torque to contact pair (i,j)
subroutine BCF_Add_Torque( this, pari, parj , Mij, Mji )
    implicit none
    class(base_ContactForce)   this
    integer(IK),intent(in)  :: pari, parj
    type(real3),intent(in)  :: Mij, Mji
    
    this%torque(pari) = this%torque(pari) + Mij
    this%torque(parj) = this%torque(parj) + Mji

end subroutine

! adding the calculated force (particle-wall) to particle i
subroutine BCF_Add_Force_PW( this, pari,  Fc )
    implicit none
    class(base_ContactForce)   this
    integer(IK),intent(in)  :: pari
    type(real3),intent(in)  :: Fc
    
    this%force(pari) = this%force(pari) + Fc
    
end subroutine

! adding the calculated torque (particle-wall) to particle i
subroutine BCF_Add_Torque_PW( this, pari, Mij )
    implicit none
    class(base_ContactForce)   this
    integer(IK),intent(in)  :: pari
    type(real3),intent(in)  :: Mij
    
    this%torque(pari) = this%torque(pari) + Mij
        
end subroutine


   
end module
