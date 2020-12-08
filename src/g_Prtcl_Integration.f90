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
!  file name  : g_Prtcl_Integration.f90 
!  module name: g_Prtcl_Integration 
! 
!  Purpose: 
!    1) to provide classes to perform integration on linear and rotational
!       equation of motions for spherical particles.
! 
! 
!********************************************************************!
    
module g_Prtcl_Integration

    
    use g_Prtcl_DefaultValues
    use g_Prtcl_OneStepIntegration
    use g_Prtcl_MultiPointIntegration
    use g_Prtcl_TwoStepIntegration
    
    
    implicit none
    
    
    !// prtcl_Integration
    type prtcl_Integration
    
        integer(IK):: Int_Method = d_PI_Method ! initially gets the default value
        integer(IK)   max_nPrtcl
        integer(IK)   numPrtcl
        real(RK)   :: dt = d_prtcl_dt          ! initially gets the default value
        
        ! pointers for particles velocity, pos, acceleration and flag
        integer(IK),pointer,dimension(:)    :: prtcl_flag
        type(real4),pointer,dimension(:)    :: Pos
        type(real3),pointer,dimension(:)    :: vel
        type(real3),pointer,dimension(:)    :: acc
        
        ! objects for performing integration
        
        ! one-step methods
        type(ForwardEuler),allocatable          :: m_ForwardEuler
        type(ModifiedEuler),allocatable         :: m_ModifiedEuler
        type(Taylor2),allocatable               :: m_Taylor2
        type(Taylor3),allocatable,dimension(:)  :: m_Taylor3
        type(Taylor4),allocatable,dimension(:)  :: m_Taylor4
        
        ! multi-point methods
        type(AB2_params),allocatable,dimension(:):: m_AB2
        type(AB3_params),allocatable,dimension(:):: m_AB3
        type(AB4_params),allocatable,dimension(:):: m_AB4
        type(AB5_params),allocatable,dimension(:):: m_AB5
        
        ! predictor-corrector methods
        type(AB2AM3_params),allocatable,dimension(:):: m_AB2AM3
        type(AB3AM4_params),allocatable,dimension(:):: m_AB3AM4
        type(AB4AM5_params),allocatable,dimension(:):: m_AB4AM5
        type(Gear3_params),allocatable,dimension(:) :: m_Gear3
        type(Gear4_params),allocatable,dimension(:) :: m_Gear4
        type(Gear5_params),allocatable,dimension(:) :: m_Gear5
        
    contains
    
    ! Initializing the object 
    procedure:: Init_Integration    => PI_Init_Integration
    
    ! performing prediction step
    procedure:: predict             => PI_predict
    
    ! performing correction step
    procedure:: correct             => PI_correct
    
    ! setting new/current number of particles in the system
    procedure:: set_numPrtcl        => PI_set_numPrtcl
    
    ! final/destructor
    final:: PI_Finalize
    
    ! allocating memory 
    procedure,private:: AllocateMethod  => PI_AllocateMethod
    
    end type !// prtcl_Integration
    
    
    ! rotational motion
    !// prtcl_rot_Integration    
    type, extends(Prtcl_Integration):: Prtcl_rot_Integration
        
	integer,allocatable:: dummy3
        type(real3),pointer,dimension(:)    :: Pos3
        
    contains
        
        ! Initializing the object
        procedure:: Init_rot_Integration    => PRI_Init_Integration
        
        ! setting current/new number of particles
        procedure:: set_numPrtcl            => PRI_set_numPrtcl
            
    end type !// prtcl_rot_Integration


contains

!**********************************************************************
! Initializing the prtcl_Integration object
!**********************************************************************
subroutine PI_Init_Integration(this, Method, max_nPrtcl, numPrtcl, dt, Pos, vel, acc , flag )

    implicit none
    class(Prtcl_Integration)  this
    
    integer(IK),                     intent(in) :: Method, max_nPrtcl, numPrtcl
    real(RK),                        intent(in) :: dt
    type(real4),pointer,dimension(:),intent(in) :: Pos
    type(real3),pointer,dimension(:),intent(in) :: vel, acc
    integer(IK),pointer,dimension(:),intent(in) :: flag
       
    this%numPrtcl   = numPrtcl
    this%max_nPrtcl = max_nPrtcl
    this%Int_Method = Method
    this%dt         = dt
    this%pos        => pos
    this%vel        => vel
    this%acc        => acc
    this%Prtcl_flag => flag
    
    ! allocating memory
    call this%AllocateMethod()
    
end subroutine

!**********************************************************************
! Initializing prtcl_rot_Integration object
!**********************************************************************
subroutine PRI_Init_Integration(this, Method, max_nPrtcl, numPrtcl, dt, Pos, vel, acc, flag )
    implicit none
    class(Prtcl_rot_Integration)                   this
    integer(IK),                     intent(in) :: Method, max_nPrtcl, numPrtcl
    real(RK),                        intent(in) :: dt
    type(real3),pointer,dimension(:),intent(in) :: Pos
    type(real3),pointer,dimension(:),intent(in) :: vel, acc
    integer(IK),pointer,dimension(:),intent(in) :: flag
    
    !// body 
    this%numPrtcl   = numPrtcl
    this%max_nPrtcl = max_nPrtcl
    this%Int_Method = Method
    this%dt         = dt
    this%pos3       => pos
    this%vel        => vel
    this%acc        => acc
    this%Prtcl_flag => flag
    
    ! allocates memory
    call this%AllocateMethod()
    
end subroutine

!**********************************************************************
! prediction step of integration
!**********************************************************************
subroutine PI_predict( this )
    implicit none
    class(Prtcl_Integration)  this
    
    ! // locals
    integer(IK) i
    type(real3) pos, vel
    
    
    !// body

    select case( this%Int_Method )
        
    case( PIM_FE, PIM_ME, PIM_Taylor2 , PIM_Taylor3 , PIM_Taylor4 ) 
        
        return ! one step integration does not have a prediction step
        
    case (PIM_AB2, PIM_AB3 , PIM_AB4, PIM_AB5)
        
        return ! multi-point integration does not have a prediction step
            
    case (PIM_Gear3)
        
        do i = 1, this%numPrtcl
            
            if( this%prtcl_flag(i) >= Pflg_inDomain ) then
            
                call this%m_Gear3(i)%predict( this%dt, pos , this%vel(i) )
                select type(this)
                type is (Prtcl_Integration)
                    this%pos(i) = real4( pos%x , pos%y, pos%z, this%pos(i)%w) 
                type is (Prtcl_rot_Integration)
                    this%pos3(i) = pos
                end select
            
            end if
                
        end do
        
    case(PIM_Gear4)
        
        do i = 1, this%numPrtcl
            
            if( this%prtcl_flag(i) >= Pflg_inDomain ) then
            
                call this%m_Gear4(i)%predict( this%dt, pos , this%vel(i) )
            
                select type(this)
                type is (Prtcl_Integration)
                    this%pos(i) = real4( pos%x , pos%y, pos%z, this%pos(i)%w) 
                type is (Prtcl_rot_Integration)
                    this%pos3(i) = pos
                end select
            
            end if
            
        end do
    
    case(PIM_Gear5)
 
        do i = 1, this%numPrtcl
            
            if( this%prtcl_flag(i) >= Pflg_inDomain ) then
                
                call this%m_Gear5(i)%predict( this%dt, pos , this%vel(i) )
                select type(this)
                type is (Prtcl_Integration)
                    this%pos(i) = real4( pos%x , pos%y, pos%z, this%pos(i)%w) 
                type is (Prtcl_rot_Integration)
                    this%pos3(i) = pos
                end select
            
            end if
                
        end do
 
    
          
    case (PIM_AB2AM3)
        
        do i = 1, this%numPrtcl
            
            if( this%prtcl_flag(i) >= Pflg_inDomain ) then
            
                call this%m_AB2AM3(i)%predict( this%dt, pos , this%vel(i) )
                select type(this)
                type is (Prtcl_Integration)
                    this%pos(i) = real4( pos%x , pos%y, pos%z, this%pos(i)%w) 
                type is (Prtcl_rot_Integration)
                    this%pos3(i) = pos
                end select
                
            end if
            
        end do
        
    case (PIM_AB3AM4)
        
        do i = 1, this%numPrtcl
            
            if( this%prtcl_flag(i) >= Pflg_inDomain ) then
            
                call this%m_AB3AM4(i)%predict( this%dt, pos , this%vel(i) )
                select type(this)
                type is (Prtcl_Integration)
                    this%pos(i) = real4( pos%x , pos%y, pos%z, this%pos(i)%w) 
                type is (Prtcl_rot_Integration)
                    this%pos3(i) = pos
                
                end select
                
            end if
            
        end do        
        
    case (PIM_AB4AM5)
    
        do i = 1, this%numPrtcl
            
            if( this%prtcl_flag(i) >= Pflg_inDomain ) then
            
                call this%m_AB4AM5(i)%predict( this%dt, pos , this%vel(i) )
                select type(this)
                type is (Prtcl_Integration)
                    this%pos(i) = real4( pos%x , pos%y, pos%z, this%pos(i)%w) 
                type is (Prtcl_rot_Integration)
                    this%pos3(i) = pos
                end select
            
            end if
            
        end do        
        
    end select
    
    
end subroutine 


!**********************************************************************
! correction step of integration
!**********************************************************************
subroutine PI_correct( this )
    implicit none
    class(Prtcl_Integration)  this
    
    !// locals
    integer(IK) i
    type(real3) pos, vel
    
    !// body
    
    select case( this%Int_Method )
    
    case( PIM_FE)    
        
        do i = 1, this%numPrtcl
            
            if( this%prtcl_flag(i) >= Pflg_inDomain ) then
                
                select type(this)
                type is (Prtcl_Integration)
                    pos = this%pos(i)
                    call this%m_ForwardEuler%Integrate( this%acc(i), this%dt, pos , this%vel(i) )
                    this%pos(i) = real4( pos%x , pos%y, pos%z, this%pos(i)%w )
                type is (Prtcl_rot_Integration)
                    call this%m_ForwardEuler%Integrate( this%acc(i), this%dt, this%pos3(i) , this%vel(i) )
                end select
            
            end if
            
        end do
    
    case( PIM_ME) 
        
        do i = 1, this%numPrtcl
            
            if( this%prtcl_flag(i) >= Pflg_inDomain ) then
                
                select type(this)
                type is(Prtcl_Integration)
                    pos = this%pos(i)
                    call this%m_ModifiedEuler%Integrate( this%acc(i), this%dt, pos , this%vel(i) )
                    this%pos(i) = real4( pos%x , pos%y, pos%z, this%pos(i)%w)
                type is(Prtcl_rot_Integration)
                    
                    call this%m_ModifiedEuler%Integrate( this%acc(i), this%dt, this%pos3(i) , this%vel(i) )
		    
                end select
                
            end if
            
        end do
        
    case( PIM_Taylor2) 
        
        do i = 1, this%numPrtcl
            
            if( this%prtcl_flag(i) >= Pflg_inDomain ) then
                
                select type(this)
                type is(Prtcl_Integration)
                    pos = this%pos(i)
                    call this%m_Taylor2%Integrate( this%acc(i), this%dt, pos , this%vel(i) )
                    this%pos(i) = real4( pos%x , pos%y, pos%z, this%pos(i)%w) 
                type is(Prtcl_rot_Integration)
                    call this%m_Taylor2%Integrate( this%acc(i), this%dt, this%pos3(i) , this%vel(i) )
                end select 
                
            end if
            
        end do
        
    case (PIM_Taylor3)     
        
        do i = 1, this%numPrtcl
            
            if( this%prtcl_flag(i) >= Pflg_inDomain ) then
                
                select type(this)
                type is(Prtcl_Integration)
                    pos = this%pos(i)
                    call this%m_Taylor3(i)%Integrate( this%acc(i), this%dt, pos , this%vel(i) )
                    this%pos(i) = real4( pos%x , pos%y, pos%z, this%pos(i)%w) 
                type is(Prtcl_rot_Integration)
                    call this%m_Taylor3(i)%Integrate( this%acc(i), this%dt, this%pos3(i) , this%vel(i) )
                end select 
            
            end if
                
        end do
    
    case (PIM_Taylor4)     
        
        do i = 1, this%numPrtcl
            
            if( this%prtcl_flag(i) >= Pflg_inDomain ) then
                
                select type(this)
                type is(Prtcl_Integration)
                    pos = this%pos(i)
                    call this%m_Taylor4(i)%Integrate( this%acc(i), this%dt, pos , this%vel(i) )
                    this%pos(i) = real4( pos%x , pos%y, pos%z, this%pos(i)%w) 
                type is(Prtcl_rot_Integration)
                    call this%m_Taylor4(i)%Integrate( this%acc(i), this%dt, this%pos3(i) , this%vel(i) )
                end select     
            
            end if
            
        end do
    
    case (PIM_AB2)
        
        do i = 1, this%numPrtcl
            
            if( this%prtcl_flag(i) >= Pflg_inDomain ) then
                
                select type(this)
                type is(Prtcl_Integration)
                    pos = this%pos(i)
                    call this%m_AB2(i)%Integrate( this%acc(i), this%dt, pos , this%vel(i) )
                    this%pos(i) = real4( pos%x , pos%y, pos%z, this%pos(i)%w) 
                type is(Prtcl_rot_Integration)
                    call this%m_AB2(i)%Integrate( this%acc(i), this%dt, this%pos3(i) , this%vel(i) )
                end select
            
            end if
            
        end do        
        
    case (PIM_AB3)
        
        do i = 1, this%numPrtcl
            
            if( this%prtcl_flag(i) >= Pflg_inDomain ) then
                
                select type(this)
                type is(Prtcl_Integration)
                    pos = this%pos(i)
                    call this%m_AB3(i)%Integrate( this%acc(i), this%dt, pos , this%vel(i) )
                    this%pos(i) = real4( pos%x , pos%y, pos%z, this%pos(i)%w) 
                type is(Prtcl_rot_Integration)
                    call this%m_AB3(i)%Integrate( this%acc(i), this%dt, this%pos3(i) , this%vel(i) )
                end select
                
            end if
            
        end do
    
    case (PIM_AB4)
        
        do i = 1, this%numPrtcl
            
            if( this%prtcl_flag(i) >= Pflg_inDomain ) then
                
                select type(this)
                type is(Prtcl_Integration)
                    pos = this%pos(i)
                    call this%m_AB4(i)%Integrate( this%acc(i), this%dt, pos , this%vel(i) )
                    this%pos(i) = real4( pos%x , pos%y, pos%z, this%pos(i)%w) 
                type is(Prtcl_rot_Integration)
                    call this%m_AB4(i)%Integrate( this%acc(i), this%dt, this%pos3(i) , this%vel(i) )
                end select
                
            end if
            
        end do
    
    case (PIM_AB5)
        
        do i = 1, this%numPrtcl
            
            if( this%prtcl_flag(i) >= Pflg_inDomain ) then
                
                select type(this)
                type is(Prtcl_Integration)
                    pos = this%pos(i)
                    call this%m_AB5(i)%Integrate( this%acc(i), this%dt, pos , this%vel(i) )
                    this%pos(i) = real4( pos%x , pos%y, pos%z, this%pos(i)%w) 
                type is(Prtcl_rot_Integration)
                    call this%m_AB5(i)%Integrate( this%acc(i), this%dt, this%pos3(i) , this%vel(i) )
                end select
                
            end if
            
        end do
        
    case(PIM_Gear3)        
        
        do i = 1, this%numPrtcl
            
            if( this%prtcl_flag(i) >= Pflg_inDomain ) then
                
                select type(this)
                type is (Prtcl_Integration)
                    call this%m_Gear3(i)%correct( this%acc(i), this%dt, pos , this%vel(i) )
                    this%pos(i) = real4( pos%x , pos%y, pos%z, this%pos(i)%w) 
                type is(Prtcl_rot_Integration)
                    call this%m_Gear3(i)%correct( this%acc(i), this%dt, this%pos3(i) , this%vel(i) )
                end select
                
            end if
            
        end do
        
    case(PIM_Gear4)        
        
        do i = 1, this%numPrtcl
            
            if( this%prtcl_flag(i) >= Pflg_inDomain ) then
                
                select type(this)
                type is (Prtcl_Integration)
                    call this%m_Gear4(i)%correct( this%acc(i), this%dt, pos , this%vel(i) )
                    this%pos(i) = real4( pos%x , pos%y, pos%z, this%pos(i)%w) 
                type is(Prtcl_rot_Integration)
                    call this%m_Gear4(i)%correct( this%acc(i), this%dt, this%pos3(i) , this%vel(i) )
                end select
            
            end if
            
        end do
        
    
    case(PIM_Gear5)        
        
        do i = 1, this%numPrtcl
            
            if( this%prtcl_flag(i) >= Pflg_inDomain ) then
                
                select type(this)
                type is (Prtcl_Integration)
                    call this%m_Gear5(i)%correct( this%acc(i), this%dt, pos , this%vel(i) )
                    this%pos(i) = real4( pos%x , pos%y, pos%z, this%pos(i)%w) 
                type is(Prtcl_rot_Integration)
                    call this%m_Gear5(i)%correct( this%acc(i), this%dt, this%pos3(i) , this%vel(i) )
                end select
            
            end if
            
        end do
         
        
    case(PIM_AB2AM3)
        
        do i = 1, this%numPrtcl
            
            if( this%prtcl_flag(i) >= Pflg_inDomain ) then
            
                select type(this)
                type is (Prtcl_Integration)
                    call this%m_AB2AM3(i)%correct( this%acc(i), this%dt, pos , this%vel(i) )
                    this%pos(i) = real4( pos%x , pos%y, pos%z, this%pos(i)%w) 
                type is(Prtcl_rot_Integration)
                    call this%m_AB2AM3(i)%correct( this%acc(i), this%dt, this%pos3(i) , this%vel(i) )
                end select
                
            end if
            
            
        end do
        
    case(PIM_AB3AM4)
        
        do i = 1, this%numPrtcl
            
            if( this%prtcl_flag(i) >= Pflg_inDomain ) then
                
                select type(this)
                type is (Prtcl_Integration)
                    call this%m_AB3AM4(i)%correct( this%acc(i), this%dt, pos , this%vel(i) )
                    this%pos(i) = real4( pos%x , pos%y, pos%z, this%pos(i)%w) 
                type is(Prtcl_rot_Integration)
                    call this%m_AB3AM4(i)%correct( this%acc(i), this%dt, this%pos3(i) , this%vel(i) )
                end select
                
            end if
            
        end do
        
    case(PIM_AB4AM5)
        
        do i = 1, this%numPrtcl
            
            if( this%prtcl_flag(i) >= Pflg_inDomain ) then
                
                select type(this)
                type is (Prtcl_Integration)
                    call this%m_AB4AM5(i)%correct( this%acc(i), this%dt, pos , this%vel(i) )
                    this%pos(i) = real4( pos%x , pos%y, pos%z, this%pos(i)%w) 
                type is(Prtcl_rot_Integration)
                    call this%m_AB4AM5(i)%correct( this%acc(i), this%dt, this%pos3(i) , this%vel(i) )
                end select
                
            end if
            
        end do
        
    end select
    
end subroutine

!**********************************************************************
! setting new number of particles 
!**********************************************************************
subroutine PI_set_numPrtcl(this, new_numPrtcl )
    implicit none
    class(Prtcl_Integration)   this
    integer(IK),intent(in)  :: new_numPrtcl
    
    !// locals
    integer(IK) i
    type(Gear3_params) G3
    type(Gear4_params) G4
    type(Gear5_params) G5
    type(AB2AM3_params) AM3
    type(AB3AM4_params) AM4
    type(AB4AM5_params) AM5
    
    !//body
    
    select case( this%Int_Method )
        
    case(PIM_AB2AM3)
        
         do i = this%numPrtcl+1,new_numPrtcl
            select type(this)
                type is (Prtcl_Integration)
                    AM3%x  = real3(this%pos(i)%x, this%pos(i)%y, this%pos(i)%z)
                type is (Prtcl_rot_Integration)
                    AM3%x  = this%pos3(i)
            end select
            AM3%v  = this%vel(i)
            AM3%v1 = zero_r3
            AM3%a  = zero_r3
            AM3%a1 = zero_r3
            call this%m_AB2AM3(i)%setval(AM3) 
         end do
    
    case(PIM_AB3AM4)
        
        do i = this%numPrtcl+1,new_numPrtcl
            select type(this)
                type is (Prtcl_Integration)
                    AM4%x  = real3(this%pos(i)%x, this%pos(i)%y, this%pos(i)%z)
                type is (Prtcl_rot_Integration)
                    AM4%x  = this%pos3(i)
            end select
            AM4%v  = this%vel(i)
            AM4%v1 = zero_r3
            AM4%v2 = zero_r3
            AM4%a  = zero_r3
            AM4%a1 = zero_r3
            AM4%a2 = zero_r3
            call this%m_AB3AM4(i)%setval(AM4)
        end do
    
    case(PIM_AB4AM5)
              
        do i = this%numPrtcl+1,new_numPrtcl
            select type(this)
                type is (Prtcl_Integration)
                    AM5%x  = real3(this%pos(i)%x, this%pos(i)%y, this%pos(i)%z)
                type is (Prtcl_rot_Integration)
                    AM5%x  = this%pos3(i)
            end select
            AM5%v  = this%vel(i)
            AM5%v1 = zero_r3
            AM5%v2 = zero_r3
            AM5%v3 = zero_r3
            AM5%a  = zero_r3
            AM5%a1 = zero_r3
            AM5%a2 = zero_r3
            AM5%a3 = zero_r3
            call this%m_AB4AM5(i)%setval(AM5)        
        end do
        
    case(PIM_Gear3)
        
        do i = this%numPrtcl+1,new_numPrtcl
            select type(this)
                type is (Prtcl_Integration)
                    G3%x = real3(this%pos(i)%x, this%pos(i)%y, this%pos(i)%z)
                type is (Prtcl_rot_Integration)
                    G3%x = this%pos3(i)
            end select  
            G3%v = this%vel(i)
            G3%a = zero_r3
            G3%b = zero_r3
            call this%m_Gear3(i)%setval(G3)
        end do

    case(PIM_Gear4)
        
        do i = this%numPrtcl+1,new_numPrtcl
            select type(this)
                type is (Prtcl_Integration)                
                    G4%x = real3(this%pos(i)%x, this%pos(i)%y, this%pos(i)%z)
                type is (Prtcl_rot_Integration)
                    G4%x = this%pos3(i)
            end select                
            G4%v = this%vel(i)
            G4%a = zero_r3
            G4%b = zero_r3
            G4%c = zero_r3
            call this%m_Gear4(i)%setval(G4)
        end do        
        
    case(PIM_Gear5)

        do i = this%numPrtcl+1,new_numPrtcl
            select type(this)
                type is (Prtcl_Integration)
                    G5%x = real3(this%pos(i)%x, this%pos(i)%y, this%pos(i)%z)
                type is  (Prtcl_rot_Integration)
                    G5%x = this%pos3(i)
            end select
            G5%v = this%vel(i)
            G5%a = zero_r3
            G5%b = zero_r3
            G5%c = zero_r3
            G5%d = zero_r3
            call this%m_Gear5(i)%setval(G5)
        end do
            
    end select
    
    this%numPrtcl = new_numPrtcl
    
end subroutine

!**********************************************************************
! setting new/current number of particles
!**********************************************************************
subroutine PRI_set_numPrtcl(this, new_numprtcl )
    implicit none
    class(Prtcl_rot_Integration) this
    integer(IK),intent(in)  :: new_numPrtcl
    
    this%numPrtcl = new_numPrtcl
    
end subroutine 

!**********************************************************************
! destructor/final procedure
!**********************************************************************
subroutine PI_Finalize(this)
    type(prtcl_Integration) this
    
    if(allocated(this%m_ForwardEuler))deallocate(this%m_ForwardEuler)
    if(allocated(this%m_ModifiedEuler))deallocate(this%m_ModifiedEuler)
    if(allocated(this%m_Taylor2)) deallocate(this%m_Taylor2)
    if(allocated(this%m_Taylor3)) deallocate(this%m_Taylor3)
    if(allocated(this%m_Taylor4)) deallocate(this%m_Taylor4)
    if(allocated(this%m_AB2)) deallocate(this%m_AB2)
    if(allocated(this%m_AB3)) deallocate(this%m_AB3)
    if(allocated(this%m_AB4)) deallocate(this%m_AB4)
    if(allocated(this%m_AB5)) deallocate(this%m_AB5)
    if(allocated(this%m_AB2AM3)) deallocate(this%m_AB2AM3)
    if(allocated(this%m_AB3AM4)) deallocate(this%m_AB3AM4)
    if(allocated(this%m_AB4AM5)) deallocate(this%m_AB4AM5)
    if(allocated(this%m_Gear3)) deallocate(this%m_Gear3)
    if(allocated(this%m_Gear4)) deallocate(this%m_Gear4)
    if(allocated(this%m_Gear5)) deallocate(this%m_Gear5)
    
end subroutine 


!**********************************************************************
! allocating memory and sets the initial value of integration parameters
! if necessary
!**********************************************************************
subroutine PI_AllocateMethod(this)
    implicit none
    class(Prtcl_Integration)  this
    
    ! locals
    integer(IK) i, max_nPrtcl
    type(Gear3_params) G3
    type(Gear4_params) G4
    type(Gear5_params) G5
    type(AB2AM3_params) AM3
    type(AB3AM4_params) AM4
    type(AB4AM5_params) AM5
    
    max_nPrtcl = this%max_nPrtcl
    
    select case(this%Int_Method)
    case(PIM_FE)
        
        if(allocated(this%m_ForwardEuler))deallocate(this%m_ForwardEuler)
        allocate(this%m_ForwardEuler)
        
    case(PIM_ME)
        
        if(allocated(this%m_ModifiedEuler))deallocate(this%m_ModifiedEuler)
        allocate(this%m_ModifiedEuler)
        
    case(PIM_Taylor2)
        
        if(allocated(this%m_Taylor2)) deallocate(this%m_Taylor2)
        allocate(this%m_Taylor2)
        
    case(PIM_Taylor3)
        
        if(allocated(this%m_Taylor3)) deallocate(this%m_Taylor3)
        allocate(this%m_Taylor3(max_nPrtcl))
        do i=1, max_nPrtcl
            call this%m_Taylor3(i)%setZero()
        end do
        
    case(PIM_Taylor4)
        
        if(allocated(this%m_Taylor4)) deallocate(this%m_Taylor4)
        allocate(this%m_Taylor4(max_nPrtcl))
        do i=1, max_nPrtcl
            call this%m_Taylor4(i)%setZero()
        end do
    
   case(PIM_AB2)
        
        if(allocated(this%m_AB2)) deallocate(this%m_AB2)
        allocate( this%m_AB2( max_nPrtcl ) )
        do i =1, max_nPrtcl
            call this%m_AB2(i)%setZero()
        end do
        
    case(PIM_AB3)
        
        if(allocated(this%m_AB3)) deallocate(this%m_AB3)
        allocate( this%m_AB3( max_nPrtcl ) )
        do i=1,max_nPrtcl
            call this%m_AB3(i)%setZero()
        end do
    
    case(PIM_AB4)
        if(allocated(this%m_AB4)) deallocate(this%m_AB4)
        allocate( this%m_AB4( max_nPrtcl ) )
        do i=1,max_nPrtcl
            call this%m_AB4(i)%setZero()
        end do
    
    case(PIM_AB5)
        if(allocated(this%m_AB5)) deallocate(this%m_AB5)
        allocate( this%m_AB5( max_nPrtcl ) )
        do i=1,max_nPrtcl
            call this%m_AB5(i)%setZero()
        end do
    
    case(PIM_AB2AM3)
        
        if(allocated(this%m_AB2AM3)) deallocate(this%m_AB2AM3)
        allocate( this%m_AB2AM3( max_nPrtcl ) )
        do i=1,max_nPrtcl
            select type(this)
                type is(Prtcl_Integration)
                    AM3%x  = real3(this%pos(i)%x, this%pos(i)%y, this%pos(i)%z)
                type is(Prtcl_rot_Integration)
                    AM3%x  = this%pos3(i)
            end select
            AM3%v  = this%vel(i)
            AM3%v1 = zero_r3
            AM3%a  = zero_r3
            AM3%a1 = zero_r3
            call this%m_AB2AM3(i)%setval(AM3)
            
        end do
        
    case(PIM_AB3AM4)
        
        if(allocated(this%m_AB3AM4)) deallocate(this%m_AB3AM4)
        allocate( this%m_AB3AM4( max_nPrtcl ) )
        do i=1,max_nPrtcl
            select type(this)
                type is(Prtcl_Integration)
                    AM4%x  = real3(this%pos(i)%x, this%pos(i)%y, this%pos(i)%z)
                type is(Prtcl_rot_Integration)
                    AM4%x  = this%pos3(i)
            end select
                
            AM4%v  = this%vel(i)
            AM4%v1 = zero_r3
            AM4%v2 = zero_r3
            AM4%a  = zero_r3
            AM4%a1 = zero_r3
            AM4%a2 = zero_r3
            call this%m_AB3AM4(i)%setval(AM4)
            
        end do
        
    case(PIM_AB4AM5)
        
        if(allocated(this%m_AB4AM5)) deallocate(this%m_AB4AM5)
        allocate( this%m_AB4AM5( max_nPrtcl ) )
        do i=1,max_nPrtcl
            select type(this)
                type is (Prtcl_Integration)
                    AM5%x  = real3(this%pos(i)%x, this%pos(i)%y, this%pos(i)%z)
                type is (Prtcl_rot_Integration)
                    AM5%x  = this%pos3(i)
            end select
            AM5%v  = this%vel(i)
            AM5%v1 = zero_r3
            AM5%v2 = zero_r3
            AM5%v3 = zero_r3
            AM5%a  = zero_r3
            AM5%a1 = zero_r3
            AM5%a2 = zero_r3
            AM5%a3 = zero_r3
            call this%m_AB4AM5(i)%setval(AM5)
            
        end do
        
    case(PIM_Gear3)
        if(allocated(this%m_Gear3)) deallocate(this%m_Gear3)
        allocate( this%m_Gear3(max_nPrtcl ) )
        do i=1,max_nPrtcl
            select type(this)
                type is (Prtcl_Integration)
                    G3%x = real3(this%pos(i)%x, this%pos(i)%y, this%pos(i)%z)
                type is (Prtcl_rot_Integration)
                    G3%x = this%pos3(i)
            end select  
            G3%v = this%vel(i)
            G3%a = zero_r3
            G3%b = zero_r3
            call this%m_Gear3(i)%setval(G3)
        end do

    case(PIM_Gear4)
        
        if(allocated(this%m_Gear4)) deallocate(this%m_Gear4)
        allocate( this%m_Gear4( max_nPrtcl ) )
        do i=1,max_nPrtcl
            select type(this)
                type is (Prtcl_Integration)                
                    G4%x = real3(this%pos(i)%x, this%pos(i)%y, this%pos(i)%z)
                type is (Prtcl_rot_Integration)
                    G4%x = this%pos3(i)
            end select
                
            G4%v = this%vel(i)
            G4%a = zero_r3
            G4%b = zero_r3
            G4%c = zero_r3
            call this%m_Gear4(i)%setval(G4)
        end do        
        
    case(PIM_Gear5)
        if(allocated(this%m_Gear5)) deallocate(this%m_Gear5)
        allocate( this%m_Gear5( max_nPrtcl ) )
        do i=1,max_nPrtcl
            select type(this)
                type is (Prtcl_Integration)
                    G5%x = real3(this%pos(i)%x, this%pos(i)%y, this%pos(i)%z)
                type is  (Prtcl_rot_Integration)
                    G5%x = this%pos3(i)
            end select
            
            G5%v = this%vel(i)
            G5%a = zero_r3
            G5%b = zero_r3
            G5%c = zero_r3
            G5%d = zero_r3
            call this%m_Gear5(i)%setval(G5)
        end do                
        
        case default
        stop "Invalid integration method in the input"
    end select
    
    
end subroutine

end module
