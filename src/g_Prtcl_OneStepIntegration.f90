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
!  file name  : g_Prtcl_OneStepIntegration.f90
!  module name: g_Prtcl_OneStepIntegration 
! 
!  Purpose:                                
!  1) providing classes for performing one-step integrations they includes:
!     - ForwardEuler 
!     - ModifiedEuler
!     - Taylor2      
!     - Taylor3 
!     - Taylor4 
!   
!   Int_OneStep_Method is the base class for all one step integration methods.
!     Each derived class (child class) should provide the definitions for two
!     methods:     
!       - Integrate
!       - setZero  
! 
!   New one step integrations can be extended based on this abstract class.
! 
!------------------------------------------------------------------------------
    
module g_Prtcl_OneStepIntegration

    use g_TypeDef
    implicit none
    
    private
    public:: Int_OneStep_Method, ForwardEuler, ModifiedEuler, Taylor2, Taylor3, Taylor4
    
    ! abstract class for one step methods
    type,abstract:: Int_OneStep_Method
    
    contains
    
        procedure( AB_Integrate ),deferred  :: Integrate
        procedure( AB_setZero),deferred     :: setZero
    
    end type
    
    ! abstract interface for one-step integrations algorithms 
    abstract interface
        subroutine AB_Integrate( this , acc , dt, pos, vel )
            import Int_OneStep_Method
            import real3, RK
            class(Int_OneStep_Method) this
            type(real3),intent(in)   :: acc
            real(RK),   intent(in)   :: dt
            type(real3),intent(inout):: pos, vel
        end subroutine
    
        subroutine AB_setZero(this)
                import Int_OneStep_Method
                class(Int_OneStep_Method) this
        end subroutine 
        
    end interface
    
    
    ! Forward Euler method 
    type,extends (Int_OneStep_Method) :: ForwardEuler
        
    contains
        procedure:: Integrate   => FE_Integrate
        procedure:: setZero     => FE_setZero
    end type
    
    
    ! Modified Euler method
    type,extends (Int_OneStep_Method) :: ModifiedEuler
        
    contains
        procedure:: Integrate   => BE_Integrate
        procedure:: setZero     => BE_setZero
    end type
        
    !Taylor 2nd order 
    type, extends (Int_OneStep_Method):: Taylor2
    
    contains 
        procedure:: Integrate   => T2_Integrate
        procedure:: setZero     => T2_setZero
    end type
    
    ! Taylor 3rd order
    type, extends(Taylor2):: Taylor3
        type(real3) a0 ! to be used for taking the first derivative of acceleration
    contains
        procedure:: Integrate   => T3_Integrate
    end type
    
    ! Taylor 4th order
    type, extends(Taylor3):: Taylor4
        type(real3) b0 ! to be used for taking the second derivative of acceleration
    contains
        procedure:: Integrate   => T4_Integrate
    end type
    
contains

subroutine FE_setZero(this)
    implicit none
    class(ForwardEuler) this
end subroutine 

subroutine BE_setZero(this)
    implicit none
    class(ModifiedEuler) this
end subroutine 

subroutine T2_setZero(this)
    implicit none
    class(Taylor2) this
    
    select type(this)
    type is (Taylor3)
        this%a0 = zero_r3
    type is (Taylor4)
        this%a0 = zero_r3
        this%b0 = zero_r3
    end select
    
end subroutine 

!**********************************************************************
! Forward Euler method
!**********************************************************************
subroutine FE_Integrate(this , acc , dt , pos, vel )
    implicit none
    class(ForwardEuler) this
    type(real3),intent(in)   :: acc   ! acceleration
    real(RK),   intent(in)   :: dt    ! time step 
    type(real3),intent(inout):: pos, vel
    
    pos = pos + vel * dt ! position equation
    vel = vel + acc * dt ! velocity equation
    
end subroutine

!**********************************************************************
! Backward Euler or modified Euler method
!**********************************************************************
subroutine BE_Integrate(this , acc , dt , pos, vel )
    implicit none
    class(ModifiedEuler) this
    type(real3),intent(in)   :: acc  ! acceleration
    real(RK),   intent(in)   :: dt   ! time step 
    type(real3),intent(inout):: pos, vel
    
    vel = vel + acc * dt ! velocity equation
    pos = pos + vel * dt ! position equation
    
end subroutine

!**********************************************************************
! Taylor 2nd order
!**********************************************************************
subroutine T2_Integrate(this , acc , dt , pos, vel )
    implicit none
    class(Taylor2) this
    type(real3),intent(in)   :: acc
    real(RK),   intent(in)   :: dt
    type(real3),intent(inout):: pos, vel
    
    pos = pos + vel * dt + 0.5_RK * acc * dt * dt
    vel = vel + acc * dt
    
end subroutine

!**********************************************************************
! Taylor 3rd order
!**********************************************************************
subroutine T3_Integrate(this , acc , dt , pos, vel )
    implicit none
    class(Taylor3) this
    type(real3),intent(in)   :: acc
    real(RK),   intent(in)   :: dt
    type(real3),intent(inout):: pos, vel
    
    !// locals
    type(real3) b
    
    !// body
    b = (acc-this%a0) !b = (acc-this%a0)/dt
    pos = pos + vel * dt + 0.5_RK * acc* dt*dt + 0.1666667_RK*b*dt*dt
    vel = vel + acc * dt + 0.5_RK * b * dt
    this%a0 = acc
    
end subroutine

!**********************************************************************
! Taylor 4th order
!**********************************************************************
subroutine T4_Integrate(this , acc , dt , pos, vel )
    implicit none
    class(Taylor4) this
    type(real3),intent(in)   :: acc
    real(RK),   intent(in)   :: dt
    type(real3),intent(inout):: pos, vel
    
    !// locals
    type(real3) b, c
    
    !// body
    b = (acc-this%a0) !b = (acc-this%a0)/dt
    c = (b - this%b0) !c = (b - this%b0)/dt
    pos = pos + vel * dt + 0.5_RK * acc* dt*dt + 0.1666667_RK*b*dt*dt + 0.04166666667_RK*c*dt**3.0_RK
    vel = vel + acc * dt + 0.5_RK * b * dt + 0.1666667_RK*c*dt*dt
    this%a0 = acc
    this%b0 = b
    
end subroutine


  
end module
