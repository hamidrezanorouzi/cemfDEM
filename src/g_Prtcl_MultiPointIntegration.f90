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
!  file name  : g_Prtcl_MultiPointIntegration.f90 
!  module name: g_Prtcl_MultiPointIntegration
!  
!  Purpose: 
!   1) providing classes for performing one-step multi-point integrations. They
!      include:
!      - AB2
!      - AB3
!      - AB4
!      - AB5
!                                                                   
!       Multi-step methods also are extended based on the one-step  
!   integration method. So, they should provide definitions for     
!   two methods:                                                    
!       - setZero                                                   
!       - Integerate                                                
!                                                                   
!------------------------------------------------------------------------------

 module g_Prtcl_MultiPointIntegration
    
    use g_TypeDef
    use g_error_handling
    use g_Prtcl_OneStepIntegration 
    implicit none
    
    private 
    public:: AB2_params, AB3_params, AB4_params, AB5_params, AB2_coeff, AB3_coeff, AB4_coeff, AB5_coeff 
    
    real(RK),parameter,dimension(2):: AB2_coeff = (/ 3.0_RK/2.0_RK, 1.0_RK/2.0_RK /)
    real(RK),parameter,dimension(3):: AB3_coeff = (/ 23.0_RK/12.0_RK, 16.0_RK/12.0_RK, 5.0_RK/12.0_RK /)
    real(RK),parameter,dimension(4):: AB4_coeff = (/ 55.0_RK/24.0_RK, 59.0_RK/24.0_RK, 37.0_RK/24.0_RK, 9.0_RK/24.0_RK /)
    real(RK),parameter,dimension(5):: AB5_coeff = (/ 1901.0_RK/720.0_RK, 2774.0_RK/720.0_RK, 2616.0_RK/720.0_RK, &
                                                     1274.0_RK/720.0_RK, 251.0_RK/720.0_RK /)
    
    ! AB2 method
    type, extends(Int_OneStep_Method):: AB2_params
        type(real3):: v1= zero_r3 
        type(real3):: a1= zero_r3
    contains
        procedure:: Integrate    => AB2_Integrate
        procedure:: setZero      => AB2_setZero    
    end type
    
    ! AB3 method 
    type, extends(AB2_params):: AB3_params
        type(real3):: v2= zero_r3
        type(real3):: a2= zero_r3
    contains
        procedure:: Integrate    => AB3_Integrate
        
    end type
    
    ! AB4 method
    type, extends(AB3_params):: AB4_params
        type(real3):: v3 = zero_r3
        type(real3):: a3 = zero_r3
    contains
        procedure:: Integrate    => AB4_Integrate
    end type
    
    ! AB5 method
    type,extends(AB4_params):: AB5_params
        type(real3):: v4 = zero_r3
        type(real3):: a4 = zero_r3
    contains
        procedure::Integrate    => AB5_Integrate
    end type
    
contains

subroutine AB2_setZero( this )
    implicit none
    class(AB2_params) this
    
    !// body
    this%v1 = zero_r3
    this%a1 = zero_r3
    
    select type (this)
    type is( AB2_params)
    type is( AB3_params)
        this%v2 = zero_r3
        this%a2 = zero_r3
    type is (AB4_params)
        this%v3 = zero_r3
        this%a3 = zero_r3
    type is(AB5_params)
        this%v3 = zero_r3
        this%a3 = zero_r3
        this%v4 = zero_r3
        this%a4 = zero_r3
        class default
            call CheckForError(ErrT_Abort, "AB2_setZero" , "Wrong class type in the input")
    end select
    
end subroutine 

!**********************************************************************
! AB2 method
!**********************************************************************
subroutine AB2_Integrate(this , acc , dt , pos, vel )
    implicit none
    class(AB2_params) this
    type(real3),intent(in)   :: acc
    real(RK),   intent(in)   :: dt
    type(real3),intent(inout):: pos, vel
        
    type(real3) vn
    
    vn = vel
    
    pos = pos + ( AB2_coeff(1) * vn  - AB2_coeff(2) * this%v1 ) * dt
    vel = vel + ( AB2_coeff(1) * acc - AB2_coeff(2) * this%a1 ) * dt        
    
    ! swaping for the next time step
    
    this%v1 = vn
    this%a1 = acc
    
end subroutine

!**********************************************************************
! AB3 method
!**********************************************************************
subroutine AB3_Integrate(this , acc , dt , pos, vel )
    implicit none
    class(AB3_params) this
    type(real3),intent(in)   :: acc
    real(RK),   intent(in)   :: dt
    type(real3),intent(inout):: pos, vel
    
    ! locals
    type(real3) vn
    
    !// body
    vn = vel
    
    ! shifting vectors one time step back 
    pos = pos + ( AB3_coeff(1) * vn  - AB3_coeff(2) * this%v1 + AB3_coeff(3) * this%v2 ) * dt
    vel = vel + ( AB3_coeff(1) * acc - AB3_coeff(2) * this%a1 + AB3_coeff(3) * this%a2 ) * dt        
    
    this%v2 = this%v1
    this%v1 = vn
    
    this%a2 = this%a1
    this%a1 = acc
    
end subroutine


!**********************************************************************
! AB4 method
!**********************************************************************
subroutine AB4_Integrate(this , acc , dt , pos, vel )
    implicit none
    class(AB4_params) this
    type(real3),intent(in)   :: acc
    real(RK),   intent(in)   :: dt
    type(real3),intent(inout):: pos, vel
    
    !// locals
    type(real3) vn
    
    !// body
    vn = vel
    
    pos = pos + dt * ( AB4_coeff(1)*vn   - AB4_coeff(2)*this%v1 + AB4_coeff(3)*this%v2 - AB4_coeff(4)*this%v3 )
    vel = vel + dt * ( AB4_coeff(1)*acc  - AB4_coeff(2)*this%a1 + AB4_coeff(3)*this%a2 - AB4_coeff(4)*this%a3 )
    
    ! shifting vectors one time step back 
    this%v3 = this%v2
    this%v2 = this%v1
    this%v1 = vn
  
    this%a3 = this%a2
    this%a2 = this%a1
    this%a1 = acc 
    
end subroutine


!**********************************************************************
! AB5 method
!**********************************************************************
subroutine AB5_Integrate(this , acc , dt , pos, vel )
    implicit none
    class(AB5_params) this
    type(real3),intent(in)   :: acc
    real(RK),   intent(in)   :: dt
    type(real3),intent(inout):: pos, vel
    
    !// locals
    type(real3) vn
    
    !// body
    
    vn = vel
    pos = pos + dt * ( AB5_coeff(1)*vn  -     &
                       AB5_coeff(2)*this%v1 + &
                       AB5_coeff(3)*this%v2 - &
                       AB5_coeff(4)*this%v3 + &
                       AB5_coeff(5)*this%v4 )

    vel = vel + dt * ( AB5_coeff(1)*acc -     &
                       AB5_coeff(2)*this%a1 + &
                       AB5_coeff(3)*this%a2 - &
                       AB5_coeff(4)*this%a3 + &
                       AB5_coeff(5)*this%a4 )
    
    this%v4 = this%v3
    this%v3 = this%v2
    this%v2 = this%v1
    this%v1 = vn
      
    this%a4 = this%a3
    this%a3 = this%a2
    this%a2 = this%a1
    this%a1 = acc 
    
    
end subroutine 
    
end module
