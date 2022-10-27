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
!  file name  : g_Prtcl_TwoStepIntegration.f90
!  module name: g_Prtcl_TwoStepIntegration
! 
!  Purpose:
!   1) providing classes for performing two step integrationv(predictor-
!      corrector) methods, they are:
!        - Gear3
!        - Gear4
!        - Gear5
!        - AB2AM3
!        - AB3AM4
!        - AB4AM5
!  
!  Int_TwoStep_Method class is a base class provides interface for all two-step
!    integration methods. Every new class which is derived from this class
!    should provide definitions for three methods:
!      - setVal
!      - Predict
!      - Correct
! 
!------------------------------------------------------------------------------
    
module g_Prtcl_TwoStepIntegration
    
    use g_TypeDef
    use g_error_handling
    use g_Prtcl_MultiPointIntegration
    implicit none
    
    private
    public:: Int_TwoStep_Method,  Gear3_params, Gear4_params, Gear5_params, AB2AM3_params, AB3AM4_params, AB4AM5_params
    
    real(RK),parameter,dimension(4):: G3_coeff = (/ 1.0_RK/12.0_RK, 5.0_RK/12.0_RK, 1.0_RK, 1.0_RK /)
    real(RK),parameter,dimension(5):: G4_coeff = (/ 19.0_RK/240.0_RK, 0.375_RK, 1.0_RK, 1.5_RK, 1.0_RK /)
    real(RK),parameter,dimension(6):: G5_coeff = (/ 3.0_RK/32.0_RK, 251.0_RK/720.0_RK, 1.0_RK, 11.0_RK/6.0_RK, 2.0_RK, 1.0_RK /)
    
    real(RK),parameter,dimension(3):: AM3_coeff = (/ 5.0_RK/12.0_RK, 8.0_RK/12.0_RK, 1.0_RK/12.0_RK /)
    real(RK),parameter,dimension(4):: AM4_coeff = (/ 9.0_RK/24.0_RK, 19.0_RK/24.0_RK, 5.0_RK/24.0_RK, 1.0_RK/24.0_RK /)
    real(RK),parameter,dimension(5):: AM5_coeff = (/ 251.0_RK/720.0_RK, 646.0_RK/720.0_RK, &
                                                     264.0_RK/720.0_RK, 106.0_RK/720.0_RK , 19.0_RK/720.0_RK /)
    
    
    type, abstract::  Int_TwoStep_Method
        
    contains
        procedure(Abs_predict),deferred:: predict
        procedure(Abs_correct),deferred:: correct
        procedure(Abs_setVal) ,deferred:: setVal
    end type
    
    
    abstract interface
        subroutine Abs_predict(this, dt, pos, vel )
        
            import Int_TwoStep_Method, RK, real3
            
            class(Int_TwoStep_Method) this
            real(RK),   intent(in)  :: dt 
            type(real3),intent(out) :: pos, vel
    
        end subroutine 
        
        subroutine Abs_correct(this, acc , dt, pos ,vel )
    
            import Int_TwoStep_Method, RK, real3
            class(Int_TwoStep_Method) this
            type(real3),intent(in)    :: acc
            real(RK),   intent(in)    :: dt
            type(real3),intent(inout) :: pos, vel ! when invoking, they contain predicted values calculated by Abs_predict
                                                  ! in return procedure replaces them with corrected (final) values 
            
        end subroutine   
        
        subroutine Abs_setVal(this, val)
            import Int_TwoStep_Method
            class(Int_TwoStep_Method) this
            class(*) val
        end subroutine 
        
    end interface
    
        
    !// Gear3 method
    type, extends(Int_TwoStep_Method):: Gear3_params
        type(real3) :: x = zero_r3 ! r
        type(real3) :: v = zero_r3 ! v    
        type(real3) :: a = zero_r3 ! a
        type(real3) :: b = zero_r3 ! b
    
    contains
    
    procedure:: predict => G3_predict
    procedure:: correct => G3_correct
    procedure:: setVal  => G3_setVal    
    
    end type!// Gear3 method
    
    
    !// Gear4 method
    type,extends(Gear3_params) :: Gear4_params
        type(real3) :: c = zero_r3 

    contains
        
    procedure:: predict => G4_predict
    procedure:: correct => G4_correct
    procedure:: setVal  => G4_setVal
    
    end type !// Gear4 method
    
    
    !// Gear5 method
    type,extends(Gear4_params) :: Gear5_params
        type(real3) :: d = zero_r3
    
    contains
    
    procedure:: predict => G5_predict
    procedure:: correct => G5_correct
    procedure:: setVal  => G5_setVal
    
    end type !// Gear5 method
    
    !// AB2AM3 method
    type,extends(Int_TwoStep_Method):: AB2AM3_params
        type(real3):: x
        type(real3):: v  = zero_r3
        type(real3):: v1 = zero_r3
        type(real3):: a  = zero_r3
        type(real3):: a1 = zero_r3
        
    contains
    
    procedure:: predict => AB2AM3_predict
    procedure:: correct => AB2AM3_correct
    procedure:: setVal  => AB2AM3_setVal
    
    endtype !// AB2AM3 method
    
    ! AB3AM4 method
    type,extends(AB2AM3_params):: AB3AM4_params
        type(real3):: v2 = zero_r3
        type(real3):: a2 = zero_r3
        
    contains
    
    procedure:: predict => AB3AM4_predict
    procedure:: correct => AB3AM4_correct
    procedure:: setVal  => AB3AM4_setVal
    
    endtype !// AB3AM4 method
    
    
    !// AB4AM5 method
    type,extends(AB3AM4_params):: AB4AM5_params
        type(real3):: v3 = zero_r3
        type(real3):: a3 = zero_r3
        
    contains
    
    procedure:: predict => AB4AM5_predict
    procedure:: correct => AB4AM5_correct
    procedure:: setVal  => AB4AM5_setVal
    
    endtype !// AB4Am5 method
    
    
contains

!**********************************************************************
! AB2AM3 prediction step
!**********************************************************************
subroutine AB2AM3_predict(this, dt, pos, vel )
    implicit none
    class( AB2AM3_params) this
    real(RK),   intent(in)  :: dt
    type(real3),intent(out) :: pos, vel
     
    ! predicts the position and velocity and returns them
    pos = this%x + ( AB2_coeff(1) * this%v - AB2_coeff(2) * this%v1 ) * dt
    vel = this%v + ( AB2_coeff(1) * this%a - AB2_coeff(2) * this%a1 ) * dt  
    
end subroutine

!**********************************************************************
! AB2AM3 correction step
!**********************************************************************
subroutine AB2AM3_correct(this, acc , dt, pos ,vel )
    implicit none
    class(AB2AM3_params) this
    type(real3),intent(in)  :: acc
    real(RK),   intent(in)  :: dt 
    type(real3),intent(inout) :: pos, vel
    
    !// locals
    type(real3) c_pos, c_vel
        
    !// body
    
    ! pos, vel and acc contain predicted values of position, velocity and acceleration
    
    ! correcting position, velocity 
    c_pos  = this%x + ( AM3_coeff(1)*vel + AM3_coeff(2)*this%v - AM3_coeff(3)*this%v1 ) * dt
    c_vel  = this%v + ( AM3_coeff(1)*acc + AM3_coeff(2)*this%a - AM3_coeff(3)*this%a1 ) * dt
    
    ! swapping variables to be used in next iteration
    this%v1 = this%v
    this%a1 = this%a
    this%x = c_pos
    this%v = c_vel
    this%a = acc
    ! returning corrected values     
    pos = c_pos
    vel = c_vel

end subroutine

!**********************************************************************
! AB2AM3 setVal method
!**********************************************************************
subroutine AB2AM3_setVal(this, val)
    implicit none
    class(AB2AM3_params) this
    class(*) val
    
    select type( val )
    type is( AB2AM3_params)
        this%x = val%x
        this%v = val%v
        this%v1= val%v1
        this%a = val%a
        this%a1= val%a1
        class default
            call CheckForError( ErrT_Abort, "AB2AM3_setVal" , "Invalid val in argument" )
    end select
    
    
end subroutine

!**********************************************************************
! AB3AM4 prediction step
!**********************************************************************
subroutine AB3AM4_predict(this, dt, pos, vel )
    
    implicit none
    class( AB3AM4_params) this
    real(RK),   intent(in)  :: dt
    type(real3),intent(out) :: pos, vel
    
    !//body
    
    ! predicting the position and velocity and returns them
    pos = this%x + ( AB3_coeff(1) * this%v - AB3_coeff(2) * this%v1 + AB3_coeff(3) * this%v2 ) * dt
    vel = this%v + ( AB3_coeff(1) * this%a - AB3_coeff(2) * this%a1 + AB3_coeff(3) * this%a2 ) * dt  
    
end subroutine

!**********************************************************************
! AB3AM4 correction step
!**********************************************************************
subroutine AB3AM4_correct(this, acc , dt, pos ,vel )
    
    implicit none
    class(AB3AM4_params) this
    type(real3),intent(in)  :: acc
    real(RK),   intent(in)  :: dt 
    type(real3),intent(inout) :: pos, vel
    
    !// locals
    type(real3) c_pos, c_vel
     
    !// body
    
    ! pos, vel and acc contain predicted values of position, velocity and acceleration
    ! corrects position and velocity 
    c_pos  = this%x + ( AM4_coeff(1)*vel + AM4_coeff(2)*this%v - AM4_coeff(3)*this%v1 + AM4_coeff(4)*this%v2 ) * dt
    c_vel  = this%v + ( AM4_coeff(1)*acc + AM4_coeff(2)*this%a - AM4_coeff(3)*this%a1 + AM4_coeff(4)*this%a2 ) * dt
    
    ! shifting variables to be used in the next iteration
    this%x = c_pos
    
    this%v2 = this%v1
    this%v1 = this%v
    this%v = c_vel
    
    this%a2 = this%a1
    this%a1 = this%a
    this%a = acc 
    
    ! returning corrected values     
    pos = c_pos
    vel = c_vel

end subroutine

!**********************************************************************
! AB3Am4 setVal method
!**********************************************************************
subroutine AB3AM4_setVal(this, val)
    implicit none
    class(AB3AM4_params) this
    class(*) val
    
    !// body    
    select type( val )
    type is( AB3AM4_params)
        this%x = val%x
        this%v = val%v
        this%v1= val%v1
        this%v2= val%v2
        this%a = val%a
        this%a1= val%a1
        this%a2= val%a2
        class default
            call CheckForError( ErrT_Abort, "AB3AM4_setVal" , "Invalid val in argument" )
    end select
    
end subroutine


!**********************************************************************
! AB4AM5 prediction step
!**********************************************************************
subroutine AB4AM5_predict(this, dt, pos, vel )
    implicit none
    class( AB4AM5_params) this
    real(RK),   intent(in)  :: dt
    type(real3),intent(out) :: pos, vel
    
    !// body
    
    ! predicting the position and velocity and returning them
    pos = this%x + ( AB4_coeff(1) * this%v - AB4_coeff(2) * this%v1 + &
                     AB4_coeff(3) * this%v2 - AB4_coeff(4) * this%v3) * dt
    vel = this%v + ( AB4_coeff(1) * this%a - AB4_coeff(2) * this%a1 + &
                     AB4_coeff(3) * this%a2 - AB4_coeff(4) * this%a3) * dt  
    
end subroutine


!**********************************************************************
! AB4AM5 correction step
!**********************************************************************
subroutine AB4AM5_correct(this, acc , dt, pos ,vel )
    implicit none
    class(AB4AM5_params) this
    type(real3),intent(in)  :: acc
    real(RK),   intent(in)  :: dt 
    type(real3),intent(inout) :: pos, vel
    
    !// locals
    type(real3) c_pos, c_vel
    
    !// body 
    
    ! pos, vel and acc contain predicted values of position, velocity and acceleration
    
    ! corrects position, velocity 
    c_pos  = this%x + ( AM5_coeff(1)*vel + AM5_coeff(2)*this%v - AM5_coeff(3)*this%v1 + &
                        AM5_coeff(4)*this%v2 - AM5_coeff(5)*this%v3) * dt
    c_vel  = this%v + ( AM5_coeff(1)*acc + AM5_coeff(2)*this%a - AM5_coeff(3)*this%a1 + &
                        AM5_coeff(4)*this%a2 - AM5_coeff(5)*this%a3) * dt
    
    ! swapping variables to be used in the next iteration
    this%x = c_pos
    this%v3 = this%v2
    this%v2 = this%v1
    this%v1 = this%v
    this%v = c_vel
    
    this%a3 = this%a2
    this%a2 = this%a1
    this%a1 = this%a
    this%a = acc 
    
    ! returning corrected values     
    pos = c_pos
    vel = c_vel

end subroutine

!**********************************************************************
! AB4Am4 setVal method
!**********************************************************************
subroutine AB4AM5_setVal(this, val)
    implicit none
    class(AB4AM5_params) this
    class(*) val
    
    select type( val )
    type is( AB4AM5_params)
        this%x = val%x
        this%v = val%v
        this%v1= val%v1
        this%v2= val%v2
        this%v3= val%v3
        this%a = val%a
        this%a1= val%a1
        this%a2= val%a2
        this%a3= val%a3
        class default
            call CheckForError( ErrT_Abort, "AB4AM5_setVal" , "Invalid val in argument" )
    end select
    
    
end subroutine

!**********************************************************************
! Gear3 prediction step
!**********************************************************************
subroutine G3_predict(this, dt, pos, vel )
    implicit none
    class( Gear3_params) this
    real(RK),   intent(in)  :: dt
    type(real3),intent(out) :: pos, vel
    
    !// locals
    real(RK) dt2
    
    dt2 = dt*dt
    
    ! predicting the position and velocity and acceleration and keeping them 
    this%x = this%x + this%v * dt +  0.5_RK * dt2 * this%a + 0.1666666666666667_RK * dt2*dt * this%b
    this%v = this%v + this%a * dt +  0.5_RK * dt2 * this%b 
    this%a = this%a + this%b * dt
    
    ! returns predicted position and velocity
    pos = this%x
    vel = this%v

end subroutine


!**********************************************************************
! Gear3 correction method
!**********************************************************************
subroutine G3_correct(this, acc , dt, pos ,vel )
    implicit none
    class(Gear3_params) this
    type(real3),intent(in)  :: acc
    real(RK),   intent(in)  :: dt 
    type(real3),intent(inout) :: pos, vel
    
    !// locals
    type(real3) delta_a
    
    !// body
    
    ! evaluating the deviations from real value
    delta_a = acc - this%a
    
    ! correcting position, velocity and acceleration
    this%x = this%x + G3_coeff(1) * dt*dt * delta_a
    this%v = this%v + G3_coeff(2) * dt    * delta_a
    this%a = this%a + G3_coeff(3)         * delta_a
    this%b = this%b + G3_coeff(4) / dt    * delta_a
    
    ! returning corrected values     
    pos = this%x
    vel = this%v

end subroutine

!**********************************************************************
! Gear4 prediction step
!**********************************************************************
subroutine G4_predict(this, dt,  pos, vel )
    implicit none
    class( Gear4_params) this
    real(RK),   intent(in)  :: dt
    type(real3),intent(out) :: pos, vel
    
    !// locals
    real(RK) dt2, dt3
    
    !// body
    dt2 = dt*dt
    dt3 = dt*dt2
    
    ! predicting the position and velocity and acceleration and higher derivatives and keeping them 
    this%x = this%x + this%v * dt +  0.5_RK * dt2 * this%a + & 
                    0.1666666666666667_RK * dt3 * this%b + 0.04166666666666667_RK * dt3*dt *this%c
    this%v = this%v + this%a * dt +  0.5_RK * dt2 * this%b + 0.1666666666666667_RK * dt3 * this%c
    this%a = this%a + this%b * dt +  0.5_RK * dt2 * this%c
    this%b = this%b + this%c * dt
    
     ! returning predicted position and velocity
    pos = this%x
    vel = this%v

end subroutine

!**********************************************************************
! Gear4 correction step
!**********************************************************************
subroutine G4_correct(this, acc , dt, pos ,vel )
    implicit none
    class(Gear4_params) this
    type(real3),intent(in)  :: acc
    real(RK),   intent(in)  :: dt
    type(real3),intent(inout) :: pos, vel
    
    !// locals
    real(RK) dt2
    type(real3) delta_a
    
    !// body
    
    dt2 = dt*dt
    
    ! evaluating the deviations from real value
    delta_a = acc - this%a
    
    ! correcting the position and velocity and acceleration and higher derivatives and keeping them 
    this%x = this%x + G4_coeff(1) * dt2   * delta_a
    this%v = this%v + G4_coeff(2) * dt    * delta_a
    this%a = this%a + G4_coeff(3)         * delta_a
    this%b = this%b + G4_coeff(4) / dt    * delta_a
    this%c = this%c + G4_coeff(5) / dt2   * delta_a

    ! returning corrected position and velocity  
    pos = this%x
    vel = this%v

end subroutine

!**********************************************************************
! Gear5 prediction step
!**********************************************************************
subroutine G5_predict(this, dt, pos, vel )
    implicit none
    class( Gear5_params) this
    real(RK),   intent(in)  :: dt
    type(real3),intent(out) :: pos, vel
    
    !// locals
    real(RK) dt2, dt3
    
    !// body
    dt2 = dt*dt
    dt3 = dt2*dt
    
    ! predicting the position and velocity and acceleration and higher derivatives and keeping them 
    this%x = this%x + this%v * dt +  0.5_RK * dt2 * this%a + & 
                    0.1666666666666667_RK * dt3 * this%b + 0.04166666666666667_RK * dt3*dt *this%c  + &
                    0.00833333333333333333_RK* dt3*dt2* this%d
    
    this%v = this%v + this%a * dt +  0.5_RK * dt2 * this%b + 0.1666666666666667_RK * dt3 * this%c + &
                    0.04166666666666667_RK * dt2*dt2 * this%d
    this%a = this%a + this%b * dt + 0.5_RK * dt2 * this%c + 0.1666666666666667_RK * dt3 * this%d
    this%b = this%b + this%c * dt + 0.5_RK * dt2 * this%d
    this%c = this%c + this%d * dt 
    
    ! returning predicted position and velocity
    pos = this%x
    vel = this%v

end subroutine

!**********************************************************************
! Gear5 correction step
!**********************************************************************
subroutine G5_correct(this, acc , dt, pos ,vel )
    implicit none
    class(Gear5_params) this
    type(real3),intent(in)  :: acc
    real(RK),   intent(in)  :: dt
    type(real3),intent(inout) :: pos, vel
    
    !// locals
    real(RK) dt2
    type(real3) delta_a
    
    !// body
    
    dt2 = dt*dt
    
    ! evaluating the deviations from real value
    delta_a = acc - this%a
    
    ! correcting the position and velocity and acceleration and higher derivatives and keeping them 
    this%x = this%x + G5_coeff(1) * dt2   * delta_a
    this%v = this%v + G5_coeff(2) * dt    * delta_a
    this%a = this%a + G5_coeff(3)         * delta_a
    this%b = this%b + G5_coeff(4) / dt    * delta_a
    this%c = this%c + G5_coeff(5) / dt2   * delta_a
    this%d = this%d + G5_coeff(6)/(dt2*dt)* delta_a 
    
    ! returning corrected position and velocity    
    pos = this%x
    vel = this%v

end subroutine
 
!**********************************************************************
!setVal for Gear3
!**********************************************************************
subroutine G3_setVal(this, val)
    implicit none
    class(Gear3_params) this
    class(*) val
    
    
    select type( val )
    type is( Gear3_params)
        this%x = val%x
        this%v = val%v
        this%a = val%a
        this%b = val%b
        class default
            call CheckForError( ErrT_Abort, "G3_setVal" , "Invalid val in argument" )
    end select
    
end subroutine

!**********************************************************************
! setVal for Gear4
!**********************************************************************
subroutine G4_setVal(this, val)
    implicit none
    class(Gear4_params) this
    class(*) val
    
    select type( val )
    type is( Gear4_params)
        this%x = val%x
        this%v = val%v
        this%a = val%a
        this%b = val%b
        this%c = val%c
        class default
           call CheckForError( ErrT_Abort, "G4_setVal" , "Invalid val in argument" )
    end select
    
end subroutine

!**********************************************************************
! setVal for Gear5
!**********************************************************************
subroutine G5_setVal(this, val)
    implicit none
    class(Gear5_params) this
    class(*) val
    
    select type( val )
    type is( Gear5_params)
        this%x = val%x
        this%v = val%v
        this%a = val%a
        this%b = val%b
        this%c = val%c
        this%d = val%d
        class default
            call CheckForError( ErrT_Abort, "G5_setVal" , "Invalid val in argument" )
    end select
    
end subroutine


end module
