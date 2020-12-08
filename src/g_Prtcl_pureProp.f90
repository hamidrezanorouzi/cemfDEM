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
!  file name  : g_Prtcl_PureProperty.f90      
!  module name: g_Prtcl_PureProperty          
! 
!  Purpose:                                   
!    Object for physical properties of a  particle or wall
! 
!------------------------------------------------------------------------------

module g_Prtcl_PureProperty

    use g_Prtcl_DefaultValues
    implicit none
    
    
type, abstract:: base_PureProperty
    real(RK)  Rad
    real(RK)  Density
    real(RK)  Vol
    real(RK)  mass
    real(RK)  inertia
    real(RK)  Young_mod
    real(RK)  Shear_mod
    real(RK):: YeildStress      = d_YeildStress
    real(RK):: poissons_ratio   = d_Poissons_ratio
    real(RK):: kappa            = d_kappa
    
contains
    
    procedure(base_PP_clc_prop),deferred:: clc_prop    
    
    procedure:: set_prop    => BPP_set_prop
    procedure:: set_rad     => BPP_set_rad
    procedure:: set_density => BPP_set_density
    procedure:: set_Young   => BPP_set_Young
    procedure:: set_Shear   => BPP_set_Shear
    procedure:: set_Pratio  => BPP_set_Pratio
    procedure:: set_Yeild   => BPP_set_Yeild
    
end type


abstract interface
    subroutine base_PP_clc_prop(this) 
        import base_PureProperty
        class(base_PureProperty) this 
    end subroutine 
end interface



type,extends(base_PureProperty) :: PureProperty
    
    contains
    procedure:: clc_prop    => PP_clc_prop
    
end type

    
    interface PureProperty
        procedure:: PureProperty_const
    end interface

type BinaryProperty
    real(RK) Yeff   ! effective Young's modulus
    real(RK) Geff   ! effective shear modulus
    real(RK) Reff   ! effective radius
    real(RK) meff   ! effective mass
    real(RK):: en   = d_en
    real(RK):: et   = d_et
    real(RK) Kn
    real(RK) Kt
    real(RK):: fric = d_fric
    real(RK):: roll_fric = d_roll_fric
    real(RK):: dmp_n = 0.0
    real(RK):: dmp_t = 0.0
end type

  
contains

subroutine PP_clc_prop(this)
    implicit none
    class(PureProperty) this
    
    this%Vol     = 4.0_RK /3.0_RK * Pi * (this%Rad)**3.0_RK
    this%mass    = this%Density * this%Vol
    this%inertia = 2.0_RK/5.0_RK*this%mass*(this%Rad**2.0_RK)
    this%kappa = (1.0_RK-this%poissons_ratio) / (1.0_RK- 0.5_RK*this%poissons_ratio )
        
end subroutine 


function PureProperty_const( Density , Young_mod , Shear_mod , poisson , YeildStress ) result(this)
    implicit none
    real(RK),intent(in) :: Density , Young_mod , Shear_mod , poisson, YeildStress
    type(PureProperty) this
    
   
    
    call this%set_prop( Density , Young_mod , Shear_mod , poisson , YeildStress )
            
end function

subroutine  BPP_set_prop (this, Density , Young_mod , Shear_mod , poisson , YeildStress )
    implicit none
    class(base_PureProperty) this
    real(RK),intent(in) :: Density , Young_mod , Shear_mod , poisson, YeildStress
    
    this%Density   = Density
    this%Young_mod = Young_mod
    this%Shear_mod = Shear_mod
    this%poissons_ratio = poisson
    this%YeildStress = YeildStress
    
    call this%clc_prop()
    
end subroutine

subroutine BPP_set_rad(this , v )
    implicit none
    class(base_PureProperty) this
    real(RK),intent(in):: v
    
    this%Rad = v
end subroutine

subroutine BPP_set_Density(this , v )
    implicit none
    class(base_PureProperty) this
    real(RK),intent(in):: v
    
    this%Density = v
end subroutine

subroutine BPP_set_Young(this , v )
    implicit none
    class(base_PureProperty) this
    real(RK),intent(in):: v
    
    this%Young_mod = v
end subroutine

subroutine BPP_set_Shear(this , v )
    implicit none
    class(base_PureProperty) this
    real(RK),intent(in):: v
    
    this%Shear_mod = v
end subroutine

subroutine BPP_set_Pratio(this , v )
    implicit none
    class(base_PureProperty) this
    real(RK),intent(in):: v
    
    this%poissons_ratio = v
end subroutine


subroutine BPP_set_Yeild(this , v )
    implicit none
    class(base_PureProperty) this
    real(RK),intent(in):: v
    
    this%YeildStress = v
end subroutine


end module
    
