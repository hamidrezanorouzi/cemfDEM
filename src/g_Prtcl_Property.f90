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

!  file name  : g_Prtcl_Properties.f90
!  module name: g_Prtcl_Properties
!  
!  Purpose:
!    Stores for physical properties of particles and walls 
! 
!------------------------------------------------------------------------------  
    
module g_Prtcl_Properties
    
    use g_Prtcl_PureProperty
    use g_MakePrtcls
 
    implicit none
    
    
    type, abstract:: base_PhysicalProperty
    
        integer(IK):: numPrtcl_type = 0
        integer(IK):: numWall_type  = 0
        
        logical:: b_purePrtcl   = .false.
	    logical:: b_pureWall    = .false.
	    logical:: b_BnryPrtcl   = .false.
	    logical:: b_BnryWall    = .false.
        
    contains

        procedure:: get_numPrtclType    => bPhP_get_numPrtclType
        procedure:: get_numWallType     => bPhP_get_numWallType
        
        procedure(base_ParticleProperty),deferred:: ParticleProperty
        procedure(base_WallProperty),   deferred:: WallProperty
        procedure(base_PP_BinaryProp),  deferred:: PP_BinaryProp
        procedure(base_PW_BinaryProp),  deferred:: PW_BinaryProp
        procedure(base_get_pure_mass)     ,deferred:: get_pure_mass
        procedure(base_get_pure_inertia)  ,deferred:: get_pure_inertia
    
    end type


    abstract interface
        
        subroutine base_ParticleProperty( this, PSDP )
            import PSD_Property, base_PhysicalProperty
            class(base_PhysicalProperty) this
            class(PSD_Property),intent(in):: PSDP
        end subroutine 
        
        subroutine base_WallProperty( this, n , props )
            import base_PhysicalProperty , IK, PureProperty 
            class(base_PhysicalProperty) this
            integer(IK),intent(in)      :: n
            type(PureProperty),dimension(:),intent(in):: props
        end subroutine
        
            
        subroutine base_PP_BinaryProp(this, DEM_opt , mu , mur, en, et )
            import base_PhysicalProperty, DEMS_Options, RK
            class(base_PhysicalProperty) this
            type(DEMS_Options),intent(in):: DEM_opt
            real(RK),           intent(in):: mu, mur, en, et
        end subroutine
        
        subroutine base_PW_BinaryProp(this, DEM_opt , mu , mur, en, et )
            import base_PhysicalProperty, DEMS_Options, RK
            class(base_PhysicalProperty) this
            type(DEMS_Options),intent(in):: DEM_opt
            real(RK),           intent(in):: mu, mur, en, et
        end subroutine
        
    
        real(RK) function base_get_pure_mass(this, pt )
            import base_PhysicalProperty, RK, IK
            class(base_PhysicalProperty) this
            integer(IK), intent(in) :: pt
        end function
    
        real(RK) function base_get_pure_inertia(this, pt )
            import base_PhysicalProperty, RK, IK
            class(base_PhysicalProperty) this
            integer(IK), intent(in) :: pt
        end function
    
    end interface


    type,extends( base_PhysicalProperty):: PhysicalProperty
    
        real(RK):: rel_vel      = d_rel_vel ! for calculating kn
        type(pureProperty),allocatable,dimension(:)    :: Prtcl_PureProp
        type(pureProperty),allocatable,dimension(:)    :: Wall_PureProp
        type(BinaryProperty),allocatable,dimension(:,:):: Prtcl_BnryProp    
        type(BinaryProperty),allocatable,dimension(:,:):: PrtclWall_BnryProp    
    
    contains
        
        procedure:: ParticleProperty    => PhP_ParticleProperty
        procedure:: WallProperty        => PhP_WallProperty
        procedure:: PP_BinaryProp       => PhP_PP_BinaryProp
        procedure:: PW_BinaryProp       => PhP_PW_BinaryProp
        procedure:: get_pure_mass       => PhP_get_pure_mass
        procedure:: get_pure_inertia    => PhP_get_pure_inertia
        procedure:: pure_getR           => PhP_pure_get_R    
        
        procedure:: Bnry_getProp        => PhP_Bnry_getProp
        procedure:: Bnry_getPropPW      => PhP_Bnry_getPropPW
        procedure:: clc_BnryPrtcl_Prop  => PhP_clc_BnryPrtcl_Prop
        
    end type


    contains
    
    !*****************************************************************************    
    ! extracting the physical property of particles from  PSD_Property object (PSDP)
    !*****************************************************************************
    subroutine PhP_ParticleProperty( this, PSDP )
        implicit none
        class(PhysicalProperty) this
        class(PSD_Property),intent(in):: PSDP
        
        integer(IK) i, num_prop
        
        num_prop = PSDP%get_num_ptype()
        this%numPrtcl_type = num_prop
        if( allocated(this%Prtcl_PureProp) ) deallocate(this%Prtcl_PureProp)
        allocate(this%Prtcl_PureProp (num_prop ) )
        
        do i =1, num_prop
            
            this%Prtcl_PureProp(i) = PSDP%get_prop(i)
            call this%Prtcl_PureProp(i)%clc_prop()
            
        end do
        
        this%b_purePrtcl = .true.
        
    end subroutine 
    
    !*****************************************************************************    
    ! setting the physical property of wall 
    !*****************************************************************************
    subroutine PhP_WallProperty( this, n , props )
        implicit none
        class(PhysicalProperty) this
        integer(IK),intent(in)      :: n
        type(PureProperty),dimension(:),intent(in):: props
        
        !// locals
        integer(IK) i
        
        this%numWall_type = n
        
        if( allocated(this%Wall_PureProp) ) deallocate(this%Wall_PureProp)
        allocate( this%Wall_PureProp(n) )
                   
       do i=1,n
    
	        this%Wall_PureProp(i)%Young_mod = props(i)%Young_mod
	        this%Wall_PureProp(i)%Shear_mod = props(i)%Shear_mod
	        this%Wall_PureProp(i)%YeildStress = props(i)%YeildStress
	        this%Wall_PureProp(i)%poissons_ratio = props(i)%poissons_ratio
	        this%Wall_PureProp(i)%kappa = props(i)%kappa
	        this%Wall_PureProp(i)%Rad  = 10000000.0_RK
	        this%Wall_PureProp(i)%Density  = 10000000.0_RK
	        this%Wall_PureProp(i)%mass = 10000000.0_RK
	        this%Wall_PureProp(i)%Vol  = 10000000.0_RK
	        this%Wall_PureProp(i)%Inertia = 10000000.0_RK
        
       end do
        
        this%b_pureWall = .true.
        
    end subroutine
    
    !*****************************************************************************    
    ! setting the binary properties for particle-particle contacts 
    !*****************************************************************************
    subroutine PhP_PP_BinaryProp(this, DEM_opt , mu , mur, en, et )
        implicit none
        class(PhysicalProperty) this
        type(DEMS_Options),intent(in):: DEM_opt
        real(RK),           intent(in):: mu, mur, en, et
        
        !// locals
        integer(IK) num, i, j
        type(PureProperty) pari, parj
        this%rel_vel = DEM_opt%rel_vel_kn
        
        num = this%get_numPrtclType()
        if( .not. this%b_purePrtcl ) then
            call checkForError(ErrT_Abort , "PhP_PP_BinaryProp", "First initialize pure particle property")    
        end if
                
        if(allocated(this%Prtcl_BnryProp)) deallocate(this%Prtcl_BnryProp)
        allocate( this%Prtcl_BnryProp(num,num) )
        
        do i=1, num
            
            do j=1,num
                
                pari = this%Prtcl_PureProp(i)
                parj = this%Prtcl_PureProp(i)
                !wall = this%Wall_pureProp(j)
                this%Prtcl_BnryProp(i,j) = this%clc_BnryPrtcl_Prop( pari, parj, mu, mur, en, et, DEM_opt%CF_Type ) 
                
            end do
            
        end do
        
        this%b_BnryPrtcl = .true.
        
    end subroutine
    
    !*****************************************************************************    
    ! setting the binary properties for particle-wall contacts 
    !*****************************************************************************
    subroutine PhP_PW_BinaryProp( this, DEM_opt , mu , mur, en, et )
        implicit none
        class(PhysicalProperty) this
        type(DEMS_Options),intent(in):: DEM_opt
        real(RK),           intent(in):: mu, mur, en, et
        
        !// locals
        integer(IK) num, numw, i, j
        type(PureProperty) pari, wall
        this%rel_vel = DEM_opt%rel_vel_kn
        
        num = this%get_numPrtclType()
        numw = this%get_numWallType()
        
        if( .not. this%b_purePrtcl ) then
            call checkForError(ErrT_Abort , "PhP_PW_BinaryProp", "First initialize pure particle property")    
        end if
        
        if( .not. this%b_pureWall ) then
            call checkForError(ErrT_Abort , "PhP_PW_BinaryProp", "First initialize pure wall property")    
        end if
        
        if(allocated(this%PrtclWall_BnryProp)) deallocate(this%PrtclWall_BnryProp)
        allocate( this%PrtclWall_BnryProp(num,numw) )
        
        
        do i = 1, num
            
            do j = 1, numw
                
                pari = this%Prtcl_PureProp(i)
                wall = this%Wall_pureProp(j)
		
                this%PrtclWall_BnryProp(i,j) = this%clc_BnryPrtcl_Prop( pari, wall , mu, mur, en, et, DEM_opt%CF_Type ) 
                
            end do
            
        end do
        
        this%b_BnryWall = .true.
        
        
    end subroutine
    
    real(RK) function PhP_get_pure_mass(this, pt )
        implicit none
        class(PhysicalProperty) this
        integer(IK), intent(in) :: pt
    
        PhP_get_pure_mass = this%Prtcl_PureProp(pt)%mass
        
    end function

    real(RK) function PhP_get_pure_inertia(this, pt )
        implicit none
        class(PhysicalProperty) this
        integer(IK), intent(in) :: pt
    
        PhP_get_pure_inertia = this%Prtcl_PureProp(pt)%inertia
        
    end function

    !*****************************************************************************    
    ! returning physical properties of a pair of particle types (pt_i, pt_j) 
    !*****************************************************************************
    function PhP_Bnry_getProp(this, pt_i, pt_j ) result(res)
        implicit none
        class(PhysicalProperty) this
        integer(IK) pt_i, pt_j
        type(BinaryProperty) res
    
        res = this%Prtcl_BnryProp(pt_i, pt_j)
    
    end function

    !*****************************************************************************    
    ! returning physical properties of particle pt_i and wall w_prop 
    !*****************************************************************************
    function PhP_Bnry_getPropPW(this, pt_i, w_prop ) result(res)
        implicit none
        class(PhysicalProperty) this
        integer(IK) pt_i, w_prop
        type(BinaryProperty) res
    
        res = this%PrtclWall_BnryProp(pt_i,w_prop)
    
    end function


    real(RK) function PhP_pure_get_R(this, pt )
    
    implicit none
    
        class(PhysicalProperty) this
        integer(IK) pt
    
        PhP_pure_get_R = this%Prtcl_PureProp(pt)%Rad
    
    end function

    !*****************************************************************************    
    ! Calculating the binary contact properties
    !*****************************************************************************
    function  PhP_clc_BnryPrtcl_Prop( this, pari, parj, mu, mur, en, et, CF_Type ) result( Bnry )
        implicit none
        class(PhysicalProperty) this
        class(PureProperty),intent(in):: pari, parj
        real(RK),           intent(in):: mu, mur, en, et
        integer(IK),        intent(in):: CF_Type
        type(BinaryProperty) Bnry
        
        real(RK)    kappa
    
        Bnry%Reff = 1.0_RK/(1.0_RK/pari%Rad + 1.0_RK/parj%Rad)
            
        Bnry%meff = 1.0_RK/(1.0_RK/pari%mass + 1.0_RK/parj%mass)
            
        Bnry%Yeff = 1.0_RK/( (1.0_RK-pari%poissons_ratio**2.0_RK)/pari%Young_mod + &
		    (1.0_RK-parj%poissons_ratio**2.0_RK)/parj%Young_mod  ) 
            
        Bnry%Geff = 1.0_RK/( (2.0_RK-pari%poissons_ratio)/pari%Shear_mod +   &
                    (2.0_RK-parj%poissons_ratio)/parj%Shear_mod  ) 
            
        Bnry%kn = 1.2024*(Bnry%meff**0.5_RK * Bnry%Yeff**2 * Bnry%Reff * this%rel_vel )**(2.0_RK/5.0_RK) 
            
        kappa = ( (1.0_RK - pari%poissons_ratio)/pari%Shear_mod + (1.0_RK - parj%poissons_ratio)/parj%Shear_mod ) / &
                ( (1.0_RK - 0.5_RK*pari%poissons_ratio)/pari%Shear_mod + (1.0_RK - 0.5_RK*parj%poissons_ratio)/parj%Shear_mod )
        Bnry%kt = kappa*Bnry%kn
    
        Bnry%en = en
        Bnry%et = et
        Bnry%fric = mu
        Bnry%roll_fric = mur
        
        Bnry%dmp_n = Damping_normal(Bnry , CF_Type) 
        Bnry%dmp_t = Damping_tangential(Bnry , CF_Type) 
    	
    end function


    function Damping_normal( Bnry, CF_Type ) result (res)
        implicit none  
        type(BinaryProperty) bnry
        integer(IK) CF_Type
        real(RK) res
        
        real(RK) K_hertz
        select case( CF_Type )
        case(CFT_LSD_nl, CFT_LSD_l)
        
            res = -2.0*log(Bnry%en)*sqrt(Bnry%meff*Bnry%kn)/sqrt( (log(Bnry%en))**2 + pi**2.0_RK);
            return
        case( CFT_nLin_nl, CFT_nLin_l )
            
            K_hertz = 4.0_RK/3.0_RK*Bnry%Yeff*sqrt(Bnry%Reff); 
            res = -2.2664*log(Bnry%en)*sqrt(Bnry%meff*K_hertz)/sqrt( log(Bnry%en)**2 + 10.1354);
            
        case default
    
            res = Bnry%dmp_n
            call checkForError(ErrT_Abort , "Damping_normal", "THis contact force model is not specified in the program")
            
        end select
    
    end function

    function Damping_tangential( Bnry, CF_Type ) result (res)
        implicit none
        type(BinaryProperty) bnry
        integer(IK) CF_Type
        real(RK) res
        
        
        res = 0.0
        print*, "No equation is considered for tangential damping yet"
        
    end function


    integer(IK) function bPhP_get_numPrtclType(this)
        implicit none
        class(base_PhysicalProperty) this
    
        bPhP_get_numPrtclType = this%numPrtcl_type
    
    end function

    integer(IK) function bPhP_get_numWallType(this)
        implicit none
        class(base_PhysicalProperty) this
    
        bPhP_get_numWallType = this%numWall_type
    
    end function


end module
