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

program main

    use g_DEMSystem
    use g_stl_reader
    use g_Prtcl_PureProperty
    implicit none
    
    integer(IK) i
    real(RK)    Y, pr
    logical     res
    type(real3) minDomain, maxDomain
    type(PS_Distribution),  pointer:: PSD
    type(PSD_Property),     pointer:: PSDP
    type(PSDP_Position),    pointer:: Pos
    class(PhysicalProperty),pointer:: Property
    class(Geometry),        pointer:: geom
    type(real3)         p1, p2, p3, p4
    type(PlaneWall)     insPlane
    type(PureProperty)  prop
    type(DEMS_Options)  DEM_opt 
    type(DEMSystem)     DEM
    
!//// Initial settings 

    ! options and settings for DEM simulation
    DEM_opt%gravity   = real3( 0.0, -9.8, 0.0 )
    DEM_opt%SaveFreq  = 1000
    DEM_opt%RunName   = "mixingRotatingDrumRolling"
    DEM_opt%Res_Dir   = "./Results/"
    DEM_opt%OutputFileType = OP_Type_VTK
    DEM_opt%CT_model  = CTM_ConstantTorque
    DEM_opt%CF_Type   = CFT_nLin_nl
    DEM_opt%CS_Method = CSM_NBS_Munjiza !_Hrchl !
    DEM_opt%PI_Method = PIM_AB4AM5
    DEM_opt%PRI_Method= PIM_AB4AM5
    
    !//// properties of particles and walls 
    Y = 1000000.0_RK
    pr = 0.23_RK
    call prop%set_prop(2500.0_RK, Y , Y/(2*(1+pr)), pr , 0.0_RK)
    
    !//// particles
    allocate(PSD)
    
    ! 19000 particles with diameter of 3 mm
    PSD  = PS_Distribution( 19000 , PSD_Uniform, 1 , 0.003_RK, 0.003_RK)
    
    ! associates particle size distribution with property
    allocate(PSDP)
    PSDP = PSD_Property( PSD , prop )
    deallocate(PSD)
   
    ! positions particles inside the drum
    allocate(Pos)
    Pos = PSDP_Position( PSDP )
    call Pos%PositionOrdered( real3(-0.071 , -0.071, 0.0), real3(0.071, 0.071, 0.03) )
    
!// main components of DEM system    
    !//// geometry
    allocate(geom)
    call ProgramDefinedGeometry(geom)
            
    !//// Property
    allocate( Property )
    call Property%ParticleProperty( PSDP )    ! particles 
    call Property%WallProperty(1 , (/prop/) ) ! walls 
    call Property%PP_BinaryProp(DEM_opt, 0.5_RK, 0.1_RK, 0.8_RK, 0.8_RK) ! binary pp
    call Property%PW_BinaryProp(DEM_opt, 0.7_RK, 0.1_RK, 0.8_RK, 0.8_RK) ! binary pw
    
    !//// DEM system with particle insertion
    ! simulation domain
    minDomain = real3(-0.11, -0.11,  -0.01)
    maxDomain = real3( 0.11,  0.11 ,  0.04)
    
    ! initializes DEM system, dt = 0.00002 sec.  
    call DEM%Initialize( 0.00002_RK, Pos ,geom, Property, minDomain, maxDomain, DEM_opt ) 
    
!//// iteration loop for 0.2 seconds. This lets particles settle under the gravity   
    call DEM%iterate(10000)
    
   
    !// sets the rotating velocity of drum, it is 10 RPM
    call geom%setWallVelocity(1, real3(0,0,0), &
                              RPMtoRAD_S(10.0_RK), &
                              p_line( real3(0,0,0) , real3(0,0,0.2) ) )
    
 !//// iteration for 0.8 more seconds to reach steady condition
    call DEM%iterate(39999)
    !// marks three groups of particles as tracers  
    call DEM%User_prtclMark()
    call DEM%iterate(1)
    
    
    !/// iteration for 20 more seconds 
    call DEM%iterate(1000000)
    
end program 
