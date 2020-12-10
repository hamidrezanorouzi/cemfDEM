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
    real(RK)    Y, pr, G
    logical     res
    type(real3) minDomain, maxDomain
    type(PS_Distribution),  pointer:: PSD
    type(PSD_Property),     pointer:: PSDP
    class(PhysicalProperty),pointer:: Property
    class(Geometry),        pointer:: geom
    type(real3)         p1, p2, p3, p4
    type(PlaneWall)     insPlane
    type(PureProperty)  prop
    type(DEMS_Options)  DEM_opt 
    type(DEMSystem)     DEM
    
!//// Initial settings 

    ! options and settings for DEM simulation
    DEM_opt%gravity   = real3( 0.0, -9.8 , 0.0 )
    DEM_opt%SaveFreq  = 1000
    DEM_opt%RunName   = "PackingBinary"
    DEM_opt%Res_Dir   = "./Results/"
    DEM_opt%OutputFileType = OP_Type_VTK !+  OP_Type_Tec
    DEM_opt%CT_model  = CTM_ConstantTorque
    DEM_opt%CF_Type   = CFT_LSD_l
    DEM_opt%CS_Method = CSM_NBS_Munjiza_Hrchl ! for a system with size distribution
    DEM_opt%PI_Method = PIM_AB3AM4
    DEM_opt%PRI_Method= PIM_AB3AM4
    
    !//// properties of particles and walls 
    Y = 1000000.0_RK ! Young's modulus 
    pr = 0.23_RK     ! Poisson's ratio
    G = Y/(2*(1+pr)) ! shear modulus 
    call prop%set_prop(2500.0_RK, Y , G , pr , 0.0_RK)
    
    !//// particles
    allocate(PSD)
    
    ! bimodal particles with d1 = 7 mm and d2 = 3.5 mm
    ! total number of particles is (7+8)*1000 = 15K
    !PSD  = (PS_Distribution( 7 , PSD_Uniform, 1 , 0.007_RK, 0.007_RK ) +  &
    !        PS_Distribution( 8 , PSD_Uniform, 1 , 0.0035_RK, 0.0035_RK ) ) * 1000

    !particles with uniform size distribution
    !Size range of particles is between 3 and 7.5 mm which are classified into 5 bins
    PSD  = PS_Distribution( 15000 , PSD_Uniform, 5 , 0.003_RK, 0.0075_RK )

    ! associates particle size distribution with a property set
    allocate(PSDP)
    PSDP = PSD_Property( PSD , prop )
    deallocate(PSD)
    
!// main components of DEM system    
    !//// geometry
    allocate(geom)
    call ProgramDefinedGeometry(geom)
            
    !//// Property
    allocate( Property )
    call Property%ParticleProperty( PSDP )    ! particles 
    call Property%WallProperty(1 , (/prop/) ) ! walls 
    call Property%PP_BinaryProp(DEM_opt, 0.2_RK, 0.1_RK, 0.8_RK, 0.8_RK) ! binary pp
    call Property%PW_BinaryProp(DEM_opt, 0.3_RK, 0.1_RK, 0.8_RK, 0.8_RK) ! binary pw
    
    !//// DEM system with particle insertion
    ! simulation domain
    minDomain = real3(-0.077, -0.01,  -0.077)
    maxDomain = real3( 0.077 , 0.21 ,  0.077)
    
    ! insertion plane
    p1 = real3(-0.05, 0.18, -0.05)
    p2 = real3(-0.05, 0.18,  0.05)
    p3 = real3( 0.05, 0.18,  0.05)
    p4 = real3( 0.05, 0.18, -0.05)
    res = insPlane%CreateWall_nv(p4,p3,p2,p1)
    
    ! initializes DEM system, dt = 0.00001 sec, particles are inserted in 150k iterations with initial velocity of 0.3 m
    call DEM%Initialize( 0.00001_RK, PSDP, insPlane, 150000 , 0.3_RK ,geom, Property, minDomain, maxDomain, DEM_opt ) 
    
    
!//// iteration loop for 2.5 seconds  
    call DEM%iterate(250000)    
    
   
end program 
