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
    DEM_opt%gravity   = real3( 0.0, 0.0, -9.8 )
    DEM_opt%SaveFreq  = 1000
    DEM_opt%RunName   = "Continuous_Mixer"
    DEM_opt%Res_Dir   = "./Results/"
    DEM_opt%OutputFileType = OP_Type_VTK   ! vtk output
    DEM_opt%CT_model  = CTM_ConstantTorque ! constant torque model
    DEM_opt%CF_Type   = CFT_nLin_nl        ! non-linear model
    DEM_opt%CS_Method = CSM_NBS_Munjiza_Hrchl !for polydisperse system
    DEM_opt%PI_Method = PIM_AB4AM5
    DEM_opt%PRI_Method= PIM_AB4AM5
    
    !//// properties of particles and walls 
    Y = 1000000.0_RK
    pr = 0.23_RK
    call prop%set_prop(2500.0_RK, Y , Y/(2*(1+pr)), pr , 0.0_RK)
    
    !//// particles
    allocate(PSD)
    
    ! tri-modal distribution with size range between 3 to 6 mm                          
    PSD  = PS_Distribution( 600000 , PSD_Uniform, 1 , 0.006_RK, 0.006_RK)
    
    ! associates particle size distribution with property
    allocate(PSDP)
    PSDP = PSD_Property( PSD , prop )
    deallocate(PSD)
   
        
!// main components of DEM system    
    !//// geometry
    allocate(geom)
    ! reads geometry from stl files                                                    
    call geom%add_stl_file("./body.stl", 1, 1, .true. ) ! shell, user_id = 1
    call geom%add_stl_file("./impeller.stl", 2, 1, .true. ) ! screw, user_id = 2
    call geom%add_stl_file("./cylender.stl", 3, 1, .true. ) ! cylender, user_id = 3
 
    ! separated output for geometry
    ! user_id = 1: separate files with a name containing "body"
    ! user_id = 2: separate files with a name containing "impeller"
    call geom%setWallOutputName( (/1,2,3/) , (/"body","impe","cyle"/) )
            
    !//// Property 
      !!!DEM_Options ,dynamic friction, rolling friction, coefficient of normal restitution, coefficient of tangential restitution 
    allocate( Property )
    call Property%ParticleProperty( PSDP )    ! particles 
    call Property%WallProperty(1 , (/prop/) ) ! walls 
    call Property%PP_BinaryProp(DEM_opt, 0.3_RK, 0.1_RK, 0.7_RK, 0.7_RK) ! binary pp       !!! inaro chi bezaram??
    call Property%PW_BinaryProp(DEM_opt, 0.4_RK, 0.1_RK, 0.7_RK, 0.7_RK) ! binary pw
    
    !//// DEM system with particle insertion
    ! simulation domain
    minDomain = real3(-0.11 , -0.01,  -0.13)
    maxDomain = real3( 0.11,  0.65 ,  0.12)
    
    ! insertion plane near the inlet gate of conveyor
    p1 = real3( -0.050 , 0.065, 0.101)
    p2 = real3(  0.050, 0.065,0.101)
    p3 = real3(  0.050 , 0.010,0.101)
    p4 = real3( -0.050 , 0.010, 0.101)
    res = insPlane%CreateWall_nv(p1,p2,p3,p4)
  !!!----------------------------------------------------------------------------------------------------------------------------
    ! initializes DEM system, dt = 0.00001 s, particles should be inserted during 2              !!!in chie?
    ! million iterations (equivalent to 20 seconds)
    call DEM%Initialize( 0.00001_RK, PSDP, InsPlane, 6000000, 0.3_RK ,geom,  Property, minDomain, maxDomain, DEM_opt )
                                                            !initial velocity m/s
 
!//// iteration loop for 0.1 seconds  
    call DEM%iterate(1)
        
    ! sets the velocity of the impeller (user_id = 2)
    ! rotation speed is 60 RPM and rotation axis is the central shaft    
    call geom%setWallVelocity( 2, real3(0,0,0), &
                             RPMtoRAD_S(-120.0_RK), &
                            p_line( real3(0,0,0) , real3(0,0.75,0.0) ) )
    
!/// iteration for 20 seconds
	call DEM%User_prtclMark(); 
	
	do i = 1,20000000
		call DEM%User_prtclMark(); 
	    call DEM%iterate(10)
	end do
        
    
end program 
