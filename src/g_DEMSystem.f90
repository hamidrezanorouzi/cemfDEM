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
!  file name  : g_DEMSystem.f90                      
!  module name: g_DEMSystem                          
!  
!  Purpose:                                          
!   Containing all objects which are necessary to do iterative DEM calculations
!    for non-cohesive particles
!   Performing iteration over time                   
!   It contains methods to initialize simulation     
!   It contains methods to write results in an external file
!                                                           
!  Note:                                                    
!   Most of important data members are public, but modifying them outside the
!    scope may have side effects.
! 
!------------------------------------------------------------------------------    
    
module g_DEMSystem

    use g_Prtcl_Integration
    use g_ContactSearch
    use g_ContactSearchPW
    use g_Prtcl_LSD_Model
    use g_Prtcl_NonLin_Model
    use g_WallOutput
    use g_timer
    use g_Prtcl_Insertion
   
    
    
    implicit none
    
    real(RK),parameter:: max_omega = 500*pi
    
    
    interface 
        subroutine ProgramDefinedGeometry( Geom )
            import Geometry
            class(Geometry),intent(in):: Geom
        end subroutine
    end interface
    
    interface 
        function User_Mark( id, flag, ptype, dpos, oldMark ) result (mark)
            import IK, real4
            integer(IK),intent(in):: id, flag, ptype, oldMark
            type(real4),intent(in):: dpos 
            integer(IK) mark
        end function
    end interface 
    
    !// DEMSystem class 
    type DEMSystem
        
        ! iteration number
        integer(IK)     :: iterNumber        = 0
        
        ! time step of integration
        real(RK)        :: prtcl_dt          = d_prtcl_dt    
        
        ! options which can be modified to choose different
        ! components of program
        type(DEMS_Options) m_DEM_Options
        
        ! determining the particle insertion method
        integer(IK) Prtcl_InsertionMethod
        
        ! object which manages particle insertion, if any
        type(Prtcl_Insertion) Prtcl_Ins
        
        
        !//  all data required to represent spherical particles
        logical:: Prtcls_init = .false.
        integer(IK) max_nPrtcl  ! maximum number of particles
        integer(IK) numPrtcl    ! current number of particles
        integer(IK),pointer,dimension(:):: prtcl_ids    ! ids of particles 
        integer(IK),pointer,dimension(:):: prtcl_flag   ! flags of particles: in domain, deleted or not inserted 
        integer(IK),pointer,dimension(:):: prtcl_type   ! property type of particles 
        integer(IK),pointer,dimension(:):: prtcl_usr_mark ! is used for marking particles in the system 
        type(real4),pointer,dimension(:):: prtcl_dPos   ! position and diameter of particles
        type(real3),pointer,dimension(:):: prtcl_vel    ! linear velocity of particles
        type(real3),pointer,dimension(:):: prtcl_teta   ! angular position of particles
        type(real3),pointer,dimension(:):: prtcl_rvel   ! angular velocity of particles
        type(real3),pointer,dimension(:):: prtcl_cntct_force    ! sum of contact forces on particles
        type(real3),pointer,dimension(:):: prtcl_torque         ! sum of torques acting on particles
        type(real3),pointer,dimension(:):: prtcl_Lin_acc        ! linear acceleration of particles
        type(real3),pointer,dimension(:):: prtcl_rot_acc        ! rotational acceleration of particles
        
        
        !// Geometry object 
        logical:: Geometry_init = .false.
        class(Geometry),pointer         :: m_Geometry
        
        ! Physical property object 
        logical:: PhysProp_init = .false.
        class(base_PhysicalProperty),pointer  :: m_PhysProp
        
        !// Contact force object 
        logical:: ContForce_init = .false.
        class(base_ContactForce),pointer:: m_Cont_force
        
        ! contact lists objects 
        logical:: ContList_init = .false.
        class(base_ContactList),pointer :: PP_Cont_list ! particle-particle
        class(base_ContactList),pointer :: PW_Cont_list ! particle-wall
        
        ! // Contact search objects 
        logical:: ContSearch_Init = .false.
        class(ContactSearch),pointer    :: Cont_Search      ! particle-particle
        class(ContactSearchPW),pointer  :: Cont_Search_PW   ! particle-wall
        
        !// Integration objects 
        logical:: Integration_init = .false.
        class(Prtcl_Integration),pointer    :: Lin_Integration  ! linear integration
        class(Prtcl_rot_Integration),pointer:: rot_Integration  ! rotational integration
        
        
        !// timers
        
        type(timer) m_total_timer
        type(timer) m_pre_iter_timer
        type(timer) m_prediction_timer
        type(timer) m_ContSearchPP_timer
        type(timer) m_ContSearchPW_timer
        type(timer) m_ForcePP_timer
        type(timer) m_ForcePW_timer
        type(timer) m_Acceleration_timer
        type(timer) m_integration_timer
        type(timer) m_write_prtcl_timer

        
    contains
    
    ! initializing DEMSystem object based on two methods
    procedure:: DEMS_Initialize1 ! inserting particles from a plane
    procedure:: DEMS_Initialize2 ! position of particles are known and all are inserted at once
    generic  :: Initialize => DEMS_Initialize2, DEMS_Initialize1
    
    ! iterating simulation for n time steps
    procedure:: iterate     => DEMS_iterate
    
    ! performing pre-iterations 
    procedure:: preIteration    => DEMS_preIteration
    
    ! TecPlot file output for geometry and particels
    procedure:: write_prtcl     => DEMS_write_prtcl
    procedure:: write_geom      => DEMS_write_geom
    
    procedure:: writeRawData    => DEMS_WriteRawData
    
    ! VTK file output for particles
    procedure:: PrtclToVTK      => DEMS_PrtclToVTK
    
    ! writing to output file 
    procedure:: outputfile     => DEMS_outputfile
    
    ! terminal output 
    procedure:: terminalOutput => DEMS_terminalOutput 
    
    ! writing packed bed condition to a file
    procedure:: WritePackedBed  => DEMS_WritePackedBed
    
    !calculating acceleration of particles - linear and angular
    procedure:: clc_Acceleration    => DEMS_clc_Acceleration
    
    
    ! setting new number of particles: is related to particle insertion
    procedure:: set_numPrtcl    => DEMS_set_numPrtcl
    
    
    ! setting geometry object
    procedure::DEMS_SetGeometry
    generic:: SetGeometry => DEMS_SetGeometry
    
    
    ! Initializing contact force object 
    procedure:: DEMS_SetContactForce 
    generic:: SetContactForce => DEMS_SetContactForce
   
    ! Initializing contact search object
    procedure:: DEMS_setContactSearch
    generic:: setContactSearch  => DEMS_setContactSearch
    
    ! Initializing integeration object
    procedure:: DEMS_setIntegration
    generic:: setIntegration    => DEMS_setIntegration
    
    ! returning number of particles in the simulation domain
    procedure:: numInDomain     => DEMS_numInDomain
    
    ! a procedure to mark particles based on the user defined method 
    procedure:: User_prtclMark  => DEMS_User_prtclMark
    
    ! writing 
    procedure:: WriteParticleData => DEMS_WriteParticleData
    
    ! final procedure for DEMsystem
    final:: Finilize_DEMSystem
    
    
    ! //// private components
    procedure,private:: DEMS_Init_rest
    procedure,private:: DEMS_timers_reset
    procedure,private:: DEMS_timers_output
    procedure,private:: DEM_GetParticles
    procedure,private:: DEM_reallocateParticles
    procedure,private:: DEMS_deallocateAll
    ! reports contact forces 
    procedure,private:: DEMS_report_contacts
    
    end type
    !// DEMSystem class
    
    
contains

!********************************************************************************
!   Initializing DEMSystem object with particles which are inserted from a 
!   predefined plane 
!********************************************************************************
subroutine DEMS_Initialize1(this, dt , Particles, Ins_Plane, num_steps, ins_vel ,Geom, Property, minDomain, maxDomain , DEM_opt )
    implicit none
    class(DEMSystem)                     this
    real(RK),               intent(in):: dt         ! time step
    class(PSD_Property),    intent(in):: Particles  ! Particles to be inserted
    class(PlaneWall),       intent(in):: Ins_Plane  ! plane wall through which particels should be inserted
    integer(IK),            intent(in):: num_steps  ! number of steps that particles should be inserted
    real(RK),               intent(in):: ins_vel    ! insertion velocity of particles
    class(Geometry),pointer,intent(in):: Geom       ! Geometry object                
    class(PhysicalProperty),pointer,intent(in):: Property  ! Property object
    type(real3),            intent(in):: minDomain, maxDomain   ! minimum and maximum points of simulation domain    
    type(DEMS_Options),     intent(in):: DEM_opt                ! DEM options which are used to choose among different components of program
    
    
    ! getting time step and DEM_Options
    this%prtcl_dt       = dt
    this%m_DEM_Options  = DEM_opt
    
    
    !// Initializing log info 
    MainLogInfo = LogInfo(d_MainLogInfo_unit, this%m_DEM_Options%Res_Dir, &
                          this%m_DEM_Options%RunName, this%m_DEM_Options%LF_file_lvl, &
                          this%m_DEM_Options%LF_cmdw_lvl ) 
    
    
    ! // starting of initializing DEMSystem 
    
    !**************************************    
    ! Step1: Initializing particles and the method they are introduced into the simulation
    
    call MainLogInfo%OutInfo("Step1: Particles are reading into DEMSystem ...", 1 )
    
    ! maximum number of particles
    this%max_nPrtcl = Particles%get_numPrtcl()
    ! number of particles. No particle is inserted yet
    this%numPrtcl = 0
    
    ! allocating memory for particles
    call this%DEM_reallocateParticles()
    
    !transferring data from "Particles" to this object
    call this%DEM_GetParticles(Particles)
    
    ! initializing insertion of particles
    call this%Prtcl_Ins%InitInsertion(Ins_Plane, num_steps , this%max_nPrtcl, ins_vel)
    
    ! this means that there is no particle inserted into the simulation
    this%Prtcl_InsertionMethod = PInsM_Box
    this%Prtcl_flag = Pflg_notInserted  
     
    !>>> log info 
     call MainLogInfo%OutInfo( "Particle insertion method is set to insertion from a plane", 2 )
     call MainLogInfo%OutInfo( "Number of particles available in the system:"// trim(num2str(this%max_nPrtcl)) ,2 )
    !>>> log info
     
     this%Prtcls_init = .true.
     
        
    ! does the rest of initialization 
    call this%DEMS_Init_rest( Geom, Property, minDomain, maxDomain , DEM_opt )
    
    
    return
    
end subroutine


!********************************************************************************
!   Initializing DEMSystem object with particles which are inserted from a 
!   predefined plane 
!********************************************************************************
subroutine  DEMS_Initialize2(this, dt , Particles, Geom, Property, minDomain, maxDomain , DEM_opt )
    implicit none
    class(DEMSystem)                     this
    real(RK),               intent(in):: dt         ! time step
    class(PSDP_Position),   intent(in):: Particles  ! particles with their positions
    class(Geometry),pointer,intent(in):: Geom       ! geometry object
    class(PhysicalProperty),pointer,intent(in):: Property  ! property object
    type(real3),            intent(in):: minDomain, maxDomain   ! minimum and maximum points of simulation domain
    type(DEMS_Options),     intent(in):: DEM_opt                ! DEM options which are used to choose among different components of program
    
    
    
    
    this%prtcl_dt = dt
    this%m_DEM_Options = DEM_opt
    
    !// Initializing main log info
    MainLogInfo = LogInfo(d_MainLogInfo_unit, this%m_DEM_Options%Res_Dir, &
                          this%m_DEM_Options%RunName, this%m_DEM_Options%LF_file_lvl, &
                          this%m_DEM_Options%LF_cmdw_lvl ) 
    
    !**************************************    
    ! Step1: Initializing particles and the method they are introduced into the simulation
    
    call MainLogInfo%OutInfo("Step1: Particles are reading into DEMSystem ...", 1 )
    
    ! maximum number of particles
    this%max_nPrtcl = Particles%get_numPrtcl()
    ! number of particles
    this%numPrtcl = this%max_nPrtcl
    
    ! allocating memory for particles
    call this%DEM_reallocateParticles()
    
    !transferring data from "Particles" to this object
    call this%DEM_GetParticles(Particles)
    
     this%Prtcl_InsertionMethod = PInsM_File
     this%Prtcl_flag = Pflg_inDomain  ! all particles are assumed to be in domain
     
     !>>> log info 
     call MainLogInfo%OutInfo( "Particle insertion method is set to normal state", 2 )
     call MainLogInfo%OutInfo( "Number of particles avaiable in the system:"// trim(num2str(this%numPrtcl)) ,2 )
     !>>> log info
     
     this%Prtcls_init = .true.
     
     ! does the rest of initialization 
     call this%DEMS_Init_rest( Geom, Property, minDomain, maxDomain , DEM_opt )
    
    return
    
end subroutine

!******************************************************************************
! rest of initialization of DEMSystem
!******************************************************************************
subroutine DEMS_Init_rest(this, Geom, Property, minDomain, maxDomain , DEM_opt )
    implicit none
    class(DEMSystem)                     this
    class(Geometry),pointer,intent(in):: Geom
    class(PhysicalProperty),pointer,intent(in):: Property
    type(real3),            intent(in):: minDomain, maxDomain
    type(DEMS_Options),     intent(in):: DEM_opt
    
    !// locals
    
    !// body
    this%m_DEM_Options = DEM_Opt
    this%prtcl_usr_mark = 1
    
    !Step2: Geometry
    !>>> log info 
     call MainLogInfo%OutInfo("Step2: Geometry is set", 1 )
    !>>> log info 
    
    this%m_Geometry=>geom
    this%Geometry_init = .true.
    
    !>>> log info 
    call MainLogInfo%OutInfo( "Geometry Contains "// trim(num2str( this%m_Geometry%get_num_pWall() )) // &
                              " Plane walls.", 2 )
    !>>> log info
    
    
    !Step3: Physical property
    
    this%m_PhysProp => Property
    this%PhysProp_init = .true.
    
    !>>> log info 
    call MainLogInfo%OutInfo("Step3: Physical properties of particels and walls are set.",1)
    
    call MainLogInfo%OutInfo("Physical properties contains "// &
                             trim( num2str(this%m_PhysProp%get_numPrtclType() ) ) // &
                             " particle types and "//trim( num2str(this%m_PhysProp%get_numWallType() ) )// " wall types.", 2)
    
    !<<< log info
    
    ! Step4: Initializing contact force
   
    
    !>>> log info 
    call MainLogInfo%OutInfo("Step4: Initializing contact force models . . . ", 1 )
    !<<< log info 
    
    call this%setContactForce()
    
    
    ! Step5: Initializing contact search method 
    call MainLogInfo%OutInfo("Step5: Initializing contact search method . . . ", 1 )
    
    call this%setContactSearch( minDomain, maxDomain )
    
    
    ! Step6: Initializing integration methods
    call MainLogInfo%OutInfo("Step6: Initializing integration method . . . ", 1 )
    
    
    call this%setIntegration()
    
    ! Step7: timers for recording the execution time of different parts of program
    call MainLogInfo%OutInfo("Step7: Initializing timers . . . ", 1 )
    call this%DEMS_timers_reset()
    
end subroutine


!********************************************************************************
!   iterating over time 
!   calls all the required methods to do numIter iterations in the DEM system
!********************************************************************************
subroutine DEMS_iterate( this, numiter )
    implicit none
    class(DEMSystem) this
    integer(IK),intent(in):: numiter ! number of time steps 
    
    !//locals
    integer(IK) i, j
    integer(IK),save:: k = 0
    integer(IK) nfrom, nto
    
    
    
    do i=1,numIter
        
        call this%m_total_timer%start()
        
        ! pre-iteration adjustments 
        call this%m_pre_iter_timer%start()
            call this%preIteration()
        call this%m_pre_iter_timer%finish()
        
        ! first predicting the position and velocity of particles based on the previous time step 
        call this%m_prediction_timer%start()
            call this%Lin_Integration%predict()
            call this%rot_Integration%predict()
        call this%m_prediction_timer%finish()
    	
	
        ! finding contacts between particles 
        call this%m_ContSearchPP_timer%start()
            call this%Cont_Search%FindContacts()
        call this%m_ContSearchPP_timer%finish()
        
        ! finding contacts between particles and walls
        call this%m_ContSearchPW_timer%start()
            call this%Cont_Search_PW%FindContacts(this%numPrtcl , this%prtcl_dpos , this%prtcl_ids)
        call this%m_ContSearchPW_timer%finish()            
        
        ! calculating contact forces and torques between particles 
        call this%m_ForcePP_timer%start()
            call this%m_cont_force%AllContactForce_PP()
        call this%m_ForcePP_timer%finish()
        
        ! calculating contact forces and torques between particles and walls
        call this%m_ForcePW_timer%start()
            call this%m_cont_force%AllContactForce_PW()
        call this%m_ForcePW_timer%finish()

        ! calculating linear and angular accelerations 
        call this%m_Acceleration_timer%start()
            call this%clc_Acceleration()
        call this%m_Acceleration_timer%finish()
    	

        ! correcting position and velocities 
        call this%m_integration_timer%start()
            call this%Lin_Integration%correct()
            call this%rot_Integration%correct()
        call this%m_integration_timer%finish()
	 
	
        ! moving walls if any wall is moving 
        call this%m_Geometry%move_walls(this%prtcl_dt)
        
        this%iterNumber = this%iterNumber + 1
        
        ! writing results to the output file 
        call this%m_write_prtcl_timer%start()    
            call this%outputFile()         
        call this%m_write_prtcl_timer%finish()
            
       
        
        call this%m_total_timer%finish()
        
        ! output to log file and terminal/command window 
        call this%terminalOutput()   
       
    end do
    
end subroutine 


subroutine DEMS_preIteration(this)
    implicit none
    class(DEMSystem) this
      
    integer(IK) numIter
    integer(IK) old_nPrtcl
    
    ! first checking if there are any particles to be inserted in the simulation domain
    
    if( this%Prtcl_insertionMethod == PInsM_Box) then
        
        old_nPrtcl = this%numPrtcl
        
        if( this%Prtcl_Ins%InsertParticle( this%numPrtcl, this%prtcl_dt, &
            this%prtcl_dpos , this%prtcl_flag , this%prtcl_vel , this%iterNumber ) ) then
            
	    
            call MainLogInfo%OutInfo("New particles have been inserted : "// &
                                      trim(num2str(this%numPrtcl-old_nPrtcl))//" particles", 2 )
            call MainLogInfo%OutInfo("Total number of particles is     : "//trim(num2str(this%numPrtcl)) , 2 )
            
            ! forcing the particle-wall contact search method to update the neighbor list of walls
            call this%Cont_Search_PW%Reset_update()
            
            ! setting new number of particles 
            call this%set_numPrtcl()
            
        end if
        
    end if
    
            
    ! wall neighbor list
    if( this%Cont_Search_PW%FindNearPrtcls( this%prtcl_dt, this%iterNumber, &
        this%numPrtcl, this%prtcl_dPos, this%prtcl_flag, this%prtcl_vel, this%prtcl_Lin_acc) ) then
            
        call MainLogInfo%OutInfo("Neighbor particles of all walls are updated", 4 )
        call MainLogInfo%OutInfo("This process is repeated in Iteration :"// &
                                  trim(num2str( this%Cont_Search_PW%Next_update() )),4)
        
    end if
    
    
    this%prtcl_cntct_force  = zero_r3
    this%prtcl_torque       = zero_r3
    
    call this%PP_Cont_List%setNumCntcts(0)
    call this%PW_Cont_List%setNumCntcts(0)
    
end subroutine 




subroutine DEMS_set_numPrtcl(this)
    
    implicit none
    class(DEMSystem) this
    
        
    ! integrating
    call this%Lin_Integration%set_numPrtcl(this%numPrtcl)
    call this%rot_Integration%set_numPrtcl(this%numPrtcl)
    
    
    ! Contact search
    call this%Cont_Search%set_numPrtcl(this%numPrtcl)
    
          
    ! Contact force 
    call this%m_Cont_force%set_numPrtcl(this%numPrtcl)
    
    
end subroutine 



!**********************************************************************
! writing particle data in a TecPlot data file
!**********************************************************************
subroutine DEMS_write_prtcl( this , file_num )
    implicit none
    class(DEMSystem) this
    integer(IK), intent(in):: file_num
    
    !// locals
    integer(IK) i, numPrtcl, num_inDomain
    integer(IK) nUnitFile
    integer(IK),save:: strnd_id = 0
    real(RK)     sol_time
    character(128) chFile, ch
    
    ! steps to produce file name
    write(ch,*) this%m_DEM_Options%base_number+file_num
       
    nUnitFile = d_tec_prtcl_unit
    
    chFile = trim(this%m_DEM_Options%Res_Dir)//'PrtclData- '// &
             trim(adjustl(ch))//"for"//trim(this%m_DEM_Options%RunName)//".plt"
    open( unit =  nUnitFile , file = chfile )
    
    
    ! counting number of particles in the domain
    numPrtcl = this%numPrtcl
    num_inDomain = this%numInDomain() 
       
    ! writing header of the file
    sol_time = this%iterNumber * this%prtcl_dt
    strnd_id =   1
    write( nUnitFile, *) "Variables = x, y, z, d, id, usr_mark , V , Vx, Vy, Vz, w,  wx, wy, wz, fc, fcx, fcy, fcz" 
    write( nUnitFile, *)"zone I = ", num_inDomain ,", SOLUTIONTIME = ", sol_time , ", STRANDID = ", strnd_id 
    
    ! writing particle data in the file 
    do i = 1, numPrtcl
        if( this%prtcl_flag(i) >= Pflg_inDomain ) then
            write( nUnitFile,*) this%prtcl_dPos(i), this%prtcl_ids(i), this%prtcl_usr_mark(i), &
                                norm(this%prtcl_vel(i)),this%prtcl_vel(i), &
                                norm(this%prtcl_rvel(i)), this%prtcl_rvel(i), &
                                norm(this%prtcl_cntct_force(i)), this%prtcl_cntct_force(i)
        end if
    end do
    
    close( nUnitFile)
    
end subroutine 


!**********************************************************************
! writing geometry data in a TecPlot data file
!**********************************************************************
subroutine DEMS_write_geom( this, num )
    implicit none
    class(DEMSystem) this
    integer(IK),intent(in):: num
    
    !// locals 
    integer(IK) i, nw
    integer(IK) nUnitFile
    integer(IK),save:: strnd_id=0
    real(RK) sol_time
    character(256) chFile, ch
    type(WallOutput), save:: TecWall
    
    
    call TecWall%reset()
    
    ! plane walls
    nw = this%m_Geometry%get_num_pWall ()
    do i = 1, nw

        call TecWall%Add_PlaneWall( this%m_Geometry%pWall(i) )
    
    end do

    
    write(ch,*) this%m_DEM_Options%base_Number + num
    chFile = trim(this%m_DEM_Options%Res_Dir)//'Walls'//trim(adjustl(ch))//"for"//trim(this%m_DEM_Options%RunName)//".plt"
    sol_time = this%iterNumber * this%prtcl_dt
    
    nUnitFile = d_tec_geom_unit
    open( unit =  nUnitFile , file = chfile )
    
    strnd_id = 2     
    call TecWall%TecPlot_Out(nUnitFile , sol_time, strnd_id )
    
    close(nUnitFile) 
    
    
end subroutine 


!**********************************************************************
! writing rawData bed condition of particles in an external file
!**********************************************************************
subroutine DEMS_WriteRawData(this, num)
    implicit none
    class(DEMSystem) this
    integer(IK),intent(in):: num
    
    !// locals 
    integer(IK):: nUnit = d_wrt_prtc_unit
    integer(IK) i, num_InDomain
    character(256) chFile, ch

    !// body
    
    ! openning file
    write(ch,*) this%m_DEM_Options%base_Number + num
    chFile = trim(this%m_DEM_Options%Res_Dir)//'RawData'//trim(adjustl(ch))//"for"//trim(this%m_DEM_Options%RunName)//".dat"
    open( file = chFile, unit = nUnit )
    
    num_InDomain = this%numInDomain()
    
    ! first writing number of particles
    write( nUnit , * ) "numberOfParticles" , num_inDomain   
    write( nUnit , * ) "id	type	posision	velocity	mark"
    ! writing particle data 
    do i=1,this%numPrtcl
        
        if( this%prtcl_flag(i) >= Pflg_inDomain ) then
            
            write( nUnit , * ) this%prtcl_ids(i) , &
                               this%prtcl_type(i), &
                               this%prtcl_dpos(i), &
                               this%prtcl_vel(i) , &
                               this%prtcl_usr_mark(i)
        end if
        
    end do
    
    close( nUnit ) 

end subroutine 

!**********************************************************************
! writing particle data to a VTK file
!**********************************************************************

subroutine DEMS_PrtclToVTK(this, file_num)
    implicit none
    class(DEMSystem) this
    integer(IK),intent(in):: file_num

    !// locals
    integer(IK) i, j, nw, nUnitFile
    integer(IK) numPrtcl, num_inDomain
    type(real3) pos
    character(256) ch, chFile, chMsg
    
    ! file name for VTK file
    write(ch,*) this%m_DEM_Options%base_Number + file_num
    chFile = trim(this%m_DEM_Options%Res_Dir)//"ParticlesFor"//trim(this%m_DEM_Options%RunName)//trim(adjustl(ch))//".vtk"
    
    ! Opens VTK file
    nUnitFile = d_tec_prtcl_unit
    open( unit =  nUnitFile , file = chfile )
    
    ! Header for VTK file
    chMsg = "Particles for run name: " // trim(this%m_DEM_Options%RunName)   
    call vtk_header( nUnitFile, chMsg )
    
    ! counting number of particles in the domain
    numPrtcl = this%numPrtcl
    num_inDomain = this%numInDomain() 
    
    write(nUnitFile, "(A)")"DATASET UNSTRUCTURED_GRID"
    write(nUnitFile, "(A,I7,A)" )"POINTS ", num_inDomain , " float"
    
    do i= 1, numPrtcl
        if( this%prtcl_flag(i) >= Pflg_inDomain ) then
            pos = this%prtcl_dPos(i)
            write(nUnitFile,*) pos
            
        end if
    end do
    
    ! creating cell types. Cell types are finite element point. Each
    ! point represents the center of the particle. 
    write( nUnitFile, "(A,I8,I8)")"CELLS", num_inDomain, 2*num_inDomain
    do i=1,num_inDomain
        write(nUnitFile, *) 1, i-1    
    end do
    
    write(nUnitFile, "(A,I8)")"CELL_TYPES ", num_inDomain
    do i=1, num_inDomain
        write(nUnitFile,*) 1
    end do
    
    write(nUnitFile, "(A,I8)") "POINT_DATA ", num_inDomain
	write(nUnitFile, "(A)" ) "SCALARS diam float 1 "
	write(nUnitFile, "(A)" ) "LOOKUP_TABLE default "
	
    do i=1,numPrtcl
        if( this%prtcl_flag(i) >= Pflg_inDomain ) write(nUnitFile, *) this%prtcl_dpos(i)%w
    end do
    
	
    !// id of particles 
    write(nUnitFile, "(A)" ) "SCALARS ID int 1 "
	write(nUnitFile, "(A)" ) "LOOKUP_TABLE default "
	
	do i=1,numPrtcl
        if( this%prtcl_flag(i) >= Pflg_inDomain ) write(nUnitFile, *) this%prtcl_ids(i)
    end do
    
    
    !// user mark of particles 
    write(nUnitFile, "(A)" ) "SCALARS usr_mark int 1 "
	write(nUnitFile, "(A)" ) "LOOKUP_TABLE default "
	
	do i=1,numPrtcl
        if( this%prtcl_flag(i) >= Pflg_inDomain ) write(nUnitFile, *) this%prtcl_usr_mark(i)
    end do
    
    
    !// field data -> velocity
    write(nUnitFile, "(A)") "FIELD FieldData 3"
    write(nUnitFile, "(A,I8,A)") "Velocity 3 ", num_inDomain, " float "
    do i=1,numPrtcl
        if( this%prtcl_flag(i) >= Pflg_inDomain ) write(nUnitFile, *) this%prtcl_vel(i)
    end do
    
 
    write(nUnitFile, "(A,I8,A)") "Omega 3 ", num_inDomain, " float "
    do i=1,numPrtcl
        if( this%prtcl_flag(i) >= Pflg_inDomain ) write(nUnitFile, *) this%prtcl_rvel(i)
    end do
    
    write(nUnitFile, "(A,I8,A)") "Fc 3 ", num_inDomain, " float "
    do i=1,numPrtcl
        if( this%prtcl_flag(i) >= Pflg_inDomain ) write(nUnitFile, *) this%prtcl_cntct_force(i)
    end do
    
    
    
end subroutine

!****
subroutine DEMS_outputFile(this )
    implicit none
    class(DEMSystem) this
    
    

    if( this%IterNumber == 1 .or. mod(this%IterNumber, this%m_DEM_Options%SaveFreq) == 0 ) then
              
        if( iAND(this%m_DEM_options%OutputFileType,OP_Type_VTK) == OP_Type_VTK ) then
            call this%m_Geometry%GeomToVTK(this%iterNumber, this%m_DEM_Options )
            call this%PrtclToVTK(this%iterNumber)
        end if
               
        if( iAND(this%m_DEM_options%OutputFileType,OP_Type_Tec) == OP_Type_Tec ) then
            call this%write_geom(this%iterNumber)
            call this%write_prtcl(this%iterNumber)
        end if
                
        if( iAND(this%m_DEM_options%OutputFileType,OP_Type_Bin) == OP_Type_Bin ) then
            !call this%write_geom(this%iterNumber)
            !call this%write_prtcl(this%iterNumber)
            print*, "*********** Binary file is not specified *****************"
        end if
        
        if( iAND(this%m_DEM_options%OutputFileType,OP_Type_Bin) == OP_Type_Raw ) then
            call this%writeRawData(this%iterNumber)
        end if
            
    endif

end subroutine 


!*******
subroutine DEMS_terminalOutput(this)
    implicit none
    class(DEMSystem) this
    
    !// locals
    integer(IK) nto
    
    ! output to log file and terminal/command window 
    if( this%IterNumber == 1 .or. mod(this%IterNumber, this%m_DEM_Options%Cmd_LFile_Freq) == 0 ) then
            
        ! command window and log file output
        nto = this%IterNumber 
        call MainLogInfo%OutInfo("Program performed " //trim(num2str(nto)) // " iterations up to here!" , 1 )
            
        call this%DEMS_timers_output()
        call this%DEMS_report_contacts()
    endif  
        
end subroutine 

!**********************************************************************
! writing packed bed condition of particles in an external file
!**********************************************************************
subroutine DEMS_WritePackedBed(this, chFile)
    implicit none
    class(DEMSystem) this
    character(*),intent(in):: chFile ! file name
    
    !// locals 
    integer(IK):: nUnit = d_wrt_prtc_unit
    integer(IK) i, num_InDomain, numProps
    class(base_PhysicalProperty),pointer:: phProp
    !// body
    
    ! openning file
    open( file = chFile, unit = nUnit )
    
    num_InDomain = this%numInDomain()
    
    ! first writing number of particles
    write( nUnit , * ) num_inDomain
    
    ! writing number of property tyeps 
    numProps = this%m_PhysProp%get_numPrtclType()
    write( nUnit , * ) numProps
    
    do i= 1, numProps
        select type (phProp => this%m_PhysProp )
            type is (PhysicalProperty)
            write( nUnit, * ) phProp%Prtcl_PureProp(i)    
            
        end select
            
    end do
    
    
    
    
    ! writing particle data 
    do i=1,this%numPrtcl
        
        if( this%prtcl_flag(i) >= Pflg_inDomain ) then
            
            write( nUnit , * ) this%prtcl_ids(i) , &
                               this%prtcl_type(i), &
                               this%prtcl_dpos(i), &
                               this%prtcl_vel(i)
        end if
        
    end do
    
    close( nUnit ) 

end subroutine 

!**********************************************************************
! calculating acceleration of particles (linear and angular) 
!**********************************************************************
subroutine DEMS_clc_Acceleration(this)
    implicit none
    class(DEMSystem) this
    
    !// locals
    integer(IK) i, pt, numPrtcl
    real(RK)    mass, inertia
    type(real3) gravity
    
    
    numPrtcl = this%numPrtcl
    gravity = this%m_DEM_Options%gravity
    
    do i = 1, numPrtcl
        
        if(this%prtcl_flag(i) >= Pflg_inDomain  )then
            
            mass    = this%m_PhysProp%get_pure_mass( this%Prtcl_Type(i) )
            inertia = this%m_PhysProp%get_pure_inertia( this%Prtcl_type(i) ) 
            
            this%prtcl_Lin_acc(i) = ( ( this%prtcl_cntct_force(i) )/mass ) + gravity
            this%prtcl_rot_acc(i) = ( ( this%prtcl_torque(i)/inertia))
             
        end if
        
    end do
    
    
end subroutine


!**********************************************************************
! setting integration method 
!**********************************************************************
subroutine DEMS_setIntegration( this )
    implicit none
    class(DEMSystem) this
    
    
    integer(IK) PI_Method, PRI_Method
    
    
    ! getting integration method
    PI_Method = this%m_DEM_Options%PI_Method
    PRI_Method = this%m_DEM_Options%PRI_Method

    ! linear integration
    if( associated(this%Lin_Integration) ) deallocate(this%Lin_Integration)
    allocate( this%Lin_Integration )

    call this%Lin_Integration%Init_Integration( PI_Method, this%max_nPrtcl, this%numPrtcl, this%prtcl_dt,  &
                                               this%prtcl_dPos, this%prtcl_vel,  this%prtcl_Lin_acc , this%prtcl_flag )
                                           
    
    ! rotational integration
    if( associated(this%rot_Integration) ) deallocate(this%rot_Integration)
    allocate( this%rot_Integration )

    call this%rot_Integration%Init_rot_Integration( PRI_Method, this%max_nPrtcl, this%numPrtcl, this%prtcl_dt,&
                                               this%prtcl_teta, this%prtcl_rvel,  this%prtcl_rot_acc, this%prtcl_flag )
                                           
    
    this%Integration_init = .true.
    
end subroutine

!**********************************************************************
! setting contact search method 
!**********************************************************************
subroutine DEMS_setContactSearch( this, minDomain, maxDomain  )
    implicit none
    class(DEMSystem)         this
    type(real3),intent(in):: minDomain, maxDomain
    
    !// locals
    integer(IK) lnumlvl, CS_Method
    real(RK)    lratio

    if( (.not.this%Prtcls_init)  .or.  (.not.this%Geometry_init) .or. &
        (.not.this%PhysProp_init) .or. (.not.this%ContForce_init) ) then
        stop "Previous components of contact search are not initialized yet!"    
    end if

    
    if(associated( this%Cont_Search) ) deallocate(this%Cont_Search)
    allocate( this%Cont_Search )
    
    
    
    ! Particle-particle contact search 
    call MainLogInfo%OutInfo("Particle-Particle contact search intialization...",2)
    
    ! getting required info for initialization
    lnumlvl     = this%m_DEM_Options%CS_numlvls 
    CS_Method   = this%m_DEM_Options%CS_Method
    lratio      = this%m_DEM_Options%CS_ratio
    
    
    if( lnumlvl /= 0 ) then
        
        ! initializing particle-particle contact search with pre-defined number of levels   
        call this%Cont_Search%InitContactSearch(CS_Method, minDomain, maxDomain, &
                                                this%max_nPrtcl, this%numPrtcl, this%prtcl_ids, this%Prtcl_flag , &
                                                this%prtcl_dPos, this%PP_cont_List, ratio = lratio  ,num_lvls = lnumlvl )
    
        
    else
        ! initializing particle-particle contact search with default number of levels  
        call this%Cont_Search%InitContactSearch(CS_Method, minDomain, maxDomain, &
                                                this%max_nPrtcl, this%numPrtcl, this%prtcl_ids, this%Prtcl_flag, &
                                                this%prtcl_dPos, this%PP_cont_List , ratio = lratio )
        
    end if
    
    
    ! initializing particle-wall contact search   
    call MainLogInfo%OutInfo("Particle-Wall contact search intialization...",2)
    
    if(associated( this%Cont_Search_PW) ) deallocate(this%Cont_Search_PW)
    allocate( this%Cont_Search_PW )
    
    call this%Cont_Search_PW%InitContactSearch( this%PW_cont_list, this%m_Geometry )
    
    
    this%ContSearch_Init = .true.
    
end subroutine

!**********************************************************************
!   setting contact force object 
!**********************************************************************
subroutine DEMS_SetContactForce( this )
    implicit none
    class(DEMSystem) this
    
    !// locals 
    integer(IK)     CF_Type, CT_Model
    integer(IK)     max_pp, max_pw
    character(62)   ch
        
    ! getting type of contact force model
    CF_Type  = this%m_DEM_Options%CF_Type
    CT_Model = this%m_DEM_Options%CT_Model 
    
    ! first checking whether all required object are initialized or not!
    if( .not. this%Prtcls_init .or. .not. this%Geometry_init .or. .not.this%PhysProp_init ) then
        stop "Either particls or geometry or physical properties are not set yet"    
    end if
    
    
    ! The contact force type 
    if( CF_Type == CFT_LSD_nl .or. CF_Type == CFT_LSD_l .or. &
        CF_Type == CFT_nLin_l .or. CF_Type == CFT_nLin_nl ) then
        
    	
        if( CF_Type == CFT_LSD_nl .or. CF_Type == CFT_LSD_l ) then
            allocate( LSD_ContactForce::this%m_Cont_force )
        else
            allocate( nonLin_ContactForce::this%m_Cont_force )
        end if
        
        

        if( CF_Type == CFT_LSD_nl ) then
            ch = "linear spring-dashpot with non-limited tangential displacement"
        elseif( CF_Type == CFT_LSD_l ) then
            ch = "linear spring-dashpot with limited tangential displacement"
        elseif( CF_Type == CFT_nLin_nl ) then
            ch = "non-linear visco-elastic model with non-limited tangential displacement"
        else
            ch = "non-linear visco-elastic model with limited tangential displacement"
        end if
        
        call MainLogInfo%OutInfo("Contact force model is "//trim(ch), 2 )
        
        call this%m_Cont_force%InitContactForce( CF_Type, CT_Model, &
                                this%max_nPrtcl, this%numPrtcl, this%prtcl_dt, &
                               this%prtcl_type, this%prtcl_dPos, this%prtcl_vel , &
                               this%prtcl_rvel, this%prtcl_cntct_force, this%prtcl_torque )
        
        
	
		
        
	! initializing particle-particle contact list
        call MainLogInfo%OutInfo("Initializing particle-particle contact list object:",2)        
      	
            allocate( ContactList:: this%PP_Cont_List )
	
            max_pp = int(PP_cntctList_Size * this%max_nPrtcl) + 1
            call this%PP_Cont_list%InitContactList( max_pp , this%max_nPrtcl )
	
        
        call MainLogInfo%OutInfo("The container size is : "//trim(num2str(max_pp)),3)
        
        
        ! initializing particle-wall contact list
        call MainLogInfo%OutInfo("Initializing particle-wall contact list object:",2)        
        
            allocate( ContactList:: this%PW_Cont_List )
            max_pw = int(Pw_cntctList_Size * this%max_nPrtcl) + 1
            call this%PW_Cont_list%InitContactList( max_pw , this%max_nPrtcl )
            
        call MainLogInfo%OutInfo("The container size is : "//trim(num2str(max_pw)),3)
            
               
        ! setting other required components for contact force object 
        call this%m_Cont_force%set_ContactLists(this%PP_Cont_list,this%PW_Cont_list )
        call this%m_cont_force%set_PysicalProperty(this%m_PhysProp)
        call this%m_cont_force%set_Geometry(this%m_Geometry)
        
    else
        stop "Wrong contact force type in input"
    end if
    
    
    this%ContForce_init = .true.
    this%ContList_init = .true.
    
end subroutine 

!**********************************************************************
! setting geometry object 
!**********************************************************************
subroutine DEMS_SetGeometry(this, geom )

    implicit none
    class(DEMSystem) this
    class(Geometry),pointer:: geom
    
    if( .not. this%Prtcls_init ) then
        stop "Particles are not initialized yet!"    
    end if
    
    this%m_Geometry => geom
    this%Geometry_init = .true.
    
end subroutine 

!**********************************************************************
! returning number of particles in the simulation domain
!**********************************************************************
function DEMS_numInDomain(this) result(res)
    implicit none
    class(DEMSystem) this
    integer(IK) res
    
    !// locals 
    integer(IK) i
    
    res = 0
    
    do i = 1, this%numPrtcl
        if( this%prtcl_flag(i) == Pflg_inDomain )  res = res + 1  
	
    end do
    
    
end function



subroutine DEMS_User_prtclMark( this )
    implicit none
    class(DEMSystem) this
    
    integer(IK) i
    
    do i=1, this%max_nPrtcl
        
        this%prtcl_usr_mark(i) = User_Mark  &
        			( this%prtcl_ids(i), 	&
        			  this%prtcl_flag(i), 	&
        			  this%prtcl_type(i), 	&
        			  this%prtcl_dpos(i),	&
        			  this%prtcl_usr_mark(i))
        
    end do
    
end subroutine 


subroutine DEMS_WriteParticleData( this, chFile )
    
    implicit none
    class(DEMSystem) this
    character(*),intent(in) :: chFile 
    
    !// locals
    integer(IK) numIn, i
    
    open( unit = d_wrt_prtc_unit , file = chFile ) 
    
    numIn = this%numInDomain();
    
    write( d_wrt_prtc_unit , * ) numIn
    
    do i=1, numIn
       if( this%prtcl_flag(i) == Pflg_inDomain ) then
           write( d_wrt_prtc_unit , * ) this%prtcl_ids(i), this%prtcl_dPos(i), this%prtcl_vel(i)
       end if
    end do
    
    
end subroutine 

!***********************************************************************
!* final procedure for DEMsystem
!***********************************************************************
subroutine Finilize_DEMSystem(this)
    
    implicit none
    type(DEMSystem) this
    
    call this%DEMS_deallocateAll()
    
end subroutine


!***********************************************************************
!resetting all timers in DEMsystem
!***********************************************************************
subroutine DEMS_timers_reset( this )
    implicit none
    class(DEMSystem) this
    
    
    call this%m_total_timer%reset()
    call this%m_pre_iter_timer%reset()
    call this%m_prediction_timer%reset()
    call this%m_ContSearchPP_timer%reset()
    call this%m_ContSearchPW_timer%reset()
    call this%m_ForcePP_timer%reset()
    call this%m_ForcePW_timer%reset()
    call this%m_Acceleration_timer%reset()
    call this%m_integration_timer%reset()

end subroutine 


subroutine DEMS_timers_output(this)
    implicit none
    class(DEMSystem) this
    
    call MainLogInfo%OutInfo( "Execution time [tot, last, ave] [sec]: "// trim( num2str(this%m_total_timer%total() ) )// &
                          ", "// trim(num2str(this%m_total_timer%last() ))//", "//trim(num2str(this%m_total_timer%average() )) , 1 )
    
    call MainLogInfo%OutInfo( "PreItertion time [tot, ave]       : "//trim(num2str(this%m_pre_iter_timer%total()))//", "// &
                              trim(num2str(this%m_pre_iter_timer%average())), 3 )
    
    call MainLogInfo%OutInfo( "Prediction time [tot, ave]        : "//trim(num2str(this%m_prediction_timer%total()))//", "// &
                               trim(num2str(this%m_prediction_timer%average())), 3 )
    call MainLogInfo%OutInfo( "Contact search P-P time [tot, ave]: "//trim(num2str(this%m_ContSearchPP_timer%total()))//", "// &
                               trim(num2str(this%m_ContSearchPP_timer%average())), 3 )
    call MainLogInfo%OutInfo( "Contact search P-W time [tot, ave]: "//trim(num2str(this%m_ContSearchPW_timer%total()))//", "// &
                               trim(num2str(this%m_ContSearchPW_timer%average())), 3 )
    call MainLogInfo%OutInfo( "Contact force P-P time [tot, ave] : "//trim(num2str(this%m_ForcePP_timer%total()))//", "// &
                               trim(num2str(this%m_ForcePP_timer%average())), 3 )
    call MainLogInfo%OutInfo( "Contact force P-W time [tot, ave] : "//trim(num2str(this%m_ForcePW_timer%total()))//", "// &
                               trim(num2str(this%m_ForcePW_timer%average())), 3 )
    call MainLogInfo%OutInfo( "Acceleration time [tot, ave]      : "//trim(num2str(this%m_Acceleration_timer%total()))//", "// &
                              trim(num2str(this%m_Acceleration_timer%average())), 3 )
    call MainLogInfo%OutInfo( "Integration time [tot, ave]       : "//trim(num2str(this%m_integration_timer%total()))//", "// &
                              trim(num2str(this%m_integration_timer%average())), 3 )
    call MainLogInfo%OutInfo( "Write to file time [tot, ave]     : "//trim(num2str(this%m_write_prtcl_timer%total()))//", "// &
                              trim(num2str(this%m_write_prtcl_timer%average())), 3 )
    
    
end subroutine


subroutine DEM_GetParticles( this , Particles )
    implicit none
    class(DEMSystem)    this
    class(PSD_Property),intent(in):: Particles
    
    
    !// locals
    integer(IK) i, max_Prtcl
    
    max_Prtcl = this%max_nPrtcl
    
    do i=1, max_Prtcl
        
        this%prtcl_ids(i)  = Particles%get_ids(i)
        this%prtcl_type(i) = Particles%get_type(i)
        this%prtcl_dpos(i)  = Particles%get_dpos(i)
        this%prtcl_teta(i)  = zero_r3
        this%prtcl_rvel(i)  = zero_r3
        
	        
        select type( Particles )
        type is(PSD_Property)
            this%prtcl_vel(i)   = zero_r3
        type is(PSDP_Position)
            this%prtcl_vel(i)   = Particles%get_vel(i)
        
        end select
        
    end do
    

    this%prtcl_cntct_force = zero_r3
    this%prtcl_torque = zero_r3
    this%prtcl_Lin_acc = zero_r3
    this%prtcl_rot_acc = zero_r3
    
    
end subroutine


subroutine DEM_reallocateParticles(this)
    implicit none
    class(DEMSystem) this
    
    
    ! deallocting all varaibales 
    call this%DEMS_deallocateAll()
        
    allocate( this%prtcl_ids(this%max_nPrtcl)  )
    allocate( this%prtcl_flag(this%max_nPrtcl) )
    allocate( this%prtcl_type(this%max_nPrtcl) )
    allocate( this%prtcl_usr_mark(this%max_nPrtcl) )
    allocate( this%prtcl_dPos(this%max_nPrtcl) )
    allocate( this%prtcl_vel (this%max_nPrtcl) )
    allocate( this% prtcl_teta(this%max_nPrtcl) )
    allocate( this%prtcl_rvel( this%max_nPrtcl) )
    allocate( this%prtcl_cntct_force( this%max_nPrtcl) )
    allocate( this%prtcl_torque( this%max_nPrtcl) )
    allocate( this%prtcl_Lin_acc( this%max_nPrtcl) )
    allocate( this%prtcl_rot_acc( this%max_nPrtcl) )
    
end subroutine


subroutine DEMS_deallocateAll(this)
    implicit none
    class(DEMSystem) this
    
    if( associated(this%prtcl_ids)  ) deallocate(this%prtcl_ids)
    if( associated(this%prtcl_flag) ) deallocate(this%prtcl_flag)
    if( associated(this%prtcl_type) ) deallocate(this%prtcl_type)
    if( associated(this%prtcl_usr_mark) ) deallocate(this%prtcl_usr_mark)
    if( associated(this%prtcl_dpos) ) deallocate(this%prtcl_dpos)
    if( associated(this%prtcl_vel)  ) deallocate(this%prtcl_vel)
    if( associated(this%prtcl_teta) ) deallocate(this%prtcl_teta)
    if( associated(this%prtcl_rvel) ) deallocate(this%prtcl_rvel)
    if( associated(this%prtcl_cntct_force) ) deallocate(this%prtcl_cntct_force)
    if( associated(this%prtcl_torque) ) deallocate(this%prtcl_torque)
    if( associated(this%prtcl_Lin_acc)) deallocate(this%prtcl_Lin_acc)
    if( associated(this%prtcl_rot_acc)) deallocate(this%prtcl_rot_acc)
    
end subroutine 



subroutine DEMS_report_contacts( this )
    implicit none
    class(DEMSystem) this
    
    !// locals 
    integer(IK):: PP_Cont, PW_Cont, Consv_Cont(2)
    character(256):: chLine
    
    ! exact contacts 
    PP_Cont = this%PP_Cont_List%getNumCntcts()
    PW_Cont = this%PW_Cont_List%getNumCntcts()
    
    ! conservative contacts 
    Consv_Cont = this%Cont_Search%get_numContact()
    
    call MainLogInfo%OutInfo("Contact information", 2)
    chLine = "No. consrvtv. contacts, same level | cross level: "// &
             trim(num2str(Consv_Cont(1)))//" | "//trim(num2str(Consv_Cont(2)))
    call MainLogInfo%OutInfo(chLine, 3)
    
    chLine = "No. exact contacts P-P | P-W : "// trim( num2str(PP_Cont) ) //" | "//trim(num2str(PW_Cont))
    call MainLogInfo%OutInfo(chLine, 3)
end subroutine




end module


