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
!  file name  : g_Prtcl_DefaultValues.f90 
!  module name: g_Prtcl_DefaultValues
!  
!  Purpose: 
!   Default values and settings for the program 
! 
!------------------------------------------------------------------------------
    
module g_Prtcl_DefaultValues
    
    use g_TypeDef
    use g_LogInfo
    
    implicit none
    
    
   
    integer(IK),parameter :: CFT_LSD_nl   = 10_IK
    integer(IK),parameter :: CFT_LSD_l    = 11_IK
    integer(IK),parameter :: CFT_nLin_nl  = 12_IK
    integer(IK),parameter :: CFT_nLin_l   = 13_IK
    
    integer(IK),parameter :: CTM_ConstantTorque = 10_IK
    integer(IK),parameter :: CTM_VariableTorque = 11_IK
    
    
    
    integer(IK),parameter:: CSM_NBS                 = 10_IK
    integer(IK),parameter:: CSM_NBS_Munjiza         = 11_IK
    integer(IK),parameter:: CSM_NBS_Hrchl           = 12_IK
    integer(IK),parameter:: CSM_NBS_Munjiza_Hrchl   = 13_IK
    integer(IK),parameter:: CSM_DESS                = 20_IK
    
    
    integer(IK),parameter:: PIM_FE          = 10_IK
    integer(IK),parameter:: PIM_ME          = 11_IK
    integer(IK),parameter:: PIM_Taylor2     = 12_IK
    integer(IK),parameter:: PIM_Taylor3     = 13_IK
    integer(IK),parameter:: PIM_Taylor4     = 14_IK
    integer(IK),parameter:: PIM_AB2         = 20_IK
    integer(IK),parameter:: PIM_AB3         = 21_IK
    integer(IK),parameter:: PIM_AB4         = 22_IK
    integer(IK),parameter:: PIM_AB5         = 23_IK
    integer(IK),parameter:: PIM_AB2AM3      = 30_IK
    integer(IK),parameter:: PIM_AB3AM4      = 31_IK
    integer(IK),parameter:: PIM_AB4AM5      = 32_IK
    integer(IK),parameter:: PIM_Gear3       = 33_IK
    integer(IK),parameter:: PIM_Gear4       = 34_IK
    integer(IK),parameter:: PIM_Gear5       = 35_IK
    
    character(26),parameter,dimension(15) ::chPIM = (/"Forward Euler             ",      &
                                                      "Modified Euler            " ,    &
                                                      "Taylor: 2nd Order         ",  &
                                                      "Taylor: 3rd Order         ",  &
                                                      "Taylor: 4th Order         ",  &
                                                      "Adams Bashforth: 2nd Order", &
                                                      "Adams Bashforth: 3rd Order", &
                                                      "Adams Bashforth: 4th Order", &
                                                      "Adams Bashforth: 5th Order", &
                                                      "AB-Moultun: 3rd Order     ",      &
                                                      "AB-Moultun: 4th Order     ",      &
                                                      "AB-Moultun: 5th Order     ",      &
                                                      "Gear: 3rd Order           ",            &
                                                      "Gear: 4th Order           ",            &
                                                      "Gear: 5th Order           " /)
    
    
    
    integer(IK),parameter:: PInsM_File      = 10_IK
    integer(IK),parameter:: PInsM_Box       = 11_IK
    
    integer(IK),parameter:: PSD_Uniform     = 10_IK
    integer(IK),parameter:: PSD_Normal      = 11_IK
    integer(IK),parameter:: PSD_LogNormal   = 12_IK
   
    
    
    integer(IK),parameter:: Pflg_inDomain   =  0_IK
    integer(IK),parameter:: pflg_deleted    = -2_IK
    integer(IK),parameter:: pflg_notInserted= -1_IK
    
    integer(IK),parameter:: x_axis          = 1_IK
    integer(IK),parameter:: y_axis          = 2_IK
    integer(IK),parameter:: z_axis          = 3_IK
    
    
    integer(IK),parameter:: OP_Type_VTK    = 1_IK
    integer(IK),parameter:: OP_Type_Tec    = 2_IK
    integer(IK),parameter:: OP_Type_Bin    = 4_IK
    integer(IK),parameter:: OP_Type_Raw    = 8_IK
    
    
    ! default values 
    
    real(RK),parameter:: d_Young_mod        = 1.0E6_RK
    real(RK),parameter:: d_Shear_mod        = 4.0E5_RK
    real(RK),parameter:: d_Poissons_ratio   = 0.25_RK
    real(RK),parameter:: d_kappa            = (1-d_Poissons_ratio)/(1-0.5*d_Poissons_ratio)
    real(RK),parameter:: d_Density          = 2500.0_RK
    real(RK),parameter:: d_YeildStress      = 5.0E5_RK
    real(RK),parameter:: d_rel_vel          = 1.0_RK
    
    real(RK),parameter:: d_fric             = 0.3_RK
    real(RK),parameter:: d_roll_fric        = 0.1_RK
    real(RK),parameter:: d_en               = 0.75_RK
    real(RK),parameter:: d_et               = 1.0_RK
    
    ! contact Force
    integer(IK),parameter:: d_CF_type  = CFT_LSD_nl
    integer(IK),parameter:: d_CT_Model = CTM_ConstantTorque 
    
    ! contact search
    integer(IK),parameter:: d_CS_Method = CSM_NBS_Munjiza
    
    ! integrtaion method
    integer(IK),parameter:: d_PI_Method = PIM_FE
    integer(IK),parameter:: d_PRI_Method =PIM_FE
    
    ! Contact List
    real(IK),parameter:: PP_cntctList_Size = 5.0
    real(IK),parameter:: PW_cntctList_Size = 3.0
    
    
    real(RK),parameter:: d_prtcl_dt = 1.0E-5_RK
    
    type(real3),parameter:: d_gravity_acc = real3(0.0_RK, -9.8_RK, 0.0_RK)
    
    ! properties
    
    
    integer(IK),parameter:: d_Base_wall_id    = 10000000 ! this preserve 10000000 particles to be simulated in the program
                                                         ! change this if number of particles is greater than 10000000. 
    integer(IK),parameter:: d_Base_Wall_id_cyl= 10000
    
    
    integer(IK),parameter:: d_Wall_max_update_iter= 100
    real(RK),parameter::    d_Wall_neighbor_ratio = 2.0_RK
    
    integer(IK),parameter:: d_PIns_Mthd = PInsM_File
    
    ! file units for diffrent routines
    integer(IK),parameter:: d_MakePrtcls_unit  = 10001
    integer(IK),parameter:: d_ReadPrtcls_unit  = 10002
    integer(IK),parameter:: d_MainLogInfo_unit = 10003
    integer(IK),parameter:: d_stl_reader_unit   = 10004
    integer(IK),parameter:: d_tec_prtcl_unit    = 10005
    integer(IK),parameter:: d_tec_geom_unit     = 10006
    integer(IK),parameter:: d_wrt_prtc_unit     = 10007
    
    
    type(LogInfo) :: MainLogInfo
    character(*),parameter:: d_RunName = "TestRun"
    
    ! Is used for calculating wall velocity
    real(RK),parameter:: d_dt_wallVel = 0.0000001_RK
    
    
    type DEMS_Options
        character(256)::RunName           = d_RunName      ! run name
        character(256)::Res_Dir           = "."            ! result directory 
        type(real3)::   gravity           = d_gravity_acc  ! gravity
        integer(IK)::   CF_Type           = d_CF_Type      ! contact force type
        integer(IK)::   CT_Model          = d_CT_Model     ! contact torque type
        integer(IK)::   CS_Method         = d_CS_Method    ! contact search method 
        integer(IK)::   CS_numlvls        = 0              ! means default behavior, number of levels in multi-level contact search
        real(RK)   ::   CS_ratio          = 1.0_RK         ! particle to cell size ratio for cell-based methods 
        integer(IK)::   PI_Method         = d_PI_Method    ! integration scheme for translational motion
        integer(IK)::   PRI_Method        = d_PRI_Method   ! integration scheme for rotational motion
        integer(IK)::   SaveFreq          = 1000_IK        ! save frequency for output file 
        integer(IK)::   OutputFileType    = OP_Type_VTK    ! save results as vtk file
        integer(IK)::   Cmd_LFile_Freq    = 500_IK         ! report frequency in the terminal 
        integer(IK)::   base_number       = 1000000_IK     ! base number for generating numbered output files 
        real(RK)   ::   rel_vel_kn        = d_rel_vel      ! relative particle velocity for calculating spring stiffness
        integer(IK)::   LF_file_lvl       = 5              ! logfile report level
        integer(IK)::   LF_cmdw_lvl       = 3              ! terminal report level 
    end type
    
    
end module
