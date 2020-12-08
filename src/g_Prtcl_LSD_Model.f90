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
!  file name  : g_Prtcl_LSD_Model.f90
!  module name: g_Prtcl_LSD_Model
! 
!  Purpose: 
!    1) To extend the base class base_ContactForce to calculate particle-particle
!       and particle-wall contact forces based on the Linear Spring-Dashpot (LSD)
!       model
! 
!   This class works in conjunction with ContactList class where the history of
!     all contact pairs are stored. ContactList class stores the tangential
!     overlap and tangential vector of contact pairs. If you are willing to
!     develop models that requires you to store more data of contact pair, you
!     should do the followings:
!      1) Extend the class base_ContactInfo or ContactInfo1 and provide methods
!         if required.
!      2) Extend the base_ContactList or ContactList and provide/modify methods
!         for it.
!      3) Extend the class base_ContactForce or LSD_ContactForce classes and 
!         develop your model there.
!                                  
!------------------------------------------------------------------------------
    module g_Prtcl_LSD_Model

    use g_Prtcl_ContactForce
    implicit none
    
    type, extends(base_ContactForce):: LSD_ContactForce
        integer,allocatable::dummy
        type(PhysicalProperty),pointer:: Ph_Props
        type(ContactList),pointer:: PP_CntctList
        type(ContactList),pointer:: PW_CntctList
        
    contains
    
    
    ! overriding of deferred bindings 

    procedure:: ContactForce_PP     => LSD_ContactForce_PP
    procedure:: ContactForce_PW     => LSD_ContactForce_PW
    procedure:: set_ContactLists    => LSD_set_ContactLists
    procedure:: set_PysicalProperty => LSD_set_PysicalProperty
    
    procedure:: clc_ContactVelocity     => LSD_clc_ContactVelocity
    procedure:: clc_ContactVelocity_PW  => LSD_clc_ContactVelocity_PW
    procedure:: clc_Torque              => LSD_clc_Torque
    
    end type
    
contains

!**********************************************************************    
! setting physical property object
!**********************************************************************
subroutine LSD_set_PysicalProperty(this, PhProps )
    implicit none
    
    class(LSD_ContactForce) this
    class(base_PhysicalProperty),pointer,intent(in) :: PhProps
    
    select type( PhProps )
    class is (PhysicalProperty)
        this%Ph_Props => PhProps
        
    end select
    
end subroutine

!**********************************************************************
! setting particle-particle and particle-wall contact lists 
!**********************************************************************
subroutine LSD_set_ContactLists(this , PP_CntctList, PW_CntctList)

    class(LSD_ContactForce) this
    class(base_ContactList), pointer, intent(in) :: PP_CntctList, PW_CntctList
    
    select type (PP_CntctList)
    class is( ContactList )
        
        this%PP_CntctList => PP_CntctList
        this%b_PP_ContList    => PP_CntctList
        
        class default
            call CheckForError( ErrT_Abort, "LSD_set_ContactLists" , "wrong type of contact list in the input" )
            
    end select
    
    select type (PW_CntctList)
    class is(ContactList)    
        this%PW_CntctList => PW_CntctList
        this%b_PW_ContList=> PW_CntctList 
        
        class default
            call CheckForError( ErrT_Abort, "LSD_set_ContactLists" , "wrong type of contact list in the input" )
            
    end select
    
            
end subroutine

!**********************************************************************
!   calculating the contact force linear spring-dashpot model and 
! tangential and rolling torques. (Particles-Particle) 
!**********************************************************************
subroutine LSD_ContactForce_PP(this , ind , lnew )
    implicit none
    class(LSD_ContactForce) this
    integer(IK),intent(in):: ind
    logical,intent(in)    :: lnew
    
    !// locals
    integer(IK) i, j
    integer(IK) pt_i , pt_j
    real(RK)    Ri, Rj, fn, ft , ft_f
    real(RK)    vrn, ovrlp, ft_fric
    type(real3) rveli, rvelj
    type(real3) norm_v
    type(real3) Vrij, Vij_n, Vij_t, ovlp_t
    type(real3) fnij, ftij, Mij, Mji, Mri, Mrj 
    type(ContactInfo1) ContPair_i
    type(BinaryProperty) Prop_ij
    
    !// body   
    ! getting the contacting pairs from contact list
    ContPair_i = this%PP_CntctList%getItem(ind)
        
    ! particle indeces 
    i = ContPair_i%pari
    j = ContPair_i%parj
    
    ! type of particle to get physical properties
    pt_i = this%get_pType(i)
    pt_j = this%get_pType(j)
    ! getting the binary properties of particles i and j 
    prop_ij = this%Ph_Props%Bnry_getProp( pt_i, pt_j)
    
    
    call this%clc_ContactVelocity(i,j, Vrij, rveli, rvelj, norm_v, ovrlp , Ri, Rj )
    
    if( lnew) then
        ! if this is a new contact, the tangential overlap should be zero 
        ContPair_i%tang_del = zero_r3
    end if
       
    ! relative velocity at contact point 
    call this%clc_ContactVelocity(i,j, Vrij, rveli, rvelj, norm_v, ovrlp , Ri, Rj )  
    
    ! normal and tangential velocities vectors 
    Vij_n = (Vrij .dot. norm_v)*norm_v
    Vij_t = Vrij - Vij_n
    ! normal relative velocity 
    vrn = Vrij .dot. norm_v
    
    ! tangential overlap vector
    ovlp_t = (Vij_t * this%dt) + ContPair_i%tang_del
    ovlp_t = this%rotate_shearPlane (ovlp_t, norm_v) 
   
    ! computing the normal and tangential contact forces 
    fn =  -prop_ij%kn * ovrlp - prop_ij%dmp_n* vrn
    fnij = fn * norm_v
    ftij = ( (-prop_ij%kt) * ovlp_t) - (prop_ij%dmp_t * Vij_t)  
    
    ! Coulomb’s friction law
    ft = norm(ftij)
    ft_fric = prop_ij%fric * abs(fn)
    if( abs(ft) .gt. ft_fric ) then
        
        if( norm(ovlp_t) .gt. 0.0_RK ) then
        
            if( this%CF_Type == CFT_nLin_l ) then
                ! limited overlap
                ftij = ftij*(ft_fric/ft)
                ovlp_t = (-1.0_RK) * (ftij/prop_ij%kt)
            else
                !non-limited overlap 
                ftij = (ftij/ft)*ft_fric
            end if
        
        else
            
            ftij = zero_r3    
            
        end if
        
    end if
    
    ! computing torque acting on spheres it also includes rolling torque 
    call this%clc_Torque( fn, ftij, prop_ij%roll_fric , norm_v , rveli , rvelj, Ri, Rj, Mij, Mji, Mri, Mrj )
    
    ! updating the contact force and torques of particles i and j
    call this%Add_Force( i, j,  fnij+ftij )
    call this%Add_Torque( i, j, Mij+Mri, Mji+Mrj )
    
    ! setting the updated contact info pair into the contact list 
    ContPair_i%tang_del = ovlp_t
    call this%PP_CntctList%setItem(ind , ContPair_i)
    
end subroutine 
    
!**********************************************************************
!   calculating the contact force linear spring-dashpot model and 
! tangential and rolling torques. (Particles-Wall) 
!**********************************************************************
subroutine LSD_ContactForce_PW(this , ind , lnew )
    implicit none
    class(LSD_ContactForce)   this
    integer(IK),intent(in) :: ind
    logical,intent(in)     :: lnew
    
    ! locals 
    integer(IK) i, j
    integer(IK) pt_i, wallID, w_pt
    real(RK)    Ri, Rj, fn, ft , ft_fric
    real(RK)    vrn, ovrlp
    type(real3) rveli, rvelj
    type(real3) norm_v
    type(real3) Vrij, Vij_n, Vij_t, ovlp_t
    type(real3) fnij, ftij, Mij, Mji, Mri, Mrj 
    type(ContactInfo1) ContPair_i
    type(BinaryProperty) Prop_ij
        
    !// body
    ! getting the contacting pairs from contact list 
    ContPair_i = this%PW_CntctList%getItem(ind)
    
    ! i = particle
    ! j = wall
    i = ContPair_i%pari
    wallID = ContPair_i%id_j 
    
    ! particle property type
    pt_i = this%get_pType(i)
    
    if( lnew) then
       ! if this is new contact, then the normal and tangential vectors should be calculated 
        ContPair_i%tang_del = zero_r3
    end if
    
    ! relative velocity at contact point 
    call this%clc_ContactVelocity_PW( i, wallID, Vrij, rveli, norm_v, ovrlp, Ri, w_pt)  
    
    ! normal and tangential velocities vectors 
    Vij_n = (Vrij .dot. norm_v)*norm_v
    Vij_t = Vrij - Vij_n
    vrn = Vrij .dot. norm_v  
     
    ! getting particle-wall property 
    prop_ij = this%Ph_Props%Bnry_getPropPW( pt_i , w_pt )

    
    ! tangential overlap vector
    ovlp_t = (Vij_t * this%dt) + ContPair_i%tang_del
    ovlp_t = this%rotate_shearPlane (ovlp_t, norm_v) 
          
    ! computing the normal and tangential contact forces 
    fn =  -prop_ij%kn * ovrlp - prop_ij%dmp_n* vrn
    fnij = fn * norm_v
    ftij = ( (-prop_ij%kt) * ovlp_t) - (prop_ij%dmp_t * Vij_t)  
    
    ! Coulomb’s friction law
    ft = norm(ftij)
    ft_fric = prop_ij%fric * abs(fn)
    if( abs(ft) .gt. ft_fric ) then
        
        if( norm(ovlp_t) .gt. 0.0_RK ) then
        
            if( this%CF_Type == CFT_nLin_l ) then
                ! limited overlap
                ftij = ftij*(ft_fric/ft)
                ovlp_t = (-1.0_RK) * (ftij/prop_ij%kt)
            else
                !non-limited overlap 
                ftij = (ftij/ft)*ft_fric
            end if
        
        else
            
            ftij = zero_r3    
            
        end if
        
    end if
    
    
    ! computing torque acting on spheres it also includes rolling torque 
    Rj = 100000000000.0_RK
    call this%clc_Torque( fn , ftij, prop_ij%roll_fric, norm_v, rveli, zero_r3, Ri, Rj , Mij, Mji, Mri, Mrj )
    
    ! updating the contact force and torques of particle i   
    call this%Add_Force_PW( i, fnij+ftij )
    call this%Add_Torque_PW( i, (Mij+Mri) )
    
    
    ! setting the updated contact pair into the contact list 
    ContPair_i%tang_del = ovlp_t
    call this%PW_CntctList%setItem(ind , ContPair_i)
    
end subroutine 



subroutine LSD_clc_ContactVelocity(this, i, j , Vrij, wi, wj, norm_v, ovrlp, Ri, Rj )
    implicit none
    class(LSD_ContactForce) this
    integer(IK),intent(in)  :: i, j
    type(real3),intent(out) :: Vrij, wi, wj, norm_v
    real(RK),   intent(out) :: ovrlp
   
    !// locals
    real(RK)    Ri, Rj
    
    type(real3) veli, velj
    type(real4) posi, posj
     ! translational velocity of particles
    veli = this%get_vel(i)
    velj = this%get_vel(j)
    
    ! rotational velocity of particles
    wi = this%get_rvel(i)
    wj = this%get_rvel(j)
                
    ! geting the center position of spherical particles
    posi = this%get_pos(i)
    posj = this%get_pos(j)
    
    ! assuming that particles are spherical 
    ! if this is the case, the bounding box info could be used too 
    Ri = 0.5* posi%w  
    Rj = 0.5* posj%w  
    
    ! computing the normal overlap of two spherical particles
    ! this can be different for spapes rather than spheres 
    ovrlp = posi .ovlp. posj    
    
    ! normal vector
    norm_v = posj .nv. posi  ! posj-posi
        
    !(vi - vj ) + cross( (Ri*wi + Rj*wj), nij) ; 
    Vrij = veli - velj  +  ( (Ri*wi + Rj*wj) .cross. norm_v )
    
end subroutine 


subroutine LSD_clc_ContactVelocity_PW(this, i, wall_id, Vrij, wi, norm_v, ovrlp, Ri, w_pt )
    implicit none
    class(LSD_ContactForce) this
    integer(IK),intent(in)  :: i, wall_id
    type(real3),intent(out) :: Vrij, wi, norm_v
    real(RK),   intent(out) :: ovrlp, Ri 
    integer(IK),intent(out) :: w_pt

    !// locals
    real(RK) dist
    type(real3) veli, velj
    type(real4) posi      
    ! getting particle translational velocity
    veli = this%get_vel(i)
    
    ! getting particle rotational velocity
    wi = this%get_rvel(i)
       
    ! particle center coordinates
    posi = this%get_pos(i)
    
    
    ! getting normal vector, normal distance, velocity at contact point and wall property type
    call this%m_Geometry%normal_dist_vel_type(posi, wall_id, norm_v, dist, velj , w_pt )
    
    ! since the normal vector points from particle i to j we must negate the normal vector
    norm_v = (-1.0_RK) * norm_v
    
    Ri =  0.5* posi%w      
    ovrlp = Ri - dist; 
    
    !(vi - vj ) + cross( (Ri*wi + Rj*wj), nij) ; 
    Vrij = veli - (velj)  +  ( (Ri*wi) .cross. norm_v )
    

end subroutine 


subroutine LSD_clc_Torque( this , fn, ftij,  roll_fric, nij, wi , wj , Ri, Rj, Mij, Mji , Mri, Mrj)
    implicit none
    class(LSD_ContactForce) this
    real(RK),   intent(in)  :: fn, roll_fric
    type(real3),intent(in)  :: ftij, nij ,wi, wj 
    real(RK),   intent(in)  :: Ri, Rj
    type(real3),intent(out) :: Mij, Mji , Mri, Mrj
    
    !// lcoals
    real(RK) Reff, w_hat_mag
    type(real3) M, w_hat
    
    
    ! tangential torque
    M = nij .cross. ftij
    Mij = Ri*M
    Mji = Rj*M
    
    ! rolling torque
    w_hat = wi-wj
    w_hat_mag = norm(w_hat)
    if( w_hat_mag .gt. 0.000001_RK ) then
        w_hat = w_hat/w_hat_mag
    else
        w_hat = zero_r3
    end if
    
    Reff = 1.0_RK /( (1.0_RK/Ri) + (1.0_RK/Rj) )
    
    if( this%CT_Model == CTM_ConstantTorque ) then    
        
        Mri = (-roll_fric*abs(fn)*Reff) * w_hat 
       
       ! removing the normal part 
       ! Mri = Mri - ( (Mri .dot. nij)*nij )
    
    else
        
       Mri = (-roll_fric * abs(fn) * Reff * norm( (Ri*wi+Rj*wj) .cross. nij ) ) * w_hat
       
       ! removing the normal part
       ! Mri = Mri - ( (Mri .dot. nij)*nij )
        
    end if
    
    Mrj = (-1.0_RK)*Mri

end subroutine 


end module
