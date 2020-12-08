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
!  file name  : g_Prtcl_NonLin_Model.f90 
!  module name: g_Prtcl_NonLin_Model 
!        
!  Purpose: 
!    1) Extending LSD_ContactForce to calculate particle-particle and particle-
!       wall contact forces based on non-linear model
! 
!------------------------------------------------------------------------------
    
module g_Prtcl_NonLin_Model

    use g_Prtcl_LSD_Model
    implicit none
    
    
    type, extends (LSD_ContactForce):: NonLin_ContactForce
    
    contains
    
        procedure:: ContactForce_PP     => NLin_ContactForce_PP
        procedure:: ContactForce_PW     => NLin_ContactForce_PW
    
    end type
    
    
contains

! particle-particle contact force 
subroutine NLin_ContactForce_PP(this , ind , lnew )
    implicit none
    class(NonLin_ContactForce) this
    integer(IK),intent(in):: ind
    logical,intent(in)    :: lnew
    
    !// locals
    integer(IK) i, j
    integer(IK) pt_i , pt_j
    real(RK)    Ri, Rj, fn, ft , ft_fric, kt
    real(RK)    vrn , ovrlp
    type(real3) rveli, rvelj, dtij, dftij, ovlp_t
    type(real3) norm_v 
    type(real3) Vrij, Vij_n, Vij_t
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
    prop_ij = this%Ph_Props%Bnry_getProp( pt_i , pt_j )
    
    if( lnew) then
        ! if this is new contact, the tangential overlap should be zero 
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
    ovlp_t = (Vij_t * this%dt) + ContPair_i%tang_del;
    ovlp_t = this%rotate_shearPlane (ovlp_t, norm_v) 
        
    ! computing the normal and tangential contact forces 
    fn = - (4.0_RK/3.0_RK   * prop_ij%Yeff * sqrt(prop_ij%Reff) * ovrlp**1.5_RK) - (prop_ij%dmp_n*ovrlp**0.25_RK * vrn)
    fnij = fn * norm_v
    ftij = (- 16.0_RK/3.0_RK * prop_ij%Geff * sqrt(prop_ij%Reff*ovrlp) ) * ovlp_t
    
    ! Coulomb's friction law
    ft = norm(ftij)
    ft_fric = prop_ij%fric * abs(fn)
    if( abs(ft) .gt. ft_fric ) then
        
        if( norm(ovlp_t) .gt. 0.0_RK ) then
        
            if( this%CF_Type == CFT_nLin_l ) then
                ! limited overlap
                ftij = ftij*(ft_fric/ft)
                kt = 16.0_RK/3.0_RK * prop_ij%Geff * sqrt(prop_ij%Reff*ovrlp)
                ovlp_t = (-1.0_RK) * (ftij/kt)
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


! particle-wall contact force 
subroutine NLin_ContactForce_PW(this , ind , lnew )
    implicit none
    class(NonLin_ContactForce)   this
    integer(IK),intent(in) :: ind
    logical,intent(in)     :: lnew
    
    !//locals 
    integer(IK) i, j
    integer(IK) pt_i , pt_j
    integer(IK) wallID, w_pt
    real(RK) Ri, Rj, fn, ft , ft_fric, kt
    real(RK) vrn, dist, ovrlp
    real(RK) mi
    type(real3) g, gg   
    type(real3) rveli, rvelj, dtij, dftij, ovlp_t
    type(real3) norm_v
    type(real3) Vrij, Vij_n, Vij_t
    type(real3) fnij, ftij, Mij, Mji, Mri, Mrj , wi
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
   ovlp_t = (Vij_t * this%dt) + ContPair_i%tang_del;
   ovlp_t = this%rotate_shearPlane (ovlp_t, norm_v) 
   
    
    ! computing the normal and tangential contact forces 
    fn = - (4.0_RK/3.0_RK   * prop_ij%Yeff * sqrt(prop_ij%Reff) * ovrlp**1.5_RK) - (prop_ij%dmp_n*ovrlp**0.25_RK * vrn)
    fnij = fn * norm_v
    ftij = (- 16.0_RK/3.0_RK * prop_ij%Geff * sqrt(prop_ij%Reff*ovrlp) ) * ovlp_t
    
    
    ! Coulomb's friction law
    ft = norm(ftij)
    ft_fric = prop_ij%fric * abs(fn)
    if( abs(ft) .gt. ft_fric ) then
        
        if( norm(ovlp_t) .gt. 0.0_RK ) then
        
            if( this%CF_Type == CFT_nLin_l ) then
                ! limited overlap
                ftij = ftij*(ft_fric/ft)
                kt = 16.0_RK/3.0_RK * prop_ij%Geff * sqrt(prop_ij%Reff*ovrlp)
                ovlp_t = (-1.0_RK) * (ftij/kt)
            else
                !non-limited overlap 
                ftij = (ftij/ft)*ft_fric
            end if
        
        else
            
            ftij = zero_r3    
            
        end if
        
    end if
    
    

    ! computing torque acting on spheres it also includes rolling torque 
    Rj = 1000000000.0_RK
    call this%clc_Torque( fn , ftij, prop_ij%roll_fric, norm_v, rveli, zero_r3, Ri, Rj , Mij, Mji, Mri, Mrj )
    
    ! updating the contact force and torques of particle i   
    call this%Add_Force_PW( i, fnij+ftij )
    call this%Add_Torque_PW( i, (Mij+Mri) )
    
    
    ! setting the updated contact pair into the contact list 
    ContPair_i%tang_del = ovlp_t
    call this%PW_CntctList%setItem(ind , ContPair_i)
    
end subroutine 

end module
