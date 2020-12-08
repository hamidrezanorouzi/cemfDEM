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
!  file name  : g_MakePrtcls.f90
!  module name: g_MakePrtcls    
!                               
!  Purpose:                     
!   A module for defining particles with size distribution and assigning
!    properties to them and positioning them.
!                                            
!                                            
!------------------------------------------------------------------------------

module g_MakePrtcls

    
    use g_Prtcl_PureProperty
    use g_error_handling
    use g_RandomNum
   
    
    implicit none
    
    
    ! size distribution of particles 
type PS_Distribution
        
        integer(IK) numPrtcl ! number of particles 
        integer(IK) numBins  ! number of bins 
        integer(IK) PSD      ! particle size distribution type
        real(RK):: minD = 0.001_RK, maxD = 0.001_RK
        integer(IK),dimension(:),allocatable:: ids          ! id vector of particles 
        type(real4),dimension(:),allocatable:: prtcls       ! dpos vector of particles 
        integer(IK),dimension(:),allocatable:: Bin_idx      ! bin index of particles 
        real(RK),dimension(:),allocatable   :: Bin_meanSize ! mean size of particles in each bin
        
        
    contains
    

    !
    procedure:: MakeDistribution        => PS_MakeDistribution
    
    ! FORTRAN takes care of copy and assignment operator with deep copying mechanism
    ! so, no need to define assignment operator
    
    ! adding two size distributions
    procedure PS_add
    generic   :: operator(+)    => PS_add
    
    procedure  PS_mult
    generic   :: operator(*)    => PS_mult
    
    ! get methods
    procedure:: get_numPrtcl    => PS_get_numPrtcl
    procedure:: get_ids         => PS_get_ids
    procedure:: get_dpos        => PS_get_dpos
    
    ! uniform size distribution
    procedure,private:: make_uniform    => PSD_make_uniform   
    
    ! normal size distribution
    procedure,private:: make_normal     => PSD_make_normal
    
    ! log-normal size distribution 
    procedure,private:: make_lognormal  => PSD_make_lognormal
    
    
    procedure,private:: PS_assign_Bins
    procedure,private:: PS_reallocate
    
    ! final procedure 
    final:: PS_final
    
end type
    
    ! constructor 
    interface PS_Distribution
        procedure PS_Distribution_const
    end interface
    
    ! particles with property
    type, extends(PS_Distribution):: PSD_Property
        
        integer(IK) num_props
	    integer,allocatable:: dummy1
        type(PureProperty),dimension(:),allocatable:: props
        
        integer(IK),dimension(:),allocatable:: p_type   ! property type
        real(RK),dimension(:),allocatable   :: mass     ! mass of particle
        real(RK),dimension(:),allocatable   :: Inertia  ! inertia of particle
        
    contains
    
        procedure:: PSD_Initialize
        
        ! adding two particle size distributions with property
        procedure:: PSD_add
         generic   :: operator(+)    => PSD_add
        
        ! get methods
        procedure:: get_prop        => PSD_get_prop
        procedure:: get_num_ptype   => PSD_get_num_ptype
        procedure:: get_mass        => PSD_get_mass
        procedure:: get_Inertia     => PSD_get_Inertia
        procedure:: get_type        => PSD_get_type
        
        ! memory allocation
        procedure,private:: PSD_reallocate
        
        ! final procedure 
        final:: PSD_final
    end type
    
    ! constructor 
    interface PSD_Property
        procedure PSD_Property_const
    end interface
    
    !// positioning particles  
    type, extends(PSD_Property):: PSDP_Position
        type(real3),dimension(:),allocatable:: vel   ! linear velocity of particles 
    
    contains
        procedure:: PositionOrdered   => PSDP_PositionOrdered
        procedure:: PositionRandom    => PSDP_PositionRandom
        procedure:: get_vel           => PSDP_get_vel
        
        procedure :: PSDP_add
        generic   :: operator(+)    => PSDP_add
        
        procedure,private:: PSDP_reallocate
    end type
    
    !// constructors 
    interface PSDP_Position
        procedure PSDP_Position_const, PSDP_Position_const2
    end interface
    
    contains


!*******************************************************************************    
! constructor for size distribution
!*******************************************************************************
function PS_Distribution_const( numPrtcl, PSD, numBin, minD, maxD ) result( this )
    implicit none
    integer(IK),intent(in):: numPrtcl, PSD, numBin
    real(RK),   intent(in):: minD, maxD
    type(PS_Distribution) this
    
    call this%MakeDistribution( numPrtcl, PSD, numBin, minD, maxD )
    
end function

!******************************************************************************
! adding size distribution (op2) to this size distribution
!******************************************************************************
function PS_add(this , op2 ) result (res)
    implicit none
    class(PS_Distribution),intent(in):: this
    type(PS_Distribution),intent(in):: op2
    type(PS_Distribution) res
    
    !// locals
    integer(IK) i
    
    res%numPrtcl = this%numPrtcl + op2%numPrtcl
    res%numBins  = this%numBins  + op2%numBins
    res%PSD      = IOR(this%PSD , op2%PSD)
    res%minD     = min( this%minD, min(this%maxD, min( op2%minD , op2%maxD) ) )
    res%maxD     = max( this%minD, max(this%maxD, max( op2%minD , op2%maxD) ) )
    
    call res%PS_reallocate()
    
    do i=1,this%numPrtcl
        res%ids(i) = this%ids(i)
        res%prtcls(i) = this%prtcls(i)    
        res%Bin_idx(i) = this%Bin_idx(i)
    end do
    
    
    do i = 1,op2%numPrtcl
        res%ids(this%numPrtcl+i)= op2%ids(i) + this%numPrtcl
        res%prtcls(this%numPrtcl+i)= op2%prtcls(i)
        res%Bin_idx(this%numPrtcl+i)= op2%Bin_idx(i) + this%numBins
    end do
    
    
    res%Bin_meanSize(1:this%numBins) = this%Bin_meanSize
    res%Bin_meanSize(this%numBins+1:)= op2%Bin_meanSize
    
       
end function


!********************************************************************************
!* multiplying the size distribution by an integer number (num)
!********************************************************************************
function PS_mult(this , num ) result (res)
    implicit none
    class(PS_Distribution),intent(in):: this
    integer(IK),intent(in):: num
    type(PS_Distribution) res
    
    !// locals
    integer(IK) i, n
    
    res%numPrtcl = this%numPrtcl * num
    res%numBins  = this%numBins  
    res%PSD      = this%PSD
    res%minD     = this%minD
    res%maxD     = this%minD
    
    call res%PS_reallocate()
    
    do n = 1, num
        do i=1,this%numPrtcl
            res%ids    (i+(n-1)*this%numPrtcl) = i+(n-1)*this%numPrtcl
            res%prtcls (i+(n-1)*this%numPrtcl) = this%prtcls(i)    
            res%Bin_idx(i+(n-1)*this%numPrtcl) = this%Bin_idx(i)
        end do
    end do
    
    res%Bin_meanSize(1:this%numBins) = this%Bin_meanSize
    
end function

!*******************************************************
! returning number of particles 
!*******************************************************
integer(IK) function PS_get_numPrtcl(this)
    implicit none
    class(PS_Distribution) this
    PS_get_numPrtcl = this%numPrtcl
end function

!*******************************************************
! returning id of particles 
!*******************************************************
integer(IK) function PS_get_ids(this , i)
    implicit none
    class(PS_Distribution) this
    integer(IK),intent(in):: i
    
    PS_get_ids = this%ids(i)
    
end function

!*******************************************************
! returning the position of particle 
!*******************************************************
type(real4) function PS_get_dpos(this, i)
    implicit none
    class(PS_Distribution) this
    integer(IK),intent(in):: i
    
    PS_get_dpos = this%prtcls(i)
    
    end function

!*******************************************************
! making a size distribution
!*******************************************************
subroutine PS_MakeDistribution( this, numPrtcl, PSD, numBins, minD, maxD )
    implicit none
    class(PS_Distribution) this
    integer(IK),intent(in):: numPrtcl   ! number of particles 
    integer(IK),intent(in):: PSD        ! type of size distribution
    integer(IK),intent(in):: numBins    ! number of bins 
    real(RK),   intent(in):: minD, maxD ! minimum and maximum size of particles
    
    !// locals
    integer(IK) i
    real(RK) mu, std
    
    !// body    
    this%numPrtcl = numPrtcl
    this%numBins  = numBins
    this%PSD      = PSD
    this%maxD     = maxD
    this%minD     = minD
    
    mu = minD
    std= maxD
    
    !// allcoating memory
    call this%PS_reallocate()
    
    do i= 1,this%numPrtcl
        this%ids(i) = i    
    end do
    
    
    select case ( this%PSD ) 
        
    case(PSD_Uniform)
        call this%make_uniform()
        return
    case(PSD_Normal)
        this%minD = mu-2*std
        this%maxD = mu+2*std
        call this%make_normal(mu, std)
        return
        
    case(PSD_logNormal)
        this%minD = max(0.0_RK, mu-2*std)
        this%maxD = mu+4*std
        call this%make_lognormal(mu,std)
        return
        
    case default 
        call CheckForError( ErrT_Pass, "PS_MakeDistribution", "Wrong input distribution type, set to uniform with one bin")
        this%numBins = 1
        this%PSD = PSD_Uniform
        call this%make_uniform()
        
        return
        
    end select
    
    
    end subroutine 


!*******************************************************
! creating a uniform size distribution 
!*******************************************************
subroutine PSD_make_uniform(this)
    implicit none
    class(PS_Distribution) this
    
    !// locals
    integer(IK) i
    integer(IK) Bin_num
    real(RK) size_dx, r
    
       
    call this%PS_assign_Bins(size_dx)  
    
    
    do i = 1, this%numPrtcl
    
        call RANDOM_NUMBER(r)
        r = this%minD + r*(this%maxD-this%minD)
        Bin_num = int( (r-this%minD)/size_dx ) + 1
        Bin_num = min( this%numBins , max(1_IK, Bin_num) )
        this%prtcls(i)%w = this%Bin_meanSize(Bin_num)
        this%Bin_idx(i) = Bin_num
        
    end do
    
end subroutine

!**********************************************************************
! creating a normal size distribtion
!**********************************************************************
subroutine PSD_make_normal( this , mu, std )
    implicit none
    class(PS_Distribution) this
    real(RK),   intent(in):: mu, std ! mean and standard deviation of distribution
    
    !// locals
    real(RK) size_dx, r
    integer(IK) i, Bin_num
    
    call this%PS_assign_Bins(size_dx) 
    
    
    do i = 1, this%numPrtcl
    
        r = random_normal_dis( mu, std, this%minD, this%maxD )
        Bin_num = int( (r-this%minD)/size_dx ) + 1
        Bin_num = min( this%numBins , max(1_IK, Bin_num) )
        this%prtcls(i)%w = this%Bin_meanSize(Bin_num)
        this%Bin_idx(i) = Bin_num
        
    end do
    
    
        
    end subroutine

!*******************************************************
! creating a log-normal size distribution
!*******************************************************
subroutine PSD_make_lognormal( this , mu, std )
    implicit none
    class(PS_Distribution) this
    real(RK), intent(in)::mu, std  ! mean and standard deviation 
    
    !// locals
    real(RK) size_dx, r
    integer(IK) i, Bin_num
    
    call this%PS_assign_Bins(size_dx) 
    
    
    do i = 1, this%numPrtcl
    
        r = random_logNormal( mu, std, this%minD, this%maxD )
        
        Bin_num = int( (r-this%minD)/size_dx ) + 1
        Bin_num = min( this%numBins , max(1_IK, Bin_num) )
        this%prtcls(i)%w = this%Bin_meanSize(Bin_num)
        this%Bin_idx(i) = Bin_num
        
    end do
    
        
    end subroutine 

!*******************************************************
! assigning bins
!*******************************************************
subroutine PS_assign_Bins(this, size_dx )
    implicit none
    class(PS_Distribution)     this
    real(RK),intent(inout)  :: size_dx
    
    !// locals
    integer(IK) i, nBins
    
    !// body
    
    nBins = this%numBins
    
    size_dx = (this%maxD-this%minD)/nBins
    if( nBins == 1 .and. size_dx .lt.  (0.001_RK*this%maxD) ) size_dx = 10*this%maxD
    
    if(nBins == 1 ) then
        
        this%Bin_meanSize(1) = 0.5_RK * (this%maxD+this%minD) 
        
    else
        
        do i=1,nBins
            this%Bin_meanSize(i) = this%minD + (i-1)*size_dx + 0.5_RK*size_dx    
        end do
        
    end if
    

    end subroutine

!*******************************************************
! memory allocation 
!*******************************************************
subroutine PS_reallocate(this)
    implicit none
    class(PS_Distribution) this
    
    if(allocated( this%ids) ) deallocate(this%ids)
    allocate( this%ids( this%numPrtcl) )
    
    if(allocated( this%prtcls ) ) deallocate(this%prtcls)
    allocate( this%prtcls( this%numPrtcl) )
    
    if(allocated( this%Bin_idx ) ) deallocate(this%Bin_idx)
    allocate( this%Bin_idx( this%numPrtcl) )
    
    if(allocated( this%Bin_meanSize) ) deallocate( this%Bin_meanSize) 
    allocate( this%Bin_meanSize( this%numBins) )
    
    
end subroutine 

!*******************************************************
! final procedure 
!*******************************************************
subroutine PS_final(this)
    implicit none
    type(PS_Distribution) this
    
    if(allocated( this%ids) ) deallocate(this%ids)
    
    if(allocated( this%prtcls ) ) deallocate(this%prtcls)
    
    if(allocated( this%Bin_idx ) ) deallocate(this%Bin_idx)
        
    if(allocated( this%Bin_meanSize) ) deallocate( this%Bin_meanSize) 
    
    end subroutine
    
    
!///////////////////////////////////////////////////////////////////////////////////
!///////////////////////////////////////////////////////////////////////////////////
    
!*******************************************************
! constructor for size distribution with property 
!*******************************************************
function PSD_Property_const( PS_dist , prop ) result(res)
    implicit none
    type(PS_Distribution), intent(in):: PS_dist
    type(PureProperty),    intent(in):: prop
    type(PSD_Property) res
    
    !print*, "PSD_Property , prop" , prop
    call res%PSD_Initialize( PS_dist, prop )
    
    
end function

!*******************************************************
! initializing the size distribution with property 
!*******************************************************   
subroutine PSD_Initialize( this, PS_dist , prop ) 
    implicit none
    class(PSD_Property)                  this
    type(PS_Distribution), intent(in):: PS_dist
    type(PureProperty),    intent(in):: prop
    
    !// locals
    integer(IK) i
    
    this%PS_Distribution = PS_dist
    this%num_props = this%numBins
    call this%PSD_reallocate()
    
    do i = 1, this%num_props
        
        this%props(i) = prop
	    this%props(i)%rad = 0.5_RK * this%Bin_meanSize(i)
        call this%props(i)%clc_prop()
         
    end do
    
    do i=1, this%numPrtcl
        this%p_type(i) = this%Bin_idx(i)
        this%mass(i)   = this%props( this%Bin_idx(i) )%mass
        this%Inertia(i)= this%props( this%Bin_idx(i) )%inertia
    end do
    
end subroutine

!****************************************************************
! adding the input size distribution with property (op2) to this 
!****************************************************************
function PSD_add(this, op2 ) result(res)
    implicit none
    class(PSD_Property),intent(in):: this
    type(PSD_Property),intent(in):: op2
    type(PSD_Property) res
    
    !// locals
    integer(IK) i
    
    
    res%PS_Distribution = this%PS_Distribution%PS_add(op2%PS_Distribution)
    
    res%num_props = this%num_props + op2%num_props
    
    call res%PSD_reallocate()
    
    
    res%props(1:this%num_props) = this%props
    res%props(this%num_props+1:)= op2%props
    
    do i =1, this%numPrtcl
        res%p_type(i)   = this%p_type(i)
        res%mass(i)     = this%mass(i)
        res%Inertia(i)  = this%Inertia(i)
    end do
    
    do i=1,op2%numPrtcl
        res%p_type(this%numPrtcl+i) = op2%p_type(i)+this%num_props
        res%mass(this%numPrtcl+i)   = op2%mass(i)
        res%Inertia(this%numPrtcl+i)= op2%Inertia(i)
    end do
    
    
end function

!*******************************************************
! returning property
!*******************************************************
type(PureProperty) function PSD_get_prop(this , i)
    implicit none
    class(PSD_Property) this
    integer(IK),intent(in):: i
    
    PSD_get_prop = this%props(i)
    
end function

!*******************************************************
! returning number of property types 
!*******************************************************
integer(IK) function PSD_get_num_ptype(this)
    implicit none
    class(PSD_Property) this
    
    PSD_get_num_ptype = this%num_props
    
end function

!*******************************************************
! returning mass of particle
!*******************************************************
real(RK) function PSD_get_mass(this, i)
    implicit none
    class(PSD_Property) this
    integer(IK),intent(in):: i
    
    PSD_get_mass = this%mass(i)
    
end function

!*******************************************************
! returning inertia of particle 
!*******************************************************
real(RK) function PSD_get_Inertia(this, i)
    implicit none
    class(PSD_Property) this
    integer(IK),intent(in):: i
    
    PSD_get_Inertia = this%Inertia(i)
    
end function

!*******************************************************
! returning property type of particle 
!*******************************************************
Integer(IK) function PSD_get_type(this, i)
    implicit none
    class(PSD_Property) this
    integer(IK),intent(in):: i
    
    PSD_get_type = this%p_type(i)
    
end function

subroutine PSD_reallocate(this)
    implicit none
    class(PSD_Property) this
    
    if( allocated(this%props) ) deallocate( this%props )
    allocate( this%props( this%num_props ) )
    
    if( allocated(this%p_type) ) deallocate( this%p_type )
    allocate( this%p_type( this%get_numPrtcl() ) )
    
    if( allocated(this%mass) ) deallocate( this%mass )
    allocate( this%mass( this%get_numPrtcl() ) )
    
    if( allocated(this%Inertia) ) deallocate( this%Inertia )
    allocate( this%Inertia( this%get_numPrtcl() ) )
    
end subroutine

subroutine PSD_final(this)
    implicit none
    type(PSD_Property) this
    
    if( allocated(this%props) ) deallocate( this%props )
    
    if( allocated(this%p_type) ) deallocate( this%p_type )
        
    if( allocated(this%mass) ) deallocate( this%mass )
        
    if( allocated(this%Inertia) ) deallocate( this%Inertia )
    
    end subroutine


!////////////////////////////////////////////////////////////////////////////////////////
!////////////////////////////////////////////////////////////////////////////////////////

!*******************************************************
! constructor for PSDP_Position based on PSDP as input
!*******************************************************
function PSDP_Position_const(  PSDP ) result (res)
    implicit none
    type(PSD_Property),intent(in):: PSDP 
    type( PSDP_Position) res
    res%PSD_Property = PSDP
   
    call res%PSDP_reallocate()

end function

!**************************************************************************
! constructor for PSDP_Position based on the input file (chFile)
!**************************************************************************
function PSDP_Position_const2(chFile) result(res)
    implicit none
    character(*),intent(in):: chFile
    type(PSDP_Position) res
    
    !// locals 
    integer(IK) iStat
    integer(IK) nUnit, numPrtcl, numProps, i, id
    
    !// body
    nUnit = d_wrt_prtc_unit
    open( unit = nUnit , file = chFile ) 
    
    !// first reading number of particles 
    read( nUnit, *, IOSTAT = iStat ) numPrtcl
    if( iStat .ne. 0 ) then
        call CheckForError(ErrT_Abort, "PSDP_Position_const2" , "Error in reading file :"//trim(chfile) )
        return
    end if
    
    !// reading number of properties
    read( nUnit, *, IOSTAT = iStat) numProps
    if( iStat .ne. 0 ) then
        call CheckForError(ErrT_Abort, "PSDP_Position_const2" , "Error in reading file :"//trim(chfile) )
        return
    end if
    
    res%numPrtcl = numPrtcl
    res%numBins  = numProps
    res%num_props= numProps
    call res%PSD_Property%PS_Distribution%PS_reallocate()
    call res%PSD_Property%PSD_reallocate()
    call res%PSDP_reallocate()
    
    !// reading properties
    do i= 1, numProps
        
        read( nUnit , *, IOSTAT = iStat ) res%props(i)
        if( iStat .ne. 0 ) then
            call CheckForError(ErrT_Abort, "PSDP_Position_const2" , "Error in reading property data from file :"//trim(chfile) )
            return
        end if
        
    end do
    
    !// reading particle data
    do i=1, numPrtcl
        
        read(nUnit , *, IOSTAT = iStat ) id , res%p_type(i), res%prtcls(i), res%vel(i)
        if( iStat .ne. 0 ) then
            call CheckForError(ErrT_Abort, "PSDP_Position_const2" , "Error in reading particle data from file :"//trim(chfile) )
            return
        end if
        res%ids(i) = i
        res%Bin_idx(i) = res%p_type(i)
        res%mass(i)    = res%props( res%p_type(i) )%mass
        res%Inertia(i) = res%props( res%p_type(i) )%inertia
        
    end do
    
    do i = 1, numProps
        res%Bin_meanSize(i) = 2.0_RK * res%props(i)%Rad
    end do
    
    
    close(nUnit)
    end function

!***********************************************************************************
! randomly positioning particles in a box with corner points min_p and max_p and 
!    the initial velocity of vel 
!***********************************************************************************
subroutine PSDP_PositionRandom( this , min_p, max_p , vel )
    implicit none
    class(PSDP_Position) this
    type(real3),    intent(in):: min_p, max_p
    type(real3),    intent(in):: vel
    optional:: vel
    
    integer(IK) n, numPrtcl, iter
    real(RK)    r, x, y, z, maxD
    type(real3) l_vel, lmin_p, ddx
    
    
    if( present(vel) )then
        l_vel = vel
    else
        l_vel = real3( 0.0, 0.0, 0.0)
    end if
    
    numPrtcl = this%get_numPrtcl()        
    maxD     = this%maxD
    
    lmin_p = min_p + (0.5_RK* real3(maxD,maxD,maxD) )
    ddx = (max_p - min_p) - real3(maxD,maxD,maxD)
    
    iter = 0      
    n = 0
    do while( n < this%numPrtcl )
        
        iter = iter+1
        call random_number(r)
        x = lmin_p%x + r*ddx%x 

        call random_number(r)
        y = lmin_p%y + r*ddx%y

        call random_number(r)
        z = lmin_p%z + r*ddx%z

        if( .not. inCont(n , this%prtcls , real4(x,y,z, this%prtcls(n+1)%w ) ) ) then
            
            n = n + 1
        
            this%Prtcls(n)%x = x
            this%Prtcls(n)%y = y
            this%Prtcls(n)%z = z
            this%vel(n)      = l_vel
        end if
        
        if( iter > 10*numPrtcl ) then
            call CheckForError( ErrT_Abort, "PSDP_PositionRandom" ,&
                                "Cannot positoin particles randomly in the specified box" )
        end if
        
    end do
    
contains
    
    logical function inCont( n , prtcls, newp )
        integer(IK) n
        type(real4),dimension(:):: prtcls
        type(real4) newp
        
        integer(IK) i
        
        inCont = .true.
        do i = 1, n
        
            if( ( prtcls(i) .ovlp. newp ) > 0.0000001_RK ) return
                
        end do
        
        inCont = .false.
        
    end function

    end subroutine

!*******************************************************************************************
! positioning particles in a box with corner points min_p and max_p and the initial 
!     velocity of vel. Positioning is done along main axes in order that is defined with
!     fill_order 
!******************************************************************************************
subroutine PSDP_PositionOrdered( this,  min_p, max_p , fill_order , vel )
    implicit none
    class(PSDP_Position) this
    type(real3),    intent(in):: min_p, max_p
    type(integer3), intent(in):: fill_order
    type(real3),    intent(in):: vel
    optional:: fill_order, vel
    
    integer(IK) n, numPrtcl
    real(RK) maxD, k
    type(integer3) l_fill_order
    type(real3) l_vel, cntr
    type(real3) l1_vec, l2_vec, l3_vec
    type(real3) lmin_p, lmax_p, rad_dx , dx, eps
    
    ! default or user defined velocity
    if( present(vel) )then
        l_vel = vel
    else
        l_vel = real3( 0.0, 0.0, 0.0)
    end if
    
    ! default or user defined fill order 
    if( present(fill_order) ) then
        l_fill_order = fill_order     
    else
        l_fill_order = integer3(x_axis,y_axis,z_axis)
    endif  
    
    
    
    maxD = maxval(this%Bin_meanSize)
    eps = 0.001_RK * real3(maxD, maxD, maxD)
    k = 1.0_RK
    
    call Vectors(l_fill_order, l1_vec, l2_vec, l3_vec)
    dx = real3(maxD, maxD, maxD)
    rad_dx = 0.5_RK * real3(maxD, maxD, maxD)
      
    ! 
    lmin_p = min_p + rad_dx 
    lmax_p = max_p - rad_dx
    
     ! start point
    cntr = lmin_p 
    
    numPrtcl = this%get_numPrtcl()
    n = 0
    
    do while( n < numPrtcl )
        
        n = n + 1
        this%prtcls(n) = real4( cntr%x, cntr%y, cntr%z, this%prtcls(n)%w )
        this%vel(n)    = l_vel
                      ! adds a dx to l1
        cntr = cntr + (l1_vec * dx) + (k*eps)
        k = -k
        ! comparing l1 values
        if( ( l1_vec .dot. cntr) .ge. (l1_vec .dot. lmax_p ) )then
            !       rewinding l1       adding a dx to l2  keeping previous value of l3    
            cntr = (lmin_p*l1_vec) + ( (cntr+dx) * l2_vec) + (cntr*l3_vec)
            
            ! compares l2 values
            if( (l2_vec .dot. cntr) .ge. (l2_vec .dot. lmax_p) ) then
                !   keeping value of l1   rewinding l2    adding a dx to l3
                cntr = (cntr*l1_vec) + (lmin_p*l2_vec) + ((cntr+dx)*l3_vec)
                
                ! comparing l3 values
                if( (l3_vec .dot. cntr) .ge. (l3_vec .dot. lmax_p) ) then
                   call CheckForError( ErrT_Abort, "PSDP_PositionOrdered" ,"Not enough space for positioning particles" ) 
                end if
                
            end if
            
        end if
        
    end do

contains

subroutine Vectors( fill_order , l1 , l2, l3 )
    implicit none
    type(integer3), intent(in) :: fill_order
    type(real3),    intent(out):: l1, l2, l3
    
    if(fill_order%x == z_axis ) then
        l1 = real3(0.0,0.0,1.0)
    elseif(fill_order%x == y_axis) then
        l1 = real3(0.0,1.0,0.0)
    else
        l1 = real3(1.0,0.0,0.0)
    endif
    
    if(fill_order%y == x_axis ) then
        l2 = real3(1.0,0.0,0.0)
    elseif(fill_order%y==z_axis) then
        l2 = real3(0.0,0.0,1.0)
    else
        l2 = real3(0.0,1.0,0.0)
    endif
    
    if(fill_order%z == x_axis) then
        l3 = real3(1.0,0.0,0.0)
    elseif(fill_order%z == y_axis) then
        l3 = real3(0.0,1.0,0.0)
    else
        l3 = real3(0.0,0.0,1.0)
    endif
    
end subroutine
    
    end subroutine

!*******************************************************
! returning velocity of particle i
!*******************************************************
type(real3) function PSDP_get_vel(this, i)
    implicit none
    class(PSDP_Position) this
    integer(IK),intent(in):: i
    
    PSDP_get_vel = this%vel(i)
    
end function
    
!*******************************************************
! adding the input PSDP_Position (op2) to this
!*******************************************************
function PSDP_add(this, op2 ) result(res)
    implicit none
    class(PSDP_Position),intent(in):: this
    type(PSDP_Position),intent(in):: op2
    type(PSDP_Position) res
    
    !// locals
    integer(IK) i
    
    
    res%PSD_Property = this%PSD_Property%PSD_add(op2%PSD_Property)
    

    res%num_props = this%num_props + op2%num_props
    
    call res%PSDP_reallocate() 
    
    do i =1, this%numPrtcl
        res%vel(i)   = this%vel(i)
    end do
    
    do i=1,op2%numPrtcl
        res%vel(this%numPrtcl+i) = op2%vel(i)
    end do
    
    
end function

    
    
subroutine PSDP_reallocate(this)
    implicit none
    class(PSDP_Position) this
    
    if( allocated(this%vel) ) deallocate( this%vel) 
    allocate( this%vel( this%get_numPrtcl() ) )
    
end subroutine 


end module
