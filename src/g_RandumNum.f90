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
!  file name  : g_RandomNum.f90
!  module name: g_RandomNum
!                          
!  Purpose:
!    generating the random number with normal and log-normal distributions
! 
!  Adapted from the following Fortran 77 code ALGORITHM 712, COLLECTED 
!    ALGORITHMS FROM ACM.THIS WORK PUBLISHED IN TRANSACTIONS ON MATHEMATICAL
!    SOFTWARE, VOL. 18, NO. 4, DECEMBER, 1992, PP. 434-435.
!
!  The function random_normal() returns a normally distributed pseudo-random
!    number with zero mean and unit variance.
! 
!  The algorithm uses the ratio of uniforms method of A.J. Kinderman and J.F.
!    Monahan augmented with quadratic bounding curves.
! 
!------------------------------------------------------------------------------
    
module g_RandomNum
    
    use g_TypeDef
    implicit none
        
contains

   
recursive function random_normal_dis( mu , std , lowcut, highcut ) result(resf)

implicit none

real(RK) mu, std 
real(RK) lowcut,highcut
real(RK) resf

!     Local variables
real(RK):: s = 0.449871_RK , t = -0.386595_RK, a = 0.19600_RK, b = 0.25472_RK,    &
            r1 = 0.27597_RK, r2 = 0.27846_RK
real(RK)   u, v, x, y, q, res


! Generate P = (u,v) uniform in rectangle enclosing acceptance region

do
  call RANDOM_NUMBER(u)
  call RANDOM_NUMBER(v)
  v = 1.7156_RK  * (v - 0.5_RK)

!     Evaluating the quadratic form
  x = u - s
  y = ABS(v) - t
  q = x**2 + y*(a*y - b*x)

!     Accept P if inside inner ellipse
  if (q < r1) exit
!     Reject P if outside outer ellipse
  if (q > r2) cycle
!     Reject P if outside acceptance region
  if (v**2 < -4.0*log(u)*u**2) exit
end do



!     Returning ratio of P's coordinates as the normal deviate
res = (v/u)*std + mu

if(res.lt. lowcut .or. res .gt. highcut ) then
    
    resf = random_normal_dis( mu , std , lowcut, highcut )
else
    resf = res

end if


    
end function


recursive function random_logNormal( mu, std, lowcut, highcut ) result (resf)

    real(RK) mu, std
    real(RK) lowcut, highcut
    real(RK) resf
    
    real(RK) variance, mean, res
    
    
    variance = log(1+std**2/mu**2)
    mean = log( mu - 0.5*variance)
    
    
        
    res = exp( random_normal_dis(mean, sqrt(variance) , -1000000000.0_RK , 1000000000.0_RK ) )
    
    if(res.lt. lowcut .or. res .gt. highcut ) then
        
        resf = random_logNormal( mu, std, lowcut, highcut )
    else
        resf = res
    end if
    
end function 
    
    
end module
