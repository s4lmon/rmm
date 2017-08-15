module backdiamond
  use constants
  implicit none
real(kind=dp),parameter :: source = 1.0_dp
contains

  subroutine backsolve(rhs,sigma,disc,length,q_array,psi_array,x_array,lhs)
    real(kind=dp) :: psi,sigma,delta,mu,psi_half,length,rhs,excess,lhs


    integer :: disc,ii,temp
    real(kind=dp),dimension(:) :: psi_array,x_array,q_array

    delta=length/dble(disc)
    mu=1


    psi_half=rhs




    x_array=delta/2.0_dp
      do ii = 2,disc
        x_array(ii)=x_array(ii-1)+delta

      end do

      do  ii=disc,1,-1


        psi_array(ii)=((1.0_dp+(sigma*delta)/(2.0_dp*abs(mu)))**(-1))*(psi_half+((delta*q_array(ii))/(2.0_dp*mu)))
!       unneeded variable
        !psi_array(ii)=psi

!print*,psi_half
        psi_half=2*psi_array(ii)-psi_half

      end do


  lhs=psi_half



  end subroutine


end module
