module diamondsolve
  use constants
  implicit none
real(kind=dp),parameter :: source = 1.0_dp
contains

  subroutine dsolve(lhs,sigma,disc,length,q_array,psi_array,x_array,rhs)
    real(kind=dp) :: psi,sigma,delta,mu,psi_half,length,rhs,excess,lhs


    integer :: disc,ii,temp
    real(kind=dp),dimension(:) :: psi_array,x_array,q_array

    delta=length/dble(disc)
    mu=1


    psi_half=lhs




    x_array=delta/2.0_dp
      do ii = 2,disc
        x_array(ii)=x_array(ii-1)+delta

      end do

      do  ii=1,disc


        psi=((1.0_dp+(sigma*delta)/(2.0_dp*abs(mu)))**(-1))*(psi_half+((delta*q_array(ii))/(2.0_dp*mu)))
        psi_array(ii)=psi


        psi_half=2*psi-psi_half

      end do


  rhs=psi_half



  end subroutine


end module
