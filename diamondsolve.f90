!x_array firxed
!todo psi correction.

module diamondsolve
  use constants

  implicit none

contains
  subroutine dsolve(sigma,disc,length,source,psi_array,x_array)
    real(kind=dp) :: psi,sigma,delta,mu,psi_half,length,source
    integer :: disc,ii,temp
    real(kind=dp),dimension(:),allocatable :: psi_array,x_array,q_array

    delta=length/(disc)
    mu=1

    allocate(psi_array(disc))
    allocate(x_array(disc))
    allocate(q_array(disc))

    q_array=0

    psi_half=source


    x_array=delta/2.0_dp
      do ii = 2,disc
        x_array(ii)=x_array(ii-1)+delta
      end do

      do  ii=1,disc
        !psi_half is half step BEFORE current evaluation of psi

        psi=((1.0_dp+(sigma*delta)/(2.0_dp*abs(mu)))**(-1))*(psi_half+((delta*q_array(ii))/(2.0_dp*mu)))
        psi_array(ii)=psi

        print*,'psi::x',psi,x_array(ii)
        psi_half=2*psi-psi_half
        !print*,'psi_half',psi_half

      end do

      ! output data to a file
      open(100, file="data1.dat", status="unknown")
        do ii=1,disc
          write(100,*) x_array(ii), psi_array(ii)
        end do
      close(100)
  end subroutine


end module
