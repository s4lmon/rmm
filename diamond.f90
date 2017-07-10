program diamond
  use diamondsolve
  implicit none
  real(kind=dp) :: psi,sigma,delta,mu,psi_half,qq,length,source
  integer :: disc,ii
  real(kind=dp),dimension(:),allocatable :: psi_array,x_array,analytical
  !insert variables
  sigma = 1.0_dp
  disc = 10000
  length = 100.0_dp
  source = 10.0_dp
  call dsolve(sigma,disc,length,source,psi_array,x_array)
  allocate(analytical(disc))
  !analytical solution
  do ii=1,disc
    
  end do



end program
