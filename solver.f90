program solver
  use rmm
  implicit none
    real(kind=dp) :: domainLength,subLength
    integer :: numDisc, numSub, ii
      real(kind=dp),dimension(4) :: sources,outputs
      real(kind=dp),dimension(4) :: captureCS
      real(kind=dp),dimension(4,4) :: particleCS,factors
      real(kind=dp) :: sourceFN,sourceTN, sourcePG, sourceSG
      real(kind=dp),dimension(:),allocatable :: psiArray,xArray
      real(kind=dp) :: fn_tn,fn_sg,tn_sg,pg_sg
      real(kind=dp),dimension(4,4) :: decays


!INPUT VALUES
    domainLength=100.0_dp
    numDisc=1000
    numSub=100
    !call discretiser(domainLength,numSub,numDisc,xArray)
    !print*,xArray
    fn_tn=0.5_dp
    fn_sg=0.5_dp
    tn_sg=0.5_dp
    pg_sg=0.5_dp



!===building for computation
    factors=0.0_dp
     factors(2,1)=fn_tn
     factors(4,1)=fn_sg
     factors(4,2)=tn_sg
     factors(4,3)=pg_sg
    sourceFN=0.0_dp
    sourceTN=0.0_dp
    sourcePG=0.0_dp
    sourceSG=1.0_dp
      sources(1)=sourceFN
      sources(2)=sourceTN
      sources(3)=sourcePG
      sources(4)=sourceSG

    subLength=domainLength/numSub
    particleCS=0.0_dp
    particleCS(1,2)=1.0_dp
    particleCS(1,4)=1.0_dp
    particleCS(2,4)=1.0_dp
    particleCS(3,4)=1.0_dp
    captureCS=(/1.0_dp,1.0_dp,1.0_dp,1.0_dp/)
    call particleSolve(sources,numDisc,captureCS,particleCS,factors,subLength,outputs)
    !call plot()
    !building alpha and beta matrix
    sources=0.0_dp
    sources(1)=1.0_dp
!=====================================================================================
!building decay constants (alpha,beta)
!switching on and off different particles and seeing what scattering they cause
    decays=0.0_dp
    call particleSolve(sources,numDisc,captureCS,particleCS,factors,subLength,outputs)
    decays(1,1)=outputs(1)
    decays(1,2)=outputs(2)
    decays(1,4)=outputs(4)
    sources=0.0_dp
    sources(2)=1.0_dp
    call particleSolve(sources,numDisc,captureCS,particleCS,factors,subLength,outputs)
    decays(2,2)=outputs(2)
    decays(2,4)=outputs(4)
    sources=0.0_dp
    sources(3)=1.0_dp
    call particleSolve(sources,numDisc,captureCS,particleCS,factors,subLength,outputs)
    decays(3,3)=outputs(3)
    decays(3,4)=outputs(4)
    sources=0.0_dp
    sources(1)=1.0_dp
    call particleSolve(sources,numDisc,captureCS,particleCS,factors,subLength,outputs)
    decays(4,4)=outputs(4)
    ! ! do ii=1,4
 !   print*, decays(ii,:)
 ! end do
end program
