  !todo fix q array

program rmm
  use diamondsolve


  implicit none

  real(kind=dp) :: domainLength,subLength,delta
  real(kind=dp) :: sMultiply !for scaling source parameter (fixed at 1)
  real(kind=dp) :: sourceFN,sourcePG,sourceTN
  real(kind=dp) :: sourceSG = 0.0_dp!unchanging, zero at start

!macroscopic cross section for each kind of particle
  real(kind=dp) :: sigma,sigmaFN,sigmaPG,sigmaTN,sigmaSG !currently sigma_T
!decay constants from source for each kind of particle
  real(kind=dp) :: alpha,alphaFN,alphaPG,alphaTN,alphaSG

!decay constants from other particles for each kind of particle
  real(kind=dp) :: betaFT,betaFS,betaTS,betaPS

!factors for q arrays
  real(kind=dp) :: factorFN,factorPG,factorTN,factorSG
!timing
  real(kind=dp) :: start,finish

!temporary arrays for calculation of discretisations
  real(kind=dp),dimension(:),allocatable :: psiArray,xArray,subPsi,subX,subQ, subQ_FN,subQ_PG,subQ_TN,subQ_SG

!data arrays for each particle (OUTPUT)
  real(kind=dp),dimension(:),allocatable :: psiFN,psiPG,psiTN,psiSG

  integer :: numDisc,numSub,total,ii !how many discretisations necessary
  call cpu_time(start)
  !==============================================================================

  sigma = 1.0_dp
  sigmaFN = 1.0_dp
  sigmaTN = 1.0_dp
  sigmaPG = 1.0_dp
  sigmaSG = 1.0_dp
  sourceFN=10.0_dp
  sourcePG=10.0_dp
  sourceTN=5.0_dp

  numDisc = 10
  numSub=10000

  domainLength = 10.0_dp
  sMultiply=10.0_dp
  total=numDisc*numSub

  !subQ=0.0_dp !!!!!!!!not working???

  !==============================================================================
  subLength=domainLength/numSub
    delta=subLength/numDisc
  allocate(psiArray(total))
  allocate(xArray(total))

  allocate(psiFN(total))
  allocate(psiPG(total))

  allocate(subPsi(numDisc))
  allocate(subX(numDisc))
  allocate(subQ(numDisc))
    allocate(subQ_FN(numDisc))
    allocate(subQ_PG(numDisc))
    allocate(subQ_TN(numDisc))
    allocate(subQ_SG(numDisc))
  subQ=0.0_dp

!building q matrices
  subQ_FN=factorFN*subQ
  subQ_PG=factorPG*subQ
  !solve that subdomain
  !GENERAL NEUTRONS
  !========================================================================
  call dsolve(sigma,numDisc,subLength,smultiply,subQ,subPsi,subX,alpha)
  !print*,subX


    do ii=1,numSub
      psiArray((ii-1)*numDisc+1:ii*numDisc)=subPsi*sMultiply

      xArray((ii-1)*numDisc+1:ii*numDisc)=subX+subLength*(ii-1)
      sMultiply=sMultiply*alpha

    end do

  !FAST NEUTRONS
  !=======================================================
  call dsolve(sigmaFN,numDisc,subLength,sourceFN,subQ,subPsi,subX,alphaFN)
    do ii=1,numSub
      psiFN((ii-1)*numDisc+1:ii*numDisc)=subPsi*sMultiply
      sourceFN=sourceFN*alphaFN

    end do
  !PRIMARY GAMMA
!=========================================================
  call dsolve(sigmaPG,numDisc,subLength,sourcePG,subQ_PG,subPsi,subX,alphaPG)
    do ii=1,numSub
      psiPG((ii-1)*numDisc+1:ii*numDisc)=subPsi*sMultiply
      sourcePG=sourcePG*alphaPG

    end do
  !THERMAL NEUTRONS
!=========================================================
  !calculate subQ_TN
  !call dsolve TN
 call dsolve(sigmaTN,numDisc,subLength,sourceTN,subQ_TN,subPsi,subX,alphaTN)
    do ii=1,numSub
      psiTN((ii-1)*numDisc+1:ii*numDisc)=subPsi*sMultiply
      sourceTN=sourceTN*alphaTN

    end do
  !SECONDARY GAMMA
!=========================================================
  !calcualte subQ_SG
  subQ_SG=sigmaSG*(psiFN(1:numDisc)*factorFN+psiTN(1:numDisc)*factorTN)
  !call dsolve SG
 call dsolve(sigmaSG,numDisc,subLength,sourceSG,subQ_SG,subPsi,subX,alphaSG)
    do ii=1,numSub
      psiSG((ii-1)*numDisc+1:ii*numDisc)=subPsi*sMultiply
      sourceSG=sourceSG*alphaSG

    end do


    !unneeded
    ! do ii=1,numDisc
    !   xArray(ii)=subX(ii)
    ! end do
    !
    ! do ii=numDisc+1,total
    !   xArray(ii)=xArray(numDisc)+(ii-numDisc)*delta
    !   !print*,xArray(ii)
    ! end do
    ! do ii = 1,size(subX)
    !   !print*,subX(ii)
    ! end do
    ! do ii = 1,total
    !   !print*,xArray(ii)
    ! end do

    ! call fullsolve(numSub,alpha,xArray)
    !/endunneeded


!==================================
  call cpu_time(finish)
    print*,"Time =", finish-start, "seconds."


    open(100, file="data1.dat", status="unknown")
      do ii=1,total
       write(100,*) xArray(ii), psiArray(ii)

       end do
   close(100)

    open(100, file="dataFN.dat", status="unknown")
     do ii=1,total
        write(100,*) xArray(ii), psiFN(ii)
      end do
    close(100)

    open(100, file="dataPG.dat", status="unknown")
    do ii=1,total
       write(100,*) xArray(ii), psiPG(ii)
     end do
    close(100)

    open(100, file="dataTN.dat", status="unknown")
     do ii=1,total
        write(100,*) xArray(ii), psiTN(ii)
      end do
    close(100)

    open(100, file="dataSG.dat", status="unknown")
    do ii=1,total
       write(100,*) xArray(ii), psiSG(ii)
     end do
    close(100)



end program
