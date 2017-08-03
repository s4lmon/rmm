program solver
  use rmm
  use extrapolate
  implicit none
    real(kind=dp) :: domainLength,subLength
    integer :: numDisc, numSub, ii,total
      real(kind=dp),dimension(4) :: sources,outputs
      real(kind=dp),dimension(4) :: captureCS
      real(kind=dp),dimension(4,4) :: particleCS,factors
      real(kind=dp) :: sourceFN,sourceTN, sourcePG, sourceSG
      real(kind=dp),dimension(:),allocatable :: xArray,psiFN,psiTN,psiPG,psiSG
      real(kind=dp),dimension(:),allocatable :: subX,subPsiFN,subPsiTN,subPsiPG,subPsiSG
      real(kind=dp) :: fn_tn,fn_sg,tn_sg,pg_sg
      real(kind=dp),dimension(4,4) :: decays

      real(kind=dp),dimension(:),allocatable :: intervalFN,intervalTN,intervalPG,intervalSG

    !storing subPsiParticle production for each input particle
    !e.g. causeFN is a four dimensional array where dimension 1 represents the FN it causes, dimension 2 is the TN etc
    real(kind=dp),dimension(:,:),allocatable :: causeFN, causeTN, causePG, causeSG
real(kind=dp) :: start,finish,omg
!numDisc=1
!do while(numDisc<=100000000)    
!call cpu_time(start)
!INPUT VALUES
 
    domainLength=10_dp
    numDisc=10000
    numSub=10
!    numDisc=1000000
!    numSub=100000000/numDisc
    print*,numDisc,numSub
    total=numSub*numDisc
    !cal discretiser(domainLength,numSub,numDisc,xArray)
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
    sourceFN=1.0_dp
    sourceTN=0.0_dp
    sourcePG=1.0_dp
    sourceSG=0.0_dp
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
    allocate(subPsiFN(numDisc))
    allocate(subPsiTN(numDisc))
    allocate(subPsiPG(numDisc))
    allocate(subPsiSG(numDisc))
    allocate(subX(numDisc))

    allocate(causeFN(4,numDisc))
    allocate(causeTN(4,numDisc))
    allocate(causePG(4,numDisc))
    allocate(causeSG(4,numDisc))
    allocate(intervalFN(numSub))
    allocate(intervalTN(numSub))
    allocate(intervalPG(numSub))
    allocate(intervalSG(numSub))
    !call plot()
    !print*,subX
    !building alpha and beta matrix
    sources=0.0_dp
    sources(1)=1.0_dp
!=====================================================================================
!building decay constants (alpha,beta)
!switching on and off different particles and seeing what scattering they cause
    decays=0.0_dp


    call particleSolve(sources,numDisc,captureCS,particleCS,factors,subLength,subX,subPsiFN,subPsiTN,subPsiPG,subPsiSG,outputs)
    decays(1,1)=outputs(1)
    decays(1,2)=outputs(2)
    decays(1,4)=outputs(4)

    causeFN(1,:)=subPsiFN
    causeFN(2,:)=subPsiTN
    causeFN(3,:)=subPsiPG
    causeFN(4,:)=subPsiSG

    sources=0.0_dp
    sources(2)=1.0_dp
    call particleSolve(sources,numDisc,captureCS,particleCS,factors,subLength,subX,subPsiFN,subPsiTN,subPsiPG,subPsiSG,outputs)
    decays(2,2)=outputs(2)
    decays(2,4)=outputs(4)

    causeTN(1,:)=subPsiFN
    causeTN(2,:)=subPsiTN
    causeTN(3,:)=subPsiPG
    causeTN(4,:)=subPsiSG

    sources=0.0_dp
    sources(3)=1.0_dp
    call particleSolve(sources,numDisc,captureCS,particleCS,factors,subLength,subX,subPsiFN,subPsiTN,subPsiPG,subPsiSG,outputs)
    decays(3,3)=outputs(3)
    decays(3,4)=outputs(4)

    causePG(1,:)=subPsiFN
    causePG(2,:)=subPsiTN
    causePG(3,:)=subPsiPG
    causePG(4,:)=subPsiSG

    sources=0.0_dp
    sources(4)=1.0_dp
    call particleSolve(sources,numDisc,captureCS,particleCS,factors,subLength,subX,subPsiFN,subPsiTN,subPsiPG,subPsiSG,outputs)
    decays(4,4)=outputs(4)

    causeSG(1,:)=subPsiFN
    causeSG(2,:)=subPsiTN
    causeSG(3,:)=subPsiPG
    causeSG(4,:)=subPsiSG
    !==========================================
    sources(1)=sourceFN
    sources(2)=sourceTN
    sources(3)=sourcePG
    sources(4)=sourceSG

    call particleSolve(sources,numDisc,captureCS,particleCS,factors,subLength,subX,subPsiFN,subPsiTN,subPsiPG,subPsiSG,outputs)


  allocate(xArray(total))


    do ii=1,numSub

          xArray((ii-1)*numDisc+1:ii*numDisc)=subX+subLength*(ii-1)
    end do

    do ii=1,numDisc
          !print*,'hello',causeFN(:,ii)
    end do
  call extra_bound(domainLength,numDisc,numSub,sources,decays,intervalFN,intervalTN,intervalPG,intervalSG)

  call extra_psi(domainLength,numDisc,numSub,xArray,intervalFN,intervalTN,intervalPG,intervalSG,causeFN,causeTN,causePG,causeSG,psiFN,psiTN,psiPG,psiSG)
    open(100, file="psiFNdata.dat", status="unknown")
      do ii=1,total
        write(100,*) xArray(ii), psiFN(ii)
      end do
    close(100)
    open(100, file="psiTNdata.dat", status="unknown")
      do ii=1,total
        write(100,*) xArray(ii), psiTN(ii)
      end do
    close(100) 
    open(100, file="psiPGdata.dat", status="unknown")
      do ii=1,total
        write(100,*) xArray(ii), psiPG(ii)
      end do
    close(100)
    open(100, file="psiSGdata.dat", status="unknown")
      do ii=1,total
        write(100,*) xArray(ii), psiSG(ii)
      end do
    close(100) 



















!      call cpu_time(finish)
!          print*,"Time =", finish-start, "seconds."
  !  open(100, file="runtime.dat", status="unknown")
  !      write(100,*) numDisc,finish-start
  !  close(100)
!    numDisc=numDisc*10
!    deallocate(subPsiFN)
!    deallocate(subPsiTN)
!    deallocate(subPsiPG)
!    deallocate(subPsiSG)
!    deallocate(subX)
!
!    deallocate(causeFN)
!    deallocate(causeTN)
!    deallocate(causePG)
!    deallocate(causeSG)
!    deallocate(intervalFN)
!    deallocate(intervalTN)
!    deallocate(intervalPG)
!    deallocate(intervalSG)
!    deallocate(xArray)
!    deallocate(psiFN)
!    deallocate(psiTN)
!    deallocate(psiPG)
!    deallocate(psiSG)
!! end do  
end program
