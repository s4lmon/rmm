program solver
  use rmm
  use extrapolate
  implicit none
    real(kind=dp) :: domainLength,subLength
    integer :: numDisc, numSub, ii,total
    !sourcesR is for right hand boundary
      real(kind=dp),dimension(4) :: sourcesL,sourcesR,outputsR,outputsL
      real(kind=dp),dimension(4) :: captureCS
      real(kind=dp),dimension(4,4) :: particleCS,factors
      real(kind=dp) :: sourceFN,sourceTN, sourcePG, sourceSG
      real(kind=dp),dimension(:),allocatable :: xArray,psiFN,psiTN_LR,psiTN_RL,psiPG,psiSG_LR,psiSG_RL
      real(kind=dp),dimension(:),allocatable :: subX,subPsiFN,subPsiTN_LR,subPsiTN_RL,subPsiPG,subPsiSG_LR,subPsiSG_RL
      real(kind=dp) :: fn_tn,fn_sg,tn_sg,pg_sg
      real(kind=dp),dimension(4,4,2) :: decays
      real(kind=dp),dimension(:),allocatable :: intervalFN,intervalTN_LR,intervalTN_RL,intervalPG,intervalSG_LR,intervalSG_RL
    !storing subPsiParticle production for each input particle
    !e.g. causeFN is a four dimensional array where dimension 1 represents the FN it causes, dimension 2 is the TN etc
    real(kind=dp),dimension(:,:),allocatable :: causeFN, causeTN, causePG, causeSG
real(kind=dp) :: start,finish,omg
!numDisc=1
!do while(numDisc<=100000000)
!call cpu_time(start)
!INPUT VALUES

    domainLength=10_dp
    numDisc=10
    numSub=1
!    numDisc=1000000
!    numSub=100000000/numDisc
    !print*,numDisc,numSub
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
    sourceTN=1.0_dp
    sourcePG=1.0_dp
    sourceSG=0.0_dp
      sourcesL(1)=sourceFN
      sourcesL(2)=sourceTN
      sourcesL(3)=sourcePG
      sourcesL(4)=sourceSG
    sourcesR(1)=0.0_dp
    sourcesR(2)=0.0_dp
    sourcesR(3)=0.0_dp
    sourcesR(4)=0.0_dp
    subLength=domainLength/numSub
    particleCS=0.0_dp
    particleCS(1,2)=1.0_dp
    particleCS(1,4)=1.0_dp
    particleCS(2,4)=1.0_dp
    particleCS(3,4)=1.0_dp
    captureCS=(/1.0_dp,1.0_dp,1.0_dp,1.0_dp/)
    allocate(subPsiFN(numDisc))
    allocate(subPsiTN_LR(numDisc))
    allocate(subPsiTN_RL(numDisc))
    allocate(subPsiPG(numDisc))
    allocate(subPsiSG_LR(numDisc))
    allocate(subPsiSG_RL(numDisc))
    allocate(subX(numDisc))
!causeFN is array of all particles caused by FN only (can cause FN, TN, PG, SG, TN_RL,SG_RL
    allocate(causeFN(6,numDisc))
    allocate(causeTN(6,numDisc))
    allocate(causePG(6,numDisc))
    allocate(causeSG(6,numDisc))
    allocate(intervalFN(numSub))
    allocate(intervalTN_LR(numSub))
    allocate(intervalTN_RL(numSub+1))
    allocate(intervalPG(numSub))
    allocate(intervalSG_LR(numSub))
    allocate(intervalSG_RL(numSub+1))
    !call plot()
    !print*,subX
    !building alpha and beta matrix
    sourcesL=0.0_dp
    sourcesR=0.0_dp
    sourcesL(1)=1.0_dp
!=====================================================================================
!building decay constants (alpha,beta)
!switching on and off different particles and seeing what scattering they cause
    decays=0.0_dp


   call particleSolve(sourcesL,sourcesR,numDisc,captureCS,particleCS,factors,subLength,subX,subPsiFN,subPsiTN_LR,subPsiTN_RL,subPsiPG,subPsiSG_LR,subPsiSG_RL,outputsR,outputsL)
    decays(1,1,1)=outputsR(1)
    decays(1,2,1)=outputsR(2)
    decays(1,4,1)=outputsR(4)

    decays(1,2,2)=outputsL(2)
    decays(1,4,2)=outputsL(4)

    causeFN(1,:)=subPsiFN
    causeFN(2,:)=subPsiTN_LR
    causeFN(3,:)=subPsiPG
    causeFN(4,:)=subPsiSG_LR
    causeFN(5,:)=subPsiTN_RL
    causeFN(6,:)=subPsiSG_RL
    sourcesL=0.0_dp
    sourcesL(2)=1.0_dp
   call particleSolve(sourcesL,sourcesR,numDisc,captureCS,particleCS,factors,subLength,subX,subPsiFN,subPsiTN_LR,subPsiTN_RL,subPsiPG,subPsiSG_LR,subPsiSG_RL,outputsR,outputsL)
    decays(2,2,1)=outputsR(2)
    decays(2,4,1)=outputsR(4)

    decays(2,4,2)=outputsL(4)

    causeTN(1,:)=subPsiFN
    causeTN(2,:)=subPsiTN_LR
    causeTN(3,:)=subPsiPG
    causeTN(4,:)=subPsiSG_LR
    causeTN(5,:)=subPsiTN_RL
    causeTN(6,:)=subPsiSG_RL

    sourcesL=0.0_dp
    sourcesL(3)=1.0_dp
   call particleSolve(sourcesL,sourcesR,numDisc,captureCS,particleCS,factors,subLength,subX,subPsiFN,subPsiTN_LR,subPsiTN_RL,subPsiPG,subPsiSG_LR,subPsiSG_RL,outputsR,outputsL)
    decays(3,3,1)=outputsR(3)
    decays(3,4,1)=outputsR(4)

    decays(3,4,2)=outputsL(4)

    causePG(1,:)=subPsiFN
    causePG(2,:)=subPsiTN_LR
    causePG(3,:)=subPsiPG
    causePG(4,:)=subPsiSG_LR
    causePG(5,:)=subPsiTN_RL
    causePG(6,:)=subPsiSG_RL

    sourcesL=0.0_dp
    sourcesL(4)=1.0_dp
   call particleSolve(sourcesL,sourcesR,numDisc,captureCS,particleCS,factors,subLength,subX,subPsiFN,subPsiTN_LR,subPsiTN_RL,subPsiPG,subPsiSG_LR,subPsiSG_RL,outputsR,outputsL)
    decays(4,4,1)=outputsR(4)


    causeSG(1,:)=subPsiFN
    causeSG(2,:)=subPsiTN_LR
    causeSG(3,:)=subPsiPG
    causeSG(4,:)=subPsiSG_LR
    causeSG(5,:)=subPsiTN_RL
    causeSG(6,:)=subPsiSG_RL
    !==========================================
    sourcesL(1)=sourceFN
    sourcesL(2)=sourceTN
    sourcesL(3)=sourcePG
    sourcesL(4)=sourceSG
   call particleSolve(sourcesL,sourcesR,numDisc,captureCS,particleCS,factors,subLength,subX,subPsiFN,subPsiTN_LR,subPsiTN_RL,subPsiPG,subPsiSG_LR,subPsiSG_RL,outputsR,outputsL)

    !todo uncomment till end
  allocate(xArray(total))


    do ii=1,numSub

          xArray((ii-1)*numDisc+1:ii*numDisc)=subX+subLength*(ii-1)
    end do

    do ii=1,numDisc
          !print*,'hello',causeFN(:,ii)
    end do
  call extra_bound(domainLength,numDisc,numSub,sourcesL,decays,intervalFN,intervalTN_LR,intervalTN_RL,intervalPG,intervalSG_LR,intervalSG_RL)

  call extra_psi(domainLength,numDisc,numSub,xArray,intervalFN,intervalTN_LR,intervalTN_RL,intervalPG,intervalSG_LR,intervalSG_RL,causeFN,causeTN,causePG,causeSG,psiFN,psiTN_LR,psiTN_RL,psiPG,psiSG_LR,psiSG_RL)
!    open(100, file="psiFNdata.dat", status="unknown")
!      do ii=1,total
!        write(100,*) xArray(ii), psiFN(ii)
!      end do
!    close(100)
!    open(100, file="psiTNdata.dat", status="unknown")
!      do ii=1,total
!        write(100,*) xArray(ii), psiTN(ii)
!      end do
!    close(100)
!    open(100, file="psiPGdata.dat", status="unknown")
!      do ii=1,total
!        write(100,*) xArray(ii), psiPG(ii)
!      end do
!    close(100)
!    open(100, file="psiSGdata.dat", status="unknown")
!      do ii=1,total
!        write(100,*) xArray(ii), psiSG(ii)
!      end do
!    close(100)



















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
