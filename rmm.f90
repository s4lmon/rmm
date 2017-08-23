
module rmm
  use constants
  use diamondsolve
  use backdiamond
  implicit none
contains
  subroutine discretiser(domainLength,numSub,numDisc,xArray)

    real(kind=dp) :: subLength,domainLength
    real(kind=dp) :: delta
    integer :: numDisc,numSub,total
    real(kind=dp),dimension(:),allocatable :: xArray
    subLength=domainLength/numSub
    delta=subLength/numDisc
    total=numDisc*numSub
    allocate(xArray(total))

    !print*,numDisc
  end subroutine
!fast neutrons and primary gamma cannot be created so only travel left to right
  subroutine particleSolve(sourcesL,sourcesR,numDisc,captureCS,particleCS,factors,length,subX,subPsiFN,subPsiTN_LR,subPsiTN_RL,subPsiPG,subPsiSG_LR,subPsiSG_RL,outputsR,outputsL)
!  input
    real(kind=dp),dimension(4) :: sourcesL,sourcesR
    real(kind=dp),dimension(4) :: captureCS,totalCS
    real(kind=dp),dimension(4,4) :: particleCS,factors
    real(kind=dp) :: length
    integer :: numDisc,numSub
!output
    real(kind=dp),dimension(4,4,2) :: decayConstants
    real(kind=dp),dimension(4) :: outputsR,outputsL
    !temporary arrays for calculation of discretisations
    real(kind=dp),dimension(:),allocatable :: subPsiFN,subPsiTN_LR,subPsiTN_RL,subPsiPG,subPsiSG_LR,subPsiSG_RL,subX,subQ_FN,subQ_PG,subQ_TN_LR,subQ_TN_RL,subQ_SG_LR,subQ_SG_RL
!split for particles going left vs right
    real(kind=dp) :: TNsplit,SGsplit
    !data arrays for each particle (OUTPUT)
    !output values: note only lhs for TN and SG, they travel right to left
    real(kind=dp) :: rhsFN,rhsTN,rhsPG,rhsSG,lhsTN,lhsSG

    real(kind=dp) :: subLength,delta

    integer :: ii
!allocation
!print*,'crosssec ',captureCS(1),particleCS(1,:)
      subLength=length/numSub
      delta=subLength/numDisc

      !creating a matrix which factors in all removals of particles from system
        do ii=1,size(captureCS)
          totalCS(ii)=sum(particleCS(ii,:))+captureCS(ii)
        end do
        !TESTING PURPOSES todo remove
        !totalCS=3


      allocate(subQ_FN(numDisc))
      allocate(subQ_PG(numDisc))
      allocate(subQ_TN_LR(numDisc))
      allocate(subQ_TN_RL(numDisc))
      allocate(subQ_SG_LR(numDisc))
      allocate(subQ_SG_RL(numDisc))
      !split is particles going right, so particles going left is (1-split)
        TNsplit=0.5_dp
        SGsplit=0.5_dp
!===================================
    !fast neutrons
    subQ_FN=0.0_dp
    !want to sum first row of cross section area
    !print*,'look',sourcesL(1),totalCS(1),numDisc,length,subQ_FN,subPsiFN,subX,rhsFN
    call dsolve(sourcesL(1),totalCS(1),numDisc,length,subQ_FN,subPsiFN,subX,rhsFN)

    ! do ii=1,numDisc
    !    print*,subPsiFN(ii),subX(ii)
    ! end do
!==================================
    !thermal neutrons
    subQ_TN_LR=TNsplit*particleCS(1,2)*subPsiFN*factors(2,1)
    subQ_TN_RL=(1-TNsplit)*particleCS(1,2)*subPsiFN*factors(2,1)

    call dsolve(sourcesL(2),totalCS(2),numDisc,length,subQ_TN_LR,subPsiTN_LR,subX,rhsTN)
    call backsolve(sourcesR(2),totalCS(2),numDisc,length,subQ_TN_RL,subPsiTN_RL,subX,lhsTN)
!=================================
    !prompt gamma
    call dsolve(sourcesL(3),totalCS(3),numDisc,length,subQ_PG,subPsiPG,subX,rhsPG)


!================================
    !secondary gamma
    subQ_SG_LR=SGsplit*particleCS(1,4)*subPsiFN*factors(4,1)+particleCS(2,4)*(subPsiTN_LR+subPsiTN_RL)*factors(4,2)+particleCS(3,4)*subPsiPG*factors(4,3)
    subQ_SG_RL=(1-SGsplit)*particleCS(1,4)*subPsiFN*factors(4,1)+particleCS(2,4)*(subPsiTN_LR+subPsiTN_RL)*factors(4,2)+particleCS(3,4)*subPsiPG*factors(4,3)
    !print*,particleCS(3,4)
    call dsolve(sourcesL(4),totalCS(4),numDisc,length,subQ_SG_LR,subPsiSG_LR,subX,rhsSG)
    call backsolve(sourcesR(4),totalCS(4),numDisc,length,subQ_SG_LR,subPsiSG_RL,subX,lhsSG)
  !  print*,particleCS(2,4)*subPsiTN*factors(4,2)
    !real(kind=dp),dimension(:) :: xArray,subPsiFN,subPsiTN,subPsiPG,subPsiSG,subX,subQ_FN,subQ_PG,subQ_TN,subQ_SG
         ! open(100, file="dataFN.dat", status="unknown")
         !  do ii=1,numDisc
         !     write(100,*) subX(ii), subPsiFN(ii)
         !   end do
         ! close(100)

         ! open(100, file="dataPG.dat", status="unknown")
         ! do ii=1,numDisc
         !    write(100,*) subX(ii), subPsiPG(ii)
         !  end do
         ! close(100)

         ! open(100, file="dataTN_LR.dat", status="unknown")
         !  do ii=1,numDisc
         !     write(100,*) subX(ii), subPsiTN_LR(ii)
         !   end do
         ! close(100)
         ! open(100, file="dataTN_RL.dat", status="unknown")
         !  do ii=1,numDisc
         !     write(100,*) subX(ii), subPsiTN_RL(ii)
         !   end do
         ! close(100)

         ! open(100, file="dataSG_LR.dat", status="unknown")
         !  do ii=1,numDisc
         !     write(100,*) subX(ii), subPsiSG_LR(ii)
         !     !print*,subX(ii)
         !   end do
         ! close(100)
         ! open(100, file="dataSG_RL.dat", status="unknown")
         !  do ii=1,numDisc
         !     write(100,*) subX(ii), subPsiSG_RL(ii)
         !   end do
         ! close(100)
         outputsR=0_dp
        ! print*,'look',rhsFN,outputsR
         outputsR(1)=rhsFN
         outputsR(2)=rhsTN
         outputsR(3)=rhsPG
         outputsR(4)=rhsSG

        !outputs for PG and FN explicit, no need for calculation
         outputsL(1)=0.0_dp
         outputsL(2)=lhsTN
         outputsL(3)=0.0_dp
         outputsL(4)=lhsSG
         !print*,outputsR
  end subroutine

!   subroutine plot(xArray,subPsiFN,subPsiTN,subPsiPG,subPsiSG)

!
!   end subroutine
end module

  !todo fix q array
