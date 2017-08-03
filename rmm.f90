
module rmm
  use constants
  use diamondsolve
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

  subroutine particleSolve(sources,numDisc,captureCS,particleCS,factors,length,subX,subPsiFN,subPsiTN,subPsiPG,subPsiSG,outputs)
!  input
    real(kind=dp),dimension(4) :: sources
    real(kind=dp),dimension(4) :: captureCS,totalCS
    real(kind=dp),dimension(4,4) :: particleCS,factors
    real(kind=dp) :: length
    integer :: numDisc,numSub
!output
    real(kind=dp),dimension(4,4) :: decayConstants
    real(kind=dp),dimension(4) :: outputs
    !temporary arrays for calculation of discretisations
    real(kind=dp),dimension(:),allocatable :: subPsiFN,subPsiTN,subPsiPG,subPsiSG,subX,subQ_FN,subQ_PG,subQ_TN,subQ_SG

    !data arrays for each particle (OUTPUT)

    real(kind=dp) :: rhsFN,rhsTN,rhsPG,rhsSG

    real(kind=dp) :: subLength,delta

    integer :: ii
!allocation
      subLength=length/numSub
      delta=subLength/numDisc

      !creating a matrix which factors in all removals of particles from system
        do ii=1,size(captureCS)
          totalCS(ii)=sum(particleCS(ii,:))+captureCS(ii)
        end do
        !TESTING PURPOSES todo remove
        !totalCS=3
    !fast neutrons
  

      allocate(subQ_FN(numDisc))
      allocate(subQ_PG(numDisc))
      allocate(subQ_TN(numDisc))
      allocate(subQ_SG(numDisc))


    subQ_FN=0.0_dp
    !want to sum first row of cross section area

    call dsolve(sources(1),totalCS(1),numDisc,length,subQ_FN,subPsiFN,subX,rhsFN)
    ! do ii=1,numDisc
    !   !print*,subPsiFN(ii),subX(ii)
    ! end do

    subQ_TN=particleCS(1,2)*subPsiFN*factors(2,1)
    call dsolve(sources(2),totalCS(2),numDisc,length,subQ_TN,subPsiTN,subX,rhsTN)
    call dsolve(sources(3),totalCS(3),numDisc,length,subQ_PG,subPsiPG,subX,rhsPG)

    subQ_SG=particleCS(1,4)*subPsiFN*factors(4,1)+particleCS(2,4)*subPsiTN*factors(4,2)+particleCS(3,4)*subPsiPG*factors(4,3)
    !print*,particleCS(3,4)
    call dsolve(sources(4),totalCS(4),numDisc,length,subQ_SG,subPsiSG,subX,rhsSG)
  !  print*,particleCS(2,4)*subPsiTN*factors(4,2)
    !real(kind=dp),dimension(:) :: xArray,subPsiFN,subPsiTN,subPsiPG,subPsiSG,subX,subQ_FN,subQ_PG,subQ_TN,subQ_SG
!         open(100, file="dataFN.dat", status="unknown")
!          do ii=1,numDisc
!             write(100,*) subX(ii), subPsiFN(ii)
!           end do
!         close(100)
!
!         open(100, file="dataPG.dat", status="unknown")
!         do ii=1,numDisc
!            write(100,*) subX(ii), subPsiPG(ii)
!          end do
!         close(100)
!
!         open(100, file="dataTN.dat", status="unknown")
!          do ii=1,numDisc
!             write(100,*) subX(ii), subPsiTN(ii)
!           end do
!         close(100)
!
!         open(100, file="dataSG.dat", status="unknown")
!          do ii=1,numDisc
!             write(100,*) subX(ii), subPsiSG(ii)
!           end do
!         close(100)
         outputs=0_dp
        ! print*,'look',rhsFN,outputs
         outputs(1)=rhsFN
         outputs(2)=rhsTN
         outputs(3)=rhsPG
         outputs(4)=rhsSG
  end subroutine

!   subroutine plot(xArray,subPsiFN,subPsiTN,subPsiPG,subPsiSG)

!
!   end subroutine
end module

  !todo fix q array
