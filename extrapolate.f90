module extrapolate
  use rmm
  implicit none



contains
  subroutine  extra_bound(domainLength,numDisc,numSub,sources,decays,intervalFN,intervalTN,intervalPG,intervalSG)

    integer :: numSub,numDisc,total,ii
    real(kind=dp) :: domainLength,writeLoop
    real(kind=dp),dimension(4,4) :: decays
    real(kind=dp),dimension(:) :: sources
    real(kind=dp),dimension(4) :: multiplier !multiplier used to represent decay at each subinterval

    !interval arrays for calculating values on boundaries
    real(kind=dp),dimension(:) :: intervalFN,intervalTN,intervalPG,intervalSG


multiplier=1.0_dp
    !starts at 1 as assigning subPsi into Psi for first interval represents no decay
    !could omit variable by assigning subPsi into Psi before calculating decay

    !psiFN(1:numDisc)=subPsiFN

    !psiPG(1:numDisc)=subPsiPG
    !psiSG(1:numDisc)=subPsiSG
    !build second stage of FN
    !print*,'dec',decays(1,2)
      intervalFN=0.0_dp
      intervalTN=0.0_dp
      intervalPG=0.0_dp
      intervalSG=0.0_dp
!      intervalFN(1,1)=sources(1)
!      intervalTN(2,1)=sources(2)
!      intervalPG(3,1)=sources(3)
!      intervalSG(4,1)=sources(4)
       intervalFN(1)=sources(1)
       intervalTN(1)=sources(2)
       intervalPG(1)=sources(3)
       intervalSG(1)=sources(4)
      !interval__ calculates the values on the boundary for each subinterval as a result of scattering into that __
        !e.g. intervalSG(1,:) is the secondary gamma left hand boundaries caused by fast neutrons, intervalSG(2,:) is caused by
        !thermal neutrons etc)
      do ii=2,numSub
      !  intervalFN(1,ii)=intervalFN(1,ii-1)*decays(1,1)
         intervalFN(ii)=intervalFN(ii-1)*decays(1,1)
      !  intervalPG(3,ii)=intervalPG(3,ii-1)*decays(3,3)
         intervalPG(ii)=intervalPG(ii-1)*decays(3,3)
      !  intervalTN(2,ii)=intervalTN(2,ii-1)*decays(2,2)
      !  intervalTN(1,ii)=intervalFN(1,ii-1)*decays(1,2)
         intervalTN(ii)=intervalTN(ii-1)*decays(2,2)+intervalFN(ii-1)*decays(1,2)
      !  intervalSG(4,ii)=intervalSG(4,ii-1)*decays(4,4)
      !  intervalSG(1,ii)=intervalFN(1,ii-1)*decays(1,4)
      !  intervalSG(2,ii)=intervalTN(2,ii-1)*decays(2,4)
      !  intervalSG(3,ii)=intervalPG(3,ii-1)*decays(3,4)
         intervalSG(ii)=intervalSG(ii-1)*decays(4,4)+intervalFN(ii-1)*decays(1,4)+intervalTN(ii-1)*decays(2,4)+intervalPG(ii-1)*decays(3,4)
      end do
!    open(100, file="intFN.dat", status="unknown")
!        writeLoop=0.0_dp
!      do ii=1,(numSub)
!
!        write(100,*) writeLoop,intervalSG(ii)
!        writeLoop=writeLoop+domainLength/numSub
!      end do
!    close(100)
!

!    do ii=1,numSub-1
!      !calculating FN, simple as no source term through domain
!      psiFN((ii-1)*numDisc+1:ii*numDisc)=subPsiFN*multiplier(1)
!      !print*,psiFN((ii-1)*numDisc+1:ii*numDisc)
!      multiplier(1)=sources(1)*decays(1,1)*multiplier(1)
!!=====================================================================
!
!      !calculating PG, simple as no source term through domain
!      psiPG((ii-1)*numDisc+1:ii*numDisc)=subPsiPG*multiplier(3)
!      multiplier(3)=sources(1)*decays(3,3)*multiplier(3)
!!=====================================================================
!
!      !calculating TN, need to take into account production
!      !psiTN((ii)*numDisc+1)=decays(1,2)*psiFN((ii-1)*numDisc+1)+decays(2,2)*psiTN((ii-1)*numDisc+1)
!
!      print*,'x',psiTN(ii)
!    end do
!
!
!    open(100, file="test1.dat", status="unknown")
!      do ii=1,total
!        write(100,*) xArray(ii), psiFN(ii)
!      end do
!    close(100)
!
!    open(100, file="test2.dat", status="unknown")
!      do ii=1,total
!        write(100,*) xArray(ii), psiPG(ii)
!      end do
!    close(100)
!    !build second stage of PG
!
!    !build second stage of Tn
!
!    !build second stage of SG
!    open(100, file="test3.dat", status="unknown")
!      do ii=1,total
!        write(100,*) xArray(ii), psiTN(ii)
!      end do
!    close(100)
!
!
!    do ii=1,numSub
!      !print*,psiFN(ii)
!    end do
  end subroutine

  subroutine extra_psi(domainLength,numDisc,numSub,xArray,intervalFN,intervalTN,intervalPG,intervalSG,causeFN,causeTN,causePG,causeSG,psiFN,psiTN,psiPG,psiSG)

    real(kind=dp) :: domainLength
    real(kind=dp),dimension(:,:) :: causeFN, causeTN, causePG, causeSG
    integer :: numSub,numDisc,total,ii
    real(kind=dp),dimension(:),allocatable :: xArray,psiFN,psiTN,psiPG,psiSG
!    real(kind=dp),dimension(:) :: subPsiFN,subPsiTN,subPsiPG,subPsiSG
    real(kind=dp),dimension(:) :: intervalFN,intervalTN,intervalPG,intervalSG
    total = numSub*numDisc

    allocate(psiFN(total))
    allocate(psiTN(total))
    allocate(psiPG(total))
    allocate(psiSG(total))
      psiFN(1:numDisc)=causeFN(1,:)*intervalFN(1)
      psiPG(1:numDisc)=causePG(3,:)*intervalPG(1)
      psiTN(1:numDisc)=causeTN(2,:)*intervalTN(1)+causeFN(2,:)*intervalFN(1)
      psiSG(1:numDisc)=causeSG(4,:)*intervalSG(1)+causeFN(4,:)*intervalFN(1)+causeTN(4,:)*intervalTN(1)+causePG(4,:)*intervalPG(1)
    !print*, numDisc-1
    do ii=2,numSub
      psiFN((ii-1)*numDisc+1:ii*numDisc)=causeFN(1,:)*intervalFN(ii)
      psiPG((ii-1)*numDisc+1:ii*numDisc)=causePG(3,:)*intervalPG(ii)
      psiTN((ii-1)*numDisc+1:ii*numDisc)=causeTN(2,:)*intervalTN(ii)+causeFN(2,:)*intervalFN(ii)
      psiSG((ii-1)*numDisc+1:ii*numDisc)=causeSG(4,:)*intervalSG(ii)+causeFN(4,:)*intervalFN(ii)+causeTN(4,:)*intervalTN(ii)+causePG(4,:)*intervalPG(ii)
    end do
    do ii=1,total
      !print*,'x',psiFN(ii)
    end do 
 !   open(100, file="psiFNtest.dat", status="unknown")
 !     do ii=1,total
 !       write(100,*) xArray(ii), psiTN(ii)
 !     end do
 !   close(100)
  end subroutine
end module
