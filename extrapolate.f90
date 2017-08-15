module extrapolate
  use rmm
  implicit none



contains
  subroutine extra_bound(domainLength,numDisc,numSub,sources,decays,intervalFN,intervalTN_LR,intervalTN_RL,intervalPG,intervalSG_LR,intervalSG_RL)

    integer :: numSub,numDisc,total,ii
    real(kind=dp) :: domainLength,writeLoop
    real(kind=dp),dimension(4,4,2) :: decays
    real(kind=dp),dimension(:) :: sources
    real(kind=dp),dimension(4) :: multiplier !multiplier used to represent decay at each subinterval

    !interval arrays for calculating values on boundaries
    real(kind=dp),dimension(:) :: intervalFN,intervalTN_LR,intervalTN_RL,intervalPG,intervalSG_LR,intervalSG_RL


multiplier=1.0_dp
    !starts at 1 as assigning subPsi into Psi for first interval represents no decay
    !could omit variable by assigning subPsi into Psi before calculating decay

    !psiFN(1:numDisc)=subPsiFN

    !psiPG(1:numDisc)=subPsiPG
    !psiSG(1:numDisc)=subPsiSG
    !build second stage of FN
    !print*,'dec',decays(1,2)
      intervalFN=0.0_dp
      intervalTN_LR=0.0_dp
      intervalTN_RL=0.0_dp
      intervalPG=0.0_dp
      intervalSG_LR=0.0_dp
      intervalSG_RL=0.0_dp
       intervalFN(1)=sources(1)
       intervalTN_LR(1)=sources(2)
       intervalPG(1)=sources(3)
       intervalSG_LR(1)=sources(4)
      !interval__ calculates the values on the boundary for each subinterval as a result of scattering into that __
        !e.g. intervalSG(1,:) is the secondary gamma left hand boundaries caused by fast neutrons, intervalSG(2,:) is caused by
        !thermal neutrons etc)
      do ii=2,numSub
         intervalFN(ii)=intervalFN(ii-1)*decays(1,1,1)
         intervalPG(ii)=intervalPG(ii-1)*decays(3,3,1)
         intervalTN_LR(ii)=intervalTN_LR(ii-1)*decays(2,2,1)+intervalFN(ii-1)*decays(1,2,1)
       end do

      intervalTN_RL(numSub+1)=0.0_dp
      intervalSG_RL(numSub+1)=0.0_dp

      do ii=numSub,2,-1
         intervalTN_RL(ii)=decays(2,2,1)*intervalTN_RL(ii+1)+decays(1,2,2)*intervalFN(ii)
         intervalSG_RL(ii)=decays(1,4,2)*intervalFN(ii)+decays(2,4,2)*intervalTN_LR(ii)+decays(3,4,2)*intervalPG(ii)+decays(4,4,1)*intervalSG_RL(ii+1)
         !intervalTN_RL(ii)=
!         print*,intervalTN_LR(ii),intervalTN_RL(ii)
!         print*,decays(1,2,1)
      end do

      do ii=2,numSub
        intervalSG_LR(ii)=intervalSG_LR(ii-1)*decays(4,4,1)+intervalFN(ii-1)*decays(1,4,1)+intervalTN_LR(ii-1)*decays(2,4,1)+intervalPG(ii-1)*decays(3,4,1)+intervalTN_RL(ii)*decays(2,4,2)
      end do

  !=============================
      print*,'look',decays(2,4,2),decays(2,4,1)
   open(100, file="intFN.dat", status="unknown")
       writeLoop=0.0_dp
     do ii=1,(numSub)

       write(100,*) writeLoop,intervalFN(ii)
       writeLoop=writeLoop+domainLength/numSub
     end do
   close(100)

   open(100, file="intTN_LR.dat", status="unknown")
       writeLoop=0.0_dp
     do ii=1,(numSub)

       write(100,*) writeLoop,intervalTN_LR(ii)
       writeLoop=writeLoop+domainLength/numSub
     end do
   close(100)

   open(100, file="intTN_RL.dat", status="unknown")
       writeLoop=0.0_dp
     do ii=1,(numSub)

       write(100,*) writeLoop,intervalTN_RL(ii)
       writeLoop=writeLoop+domainLength/numSub
     end do
   close(100)

   open(100, file="intSG_LR.dat", status="unknown")
       writeLoop=0.0_dp
     do ii=1,(numSub)

       write(100,*) writeLoop,intervalSG_LR(ii)
       writeLoop=writeLoop+domainLength/numSub
     end do
   close(100)

   open(100, file="intSG_RL.dat", status="unknown")
       writeLoop=0.0_dp
     do ii=1,(numSub)

       write(100,*) writeLoop,intervalSG_RL(ii)
       writeLoop=writeLoop+domainLength/numSub
     end do
   close(100)

   open(100, file="intPG.dat", status="unknown")
       writeLoop=0.0_dp
     do ii=1,(numSub)

       write(100,*) writeLoop,intervalPG(ii)
       writeLoop=writeLoop+domainLength/numSub
     end do
   close(100)

!===============================================================

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

  subroutine extra_psi(domainLength,numDisc,numSub,xArray,intervalFN,intervalTN_LR,intervalTN_RL,intervalPG,intervalSG_LR,intervalSG_RL,causeFN,causeTN,causePG,causeSG,psiFN,psiTN_LR,psiTN_RL,psiPG,psiSG_LR,psiSG_RL)

    real(kind=dp) :: domainLength
    real(kind=dp),dimension(:,:) :: causeFN, causeTN, causePG, causeSG
    real(kind=dp),dimension(:,:),allocatable :: tempFN,tempTN,tempPG,tempSG
    integer :: numSub,numDisc,total,ii,jj
    real(kind=dp),dimension(:),allocatable :: xArray,psiFN,psiTN_LR,psiTN_RL,psiPG,psiSG_LR,psiSG_RL
!    real(kind=dp),dimension(:) :: subPsiFN,subPsiTN,subPsiPG,subPsiSG
    real(kind=dp),dimension(:) :: intervalFN,intervalTN_LR,intervalTN_RL,intervalPG,intervalSG_LR,intervalSG_RL
    total = numSub*numDisc

    allocate(tempTN(6,numDisc))
    allocate(tempFN(6,numDisc))
    allocate(tempPG(6,numDisc))
    allocate(tempSG(6,numDisc))
    allocate(psiFN(total))
    allocate(psiTN_LR(total))
    allocate(psiTN_RL(total))
    allocate(psiPG(total))
    allocate(psiSG_LR(total))
    allocate(psiSG_RL(total))

      !psiFN(1:numDisc)=causeFN(1,:)*intervalFN(1)
      !====
      ! psiPG(1:numDisc)=causePG(3,:)*intervalPG(1)
      ! psiTN_LR(1:numDisc)=causeTN(2,:)*intervalTN_LR(1)+causeFN(2,:)*intervalFN(1)
      ! psiSG_LR(1:numDisc)=causeSG(4,:)*intervalSG_LR(1)+causeFN(4,:)*intervalFN(1)+causeTN(4,:)*intervalTN_LR(1)+causePG(4,:)*intervalPG(1)+causeTN(6,:)*intervalTN_RL(2)
      !====
 !next two do loops build reverse TN then reverse SG

    !print*,tempTN
  !  print*,size(psiTN_RL(1:numDisc)),size(intervalTN_RL)


    do ii=1,numDisc
        tempFN(5,((numDisc+1)-ii))=causeFN(5,ii)
        tempTN(2,((numDisc+1)-ii))=causeTN(2,ii)
        tempFN(6,((numDisc+1)-ii))=causeFN(6,ii)
        !TN causing SG going right to left
        tempTN(6,((numDisc+1)-ii))=causeTN(6,ii)
        tempPG(6,((numDisc+1)-ii))=causePG(6,ii)
        tempSG(4,((numDisc+1)-ii))=causeSG(4,ii)
        !TN causing SG going left to right
        tempTN(4,((numDisc+1)-ii))=causeTN(4,ii)
    end do
    print*,causeFN(5,:)
    open(100, file="tempFN.dat", status="unknown")
     do ii=1,numDisc
       write(100,*) xArray(ii), tempFN(5,ii)
     end do
    close(100)

          ! psiTN_RL(1:numDisc)=tempTN(2,:)*intervalTN_RL(2)+tempFN(5,:)*intervalFN(2)
          !
          ! psiSG_RL(1:numDisc)=tempFN(6,:)*intervalFN(1)+tempTN(6,:)*intervalTN_LR(2)+tempTN(4,:)*intervalTN_RL(1)+tempSG(4,:)*intervalSG_RL(1)
      !print*,psiTN_RL
    !  psiTN_RL(1:numDisc)=intervalTN_RL(1)*causes
    !print*, numDisc-1
    do ii=1,numSub

      psiFN((ii-1)*numDisc+1:ii*numDisc)=causeFN(1,:)*intervalFN(ii)
      psiPG((ii-1)*numDisc+1:ii*numDisc)=causePG(3,:)*intervalPG(ii)
      psiSG_LR((ii-1)*numDisc+1:ii*numDisc)=causeSG(4,:)*intervalSG_LR(ii)+causeFN(4,:)*intervalFN(ii)+causeTN(4,:)*intervalTN_LR(ii)+causePG(4,:)*intervalPG(ii)+causeTN(6,:)*intervalTN_RL(ii+1)

      psiTN_RL((ii-1)*numDisc+1:ii*numDisc)=tempTN(2,:)*intervalTN_RL(ii+1)+causeFN(5,:)*intervalFN(ii)
      psiSG_RL((ii-1)*numDisc+1:ii*numDisc)=causeFN(6,:)*intervalFN(ii)+causeTN(6,:)*intervalTN_LR(ii)+tempTN(4,:)*intervalTN_RL(ii+1)+tempSG(4,:)*intervalSG_RL(ii+1)+causePG(6,:)*intervalPG(ii)

      psiTN_LR((ii-1)*numDisc+1:ii*numDisc)=causeTN(2,:)*intervalTN_LR(ii)+causeFN(2,:)*intervalFN(ii)
      !psiSG_LR((ii-1)*numDisc+1:ii*numDisc)=causeSG(4,:)*intervalSG_LR(ii)+causeFN(4,:)*intervalFN(ii)+causeTN(4,:)*intervalTN_LR(ii)+causePG(4,:)*intervalPG(ii)+

      !todo

    end do
    do ii=1,total
      !print*,psiTN_RL(ii)
    end do
    do ii=numSub-1,1
    !  psiSG_RL(ii*numDisc:(ii-1)*numDisc)=
    end do
    do ii=1,total
      !print*,'x',psiFN(ii)
    end do
    open(100, file="psiTN_RLtest2.dat", status="unknown")
     do ii=1,total
       write(100,*) xArray(ii), psiTN_RL(ii)
     end do
   close(100)
   open(100, file="psiTN_LRtest2.dat", status="unknown")
    do ii=1,total
      write(100,*) xArray(ii), psiTN_LR(ii)
    end do
  close(100)
  open(100, file="psiFNtest2.dat", status="unknown")
   do ii=1,total
     write(100,*) xArray(ii), psiFN(ii)
   end do
 close(100)
 open(100, file="psiPGtest2.dat", status="unknown")
  do ii=1,total
    write(100,*) xArray(ii), psiPG(ii)
  end do
close(100)
open(100, file="psiSG_LRtest2.dat", status="unknown")
 do ii=1,total
   write(100,*) xArray(ii), psiSG_LR(ii)
 end do
close(100)
open(100, file="psiSG_RLtest2.dat", status="unknown")
 do ii=1,total
   write(100,*) xArray(ii), psiSG_RL(ii)
 end do
close(100)



  end subroutine
end module
