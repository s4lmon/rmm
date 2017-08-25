module extrapolate
  use rmm
  use materialtypes
  implicit none



contains
  subroutine extra_bound(domainLength,numDisc,numSub,sources,matArray,intervalFN,intervalTN_LR,intervalTN_RL,intervalPG,intervalSG_LR,intervalSG_RL,indexMat)
    type(material),dimension(:),intent(in) :: matArray
    integer,dimension(:) :: indexMat 

    integer :: numSub,numDisc,total,ii,temp
    real(kind=dp) :: domainLength,writeLoop
    real(kind=dp),dimension(4,4,2) :: decays
    real(kind=dp),dimension(:) :: sources
    real(kind=dp),dimension(4) :: multiplier !multiplier used to represent decay at each subinterval

    !interval arrays for calculating values on boundaries
    real(kind=dp),dimension(:) :: intervalFN,intervalTN_LR,intervalTN_RL,intervalPG,intervalSG_LR,intervalSG_RL

    ! Referencing stuff in material array:
    ! bla bla = m_array(ii)%decays 

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

         temp=indexMat((ii-1)*numDisc+1)
         
         !print*,'temp',temp
         intervalFN(ii)=intervalFN(ii-1)*matArray(temp)%decays(1,1,1)
         !print*,intervalFN(ii)
         intervalPG(ii)=intervalPG(ii-1)*matArray(temp)%decays(3,3,1)
         intervalTN_LR(ii)=intervalTN_LR(ii-1)*matArray(temp)%decays(2,2,1)+intervalFN(ii-1)*matArray(temp)%decays(1,2,1)
         
       end do

      intervalTN_RL(numSub+1)=0.0_dp
      intervalSG_RL(numSub+1)=0.0_dp

      do ii=numSub,2,-1
       temp=indexMat((ii-1)*numDisc+1)
         intervalTN_RL(ii)=matArray(temp)%decays(2,2,1)*intervalTN_RL(ii+1)+matArray(temp)%decays(1,2,2)*intervalFN(ii)
         intervalSG_RL(ii)=matArray(temp)%decays(1,4,2)*intervalFN(ii)+matArray(temp)%decays(2,4,2)*intervalTN_LR(ii)+matArray(temp)%decays(3,4,2)*intervalPG(ii)+matArray(temp)%decays(4,4,1)*intervalSG_RL(ii+1)
         !print*,temp,'mat',matArray(temp)%decays(2,2,1),intervalSG_RL(ii),intervalSG_RL(ii+1)-intervalSG_LR(ii)
         !intervalTN_RL(ii)=
!         print*,intervalTN_LR(ii),intervalTN_RL(ii)
!         print*,decays(1,2,1)
      end do

      do ii=2,numSub
        temp=indexMat((ii-1)*numDisc+1)
        intervalSG_LR(ii)=intervalSG_LR(ii-1)*matArray(temp)%decays(4,4,1)+intervalFN(ii-1)*matArray(temp)%decays(1,4,1)+intervalTN_LR(ii-1)*matArray(temp)%decays(2,4,1)+intervalPG(ii-1)*matArray(temp)%decays(3,4,1)+intervalTN_RL(ii)*matArray(temp)%decays(2,4,2)
       ! print*,temp,'mat',matArray(temp)%decays(4,4,1),intervalSG_LR(ii),intervalSG_LR(ii)-intervalSG_LR(ii-1)
      end do

  !=============================uncomment
  !     print*,'look',decays(2,4,2),decays(2,4,1)
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

  subroutine extra_psi(domainLength,numDisc,numSub,xArray,intervalFN,intervalTN_LR,intervalTN_RL,intervalPG,intervalSG_LR,intervalSG_RL,matArray,psiFN,psiTN_LR,psiTN_RL,psiPG,psiSG_LR,psiSG_RL,numMat,indexMat)
    type(material),dimension(:),intent(inout) :: matArray
    integer,dimension(:) :: indexMat 
    integer :: hold

    real(kind=dp) :: domainLength

    integer :: numSub,numDisc, numMat,total,ii,jj,psiIndex

    real(kind=dp),dimension(:) :: xArray
    real(kind=dp),dimension(:),allocatable :: psiFN,psiTN_LR,psiTN_RL,psiPG,psiSG_LR,psiSG_RL
!    real(kind=dp),dimension(:) :: subPsiFN,subPsiTN,subPsiPG,subPsiSG
    real(kind=dp),dimension(:) :: intervalFN,intervalTN_LR,intervalTN_RL,intervalPG,intervalSG_LR,intervalSG_RL
    total = numSub*numDisc
  
    allocate(psiFN(total))
    allocate(psiTN_LR(total))
    allocate(psiTN_RL(total))
    allocate(psiPG(total))
    allocate(psiSG_LR(total))
    allocate(psiSG_RL(total))

    !print*,matArray(1)%tempFN(1,1)

    do jj=1,numMat
      do ii=1,numDisc
        matArray(jj)%tempFN(5,((numDisc+1)-ii))=matArray(jj)%causeFN(5,ii)
        matArray(jj)%tempTN(2,((numDisc+1)-ii))=matArray(jj)%causeTN(2,ii)
        matArray(jj)%tempFN(6,((numDisc+1)-ii))=matArray(jj)%causeFN(6,ii)
        !TN causing SG going right to left
        matArray(jj)%tempTN(6,((numDisc+1)-ii))=matArray(jj)%causeTN(6,ii)
        matArray(jj)%tempPG(6,((numDisc+1)-ii))=matArray(jj)%causePG(6,ii)
        matArray(jj)%tempSG(4,((numDisc+1)-ii))=matArray(jj)%causeSG(4,ii)
        !TN causing SG going left to right
        matArray(jj)%tempTN(4,((numDisc+1)-ii))=matArray(jj)%causeTN(4,ii)
      end do
    end do 
    !print*,matArray(1)%tempTN(4,:)

    do ii=1,numSub
      hold=indexMat((ii-1)*numDisc+1)
      psiFN((ii-1)*numDisc+1:ii*numDisc)=matArray(hold)%causeFN(1,:)*intervalFN(ii)
      psiPG((ii-1)*numDisc+1:ii*numDisc)=matArray(hold)%causePG(3,:)*intervalPG(ii)
      psiSG_LR((ii-1)*numDisc+1:ii*numDisc)=matArray(hold)%causeSG(4,:)*intervalSG_LR(ii)+matArray(hold)%causeFN(4,:)*intervalFN(ii)+matArray(hold)%causeTN(4,:)*intervalTN_LR(ii)+matArray(hold)%causePG(4,:)*intervalPG(ii)+matArray(hold)%causeTN(6,:)*intervalTN_RL(ii+1)

      psiTN_RL((ii-1)*numDisc+1:ii*numDisc)=matArray(hold)%tempTN(2,:)*intervalTN_RL(ii+1)+matArray(hold)%causeFN(5,:)*intervalFN(ii)
      psiSG_RL((ii-1)*numDisc+1:ii*numDisc)=matArray(hold)%causeFN(6,:)*intervalFN(ii)+matArray(hold)%causeTN(6,:)*intervalTN_LR(ii)+matArray(hold)%tempTN(4,:)*intervalTN_RL(ii+1)+matArray(hold)%tempSG(4,:)*intervalSG_RL(ii+1)+matArray(hold)%causePG(6,:)*intervalPG(ii)

      psiTN_LR((ii-1)*numDisc+1:ii*numDisc)=matArray(hold)%causeTN(2,:)*intervalTN_LR(ii)+matArray(hold)%causeFN(2,:)*intervalFN(ii)
     
      print*,hold,psiFN(ii*numDisc),psiFN(ii*numDisc)-psiFN(((ii-1)*numDisc)+1),ii
      
      !todo
      print*,xArray(ii*numDisc)
    end do

    psiIndex = 1
    do ii=1,numSub
      do jj=1,numDisc
        hold=indexMat(jj*ii)
        psiFN(psiIndex)=matArray(hold)%causeFN(1,jj)*intervalFN(ii)
        psiIndex = psiIndex+1
      end do
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
