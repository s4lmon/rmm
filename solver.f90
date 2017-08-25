program solver
  use rmm
  use extrapolate
  use materialtypes
  implicit none
    real(kind=dp) :: domainLength,subLength
    integer :: numDisc, numSub, ii,total
    integer :: numMat
    real(kind=dp),dimension(:),allocatable :: lengthMat
    integer,dimension(:),allocatable :: indexMat 
    !sourcesR is for right hand boundary
      real(kind=dp),dimension(4) :: sourcesL,sourcesR,outputsR,outputsL
      real(kind=dp),dimension(4) :: captureCS
      real(kind=dp),dimension(4,4) :: particleCS,factors
      real(kind=dp) :: sourceFN,sourceTN, sourcePG, sourceSG
      real(kind=dp),dimension(:),allocatable :: xArray,psiFN,psiTN_LR,psiTN_RL,psiPG,psiSG_LR,psiSG_RL
      real(kind=dp),dimension(:),allocatable :: subX,subPsiFN,subPsiTN_LR,subPsiTN_RL,subPsiPG,subPsiSG_LR,subPsiSG_RL
      ! real(kind=dp) :: fn_tn,fn_sg,tn_sg,pg_sg
      ! real(kind=dp),dimension(4,4,2) :: decays
      real(kind=dp),dimension(:),allocatable :: intervalFN,intervalTN_LR,intervalTN_RL,intervalPG,intervalSG_LR,intervalSG_RL
      type(material),dimension(:),allocatable :: material_properties
    !storing subPsiParticle production for each input particle
    !e.g. causeFN is a four dimensional array where dimension 1 represents the FN it causes, dimension 2 is the TN etc
    real(kind=dp),dimension(:,:),allocatable :: causeFN, causeTN, causePG, causeSG
! real(kind=dp) :: start,finish,omg
! numDisc=1
! do while(numDisc<=100000)
! call cpu_time(start)
!INPUT VALUES
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
 
    numDisc=50
    numSub=30
    ! numDisc=1
    ! numSub=100000/numDisc
    !print*,numDisc,numSub
    total=numSub*numDisc

    numMat=3
    allocate(lengthMat(numMat))
    allocate(material_properties(numMat))
    lengthMat(1)=0.5_dp
    lengthMat(2)=0.5_dp
    lengthMat(3)=0.5_dp
    domainLength=sum(lengthMat)
    !print*,lengthMat(1)/domainLength

    if (mod(real(total),domainLength)==0.0_dp) then
    allocate(indexMat((total)))
      do ii=1,numMat
        indexMat((ii-1)*int((lengthMat(ii)/domainLength)*total)+1:(ii)*int((lengthMat(ii)/domainLength)*total))=ii
      end do 
      else 
      print*,"mismatch"
      stop
    end if 
      !print*,indexMat


subLength=domainLength/numSub



    allocate(intervalFN(numSub))
    allocate(intervalTN_LR(numSub))
    allocate(intervalTN_RL(numSub+1))
    allocate(intervalPG(numSub))
    allocate(intervalSG_LR(numSub))
    allocate(intervalSG_RL(numSub+1))
    allocate(xArray(total))


    
    !call plot()

    !!!!!!!!!!!!!!!!!!
    ! Material Setup !
    !!!!!!!!!!!!!!!!!!

    do ii=1,numMat
      call material_setup(material_properties(ii),ii,numDisc,subLength)
      call material_properties(ii)%psi_interface
    end do
    !print*,'test',material_properties(1)%factors(2,1),material_properties(2)%factors(2,1)
    !print*,'test',material_properties%decays

    !print*,material_properties(1)%decays
    ! do ii=1,size(material_properties(1)%subX)
    !   print*,material_properties(1)%subX(ii)
    ! end do 
    !pprint*,material_properties(1)%subX
    do ii=1,numSub

           xArray((ii-1)*numDisc+1:ii*numDisc)=material_properties(1)%subX+subLength*(ii-1)
    end do
    !print*,xArray
!(material_properties,particleCS,captureCS,decays)
! print*, "Decays test: ",
! do ii=1,
!print*,material_properties(1)%decays
!     !todo uncomment till end


!     do ii=1,numDisc
!           !print*,'hello',causeFN(:,ii)
!     end do
call extra_bound(domainLength,numDisc,numSub,sourcesL,material_properties,intervalFN,intervalTN_LR,intervalTN_RL,intervalPG,intervalSG_LR,intervalSG_RL,indexMat)

! call extra_bound(material_properties, material_index_array, [rest...])
 
print*,material_properties(1)%causeFN(5,2)
ii=1
!print*,indexMat((ii-1)*numDisc+1)
!print*,indexMat((ii-1)*numDisc+1)
!rint*,material_properties(1)%tempFN(1,5)
!rint*,material_properties(1)%tempFN(1,1)
call extra_psi(domainLength,numDisc,numSub,xArray,intervalFN,intervalTN_LR,intervalTN_RL,intervalPG,intervalSG_LR,intervalSG_RL,material_properties,psiFN,psiTN_LR,psiTN_RL,psiPG,psiSG_LR,psiSG_RL,numMat,indexMat)
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



















   !   call cpu_time(finish)
   !       print*,"Time =", finish-start, "seconds."
   ! open(100, file="runtime.dat", status="unknown")
   !     write(100,*) numDisc,finish-start
   ! close(100)
   ! numDisc=numDisc*10
  !  deallocate(subPsiFN)
  !  deallocate(subPsiTN_LR)
  !  deallocate(subPsiTN_RL)
  !  deallocate(subPsiPG)
  !  deallocate(subPsiSG_LR)
  !  deallocate(subPsiSG_RL)
  !  deallocate(subX)
 
  !  deallocate(causeFN)
  !  deallocate(causeTN)
  !  deallocate(causePG)
  !  deallocate(causeSG)
  ! !  deallocate(tempFN)
  ! !  deallocate(tempTN)
  ! !  deallocate(tempPG)
  ! !  deallocate(tempSG)
  !  deallocate(intervalFN)
  !  deallocate(intervalTN_LR)
  !  deallocate(intervalTN_RL)
  !  deallocate(intervalPG)
  !  deallocate(intervalSG_LR)
  !  deallocate(intervalSG_RL)
  !  deallocate(xArray)
  !  deallocate(psiFN)
  !  deallocate(psiTN_LR)
  !  deallocate(psiTN_RL)
  !  deallocate(psiPG)
  !  deallocate(psiSG_LR)
  !  deallocate(psiSG_RL)
! end do
end program
