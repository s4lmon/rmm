module materialtypes
	use constants
  use diamondsolve
  use backdiamond
  use rmm
  use extrapolate
	implicit none

	 type :: material
    real(kind=dp),dimension(:),allocatable :: todo
    real(kind=dp) :: domainLength,subLength

    !Alphas and betas here
    integer :: numDisc, numSub, ii,total
    !sourcesR is for right hand boundary
      real(kind=dp),dimension(4) :: sourcesL,sourcesR,outputsR,outputsL
      real(kind=dp),dimension(4) :: captureCS
      real(kind=dp),dimension(4,4) :: particleCS,factors
      real(kind=dp) :: sourceFN,sourceTN, sourcePG, sourceSG
      !real(kind=dp),dimension(:),allocatable :: xArray,psiFN,psiTN_LR,psiTN_RL,psiPG,psiSG_LR,psiSG_RL
      real(kind=dp),dimension(:),allocatable :: subX,subPsiFN,subPsiTN_LR,subPsiTN_RL,subPsiPG,subPsiSG_LR,subPsiSG_RL
      real(kind=dp) :: fn_tn,fn_sg,tn_sg,pg_sg
      real(kind=dp),dimension(4,4,2) :: decays
      !real(kind=dp),dimension(:),allocatable :: intervalFN,intervalTN_LR,intervalTN_RL,intervalPG,intervalSG_LR,intervalSG_RL
    !storing subPsiParticle production for each input particle
    !e.g. causeFN is a four dimensional array where dimension 1 represents the FN it causes, dimension 2 is the TN etc
    real(kind=dp),dimension(:,:),allocatable :: causeFN, causeTN, causePG, causeSG
  ! real(kind=dp) :: start,finish,omg
    contains
  	 procedure :: psi_interface=>psi_interface_material
    end type


  contains
    subroutine material_setup(m, countMat,numDisc,subLength)
      type(material), intent(inout) :: m
      integer, intent(in) :: countMat,numDisc
      real(kind=dp),intent(in) :: subLength

      m%subLength = subLength
      m%numDisc = numDisc

      allocate(m%subPsiFN(numDisc))
      allocate(m%subPsiTN_LR(numDisc))
      allocate(m%subPsiTN_RL(numDisc))
      allocate(m%subPsiPG(numDisc))
      allocate(m%subPsiSG_LR(numDisc))
      allocate(m%subPsiSG_RL(numDisc))
      allocate(m%subX(numDisc))
      allocate(m%causeFN(6,numDisc))
      allocate(m%causeTN(6,numDisc))
      allocate(m%causePG(6,numDisc))
      allocate(m%causeSG(6,numDisc))

      if (countMat .eq. 1) then
         m%fn_tn=0.5_dp
          m%fn_sg=0.5_dp
          m%tn_sg=0.5_dp
          m%pg_sg=0.5_dp

!===building for computation
          m%factors=0.0_dp
           m%factors(2,1)=m%fn_tn
           m%factors(4,1)=m%fn_sg
           m%factors(4,2)=m%tn_sg
           m%factors(4,3)=m%pg_sg
          m%sourceFN=1.0_dp
          m%sourceTN=1.0_dp
          m%sourcePG=1.0_dp
          m%sourceSG=0.0_dp
           m%sourcesL(1)=m%sourceFN
            m%sourcesL(2)=m%sourceTN
            m%sourcesL(3)=m%sourcePG
            m%sourcesL(4)=m%sourceSG
          m%sourcesR(1)=0.0_dp
          m%sourcesR(2)=0.0_dp
         m%sourcesR(3)=0.0_dp
          m%sourcesR(4)=0.0_dp
                  m%particleCS=0.0_dp
          m%particleCS(1,2)=1.0_dp
          m%particleCS(1,4)=1.0_dp
          m%particleCS(2,4)=1.0_dp
          m%particleCS(3,4)=1.0_dp
          m%captureCS=(/1.0_dp,1.0_dp,1.0_dp,1.0_dp/)
      ! else if (countMat .eq. 2) then
      ! ...
      end if 
      end subroutine
    subroutine psi_interface_material(self)!,factors,particleCS,captureCS,decays)!,CS)
      class(material),intent(inout) :: self



            ! real (kind = dp), intent(in) :: factors
      ! real (kind = dp), intent(in) :: particleCS
      ! real (kind = dp), intent(in) :: captureCS
      ! real (kind = dp), intent(inout) :: decays
     ! type(material) :: material
!       real(kind=dp),dimension(:),allocatable :: todo
!       real(kind=dp) :: domainLength,subLength
!       integer :: numDisc, numSub, ii,total
!       !sourcesR is for right hand boundary
!       real(kind=dp),dimension(4) :: sourcesL, sourcesR, outputsR, outputsL
!       real(kind=dp),dimension(4) :: captureCS
!       ! real(kind=dp),dimension(4,4) :: particleCS,factors
!       real(kind=dp) :: sourceFN,sourceTN, sourcePG, sourceSG
!       real(kind=dp),dimension(:),allocatable :: xArray,psiFN,psiTN_LR,psiTN_RL,psiPG,psiSG_LR,psiSG_RL
!       real(kind=dp),dimension(:),allocatable :: subX,subPsiFN,subPsiTN_LR,subPsiTN_RL,subPsiPG,subPsiSG_LR,subPsiSG_RL
!       real(kind=dp) :: fn_tn,fn_sg,tn_sg,pg_sg
!       !real(kind=dp),dimension(4,4,2) :: decays
!       real(kind=dp),dimension(:),allocatable :: intervalFN,intervalTN_LR,intervalTN_RL,intervalPG,intervalSG_LR,intervalSG_RL
!       !storing subPsiParticle production for each input particle
!       !e.g. causeFN is a four dimensional array where dimension 1 represents the FN it causes, dimension 2 is the TN etc
!       real(kind=dp),dimension(:,:),allocatable :: causeFN, causeTN, causePG, causeSG
!     ! real(kind=dp) :: start,finish,omg 
!           !print*,subX


      !building alpha and beta matrix
      self%sourcesL=0.0_dp
      self%sourcesR=0.0_dp
      self%sourcesL(1)=1.0_dp
!=====================================================================================
!building decay constants (alpha,beta)
!switching on and off different particles and seeing what scattering they cause
    self%decays=0.0_dp
  
    !print*,'ttt',self%particleCS,'end'
   call particleSolve(self%sourcesL,self%sourcesR,self%numDisc,self%captureCS,self%particleCS,self%factors,self%subLength,self%subX,self%subPsiFN,self%subPsiTN_LR,self%subPsiTN_RL,self%subPsiPG,self%subPsiSG_LR,self%subPsiSG_RL,self%outputsR,self%outputsL)
   
    self%decays(1,1,1)=self%outputsR(1)
    self%decays(1,2,1)=self%outputsR(2)
    self%decays(1,4,1)=self%outputsR(4)

    self%decays(1,2,2)=self%outputsL(2)
    self%decays(1,4,2)=self%outputsL(4)
    !prin*t,'xxx',self%decays(1,1,1)

    self%causeFN(1,:)=self%subPsiFN
    self%causeFN(2,:)=self%subPsiTN_LR
    self%causeFN(3,:)=self%subPsiPG
    self%causeFN(4,:)=self%subPsiSG_LR
    self%causeFN(5,:)=self%subPsiTN_RL
    self%causeFN(6,:)=self%subPsiSG_RL
    self%sourcesL=0.0_dp
    self%sourcesL(2)=1.0_dp
   call particleSolve(self%sourcesL,self%sourcesR,self%numDisc,self%captureCS,self%particleCS,self%factors,self%subLength,self%subX,self%subPsiFN,self%subPsiTN_LR,self%subPsiTN_RL,self%subPsiPG,self%subPsiSG_LR,self%subPsiSG_RL,self%outputsR,self%outputsL)
    self%decays(2,2,1)=self%outputsR(2)
    self%decays(2,4,1)=self%outputsR(4)

    self%decays(2,4,2)=self%outputsL(4)

    self%causeTN(1,:)=self%subPsiFN
    self%causeTN(2,:)=self%subPsiTN_LR
    self%causeTN(3,:)=self%subPsiPG
    self%causeTN(4,:)=self%subPsiSG_LR
    self%causeTN(5,:)=self%subPsiTN_RL
    self%causeTN(6,:)=self%subPsiSG_RL

    self%sourcesL=0.0_dp
    self%sourcesL(3)=1.0_dp
   call particleSolve(self%sourcesL,self%sourcesR,self%numDisc,self%captureCS,self%particleCS,self%factors,self%subLength,self%subX,self%subPsiFN,self%subPsiTN_LR,self%subPsiTN_RL,self%subPsiPG,self%subPsiSG_LR,self%subPsiSG_RL,self%outputsR,self%outputsL)
    self%decays(3,3,1)=self%outputsR(3)
    self%decays(3,4,1)=self%outputsR(4)

    self%decays(3,4,2)=self%outputsL(4)

    self%causePG(1,:)=self%subPsiFN
    self%causePG(2,:)=self%subPsiTN_LR
    self%causePG(3,:)=self%subPsiPG
    self%causePG(4,:)=self%subPsiSG_LR
    self%causePG(5,:)=self%subPsiTN_RL
    self%causePG(6,:)=self%subPsiSG_RL

    self%sourcesL=0.0_dp
    self%sourcesL(4)=1.0_dp
   call particleSolve(self%sourcesL,self%sourcesR,self%numDisc,self%captureCS,self%particleCS,self%factors,self%subLength,self%subX,self%subPsiFN,self%subPsiTN_LR,self%subPsiTN_RL,self%subPsiPG,self%subPsiSG_LR,self%subPsiSG_RL,self%outputsR,self%outputsL)
    self%decays(4,4,1)=self%outputsR(4)


    self%causeSG(1,:)=self%subPsiFN
    self%causeSG(2,:)=self%subPsiTN_LR
    self%causeSG(3,:)=self%subPsiPG
    self%causeSG(4,:)=self%subPsiSG_LR
    self%causeSG(5,:)=self%subPsiTN_RL
    self%causeSG(6,:)=self%subPsiSG_RL
    !==========================================
    self%sourcesL(1)=self%sourceFN
    self%sourcesL(2)=self%sourceTN
    self%sourcesL(3)=self%sourcePG
    self%sourcesL(4)=self%sourceSG
   call particleSolve(self%sourcesL,self%sourcesR,self%numDisc,self%captureCS,self%particleCS,self%factors,self%subLength,self%subX,self%subPsiFN,self%subPsiTN_LR,self%subPsiTN_RL,self%subPsiPG,self%subPsiSG_LR,self%subPsiSG_RL,self%outputsR,self%outputsL)

    end subroutine
  end module