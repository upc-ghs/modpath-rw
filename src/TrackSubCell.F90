module TrackSubCellModule
  use ParticleLocationModule,only : ParticleLocationType
  use TrackSubCellResultModule,only : TrackSubCellResultType
  use ModpathSubCellDataModule,only : ModpathSubCellDataType
  use ParticleTrackingOptionsModule,only : ParticleTrackingOptionsType
  use ModpathCellDataModule,only : ModpathCellDataType
  use rng_par_zig, only : rng_uni, rng_init ! RGN parallel ifort
#ifdef _OPENMP
  use omp_lib
#endif
  implicit none
  
! Set default access status to private
  private

  type,public :: TrackSubCellType
    type(ModpathSubCellDataType) :: SubCellData

    ! RWPT
    ! qx, qy, qz, qnorm
    doubleprecision, dimension(4) :: qCorner000
    doubleprecision, dimension(4) :: qCorner100
    doubleprecision, dimension(4) :: qCorner010
    doubleprecision, dimension(4) :: qCorner110
    doubleprecision, dimension(4) :: qCorner001
    doubleprecision, dimension(4) :: qCorner101
    doubleprecision, dimension(4) :: qCorner011
    doubleprecision, dimension(4) :: qCorner111

    ! corner porosities
    doubleprecision :: porosity000
    doubleprecision :: porosity100
    doubleprecision :: porosity010
    doubleprecision :: porosity110
    doubleprecision :: porosity001
    doubleprecision :: porosity101
    doubleprecision :: porosity011
    doubleprecision :: porosity111

    ! These arrays will contain the dispersion terms 
    ! at the corners of each cell employed to calculate
    ! the divergence of dispersion The naming convention
    ! follows the rule: 
    !
    !   - Dxyx: evaluate Dxy terms for derivative x
    !   - Dyzz: evaluate Dzy terms for derivative z
    !   - and so on ...
    !
    ! (1): 000
    ! (2): 100
    ! (3): 010
    ! (4): 110
    ! (5): 001
    ! (6): 101
    ! (7): 011
    ! (8): 111
    doubleprecision, dimension(8) :: Dxxx
    doubleprecision, dimension(8) :: Dyyy
    doubleprecision, dimension(8) :: Dzzz
    doubleprecision, dimension(8) :: Dxyx
    doubleprecision, dimension(8) :: Dxzx
    doubleprecision, dimension(8) :: Dxyy
    doubleprecision, dimension(8) :: Dyzy
    doubleprecision, dimension(8) :: Dxzz
    doubleprecision, dimension(8) :: Dyzz
                                    
    ! Interpolation indexes !
    ! Corner discharge components indexes
    ! 8 corners, 3 subcell indexes, 1 face id
    integer, dimension(8,4) :: cornerXComponentIndexes, &
                               cornerYComponentIndexes, & 
                               cornerZComponentIndexes
    ! Corner porosity sub cell indexes
    ! These indexes also work for other quantities
    ! that could be averaged weighted by volume.
    integer, dimension(8,6)  :: cornerPorosityIndexes ! 8 corners, 6 neighbor subcells

    ! RWPT Pointers
    procedure(Advection)            , pass, pointer :: AdvectionDisplacement=>null()
    procedure(ExitFaceAndTimeStep)  , pass, pointer :: ExitFaceAndUpdateTimeStep=>null()
    procedure(DispersionModel)      , pass, pointer :: ComputeRWPTDisplacements=>null()
    procedure(RandomDisplacement)   , pass, pointer :: DisplacementRandomDischarge=>null()
    procedure(RandomDisplacementAxi), pass, pointer :: DisplacementRandomDischargeAxisymmetric=>null()
    procedure(RandomGenerator)      , pass, pointer :: GenerateStandardNormalRandom=>null()
    procedure(CornerPorosity)       , pass, pointer :: ComputeCornerPorosity=>null()
    procedure(CornerDispersion)     , pass, pointer :: ComputeCornerDispersion=>null()

    ! Needed for OBS (?)
    type(TrackSubCellResultType) :: TrackSubCellResult

    ! Point towards trackingOptions
    integer, dimension(:), pointer :: dimensions
    integer, pointer :: nDim

  contains
    procedure,private :: CalculateDT=>pr_CalculateDT
    procedure,private :: NewXYZ=>pr_NewXYZ
    procedure :: ExecuteTracking=>pr_ExecuteTracking

    ! RWPT
    procedure :: InitializeRandomWalk=>pr_InitializeRandomWalk
    procedure :: ExecuteRandomWalkParticleTracking=>pr_ExecuteRandomWalkParticleTracking
    procedure :: LinearInterpolationVelocities=>pr_LinearInterpolationVelocities
    procedure :: ComputeRandomWalkTimeStep=>pr_ComputeRandomWalkTimeStep
    procedure :: ComputeCornerVariables=>pr_ComputeCornerVariables
    procedure :: ComputeCornerDischarge=>pr_ComputeCornerDischarge
    procedure :: GetInterpolatedCornerDischarge=>pr_GetInterpolatedCornerDischarge
    procedure :: SetCornerComponentsIndexes=>pr_SetCornerComponentsIndexes
    procedure :: GetInterpolatedCornerPorosity=>pr_GetInterpolatedCornerPorosity
    procedure :: SetCornerPorosityIndexes=> pr_SetCornerPorosityIndexes
    procedure :: Trilinear=>pr_Trilinear
    procedure :: TrilinearDerivative=>pr_TrilinearDerivative
    procedure :: TrilinearDerivativeX=>pr_TrilinearDerivativeX
    procedure :: TrilinearDerivativeY=>pr_TrilinearDerivativeY
    procedure :: TrilinearDerivativeZ=>pr_TrilinearDerivativeZ
    procedure :: SetDispersionDisplacement=>pr_SetDispersionDisplacement
    procedure :: DispersionDivergence=>pr_DispersionDivergence
    procedure :: AdvectionDisplacementExponential=>pr_AdvectionDisplacementExponential
    procedure :: AdvectionDisplacementEulerian=>pr_AdvectionDisplacementEulerian
    procedure :: NewtonRaphsonTimeStep=>NewtonRaphsonTimeStepExponentialAdvection

  end type


  ! Interfaces
  abstract interface

    ! Advection model
    subroutine Advection( this, x, y, z, dt, vx, vy, vz, &
                          dAdvx, dAdvy, dAdvz )
      import TrackSubCellType
      ! this
      class(TrackSubCellType) :: this 
      !input
      doubleprecision, intent(in) :: x, y, z, dt
      doubleprecision, intent(in) :: vx, vy, vz
      ! output
      doubleprecision, intent(inout) :: dAdvx, dAdvy, dAdvz
    end subroutine Advection

    ! Update time step model,
    ! determined by advection model
    subroutine ExitFaceAndTimeStep( this, x, y, z, nx, ny, nz, & 
               vx, vy, vz, divDx, divDy, divDz, dBx, dBy, dBz, &
               t, dt, dtxyz, exitFace )
      import TrackSubCellType
      ! this
      class(TrackSubCellType) :: this 
      ! input
      doubleprecision, intent(in) :: x, y, z
      doubleprecision, intent(in) :: nx, ny, nz
      doubleprecision, intent(in) :: vx, vy, vz
      doubleprecision, intent(in) :: divDx, divDy, divDz
      doubleprecision, intent(in) :: dBx, dBy, dBz
      ! output
      doubleprecision, intent(inout) :: t, dt
      integer, intent(inout)         :: exitFace
      doubleprecision, dimension(3), intent(inout) :: dtxyz
    end subroutine ExitFaceAndTimeStep

    ! Dispersion model 
    ! Interfaces between linear and nonlinear dispersion
    subroutine DispersionModel( this, x, y, z, vx, vy, vz, &
                                      dt, trackingOptions, &
                                      dAdvx, dAdvy, dAdvz, &
                                            dBx, dBy, dBz, &
                                      divDx, divDy, divDz  )
      import TrackSubCellType
      import ParticleTrackingOptionsType
      !----------------------------------------------------------------
      ! this
      class(TrackSubCellType) :: this
      ! input
      type(ParticleTrackingOptionsType),intent(in) :: trackingOptions
      doubleprecision, intent(in)    :: x, y, z, dt
      ! output
      doubleprecision, intent(inout) :: vx, vy, vz 
      doubleprecision, intent(inout) :: dAdvx, dAdvy, dAdvz 
      doubleprecision, intent(inout) :: dBx, dBy, dBz 
      doubleprecision, intent(inout) :: divDx, divDy, divDz  
    end subroutine DispersionModel 

    ! Random displacement function
    ! Determined by dimensions 
    subroutine RandomDisplacement( this, x, y, z, alphaL, alphaT, &
                                             dMEff, dBx, dBy, dBz )
      import TrackSubCellType
      !------------------------------------------------------------
      class (TrackSubCellType) :: this
      ! input
      doubleprecision, intent(in)    :: x, y, z
      doubleprecision, intent(in)    :: alphaL, alphaT, dMEff
      ! output
      doubleprecision, intent(out) :: dBx, dBy, dBz
    end subroutine RandomDisplacement

    ! Random displacement for axisymmetric dispersion 
    ! Determined by dimensions 
    subroutine RandomDisplacementAxi( this, x, y, z, alphaL, alphaT, alphaTH, &
                                                         dMEff, dBx, dBy, dBz )
      import TrackSubCellType
      !------------------------------------------------------------
      class (TrackSubCellType) :: this
      ! input
      doubleprecision, intent(in)    :: x, y, z
      doubleprecision, intent(in)    :: alphaL, alphaT, alphaTH, dMEff
      ! output
      doubleprecision, intent(out) :: dBx, dBy, dBz
    end subroutine RandomDisplacementAxi

    ! Random generator function
    subroutine RandomGenerator( this, random_value )
      import TrackSubCellType
      !---------------------------------------------
      class (TrackSubCellType) :: this
      ! input/output
      doubleprecision, intent(out) :: random_value
    end subroutine RandomGenerator

    ! Compute corner porosity function 
    subroutine CornerPorosity( this, neighborSubCellVolume,&
                                   neighborSubCellPorosity )
      import TrackSubCellType
      !-------------------------------------------------------------------
      class (TrackSubCellType) :: this
      ! input
      doubleprecision, dimension(18), intent(in) :: neighborSubCellVolume   
      doubleprecision, dimension(18), intent(in) :: neighborSubCellPorosity
    end subroutine CornerPorosity 

    ! Compute corner dispersion function 
    subroutine CornerDispersion( this, &
                             qprod000, & 
                             qprod100, &
                             qprod010, &
                             qprod110, &
                             qprod001, &
                             qprod101, &
                             qprod011, &
                             qprod111  )
      import TrackSubCellType
      !-------------------------------------------------------------------
      class (TrackSubCellType) :: this
      ! input
      doubleprecision, dimension(8), intent(in) :: qprod000
      doubleprecision, dimension(8), intent(in) :: qprod100
      doubleprecision, dimension(8), intent(in) :: qprod010
      doubleprecision, dimension(8), intent(in) :: qprod110
      doubleprecision, dimension(8), intent(in) :: qprod001
      doubleprecision, dimension(8), intent(in) :: qprod101
      doubleprecision, dimension(8), intent(in) :: qprod011
      doubleprecision, dimension(8), intent(in) :: qprod111
    end subroutine CornerDispersion

  end interface


contains


  subroutine pr_ExecuteTracking(this,stopIfNoExit,initialLocation,maximumTime, trackingResult)
    !-------------------------------------------------------------------
    !
    !-------------------------------------------------------------------
    ! Specifications
    !-------------------------------------------------------------------
    implicit none
    class(TrackSubCellType) :: this
    logical,intent(in) :: stopIfNoExit
    type(ParticleLocationType),intent(in) :: initialLocation
    doubleprecision,intent(in) :: maximumTime
    type(TrackSubCellResultType),intent(inout) :: trackingResult
    integer :: cellNumber
    doubleprecision :: initialX,initialY,initialZ,initialTime
    doubleprecision :: vx1,vx2,vy1,vy2,vz1,vz2
    doubleprecision :: vx,dvxdx,dtx,vy,dvydy,dty,vz,dvzdz,dtz,dt
    doubleprecision :: t,x,y,z
    integer :: exitFace
    integer :: statusVX,statusVY,statusVZ
    !-------------------------------------------------------------------

    call trackingResult%Reset()
    
    cellNumber = initialLocation%CellNumber
    initialX = initialLocation%LocalX
    initialY = initialLocation%LocalY
    initialZ = initialLocation%LocalZ
    initialTime = initialLocation%TrackingTime
    
    trackingResult%CellNumber = cellNumber
    trackingResult%Row = this%SubCellData%Row
    trackingResult%Column = this%SubCellData%Column
    trackingResult%InitialLocation%CellNumber = cellNumber
    trackingResult%InitialLocation%LocalX = initialX
    trackingResult%InitialLocation%LocalY = initialY
    trackingResult%InitialLocation%LocalZ = initialZ
    trackingResult%InitialLocation%TrackingTime = initialTime
    trackingResult%FinalLocation%LocalX = initialX
    trackingResult%FinalLocation%LocalY = initialY
    trackingResult%FinalLocation%LocalZ = initialZ
    trackingResult%FinalLocation%TrackingTime = initialTime
    trackingResult%MaximumTime = maximumTime
    trackingResult%Status = trackingResult%Status_Undefined()
    
    if(stopIfNoExit) then
      if(.not. this%SubCellData%HasExitFace()) then
        trackingResult%Status = trackingResult%Status_NoExitPossible()
        return
      end if
    end if
    
    ! Make local copies of face velocities for convenience
    vx1 = this%SubCellData%VX1
    vx2 = this%SubCellData%VX2
    vy1 = this%SubCellData%VY1
    vy2 = this%SubCellData%VY2
    vz1 = this%SubCellData%VZ1
    vz2 = this%SubCellData%VZ2

    ! Compute time of travel to each possible exit face
    statusVX = this%CalculateDT(vx1, vx2, this%SubCellData%DX, initialX, vx, dvxdx, dtx)
    statusVY = this%CalculateDT(vy1, vy2, this%SubCellData%DY, initialY, vy, dvydy, dty)
    statusVZ = this%CalculateDT(vz1, vz2, this%SubCellData%DZ, initialZ, vz, dvzdz, dtz)

    ! This snippet selects dt and determines the exitFace
    exitFace = 0
    dt = 1.0d+30
    if((statusVX .lt. 2) .or. (statusVY .lt. 2) .or. (statusVZ .lt. 2)) then
      dt = dtx
      if(vx .lt. 0d0) then
        exitFace = 1
      else if(vx .gt. 0) then
        exitFace = 2
      end if
      
      if(dty .lt. dt) then
        dt = dty
        if(vy .lt. 0d0) then
          exitFace = 3
        else if(vy .gt. 0d0) then
          exitFace = 4
        end if
      end if
      
      if(dtz .lt. dt) then
        dt = dtz
        if(vz .lt. 0d0) then
          exitFace = 5
        else if(vz .gt. 0d0) then
          exitFace = 6
        end if
      end if
    else
    end if

    ! Increment t
    t = initialTime + dt
    
    ! If the maximum time is less than the computed exit time, then calculate the
    ! particle location at the maximum time and set the final time equal to
    ! the maximum time. Set status to 2 (MaximumTimeReached) and return.
    if(maximumTime .lt. t) then
      dt = maximumTime - initialTime
      t = maximumTime
      x = this%NewXYZ(vx, dvxdx, vx1, vx2, dt, initialX, this%SubCellData%DX, statusVX)
      y = this%NewXYZ(vy, dvydy, vy1, vy2, dt, initialY, this%SubCellData%DY, statusVY)
      z = this%NewXYZ(vz, dvzdz, vz1, vz2, dt, initialZ, this%SubCellData%DZ, statusVZ)
      trackingResult%ExitFace = 0
      trackingResult%FinalLocation%CellNumber = cellNumber
      trackingResult%FinalLocation%LocalX = x
      trackingResult%FinalLocation%LocalY = y
      trackingResult%FinalLocation%LocalZ = z
      trackingResult%FinalLocation%TrackingTime = t
      trackingResult%Status = trackingResult%Status_ReachedMaximumTime()
    else
        ! Otherwise, if the computed exit time is less than or equal to the maximum time,
        ! then calculate the exit location and set the final time equal to the computed
        ! exit time.
        if((exitFace .eq. 1) .or. (exitFace .eq.2)) then
            x = 0d0
            y = this%NewXYZ(vy, dvydy, vy1, vy2, dt, initialY, this%SubCellData%DY, statusVY)
            z = this%NewXYZ(vz, dvzdz, vz1, vz2, dt, initialZ, this%SubCellData%DZ, statusVZ)
            if(exitFace .eq. 2) x = 1.0d0
        else if((exitFace .eq. 3) .or. (exitFace .eq.4)) then
            x = this%NewXYZ(vx, dvxdx, vx1, vx2, dt, initialX, this%SubCellData%DX, statusVX)
            y = 0d0
            z = this%NewXYZ(vz, dvzdz, vz1, vz2, dt, initialZ, this%SubCellData%DZ, statusVZ)
            if(exitFace .eq. 4) y = 1.0d0
        else if((exitFace .eq. 5) .or. (exitFace .eq.6)) then
            x = this%NewXYZ(vx, dvxdx, vx1, vx2, dt, initialX, this%SubCellData%DX, statusVX)
            y = this%NewXYZ(vy, dvydy, vy1, vy2, dt, initialY, this%SubCellData%DY, statusVY)
            z = 0d0
            if(exitFace .eq. 6) z = 1.0d0
        else
            ! If it gets this far, something went wrong. Signal an error condition by
            ! setting Status = 0 (undefined) and then return.
            trackingResult%ExitFace = exitFace
            trackingResult%Status = trackingResult%Status_Undefined()
            return
        end if
        
        ! Assign tracking result data
        trackingResult%ExitFaceConnection = this%SubCellData%Connection(exitFace)
        if(this%SubCellData%Connection(exitFace) .lt. 0) then
            trackingResult%Status = trackingResult%Status_ExitAtInternalFace()
        else
            trackingResult%Status = trackingResult%Status_ExitAtCellFace()
        end if
        trackingResult%ExitFace = exitFace
        trackingResult%FinalLocation%CellNumber = cellNumber
        trackingResult%FinalLocation%LocalX = x
        trackingResult%FinalLocation%LocalY = y
        trackingResult%FinalLocation%LocalZ = z
        trackingResult%FinalLocation%TrackingTime = t
    end if

  end subroutine pr_ExecuteTracking
  
!-------------------------------------------------------------------
  function pr_CalculateDT(this,v1,v2,dx,xL,v,dvdx,dt) result(status)
  implicit none
  class(TrackSubCellType) :: this
  doubleprecision,intent(in) :: v1,v2,dx,xL
  doubleprecision,intent(inout) :: v,dvdx,dt
  doubleprecision :: v2a,v1a,dv,dva,vv,vvv,zro,zrom,x,tol
  doubleprecision :: vr1,vr2,vr,v1v2
  integer :: status
  logical :: noOutflow
  
  ! Initialize variables
  status = -1
  dt = 1.0d+20
  v2a = v2
  if(v2a .lt. 0d0) v2a = -v2a
  v1a = v1
  if(v1a .lt. 0d0) v1a = -v1a
  dv = v2 - v1
  dva = dv
  if(dva .lt. 0d0) dva = -dva
  
  ! Check for a uniform zero velocity in this direction.
  ! If so, set status = 2 and return (dt = 1.0d+20).
  tol = 1.0d-15
  if((v2a .lt. tol) .and. (v1a .lt. tol)) then
    v = 0d0
    dvdx = 0d0
    status = 2
    return
  end if
  
  ! Check for uniform non-zero velocity in this direction. 
  ! If so, set compute dt using the constant velocity, 
  ! set status = 1 and return.
  vv = v1a
  if(v2a .gt. vv) vv = v2a
  vvv = dva / vv
  if(vvv .lt. 1.0d-4) then
    zro = tol
    zrom = -zro
    v = v1
    x = xL * dx
    if(v1 .gt. zro) dt = (dx - x) / v1
    if(v1 .lt. zrom) dt = -x / v1
    dvdx = 0d0
    status = 1
    return
  end if
  
  ! Velocity has a linear variation.
  ! Compute velocity corresponding to particle position
  dvdx = dv / dx
  v = (1.0d0 - xL)*v1 + xL*v2
  
  ! If flow is into the cell from both sides there is no outflow.
  ! In that case, set status = 3 and return
  noOutflow = .true.
  if(v1 .lt. 0d0) noOutflow = .false.
  if(v2 .gt. 0d0) noOutflow = .false.
  if(noOutflow) then
    status = 3
    return
  end if
  
  ! If there is a divide in the cell for this flow direction, check to see if the
  ! particle is located exactly on the divide. If it is, move it very slightly to
  ! get it off the divide. This avoids possible numerical problems related to 
  ! stagnation points.
  if((v1 .le. 0d0) .and. (v2 .ge. 0d0)) then
    if(abs(v)  .le. 0d0) then
      v = 1.0d-20
      if(v2 .le. 0d0) v = -v
    end if
  end if

  ! If there is a flow divide, this check finds out what side of the divide the particle
  ! is on and sets the value of vr appropriately to reflect that location.
  vr1 = v1 / v
  vr2 = v2 / v
  vr = vr1
  if(vr .le. 0d0) vr = vr2
  
  ! Check to see if the velocity is in the same direction throughout the cell (i.e. no flow divide).
  ! Check to see if the product v1*v2 > 0 then the velocity is in the same direction throughout 
  ! the cell (i.e. no flow divide). If so, set the value of vr to reflect the appropriate direction.
  v1v2 = v1*v2
  if(v1v2 .gt. 0d0) then
    if(v .gt. 0d0) vr = vr2
    if(v .lt. 0d0) vr = vr1
  end if
  
  ! Compute travel time to exit face. Return with status = 0
  dt = log(vr) / dvdx
  status = 0
    
  end function pr_CalculateDT
  
!-------------------------------------------------------------------
  function pr_NewXYZ(this,v,dvdx,v1,v2,dt,x,dx,velocityProfileStatus) result(newX)
  implicit none
  class(TrackSubCellType) :: this
  integer,intent(in) :: velocityProfileStatus
  doubleprecision,intent(in) :: v,dvdx,v1,v2,dt,x,dx
  doubleprecision :: newX

  newX = x
  select case (velocityProfileStatus)
    case (1)
      newX = newX + (v1*dt/dx)
    case default
      if(v .ne. 0d0) then
        newX = newX + (v*(exp(dvdx*dt) - 1.0d0)/dvdx/dx)
      end if
  end select
  if(newX .lt. 0d0) newX = 0d0
  if(newX .gt. 1.0d0) newX = 1.0d0
  
  end function pr_NewXYZ


! RWPT
!-------------------------------------------------------------------
  subroutine pr_InitializeRandomWalk(this, trackingOptions)
    use iso_fortran_env, only: int64
    !------------------------------------------------------------
    ! Called in trackingEngine%Initialize
    !------------------------------------------------------------
    ! Specifications
    !------------------------------------------------------------
    implicit none
    class(TrackSubCellType) :: this
    type(ParticleTrackingOptionsType),intent(in), target :: trackingOptions
    ! RNG parallel ifort
    integer(int64), dimension(2) :: seedzigrng
    integer(int64), parameter :: baseseed(2) = [int(Z'1DADBEEFBAADD0D0', int64), &
                                                int(Z'5BADD0D0DEADBEEF', int64)]
    integer :: ompNumThreads
    !------------------------------------------------------------

    ! Assign displacement pointers
    if ( trackingOptions%advectionKind .eq. 1 ) then 
      this%AdvectionDisplacement=>pr_AdvectionDisplacementExponential
      this%ExitFaceAndUpdateTimeStep=>pr_DetectExitFaceAndUpdateTimeStepNewton
    else if ( trackingOptions%advectionKind .eq. 2 ) then 
      this%AdvectionDisplacement=>pr_AdvectionDisplacementEulerian
      this%ExitFaceAndUpdateTimeStep=>pr_DetectExitFaceAndUpdateTimeStepQuadratic
    end if

    ! Set indexes for corner x,y,z components
    call this%SetCornerComponentsIndexes() 

    ! interface for porosities    
    if ( trackingOptions%isUniformPorosity ) then 
      this%ComputeCornerPorosity => pr_ComputeCornerPorosityUniform
    else
      ! Set indexes for corner porosities
      call this%SetCornerPorosityIndexes()
      this%ComputeCornerPorosity => pr_ComputeCornerPorosityInterpolated 
    end if 

    ! Pointers to dimensionality
    this%dimensions => trackingOptions%dimensions
    this%nDim => trackingOptions%nDim
    
    ! Assign random displacement function
    select case(this%nDim) 
      case(1)
        this%DisplacementRandomDischarge => pr_DisplacementRandomDischarge1D
      case(2)
        this%DisplacementRandomDischarge => pr_DisplacementRandomDischarge2D
      case(3)
        this%DisplacementRandomDischarge => pr_DisplacementRandomDischarge
    end select

    select case(this%nDim)
      ! For 1D or 2D problems these are the same than the isotropic model, 
      ! only for interface purposes. 
      case(1)
        this%DisplacementRandomDischargeAxisymmetric => pr_DisplacementRandomDischargeAxisymmetric1D
      case(2)
        this%DisplacementRandomDischargeAxisymmetric => pr_DisplacementRandomDischargeAxisymmetric2D
      case(3)
        this%DisplacementRandomDischargeAxisymmetric => pr_DisplacementRandomDischargeAxisymmetric
    end select

    ! Dispersion displacement function is set in particletrackingengine
    !call this%SetDispersionDisplacement( trackingOptions%dispersionModel )

    ! Assign random generator
    select case(trackingOptions%randomGenFunction)
    case(1)
      this%GenerateStandardNormalRandom => pr_GenerateStandardNormalRandom
    case(2)
#ifdef _OPENMP
      ompNumThreads = omp_get_max_threads()
#else
      ompNumThreads = 1
#endif
      seedzigrng    = baseseed
      call rng_init(ompNumThreads, seedzigrng) 
      this%GenerateStandardNormalRandom => pr_GenerateStandardNormalRandomZig
    case default
      write(*,*)'Selected random generator function not implemented. Stop.'
      stop
    end select


    ! Done
    return


  end subroutine pr_InitializeRandomWalk


  subroutine pr_SetDispersionDisplacement(this, dispersionModel)
  !------------------------------------------------------------
  !
  !------------------------------------------------------------
  ! Specifications
  !------------------------------------------------------------
  implicit none
  class(TrackSubCellType) :: this
  integer :: dispersionModel 
  !------------------------------------------------------------

    ! Assign displacement pointers
    select case(dispersionModel)
    case(0)
      !  Linear isotropic
      this%ComputeRWPTDisplacements => pr_RWPTDisplacementsLinear
      this%ComputeCornerDispersion => pr_CornerDispersionLinearIsotropic
    case(1)
      ! Linear axisymmetric 
      this%ComputeRWPTDisplacements => pr_RWPTDisplacementsAxisymmetric
      this%ComputeCornerDispersion => pr_CornerDispersionAxisymmetric
    !case(2)else if ( dispersionModel .eq.2 ) then
    ! ! Non linear 
    ! this%ComputeRWPTDisplacements => pr_RWPTDisplacementsNonlinear
    case default
      ! Not set !
      ! Some kind of error handling
      write(*,*)&
        'TrackSubCell:SetDispersionDisplacement: dispersionModel ' , dispersionModel ,' NOT implemented. Stop.'
      stop
    end select

    ! Done
    return


  end subroutine pr_SetDispersionDisplacement


  subroutine pr_RWPTDisplacementsLinear(this, x, y, z, vx, vy, vz, &
                                              dt, trackingOptions, &
                                              dAdvx, dAdvy, dAdvz, &
                                                    dBx, dBy, dBz, &
                                              divDx, divDy, divDz  )
      !------------------------------------------------------------
      !------------------------------------------------------------
      implicit none
      class(TrackSubCellType) :: this
      type(ParticleTrackingOptionsType),intent(in) :: trackingOptions
      doubleprecision, intent(in)    :: x, y, z, dt
      doubleprecision, intent(inout) :: vx, vy, vz 
      doubleprecision, intent(inout) :: dAdvx, dAdvy, dAdvz 
      doubleprecision, intent(inout) :: dBx, dBy, dBz 
      doubleprecision, intent(inout) :: divDx, divDy, divDz  
      ! local
      doubleprecision :: alphaL, alphaT, dMEff
      !------------------------------------------------------------
      ! Specifications
      !------------------------------------------------------------

      dMEff  = this%SubCellData%dMEff 
      alphaL = this%SubCellData%alphaLH
      alphaT = this%SubCellData%alphaTH

      call this%LinearInterpolationVelocities( x, y, z, vx, vy, vz )
      call this%DispersionDivergence( x, y, z, divDx, divDy, divDz )
      call this%DisplacementRandomDischarge( x, y, z, alphaL, alphaT, dMEff, dBx, dBy, dBz )
      call this%AdvectionDisplacement( x, y, z, dt, vx, vy, vz, dAdvx, dAdvy, dAdvz )

      return


  end subroutine pr_RWPTDisplacementsLinear


  subroutine pr_RWPTDisplacementsNonlinear(this, x, y, z, vx, vy, vz, &
                                                 dt, trackingOptions, &
                                                 dAdvx, dAdvy, dAdvz, &
                                                       dBx, dBy, dBz, &
                                                 divDx, divDy, divDz  )
      !------------------------------------------------------------
      !------------------------------------------------------------
      implicit none
      class(TrackSubCellType) :: this
      type(ParticleTrackingOptionsType),intent(in) :: trackingOptions
      doubleprecision, intent(in)    :: x, y, z, dt
      doubleprecision, intent(inout) :: vx, vy, vz 
      doubleprecision, intent(inout) :: dAdvx, dAdvy, dAdvz 
      doubleprecision, intent(inout) :: dBx, dBy, dBz 
      doubleprecision, intent(inout) :: divDx, divDy, divDz  
      ! local
      doubleprecision :: Daqueous, Dmol, betaL, betaT
      doubleprecision :: mediumDistance, mediumDelta
      doubleprecision :: alphaL, alphaT
      !------------------------------------------------------------
      ! Specifications
      !------------------------------------------------------------
        
      ! Will consider by convention that molecular diffusion 
      ! specified at configuration file is the aqueous

      ! THIS IS TEMPORARY:
      Daqueous       = trackingOptions%Dmol 
      Dmol           = Daqueous*this%SubCellData%Porosity ! Pore diffusion approx Daq*phi
      mediumDistance = trackingOptions%mediumDistance
      mediumDelta    = trackingOptions%mediumDelta
      betaL          = trackingOptions%betaLong 
      betaT          = trackingOptions%betaTrans


      alphaL = 0d0
      alphaT = 0d0


      ! Nonlinear dispersion
      call this%LinearInterpolationVelocities( x, y, z, vx, vy, vz )
      ! Compute dispersivities
      call pr_ComputeNonlinearDispersivities( this, vx, vy, vz, Daqueous, &
                mediumDistance, mediumDelta, betaL, betaT, alphaL, alphaT )
      call this%DispersionDivergence( x, y, z, divDx, divDy, divDz )
      call this%DisplacementRandomDischarge( x, y, z, alphaL, alphaT, Dmol, dBx, dBy, dBz )
      call this%AdvectionDisplacement( x, y, z, dt, vx, vy, vz, dAdvx, dAdvy, dAdvz )


      return


  end subroutine pr_RWPTDisplacementsNonlinear


  subroutine pr_RWPTDisplacementsAxisymmetric(this, x, y, z, vx, vy, vz, &
                                                    dt, trackingOptions, &
                                                    dAdvx, dAdvy, dAdvz, &
                                                          dBx, dBy, dBz, &
                                                    divDx, divDy, divDz  )
      !------------------------------------------------------------
      ! The dispersion model from Lichtner et al. 2002
      !------------------------------------------------------------
      implicit none
      class(TrackSubCellType) :: this
      type(ParticleTrackingOptionsType),intent(in) :: trackingOptions
      doubleprecision, intent(in)    :: x, y, z, dt
      doubleprecision, intent(inout) :: vx, vy, vz 
      doubleprecision, intent(inout) :: dAdvx, dAdvy, dAdvz 
      doubleprecision, intent(inout) :: dBx, dBy, dBz 
      doubleprecision, intent(inout) :: divDx, divDy, divDz  
      ! local
      doubleprecision :: alphaL, alphaT, dMEff
      doubleprecision :: cosinesq, vnorm
      !------------------------------------------------------------
      ! Specifications
      !------------------------------------------------------------

      dMEff  = this%SubCellData%dMEff 
      
      ! velocities and dispersivities  
      call this%LinearInterpolationVelocities( x, y, z, vx, vy, vz )
      vnorm    = sqrt(vx**2d0+vy**2d0+vz**2d0)
      cosinesq = 0d0
      if ( vnorm .gt. 0d0 ) cosinesq = vz**2d0/vnorm 
      alphaL = &
        this%SubCellData%alphaLH + cosinesq*(this%SubCellData%alphaLV-this%SubCellData%alphaLH)
      alphaT = &
        this%SubCellData%alphaTV + cosinesq*(this%SubCellData%alphaTH-this%SubCellData%alphaTV)

      call this%DispersionDivergence( x, y, z, divDx, divDy, divDz )
      call this%DisplacementRandomDischargeAxisymmetric( x, y, z, & 
        alphaL, alphaT, this%SubCellData%alphaTH, dMEff, dBx, dBy, dBz )
      call this%AdvectionDisplacement( x, y, z, dt, vx, vy, vz, dAdvx, dAdvy, dAdvz )

      return


  end subroutine pr_RWPTDisplacementsAxisymmetric



!-------------------------------------------------------------------
  subroutine pr_ExecuteRandomWalkParticleTracking(this,stopIfNoExit, &
          initialLocation,maximumTime,trackingResult,trackingOptions)
  !------------------------------------------------------------
  ! Function moves particles following RWPT protocol.
  ! It has been adapted from ExecuteTracking but logic
  ! differs significantly as particles are moved numerically 
  ! instead of semi-analytically
  !------------------------------------------------------------
  ! Specifications
  !------------------------------------------------------------
  implicit none
  class(TrackSubCellType) :: this
  logical,intent(in) :: stopIfNoExit
  type(ParticleLocationType),intent(in) :: initialLocation
  doubleprecision,intent(in) :: maximumTime
  type(TrackSubCellResultType),intent(inout) :: trackingResult
  integer :: cellNumber
  doubleprecision :: initialX,initialY,initialZ,initialTime
  doubleprecision :: vx,vy,vz,dt
  doubleprecision :: t,x,y,z
  integer :: exitFace
  ! RWPT
  type(ParticleTrackingOptionsType),intent(in) :: trackingOptions
  doubleprecision :: divDx, divDy, divDz
  doubleprecision :: dAdvx, dAdvy, dAdvz
  doubleprecision :: dBx, dBy, dBz
  doubleprecision :: dx, dy, dz
  doubleprecision :: nx, ny, nz
  doubleprecision :: xi, yi, zi
  doubleprecision :: dxrw, dyrw, dzrw
  logical         :: continueTimeLoop
  logical         :: reachedMaximumTime
  doubleprecision :: tinit, dtold, dtcell
  doubleprecision, dimension(3) :: dtxyz
  integer :: dtLoopCounter, posRestartCounter
  integer :: reboundCounter, intLoopCounter
  integer, parameter :: maxInterfaceLoopCounter   = 5
  integer, parameter :: maxRestartPositionCounter = 10
  doubleprecision, parameter :: maxRelativeJump   = 0.5
  doubleprecision, parameter :: dtReductionFactor = 0.1
  !------------------------------------------------------------

    ! Initialize trackingResult
    call trackingResult%Reset()
    
    cellNumber = initialLocation%CellNumber
    initialX = initialLocation%LocalX
    initialY = initialLocation%LocalY
    initialZ = initialLocation%LocalZ
    initialTime = initialLocation%TrackingTime
    
    trackingResult%CellNumber = cellNumber
    trackingResult%Row = this%SubCellData%Row
    trackingResult%Column = this%SubCellData%Column
    trackingResult%InitialLocation%CellNumber = cellNumber
    trackingResult%InitialLocation%LocalX = initialX
    trackingResult%InitialLocation%LocalY = initialY
    trackingResult%InitialLocation%LocalZ = initialZ
    trackingResult%InitialLocation%TrackingTime = initialTime
    trackingResult%FinalLocation%LocalX = initialX
    trackingResult%FinalLocation%LocalY = initialY
    trackingResult%FinalLocation%LocalZ = initialZ
    trackingResult%FinalLocation%TrackingTime = initialTime
    trackingResult%MaximumTime = maximumTime
    trackingResult%Status = trackingResult%Status_Undefined()
    
    if(stopIfNoExit) then
      if(.not. this%SubCellData%HasExitFace()) then
        trackingResult%Status = trackingResult%Status_NoExitPossible()
        return
      end if
    end if

    ! Advection model pointers are assigned in InitializeRandomWalk 
    ! called at initialize tracking engine
    
    ! Local copies of cell size
    dx = this%SubCellData%DX
    dy = this%SubCellData%DY
    dz = this%SubCellData%DZ

    ! Initialize positions
    x  = initialLocation%LocalX
    y  = initialLocation%LocalY
    z  = initialLocation%LocalZ
    nx = initialLocation%LocalX
    ny = initialLocation%LocalY
    nz = initialLocation%LocalZ


    if ( this%SubCellData%dry ) then
      ! If cell is completely dry, then particle is not displaced
      ! Q: If particle is set to InactiveCell, 
      ! can be displaced in a later cycle if cell is rewetted ?
      ! A: in MPathRW.f90 there is a verification. If cell 
      ! is partially dry, particle status is set to active for tracking
      trackingResult%Status = trackingResult%Status_InactiveCell()
      trackingResult%FinalLocation%CellNumber = cellNumber
      trackingResult%FinalLocation%LocalX = x
      trackingResult%FinalLocation%LocalY = y
      trackingResult%FinalLocation%LocalZ = z
      trackingResult%FinalLocation%TrackingTime = t
      return
    end if 

    ! In case something needs to be done for partially dry cells
    if ( this%SubCellData%partiallyDry ) then
      continue
    end if 

    ! Initialize displacements
    dxrw = 0d0
    dyrw = 0d0
    dzrw = 0d0

    ! Compute time step for RWPT
    call this%ComputeRandomWalkTimeStep( trackingOptions, dtcell )

    ! Initializes current time
    t     = initialTime
    tinit = initialTime
    dtold = dtcell

    ! Something wrong, leave
    if ( dtcell .le. 0d0 ) then 
      trackingResult%ExitFace = exitFace
      trackingResult%Status = trackingResult%Status_Undefined()
      trackingResult%FinalLocation%CellNumber = cellNumber
      trackingResult%FinalLocation%LocalX = x
      trackingResult%FinalLocation%LocalY = y
      trackingResult%FinalLocation%LocalZ = z
      trackingResult%FinalLocation%TrackingTime = t
      return
    end if

    ! Sanity counters
    dtLoopCounter     = 0
    intLoopCounter    = 0
    posRestartCounter = 0

    ! Local cell time loop 
    exitFace = 0
    continueTimeLoop = .true.
    reachedMaximumTime = .false.
    do while( continueTimeLoop )

      ! Time loop counter
      dtLoopCounter = dtLoopCounter + 1

      ! Initialize dt with dtcell
      dt = dtcell

      ! Update current time
      t = t + dt

      ! Recompute dt for maximumTime 
      if (maximumTime .lt. t) then
        t  = t - dt
        dt = max( maximumTime - t, 0d0 ) ! avoids numerical error 
        t  = maximumTime
        reachedMaximumTime = .true.
      end if

      ! Compute RWPT terms
      call this%ComputeRWPTDisplacements( &
                     x, y, z, vx, vy, vz, &
                     dt, trackingOptions, &
                     dAdvx, dAdvy, dAdvz, &
                           dBx, dBy, dBz, &
                     divDx, divDy, divDz  )

      ! RW displacements
      dxrw = dAdvx + divDx*dt + dBx*sqrt( dt )
      dyrw = dAdvy + divDy*dt + dBy*sqrt( dt )
      dzrw = dAdvz + divDz*dt + dBz*sqrt( dt )

      ! Reduce dt if large relative jumps
      do while ( & 
        ( abs(dxrw/dx) .gt. maxRelativeJump ) .or. &
        ( abs(dyrw/dy) .gt. maxRelativeJump ) .or. &
        ( abs(dzrw/dz) .gt. maxRelativeJump ) )

        ! Rollback time and reduce time step
        t  = t - dt
        dt = dtReductionFactor*dt

        ! Given new dt, recompute RW displacements
        call this%AdvectionDisplacement( x, y, z, dt, vx, vy, vz, & 
                                              dAdvx, dAdvy, dAdvz )
        dxrw = dAdvx + divDx*dt + dBx*sqrt( dt )
        dyrw = dAdvy + divDy*dt + dBy*sqrt( dt )
        dzrw = dAdvz + divDz*dt + dBz*sqrt( dt )
        
        ! and update time 
        t = t + dt

        ! Disable the flag for maximum time if it was set
        if ( reachedMaximumTime ) then
          reachedMaximumTime = .false.
        end if

      end do

      ! new positions
      nx   = x + dxrw/dx
      ny   = y + dyrw/dy
      nz   = z + dzrw/dz

      ! particleLeavingCell:
      ! Detect if particle leaving the cell
      ! and force the particle into exactly one
      ! interface by computing the required dt
      intLoopCounter = 0
      do while (                                       &
          ( nx .gt. 1.0d0 ) .or. ( nx .lt. 0d0 )  .or. &
          ( ny .gt. 1.0d0 ) .or. ( ny .lt. 0d0 )  .or. &
          ( nz .gt. 1.0d0 ) .or. ( nz .lt. 0d0 )       & 
      )

          intLoopCounter = intLoopCounter + 1
          dtxyz(:) = 0d0

          ! Recompute dt for exact interface
          call this%ExitFaceAndUpdateTimeStep( x, y, z, nx, ny, nz, &
                                   vx, vy, vz, divDx, divDy, divDz, &
                              dBx, dBy, dBz, t, dt, dtxyz, exitFace )

          ! Health control
          if ( intLoopCounter .gt. maxInterfaceLoopCounter ) then
            ! Restart coordinates 
            nx = initialLocation%LocalX
            ny = initialLocation%LocalY
            nz = initialLocation%LocalZ
            t  = tinit
            dt = dtcell
            posRestartCounter = posRestartCounter + 1
            if ( posRestartCounter .gt. maxRestartPositionCounter ) then 
              ! Something wrong, leave
              trackingResult%ExitFace = 0
              trackingResult%Status = trackingResult%Status_Undefined()
              trackingResult%FinalLocation%CellNumber = cellNumber
              trackingResult%FinalLocation%LocalX = x
              trackingResult%FinalLocation%LocalY = y
              trackingResult%FinalLocation%LocalZ = z
              trackingResult%FinalLocation%TrackingTime = t
              return
            end if
            ! exit interface loop
            exit
          end if

          ! Given new dt, recompute advection displacements
          call this%AdvectionDisplacement( x, y, z, dt, vx, vy, vz, & 
                                                dAdvx, dAdvy, dAdvz )

          ! If maximumTime was reached, but particle left
          ! the cell, then the condition is reset
          if (reachedMaximumTime) then
            reachedMaximumTime = .false.
          end if

          ! Find new RWPT displacements
          select case (exitFace)
            ! X Face
            case(1,2)
              nx   = 1.0d0
              if ( exitFace .eq. 1 ) nx=0d0
              dyrw = dAdvy + divDy*dt + dBy*sqrt( dt )
              ny   = y + dyrw/dy
              dzrw = dAdvz + divDz*dt + dBz*sqrt( dt )
              nz   = z + dzrw/dz
            ! Y Face
            case(3,4)
              dxrw = dAdvx + divDx*dt + dBx*sqrt( dt )
              nx   = x + dxrw/dx
              ny   = 1.0d0
              if ( exitFace .eq. 3 ) ny=0d0
              dzrw = dAdvz + divDz*dt + dBz*sqrt( dt )
              nz   = z + dzrw/dz
            ! Z Face
            case(5,6)
              dxrw = dAdvx + divDx*dt + dBx*sqrt( dt )
              nx   = x + dxrw/dx
              dyrw = dAdvy + divDy*dt + dBy*sqrt( dt )
              ny   = y + dyrw/dy
              nz   = 1.0d0
              if ( exitFace .eq. 5 ) nz=0d0
            ! No exit
            case(0)
              ! Fallback
              ! If exitFace .eq. 0, restart coordinates 
              nx = initialLocation%LocalX
              ny = initialLocation%LocalY
              nz = initialLocation%LocalZ
              t  = tinit
              dt = dtcell
              posRestartCounter = posRestartCounter + 1
              if ( posRestartCounter .gt. maxRestartPositionCounter ) then 
                ! Something wrong, leave
                trackingResult%ExitFace = 0
                trackingResult%Status = trackingResult%Status_Undefined()
                trackingResult%FinalLocation%CellNumber = cellNumber
                trackingResult%FinalLocation%LocalX = x
                trackingResult%FinalLocation%LocalY = y
                trackingResult%FinalLocation%LocalZ = z
                trackingResult%FinalLocation%TrackingTime = t
                return
              end if
              ! Exit interface loop and try again,
              ! new random displacements
              exit
          end select

          ! Think how to integrate these conditions into single 
          ! interface loop.  

          ! Found proper interface ?
          !
          ! It is possible that nx,ny,nz are not consistent values
          ! after previous displacement "to the interface".
          ! This effect is maybe originated due to changing direction
          ! of the displacement vector for different dts. So, even after 
          ! finding a time step for exitFace, this new time step may lead
          ! to landing outside the cell from an orthogonal direction, 
          ! typical case of cell corners. If that is the case, then 
          ! interface loop continues, but now finding a smaller time step
          ! for the "real" crossing.
          if (                                               &
                ( nx .gt. 1.0d0 ) .or. ( nx .lt. 0d0 )  .or. &
                ( ny .gt. 1.0d0 ) .or. ( ny .lt. 0d0 )  .or. &
                ( nz .gt. 1.0d0 ) .or. ( nz .lt. 0d0 )       & 
            ) then 
            ! Continue looking exact interface
            continue
          else
            ! Once the interface is consistent apply boundary conditions

            ! Consider a flag indicating whether the center cell is 
            ! connected to any boundary 

            ! elasticRebound:
            ! Verify if particle has to rebound against boundary face.
            ! Note: one of the possible outputs from a rebound boundary relies on 
            ! the boundary processing being inside the interface detection loop
            reboundCounter = 0
            do while( ( exitFace .gt. 0 ) )

              ! If not connected to rebound boundary cell, leave
              if ( this%SubCellData%MassBoundary(exitFace) .ne. 1 ) then 
                exit
              end if

              ! reboundCounter and a catch for unexpected cases
              reboundCounter = reboundCounter + 1
              if ( reboundCounter .gt. maxRestartPositionCounter ) then
                ! If particle has been rebounding for a long time, stop
                trackingResult%ExitFace = 0
                trackingResult%Status   = trackingResult%Status_Undefined()
                trackingResult%FinalLocation%CellNumber = cellNumber
                trackingResult%FinalLocation%LocalX = x
                trackingResult%FinalLocation%LocalY = y
                trackingResult%FinalLocation%LocalZ = z
                trackingResult%FinalLocation%TrackingTime = t
                return
              end if 

              ! Save interface positions 
              xi = nx
              yi = ny
              zi = nz
              
              ! If dt .eq. 0d0 then the particle is exactly at the interface,
              ! and is a rebound interface. In the meantime, restart.
              if ( dt .eq. 0d0 ) then
                ! Restart coordinates 
                nx = initialLocation%LocalX
                ny = initialLocation%LocalY
                nz = initialLocation%LocalZ
                t  = t - dt
                dt = dtcell
                exitFace = 0
                posRestartCounter = posRestartCounter + 1
                ! With nx, ny, nz adopting the starting values, 
                ! interface loop is also broken
                ! Exit rebound loop 
                exit
              end if
              
              ! Update current time with same time step 
              ! determined for interface displacement 
              t = t + dt

              ! In case the new time for an elastic rebound 
              ! is higher than maximum time, updates dt
              if (maximumTime .lt. t) then
                ! Note: if maximumTime is reached in second half of 
                ! rebound, rebound is shortened

                ! Recompute time step 
                t  = t - dt
                dt = max( maximumTime - t, 0d0 ) ! avoids numerical error
                t  = maximumTime
                reachedMaximumTime = .true.

                ! Recompute advection displacement for the new time step, 
                ! as if starting from original position
                call this%AdvectionDisplacement( x, y, z, dt, vx, vy, vz, dAdvx, dAdvy, dAdvz )

                ! Recompute RW displacements for new dt
                dxrw = dAdvx + divDx*dt + dBx*sqrt( dt )
                dyrw = dAdvy + divDy*dt + dBy*sqrt( dt )
                dzrw = dAdvz + divDz*dt + dBz*sqrt( dt )

                ! If particle lands outside cell, entering interface loop will modify
                ! reachedMaximumTime back to .false.
                !exitFace = 0 ! It is used later so don't reset yet
              end if ! maximumTime 

              ! Particle rebounds with elastic reflection
              select case ( exitFace )
                ! X Face
                case(1,2)
                  nx = nx - dxrw/dx
                  ny = ny + dyrw/dy
                  nz = nz + dzrw/dz
                ! Y Face
                case(3,4)
                  nx = nx + dxrw/dx
                  ny = ny - dyrw/dy
                  nz = nz + dzrw/dz
                ! Z Face
                case(5,6)
                  nx = nx + dxrw/dx
                  ny = ny + dyrw/dy
                  nz = nz - dzrw/dz
              end select

              ! If nx, ny or nz are outside cell interfaces, 
              ! then the interface loop will continue.
              ! Set initial particle position to rebound interface
              if (                                             &
                  ( nx .gt. 1.0d0 ) .or. ( nx .lt. 0d0 )  .or. &
                  ( ny .gt. 1.0d0 ) .or. ( ny .lt. 0d0 )  .or. &
                  ( nz .gt. 1.0d0 ) .or. ( nz .lt. 0d0 )       & 
              ) then

                ! Set starting position to rebound interface
                x = xi
                y = yi
                z = zi

                ! As particle is at the rebound interface, reverts
                ! displacements of the rebound direction. 
                ! These are going to be used in determining 
                ! time step and exact exit position at interface loop
                select case ( exitFace )
                  ! X Face
                  case(1,2)
                    vx    = -vx
                    divDx = -divDx
                    dBx   = -dBx
                  ! Y Face
                  case(3,4)
                    vy    = -vy
                    divDy = -divDy
                    dBy   = -dBy
                  ! Z Face
                  case(5,6)
                    vz    = -vz
                    divDz = -divDz
                    dBz   = -dBz
                end select 

                ! Go to: particleLeavingCell
                exitFace = 0

              else
                ! If nx, ny and nz are inside the cell, then no problem
                ! particleLeavingCell loop is broken and will update 
                ! particle position to rebound position and continue
                ! time loop with cell characteristic time step
                dt =  dtcell
                exitFace = 0

                ! If one of the positions is exactly an interface
                if (                                          & 
                    ( nx .eq. 0d0 ) .or. ( nx .eq. 1d0 ) .or. & 
                    ( ny .eq. 0d0 ) .or. ( ny .eq. 1d0 ) .or. & 
                    ( nz .eq. 0d0 ) .or. ( nz .eq. 1d0 )      & 
                ) then
                    if ( nx .eq. 0d0 ) exitFace = 1 
                    if ( nx .eq. 1d0 ) exitFace = 2
                    if ( ny .eq. 0d0 ) exitFace = 3
                    if ( ny .eq. 1d0 ) exitFace = 4
                    if ( nz .eq. 0d0 ) exitFace = 5
                    if ( nz .eq. 1d0 ) exitFace = 6
                end if

              end if ! outsideInterfaces  

            end do ! elasticRebound

          end if ! If found proper interface

      end do ! particleLeavingCell


      ! Report and leave

      ! Time is up
      if (reachedMaximumTime) then
        ! In this case it is allowed for a particle
        ! to reach the maximum time and eventually have an exit face
        trackingResult%Status = trackingResult%Status_ReachedMaximumTime()
        trackingResult%ExitFace = exitFace
        trackingResult%FinalLocation%CellNumber = cellNumber
        trackingResult%FinalLocation%LocalX = nx
        trackingResult%FinalLocation%LocalY = ny
        trackingResult%FinalLocation%LocalZ = nz
        trackingResult%FinalLocation%TrackingTime = t
        continueTimeLoop = .false.
        ! Done
        return
      end if

      ! Particle left cell
      if( exitFace .gt. 0 ) then 
        ! Based on connected cell id, it is determined 
        ! what happens to the particle (see: source/ModpathCellData.f90:FillSubCellDataBuffer ) 
        ! if ExitFaceConnection = 0, domainBoundary
        ! if ExitFaceConnection > 0, anotherCell
        ! if ExitFaceConnection < 0, subCell, usg-grid
        trackingResult%ExitFaceConnection = this%SubCellData%Connection(exitFace)
        if( trackingResult%ExitFaceConnection .lt. 0 ) then
          ! Internal transfer. This is the case for unstructured grid
          trackingResult%Status = trackingResult%Status_ExitAtInternalFace()
        else
          ! Transfer to another cell
          trackingResult%Status = trackingResult%Status_ExitAtCellFace()
        end if
        trackingResult%ExitFace                   = exitFace
        trackingResult%FinalLocation%CellNumber   = cellNumber
        trackingResult%FinalLocation%LocalX       = nx
        trackingResult%FinalLocation%LocalY       = ny
        trackingResult%FinalLocation%LocalZ       = nz
        trackingResult%FinalLocation%TrackingTime = t
        continueTimeLoop = .false.
        ! Done
        return
      end if 

      ! Update particle positions
      x = nx
      y = ny
      z = nz

    end do ! continueTimeLoop 

    ! Done
    return

  end subroutine pr_ExecuteRandomWalkParticleTracking

  ! RWPT
  subroutine pr_ComputeRandomWalkTimeStep( this, trackingOptions, dt )
  !----------------------------------------------------------------
  ! From user defined options and cell properties, 
  ! compute time step  
  ! 
  ! Params:
  !     - trackingOptions : simulation options
  !     - dt              : time step, output
  !----------------------------------------------------------------
  ! Specifications
  !----------------------------------------------------------------
  implicit none
  class (TrackSubCellType) :: this
  ! input
  type(ParticleTrackingOptionsType),intent(in) :: trackingOptions
  ! output
  doubleprecision, intent(inout) :: dt
  ! local
  doubleprecision :: vx1, vx2, vy1, vy2, vz1, vz2 
  doubleprecision :: dx, dy, dz
  doubleprecision :: alphaL, alphaT
  doubleprecision, dimension(2) :: dts
  !----------------------------------------------------------------

    ! Initialize
    dts(:) = 0d0

    ! Make local copies of face velocities for convenience
    vx1 = this%SubCellData%VX1
    vx2 = this%SubCellData%VX2
    vy1 = this%SubCellData%VY1
    vy2 = this%SubCellData%VY2
    vz1 = this%SubCellData%VZ1
    vz2 = this%SubCellData%VZ2
   
    ! Local copies of cell size
    dx = this%SubCellData%DX
    dy = this%SubCellData%DY
    dz = this%SubCellData%DZ

    ! Local copies of dispersivities
    ! Note: simplified form taking only H values 
    alphaL = this%SubCellData%alphaLH
    alphaT = this%SubCellData%alphaTH

    ! Missing diffusion

    ! Compute time step
    select case (trackingOptions%timeStepKind)
      case (1)
        ! Advection criteria
        dt = trackingOptions%timeStepParameters(1)/( &
            max(abs(vx1), abs(vx2))/dx +             &
            max(abs(vy1), abs(vy2))/dy +             &
            max(abs(vz1), abs(vz2))/dz )
      case (2)
        ! Dispersion criteria
        ! dt = c_T dx**2/D
        dt = trackingOptions%timeStepParameters(2)/(        &
                 alphaL*max(abs(vx1), abs(vx2))/( dx**2 ) + & 
                 alphaT*max(abs(vy1), abs(vy2))/( dy**2 ) + &
                 alphaT*max(abs(vz1), abs(vz2))/( dz**2 ) )
      case (3)
        ! Advection condition
        ! dt = CFL dx / v 
        dts(1) = trackingOptions%timeStepParameters(1)/( & 
            max(abs(vx1), abs(vx2))/dx +                 &
            max(abs(vy1), abs(vy2))/dy +                 &
            max(abs(vz1), abs(vz2))/dz )
        ! Dispersion condition
        ! dt = c_T dx**2/D
        dts(2) = trackingOptions%timeStepParameters(2)/(    &
                 alphaL*max(abs(vx1), abs(vx2))/( dx**2 ) + & 
                 alphaT*max(abs(vy1), abs(vy2))/( dy**2 ) + &
                 alphaT*max(abs(vz1), abs(vz2))/( dz**2 ) )
        ! Compute minimum
        dt     = minval( dts, dts > 0 )
      case (4)
        ! Fixed
        dt = trackingOptions%timeStepParameters(1)
    end select

  end subroutine pr_ComputeRandomWalkTimeStep

  ! RWPT
  subroutine pr_LinearInterpolationVelocities( this, x, y, z, vx, vy, vz )
  !----------------------------------------------------------------
  ! Computes velocity at given location using cell face velocities
  ! 
  ! Params:
  !     - x, y, z    : local cell coordinates
  !     - vx, vy, vz : holders for velocities, output 
  !                                       
  !----------------------------------------------------------------
  ! Specifications
  !----------------------------------------------------------------
  implicit none
  class(TrackSubCellType) :: this
  ! input
  doubleprecision, intent(in) :: x, y, z
  ! output
  doubleprecision, intent(out) :: vx, vy, vz
  !----------------------------------------------------------------

    vx = ( 1.0d0 - x )*this%SubCellData%vx1 + x*this%SubCellData%vx2
    vy = ( 1.0d0 - y )*this%SubCellData%vy1 + y*this%SubCellData%vy2
    vz = ( 1.0d0 - z )*this%SubCellData%vz1 + z*this%SubCellData%vz2

  end subroutine pr_LinearInterpolationVelocities

  ! RWPT
  subroutine pr_AdvectionDisplacementExponential( this, x, y, z, dt, vx, vy, vz, dAdvx, dAdvy, dAdvz )
  !----------------------------------------------------------------
  ! Compute advection displacement using exponential analytical
  ! expression
  ! 
  ! Params: 
  !     - x, y, z             : local cell coordinates
  !     - dt                  : time step
  !     - vx, vy, vz          : interpolated velocities 
  !     - dAdvx, dAdvy, dAdvz : holders for displacements
  !
  ! Interface:
  !     - Advection
  !----------------------------------------------------------------
  ! Specifications
  !----------------------------------------------------------------
  implicit none
  class(TrackSubCellType) :: this
  ! input
  doubleprecision, intent(in) :: x, y, z, dt
  doubleprecision, intent(in) :: vx, vy, vz
  ! output
  doubleprecision, intent(inout) :: dAdvx, dAdvy, dAdvz
  ! local
  doubleprecision :: dvxdx, dvydy, dvzdz
  doubleprecision :: dvtol = 1.0d-10
  !----------------------------------------------------------------

    ! Compute displacement
    ! x
    dvxdx = ( this%SubCellData%vx2 - this%SubCellData%vx1 )/this%SubCellData%dx
    if ( ( abs(dvxdx) .gt. dvtol ) ) then
      dAdvx = vx*( exp(dvxdx*dt) - 1.0d0 )/dvxdx
    else
      dAdvx = vx*dt
    end if
    ! y
    dvydy = ( this%SubCellData%vy2 - this%SubCellData%vy1 )/this%SubCellData%dy
    if ( ( abs(dvydy) .gt. dvtol ) ) then
      dAdvy = vy*( exp(dvydy*dt) - 1.0d0 )/dvydy
    else
      dAdvy = vy*dt
    end if
    ! z
    dvzdz = ( this%SubCellData%vz2 - this%SubCellData%vz1 )/this%SubCellData%dz
    if ( ( abs(dvzdz) .gt. dvtol ) ) then
      dAdvz = vz*( exp(dvzdz*dt) - 1.0d0 )/dvzdz
    else
      dAdvz = vz*dt
    end if

  end subroutine pr_AdvectionDisplacementExponential

  ! RWPT
  subroutine pr_AdvectionDisplacementEulerian( this, x, y, z, dt, vx, vy, vz, dAdvx, dAdvy, dAdvz )
  !----------------------------------------------------------------
  ! Compute advection displacement using eulerian approximation 
  ! 
  ! Params: 
  !     - x, y, z             : local cell coordinates
  !     - dt                  : time step
  !     - vx, vy, vz          : interpolated velocities 
  !     - dAdvx, dAdvy, dAdvz : holders for displacements
  !
  ! Interface:
  !     - Advection
  !----------------------------------------------------------------
  ! Specifications
  !----------------------------------------------------------------
  implicit none
  class(TrackSubCellType) :: this
  ! input
  doubleprecision, intent(in) :: x, y, z, dt
  doubleprecision, intent(in) :: vx, vy, vz
  ! output
  doubleprecision, intent(inout) :: dAdvx, dAdvy, dAdvz
  !----------------------------------------------------------------
      
    ! Compute displacement
    ! x
    dAdvx = vx*dt
    ! y
    dAdvy = vy*dt
    ! z
    dAdvz = vz*dt

  end subroutine pr_AdvectionDisplacementEulerian

  ! RWPT
  subroutine pr_DetectExitFaceAndUpdateTimeStepQuadratic( this, x, y, z, nx, ny, nz, & 
               vx, vy, vz, divDx, divDy, divDz, dBx, dBy, dBz, t, dt, dtxyz, exitFace )
  !----------------------------------------------------------------
  ! Detects exit face and computes time step to force particle
  ! into the interface, using analytical quadratic solution method.
  ! Related to eulerian advection model
  ! 
  ! Params:
  !     - x, y, z             : local cell coordinates
  !     - nx, ny, nz          : updated coordinates after RWPT
  !     - vx, vy, vz          : interpolated velocities
  !     - divDx, divDy, divDz : dispersion divergence
  !     - dBx, dBy, dBz       : random dispersion displacements
  !     - t                   : current time
  !     - dt                  : time step
  !     - dtxyz               : array for saving each dt
  !     - exitFace            : holder for exit face
  !  
  ! Interface:
  !     - ExitFaceAndTimeStep
  !----------------------------------------------------------------
  ! Specifications
  !----------------------------------------------------------------
  implicit none
  class (TrackSubCellType) :: this
  ! input
  doubleprecision, intent(in) :: x, y, z
  doubleprecision, intent(in) :: nx, ny, nz
  doubleprecision, intent(in) :: vx, vy, vz
  doubleprecision, intent(in) :: divDx, divDy, divDz
  doubleprecision, intent(in) :: dBx, dBy, dBz
  ! output
  doubleprecision, intent(inout) :: t, dt
  integer, intent(inout)         :: exitFace
  doubleprecision, dimension(3), intent(inout) :: dtxyz
  ! local
  integer :: exitFaceX, exitFaceY, exitFaceZ
  integer :: imindt
  doubleprecision :: dx, dy, dz
  doubleprecision :: dInterface
  doubleprecision :: AFace, BFace, z1, z2, zsqrt, zsqrtarg
  !----------------------------------------------------------------
        
    ! Initialize
    exitFaceX  = 0
    exitFaceY  = 0
    exitFaceZ  = 0

    ! Local copies of cell size 
    dx = this%SubCellData%DX
    dy = this%SubCellData%DY
    dz = this%SubCellData%DZ

    !Reset t, dt will be replaced
    t = t - dt

    ! Leaving through x face
    dInterface = 0d0
    AFace      = 0d0
    BFace      = 0d0
    z1         = 0d0
    z2         = 0d0
    zsqrt      = 0d0 
    zsqrtarg   = 0d0 
    if ( ( nx .gt. 1.0d0 ) .or. ( nx .lt. 0d0 )  ) then

      ! Compute dInterface
      if ( nx .gt. 1.0d0 ) then 
        dInterface = dx*( 1.0d0 - x )
        exitFaceX  = 2
      else 
        dInterface = -dx*x
        exitFaceX  = 1
      end if

      ! Exactly at interface, force new cell
      if ( dInterface .eq. 0d0 ) then
        exitFace = exitFaceX
        dt = 0.0 
        return
      end if

      ! Coefficients
      AFace = vx + divDx
      BFace = dBx

      ! Given dInterface, compute new dt
      zsqrtarg = BFace**2 + 4*dInterface*AFace

      if ( ( zsqrtarg .ge. 0d0 ) .and. ( AFace .ne. 0d0 ) ) then 
        zsqrt = sqrt( zsqrtarg )
        z1    = (-BFace + zsqrt )/( 2*AFace )
        z2    = (-BFace - zsqrt )/( 2*AFace )

        if ( any( (/ z1, z2 /) .gt. 0d0) ) then 
          ! Compute new dt
          dtxyz(1) = minval( (/z1, z2/), mask=(/z1, z2/)>0d0 )**2
          ! If computed dt is zero, then 
          ! particle is at the interface
          if ( dtxyz(1) .eq. 0d0 ) then
            exitFace = exitFaceX
            dt = 0.0 
            return
          end if
        else if ( any( (/ z1, z2 /) .eq. 0d0) ) then 
          ! This is equivalent to being at the interface, 
          ! It didn't detected any solution gt 0, but one zero,
          ! which means interface
          exitFace = exitFaceX
          dt       = 0.0
          return
        end if

      else if ( AFace .eq. 0d0 ) then 
        ! If A is zero, then the equation is linear
        ! At this point, it is important to note
        ! that dInterface is non zero, but could be really small
        if ( BFace .ne. 0d0 ) then 
          dtxyz(1) = ( dInterface/BFace )**2
          ! If computed dt is zero, then 
          ! particle is at the interface
          if ( dtxyz(1) .eq. 0d0 ) then
            exitFace = exitFaceX
            dt = 0.0 
            return
          end if
        end if 
      end if

      ! If computed dt is higher than current
      ! is not valid, set to zero
      if ( dtxyz(1) .gt. dt ) dtxyz(1) = 0d0

    end if

    ! Leaving through y face
    dInterface = 0d0
    AFace      = 0d0
    BFace      = 0d0
    z1         = 0d0
    z2         = 0d0
    zsqrt      = 0d0 
    zsqrtarg   = 0d0 
    if ( ( ny .gt. 1.0d0 ) .or. ( ny .lt. 0d0 )  ) then

      ! Compute dInterface
      if ( ny .gt. 1.0d0 ) then 
        dInterface = dy*( 1.0d0 - y )
        exitFaceY  = 4
      else 
        dInterface = -dy*y
        exitFaceY  = 3
      end if

      ! Exactly at interface, force new cell
      if ( dInterface .eq. 0.0 ) then
        exitFace = exitFaceY
        dt = 0.0
        return
      end if

      ! Coefficients
      AFace = vy + divDy
      BFace = dBy

      ! Given dInterface, compute new dt
      zsqrtarg = BFace**2 + 4*dInterface*AFace

      if ( ( zsqrtarg .ge. 0d0 ) .and. ( AFace .ne. 0d0 ) ) then
        zsqrt = sqrt( zsqrtarg )
        z1    = (-BFace + zsqrt )/( 2*AFace )
        z2    = (-BFace - zsqrt )/( 2*AFace )

        if ( any( (/ z1, z2 /) .gt. 0d0) ) then 
          ! Compute new dt
          dtxyz(2) = minval( (/z1, z2/), mask=(/z1, z2/)>0d0 )**2
          ! If computed dt is zero, then 
          ! particle is at the interface
          if ( dtxyz(2) .eq. 0d0 ) then
            exitFace = exitFaceY
            dt = 0.0 
            return
          end if
        else if ( any( (/ z1, z2 /) .eq. 0d0) ) then 
          ! This is equivalent to being at the interface, 
          ! It didn't detected any solution gt 0, but one zero,
          ! which means interface
          exitFace = exitFaceY
          dt       = 0.0
          return
        end if

      else if ( AFace .eq. 0d0 ) then 
        ! If A is zero, then the equation is linear
        ! At this point, it is important to note
        ! that dInterface is non zero, but could be really small
        if ( BFace .ne. 0d0 ) then 
          dtxyz(2) = ( dInterface/BFace )**2
          ! If computed dt is zero, then 
          ! particle is at the interface
          if ( dtxyz(2) .eq. 0d0 ) then
            exitFace = exitFaceY
            dt = 0.0 
            return
          end if
        end if 
      end if

      ! If computed dt is higher than current
      ! is not valid, set to zero
      if ( dtxyz(2) .gt. dt ) dtxyz(2) = 0d0

    end if

    ! Leaving through z face
    dInterface = 0d0
    AFace      = 0d0
    BFace      = 0d0
    z1         = 0d0
    z2         = 0d0
    zsqrt      = 0d0 
    zsqrtarg   = 0d0 
    if ( ( nz .gt. 1.0d0 ) .or. ( nz .lt. 0d0 )  ) then

      ! Compute dInterface
      if ( nz .gt. 1.0d0 ) then
        dInterface = dz*( 1.0d0 - z )
        exitFaceZ  = 6
      else
        dInterface = -dz*z
        exitFaceZ  = 5
      end if

      ! Exactly at interface, force new cell
      if ( dInterface .eq. 0.0 ) then
        exitFace = exitFaceZ
        dt = 0.0
        return
      end if

      ! Coefficients
      AFace = vz + divDz
      BFace = dBz
      
      ! Given dInterface, compute new dt
      zsqrtarg = BFace**2 + 4*dInterface*AFace

      if ( ( zsqrtarg .ge. 0d0 ) .and. ( AFace .ne. 0d0 ) )  then
        zsqrt = sqrt( zsqrtarg )
        z1    = (-BFace + zsqrt )/( 2*AFace )
        z2    = (-BFace - zsqrt )/( 2*AFace )

        if ( any( (/ z1, z2 /) .gt. 0d0 ) ) then 
          ! Compute new dt
          dtxyz(3) = minval( (/z1, z2/), mask=(/z1, z2/)>0d0 )**2
          ! If computed dt is zero, then 
          ! particle is at the interface
          if ( dtxyz(3) .eq. 0d0 ) then
            exitFace = exitFaceZ
            dt = 0.0 
            return
          end if
        else if ( any( (/ z1, z2 /) .eq. 0d0) ) then 
          ! This is equivalent to being at the interface, 
          ! It didn't detected any solution gt 0, but one zero,
          ! which means interface
          exitFace = exitFaceZ
          dt       = 0.0
          return
        end if

      else if ( AFace .eq. 0d0 ) then 
        ! If A is zero, then the equation is linear
        ! At this point, it is important to note
        ! that dInterface is non zero
        if ( BFace .ne. 0d0 ) then 
          dtxyz(3) = ( dInterface/BFace )**2
          ! If computed dt is zero, then 
          ! particle is at the interface
          if ( dtxyz(3) .eq. 0d0 ) then
            exitFace = exitFaceZ
            dt = 0.0 
            return
          end if
        end if 
      end if

      ! If computed dt is higher than current
      ! is not valid, set to zero
      if ( dtxyz(3) .gt. dt ) dtxyz(3) = 0d0 

    end if

    ! If all dts where zero, restore t + dt
    ! and leave
    if ( all( dtxyz .eq. 0d0 ) ) then 
      t = t + dt
      return
    end if

    ! Find minimum dt and 
    ! assign values accordingly
    imindt = minloc( dtxyz, dim=1, mask=(dtxyz > 0d0) )
    dt     = dtxyz(imindt)

    if ( imindt .eq. 1 ) then
      exitFace = exitFaceX
    else if ( imindt .eq. 2 ) then
      exitFace = exitFaceY
    else if ( imindt .eq. 3 ) then
      exitFace = exitFaceZ
    end if

    ! Update time
    t = t + dt

  end subroutine pr_DetectExitFaceAndUpdateTimeStepQuadratic


  subroutine pr_DetectExitFaceAndUpdateTimeStepNewton( this, x, y, z, nx, ny, nz, & 
            vx, vy, vz, divDx, divDy, divDz, dBx, dBy, dBz, t, dt, dtxyz, exitFace )
  !----------------------------------------------------------------
  ! Detects exit face and computes time step to force particle
  ! into the interface, using newton raphson method. 
  ! Related to exponential advection model 
  ! 
  ! Params:
  !     - x, y, z             : local cell coordinates
  !     - nx, ny, nz          : updated coordinates after RWPT
  !     - vx, vy, vz          : interpolated velocities
  !     - divDx, divDy, divDz : dispersion divergence
  !     - dBx, dBy, dBz       : random dispersion displacements
  !     - t                   : current time
  !     - dt                  : time step
  !     - dtxyz               : array for saving each dt
  !     - exitFace            : holder for exit face
  !  
  ! Interface:
  !     - ExitFaceAndTimeStep
  !----------------------------------------------------------------
  ! Specifications
  !----------------------------------------------------------------
  implicit none
  class (TrackSubCellType) :: this
  ! input
  doubleprecision, intent(in) :: x, y, z
  doubleprecision, intent(in) :: nx, ny, nz
  doubleprecision, intent(in) :: vx, vy, vz
  doubleprecision, intent(in) :: divDx, divDy, divDz
  doubleprecision, intent(in) :: dBx, dBy, dBz
  ! output
  doubleprecision, intent(inout) :: t, dt
  integer, intent(inout)         :: exitFace
  doubleprecision, dimension(3), intent(inout) :: dtxyz
  ! local
  integer :: exitFaceX, exitFaceY, exitFaceZ
  integer :: imindt
  doubleprecision :: dx, dy, dz
  doubleprecision :: dInterface
  !----------------------------------------------------------------

    ! Initialize
    exitFaceX  = 0
    exitFaceY  = 0
    exitFaceZ  = 0

    ! Local copies of cell size 
    dx = this%SubCellData%DX
    dy = this%SubCellData%DY
    dz = this%SubCellData%DZ

    !Reset t, dt will be replaced
    t = t - dt

    ! Leaving through x face
    dInterface = 0d0
    if ( ( nx .gt. 1.0d0 ) .or. ( nx .lt. 0d0 )  ) then

      ! Compute dInterface
      if ( nx .gt. 1.0d0 ) then 
        dInterface = dx*( 1.0d0 - x )
        exitFaceX  = 2
      else 
        dInterface = -dx*x
        exitFaceX  = 1
      end if

      ! Exactly at interface, force new cell
      if ( dInterface .eq. 0d0 ) then
        exitFace = exitFaceX
        dt = 0d0
        return
      end if

      ! Solve
      call this%NewtonRaphsonTimeStep( dt, vx,&
                         this%SubCellData%vx1,&
                         this%SubCellData%vx2,&
                         this%SubCellData%dx, &
                         dInterface, divDx,   & 
                         dBx, dtxyz(1) )

      ! If computed dt is higher than current
      ! is not valid, set to zero
      if ( dtxyz(1) .gt. dt ) dtxyz(1) = 0d0

    end if


    ! Leaving through y face
    dInterface = 0d0
    if ( ( ny .gt. 1.0d0 ) .or. ( ny .lt. 0d0 )  ) then

      ! Compute dInterface
      if ( ny .gt. 1.0d0 ) then 
        dInterface = dy*( 1.0d0 - y )
        exitFaceY  = 4
      else 
        dInterface = -dy*y
        exitFaceY  = 3
      end if

      ! Exactly at interface, force new cell
      if ( dInterface .eq. 0d0 ) then
        exitFace = exitFaceY
        dt = 0d0
        return
      end if

      ! Solve
      call this%NewtonRaphsonTimeStep( dt, vy,&
                         this%SubCellData%vy1,&
                         this%SubCellData%vy2,&
                         this%SubCellData%dy, &
                         dInterface, divDy,   & 
                         dBy, dtxyz(2) )

      ! If computed dt is higher than current
      ! is not valid, set to zero
      if ( dtxyz(2) .gt. dt ) dtxyz(2) = 0d0

    end if


    ! Leaving through z face
    dInterface = 0d0
    if ( ( nz .gt. 1.0d0 ) .or. ( nz .lt. 0d0 )  ) then

      ! Compute dInterface
      if ( nz .gt. 1.0d0 ) then
        dInterface = dz*( 1.0d0 - z )
        exitFaceZ  = 6
      else
        dInterface = -dz*z
        exitFaceZ  = 5
      end if

      ! Exactly at interface, force new cell
      if ( dInterface .eq. 0d0 ) then
        exitFace = exitFaceZ
        dt = 0d0
        return
      end if

      ! Solve
      call this%NewtonRaphsonTimeStep( dt, vz,&
                         this%SubCellData%vz1,&
                         this%SubCellData%vz2,&
                         this%SubCellData%dz, &
                         dInterface, divDz,   & 
                         dBz, dtxyz(3) )

      ! If computed dt is higher than current
      ! is not valid, set to zero
      if ( dtxyz(3) .gt. dt ) dtxyz(3) = 0d0

    end if


    ! If all dts where zero, restore t + dt
    ! and leave
    if ( all( dtxyz .eq. 0d0 ) ) then 
      t = t + dt
      return
    end if


    ! Find minimum dt and 
    ! assign values accordingly
    imindt = minloc( dtxyz, dim=1, mask=(dtxyz > 0d0) )
    dt     = dtxyz(imindt)
    if ( imindt .eq. 1 ) then
      exitFace = exitFaceX
    else if ( imindt .eq. 2 ) then
      exitFace = exitFaceY
    else if ( imindt .eq. 3 ) then
      exitFace = exitFaceZ
    end if


    ! Update time
    t = t + dt

  end subroutine pr_DetectExitFaceAndUpdateTimeStepNewton


  subroutine pr_NewtonRaphsonVariablesExponential( this, dt, v, v1, v2, &
                                              dx, dInterface, divD, dB, & 
                                      nrf0, nrfprim, nrf2prim, nrf3prim )
  !-----------------------------------------------------------------------------
  ! Specifications
  !-----------------------------------------------------------------------------
  implicit none
  class(TrackSubCellType) :: this
  doubleprecision, intent(in)    :: dt
  doubleprecision, intent(in)    :: v, v1, v2, dx, dInterface, divD, dB
  doubleprecision, intent(inout) :: nrf0, nrfprim, nrf2prim, nrf3prim
  doubleprecision :: dvdx, dAdv
  doubleprecision :: dvtol = 1.0d-10 ! a parameter declaration at module level ? 
  !-----------------------------------------------------------------------------

    if ( dt.eq.0d0 ) return

    ! Initialize
    nrf0      = 0d0
    nrfprim   = 0d0
    nrf2prim  = 0d0
    nrf3prim  = 0d0
    dvdx      = ( v2 - v1 )/dx

    ! Advection
    if ( abs(dvdx) > dvtol ) then 
      dAdv = v*( exp( dvdx*dt ) - 1.0 )/dvdx
    else
      dAdv = v*dt
    end if

    ! Derivatives and iteration terms
    nrf0     = dAdv + divD*dt + dB*sqrt( dt ) - dInterface
    nrfprim  = v*exp( dvdx*dt ) + divD + 0.5*dB/sqrt( dt )
    nrf2prim = dvdx*v*exp( dvdx*dt ) - 0.25*dB*dt**(-1.5)
    nrf3prim = dvdx*dvdx*v*exp( dvdx*dt ) + 0.375*dB*dt**(-2.5)

    ! Done
    return

  end subroutine pr_NewtonRaphsonVariablesExponential


  function pr_GetConvergenceFunction( this, nrf0, nrfprim, nrf2prim ) result( gprim )
  !----------------------------------------------------------------
  ! Specifications
  !----------------------------------------------------------------
  implicit none
  class(TrackSubCellType) :: this
  doubleprecision, intent(in) :: nrf0, nrfprim, nrf2prim
  doubleprecision :: gprim 
  !----------------------------------------------------------------
  
    gprim    = abs( nrf2prim*nrf0 / nrfprim**2 )

    ! Done
    return

  end function pr_GetConvergenceFunction


  function pr_GetConvergenceFunctionDerivative( this, nrf0, nrfprim, nrf2prim, nrf3prim ) result( gprimder )
  !----------------------------------------------------------------
  !----------------------------------------------------------------
  implicit none
  class(TrackSubCellType) :: this
  doubleprecision, intent(in) :: nrf0, nrfprim, nrf2prim, nrf3prim
  doubleprecision :: gprimnoabs, gprimder, gprimprim
  !----------------------------------------------------------------

    gprimnoabs = nrf2prim*nrf0/nrfprim**2
    gprimprim  = -2d0*( nrf0*nrf2prim**2 )/( nrfprim**3 ) + ( nrfprim*nrf2prim + nrf0*nrf3prim )/( nrfprim**2 )
    gprimder   = gprimnoabs/abs( gprimnoabs )*gprimprim

    ! Done      
    return

  end function pr_GetConvergenceFunctionDerivative


  function pr_SetInitialGuess( this, dt, v, v1, v2, dx, dInterface, divD, dB ) result( dtnew )
  !----------------------------------------------------------------
  ! Specifications
  !----------------------------------------------------------------
  implicit none
  class(TrackSubCellType) :: this
  doubleprecision, intent(in) :: dt
  doubleprecision, intent(in) :: v, v1, v2, dx, dInterface, divD, dB
  doubleprecision :: dtnew, dtold, dt0 
  doubleprecision :: nrf0, nrfprim
  doubleprecision :: nrf2prim, nrf3prim
  doubleprecision :: gprimdernew, gprimnew
  doubleprecision :: gprimder, gprim, gprimmin
  doubleprecision :: gamman
  doubleprecision :: dtfraction, dtoldfraction
  doubleprecision :: mindtfrac
  doubleprecision :: gdtol  = 1e-6
  integer         :: maxter = 50
  integer         :: counter
  integer         :: ccount
  integer         :: bcount
  integer         :: zcount
  integer         :: rcount
  integer         :: gcount
  !----------------------------------------------------------------

    ! Check 
    dtnew = 0d0
    if ( dt.eq.0d0 ) return

    ! Initialize
    gprimmin  = 1e10
    gprim     = 1e10

    ! Estimate order of initial dt fraction
    ! Run over very small fractions
    gcount = 0
    do while ( ( gprim .gt. 1 ) .and. ( gcount .lt. 7 ) )
      gcount     = gcount + 1
      dtfraction = 10d0**( -8 + gcount )
      dt0        = dt*dtfraction
      call pr_NewtonRaphsonVariablesExponential( this, dt0, v, v1, v2, &
                                             dx, dInterface, divD, dB, & 
                                     nrf0, nrfprim, nrf2prim, nrf3prim )
      gprim = pr_GetConvergenceFunction( this, nrf0, nrfprim, nrf2prim )
      if ( gprim .lt. gprimmin ) then 
        gprimmin  = gprim
        mindtfrac = dtfraction
      end if 
    end do


    ! Set smaller fraction for gradient descent 
    if ( ( gprim .lt. 1 ) .and. ( gcount .eq. 1 ) ) then  
      dtoldfraction = 1d-8 
    else if ( gprim .lt. 1 ) then 
      dtoldfraction =  10d0**( -8 + gcount - 1 )
    else
      dtfraction    = mindtfrac
      dtoldfraction = 0.5*dtfraction/10d0      
    end if 
   
    ! Set initial dt, dtold
    dtnew  = dtfraction*dt
    dtold  = dtoldfraction*dt 

    ! Do it or leave ? 
    call pr_NewtonRaphsonVariablesExponential( this, dtnew, v, v1, v2, &
                                           dx, dInterface, divD, dB, & 
                                   nrf0, nrfprim, nrf2prim, nrf3prim )
    gprim = pr_GetConvergenceFunction( this, nrf0, nrfprim, nrf2prim )

    if ( gprim .lt. 0.5 ) then 
        ! Done 
        return
    end if

    
    ! Continue to gradient descent  
    gprimnew = 10d0
    gprim    = 10d0
    ccount  = 0
    bcount  = 0
    zcount  = 0
    rcount  = 0
    counter = 0
    do while ( counter .lt. maxter )

        counter = counter + 1 

        ! Compute quantities for gradient descent
        call pr_NewtonRaphsonVariablesExponential( this, dtold, v, v1, v2, &
                                                 dx, dInterface, divD, dB, & 
                                         nrf0, nrfprim, nrf2prim, nrf3prim )
        gprimder  = pr_GetConvergenceFunctionDerivative( this, nrf0, nrfprim, nrf2prim, nrf3prim )

        call pr_NewtonRaphsonVariablesExponential( this, dtnew, v, v1, v2, &
                                                 dx, dInterface, divD, dB, & 
                                         nrf0, nrfprim, nrf2prim, nrf3prim )
        gprimdernew = pr_GetConvergenceFunctionDerivative( this, nrf0, nrfprim, nrf2prim, nrf3prim )
        gprim       = pr_GetConvergenceFunction( this, nrf0, nrfprim, nrf2prim )

        ! Compute gamma gradient descent and update
        if ( ( gprimdernew - gprimder ) .gt. 0d0 ) then
            gamman = abs( ( dtnew - dtold )*( gprimdernew - gprimder ) )/( ( gprimdernew - gprimder )**2 )
        else 
            gamman = 0d0
        end if 
        dtold = dtnew
        dtnew = dtnew - gamman*gprimdernew

        ! If nan, leave loop
        if ( isnan( dtnew ) ) then 
           dtnew = 0d0 
           exit
        end if 

        ! Bound dtnew to values smaller than one or leave if ready
        if ( ( gprim .lt. 1 ) .and. ( dtnew/dt .gt. 1 ) ) then

            dtnew = dtold

            ! Done
            return
      
        else if ( dtnew/dt .gt. 1 ) then 

            dtnew = 0.99*dt

        end if 

        ! Bound dtnew values higher than zero
        ! If convergence parameter already smaller than one
        ! set to dtold and leave.
        if ( ( gprim .lt. 1 ) .and. ( dtnew .lt. 0 ) .and. ( dtold .gt. 0 )  )  then

            dtnew = dtold

            ! Done
            return

        else if ( dtnew .lt. 0 ) then
            ! If smaller than zero but no convergence 
            ! yet, count how many times this happens 
            ! and leave after twice. Probably a very small value
            ! close to zero  

            dtnew = 1e-8*dt 
            zcount = zcount + 1 

            if ( zcount .gt. 1 ) then 
            
                dtnew = 0d0
                
                ! Leave loop
                exit

            end if 

        end if


        call pr_NewtonRaphsonVariablesExponential( this, dtnew, v, v1, v2, &
                                                 dx, dInterface, divD, dB, & 
                                         nrf0, nrfprim, nrf2prim, nrf3prim )
        gprimnew = pr_GetConvergenceFunction( this, nrf0, nrfprim, nrf2prim )


        ! If values appropiate for convergence,
        ! leave after two consecutive occurrences
        if ( ( gprimnew .lt. 1 ) .and. ( dtnew .gt. 0 ) .and. ( dtnew/dt .lt. 1 ) ) then 

            rcount = rcount + 1
            if ( rcount .gt. 1 ) then 
                ! Done
                return
            end if
        else 
            rcount = 0
        end if 


        ! If convergence function is already
        ! less than one and the method provides
        ! a higher value for the new 
        ! convergence function, leave
        if ( ( gprim .lt. 1 ) .and. ( gprimnew .gt. gprim ) ) then

            dtnew = dtold
            
            ! Done
            return

        end if 


        ! If max iterations, exit do loop
        if ( counter .gt. maxter ) then 
            ! Leave loop  
            exit
        end if 


        ! If changes in initial guess
        ! are small, at least twice consecutive,
        ! then leave
        if ( ( abs( (dtold - dtnew)/dtold ) < gdtol ) .and. ( gprim .lt. 1 )  ) then

            ccount = ccount + 1
            if ( ccount .gt. 1 ) then 
                ! Done 
                return 
            end if 
        else 
            ccount = 0
        end if

        
        ! If relative changes in estimate are small and 
        ! convergence function not good, leave and try 
        ! higher orders
        if ( ( abs( (dtold - dtnew)/dtold ) < gdtol ) .and. ( gprim .gt. 1 )  ) then

            bcount = bcount + 1

            if ( bcount .gt. 1 ) then
                ! Leave loop  
                exit
            end if 
    
        else 
            bcount = 0
        end if


    end do


    ! Needed ?
    if ( gprim .lt. 1 ) then 
        ! Done
        return 
    end if 


    ! Continue and try with higher order fractions
    gprim     = 1e10
    gprimmin  = 1e10
    mindtfrac = 0d0

    ! Run over higher order fractions 
    gcount = 1
    do while ( ( gprim .gt. 1 ) .and. ( gcount .lt. 9 ) )

      gcount     = gcount + 1
      dtfraction = 0.1*gcount
      dt0        = dt*dtfraction

      call pr_NewtonRaphsonVariablesExponential( this, dt0, v, v1, v2, &
                                             dx, dInterface, divD, dB, & 
                                     nrf0, nrfprim, nrf2prim, nrf3prim )
      gprim = pr_GetConvergenceFunction( this, nrf0, nrfprim, nrf2prim )
  
      if ( gprim .lt. gprimmin ) then  
          gprimmin  = gprim
          mindtfrac = dtfraction
      end if 

    end do

    ! If any was smaller than one,
    ! set those values and leave
    if ( gprimmin .lt. 1 ) then 
        dtnew = dt*mindtfrac
    end if 

    ! Set advection condition if zero and v, otherwise zero 
    if ( ( dtnew .le. 0d0 ) .and. ( abs( v ) .gt. 0d0 ) ) then 
        dtnew = abs( dInterface/v )
    end if 


    ! Done
    return
        

  end function pr_SetInitialGuess


  subroutine NewtonRaphsonTimeStepExponentialAdvection( this, dt, v, v1, v2, &
                                              dx, dInterface, divD, dB, dtnr )
  !----------------------------------------------------------------
  ! Computes time step required to reach interface by RWPT
  ! with exponential integration of advection, only one axis
  !
  ! Params:
  !     - dt         : time step that moved a particle outside the cell
  !     - v          : interpolated velocity at particle's position
  !     - v1, v2     : cell faces velocity
  !     - dx         : cell size
  !     - dInterface : distance to axis interface
  !     - divD       : dispersion divergence
  !     - dB         : random dispersion displacement
  !     - dtnr       : time step computed with NR, output
  !
  !----------------------------------------------------------------
  ! Specifications
  !----------------------------------------------------------------
  implicit none
  class(TrackSubCellType) :: this
  ! input 
  doubleprecision, intent(in)  :: dt 
  doubleprecision, intent(in)  :: v, v1, v2, dx, dInterface, divD, dB
  ! output
  doubleprecision, intent(out) :: dtnr
  ! local
  doubleprecision :: dvdx, dAdv
  doubleprecision :: dt0
  doubleprecision :: nrf0, nrfprim, nrerror
  doubleprecision :: dvtol = 1.0d-10
  doubleprecision :: nrtol = 1e-6
  integer :: countIter
  integer :: maxIter = 50
  ! Note
  ! Compute 
  !    abs(f*ftwoprim) .lt. abs(fprim**2)
  ! https://math.stackexchange.com/questions/3136446/condition-for-convergence-of-newton-raphson-method
  !----------------------------------------------------------------

      ! Initialize
      nrf0      = 0d0
      nrfprim   = 0d0
      nrerror   = 1d6 ! Something big
      dvdx      = ( v2 - v1 )/dx

      ! Initial time step estimate
      dt0 = pr_SetInitialGuess( this, dt, v, v1, v2, dx, dInterface, divD, dB )

      ! If no proper initial guess, leave
      if ( dt0 .eq. 0d0 ) then 
        dtnr = 0d0 
        return
      end if 

      ! Iteration until convergence or maxIterations
      countIter = 0
      do while( ( abs(nrerror/dt0) .gt. nrtol ) .and. ( countIter .lt. maxIter ) )

        countIter = countIter + 1

        ! Compute displacement, although initially this 
        ! should be known
        if ( ( abs(dvdx) .gt. dvtol ) ) then
          dAdv = v*( exp(dvdx*dt0) - 1.0d0 )/dvdx
        else
          dAdv = v*dt0
        end if

        ! RWPT displacement with initial guess
        nrf0  = dAdv + divD*dt0 + dB*sqrt( dt0 ) - dInterface

        ! Analytical derivative of RWPT displacement
        nrfprim  = v*exp(dvdx*dt0) + divD + 0.5*dB/sqrt(dt0) 

        ! NR error and new time step
        nrerror = -nrf0/nrfprim
        dt0     = dt0 + nrerror 

        ! It can't be smaller than zero
        if ( dt0 .lt. 0d0 ) then
          dt0 = -0.9*dt0   
        end if 

      end do

      ! Assign return value
      dtnr = dt0

      ! If new value higher than the previous, return zero 
      if ( dt0 .gt. dt ) then
        dtnr = dt
        ! Done
        return 
      end if 

      ! If no convergence, return zero
      if ( ( countIter .eq. maxIter ) .and. ( abs(nrerror/dt0) .gt. nrtol ) ) then
        ! Done
        dtnr = 0d0 
        return 
      end if

  end subroutine NewtonRaphsonTimeStepExponentialAdvection

  ! RWPT
  subroutine pr_CornerDispersionLinearIsotropic( this, &
                                             qprod000, & 
                                             qprod100, &
                                             qprod010, &
                                             qprod110, &
                                             qprod001, &
                                             qprod101, &
                                             qprod011, &
                                             qprod111  )
  !-------------------------------------------------------------------
  ! Compute dispersion at the cell corners, using the 
  ! linear isotropic dispersion model 
  !
  ! Interface:
  !   CornerDispersion
  !-------------------------------------------------------------------
  ! Specifications
  !-------------------------------------------------------------------
  implicit none
  class (TrackSubCellType) :: this
  ! input
  doubleprecision, dimension(8), intent(in) :: qprod000
  doubleprecision, dimension(8), intent(in) :: qprod100
  doubleprecision, dimension(8), intent(in) :: qprod010
  doubleprecision, dimension(8), intent(in) :: qprod110
  doubleprecision, dimension(8), intent(in) :: qprod001
  doubleprecision, dimension(8), intent(in) :: qprod101
  doubleprecision, dimension(8), intent(in) :: qprod011
  doubleprecision, dimension(8), intent(in) :: qprod111
  ! local
  doubleprecision :: alphaL, alphaT, dMEff
  doubleprecision :: aLMinusaT
  doubleprecision :: aTQPlusPDMEff
  !-------------------------------------------------------------------

    ! This is ok for uniform dispersion parameters 
    dMEff     = this%SubCellData%dMEff 
    alphaL    = this%SubCellData%alphaLH
    alphaT    = this%SubCellData%alphaTH
    aLMinusaT = alphaL - alphaT

    ! These calculations correspond to the isotropic dispersion model
    ! 000
    aTQPlusPDMeff = alphaT*this%qCorner000(4) + this%porosity000*dMEff 
    this%Dxxx(1)  = ( aTQPlusPDMEff ) + ( aLMinusaT )*qprod000(1)
    this%Dyyy(1)  = ( aTQPlusPDMEff ) + ( aLMinusaT )*qprod000(4)
    this%Dzzz(1)  = ( aTQPlusPDMEff ) + ( aLMinusaT )*qprod000(6)
    this%Dxyx(1)  = ( aLMinusaT )*qprod000(2)
    this%Dxzx(1)  = ( aLMinusaT )*qprod000(3)
    this%Dxyy(1)  = this%Dxyx(1)
    this%Dyzy(1)  = ( aLMinusaT )*qprod000(5)
    this%Dxzz(1)  = this%Dxzx(1)
    this%Dyzz(1)  = this%Dyzy(1)
    ! 100
    aTQPlusPDMeff = alphaT*this%qCorner100(4) + this%porosity100*dMEff 
    this%Dxxx(2)  = ( aTQPlusPDMEff ) + ( aLMinusaT )*qprod100(1)
    this%Dyyy(2)  = ( aTQPlusPDMEff ) + ( aLMinusaT )*qprod100(4)
    this%Dzzz(2)  = ( aTQPlusPDMEff ) + ( aLMinusaT )*qprod100(6)
    this%Dxyx(2)  = ( aLMinusaT )*qprod100(2)
    this%Dxzx(2)  = ( aLMinusaT )*qprod100(3)
    this%Dxyy(2)  = this%Dxyx(2)
    this%Dyzy(2)  = ( aLMinusaT )*qprod100(5)
    this%Dxzz(2)  = this%Dxzx(2)
    this%Dyzz(2)  = this%Dyzy(2)
    ! 010
    aTQPlusPDMeff = alphaT*this%qCorner100(4) + this%porosity010*dMEff 
    this%Dxxx(3)  = ( aTQPlusPDMEff ) + ( aLMinusaT )*qprod010(1)
    this%Dyyy(3)  = ( aTQPlusPDMEff ) + ( aLMinusaT )*qprod010(4)
    this%Dzzz(3)  = ( aTQPlusPDMEff ) + ( aLMinusaT )*qprod010(6)
    this%Dxyx(3)  = ( aLMinusaT )*qprod010(2)
    this%Dxzx(3)  = ( aLMinusaT )*qprod010(3)
    this%Dxyy(3)  = this%Dxyx(3)
    this%Dyzy(3)  = ( aLMinusaT )*qprod010(5)
    this%Dxzz(3)  = this%Dxzx(3)
    this%Dyzz(3)  = this%Dyzy(3)
    ! 110
    aTQPlusPDMeff = alphaT*this%qCorner100(4) + this%porosity110*dMEff 
    this%Dxxx(4)  = ( aTQPlusPDMEff ) + ( aLMinusaT )*qprod110(1)
    this%Dyyy(4)  = ( aTQPlusPDMEff ) + ( aLMinusaT )*qprod110(4)
    this%Dzzz(4)  = ( aTQPlusPDMEff ) + ( aLMinusaT )*qprod110(6)
    this%Dxyx(4)  = ( aLMinusaT )*qprod110(2)
    this%Dxzx(4)  = ( aLMinusaT )*qprod110(3)
    this%Dxyy(4)  = this%Dxyx(4)
    this%Dyzy(4)  = ( aLMinusaT )*qprod110(5)
    this%Dxzz(4)  = this%Dxzx(4)
    this%Dyzz(4)  = this%Dyzy(4)
    ! 001
    aTQPlusPDMeff = alphaT*this%qCorner100(4) + this%porosity001*dMEff 
    this%Dxxx(5)  = ( aTQPlusPDMEff ) + ( aLMinusaT )*qprod001(1)
    this%Dyyy(5)  = ( aTQPlusPDMEff ) + ( aLMinusaT )*qprod001(4)
    this%Dzzz(5)  = ( aTQPlusPDMEff ) + ( aLMinusaT )*qprod001(6)
    this%Dxyx(5)  = ( aLMinusaT )*qprod001(2)
    this%Dxzx(5)  = ( aLMinusaT )*qprod001(3)
    this%Dxyy(5)  = this%Dxyx(5)
    this%Dyzy(5)  = ( aLMinusaT )*qprod001(5)
    this%Dxzz(5)  = this%Dxzx(5)
    this%Dyzz(5)  = this%Dyzy(5)
    ! 101
    aTQPlusPDMeff = alphaT*this%qCorner100(4) + this%porosity101*dMEff 
    this%Dxxx(6)  = ( aTQPlusPDMEff ) + ( aLMinusaT )*qprod101(1)
    this%Dyyy(6)  = ( aTQPlusPDMEff ) + ( aLMinusaT )*qprod101(4)
    this%Dzzz(6)  = ( aTQPlusPDMEff ) + ( aLMinusaT )*qprod101(6)
    this%Dxyx(6)  = ( aLMinusaT )*qprod101(2)
    this%Dxzx(6)  = ( aLMinusaT )*qprod101(3)
    this%Dxyy(6)  = this%Dxyx(6)
    this%Dyzy(6)  = ( aLMinusaT )*qprod101(5)
    this%Dxzz(6)  = this%Dxzx(6)
    this%Dyzz(6)  = this%Dyzy(6)
    ! 011
    aTQPlusPDMeff = alphaT*this%qCorner100(4) + this%porosity011*dMEff 
    this%Dxxx(7)  = ( aTQPlusPDMEff ) + ( aLMinusaT )*qprod011(1)
    this%Dyyy(7)  = ( aTQPlusPDMEff ) + ( aLMinusaT )*qprod011(4)
    this%Dzzz(7)  = ( aTQPlusPDMEff ) + ( aLMinusaT )*qprod011(6)
    this%Dxyx(7)  = ( aLMinusaT )*qprod011(2)
    this%Dxzx(7)  = ( aLMinusaT )*qprod011(3)
    this%Dxyy(7)  = this%Dxyx(7)
    this%Dyzy(7)  = ( aLMinusaT )*qprod011(5)
    this%Dxzz(7)  = this%Dxzx(7)
    this%Dyzz(7)  = this%Dyzy(7)
    ! 111
    aTQPlusPDMeff = alphaT*this%qCorner100(4) + this%porosity111*dMEff 
    this%Dxxx(8)  = ( aTQPlusPDMEff ) + ( aLMinusaT )*qprod111(1)
    this%Dyyy(8)  = ( aTQPlusPDMEff ) + ( aLMinusaT )*qprod111(4)
    this%Dzzz(8)  = ( aTQPlusPDMEff ) + ( aLMinusaT )*qprod111(6)
    this%Dxyx(8)  = ( aLMinusaT )*qprod111(2)
    this%Dxzx(8)  = ( aLMinusaT )*qprod111(3)
    this%Dxyy(8)  = this%Dxyx(8)
    this%Dyzy(8)  = ( aLMinusaT )*qprod111(5)
    this%Dxzz(8)  = this%Dxzx(8)
    this%Dyzz(8)  = this%Dyzy(8)

    ! Done
    return

  end subroutine pr_CornerDispersionLinearIsotropic


  ! RWPT
  subroutine pr_DispersionDivergence( this, x, y, z, divDx, divDy, divDz )
      !----------------------------------------------------------------
      ! Compute dispersion divergence terms 
      ! 
      ! Params:
      !     - x, y, z             : local cell coordinates
      !     - alphaL              : longidutinal dispersivity
      !     - alphaT              : transverse dispersivity
      !     - dMEff               : effective molecular diffusion (corrected by tortuosity)
      !     - divDx, divDy, divDz : dispersion divergence, output 
      !----------------------------------------------------------------
      ! Specifications
      !----------------------------------------------------------------
      implicit none
      class (TrackSubCellType), target :: this
      ! input
      doubleprecision, intent(in) :: x, y, z
      ! output
      doubleprecision, intent(out) :: divDx, divDy, divDz
      ! local
      doubleprecision :: dDxxdx, dDxydy, dDxzdz, &
                         dDxydx, dDyydy, dDyzdz, &
                         dDxzdx, dDyzdy, dDzzdz
      !---------------------------------------------------------------- 

      ! Interpolated derivates !
      call this%TrilinearDerivativeX( x, y, z, &
                this%Dxxx(1), & 
                this%Dxxx(2), &
                this%Dxxx(3), &
                this%Dxxx(4), &
                this%Dxxx(5), &
                this%Dxxx(6), &
                this%Dxxx(7), &
                this%Dxxx(8), &
                dDxxdx )
      call this%TrilinearDerivativeY( x, y, z, &
                this%Dyyy(1), &
                this%Dyyy(2), &
                this%Dyyy(3), &
                this%Dyyy(4), &
                this%Dyyy(5), &
                this%Dyyy(6), &
                this%Dyyy(7), &
                this%Dyyy(8), &
                dDyydy )
      call this%TrilinearDerivativeZ( x, y, z, &
                this%Dzzz(1), &
                this%Dzzz(2), &
                this%Dzzz(3), &
                this%Dzzz(4), &
                this%Dzzz(5), &
                this%Dzzz(6), &
                this%Dzzz(7), &
                this%Dzzz(8), &
                dDzzdz )
      call this%TrilinearDerivativeX( x, y, z, & 
                this%Dxyx(1), & 
                this%Dxyx(2), &
                this%Dxyx(3), &
                this%Dxyx(4), &
                this%Dxyx(5), &
                this%Dxyx(6), &
                this%Dxyx(7), &
                this%Dxyx(8), &
                dDxydx )
      call this%TrilinearDerivativeX( x, y, z, &
                this%Dxzx(1), &
                this%Dxzx(2), &
                this%Dxzx(3), &
                this%Dxzx(4), &
                this%Dxzx(5), &
                this%Dxzx(6), &
                this%Dxzx(7), &
                this%Dxzx(8), &
                dDxzdx )
      call this%TrilinearDerivativeY( x, y, z, &
                this%Dxyy(1), &
                this%Dxyy(2), &
                this%Dxyy(3), &
                this%Dxyy(4), &
                this%Dxyy(5), &
                this%Dxyy(6), &
                this%Dxyy(7), &
                this%Dxyy(8), &
                dDxydy )
      call this%TrilinearDerivativeY( x, y, z, &
                this%Dyzy(1), &
                this%Dyzy(2), &
                this%Dyzy(3), &
                this%Dyzy(4), &
                this%Dyzy(5), &
                this%Dyzy(6), &
                this%Dyzy(7), &
                this%Dyzy(8), &
                dDyzdy )
      call this%TrilinearDerivativeZ( x, y, z, &
                this%Dxzz(1), &
                this%Dxzz(2), &
                this%Dxzz(3), &
                this%Dxzz(4), &
                this%Dxzz(5), &
                this%Dxzz(6), &
                this%Dxzz(7), &
                this%Dxzz(8), &
                dDxzdz )
      call this%TrilinearDerivativeZ( x, y, z, &
                this%Dyzz(1), &
                this%Dyzz(2), &
                this%Dyzz(3), &
                this%Dyzz(4), &
                this%Dyzz(5), &
                this%Dyzz(6), &
                this%Dyzz(7), &
                this%Dyzz(8), &
                dDyzdz )

      ! Notice correction by porosity and retardation 
      divDx = ( dDxxdx + dDxydy + dDxzdz )/this%SubCellData%Porosity/this%SubCellData%Retardation
      divDy = ( dDxydx + dDyydy + dDyzdz )/this%SubCellData%Porosity/this%SubCellData%Retardation
      divDz = ( dDxzdx + dDyzdy + dDzzdz )/this%SubCellData%Porosity/this%SubCellData%Retardation


  end subroutine pr_DispersionDivergence


  subroutine pr_DisplacementRandomDischarge( this, x, y, z, alphaL, alphaT, dMEff, dBx, dBy, dBz ) 
      !----------------------------------------------------------------
      ! Computes the product between displacement matrix and random 
      ! vector
      !
      ! Params:
      !     - x, y, z       : local cell coordinates
      !     - alphaL        : longidutinal dispersivity
      !     - alphaT        : transverse dispersivity
      !     - dMEff         : effective molecular diffusion (corrected by tortuosity )
      !     - dBx, dBy, dBz : random dispersion displacement, output
      !----------------------------------------------------------------
      ! Specifications
      !----------------------------------------------------------------
      implicit none
      class (TrackSubCellType) :: this
      ! input
      doubleprecision, intent(in)    :: x, y, z
      doubleprecision, intent(in)    :: alphaL, alphaT, dMEff
      ! output
      doubleprecision, intent(out) :: dBx, dBy, dBz
      ! local
      doubleprecision :: vBx, vBy, vBz, vBnorm, vBnormxy
      doubleprecision :: B11, B12, B13, B21, B22, B23, B31, B32
      doubleprecision :: rdmx, rdmy, rdmz
      doubleprecision :: RFactor
      doubleprecision, dimension(4) :: v000
      doubleprecision, dimension(4) :: v100
      doubleprecision, dimension(4) :: v010
      doubleprecision, dimension(4) :: v110
      doubleprecision, dimension(4) :: v001
      doubleprecision, dimension(4) :: v101
      doubleprecision, dimension(4) :: v011
      doubleprecision, dimension(4) :: v111
      !----------------------------------------------------------------

      ! Initialize
      dBx = 0d0
      dBy = 0d0
      dBz = 0d0

      ! Local copies of corner velocities
      ! pointers ?
      v000 = this%qCorner000 / this%porosity000 
      v100 = this%qCorner100 / this%porosity100 
      v010 = this%qCorner010 / this%porosity010 
      v110 = this%qCorner110 / this%porosity110 
      v001 = this%qCorner001 / this%porosity001 
      v101 = this%qCorner101 / this%porosity101 
      v011 = this%qCorner011 / this%porosity011 
      v111 = this%qCorner111 / this%porosity111

      ! Extract R
      RFactor = this%SubCellData%Retardation

      ! Trilinear interpolation of velocities and norm
      call this%Trilinear( x, y, z, &
                           v000(1), v100(1), v010(1), v110(1), &
                           v001(1), v101(1), v011(1), v111(1), &
                           vBx )
      call this%Trilinear( x, y, z, &
                           v000(2), v100(2), v010(2), v110(2), &
                           v001(2), v101(2), v011(2), v111(2), &
                           vBy )
      call this%Trilinear( x, y, z, &
                           v000(3), v100(3), v010(3), v110(3), &
                           v001(3), v101(3), v011(3), v111(3), &
                           vBz )
      vBnorm   = sqrt( vBx**2 + vBy**2 + vBz**2 )
      vBnormxy = sqrt( vBx**2 + vBy**2 )
   

      ! Displacement matrix terms
      ! Refs: Fernàndez-Garcia et al. 2005; Salamon et al. 2006
      ! Handles the case of zero vBnorm
      B11 = 0d0 
      B12 = 0d0
      B13 = 0d0
      B21 = 0d0
      B22 = 0d0
      B23 = 0d0
      B31 = 0d0
      B32 = 0d0
      if ( vBnorm .gt. 0d0 ) then
        B11 =       vBx*sqrt( 2*( alphaL*vBnorm + dMEff )/RFactor )/vBnorm
        B21 =       vBy*sqrt( 2*( alphaL*vBnorm + dMEff )/RFactor )/vBnorm
        B31 =       vBz*sqrt( 2*( alphaL*vBnorm + dMEff )/RFactor )/vBnorm
        B32 =  vBnormxy*sqrt( 2*( alphaT*vBnorm + dMEff )/RFactor )/vBnorm
        if ( vBnormxy .gt. 0d0 ) then
          B12 =  -vBx*vBz*sqrt( 2*( alphaT*vBnorm + dMEff )/RFactor )/vBnorm/vBnormxy
          B13 =      -vBy*sqrt( 2*( alphaT*vBnorm + dMEff )/RFactor )/vBnormxy
          B22 =  -vBy*vBz*sqrt( 2*( alphaT*vBnorm + dMEff )/RFactor )/vBnorm/vBnormxy
          B23 =       vBx*sqrt( 2*( alphaT*vBnorm + dMEff )/RFactor )/vBnormxy
        end if
      end if 

      ! Compute random numbers
      call this%GenerateStandardNormalRandom( rdmx ) 
      call this%GenerateStandardNormalRandom( rdmy ) 
      call this%GenerateStandardNormalRandom( rdmz ) 

      ! Compute displacement times random
      dBx = B11*rdmx + B12*rdmy + B13*rdmz 
      dBy = B21*rdmx + B22*rdmy + B23*rdmz 
      dBz = B31*rdmx + B32*rdmy 


  end subroutine pr_DisplacementRandomDischarge


  subroutine pr_DisplacementRandomDischarge2D( this, x, y, z, alphaL, alphaT, dMEff, dBx, dBy, dBz )
      !----------------------------------------------------------------
      ! Computes the product between displacement matrix and random 
      ! vector
      !
      ! Params:
      !     - x, y, z       : local cell coordinates
      !     - alphaL        : longidutinal dispersivity ( assumed for idDim1 )
      !     - alphaT        : transverse dispersivity   ( assumed for idDim2 )
      !     - dMEff         : effective molecular diffusion (corrected by tortuosity )
      !     - dBx, dBy, dBz : random dispersion displacement, output
      !----------------------------------------------------------------
      ! Specifications
      !----------------------------------------------------------------
      implicit none
      class (TrackSubCellType) :: this
      ! input
      doubleprecision, intent(in)    :: x, y, z
      doubleprecision, intent(in)    :: alphaL, alphaT, dMEff
      ! output
      doubleprecision, intent(out) :: dBx, dBy, dBz
      ! local
      doubleprecision :: vB1, vB2, vBnorm
      doubleprecision :: B11, B12, B21, B22
      doubleprecision :: rdm1, rdm2
      doubleprecision :: RFactor
      doubleprecision, dimension(4) :: v000
      doubleprecision, dimension(4) :: v100
      doubleprecision, dimension(4) :: v010
      doubleprecision, dimension(4) :: v110
      doubleprecision, dimension(4) :: v001
      doubleprecision, dimension(4) :: v101
      doubleprecision, dimension(4) :: v011
      doubleprecision, dimension(4) :: v111
      doubleprecision, dimension(3) :: dB
      integer :: idDim1, idDim2
      !----------------------------------------------------------------

      ! Initialize
      dBx    = 0d0
      dBy    = 0d0
      dBz    = 0d0
      dB(:)  = 0d0
      idDim1 = this%dimensions(1)
      idDim2 = this%dimensions(2)

      ! Local copies of corner velocities
      v000 = this%qCorner000 / this%porosity000 
      v100 = this%qCorner100 / this%porosity100 
      v010 = this%qCorner010 / this%porosity010 
      v110 = this%qCorner110 / this%porosity110 
      v001 = this%qCorner001 / this%porosity001 
      v101 = this%qCorner101 / this%porosity101 
      v011 = this%qCorner011 / this%porosity011 
      v111 = this%qCorner111 / this%porosity111

      ! Extract R
      RFactor = this%SubCellData%Retardation

      ! Trilinear interpolation of velocities and norm
      call this%Trilinear( x, y, z, &
                           v000(idDim1), v100(idDim1), v010(idDim1), v110(idDim1), &
                           v001(idDim1), v101(idDim1), v011(idDim1), v111(idDim1), &
                           vB1 )
      call this%Trilinear( x, y, z, &
                           v000(idDim2), v100(idDim2), v010(idDim2), v110(idDim2), &
                           v001(idDim2), v101(idDim2), v011(idDim2), v111(idDim2), &
                           vB2 )
      vBnorm   = sqrt( vB1**2 + vB2**2 )

      ! Displacement matrix terms
      ! Refs: Fernàndez-Garcia et al. 2005; Salamon et al. 2006
      ! Handles the case of zero vBnorm
      B11 = 0d0 
      B12 = 0d0
      B21 = 0d0
      B22 = 0d0
      if ( vBnorm .gt. 0d0 ) then
        B11 =  vB1*sqrt( 2*( alphaL*vBnorm + dMEff )/RFactor )/vBnorm
        B21 =  vB2*sqrt( 2*( alphaL*vBnorm + dMEff )/RFactor )/vBnorm
        B12 = -vB2*sqrt( 2*( alphaT*vBnorm + dMEff )/RFactor )/vBnorm
        B22 =  vB1*sqrt( 2*( alphaT*vBnorm + dMEff )/RFactor )/vBnorm
      end if 

      ! Compute random numbers
      call this%GenerateStandardNormalRandom( rdm1 ) 
      call this%GenerateStandardNormalRandom( rdm2 )

      ! Compute displacement times random
      dB(idDim1) = B11*rdm1 + B12*rdm2
      dB(idDim2) = B21*rdm1 + B22*rdm2

      dBx = dB(1)
      dBy = dB(2)
      dBz = dB(3)


  end subroutine pr_DisplacementRandomDischarge2D


  subroutine pr_DisplacementRandomDischarge1D( this, x, y, z, alphaL, alphaT, dMEff, dBx, dBy, dBz )
      !----------------------------------------------------------------
      ! Computes the product between displacement matrix and random 
      ! vector
      !
      ! Params:
      !     - x, y, z       : local cell coordinates
      !     - alphaL        : longidutinal dispersivity ( assumed for idDim1 )
      !     - alphaT        : transverse dispersivity   ( not used )
      !     - dMEff         : effective molecular diffusion (corrected by tortuosity )
      !     - dBx, dBy, dBz : random dispersion displacement, output
      !----------------------------------------------------------------
      ! Specifications
      !----------------------------------------------------------------
      implicit none
      class (TrackSubCellType) :: this
      ! input
      doubleprecision, intent(in)    :: x, y, z
      doubleprecision, intent(in)    :: alphaL, alphaT, dMEff
      ! output
      doubleprecision, intent(out) :: dBx, dBy, dBz
      ! local
      doubleprecision :: vB1, vBnorm
      doubleprecision :: B11
      doubleprecision :: rdm1
      doubleprecision :: RFactor
      doubleprecision, dimension(4) :: v000
      doubleprecision, dimension(4) :: v100
      doubleprecision, dimension(4) :: v010
      doubleprecision, dimension(4) :: v110
      doubleprecision, dimension(4) :: v001
      doubleprecision, dimension(4) :: v101
      doubleprecision, dimension(4) :: v011
      doubleprecision, dimension(4) :: v111
      doubleprecision, dimension(3) :: dB
      integer :: idDim1
      !----------------------------------------------------------------

      ! Initialize
      dBx    = 0d0
      dBy    = 0d0
      dBz    = 0d0
      dB(:)  = 0d0
      idDim1 = this%dimensions(1)

      ! Local copies of corner velocities
      v000 = this%qCorner000 / this%porosity000 
      v100 = this%qCorner100 / this%porosity100 
      v010 = this%qCorner010 / this%porosity010 
      v110 = this%qCorner110 / this%porosity110 
      v001 = this%qCorner001 / this%porosity001 
      v101 = this%qCorner101 / this%porosity101 
      v011 = this%qCorner011 / this%porosity011 
      v111 = this%qCorner111 / this%porosity111

      ! Extract R
      RFactor = this%SubCellData%Retardation

      ! Trilinear interpolation of velocities and norm
      call this%Trilinear( x, y, z, &
                           v000(idDim1), v100(idDim1), v010(idDim1), v110(idDim1), &
                           v001(idDim1), v101(idDim1), v011(idDim1), v111(idDim1), &
                           vB1 )
      vBnorm = sqrt(vB1**2)

      ! Displacement matrix terms
      ! Refs: Fernàndez-Garcia et al. 2005; Salamon et al. 2006
      ! Handles the case of zero vBnorm
      B11 = 0d0 
      if ( vBnorm .gt. 0d0 ) B11 = sqrt( 2*( alphaL*vBnorm + dMEff )/RFactor )

      ! Compute random numbers
      call this%GenerateStandardNormalRandom( rdm1 ) 

      ! Compute displacement times random
      dB(idDim1) = B11*rdm1

      dBx = dB(1)
      dBy = dB(2)
      dBz = dB(3)


  end subroutine pr_DisplacementRandomDischarge1D


  ! Axisymmetric dispersion model
  subroutine pr_CornerDispersionAxisymmetric( this, &
                                             qprod000, & 
                                             qprod100, &
                                             qprod010, &
                                             qprod110, &
                                             qprod001, &
                                             qprod101, &
                                             qprod011, &
                                             qprod111  )
  !-------------------------------------------------------------------
  ! Compute dispersion at the cell corners, using the 
  ! axisymmetric dispersion model in Lichtner et al. 2002 
  !
  ! Interface:
  !   CornerDispersion
  !-------------------------------------------------------------------
  ! Specifications
  !-------------------------------------------------------------------
  implicit none
  class (TrackSubCellType) :: this
  ! input
  doubleprecision, dimension(8), intent(in) :: qprod000
  doubleprecision, dimension(8), intent(in) :: qprod100
  doubleprecision, dimension(8), intent(in) :: qprod010
  doubleprecision, dimension(8), intent(in) :: qprod110
  doubleprecision, dimension(8), intent(in) :: qprod001
  doubleprecision, dimension(8), intent(in) :: qprod101
  doubleprecision, dimension(8), intent(in) :: qprod011
  doubleprecision, dimension(8), intent(in) :: qprod111
  ! local
  doubleprecision :: alphaL, alphaT, dMEff
  doubleprecision :: alphaLH, alphaLV, alphaTH, alphaTV
  doubleprecision :: aLMinusaT
  doubleprecision :: cosinesq

  doubleprecision :: pDmeff000
  doubleprecision :: pDmeff100
  doubleprecision :: pDmeff010
  doubleprecision :: pDmeff110
  doubleprecision :: pDmeff001
  doubleprecision :: pDmeff101
  doubleprecision :: pDmeff011
  doubleprecision :: pDmeff111
  !-------------------------------------------------------------------

    ! This is ok for uniform dispersion parameters
    dMEff   = this%SubCellData%dMEff
    alphaLH = this%SubCellData%alphaLH
    alphaLV = this%SubCellData%alphaLV
    alphaTH = this%SubCellData%alphaTH
    alphaTV = this%SubCellData%alphaTV
  
    ! Product between porosities and effecfive molecular diffusion 
    pDMeff000 = this%porosity000*dMEff 
    pDMeff100 = this%porosity100*dMEff
    pDMeff010 = this%porosity010*dMEff
    pDMeff110 = this%porosity110*dMEff
    pDMeff001 = this%porosity001*dMEff
    pDMeff101 = this%porosity101*dMEff
    pDMeff011 = this%porosity011*dMEff
    pDMeff111 = this%porosity111*dMEff

    ! 000
    cosinesq = 0d0
    if ( this%qCorner000(4) .gt. 0d0 ) cosinesq = this%qCorner000(3)**2d0/this%qCorner000(4)
    alphaL       = alphaLH + cosinesq*(alphaLV-alphaLH)
    alphaT       = alphaTV + cosinesq*(alphaTH-alphaTV)
    aLMinusaT    = alphaL - alphaT
    this%Dxxx(1) = pDMeff000 + ( alphaL + alphaT*qprod000(8) )*qprod000(1) + alphaTH*qprod000(4)*( 1d0 + qprod000(8) )
    this%Dyyy(1) = pDMeff000 + alphaTH*qprod000(1)*( 1d0 + qprod000(8) ) + qprod000(4)*( alphaL + alphaT*qprod000(8) )
    this%Dxyx(1) = ( alphaL - alphaTH*( 1d0 + qprod000(8) ) + alphaT*qprod000(8) )*qprod000(1)
    this%Dxyy(1) = this%Dxyx(1)
    this%Dzzz(1) = pDMeff000 + alphaT*qprod000(7) + alphaL*qprod000(6)
    this%Dxzx(1) = ( aLMinusaT )*qprod000(3)
    this%Dyzy(1) = ( aLMinusaT )*qprod000(5)
    this%Dxzz(1) = this%Dxzx(1)
    this%Dyzz(1) = this%Dyzy(1)
    ! 100
    cosinesq = 0d0
    if ( this%qCorner100(4) .gt. 0d0 ) cosinesq = this%qCorner100(3)**2d0/this%qCorner100(4)
    alphaL       = alphaLH + cosinesq*(alphaLV-alphaLH)
    alphaT       = alphaTV + cosinesq*(alphaTH-alphaTV)
    aLMinusaT    = alphaL - alphaT
    this%Dxxx(2) = pDMeff100 + ( alphaL + alphaT*qprod100(8) )*qprod100(1) + alphaTH*qprod100(4)*( 1d0 + qprod100(8) )
    this%Dyyy(2) = pDMeff100 + alphaTH*qprod100(1)*( 1d0 + qprod100(8) ) + qprod100(4)*( alphaL + alphaT*qprod100(8) )
    this%Dxyx(2) = ( alphaL - alphaTH*( 1d0 + qprod100(8) ) + alphaT*qprod100(8) )*qprod100(1)
    this%Dxyy(2) = this%Dxyx(2)
    this%Dzzz(2) = pDMeff100 + alphaT*qprod100(7) + alphaL*qprod100(6)
    this%Dxzx(2) = ( aLMinusaT )*qprod100(3)
    this%Dyzy(2) = ( aLMinusaT )*qprod100(5)
    this%Dxzz(2) = this%Dxzx(2)
    this%Dyzz(2) = this%Dyzy(2)
    ! 010
    cosinesq = 0d0
    if ( this%qCorner010(4) .gt. 0d0 ) cosinesq = this%qCorner010(3)**2d0/this%qCorner010(4)
    alphaL       = alphaLH + cosinesq*(alphaLV-alphaLH)
    alphaT       = alphaTV + cosinesq*(alphaTH-alphaTV)
    aLMinusaT    = alphaL - alphaT
    this%Dxxx(3) = pDMeff010 + ( alphaL + alphaT*qprod010(8) )*qprod010(1) + alphaTH*qprod010(4)*( 1d0 + qprod010(8) )
    this%Dyyy(3) = pDMeff010 + alphaTH*qprod010(1)*( 1d0 + qprod010(8) ) + qprod010(4)*( alphaL + alphaT*qprod010(8) )
    this%Dxyx(3) = ( alphaL - alphaTH*( 1d0 + qprod010(8) ) + alphaT*qprod010(8) )*qprod010(1)
    this%Dxyy(3) = this%Dxyx(3)
    this%Dzzz(3) = pDMeff010 + alphaT*qprod010(7) + alphaL*qprod010(6)
    this%Dxzx(3) = ( aLMinusaT )*qprod010(3)
    this%Dyzy(3) = ( aLMinusaT )*qprod010(5)
    this%Dxzz(3) = this%Dxzx(3)
    this%Dyzz(3) = this%Dyzy(3)
    ! 110
    cosinesq = 0d0
    if ( this%qCorner110(4) .gt. 0d0 ) cosinesq = this%qCorner110(3)**2d0/this%qCorner110(4)
    alphaL       = alphaLH + cosinesq*(alphaLV-alphaLH)
    alphaT       = alphaTV + cosinesq*(alphaTH-alphaTV)
    aLMinusaT    = alphaL - alphaT
    this%Dxxx(4) = pDMeff110 + ( alphaL + alphaT*qprod110(8) )*qprod110(1) + alphaTH*qprod110(4)*( 1d0 + qprod110(8) )
    this%Dyyy(4) = pDMeff110 + alphaTH*qprod110(1)*( 1d0 + qprod110(8) ) + qprod110(4)*( alphaL + alphaT*qprod110(8) )
    this%Dxyx(4) = ( alphaL - alphaTH*( 1d0 + qprod110(8) ) + alphaT*qprod110(8) )*qprod110(1)
    this%Dxyy(4) = this%Dxyx(4)
    this%Dzzz(4) = pDMeff110 + alphaT*qprod110(7) + alphaL*qprod110(6)
    this%Dxzx(4) = ( aLMinusaT )*qprod110(3)
    this%Dyzy(4) = ( aLMinusaT )*qprod110(5)
    this%Dxzz(4) = this%Dxzx(4)
    this%Dyzz(4) = this%Dyzy(4)
    ! 001
    cosinesq = 0d0
    if ( this%qCorner001(4) .gt. 0d0 ) cosinesq = this%qCorner001(3)**2d0/this%qCorner001(4)
    alphaL       = alphaLH + cosinesq*(alphaLV-alphaLH)
    alphaT       = alphaTV + cosinesq*(alphaTH-alphaTV)
    aLMinusaT    = alphaL - alphaT
    this%Dxxx(5) = pDMeff001 + ( alphaL + alphaT*qprod001(8) )*qprod001(1) + alphaTH*qprod001(4)*( 1d0 + qprod001(8) )
    this%Dyyy(5) = pDMeff001 + alphaTH*qprod001(1)*( 1d0 + qprod001(8) ) + qprod001(4)*( alphaL + alphaT*qprod001(8) )
    this%Dxyx(5) = ( alphaL - alphaTH*( 1d0 + qprod001(8) ) + alphaT*qprod001(8) )*qprod001(1)
    this%Dxyy(5) = this%Dxyx(5)
    this%Dzzz(5) = pDMeff001 + alphaT*qprod001(7) + alphaL*qprod001(6)
    this%Dxzx(5) = ( aLMinusaT )*qprod001(3)
    this%Dyzy(5) = ( aLMinusaT )*qprod001(5)
    this%Dxzz(5) = this%Dxzx(5)
    this%Dyzz(5) = this%Dyzy(5)
    ! 101
    cosinesq = 0d0
    if ( this%qCorner101(4) .gt. 0d0 ) cosinesq = this%qCorner101(3)**2d0/this%qCorner101(4)
    alphaL       = alphaLH + cosinesq*(alphaLV-alphaLH)
    alphaT       = alphaTV + cosinesq*(alphaTH-alphaTV)
    aLMinusaT    = alphaL - alphaT
    this%Dxxx(6) = pDMeff101 + ( alphaL + alphaT*qprod101(8) )*qprod101(1) + alphaTH*qprod101(4)*( 1d0 + qprod101(8) )
    this%Dyyy(6) = pDMeff101 + alphaTH*qprod101(1)*( 1d0 + qprod101(8) ) + qprod101(4)*( alphaL + alphaT*qprod101(8) )
    this%Dxyx(6) = ( alphaL - alphaTH*( 1d0 + qprod101(8) ) + alphaT*qprod101(8) )*qprod101(1)
    this%Dxyy(6) = this%Dxyx(6)
    this%Dzzz(6) = pDMeff101 + alphaT*qprod101(7) + alphaL*qprod101(6)
    this%Dxzx(6) = ( aLMinusaT )*qprod101(3)
    this%Dyzy(6) = ( aLMinusaT )*qprod101(5)
    this%Dxzz(6) = this%Dxzx(6)
    this%Dyzz(6) = this%Dyzy(6)
    ! 011
    cosinesq = 0d0
    if ( this%qCorner011(4) .gt. 0d0 ) cosinesq = this%qCorner011(3)**2d0/this%qCorner011(4)
    alphaL       = alphaLH + cosinesq*(alphaLV-alphaLH)
    alphaT       = alphaTV + cosinesq*(alphaTH-alphaTV)
    aLMinusaT    = alphaL - alphaT
    this%Dxxx(7) = pDMeff011 + ( alphaL + alphaT*qprod011(8) )*qprod011(1) + alphaTH*qprod011(4)*( 1d0 + qprod011(8) )
    this%Dyyy(7) = pDMeff011 + alphaTH*qprod011(1)*( 1d0 + qprod011(8) ) + qprod011(4)*( alphaL + alphaT*qprod011(8) )
    this%Dxyx(7) = ( alphaL - alphaTH*( 1d0 + qprod011(8) ) + alphaT*qprod011(8) )*qprod011(1)
    this%Dxyy(7) = this%Dxyx(7)
    this%Dzzz(7) = pDMeff011 + alphaT*qprod011(7) + alphaL*qprod011(6)
    this%Dxzx(7) = ( aLMinusaT )*qprod011(3)
    this%Dyzy(7) = ( aLMinusaT )*qprod011(5)
    this%Dxzz(7) = this%Dxzx(7)
    this%Dyzz(7) = this%Dyzy(7)
    ! 111
    cosinesq = 0d0
    if ( this%qCorner111(4) .gt. 0d0 ) cosinesq = this%qCorner111(3)**2d0/this%qCorner111(4)
    alphaL       = alphaLH + cosinesq*(alphaLV-alphaLH)
    alphaT       = alphaTV + cosinesq*(alphaTH-alphaTV)
    aLMinusaT    = alphaL - alphaT
    this%Dxxx(8) = pDMeff111 + ( alphaL + alphaT*qprod111(8) )*qprod111(1) + alphaTH*qprod111(4)*( 1d0 + qprod111(8) )
    this%Dyyy(8) = pDMeff111 + alphaTH*qprod111(1)*( 1d0 + qprod111(8) ) + qprod111(4)*( alphaL + alphaT*qprod111(8) )
    this%Dxyx(8) = ( alphaL - alphaTH*( 1d0 + qprod111(8) ) + alphaT*qprod111(8) )*qprod111(1)
    this%Dxyy(8) = this%Dxyx(8)
    this%Dzzz(8) = pDMeff111 + alphaT*qprod111(7) + alphaL*qprod111(6)
    this%Dxzx(8) = ( aLMinusaT )*qprod111(3)
    this%Dyzy(8) = ( aLMinusaT )*qprod111(5)
    this%Dxzz(8) = this%Dxzx(8)
    this%Dyzz(8) = this%Dyzy(8)

    ! Done
    return

  end subroutine pr_CornerDispersionAxisymmetric


  subroutine pr_DisplacementRandomDischargeAxisymmetric( this, x, y, z, & 
                          alphaL, alphaT, alphaTH, dMEff, dBx, dBy, dBz ) 
      !----------------------------------------------------------------
      ! Computes the product between displacement matrix, computed as 
      ! in Lichtner et al. 2002, and the random vector
      !
      ! Params:
      !     - x, y, z       : local cell coordinates
      !     - alphaL        : longidutinal dispersivity
      !     - alphaT        : transverse dispersivity
      !     - alphaTH       : transverse dispersivity
      !     - dMEff         : effective molecular diffusion (corrected by tortuosity )
      !     - dBx, dBy, dBz : random dispersion displacement, output
      !----------------------------------------------------------------
      ! Specifications
      !----------------------------------------------------------------
      implicit none
      class (TrackSubCellType) :: this
      ! input
      doubleprecision, intent(in)    :: x, y, z
      doubleprecision, intent(in)    :: alphaL, alphaT, alphaTH, dMEff
      ! output
      doubleprecision, intent(out) :: dBx, dBy, dBz
      ! local
      doubleprecision :: vBx, vBy, vBz, vBnorm, vBnormxy
      doubleprecision :: B11, B12, B13, B21, B22, B23, B31, B32
      doubleprecision :: rdmx, rdmy, rdmz
      doubleprecision :: RFactor
      doubleprecision, dimension(4) :: v000
      doubleprecision, dimension(4) :: v100
      doubleprecision, dimension(4) :: v010
      doubleprecision, dimension(4) :: v110
      doubleprecision, dimension(4) :: v001
      doubleprecision, dimension(4) :: v101
      doubleprecision, dimension(4) :: v011
      doubleprecision, dimension(4) :: v111
      !----------------------------------------------------------------

      ! Initialize
      dBx = 0d0
      dBy = 0d0
      dBz = 0d0

      ! Local copies of corner velocities
      ! pointers ?
      v000 = this%qCorner000 / this%porosity000 
      v100 = this%qCorner100 / this%porosity100 
      v010 = this%qCorner010 / this%porosity010 
      v110 = this%qCorner110 / this%porosity110 
      v001 = this%qCorner001 / this%porosity001 
      v101 = this%qCorner101 / this%porosity101 
      v011 = this%qCorner011 / this%porosity011 
      v111 = this%qCorner111 / this%porosity111

      ! Extract R
      RFactor = this%SubCellData%Retardation

      ! Trilinear interpolation of velocities and norm
      call this%Trilinear( x, y, z, &
                           v000(1), v100(1), v010(1), v110(1), &
                           v001(1), v101(1), v011(1), v111(1), &
                           vBx )
      call this%Trilinear( x, y, z, &
                           v000(2), v100(2), v010(2), v110(2), &
                           v001(2), v101(2), v011(2), v111(2), &
                           vBy )
      call this%Trilinear( x, y, z, &
                           v000(3), v100(3), v010(3), v110(3), &
                           v001(3), v101(3), v011(3), v111(3), &
                           vBz )
      vBnorm   = sqrt( vBx**2 + vBy**2 + vBz**2 )
      vBnormxy = sqrt( vBx**2 + vBy**2 )
   

      ! Displacement matrix terms
      ! Refs: Fernàndez-Garcia et al. 2005; Salamon et al. 2006
      ! Handles the case of zero vBnorm
      B11 = 0d0 
      B12 = 0d0
      B13 = 0d0
      B21 = 0d0
      B22 = 0d0
      B23 = 0d0
      B31 = 0d0
      B32 = 0d0
      if ( vBnorm .gt. 0d0 ) then
        B11 =       vBx*sqrt( 2*( alphaL*vBnorm + dMEff )/RFactor )/vBnorm
        B21 =       vBy*sqrt( 2*( alphaL*vBnorm + dMEff )/RFactor )/vBnorm
        B31 =       vBz*sqrt( 2*( alphaL*vBnorm + dMEff )/RFactor )/vBnorm
        B32 =  vBnormxy*sqrt( 2*( alphaT*vBnorm + dMEff )/RFactor )/vBnorm
        if ( vBnormxy .gt. 0d0 ) then
          B12 =  -vBx*vBz*sqrt( 2*( alphaT*vBnorm + dMEff )/RFactor )/vBnorm/vBnormxy
          B13 =      -vBy*sqrt( 2*( alphaTH*vBnorm + dMEff )/RFactor )/vBnormxy
          B22 =  -vBy*vBz*sqrt( 2*( alphaT*vBnorm + dMEff )/RFactor )/vBnorm/vBnormxy
          B23 =       vBx*sqrt( 2*( alphaTH*vBnorm + dMEff )/RFactor )/vBnormxy
        end if
      end if 

      ! Compute random numbers
      call this%GenerateStandardNormalRandom( rdmx ) 
      call this%GenerateStandardNormalRandom( rdmy ) 
      call this%GenerateStandardNormalRandom( rdmz ) 

      ! Compute displacement times random
      dBx = B11*rdmx + B12*rdmy + B13*rdmz 
      dBy = B21*rdmx + B22*rdmy + B23*rdmz 
      dBz = B31*rdmx + B32*rdmy 


  end subroutine pr_DisplacementRandomDischargeAxisymmetric


  subroutine pr_DisplacementRandomDischargeAxisymmetric2D( this, x, y, z, & 
                            alphaL, alphaT, alphaTH, dMEff, dBx, dBy, dBz )
      !----------------------------------------------------------------
      ! Computes the product between displacement matrix and random 
      ! vector
      !
      ! Params:
      !     - x, y, z       : local cell coordinates
      !     - alphaL        : longidutinal dispersivity ( assumed for idDim1 )
      !     - alphaT        : transverse dispersivity  (not used) 
      !     - alphaTH       : horizontal transverse dispersivity   
      !     - dMEff         : effective molecular diffusion (corrected by tortuosity )
      !     - dBx, dBy, dBz : random dispersion displacement, output
      !----------------------------------------------------------------
      ! Specifications
      !----------------------------------------------------------------
      implicit none
      class (TrackSubCellType) :: this
      ! input
      doubleprecision, intent(in)    :: x, y, z
      doubleprecision, intent(in)    :: alphaL, alphaT, alphaTH, dMEff
      ! output
      doubleprecision, intent(out) :: dBx, dBy, dBz
      ! local
      doubleprecision :: vB1, vB2, vBnorm
      doubleprecision :: B11, B12, B21, B22
      doubleprecision :: rdm1, rdm2
      doubleprecision :: RFactor
      doubleprecision, dimension(4) :: v000
      doubleprecision, dimension(4) :: v100
      doubleprecision, dimension(4) :: v010
      doubleprecision, dimension(4) :: v110
      doubleprecision, dimension(4) :: v001
      doubleprecision, dimension(4) :: v101
      doubleprecision, dimension(4) :: v011
      doubleprecision, dimension(4) :: v111
      doubleprecision, dimension(3) :: dB
      integer :: idDim1, idDim2
      !----------------------------------------------------------------

      ! Initialize
      dBx    = 0d0
      dBy    = 0d0
      dBz    = 0d0
      dB(:)  = 0d0
      idDim1 = this%dimensions(1)
      idDim2 = this%dimensions(2)

      ! Local copies of corner velocities
      v000 = this%qCorner000 / this%porosity000 
      v100 = this%qCorner100 / this%porosity100 
      v010 = this%qCorner010 / this%porosity010 
      v110 = this%qCorner110 / this%porosity110 
      v001 = this%qCorner001 / this%porosity001 
      v101 = this%qCorner101 / this%porosity101 
      v011 = this%qCorner011 / this%porosity011 
      v111 = this%qCorner111 / this%porosity111

      ! Extract R
      RFactor = this%SubCellData%Retardation

      ! Trilinear interpolation of velocities and norm
      call this%Trilinear( x, y, z, &
                           v000(idDim1), v100(idDim1), v010(idDim1), v110(idDim1), &
                           v001(idDim1), v101(idDim1), v011(idDim1), v111(idDim1), &
                           vB1 )
      call this%Trilinear( x, y, z, &
                           v000(idDim2), v100(idDim2), v010(idDim2), v110(idDim2), &
                           v001(idDim2), v101(idDim2), v011(idDim2), v111(idDim2), &
                           vB2 )
      vBnorm   = sqrt( vB1**2 + vB2**2 )

      ! Displacement matrix terms
      ! Refs: Fernàndez-Garcia et al. 2005; Salamon et al. 2006
      ! Handles the case of zero vBnorm
      B11 = 0d0 
      B12 = 0d0
      B21 = 0d0
      B22 = 0d0
      if ( vBnorm .gt. 0d0 ) then
        B11 =  vB1*sqrt( 2*( alphaL*vBnorm + dMEff )/RFactor )/vBnorm
        B21 =  vB2*sqrt( 2*( alphaL*vBnorm + dMEff )/RFactor )/vBnorm
        B12 = -vB2*sqrt( 2*( alphaTH*vBnorm + dMEff )/RFactor )/vBnorm
        B22 =  vB1*sqrt( 2*( alphaTH*vBnorm + dMEff )/RFactor )/vBnorm
      end if 

      ! Compute random numbers
      call this%GenerateStandardNormalRandom( rdm1 ) 
      call this%GenerateStandardNormalRandom( rdm2 )

      ! Compute displacement times random
      dB(idDim1) = B11*rdm1 + B12*rdm2
      dB(idDim2) = B21*rdm1 + B22*rdm2

      dBx = dB(1)
      dBy = dB(2)
      dBz = dB(3)


  end subroutine pr_DisplacementRandomDischargeAxisymmetric2D


  subroutine pr_DisplacementRandomDischargeAxisymmetric1D( this, x, y, z, &
                alphaL, alphaT, alphaTH, dMEff, dBx, dBy, dBz )
      !----------------------------------------------------------------
      ! Computes the product between displacement matrix and random 
      ! vector
      !
      ! Params:
      !     - x, y, z       : local cell coordinates
      !     - alphaL        : longidutinal dispersivity ( assumed for idDim1 )
      !     - alphaT        : transverse dispersivity   ( not used )
      !     - alphaTH       : horizontal transverse dispersivity   ( not used )
      !     - dMEff         : effective molecular diffusion (corrected by tortuosity )
      !     - dBx, dBy, dBz : random dispersion displacement, output
      !----------------------------------------------------------------
      ! Specifications
      !----------------------------------------------------------------
      implicit none
      class (TrackSubCellType) :: this
      ! input
      doubleprecision, intent(in)    :: x, y, z
      doubleprecision, intent(in)    :: alphaL, alphaT, alphaTH, dMEff
      ! output
      doubleprecision, intent(out) :: dBx, dBy, dBz
      ! local
      doubleprecision :: vB1, vBnorm
      doubleprecision :: B11
      doubleprecision :: rdm1
      doubleprecision :: RFactor
      doubleprecision, dimension(4) :: v000
      doubleprecision, dimension(4) :: v100
      doubleprecision, dimension(4) :: v010
      doubleprecision, dimension(4) :: v110
      doubleprecision, dimension(4) :: v001
      doubleprecision, dimension(4) :: v101
      doubleprecision, dimension(4) :: v011
      doubleprecision, dimension(4) :: v111
      doubleprecision, dimension(3) :: dB
      integer :: idDim1
      !----------------------------------------------------------------

      ! Initialize
      dBx    = 0d0
      dBy    = 0d0
      dBz    = 0d0
      dB(:)  = 0d0
      idDim1 = this%dimensions(1)

      ! Local copies of corner velocities
      v000 = this%qCorner000 / this%porosity000 
      v100 = this%qCorner100 / this%porosity100 
      v010 = this%qCorner010 / this%porosity010 
      v110 = this%qCorner110 / this%porosity110 
      v001 = this%qCorner001 / this%porosity001 
      v101 = this%qCorner101 / this%porosity101 
      v011 = this%qCorner011 / this%porosity011 
      v111 = this%qCorner111 / this%porosity111

      ! Extract R
      RFactor = this%SubCellData%Retardation

      ! Trilinear interpolation of velocities and norm
      call this%Trilinear( x, y, z, &
                           v000(idDim1), v100(idDim1), v010(idDim1), v110(idDim1), &
                           v001(idDim1), v101(idDim1), v011(idDim1), v111(idDim1), &
                           vB1 )
      vBnorm = sqrt(vB1**2)

      ! Displacement matrix terms
      ! Refs: Fernàndez-Garcia et al. 2005; Salamon et al. 2006
      ! Handles the case of zero vBnorm
      B11 = 0d0 
      if ( vBnorm .gt. 0d0 ) B11 = sqrt( 2*( alphaL*vBnorm + dMEff )/RFactor )

      ! Compute random numbers
      call this%GenerateStandardNormalRandom( rdm1 ) 

      ! Compute displacement times random
      dB(idDim1) = B11*rdm1

      dBx = dB(1)
      dBy = dB(2)
      dBz = dB(3)


  end subroutine pr_DisplacementRandomDischargeAxisymmetric1D


  ! RWPT
  subroutine pr_GenerateStandardNormalRandom( this, random_value )
    !----------------------------------------------------------------
    ! Generate a random number from an standard normal distribution
    ! Inherited from RW3D:library_gslib:random_normal
    !----------------------------------------------------------------
    ! Specifications
    !----------------------------------------------------------------
    implicit none
    class(TrackSubCellType) :: this 
    ! input/output
    doubleprecision, intent(out) :: random_value
    ! local
    doubleprecision :: harvest(12)
    !----------------------------------------------------------------

    call random_number (harvest)
    random_value = sum(harvest)-6.d0

    return

  end subroutine  pr_GenerateStandardNormalRandom


  subroutine pr_GenerateStandardNormalRandomZig( this, random_value )
    !----------------------------------------------------------------
    ! Generate a random number from an standard normal distribution
    ! It uses rng_par_zig library  
    !----------------------------------------------------------------
    ! Specifications
    !----------------------------------------------------------------
    implicit none
    class(TrackSubCellType) :: this 
    ! input/output
    doubleprecision, intent(out) :: random_value
    ! local
    integer, parameter :: nvals = 12
    integer :: m, threadId
    doubleprecision :: r,p 
    !----------------------------------------------------------------

#ifdef _OPENMP
    threadId = omp_get_thread_num()
#else
    threadId = 0
#endif
    
    r = 0d0
    do m=1,nvals
      call rng_uni(p,threadId)
      r = r + p
    end do
    random_value = r - 6.d0
   
    return

  end subroutine  pr_GenerateStandardNormalRandomZig


  ! RWPT
  subroutine pr_ComputeCornerVariables( this, currentCellData, neighborCellData )
  !----------------------------------------------------------------
  ! Function called in TrackCell
  !----------------------------------------------------------------
  ! Specifications
  !----------------------------------------------------------------
  implicit none
  class(TrackSubCellType) :: this
  type(ModpathCellDataType) :: currentCellData
  type(ModpathCellDataType), dimension(2,18) :: neighborCellData
  integer, dimension(3,18)         :: neighborSubCellIndexes   ! nbcell, subRow, subColumn
  doubleprecision, dimension(6,18) :: neighborSubCellFaceFlows ! nbcell, flowFaceNumber
  doubleprecision, dimension(3,18) :: neighborSubCellFaceAreas ! nbcell, faceDirection
  doubleprecision, dimension(18)   :: neighborSubCellVolume   
  doubleprecision, dimension(18)   :: neighborSubCellPorosity

  ! Calculate discharge products necessary for divergence
  ! of local dispersion
  ! xx, xy, xz, yy, yz, zz, (xx+yy), zz/(xx+yy)
  doubleprecision, dimension(8) :: qprod000
  doubleprecision, dimension(8) :: qprod100
  doubleprecision, dimension(8) :: qprod010
  doubleprecision, dimension(8) :: qprod110
  doubleprecision, dimension(8) :: qprod001
  doubleprecision, dimension(8) :: qprod101
  doubleprecision, dimension(8) :: qprod011
  doubleprecision, dimension(8) :: qprod111
  !doubleprecision :: alphaL, alphaT, dMEff
  !doubleprecision :: aLMinusaT
  !doubleprecision :: aTQPlusPDMEff
  !----------------------------------------------------------------
  
    ! Get sub cell indexes for current sub cell location
    neighborSubCellIndexes = currentCellData%GetNeighborSubCellIndexes( &
                           this%SubCellData%Row, this%SubCellData%Column )

    ! Fill neighbor cells faceFlows
    call pr_FillNeighborSubCellVariables( this, currentCellData, neighborCellData, &
       neighborSubCellIndexes, neighborSubCellFaceFlows, neighborSubCellFaceAreas, &
                                    neighborSubCellVolume, neighborSubCellPorosity )
    ! Compute discharge
    call pr_ComputeCornerDischarge( this, currentCellData, & 
        neighborSubCellFaceFlows, neighborSubCellFaceAreas )

    ! Compute porosities
    call this%ComputeCornerPorosity( neighborSubCellVolume, neighborSubCellPorosity )

    ! Calculate discharge products necessary for divergence
    ! of local dispersion
    ! xx, xy, xz, yy, yz, zz, (xx+yy), zz/(xx+yy)
    ! Note: entries 1:7 are normalized by qnorm
    ! Note: entries 7,8 are only needed for the axisymmetric model
    ! Note: these are needed for the isotropic and axisymmetric disp model
    ! 000 
    qprod000(:) = 0d0
    if( this%qCorner000(4) .gt. 0d0 ) then
      qprod000(1) = this%qCorner000(1)**2d0
      qprod000(2) = this%qCorner000(1)*this%qCorner000(2)
      qprod000(3) = this%qCorner000(1)*this%qCorner000(3)
      qprod000(4) = this%qCorner000(2)**2d0
      qprod000(5) = this%qCorner000(2)*this%qCorner000(3)
      qprod000(6) = this%qCorner000(3)**2d0
      qprod000(7) = qprod000(1) + qprod000(4)
      if ( qprod000(7).gt.0d0 ) qprod000(8) = qprod000(6)/qprod000(7)
      qprod000(1:7) = qprod000(1:7)/this%qCorner000(4)
    end if 
    ! 100 
    qprod100(:) = 0d0
    if( this%qCorner100(4) .gt. 0d0 ) then
      qprod100(1) = this%qCorner100(1)**2d0
      qprod100(2) = this%qCorner100(1)*this%qCorner100(2)
      qprod100(3) = this%qCorner100(1)*this%qCorner100(3)
      qprod100(4) = this%qCorner100(2)**2d0
      qprod100(5) = this%qCorner100(2)*this%qCorner100(3)
      qprod100(6) = this%qCorner100(3)**2d0
      qprod100(7) = qprod100(1) + qprod100(4)
      if ( qprod100(7).gt.0d0 ) qprod100(8) = qprod100(6)/qprod100(7)
      qprod100(1:7) = qprod100(1:7)/this%qCorner100(4)
    end if 
    ! 010 
    qprod010(:) = 0d0
    if( this%qCorner010(4) .gt. 0d0 ) then
      qprod010(1) = this%qCorner010(1)**2d0
      qprod010(2) = this%qCorner010(1)*this%qCorner010(2)
      qprod010(3) = this%qCorner010(1)*this%qCorner010(3)
      qprod010(4) = this%qCorner010(2)**2d0
      qprod010(5) = this%qCorner010(2)*this%qCorner010(3)
      qprod010(6) = this%qCorner010(3)**2d0
      qprod010(7) = qprod010(1) + qprod010(4)
      if ( qprod010(7).gt.0d0 ) qprod010(8) = qprod010(6)/qprod010(7)
      qprod010(1:7) = qprod010(1:7)/this%qCorner010(4)
    end if 
    ! 110 
    qprod110(:) = 0d0
    if( this%qCorner110(4) .gt. 0d0 ) then
      qprod110(1) = this%qCorner110(1)**2d0
      qprod110(2) = this%qCorner110(1)*this%qCorner110(2)
      qprod110(3) = this%qCorner110(1)*this%qCorner110(3)
      qprod110(4) = this%qCorner110(2)**2d0
      qprod110(5) = this%qCorner110(2)*this%qCorner110(3)
      qprod110(6) = this%qCorner110(3)**2d0
      qprod110(7) = qprod110(1) + qprod110(4)
      if ( qprod110(7).gt.0d0 ) qprod110(8) = qprod110(6)/qprod110(7)
      qprod110(1:7) = qprod110(1:7)/this%qCorner110(4)
    end if 
    ! 001
    qprod001(:) = 0d0
    if( this%qCorner001(4) .gt. 0d0 ) then
      qprod001(1) = this%qCorner001(1)**2d0
      qprod001(2) = this%qCorner001(1)*this%qCorner001(2)
      qprod001(3) = this%qCorner001(1)*this%qCorner001(3)
      qprod001(4) = this%qCorner001(2)**2d0
      qprod001(5) = this%qCorner001(2)*this%qCorner001(3)
      qprod001(6) = this%qCorner001(3)**2d0
      qprod001(7) = qprod001(1) + qprod001(4)
      if ( qprod001(7).gt.0d0 ) qprod001(8) = qprod001(6)/qprod001(7)
      qprod001(1:7) = qprod001(1:7)/this%qCorner001(4)
    end if 
    ! 101
    qprod101(:) = 0d0
    if( this%qCorner101(4) .gt. 0d0 ) then
      qprod101(1) = this%qCorner101(1)**2d0
      qprod101(2) = this%qCorner101(1)*this%qCorner101(2)
      qprod101(3) = this%qCorner101(1)*this%qCorner101(3)
      qprod101(4) = this%qCorner101(2)**2d0
      qprod101(5) = this%qCorner101(2)*this%qCorner101(3)
      qprod101(6) = this%qCorner101(3)**2d0
      qprod101(7) = qprod101(1) + qprod101(4)
      if ( qprod101(7).gt.0d0 ) qprod101(8) = qprod101(6)/qprod101(7)
      qprod101(1:7) = qprod101(1:7)/this%qCorner101(4)
    end if 
    ! 011
    qprod011(:) = 0d0
    if( this%qCorner011(4) .gt. 0d0 ) then
      qprod011(1) = this%qCorner011(1)**2d0
      qprod011(2) = this%qCorner011(1)*this%qCorner011(2)
      qprod011(3) = this%qCorner011(1)*this%qCorner011(3)
      qprod011(4) = this%qCorner011(2)**2d0
      qprod011(5) = this%qCorner011(2)*this%qCorner011(3)
      qprod011(6) = this%qCorner011(3)**2d0
      qprod011(7) = qprod011(1) + qprod011(4)
      if ( qprod011(7).gt.0d0 ) qprod011(8) = qprod011(6)/qprod011(7)
      qprod011(1:7) = qprod011(1:7)/this%qCorner011(4)
    end if 
    ! 111
    qprod111(:) = 0d0
    if( this%qCorner111(4) .gt. 0d0 ) then
      qprod111(1) = this%qCorner111(1)**2d0
      qprod111(2) = this%qCorner111(1)*this%qCorner111(2)
      qprod111(3) = this%qCorner111(1)*this%qCorner111(3)
      qprod111(4) = this%qCorner111(2)**2d0
      qprod111(5) = this%qCorner111(2)*this%qCorner111(3)
      qprod111(6) = this%qCorner111(3)**2d0
      qprod111(7) = qprod111(1) + qprod111(4)
      if ( qprod111(7).gt.0d0 ) qprod111(8) = qprod111(6)/qprod111(7)
      qprod111(1:7) = qprod111(1:7)/this%qCorner111(4)
    end if

    ! Calculates dispersion at the cell corners depending
    ! on the selected dispersion model.
    !   - this%Dxxx
    !   - this%Dyyy
    !   - this%Dzzz
    !   - this%Dxyx
    !   - this%Dxzx
    !   - this%Dxyy
    !   - this%Dyzy
    !   - this%Dxzz
    !   - this%Dyzz
    call this%ComputeCornerDispersion( &
        qprod000, & 
        qprod100, &
        qprod010, &
        qprod110, &
        qprod001, &
        qprod101, &
        qprod011, &
        qprod111  )


    ! Done
    return

  end subroutine pr_ComputeCornerVariables



  ! RWPT
  ! Could be differentiated for usg vs structured to avoid getsubcellcount if check 
  subroutine pr_FillNeighborSubCellVariables( this, centerCellData, neighborCellData, &
          neighborSubCellIndexes, neighborSubCellFaceFlows, neighborSubCellFaceAreas, & 
                                      neighborSubCellVolume, neighborSubCellPorosity  )
    !-----------------------------------------------------------
    !
    !-----------------------------------------------------------
    ! Specifications
    !-----------------------------------------------------------
    implicit none
    class (TrackSubCellType) :: this
    class (ModpathCellDataType), intent(in) :: centerCellData 
    class (ModpathCellDataType), dimension(2,18), intent(in) :: neighborCellData 
    integer, dimension(3,18), intent(in) ::  neighborSubCellIndexes ! nCellBuffer, subRow, subCol
    doubleprecision, dimension(6,18), intent(out) :: neighborSubCellFaceFlows ! faceFlows
    doubleprecision, dimension(3,18), intent(out) :: neighborSubCellFaceAreas ! faceAreas
    doubleprecision, dimension(18)  , intent(out) :: neighborSubCellVolume    ! volumes
    doubleprecision, dimension(18)  , intent(out) :: neighborSubCellPorosity  ! porosities
    logical :: skipSubCells = .false.
    integer :: subConnectionIndex, neighborSubRow, neighborSubColumn
    integer :: n
    !-----------------------------------------------------------


    ! Reset arrays
    neighborSubCellFaceFlows = 0d0
    neighborSubCellFaceAreas = 0d0
    neighborSubCellVolume    = 0d0
    neighborSubCellPorosity  = 0d0


    ! If current cell is not refined
    if ( centerCellData%GetSubCellCount() .eq. 1 ) then 
      ! Cell is not refined then neighbors
      ! have the same size, some maybe refined
      ! due to connections with smaller cells.
      ! Skip sub cell indexation in such case. 
      skipSubCells = .true.
      do n = 1,18

        if ( neighborCellData( 1, neighborSubCellIndexes( 1, n ) )%CellNumber .eq. 0 ) cycle

        call neighborCellData( &
            1, neighborSubCellIndexes( 1, n ) )%FillSubCellFaceFlows( 1, 1, &
                              neighborSubCellFaceFlows( :, n ), skipSubCells )
        call neighborCellData( & 
            1, neighborSubCellIndexes( 1, n ) )%FillSubCellFaceAreas( neighborSubCellFaceAreas( :, n ), skipSubCells ) 

        neighborSubCellVolume(n)   = neighborCellData( 1, neighborSubCellIndexes( 1, n ) )%GetVolume( skipSubCells )
        neighborSubCellPorosity(n) = neighborCellData( 1, neighborSubCellIndexes( 1, n ) )%Porosity

      end do

      ! Done 
      return

    end if  


    ! If current cell is refined
    do n = 1, 18
      skipSubCells = .false.

      ! Fill face flows with one of its own subcells
      if ( neighborSubCellIndexes( 1, n ) .eq. 0 ) then
        call centerCellData%FillSubCellFaceFlowsBuffer( &
                          neighborSubCellIndexes( 2, n ), neighborSubCellIndexes( 3, n ), & 
                                                         neighborSubCellFaceFlows( :, n ) )

        call centerCellData%FillSubCellFaceAreas( neighborSubCellFaceAreas( :, n ), skipSubCells )

        neighborSubCellVolume(n)   = centerCellData%GetVolume( skipSubCells )
        neighborSubCellPorosity(n) = centerCellData%Porosity

        cycle

      end if

      ! Fill face flows with data obtained from neighbors
      if ( neighborCellData( 2, neighborSubCellIndexes( 1, n ) )%CellNumber .gt. 0 ) then
        ! If double buffer 

        ! When the buffer is double, 
        ! these cells are smaller, then skip 
        ! sub cell indexation if refined. 
        if ( & 
          ( neighborCellData(1, neighborSubCellIndexes( 1, n ) )%GetSubCellCount() .gt. 1 ) .or. &
          ( neighborCellData(2, neighborSubCellIndexes( 1, n ) )%GetSubCellCount() .gt. 1 ) ) then 
          skipSubCells = .true.
        end if 

        ! Detect from which direction where requested.
        if ( neighborCellData( 2, neighborSubCellIndexes( 1, n ) )%requestedFromDirection .eq. 1 ) then 
          ! If x-direction, location in buffer is given by subcell row index 
          subConnectionIndex = neighborSubCellIndexes( 2, n )
          neighborSubRow     = 1 
          neighborSubColumn  = 1
        else if ( neighborCellData( 2, neighborSubCellIndexes( 1, n ) )%requestedFromDirection .eq. 2 ) then
          ! If y-direction, location in buffer is given by subcell column index 
          subConnectionIndex = neighborSubCellIndexes( 3, n )
          neighborSubRow     = 1 
          neighborSubColumn  = 1
        else 
          ! If z-direction, location in buffer is also given 
          ! by subcell column index. This is true because double 
          ! buffers requested from a vertical face are possible 
          ! only if horizontal buffer is double too. This means 
          ! that this case is related to indirect buffers, which 
          ! in the current protocol of neighbor cells initialization 
          ! are only possible after initialization of a direct connection 
          ! through the y-direction. 
          subConnectionIndex = neighborSubCellIndexes( 3, n )
          neighborSubRow     = 1 
          neighborSubColumn  = 1
        end if 

      else if ( &
        ( neighborCellData( 1, neighborSubCellIndexes( 1, n ) )%CellNumber .gt. 0 ) .or. & 
        ( neighborCellData( 1, neighborSubCellIndexes( 1, n ) )%fromSubCell ) ) then
        ! If single buffer 

        if ( neighborCellData(1, neighborSubCellIndexes( 1, n ) )%GetSubCellCount() .gt. 1 ) then 
          ! If its refined request the sub cell indexes
          subConnectionIndex = 1
          neighborSubRow     = neighborSubCellIndexes( 2, n )
          neighborSubColumn  = neighborSubCellIndexes( 3, n ) 
        else
          ! If its not refined, all ones
          subConnectionIndex = 1
          neighborSubRow     = 1
          neighborSubColumn  = 1
        end if

      else
        ! Buffer is empty, continue
        cycle
      end if  

      ! Fill variables 
      call neighborCellData( subConnectionIndex, neighborSubCellIndexes( 1, n ) )%FillSubCellFaceFlows( &
                      neighborSubRow, neighborSubColumn, neighborSubCellFaceFlows( :, n ), skipSubCells )

      call neighborCellData( subConnectionIndex, neighborSubCellIndexes( 1, n ) )%FillSubCellFaceAreas( &
                                                         neighborSubCellFaceAreas( :, n ), skipSubCells )  

      neighborSubCellVolume(n)   = neighborCellData(&
        subConnectionIndex, neighborSubCellIndexes( 1, n ) )%GetVolume( skipSubCells )
      neighborSubCellPorosity(n) = neighborCellData(&
        subConnectionIndex, neighborSubCellIndexes( 1, n ) )%Porosity

    end do


    ! Done
    return


  end subroutine pr_FillNeighborSubCellVariables


  ! RWPT
  subroutine pr_ComputeCornerDischarge( this, currentCellData, &
            neighborSubCellFaceFlows, neighborSubCellFaceAreas )
  !----------------------------------------------------------------
  !----------------------------------------------------------------
  ! Specifications
  !----------------------------------------------------------------
  implicit none
  class(TrackSubCellType) :: this
  ! input
  type(ModpathCellDataType) :: currentCellData
  doubleprecision, dimension(6,18), intent(in) :: neighborSubCellFaceFlows ! nbcell, flowFaceNumber
  doubleprecision, dimension(3,18), intent(in) :: neighborSubCellFaceAreas ! nbcell, faceDirection
  ! local
  doubleprecision, dimension(6) :: centerSubCellFaceFlows
  integer :: n
  !----------------------------------------------------------------


    ! Fill faceFlows from current subcell
    call currentCellData%FillSubCellFaceFlowsBuffer( &
        this%SubCellData%Row, this%SubCellData%Column, centerSubCellFaceFlows )

    ! For each dimension
    ! dimension mask ?
    do n = 1, 3
      this%qCorner000(n) = this%GetInterpolatedCornerDischarge( & 
        centerSubCellFaceFlows, neighborSubCellFaceFlows, neighborSubCellFaceAreas, 1, n ) 
      this%qCorner100(n) = this%GetInterpolatedCornerDischarge( & 
        centerSubCellFaceFlows, neighborSubCellFaceFlows, neighborSubCellFaceAreas, 2, n ) 
      this%qCorner010(n) = this%GetInterpolatedCornerDischarge( & 
        centerSubCellFaceFlows, neighborSubCellFaceFlows, neighborSubCellFaceAreas, 3, n ) 
      this%qCorner110(n) = this%GetInterpolatedCornerDischarge( & 
        centerSubCellFaceFlows, neighborSubCellFaceFlows, neighborSubCellFaceAreas, 4, n ) 
      this%qCorner001(n) = this%GetInterpolatedCornerDischarge( & 
        centerSubCellFaceFlows, neighborSubCellFaceFlows, neighborSubCellFaceAreas, 5, n ) 
      this%qCorner101(n) = this%GetInterpolatedCornerDischarge( & 
        centerSubCellFaceFlows, neighborSubCellFaceFlows, neighborSubCellFaceAreas, 6, n ) 
      this%qCorner011(n) = this%GetInterpolatedCornerDischarge( & 
        centerSubCellFaceFlows, neighborSubCellFaceFlows, neighborSubCellFaceAreas, 7, n ) 
      this%qCorner111(n) = this%GetInterpolatedCornerDischarge( & 
        centerSubCellFaceFlows, neighborSubCellFaceFlows, neighborSubCellFaceAreas, 8, n ) 
    end do 

    ! Norm
    this%qCorner000(4) = sqrt( this%qCorner000(1)**2 + this%qCorner000(2)**2 + this%qCorner000(3)**2 )
    this%qCorner100(4) = sqrt( this%qCorner100(1)**2 + this%qCorner100(2)**2 + this%qCorner100(3)**2 )
    this%qCorner010(4) = sqrt( this%qCorner010(1)**2 + this%qCorner010(2)**2 + this%qCorner010(3)**2 )
    this%qCorner110(4) = sqrt( this%qCorner110(1)**2 + this%qCorner110(2)**2 + this%qCorner110(3)**2 )
    this%qCorner001(4) = sqrt( this%qCorner001(1)**2 + this%qCorner001(2)**2 + this%qCorner001(3)**2 )
    this%qCorner101(4) = sqrt( this%qCorner101(1)**2 + this%qCorner101(2)**2 + this%qCorner101(3)**2 )
    this%qCorner011(4) = sqrt( this%qCorner011(1)**2 + this%qCorner011(2)**2 + this%qCorner011(3)**2 )
    this%qCorner111(4) = sqrt( this%qCorner111(1)**2 + this%qCorner111(2)**2 + this%qCorner111(3)**2 )


  end subroutine pr_ComputeCornerDischarge



  ! RWPT
  function pr_GetInterpolatedCornerDischarge( this, centerSubCellFaceFlows, &
                        neighborSubCellFaceFlows, neighborSubCellFaceAreas, & 
                                                cornerIndex, faceDirection  ) result( qCorner )
  !---------------------------------------------------------------------
  ! Compute interpolated corner discharge as the sum of flow rates
  ! of contributing faces divided by the sum of face areas
  !
  ! Initializes sums with contribution of center cell and 
  ! loop over remaining involved cells
  !---------------------------------------------------------------------
  implicit none
  class( TrackSubCellType ) :: this
  ! input
  doubleprecision, dimension(6)   , intent(in) :: centerSubCellFaceFlows 
  doubleprecision, dimension(6,18), intent(in) :: neighborSubCellFaceFlows
  doubleprecision, dimension(3,18), intent(in) :: neighborSubCellFaceAreas
  integer, intent(in) :: cornerIndex, faceDirection
  ! local
  doubleprecision :: sumArea, sumFlowRate = 0d0
  integer :: m 
  ! output 
  doubleprecision :: qCorner    
  !---------------------------------------------------------------------
  
    
    select case( faceDirection ) 
      case (1) 
        sumArea     = this%SubCellData%DY*this%SubCellData%DZ
        sumFlowRate = centerSubCellFaceFlows( this%cornerXComponentIndexes( cornerIndex, 4 ) )
        do m = 1, 3
          sumArea = sumArea + &
              neighborSubCellFaceAreas( 1, this%cornerXComponentIndexes( cornerIndex, m ) )
          sumFlowRate = sumFlowRate + & 
              neighborSubCellFaceFlows( this%cornerXComponentIndexes( cornerIndex, 4 ) , &
                                        this%cornerXComponentIndexes( cornerIndex, m ) )  
        end do
      case (2)
        sumArea     = this%SubCellData%DX*this%SubCellData%DZ
        sumFlowRate = centerSubCellFaceFlows( this%cornerYComponentIndexes( cornerIndex, 4 ) )
        do m = 1, 3
          sumArea = sumArea + &
              neighborSubCellFaceAreas( 2, this%cornerYComponentIndexes( cornerIndex, m ) )
          sumFlowRate = sumFlowRate + & 
              neighborSubCellFaceFlows( this%cornerYComponentIndexes( cornerIndex, 4 ) , &
                                        this%cornerYComponentIndexes( cornerIndex, m ) )  
        end do
      case (3)
        sumArea     = this%SubCellData%DX*this%SubCellData%DY
        sumFlowRate = centerSubCellFaceFlows( this%cornerZComponentIndexes( cornerIndex, 4 ) )
        do m = 1, 3
          sumArea = sumArea + &
              neighborSubCellFaceAreas( 3, this%cornerZComponentIndexes( cornerIndex, m ) )
          sumFlowRate = sumFlowRate + & 
              neighborSubCellFaceFlows( this%cornerZComponentIndexes( cornerIndex, 4 ) , &
                                        this%cornerZComponentIndexes( cornerIndex, m ) )  
        end do
    end select 
 

    ! flowContributionFactor=0.25 is cancelled implicitly    
    qCorner = sumFlowRate/sumArea
  
  
  end function pr_GetInterpolatedCornerDischarge



  subroutine pr_SetCornerComponentsIndexes( this ) 
      !-----------------------------------------------------------------
      ! Set neighbor sub cells indexes for computation of 
      ! x,y,z components of corner quantities
      ! 
      ! Convention for sub cell corners
      !     1: 000
      !     2: 100
      !     3: 010
      !     4: 110
      !     5: 001
      !     6: 101
      !     7: 011
      !     8: 111
      !
      ! dev: find a proper place for doing this only once
      !-----------------------------------------------------------------
      implicit none
      class( TrackSubCellType ) :: this
      !-----------------------------------------------------------------


      ! Cell indexes for x-component
      ! 000: 7-8-13   : fn1 : 
      ! 100: 7-8-13   : fn2 : 
      ! 010: 10-11-13 : fn1 : 
      ! 110: 10-11-13 : fn2 :  
      ! 001: 7-9-16   : fn1 :
      ! 101: 7-9-16   : fn2 :
      ! 011: 10-12-16 : fn1 :
      ! 111: 10-12-16 : fn2 :
      this%cornerXComponentIndexes(1,:) = [ 7, 8,13,1] 
      this%cornerXComponentIndexes(2,:) = [ 7, 8,13,2] 
      this%cornerXComponentIndexes(3,:) = [10,11,13,1] 
      this%cornerXComponentIndexes(4,:) = [10,11,13,2] 
      this%cornerXComponentIndexes(5,:) = [ 7, 9,16,1] 
      this%cornerXComponentIndexes(6,:) = [ 7, 9,16,2] 
      this%cornerXComponentIndexes(7,:) = [10,12,16,1] 
      this%cornerXComponentIndexes(8,:) = [10,12,16,2] 


      ! Cell indexes for y-component
      !000: 1-13-14 : fn3 :
      !100: 4-13-15 : fn3 :
      !010: 1-13-14 : fn4 :
      !110: 4-13-15 : fn4 :
      !001: 1-16-17 : fn3 :
      !101: 4-16-18 : fn3 :
      !011: 1-16-17 : fn4 :
      !111: 4-16-18 : fn4 :
      this%cornerYComponentIndexes(1,:) = [ 1,13,14,3] 
      this%cornerYComponentIndexes(2,:) = [ 4,13,15,3] 
      this%cornerYComponentIndexes(3,:) = [ 1,13,14,4] 
      this%cornerYComponentIndexes(4,:) = [ 4,13,15,4] 
      this%cornerYComponentIndexes(5,:) = [ 1,16,17,3] 
      this%cornerYComponentIndexes(6,:) = [ 4,16,18,3] 
      this%cornerYComponentIndexes(7,:) = [ 1,16,17,4] 
      this%cornerYComponentIndexes(8,:) = [ 4,16,18,4]


      ! Cell indexes for z-component
      !000: 1-2-7   : fn5 :
      !100: 4-5-7   : fn5 :
      !010: 1-3-10  : fn5 :
      !110: 4-6-10  : fn5 :
      !001: 1-2-7   : fn6 :
      !101: 4-5-7   : fn6 :
      !011: 1-3-10  : fn6 :
      !111: 4-6-10  : fn6 :
      this%cornerZComponentIndexes(1,:) = [ 1, 2, 7,5] 
      this%cornerZComponentIndexes(2,:) = [ 4, 5, 7,5] 
      this%cornerZComponentIndexes(3,:) = [ 1, 3,10,5] 
      this%cornerZComponentIndexes(4,:) = [ 4, 6,10,5] 
      this%cornerZComponentIndexes(5,:) = [ 1, 2, 7,6] 
      this%cornerZComponentIndexes(6,:) = [ 4, 5, 7,6] 
      this%cornerZComponentIndexes(7,:) = [ 1, 3,10,6] 
      this%cornerZComponentIndexes(8,:) = [ 4, 6,10,6]


      return


  end subroutine pr_SetCornerComponentsIndexes


  subroutine pr_ComputeCornerPorosityUniform( this, neighborSubCellVolume, neighborSubCellPorosity )
    !----------------------------------------------------------------
    ! Assign uniform porosity 
    !----------------------------------------------------------------
    ! Specifications
    !----------------------------------------------------------------
    implicit none
    class( TrackSubCellType ) :: this
    ! input
    doubleprecision, dimension(18), intent(in) :: neighborSubCellVolume   
    doubleprecision, dimension(18), intent(in) :: neighborSubCellPorosity
    !-----------------------------------------------------------------

    ! spatially constant porosity
    this%porosity000 = this%SubCellData%Porosity   
    this%porosity100 = this%SubCellData%Porosity
    this%porosity010 = this%SubCellData%Porosity       
    this%porosity110 = this%SubCellData%Porosity
    this%porosity001 = this%SubCellData%Porosity
    this%porosity101 = this%SubCellData%Porosity       
    this%porosity011 = this%SubCellData%Porosity
    this%porosity111 = this%SubCellData%Porosity

    return
    
  end subroutine pr_ComputeCornerPorosityUniform


  ! RWPT
  subroutine pr_ComputeCornerPorosityInterpolated( this, neighborSubCellVolume, neighborSubCellPorosity )
    !----------------------------------------------------------------
    ! Calculate porosities at the corners as the volume weighted 
    ! averaged porosity
    !----------------------------------------------------------------
    ! Specifications
    !----------------------------------------------------------------
    implicit none
    class( TrackSubCellType ) :: this
    ! input
    doubleprecision, dimension(18), intent(in) :: neighborSubCellVolume   
    doubleprecision, dimension(18), intent(in) :: neighborSubCellPorosity
    !-----------------------------------------------------------------

    ! Assign interpolated values
    this%porosity000 = this%GetInterpolatedCornerPorosity( neighborSubCellVolume, neighborSubCellPorosity, 1 )
    this%porosity100 = this%GetInterpolatedCornerPorosity( neighborSubCellVolume, neighborSubCellPorosity, 2 )
    this%porosity010 = this%GetInterpolatedCornerPorosity( neighborSubCellVolume, neighborSubCellPorosity, 3 )
    this%porosity110 = this%GetInterpolatedCornerPorosity( neighborSubCellVolume, neighborSubCellPorosity, 4 )
    this%porosity001 = this%GetInterpolatedCornerPorosity( neighborSubCellVolume, neighborSubCellPorosity, 5 )
    this%porosity101 = this%GetInterpolatedCornerPorosity( neighborSubCellVolume, neighborSubCellPorosity, 6 )
    this%porosity011 = this%GetInterpolatedCornerPorosity( neighborSubCellVolume, neighborSubCellPorosity, 7 )
    this%porosity111 = this%GetInterpolatedCornerPorosity( neighborSubCellVolume, neighborSubCellPorosity, 8 )

    ! Done  
    return
      
  end subroutine pr_ComputeCornerPorosityInterpolated


  ! RWPT-USG
  function pr_GetInterpolatedCornerPorosity( this, neighborSubCellVolume, &
            neighborSubCellPorosity, cornerIndex ) result( cornerPorosity ) 
  !-------------------------------------------------------------------------
  ! Compute equivalent corner porosity for a given cornerIndex, 
  ! using source information in neighborSubCellVolume and 
  ! neighborSubCellPorosity
  !
  ! Equivalent corner porosity is a volume weighted average
  ! between 7 cells: 6 neighbors and center cell
  !-------------------------------------------------------------------------
  implicit none
  class( TrackSubCellType ) :: this
  ! input
  doubleprecision, dimension(18), intent(in) :: neighborSubCellVolume
  doubleprecision, dimension(18), intent(in) :: neighborSubCellPorosity
  integer, intent(in) :: cornerIndex
  ! output
  doubleprecision     :: cornerPorosity
  ! local 
  doubleprecision :: sumVolume, sumWeightedPorosity
  integer :: m
  !------------------------------------------------------------------------

    ! self
    sumVolume           = this%SubCellData%DX*this%SubCellData%DY*this%SubCellData%DZ
    sumWeightedPorosity = this%SubCellData%Porosity*sumVolume
    ! neighbors
    do m = 1, 6
      sumVolume = sumVolume + &
        neighborSubCellVolume( this%cornerPorosityIndexes( cornerIndex, m ) )
      sumWeightedPorosity = sumWeightedPorosity + & 
        neighborSubCellVolume( this%cornerPorosityIndexes( cornerIndex, m ) )*&
        neighborSubCellPorosity( this%cornerPorosityIndexes( cornerIndex, m ) ) 
    end do
    cornerPorosity = sumWeightedPorosity/sumVolume

    ! Done
    return

  end function pr_GetInterpolatedCornerPorosity


  subroutine pr_SetCornerPorosityIndexes( this ) 
   !-----------------------------------------------------------------
   ! Set neighbor sub cells indexes for computation of 
   ! equivalent corner porosities
   ! 
   ! Convention for sub cell corners
   !     1: 000
   !     2: 100
   !     3: 010
   !     4: 110
   !     5: 001
   !     6: 101
   !     7: 011
   !     8: 111
   !
   ! dev: find a proper place for doing this only once
   !-----------------------------------------------------------------
   implicit none
   class( TrackSubCellType ) :: this
   !-----------------------------------------------------------------

      ! Cell indexes for corner porosity
      this%cornerPorosityIndexes(1,:) = [ 1, 2,  7,  8, 13, 14 ]
      this%cornerPorosityIndexes(2,:) = [ 4, 5,  7,  8, 13, 15 ]
      this%cornerPorosityIndexes(3,:) = [ 1, 3, 10, 11, 13, 14 ]
      this%cornerPorosityIndexes(4,:) = [ 4, 6, 10, 12, 16, 18 ]
      this%cornerPorosityIndexes(5,:) = [ 1, 2,  7,  9, 16, 17 ]
      this%cornerPorosityIndexes(6,:) = [ 4, 5,  7,  9, 16, 18 ]
      this%cornerPorosityIndexes(7,:) = [ 1, 3, 10, 12, 16, 17 ]
      this%cornerPorosityIndexes(8,:) = [ 4, 6, 10, 11, 13, 15 ]

      ! Done
      return

  end subroutine pr_SetCornerPorosityIndexes


  ! RWPT
  subroutine pr_TrilinearDerivative( this, direction, x, y, z, v000, v100, v010, v110, v001, v101, v011, v111, output )
      !-----------------------------------------------------------
      ! Compute derivative of a trilinear interpolation in a given
      ! direction
      ! 
      ! Params
      !     - x, y, z   : local cell coordinates 
      !     - direction : specifies x=1, y=2 or z=3 direction
      !     - vijk      : values at corresponding corners
      !     - output    : the output variable
      !-----------------------------------------------------------
      ! Specifications
      !-----------------------------------------------------------
      implicit none
      class (TrackSubCellType) :: this
      ! input
      integer        , intent(in):: direction
      doubleprecision, intent(in):: x, y, z
      doubleprecision, intent(in):: v000, v100, v010, v110, v001, v101, v011, v111
      ! output
      doubleprecision, intent(inout) :: output
      ! local
      doubleprecision :: v0, v1, v00, v10, v01, v11
      !-----------------------------------------------------------

      ! Initialize
      output = 0d0

      select case (direction)
          ! x direction
          case (1)
              v00    = ( v100 - v000 )/this%SubCellData%DX
              v01    = ( v101 - v001 )/this%SubCellData%DX
              v10    = ( v110 - v010 )/this%SubCellData%DX
              v11    = ( v111 - v011 )/this%SubCellData%DX
              v0     = ( 1.0d0 - y )*v00 + y*v10
              v1     = ( 1.0d0 - y )*v01 + y*v11
              output = ( 1.0d0 - z )*v0  + z*v1 
              return 
          ! y direction
          case (2)
              v00    = ( v010 - v000 )/this%SubCellData%DY
              v01    = ( v011 - v001 )/this%SubCellData%DY
              v10    = ( v110 - v100 )/this%SubCellData%DY
              v11    = ( v111 - v101 )/this%SubCellData%DY
              v0     = ( 1.0d0 - x )*v00 + x*v10
              v1     = ( 1.0d0 - x )*v01 + x*v11
              output = ( 1.0d0 - z )*v0  + z*v1 
              return 
          ! z direction
          case (3) 
              v00    = ( v001 - v000 )/this%SubCellData%DZ
              v01    = ( v011 - v010 )/this%SubCellData%DZ
              v10    = ( v101 - v100 )/this%SubCellData%DZ
              v11    = ( v111 - v110 )/this%SubCellData%DZ
              v0     = ( 1.0d0 - x )*v00 + x*v10
              v1     = ( 1.0d0 - x )*v01 + x*v11
              output = ( 1.0d0 - y )*v0  + y*v1 
              return 
      end select


  end subroutine pr_TrilinearDerivative


  subroutine pr_TrilinearDerivativeX( this, x, y, z, v000, v100, v010, v110, v001, v101, v011, v111, output )
  !-----------------------------------------------------------
  ! Compute derivative of a trilinear interpolation in a given
  ! direction
  ! 
  ! Params
  !     - x, y, z   : local cell coordinates 
  !     - vijk      : values at corresponding corners
  !     - output    : the output variable
  !-----------------------------------------------------------
  ! Specifications
  !-----------------------------------------------------------
  implicit none
  class (TrackSubCellType) :: this
  ! input
  doubleprecision, intent(in):: x, y, z
  doubleprecision, intent(in):: v000, v100, v010, v110, v001, v101, v011, v111
  ! output
  doubleprecision, intent(inout) :: output
  ! local
  doubleprecision :: v0, v1, v00, v10, v01, v11
  !-----------------------------------------------------------

      ! x direction
      v00    = ( v100 - v000 )/this%SubCellData%DX
      v01    = ( v101 - v001 )/this%SubCellData%DX
      v10    = ( v110 - v010 )/this%SubCellData%DX
      v11    = ( v111 - v011 )/this%SubCellData%DX
      v0     = ( 1.0d0 - y )*v00 + y*v10
      v1     = ( 1.0d0 - y )*v01 + y*v11
      output = ( 1.0d0 - z )*v0  + z*v1

      ! Done 
      return 

  end subroutine pr_TrilinearDerivativeX


  subroutine pr_TrilinearDerivativeY( this, x, y, z, v000, v100, v010, v110, v001, v101, v011, v111, output )
  !-----------------------------------------------------------
  ! Compute derivative of a trilinear interpolation in a given
  ! direction
  ! 
  ! Params
  !     - x, y, z   : local cell coordinates 
  !     - vijk      : values at corresponding corners
  !     - output    : the output variable
  !-----------------------------------------------------------
  ! Specifications
  !-----------------------------------------------------------
  implicit none
  class (TrackSubCellType) :: this
  ! input
  doubleprecision, intent(in):: x, y, z
  doubleprecision, intent(in):: v000, v100, v010, v110, v001, v101, v011, v111
  ! output
  doubleprecision, intent(inout) :: output
  ! local
  doubleprecision :: v0, v1, v00, v10, v01, v11
  !-----------------------------------------------------------

    ! y direction
    v00    = ( v010 - v000 )/this%SubCellData%DY
    v01    = ( v011 - v001 )/this%SubCellData%DY
    v10    = ( v110 - v100 )/this%SubCellData%DY
    v11    = ( v111 - v101 )/this%SubCellData%DY
    v0     = ( 1.0d0 - x )*v00 + x*v10
    v1     = ( 1.0d0 - x )*v01 + x*v11
    output = ( 1.0d0 - z )*v0  + z*v1
    
    ! Done 
    return 

  end subroutine pr_TrilinearDerivativeY


  subroutine pr_TrilinearDerivativeZ( this, x, y, z, v000, v100, v010, v110, v001, v101, v011, v111, output )
  !-----------------------------------------------------------
  ! Compute derivative of a trilinear interpolation in a given
  ! direction
  ! 
  ! Params
  !     - x, y, z   : local cell coordinates 
  !     - vijk      : values at corresponding corners
  !     - output    : the output variable
  !-----------------------------------------------------------
  ! Specifications
  !-----------------------------------------------------------
  implicit none
  class (TrackSubCellType) :: this
  ! input
  doubleprecision, intent(in):: x, y, z
  doubleprecision, intent(in):: v000, v100, v010, v110, v001, v101, v011, v111
  ! output
  doubleprecision, intent(inout) :: output
  ! local
  doubleprecision :: v0, v1, v00, v10, v01, v11
  !-----------------------------------------------------------

    ! z direction
    v00    = ( v001 - v000 )/this%SubCellData%DZ
    v01    = ( v011 - v010 )/this%SubCellData%DZ
    v10    = ( v101 - v100 )/this%SubCellData%DZ
    v11    = ( v111 - v110 )/this%SubCellData%DZ
    v0     = ( 1.0d0 - x )*v00 + x*v10
    v1     = ( 1.0d0 - x )*v01 + x*v11
    output = ( 1.0d0 - y )*v0  + y*v1

    ! Done 
    return 

  end subroutine pr_TrilinearDerivativeZ


  subroutine pr_Trilinear( this, x, y, z, v000, v100, v010, v110, v001, v101, v011, v111, output )
  !-----------------------------------------------------------
  ! Compute trilinear interpolation of corner values
  ! into the given coordinates 
  !  
  ! Params
  !     - x, y, z : coordinates to interpolate
  !     - vijk    : values at corresponding corners
  !     - output  : the output variable
  !-----------------------------------------------------------
  ! Specifications
  !-----------------------------------------------------------
  implicit none
  class (TrackSubCellType) :: this
  ! input
  doubleprecision, intent(in):: x, y, z
  doubleprecision, intent(in):: v000, v100, v010, v110, v001, v101, v011, v111
  ! output
  doubleprecision, intent(inout) :: output
  ! local
  doubleprecision :: v0, v1, v00, v10, v01, v11
  !-----------------------------------------------------------
      
    v00    = ( 1.0d0 - x )*v000 + x*v100
    v01    = ( 1.0d0 - x )*v001 + x*v101
    v10    = ( 1.0d0 - x )*v010 + x*v110
    v11    = ( 1.0d0 - x )*v011 + x*v111
    v0     = ( 1.0d0 - y )*v00  + y*v10
    v1     = ( 1.0d0 - y )*v01  + y*v11
    output = ( 1.0d0 - z )*v0   + z*v1

    ! Done 
    return 

  end subroutine pr_Trilinear


  ! RWPT
  subroutine pr_ComputeNonlinearDispersivities( this, vx, vy, vz, Daqueous, distance, &
                                                  delta, betaL, betaT, alphaL, alphaT )
  !----------------------------------------------------------------
  ! Compute nonlinear equivalent dispersivities
  ! establishing analogy with 
  ! model from  Chiogna et al. 2010, Rolle et al. 2013
  !  
  ! Params:
  !----------------------------------------------------------------
  ! Specifications
  !----------------------------------------------------------------
  implicit none
  class (TrackSubCellType)    :: this
  ! input
  doubleprecision, intent(in) :: vx, vy, vz, Daqueous
  doubleprecision, intent(in) :: distance, delta, betaL, betaT
  ! output
  doubleprecision, intent(inout) :: alphaL, alphaT
  ! local
  doubleprecision :: v, peclet, fbase
  !---------------------------------------------------------------- 
        
    v      = 0d0
    peclet = 0d0
    fbase  = 0d0  


    v      = sqrt( vx**2 + vy**2 + vz**2 )
    peclet = v*distance/Daqueous
    fbase  = peclet**2/( peclet + 2 + 4*delta**2 )

    alphaL = (distance/peclet)*( fbase )**betaL 
    alphaT = (distance/peclet)*( fbase )**betaT

    return

  end subroutine pr_ComputeNonlinearDispersivities


end module TrackSubCellModule
