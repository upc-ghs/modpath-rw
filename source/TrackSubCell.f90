module TrackSubCellModule
  use ParticleLocationModule,only : ParticleLocationType
  use TrackSubCellResultModule,only : TrackSubCellResultType
  use ModpathSubCellDataModule,only : ModpathSubCellDataType

  ! RWPT 
  use ParticleTrackingOptionsModule,only : ParticleTrackingOptionsType
  use ModpathCellDataModule,only : ModpathCellDataType
  use omp_lib ! mainly debugging
  implicit none
  
! Set default access status to private
  private

  type,public :: TrackSubCellType
    type(ModpathSubCellDataType) :: SubCellData

    ! RWPT
    ! vx, vy, vz, vnorm
    doubleprecision, dimension(4) :: vCorner000
    doubleprecision, dimension(4) :: vCorner100
    doubleprecision, dimension(4) :: vCorner010
    doubleprecision, dimension(4) :: vCorner110
    doubleprecision, dimension(4) :: vCorner001
    doubleprecision, dimension(4) :: vCorner101
    doubleprecision, dimension(4) :: vCorner011
    doubleprecision, dimension(4) :: vCorner111

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

    ! Corner components indexes
    ! 8 corners, 3 subcell indexes, 1 face id
    integer, dimension(8,4) :: cornerXComponentIndexes, &
                               cornerYComponentIndexes, & 
                               cornerZComponentIndexes
    

    ! Prototype that might be used in production 
    !doubleprecision, dimension(0:1,0:1,0:1) :: cornerPorosity

    ! RWPT Pointers
    procedure(Advection), pass, pointer :: AdvectionDisplacement=>null()
    procedure(ExitFaceAndTimeStep), pass, pointer :: ExitFaceAndUpdateTimeStep=>null()

  contains
    procedure,private :: CalculateDT=>pr_CalculateDT
    procedure,private :: NewXYZ=>pr_NewXYZ
    procedure :: ExecuteTracking=>pr_ExecuteTracking

    ! RWPT
    procedure :: ExecuteRandomWalkParticleTracking=>pr_ExecuteRandomWalkParticleTracking
    procedure :: LinearInterpolationVelocities=>pr_LinearInterpolationVelocities
    procedure :: ComputeRandomWalkTimeStep=>pr_ComputeRandomWalkTimeStep
    procedure :: ComputeCornerVelocities=>pr_ComputeCornerVelocities
    procedure :: ComputeCornerDischarge=>pr_ComputeCornerDischarge
    procedure :: GetInterpolatedCornerDischarge=>pr_GetInterpolatedCornerDischarge
    procedure :: SetCornerComponentsIndexes=>pr_SetCornerComponentsIndexes
    procedure :: Trilinear=>pr_Trilinear
    procedure :: TrilinearDerivative=>pr_TrilinearDerivative
    procedure :: DispersionDivergence=>pr_DispersionDivergence
    procedure :: DispersionDivergenceDischarge=>pr_DispersionDivergenceDischarge
    procedure :: ComputeCornerPorosity=>pr_ComputeCornerPorosity
    procedure :: DisplacementRandom=>pr_DisplacementRandom
    procedure :: DisplacementRandomDischarge=>pr_DisplacementRandomDischarge
    procedure :: GenerateStandardNormalRandom=>pr_GenerateStandardNormalRandom
    procedure :: AdvectionDisplacementExponential=>pr_AdvectionDisplacementExponential
    procedure :: AdvectionDisplacementEulerian=>pr_AdvectionDisplacementEulerian
    procedure :: NewtonRaphsonTimeStep=>NewtonRaphsonTimeStepExponentialAdvection
  end type


  ! RWPT
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

  end interface


contains
!-------------------------------------------------------------------
  subroutine pr_ExecuteTracking(this,stopIfNoExit,initialLocation,maximumTime, trackingResult)
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
  integer :: exitFace,exitStatus
  integer :: statusVX,statusVY,statusVZ

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
  subroutine pr_ExecuteRandomWalkParticleTracking(this,stopIfNoExit, &
          initialLocation,maximumTime,trackingResult,trackingOptions)
      !------------------------------------------------------------
      ! Function moves particles following RWPT protocol.
      ! It has been adapted from ExecuteTracking but logic
      ! differs significantly as particles are moved numerically 
      ! instead of analytically
      ! 
      ! Params:
      !     - ToDo
      !
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
      integer :: exitFace,exitStatus
      integer :: statusVX,statusVY,statusVZ

      ! RWPT
      type(ParticleTrackingOptionsType),intent(in) :: trackingOptions
      doubleprecision :: alphaT, alphaL, Dmol
      doubleprecision :: divDx, divDy, divDz
      doubleprecision :: dAdvx, dAdvy, dAdvz
      doubleprecision :: dBx, dBy, dBz
      doubleprecision :: dx, dy, dz
      doubleprecision :: nx, ny, nz
      doubleprecision :: nnx, nny, nnz ! remove
      doubleprecision :: dxrw, dyrw, dzrw
      doubleprecision :: drwtol = 1.0d-14
      logical         :: continueTimeLoop
      logical         :: reachedMaximumTime
      logical         :: twoDimensions
      doubleprecision :: dtold
      doubleprecision, dimension(3) :: dts
      doubleprecision, dimension(3) :: dtxyz
      integer :: dtLoopCounter, posRestartCounter
      !------------------------------------------------------------

      ! Needs cleaning

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

      ! Assign pointers
      if ( trackingOptions%advectionKind .eq. 1 ) then 
          this%AdvectionDisplacement=>pr_AdvectionDisplacementExponential
          this%ExitFaceAndUpdateTimeStep=>pr_DetectExitFaceAndUpdateTimeStepNewton
      else if ( trackingOptions%advectionKind .eq. 2 ) then 
          this%AdvectionDisplacement=>pr_AdvectionDisplacementEulerian
          this%ExitFaceAndUpdateTimeStep=>pr_DetectExitFaceAndUpdateTimeStepQuadratic
      end if

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

      ! Initialize displacements
      dxrw = 0d0
      dyrw = 0d0
      dzrw = 0d0

      ! Initialize kind of domain solver
      twoDimensions = trackingOptions%twoDimensions

      ! Initialize dispersion
      alphaL = trackingOptions%alphaL
      alphaT = trackingOptions%alphaT
      Dmol   = trackingOptions%Dmol

      ! Compute time step for RWPT
      call this%ComputeRandomWalkTimeStep( trackingOptions, dt )

      ! Initializes current time
      t     = initialTime
      dtold = dt

      ! Something wrong, leave
      if ( dt .eq. 0 ) then 
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
      posRestartCounter = 0


      !!!!!!!!!!1 DEV/REMOVE !!!!!!!!!!!!!!
      !dBx = 0
      !dBy = 0
      !dBz = 0
      !print *, '********************************************************'
      !print *, '** TrackSubCell: cellNumber: ', cellNumber
      !print *, '** TrackSubCell: initialTime: ', initialTime
      !print *, '** TrackSubCell: LOOP ...'
      !!print *, '********************************************************'
      !!!!!!!!!!1 DEV/REMOVE !!!!!!!!!!!!!!


      ! Local cell time loop 
      exitFace = 0
      continueTimeLoop = .true.
      reachedMaximumTime = .false.
      do while( continueTimeLoop )

          ! Update current time
          t = t + dt

          ! Recompute dt for maximumTime 
          if (maximumTime .lt. t) then
              t  = t - dt
              dt = maximumTime - t
              t  = maximumTime
              reachedMaximumTime = .true.
          end if 

          ! Compute RWPT movement
          ! Discharge formulation
          call this%LinearInterpolationVelocities( x, y, z, vx, vy, vz )
          call this%DispersionDivergenceDischarge( x, y, z, alphaL, alphaT, Dmol, divDx, divDy, divDz )
          call this%DisplacementRandomDischarge( x, y, z, alphaL, alphaT, Dmol, dBx, dBy, dBz )
          call this%AdvectionDisplacement( x, y, z, dt, vx, vy, vz, dAdvx, dAdvy, dAdvz )
          ! Velocities formulation
          !call this%LinearInterpolationVelocities( x, y, z, vx, vy, vz )
          !call this%DispersionDivergence( x, y, z, alphaL, alphaT, Dmol, divDx, divDy, divDz )
          !call this%DisplacementRandom( x, y, z, alphaL, alphaT, Dmol, dBx, dBy, dBz )
          !call this%AdvectionDisplacement( x, y, z, dt, vx, vy, vz, dAdvx, dAdvy, dAdvz )
          dxrw = dAdvx + divDx*dt + dBx*sqrt( dt )
          nx   = x + dxrw/dx
          dyrw = dAdvy + divDy*dt + dBy*sqrt( dt )
          ny   = y + dyrw/dy
          if ( .not. twoDimensions ) then 
              dzrw = dAdvz + divDz*dt + dBz*sqrt( dt )
              nz   = z + dzrw/dz
          end if

          print *, '** TrackSubCell: TIME', t
          print *, '** TrackSubCell: DT', dt
          print *, '** TrackSubCell: dAdvx, divDx, dBx', dAdvx/dx, divDx*dt/dx, dBx*sqrt(dt)/dx
          print *, '** TrackSubCell: dAdvy, divDy, dBy', dAdvy/dy, divDy*dt/dy, dBy*sqrt(dt)/dy

          ! Detect if particle leaving the cell
          ! and force the particle into exactly one
          ! interface by computing the required dt
          do while (                                       &
              ( nx .gt. 1.0d0 ) .or. ( nx .lt. 0d0 )  .or. &
              ( ny .gt. 1.0d0 ) .or. ( ny .lt. 0d0 )  .or. &
              ( nz .gt. 1.0d0 ) .or. ( nz .lt. 0d0 )       & 
          )
              
              print *, '** TrackSubCell: LEAVING SUBCELL', nx, ny, nz

              dtxyz(:) = 0d0

              ! Recompute dt for exact interface
              call this%ExitFaceAndUpdateTimeStep( x, y, z, nx, ny, nz, &
                                       vx, vy, vz, divDx, divDy, divDz, &
                                       dBx, dBy, dBz, t, dt, dtxyz, exitFace )

              ! Given new dt, recompute advection displacements
              call this%AdvectionDisplacement( x, y, z, dt, vx, vy, vz, & 
                                                    dAdvx, dAdvy, dAdvz )

              print *, '** TrackSubCell: exitFace, t, dt', exitFace, t, dt


              ! If maximumTime was reached, but particle left
              ! the cell, then the condition is resetted
              if (reachedMaximumTime) then
                  reachedMaximumTime = .false.
              end if

              ! Find new RWPT displacements
              if ( ( exitFace .eq. 1 ) .or. ( exitFace .eq. 2 ) ) then
                  dyrw = dAdvy + divDy*dt + dBy*sqrt( dt )
                  ny   = y + dyrw/dy
                  if ( .not. twoDimensions ) then
                      dzrw = dAdvz + divDz*dt + dBz*sqrt( dt )
                      nz   = z + dzrw/dz
                  end if
                  nx = 1.0d0
                  if ( exitFace .eq. 1 ) nx=0d0
              else if ( ( exitFace .eq. 3 ) .or. ( exitFace .eq. 4 ) ) then 
                  dxrw = dAdvx + divDx*dt + dBx*sqrt( dt )
                  nx   = x + dxrw/dx
                  if ( .not. twoDimensions ) then 
                      dzrw = dAdvz + divDz*dt + dBz*sqrt( dt )
                      nz   = z + dzrw/dz
                  end if
                  ny = 1.0d0
                  if ( exitFace .eq. 3 ) ny=0d0
              else if ( ( exitFace .eq. 5 ) .or. ( exitFace .eq. 6 ) ) then
                  dxrw = dAdvx + divDx*dt + dBx*sqrt( dt )
                  dyrw = dAdvy + divDy*dt + dBy*sqrt( dt )
                  nx = x + dxrw/dx
                  ny = y + dyrw/dy
                  z = 1.0d0
                  if ( exitFace .eq. 5 ) nz=0d0
              else
                  print *, '** TrackSubCell: RESET_ODK'
                  print *, nx, ny, nz
                  ! Restart nx, ny, nz and try again
                  ! if not a valid time step and exitFace
                  nx = x
                  ny = y
                  nz = z
                  t  = t - dt
                  dt = dtold
                  dtLoopCounter = 0
                  posRestartCounter = posRestartCounter + 1
                  print *, nx, ny, nz
                  exit
              end if

              ! Restart if method did not
              ! found a valid position after two tries.
              dtLoopCounter = dtLoopCounter + 1
              if ( dtLoopCounter .eq. 2 ) then
                  ! When using exponential integration 
                  ! and rwpt, it has been observed that 
                  ! new positions can pivot between invalid 
                  ! positions which happened when not forcing NR 
                  ! to return a new dt smaller than the previous.
                  ! This block controls that if that happens, 
                  ! positions are restarted after two tries, 
                  
                  print *, '** TrackSubCell: RESET_OUT'

                  ! Restart nx, ny, nz and try again
                  nx = x
                  ny = y
                  nz = z
                  t  = t - dt
                  dt = dtold
                  exitFace = 0
                  dtLoopCounter = 0
                  posRestartCounter = posRestartCounter + 1
                  exit
              end if

          end do

          ! Report and leave
          if ( (reachedMaximumTime) .or. ( exitFace .ne. 0 ) ) then
              if (reachedMaximumTime) then
                  ! In this case it is allowed for a particle
                  ! to reach the maximum time and eventually have an exit face
                  trackingResult%Status = trackingResult%Status_ReachedMaximumTime()
              else
                  trackingResult%ExitFaceConnection = this%SubCellData%Connection(exitFace)
                  if(this%SubCellData%Connection(exitFace) .lt. 0) then
                      ! This is the case for unstructured grid
                      trackingResult%Status = trackingResult%Status_ExitAtInternalFace()
                  else
                      trackingResult%Status = trackingResult%Status_ExitAtCellFace()
                  end if
              end if
              trackingResult%ExitFace = exitFace
              trackingResult%FinalLocation%CellNumber = cellNumber
              trackingResult%FinalLocation%LocalX = nx
              trackingResult%FinalLocation%LocalY = ny
              trackingResult%FinalLocation%LocalZ = nz

              if (.not. trackingResult%FinalLocation%Valid() ) then 
                  print *, 'DEBUG: ** TrackSubCell'
                  print *, 'DEBUG: xyz', x, y, z
                  print *, 'DEBUG: nxyz', nx, ny, nz
                  print *, 'DEBUG: exitFace', exitFace
              end if

              trackingResult%FinalLocation%TrackingTime = t
              continueTimeLoop = .false.

              !print *, '** TrackSubCell: FINAL TIME', t

              return

          end if

          ! Update particle positions
          x = nx
          y = ny
          z = nz

      end do 

      return

  end subroutine pr_ExecuteRandomWalkParticleTracking


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
      doubleprecision :: x, y, z
      ! output
      doubleprecision :: vx, vy, vz
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

      ! Compute time step
      select case (trackingOptions%timeStepKind)
          case (1)
              ! Fix Courant 
              dt = trackingOptions%timeStepParameters(1)/( &
                  max(abs(vx1), abs(vx2))/dx +             &
                  max(abs(vy1), abs(vy2))/dy +             &
                  max(abs(vz1), abs(vz2))/dz )
          case (2)
              ! Fix Peclet 
              dt = 1/( trackingOptions%timeStepParameters(2)*(                  &
                  trackingOptions%alphaL*max(abs(vx1), abs(vx2))/( dx**2 ) +    &
                  trackingOptions%alphaT*max(abs(vy1), abs(vy2))/( dy**2 ) +    &
                  trackingOptions%alphaT*max(abs(vz1), abs(vz2))/( dz**2 ) ) )
          case (3)
              ! Courant condition
              dts(1) = trackingOptions%timeStepParameters(1)/( & 
                  max(abs(vx1), abs(vx2))/dx +                 &
                  max(abs(vy1), abs(vy2))/dy +                 &
                  max(abs(vz1), abs(vz2))/dz )
              ! Peclet condition
              dts(2) = 1/( trackingOptions%timeStepParameters(2)*(              &
                  trackingOptions%alphaL*max(abs(vx1), abs(vx2))/( dx**2 ) +    & 
                  trackingOptions%alphaT*max(abs(vy1), abs(vy2))/( dy**2 ) +    &
                  trackingOptions%alphaT*max(abs(vz1), abs(vz2))/( dz**2 ) ) )
              ! Compute minimum
              dt     = minval( dts, dts > 0 )
      end select


  end subroutine pr_ComputeRandomWalkTimeStep


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
      doubleprecision :: AFace, BFace, z1, z2, zsqrt
      !----------------------------------------------------------------

      ! Initialize
      AFace      = 0d0
      BFace      = 0d0
      z1         = 0d0
      z2         = 0d0
      zsqrt      = 0d0 
      dInterface = 0d0
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
          if ( dInterface .eq. 0.0 ) then
              exitFace = exitFaceX
              dt = 0.0
              return
          end if

          ! Coefficients
          AFace = vx + divDx
          BFace = dBx

          ! Given dInterface, compute new dt
          zsqrt = sqrt( BFace**2 + 4*dInterface*AFace )
          z1    = (-BFace + zsqrt )/( 2*AFace )
          z2    = (-BFace - zsqrt )/( 2*AFace )

          ! Compute new dt
          if ( ( z1 .gt. 0d0 ) .and. ( z2 .gt. 0d0 ) ) then
              dtxyz(1) = min( z1**2, z2**2 )
          else if ( z1 .gt. 0d0 ) then
              dtxyz(1) = z1**2
          else if ( z2 .gt. 0d0 ) then
              dtxyz(1) = z2**2
          end if

          ! If computed dt is higher than current
          ! is not valid, set to zero
          if ( dtxyz(1) .gt. dt ) dtxyz(1) = 0d0 

      end if


      ! Leaving through y face
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
          zsqrt = sqrt( BFace**2 + 4*dInterface*AFace )
          z1    = (-BFace + zsqrt )/( 2*AFace )
          z2    = (-BFace - zsqrt )/( 2*AFace )

          ! Determine new dt
          if ( ( z1 .gt. 0d0 ) .and. ( z2 .gt. 0d0 ) ) then
              dtxyz(2) = min( z1**2, z2**2 )
          else if ( z1 .gt. 0d0 ) then
              dtxyz(2) = z1**2
          else if ( z2 .gt. 0d0 ) then
              dtxyz(2) = z2**2
          end if

          ! If computed dt is higher than current
          ! is not valid, set to zero
          if ( dtxyz(2) .gt. dt ) dtxyz(2) = 0d0 

      end if


      ! Leaving through z face
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
          zsqrt = sqrt( BFace**2 + 4*dInterface*AFace )
          z1    = (-BFace + zsqrt )/( 2*AFace )
          z2    = (-BFace - zsqrt )/( 2*AFace )

          if ( ( z1 .gt. 0d0 ) .and. ( z2 .gt. 0d0 ) ) then
              dtxyz(3) = min( z1**2, z2**2 )
          else if ( z1 .gt. 0d0 ) then
              dtxyz(3) = z1**2
          else if ( z2 .gt. 0d0 ) then
              dtxyz(3) = z2**2
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
      imindt = minloc( dtxyz, dim=1, mask=(dtxyz > 0) )
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
      dInterface = 0d0

      ! Local copies of cell size 
      dx = this%SubCellData%DX
      dy = this%SubCellData%DY
      dz = this%SubCellData%DZ

      !Reset t, dt will be replaced
      t = t - dt

      ! Leaving through x face
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
          if ( dInterface .eq. 0.0 ) then
              exitFace = exitFaceX
              dt = .0
              return
          end if

          ! Solve
          call this%NewtonRaphsonTimeStep( dt, vx,              &
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

          ! Solve
          call this%NewtonRaphsonTimeStep( dt, vy,              &
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

          ! Solve
          call this%NewtonRaphsonTimeStep( dt, vz,              &
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
      imindt = minloc( dtxyz, dim=1, mask=(dtxyz > 0) )
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
      ! How to determine a proper tolerance ?
      ! Convergence properties display 
      ! problem dependence
      !doubleprecision :: nrtol = 1.0d-3 
      doubleprecision :: nrtol
      integer :: countIter
      integer :: maxIter = 100
      !----------------------------------------------------------------

      ! Initialize
      countIter = 0
      ! Force initial guess smaller
      ! than original value
      dt0       = 0.01*dt 
      nrf0      = 0d0
      nrfprim   = 0d0
      nrerror   = 1d6 ! Something big

      ! Define nrtol from current dt
      nrtol    = 0.1*dt

      print *, '*** TrackSubCell: NewtonRapshonEXPONENTIAL: original DT', dt

      ! Iteration until convergence or maxIterations
      do while( ( abs(nrerror/dt) .gt. 0.01 ) .and. ( countIter .lt. maxIter ) )
      !do while( ( abs(nrerror) .gt. nrtol ) .and. ( countIter .lt. maxIter ) )

          countIter = countIter + 1

          ! Compute displacement, although initially this 
          ! should be known
          dvdx = ( v2 - v1 )/dx
          if ( ( abs(dvdx) .gt. dvtol ) ) then
              dAdv = v*( exp(dvdx*dt0) - 1.0d0 )/dvdx
          else
              dAdv = v*dt
          end if

          ! RWPT displacement with initial guess
          nrf0  = dAdv + divD*dt0 + dB*sqrt( dt0 ) - dInterface

          ! Analytical derivative of RWPT displacement
          nrfprim = v*exp(dvdx*dt0) + divD + 0.5*dB/sqrt(dt0) 

          ! NR error and new time step
          nrerror = -nrf0/nrfprim
          dt0     = dt0 + nrerror 

          ! It can't be smaller than zero
          if ( dt0 .lt. 0d0 ) dt0 = -0.5*dt0 

          !print *, nrerror, dt0
          

      end do

      ! Assign return value
      dtnr = dt0

      ! If new value higher than the previous, return zero 
      if ( dt0 .gt. dt) then
          print *, '*** TrackSubCell: Inconsistency in Newton Raphson, ', &
              'new dt higher than previous, dt, dt0, nrerror', dt, dt0, nrerror 
          dtnr = 0d0 ! Is this a proper exit condition ?
      end if 

      ! If no convergence, return zero
      if ( ( countIter .eq. maxIter ) .and. ( abs(nrerror/dt) .gt. 0.01 ) ) then
      !if ( ( countIter .eq. maxIter ) .and. ( abs(nrerror) .gt. nrtol ) ) then
          print *, '*** TrackSubCell: No convergence Newton Raphson dt, dt0, nrerror', dt, dt0, nrerror 
          dtnr = 0d0 ! Is this a proper exit condition ?
      end if


  end subroutine NewtonRaphsonTimeStepExponentialAdvection


  ! RWPT
  subroutine pr_DispersionDivergence( this, x, y, z, alphaL, alphaT, Dmol, divDx, divDy, divDz )
      !----------------------------------------------------------------
      ! Compute dispersion divergence terms 
      ! 
      ! Params:
      !     - x, y, z             : local cell coordinates
      !     - alphaL              : longidutinal dispersivity
      !     - alphaT              : transverse dispersivity
      !     - Dmol                : molecular diffusion
      !     - divDx, divDy, divDz : dispersion divergence, output 
      !----------------------------------------------------------------
      ! Specifications
      !----------------------------------------------------------------
      implicit none
      class (TrackSubCellType)    :: this
      ! input
      doubleprecision, intent(in) :: x, y, z
      doubleprecision, intent(in) :: alphaL, alphaT, Dmol
      ! output
      doubleprecision, intent(inout) :: divDx, divDy, divDz
      ! local
      doubleprecision, dimension(4) :: v000
      doubleprecision, dimension(4) :: v100
      doubleprecision, dimension(4) :: v010
      doubleprecision, dimension(4) :: v110
      doubleprecision, dimension(4) :: v001
      doubleprecision, dimension(4) :: v101
      doubleprecision, dimension(4) :: v011
      doubleprecision, dimension(4) :: v111
      doubleprecision :: dDxxdx, dDxydy, dDxzdz, &
                         dDxydx, dDyydy, dDyzdz, &
                         dDxzdx, dDyzdy, dDzzdz
      !---------------------------------------------------------------- 

      ! Initialize
      divDx = 0d0
      divDy = 0d0
      divDz = 0d0

      ! Local copies of corner velocities
      ! pointers ?
      v000 = this%vCorner000 
      v100 = this%vCorner100
      v010 = this%vCorner010
      v110 = this%vCorner110
      v001 = this%vCorner001
      v101 = this%vCorner101
      v011 = this%vCorner011
      v111 = this%vCorner111
    
      ! Direction, coordinates, corner values
      call this%TrilinearDerivative( 1, x, y, z, &
                ( alphaT*v000(4) + Dmol ) + ( alphaL - alphaT )*v000(1)**2/v000(4), & 
                ( alphaT*v100(4) + Dmol ) + ( alphaL - alphaT )*v100(1)**2/v100(4), &
                ( alphaT*v010(4) + Dmol ) + ( alphaL - alphaT )*v010(1)**2/v010(4), &
                ( alphaT*v110(4) + Dmol ) + ( alphaL - alphaT )*v110(1)**2/v110(4), &
                ( alphaT*v001(4) + Dmol ) + ( alphaL - alphaT )*v001(1)**2/v001(4), &
                ( alphaT*v101(4) + Dmol ) + ( alphaL - alphaT )*v101(1)**2/v101(4), &
                ( alphaT*v011(4) + Dmol ) + ( alphaL - alphaT )*v011(1)**2/v011(4), &
                ( alphaT*v111(4) + Dmol ) + ( alphaL - alphaT )*v111(1)**2/v111(4), &
                dDxxdx )
      call this%TrilinearDerivative( 2, x, y, z, &
                ( alphaT*v000(4) + Dmol ) + ( alphaL - alphaT )*v000(2)**2/v000(4), &
                ( alphaT*v100(4) + Dmol ) + ( alphaL - alphaT )*v100(2)**2/v100(4), &
                ( alphaT*v010(4) + Dmol ) + ( alphaL - alphaT )*v010(2)**2/v010(4), &
                ( alphaT*v110(4) + Dmol ) + ( alphaL - alphaT )*v110(2)**2/v110(4), &
                ( alphaT*v001(4) + Dmol ) + ( alphaL - alphaT )*v001(2)**2/v001(4), &
                ( alphaT*v101(4) + Dmol ) + ( alphaL - alphaT )*v101(2)**2/v101(4), &
                ( alphaT*v011(4) + Dmol ) + ( alphaL - alphaT )*v011(2)**2/v011(4), &
                ( alphaT*v111(4) + Dmol ) + ( alphaL - alphaT )*v111(2)**2/v111(4), &
                dDyydy )
      call this%TrilinearDerivative( 3, x, y, z, &
                ( alphaT*v000(4) + Dmol ) + ( alphaL - alphaT )*v000(3)**2/v000(4), &
                ( alphaT*v100(4) + Dmol ) + ( alphaL - alphaT )*v100(3)**2/v100(4), &
                ( alphaT*v010(4) + Dmol ) + ( alphaL - alphaT )*v010(3)**2/v010(4), &
                ( alphaT*v110(4) + Dmol ) + ( alphaL - alphaT )*v110(3)**2/v110(4), &
                ( alphaT*v001(4) + Dmol ) + ( alphaL - alphaT )*v001(3)**2/v001(4), &
                ( alphaT*v101(4) + Dmol ) + ( alphaL - alphaT )*v101(3)**2/v101(4), &
                ( alphaT*v011(4) + Dmol ) + ( alphaL - alphaT )*v011(3)**2/v011(4), &
                ( alphaT*v111(4) + Dmol ) + ( alphaL - alphaT )*v111(3)**2/v111(4), &
                dDzzdz )
      call this%TrilinearDerivative( 1, x, y, z, & 
                ( alphaL - alphaT )*v000(1)*v000(2)/v000(4), & 
                ( alphaL - alphaT )*v100(1)*v100(2)/v100(4), &
                ( alphaL - alphaT )*v010(1)*v010(2)/v010(4), &
                ( alphaL - alphaT )*v110(1)*v110(2)/v110(4), &
                ( alphaL - alphaT )*v001(1)*v001(2)/v001(4), &
                ( alphaL - alphaT )*v101(1)*v101(2)/v101(4), &
                ( alphaL - alphaT )*v011(1)*v011(2)/v011(4), &
                ( alphaL - alphaT )*v111(1)*v111(2)/v111(4), &
                dDxydx )
      call this%TrilinearDerivative( 1, x, y, z, &
                ( alphaL - alphaT )*v000(1)*v000(3)/v000(4), &
                ( alphaL - alphaT )*v100(1)*v100(3)/v100(4), &
                ( alphaL - alphaT )*v010(1)*v010(3)/v010(4), &
                ( alphaL - alphaT )*v110(1)*v110(3)/v110(4), &
                ( alphaL - alphaT )*v001(1)*v001(3)/v001(4), &
                ( alphaL - alphaT )*v101(1)*v101(3)/v101(4), &
                ( alphaL - alphaT )*v011(1)*v011(3)/v011(4), &
                ( alphaL - alphaT )*v111(1)*v111(3)/v111(4), &
                dDxzdx )
      call this%TrilinearDerivative( 2, x, y, z, &
                ( alphaL - alphaT )*v000(1)*v000(2)/v000(4), &
                ( alphaL - alphaT )*v100(1)*v100(2)/v100(4), &
                ( alphaL - alphaT )*v010(1)*v010(2)/v010(4), &
                ( alphaL - alphaT )*v110(1)*v110(2)/v110(4), &
                ( alphaL - alphaT )*v001(1)*v001(2)/v001(4), &
                ( alphaL - alphaT )*v101(1)*v101(2)/v101(4), &
                ( alphaL - alphaT )*v011(1)*v011(2)/v011(4), &
                ( alphaL - alphaT )*v111(1)*v111(2)/v111(4), &
                dDxydy )
      call this%TrilinearDerivative( 2, x, y, z, &
                ( alphaL - alphaT )*v000(2)*v000(3)/v000(4), &
                ( alphaL - alphaT )*v100(2)*v100(3)/v100(4), &
                ( alphaL - alphaT )*v010(2)*v010(3)/v010(4), &
                ( alphaL - alphaT )*v110(2)*v110(3)/v110(4), &
                ( alphaL - alphaT )*v001(2)*v001(3)/v001(4), &
                ( alphaL - alphaT )*v101(2)*v101(3)/v101(4), &
                ( alphaL - alphaT )*v011(2)*v011(3)/v011(4), &
                ( alphaL - alphaT )*v111(2)*v111(3)/v111(4), &
                dDyzdy )
      call this%TrilinearDerivative( 3, x, y, z, &
                ( alphaL - alphaT )*v000(1)*v000(3)/v000(4), &
                ( alphaL - alphaT )*v100(1)*v100(3)/v100(4), &
                ( alphaL - alphaT )*v010(1)*v010(3)/v010(4), &
                ( alphaL - alphaT )*v110(1)*v110(3)/v110(4), &
                ( alphaL - alphaT )*v001(1)*v001(3)/v001(4), &
                ( alphaL - alphaT )*v101(1)*v101(3)/v101(4), &
                ( alphaL - alphaT )*v011(1)*v011(3)/v011(4), &
                ( alphaL - alphaT )*v111(1)*v111(3)/v111(4), &
                dDxzdz )
      call this%TrilinearDerivative( 3, x, y, z, &
                ( alphaL - alphaT )*v000(2)*v000(3)/v000(4), &
                ( alphaL - alphaT )*v100(2)*v100(3)/v100(4), &
                ( alphaL - alphaT )*v010(2)*v010(3)/v010(4), &
                ( alphaL - alphaT )*v110(2)*v110(3)/v110(4), &
                ( alphaL - alphaT )*v001(2)*v001(3)/v001(4), &
                ( alphaL - alphaT )*v101(2)*v101(3)/v101(4), &
                ( alphaL - alphaT )*v011(2)*v011(3)/v011(4), &
                ( alphaL - alphaT )*v111(2)*v111(3)/v111(4), &
                dDyzdz )

      divDx = dDxxdx + dDxydy + dDxzdz
      divDy = dDxydx + dDyydy + dDyzdz
      divDz = dDxzdx + dDyzdy + dDzzdz


  end subroutine pr_DispersionDivergence


  ! RWPT
  subroutine pr_DispersionDivergenceDischarge( this, x, y, z, alphaL, alphaT, Dmol, divDx, divDy, divDz )
      !----------------------------------------------------------------
      ! Compute dispersion divergence terms 
      ! 
      ! Params:
      !     - x, y, z             : local cell coordinates
      !     - alphaL              : longidutinal dispersivity
      !     - alphaT              : transverse dispersivity
      !     - Dmol                : molecular diffusion
      !     - divDx, divDy, divDz : dispersion divergence, output 
      !----------------------------------------------------------------
      ! Specifications
      !----------------------------------------------------------------
      implicit none
      class (TrackSubCellType)    :: this
      ! input
      doubleprecision, intent(in) :: x, y, z
      doubleprecision, intent(in) :: alphaL, alphaT, Dmol
      ! output
      doubleprecision, intent(inout) :: divDx, divDy, divDz
      ! local
      doubleprecision, dimension(4) :: v000
      doubleprecision, dimension(4) :: v100
      doubleprecision, dimension(4) :: v010
      doubleprecision, dimension(4) :: v110
      doubleprecision, dimension(4) :: v001
      doubleprecision, dimension(4) :: v101
      doubleprecision, dimension(4) :: v011
      doubleprecision, dimension(4) :: v111

      doubleprecision :: p000
      doubleprecision :: p100
      doubleprecision :: p010
      doubleprecision :: p110
      doubleprecision :: p001
      doubleprecision :: p101
      doubleprecision :: p011
      doubleprecision :: p111

      doubleprecision :: dDxxdx, dDxydy, dDxzdz, &
                         dDxydx, dDyydy, dDyzdz, &
                         dDxzdx, dDyzdy, dDzzdz
      !---------------------------------------------------------------- 

      ! Initialize
      divDx = 0d0
      divDy = 0d0
      divDz = 0d0

      ! Local copies of SPECIFIC DISCHARGE
      v000 = this%qCorner000 
      v100 = this%qCorner100
      v010 = this%qCorner010
      v110 = this%qCorner110
      v001 = this%qCorner001
      v101 = this%qCorner101
      v011 = this%qCorner011
      v111 = this%qCorner111
  
      ! Porosities
      p000 = this%porosity000 
      p100 = this%porosity100
      p010 = this%porosity010
      p110 = this%porosity110
      p001 = this%porosity001
      p101 = this%porosity101
      p011 = this%porosity011
      p111 = this%porosity111


      ! Direction, coordinates, corner values
      call this%TrilinearDerivative( 1, x, y, z, &
                ( alphaT*v000(4) + p000*Dmol ) + ( alphaL - alphaT )*v000(1)**2/v000(4), & 
                ( alphaT*v100(4) + p100*Dmol ) + ( alphaL - alphaT )*v100(1)**2/v100(4), &
                ( alphaT*v010(4) + p010*Dmol ) + ( alphaL - alphaT )*v010(1)**2/v010(4), &
                ( alphaT*v110(4) + p110*Dmol ) + ( alphaL - alphaT )*v110(1)**2/v110(4), &
                ( alphaT*v001(4) + p001*Dmol ) + ( alphaL - alphaT )*v001(1)**2/v001(4), &
                ( alphaT*v101(4) + p101*Dmol ) + ( alphaL - alphaT )*v101(1)**2/v101(4), &
                ( alphaT*v011(4) + p011*Dmol ) + ( alphaL - alphaT )*v011(1)**2/v011(4), &
                ( alphaT*v111(4) + p111*Dmol ) + ( alphaL - alphaT )*v111(1)**2/v111(4), &
                dDxxdx )
      call this%TrilinearDerivative( 2, x, y, z, &
                ( alphaT*v000(4) + p000*Dmol ) + ( alphaL - alphaT )*v000(2)**2/v000(4), &
                ( alphaT*v100(4) + p100*Dmol ) + ( alphaL - alphaT )*v100(2)**2/v100(4), &
                ( alphaT*v010(4) + p010*Dmol ) + ( alphaL - alphaT )*v010(2)**2/v010(4), &
                ( alphaT*v110(4) + p110*Dmol ) + ( alphaL - alphaT )*v110(2)**2/v110(4), &
                ( alphaT*v001(4) + p001*Dmol ) + ( alphaL - alphaT )*v001(2)**2/v001(4), &
                ( alphaT*v101(4) + p101*Dmol ) + ( alphaL - alphaT )*v101(2)**2/v101(4), &
                ( alphaT*v011(4) + p011*Dmol ) + ( alphaL - alphaT )*v011(2)**2/v011(4), &
                ( alphaT*v111(4) + p111*Dmol ) + ( alphaL - alphaT )*v111(2)**2/v111(4), &
                dDyydy )
      call this%TrilinearDerivative( 3, x, y, z, &
                ( alphaT*v000(4) + p000*Dmol ) + ( alphaL - alphaT )*v000(3)**2/v000(4), &
                ( alphaT*v100(4) + p100*Dmol ) + ( alphaL - alphaT )*v100(3)**2/v100(4), &
                ( alphaT*v010(4) + p010*Dmol ) + ( alphaL - alphaT )*v010(3)**2/v010(4), &
                ( alphaT*v110(4) + p110*Dmol ) + ( alphaL - alphaT )*v110(3)**2/v110(4), &
                ( alphaT*v001(4) + p001*Dmol ) + ( alphaL - alphaT )*v001(3)**2/v001(4), &
                ( alphaT*v101(4) + p101*Dmol ) + ( alphaL - alphaT )*v101(3)**2/v101(4), &
                ( alphaT*v011(4) + p011*Dmol ) + ( alphaL - alphaT )*v011(3)**2/v011(4), &
                ( alphaT*v111(4) + p111*Dmol ) + ( alphaL - alphaT )*v111(3)**2/v111(4), &
                dDzzdz )
      call this%TrilinearDerivative( 1, x, y, z, & 
                ( alphaL - alphaT )*v000(1)*v000(2)/v000(4), & 
                ( alphaL - alphaT )*v100(1)*v100(2)/v100(4), &
                ( alphaL - alphaT )*v010(1)*v010(2)/v010(4), &
                ( alphaL - alphaT )*v110(1)*v110(2)/v110(4), &
                ( alphaL - alphaT )*v001(1)*v001(2)/v001(4), &
                ( alphaL - alphaT )*v101(1)*v101(2)/v101(4), &
                ( alphaL - alphaT )*v011(1)*v011(2)/v011(4), &
                ( alphaL - alphaT )*v111(1)*v111(2)/v111(4), &
                dDxydx )
      call this%TrilinearDerivative( 1, x, y, z, &
                ( alphaL - alphaT )*v000(1)*v000(3)/v000(4), &
                ( alphaL - alphaT )*v100(1)*v100(3)/v100(4), &
                ( alphaL - alphaT )*v010(1)*v010(3)/v010(4), &
                ( alphaL - alphaT )*v110(1)*v110(3)/v110(4), &
                ( alphaL - alphaT )*v001(1)*v001(3)/v001(4), &
                ( alphaL - alphaT )*v101(1)*v101(3)/v101(4), &
                ( alphaL - alphaT )*v011(1)*v011(3)/v011(4), &
                ( alphaL - alphaT )*v111(1)*v111(3)/v111(4), &
                dDxzdx )
      call this%TrilinearDerivative( 2, x, y, z, &
                ( alphaL - alphaT )*v000(1)*v000(2)/v000(4), &
                ( alphaL - alphaT )*v100(1)*v100(2)/v100(4), &
                ( alphaL - alphaT )*v010(1)*v010(2)/v010(4), &
                ( alphaL - alphaT )*v110(1)*v110(2)/v110(4), &
                ( alphaL - alphaT )*v001(1)*v001(2)/v001(4), &
                ( alphaL - alphaT )*v101(1)*v101(2)/v101(4), &
                ( alphaL - alphaT )*v011(1)*v011(2)/v011(4), &
                ( alphaL - alphaT )*v111(1)*v111(2)/v111(4), &
                dDxydy )
      call this%TrilinearDerivative( 2, x, y, z, &
                ( alphaL - alphaT )*v000(2)*v000(3)/v000(4), &
                ( alphaL - alphaT )*v100(2)*v100(3)/v100(4), &
                ( alphaL - alphaT )*v010(2)*v010(3)/v010(4), &
                ( alphaL - alphaT )*v110(2)*v110(3)/v110(4), &
                ( alphaL - alphaT )*v001(2)*v001(3)/v001(4), &
                ( alphaL - alphaT )*v101(2)*v101(3)/v101(4), &
                ( alphaL - alphaT )*v011(2)*v011(3)/v011(4), &
                ( alphaL - alphaT )*v111(2)*v111(3)/v111(4), &
                dDyzdy )
      call this%TrilinearDerivative( 3, x, y, z, &
                ( alphaL - alphaT )*v000(1)*v000(3)/v000(4), &
                ( alphaL - alphaT )*v100(1)*v100(3)/v100(4), &
                ( alphaL - alphaT )*v010(1)*v010(3)/v010(4), &
                ( alphaL - alphaT )*v110(1)*v110(3)/v110(4), &
                ( alphaL - alphaT )*v001(1)*v001(3)/v001(4), &
                ( alphaL - alphaT )*v101(1)*v101(3)/v101(4), &
                ( alphaL - alphaT )*v011(1)*v011(3)/v011(4), &
                ( alphaL - alphaT )*v111(1)*v111(3)/v111(4), &
                dDxzdz )
      call this%TrilinearDerivative( 3, x, y, z, &
                ( alphaL - alphaT )*v000(2)*v000(3)/v000(4), &
                ( alphaL - alphaT )*v100(2)*v100(3)/v100(4), &
                ( alphaL - alphaT )*v010(2)*v010(3)/v010(4), &
                ( alphaL - alphaT )*v110(2)*v110(3)/v110(4), &
                ( alphaL - alphaT )*v001(2)*v001(3)/v001(4), &
                ( alphaL - alphaT )*v101(2)*v101(3)/v101(4), &
                ( alphaL - alphaT )*v011(2)*v011(3)/v011(4), &
                ( alphaL - alphaT )*v111(2)*v111(3)/v111(4), &
                dDyzdz )

      divDx = ( dDxxdx + dDxydy + dDxzdz )/this%SubCellData%Porosity/this%SubCellData%Retardation
      divDy = ( dDxydx + dDyydy + dDyzdz )/this%SubCellData%Porosity/this%SubCellData%Retardation
      divDz = ( dDxzdx + dDyzdy + dDzzdz )/this%SubCellData%Porosity/this%SubCellData%Retardation


  end subroutine pr_DispersionDivergenceDischarge


  ! RWPT
  subroutine pr_DisplacementRandom( this, x, y, z, alphaL, alphaT, Dmol, dBx, dBy, dBz ) 
      !----------------------------------------------------------------
      ! Computes the product between displacement matrix and random 
      ! vector
      !
      ! Params:
      !     - x, y, z       : local cell coordinates
      !     - alphaL        : longidutinal dispersivity
      !     - alphaT        : transverse dispersivity
      !     - Dmol          : molecular diffusion
      !     - dBx, dBy, dBz : random dispersion displacement, output
      !----------------------------------------------------------------
      ! Specifications
      !----------------------------------------------------------------
      implicit none
      class (TrackSubCellType) :: this
      ! input
      doubleprecision, intent(in)    :: x, y, z
      doubleprecision, intent(in)    :: alphaL, alphaT, Dmol
      ! output
      doubleprecision, intent(inout) :: dBx, dBy, dBz
      ! local
      doubleprecision :: vBx, vBy, vBz, vBnorm, vBnormxy
      doubleprecision :: B11, B12, B13, B21, B22, B23, B31, B32
      doubleprecision :: rdmx, rdmy, rdmz
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
      v000 = this%vCorner000 
      v100 = this%vCorner100
      v010 = this%vCorner010
      v110 = this%vCorner110
      v001 = this%vCorner001
      v101 = this%vCorner101
      v011 = this%vCorner011
      v111 = this%vCorner111

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
      ! Refs: Fernàndez-García et al. 2005; Salamon et al. 2006
      ! Requires some kind of handling for the case 
      ! of zero vBnorm
      B11 =       vBx*sqrt( 2*( alphaL*vBnorm + Dmol ) )/vBnorm
      B12 =  -vBx*vBz*sqrt( 2*( alphaT*vBnorm + Dmol ) )/vBnorm/vBnormxy
      B13 =      -vBy*sqrt( 2*( alphaT*vBnorm + Dmol ) )/vBnormxy
      B21 =       vBy*sqrt( 2*( alphaL*vBnorm + Dmol ) )/vBnorm
      B22 =  -vBy*vBz*sqrt( 2*( alphaT*vBnorm + Dmol ) )/vBnorm/vBnormxy
      B23 =       vBx*sqrt( 2*( alphaT*vBnorm + Dmol ) )/vBnormxy
      B31 =       vBz*sqrt( 2*( alphaL*vBnorm + Dmol ) )/vBnorm
      B32 =  vBnormxy*sqrt( 2*( alphaT*vBnorm + Dmol ) )/vBnorm
   

      ! Compute random numbers
      call this%GenerateStandardNormalRandom( rdmx ) 
      call this%GenerateStandardNormalRandom( rdmy ) 
      call this%GenerateStandardNormalRandom( rdmz ) 


      ! Compute displacement times random
      dBx = B11*rdmx + B12*rdmy + B13*rdmz 
      dBy = B21*rdmx + B22*rdmy + B23*rdmz 
      dBz = B31*rdmx + B32*rdmy 


  end subroutine pr_DisplacementRandom


  subroutine pr_DisplacementRandomDischarge( this, x, y, z, alphaL, alphaT, Dmol, dBx, dBy, dBz ) 
      !----------------------------------------------------------------
      ! Computes the product between displacement matrix and random 
      ! vector
      !
      ! Params:
      !     - x, y, z       : local cell coordinates
      !     - alphaL        : longidutinal dispersivity
      !     - alphaT        : transverse dispersivity
      !     - Dmol          : molecular diffusion
      !     - dBx, dBy, dBz : random dispersion displacement, output
      !----------------------------------------------------------------
      ! Specifications
      !----------------------------------------------------------------
      implicit none
      class (TrackSubCellType) :: this
      ! input
      doubleprecision, intent(in)    :: x, y, z
      doubleprecision, intent(in)    :: alphaL, alphaT, Dmol
      ! output
      doubleprecision, intent(inout) :: dBx, dBy, dBz
      ! local
      doubleprecision :: vBx, vBy, vBz, vBnorm, vBnormxy
      doubleprecision :: B11, B12, B13, B21, B22, B23, B31, B32
      doubleprecision :: rdmx, rdmy, rdmz
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
      v000 = this%qCorner000 / this%porosity000 / this%SubCellData%Retardation  
      v100 = this%qCorner100 / this%porosity100 / this%SubCellData%Retardation
      v010 = this%qCorner010 / this%porosity010 / this%SubCellData%Retardation
      v110 = this%qCorner110 / this%porosity110 / this%SubCellData%Retardation
      v001 = this%qCorner001 / this%porosity001 / this%SubCellData%Retardation
      v101 = this%qCorner101 / this%porosity101 / this%SubCellData%Retardation
      v011 = this%qCorner011 / this%porosity011 / this%SubCellData%Retardation
      v111 = this%qCorner111 / this%porosity111 / this%SubCellData%Retardation

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
      ! Refs: Fernàndez-García et al. 2005; Salamon et al. 2006
      ! Requires some kind of handling for the case 
      ! of zero vBnorm
      B11 =       vBx*sqrt( 2*( alphaL*vBnorm + Dmol ) )/vBnorm
      B12 =  -vBx*vBz*sqrt( 2*( alphaT*vBnorm + Dmol ) )/vBnorm/vBnormxy
      B13 =      -vBy*sqrt( 2*( alphaT*vBnorm + Dmol ) )/vBnormxy
      B21 =       vBy*sqrt( 2*( alphaL*vBnorm + Dmol ) )/vBnorm
      B22 =  -vBy*vBz*sqrt( 2*( alphaT*vBnorm + Dmol ) )/vBnorm/vBnormxy
      B23 =       vBx*sqrt( 2*( alphaT*vBnorm + Dmol ) )/vBnormxy
      B31 =       vBz*sqrt( 2*( alphaL*vBnorm + Dmol ) )/vBnorm
      B32 =  vBnormxy*sqrt( 2*( alphaT*vBnorm + Dmol ) )/vBnorm
   

      ! Compute random numbers
      call this%GenerateStandardNormalRandom( rdmx ) 
      call this%GenerateStandardNormalRandom( rdmy ) 
      call this%GenerateStandardNormalRandom( rdmz ) 


      ! Compute displacement times random
      dBx = B11*rdmx + B12*rdmy + B13*rdmz 
      dBy = B21*rdmx + B22*rdmy + B23*rdmz 
      dBz = B31*rdmx + B32*rdmy 


  end subroutine pr_DisplacementRandomDischarge




  ! RWPT
  subroutine pr_GenerateStandardNormalRandom( this, random_value )
      !----------------------------------------------------------------
      ! Generate a random number from an standard normal distribution
      ! Inherited from RW3D:library_gslib:random_normal
      !
      !----------------------------------------------------------------
      ! Specifications
      !----------------------------------------------------------------
      implicit none
      class(TrackSubCellType) :: this 
      ! input/output
      doubleprecision, intent(inout) :: random_value
      ! local
      doubleprecision :: harvest(12)
      !----------------------------------------------------------------

      call random_number (harvest)
      random_value = sum(harvest)-6.d0

  end subroutine


  ! RWPT
  ! DEPRECATION WARNING
  subroutine pr_ComputeCornerVelocities( this, neighborSubCellData )
      !----------------------------------------------------------------
      ! From its subCellData and neighborSubCellData array, 
      ! computes velocities at cell corners
      !
      !----------------------------------------------------------------
      ! Specifications
      !----------------------------------------------------------------
      implicit none
      class (TrackSubCellType) :: this
      type(ModpathSubCellDataType), dimension(18) :: neighborSubCellData
      logical :: twoDimensionsDomain
      real    :: nominalFactor
      !----------------------------------------------------------------

      ! Do something more elegant please

      ! Think of a counter or something for dealing with interfaces
      ! with missing connections, cells
      ! Remember indexation
      ! 1: connection 1
      ! 4: connection 2
      ! 7: connection 3
      ! 10: connection 4
      ! 13: connection 5
      ! 16: connection 6
      ! and so on...
      !
      ! The 0.25 factor is for the case
      ! in which there are four cells for computing 
      ! values at the corner, should be verified 
      ! and maybe computed by finding the valid cells

      ! In 2D, only corners with third
      ! index equal 0

      nominalFactor = 0.25 ! 3D
      !twoDimensionsDomain = .true.
      !if ( twoDimensionsDomain ) nominalFactor = 0.5 ! 2D

      this%vCorner000(1) = nominalFactor*( this%SubCellData%VX1 + neighborSubCellData(7)%VX1  + &
          neighborSubCellData(13)%VX1 + neighborSubCellData(8)%VX1 )
      this%vCorner100(1) = nominalFactor*( this%SubCellData%VX2 + neighborSubCellData(7)%VX2  + &
          neighborSubCellData(13)%VX2 + neighborSubCellData(8)%VX2 )
      this%vCorner010(1) = nominalFactor*( this%SubCellData%VX1 + neighborSubCellData(10)%VX1 + &
          neighborSubCellData(13)%VX1 + neighborSubCellData(11)%VX1 )
      this%vCorner110(1) = nominalFactor*( this%SubCellData%VX2 + neighborSubCellData(10)%VX2 + &
          neighborSubCellData(13)%VX2 + neighborSubCellData(11)%VX2 )
      this%vCorner001(1) = nominalFactor*( this%SubCellData%VX1 + neighborSubCellData(7)%VX1  + &
          neighborSubCellData(16)%VX1 + neighborSubCellData(9)%VX1 )
      this%vCorner101(1) = nominalFactor*( this%SubCellData%VX2 + neighborSubCellData(7)%VX2  + &
          neighborSubCellData(16)%VX2 + neighborSubCellData(9)%VX2 )
      this%vCorner011(1) = nominalFactor*( this%SubCellData%VX1 + neighborSubCellData(10)%VX1 + &
          neighborSubCellData(16)%VX1 + neighborSubCellData(12)%VX1 )
      this%vCorner111(1) = nominalFactor*( this%SubCellData%VX2 + neighborSubCellData(10)%VX2 + &
          neighborSubCellData(16)%VX2 + neighborSubCellData(12)%VX2 )


      this%vCorner000(2) = nominalFactor*( this%SubCellData%VY1 + neighborSubCellData(1)%VY1  + &
          neighborSubCellData(13)%VY1 + neighborSubCellData(14)%VY1 )
      this%vCorner010(2) = nominalFactor*( this%SubCellData%VY2 + neighborSubCellData(1)%VY2  + &
          neighborSubCellData(13)%VY2 + neighborSubCellData(14)%VY2 )
      this%vCorner100(2) = nominalFactor*( this%SubCellData%VY1 + neighborSubCellData(4)%VY1  + &
          neighborSubCellData(13)%VY1 + neighborSubCellData(15)%VY1 )
      this%vCorner110(2) = nominalFactor*( this%SubCellData%VY2 + neighborSubCellData(4)%VY2  + &
          neighborSubCellData(13)%VY2 + neighborSubCellData(15)%VY2 )
      this%vCorner001(2) = nominalFactor*( this%SubCellData%VY1 + neighborSubCellData(1)%VY1  + &
          neighborSubCellData(16)%VY1 + neighborSubCellData(17)%VY1 )
      this%vCorner011(2) = nominalFactor*( this%SubCellData%VY2 + neighborSubCellData(1)%VY2  + &
          neighborSubCellData(16)%VY2 + neighborSubCellData(17)%VY2 )
      this%vCorner101(2) = nominalFactor*( this%SubCellData%VY1 + neighborSubCellData(4)%VY1  + &
          neighborSubCellData(16)%VY1 + neighborSubCellData(18)%VY1 )
      this%vCorner111(2) = nominalFactor*( this%SubCellData%VY2 + neighborSubCellData(4)%VY2  + &
          neighborSubCellData(16)%VY2 + neighborSubCellData(18)%VY2 )


      this%vCorner000(3) = nominalFactor*( this%SubCellData%VZ1 + neighborSubCellData(1)%VZ1  + &
          neighborSubCellData(7)%VZ1 + neighborSubCellData(2)%VZ1 )
      this%vCorner001(3) = nominalFactor*( this%SubCellData%VZ2 + neighborSubCellData(1)%VZ2  + &
          neighborSubCellData(7)%VZ2 + neighborSubCellData(2)%VZ2 )
      this%vCorner100(3) = nominalFactor*( this%SubCellData%VZ1 + neighborSubCellData(4)%VZ1  + &
          neighborSubCellData(7)%VZ1 + neighborSubCellData(5)%VZ1 )
      this%vCorner101(3) = nominalFactor*( this%SubCellData%VZ2 + neighborSubCellData(4)%VZ2  + &
          neighborSubCellData(7)%VZ2 + neighborSubCellData(5)%VZ2 )
      this%vCorner010(3) = nominalFactor*( this%SubCellData%VZ1 + neighborSubCellData(1)%VZ1  + &
          neighborSubCellData(10)%VZ1 + neighborSubCellData(3)%VZ1 )
      this%vCorner011(3) = nominalFactor*( this%SubCellData%VZ2 + neighborSubCellData(1)%VZ2  + &
          neighborSubCellData(10)%VZ2 + neighborSubCellData(3)%VZ2 )
      this%vCorner110(3) = nominalFactor*( this%SubCellData%VZ1 + neighborSubCellData(4)%VZ1  + &
          neighborSubCellData(10)%VZ1 + neighborSubCellData(6)%VZ1 )
      this%vCorner111(3) = nominalFactor*( this%SubCellData%VZ2 + neighborSubCellData(4)%VZ2  + &
          neighborSubCellData(10)%VZ2 + neighborSubCellData(6)%VZ2 )


      this%vCorner000(4) = sqrt( this%vCorner000(1)**2 + this%vCorner000(2)**2 + this%vCorner000(3)**2 )
      this%vCorner100(4) = sqrt( this%vCorner100(1)**2 + this%vCorner100(2)**2 + this%vCorner100(3)**2 )
      this%vCorner010(4) = sqrt( this%vCorner010(1)**2 + this%vCorner010(2)**2 + this%vCorner010(3)**2 )
      this%vCorner110(4) = sqrt( this%vCorner110(1)**2 + this%vCorner110(2)**2 + this%vCorner110(3)**2 )
      this%vCorner001(4) = sqrt( this%vCorner001(1)**2 + this%vCorner001(2)**2 + this%vCorner001(3)**2 )
      this%vCorner101(4) = sqrt( this%vCorner101(1)**2 + this%vCorner101(2)**2 + this%vCorner101(3)**2 )
      this%vCorner011(4) = sqrt( this%vCorner011(1)**2 + this%vCorner011(2)**2 + this%vCorner011(3)**2 )
      this%vCorner111(4) = sqrt( this%vCorner111(1)**2 + this%vCorner111(2)**2 + this%vCorner111(3)**2 )
      

  end subroutine pr_ComputeCornerVelocities


  ! RWPT
  subroutine pr_ComputeCornerDischarge( this, currentCellData, neighborCellData )
      !----------------------------------------------------------------
      ! From its subCellData and neighborSubCellData array, 
      ! computes velocities at cell corners
      !
      !----------------------------------------------------------------
      ! Specifications
      !----------------------------------------------------------------
      implicit none
      class(TrackSubCellType) :: this
      type(ModpathCellDataType) :: currentCellData
      type(ModpathCellDataType), dimension(2,18) :: neighborCellData
      integer, dimension(18,3)         :: neighborSubCellIndexes   ! nbcell, subRow, subColumn
      doubleprecision, dimension(18,6) :: neighborSubCellFaceFlows ! nbcell, flowFaceNumber
      doubleprecision, dimension(18,3) :: neighborSubCellFaceAreas ! nbcell, faceDirection
      doubleprecision, dimension(6)    :: centerSubCellFaceFlows
      integer :: n, m
      ! deprecated
      logical :: twoDimensionsDomain
      real    :: flowContributionFactor
      integer :: cid
      doubleprecision :: dx, dy, dz, dz0, dz1
      doubleprecision :: areaFlowX0, areaFlowX1
      doubleprecision :: areaFlowY0, areaFlowY1
      doubleprecision :: areaFlowZ
      !----------------------------------------------------------------

    
      ! Get sub cell indexes for current sub cell location
      neighborSubCellIndexes = currentCellData%GetNeighborSubCellIndexes( &
                             this%SubCellData%Row, this%SubCellData%Column )

      ! Fill neighbor cells faceFlows
      call pr_FillNeighborSubCellFaceFlowsAreas( this, currentCellData, neighborCellData, &
                neighborSubCellIndexes, neighborSubCellFaceFlows, neighborSubCellFaceAreas )

      ! Fill faceFlows from current subcell
      call currentCellData%FillSubCellFaceFlowsBuffer( &
          this%SubCellData%Row, this%SubCellData%Column, centerSubCellFaceFlows )

      ! Set indexes for corner x,y,z components
      call this%SetCornerComponentsIndexes() 

      ! For each dimension
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


      print *, '**** TrackSubCell:currentCell number', currentCellData%CellNumber
      print *, '**** TrackSubCell:currentCell dx, dy, dz', currentCellData%dx, currentCellData%dy, &
          currentCellData%GetDZ()
      print *, '**** TrackSubCell:currenttracksubCell dx, dy, dz', this%SubCellData%dx, this%SubCellData%dy, &
          this%SubCellData%dz
      print *, '**** TrackSubCell:ComputeCornerDischarge: will print areas' 
        
      do n=1, 18
          print *, n, neighborSubCellFaceAreas( n, : ) 
      end do


      call exit(0)


  end subroutine pr_ComputeCornerDischarge



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



  ! RWPT
  function pr_GetInterpolatedCornerDischarge( this, centerSubCellFaceFlows, &
                        neighborSubCellFaceFlows, neighborSubCellFaceAreas, & 
                                                cornerIndex, faceDirection  ) result( qCorner )
  !---------------------------------------------------------------------
  ! Compute interpolated corner discharge as the sum of flow rates
  ! of contributing faces dividing by the sum of face areas
  !
  ! Initializes sums with contribution of center cell and 
  ! loop over remaining involved cells
  !---------------------------------------------------------------------
  class( TrackSubCellType ) :: this
  ! input
  doubleprecision, dimension(6)   , intent(in) :: centerSubCellFaceFlows 
  doubleprecision, dimension(18,6), intent(in) :: neighborSubCellFaceFlows
  doubleprecision, dimension(18,3), intent(in) :: neighborSubCellFaceAreas
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
                     neighborSubCellFaceAreas( this%cornerXComponentIndexes( cornerIndex, m ), 1 )
                 sumFlowRate = sumFlowRate + & 
                     neighborSubCellFaceFlows( this%cornerXComponentIndexes( cornerIndex, m ) , &
                                               this%cornerXComponentIndexes( cornerIndex, 4 ) )  
             end do
         case (2)
             sumArea     = this%SubCellData%DX*this%SubCellData%DZ
             sumFlowRate = centerSubCellFaceFlows( this%cornerYComponentIndexes( cornerIndex, 4 ) )
             do m = 1, 3
                 sumArea = sumArea + &
                     neighborSubCellFaceAreas( this%cornerYComponentIndexes( cornerIndex, m ), 2 )
                 sumFlowRate = sumFlowRate + & 
                     neighborSubCellFaceFlows( this%cornerYComponentIndexes( cornerIndex, m ) , &
                                               this%cornerYComponentIndexes( cornerIndex, 4 ) )  
             end do
         case (3)
             sumArea     = this%SubCellData%DX*this%SubCellData%DY
             sumFlowRate = centerSubCellFaceFlows( this%cornerZComponentIndexes( cornerIndex, 4 ) )
             do m = 1, 3
                 sumArea = sumArea + &
                     neighborSubCellFaceAreas( this%cornerZComponentIndexes( cornerIndex, m ), 3 )
                 sumFlowRate = sumFlowRate + & 
                     neighborSubCellFaceFlows( this%cornerZComponentIndexes( cornerIndex, m ) , &
                                               this%cornerZComponentIndexes( cornerIndex, 4 ) )  
             end do
     end select 
  
     ! flowContributionFactor is cancelled implicitly    
     qCorner = sumFlowRate/sumArea
  
  
  
  end function pr_GetInterpolatedCornerDischarge



  ! RWPT
  subroutine pr_ComputeCornerPorosity( this, currentCellData, neighborCellData )
      !----------------------------------------------------------------
      ! From its subCellData and neighborSubCellData array, 
      ! computes velocities at cell corners
      !
      !----------------------------------------------------------------
      ! Specifications
      !----------------------------------------------------------------
      implicit none
      class (TrackSubCellType) :: this
      !type(ModpathSubCellDataType), dimension(18) :: neighborSubCellData
      type(ModpathCellDataType) :: currentCellData
      type(ModpathCellDataType), dimension(18) :: neighborCellData
      integer, dimension(18,3) :: neighborSubCellIndexes ! nbcell, subRow, subColumn
      doubleprecision, dimension(18) :: neighborSubCellPorosity 
      doubleprecision, dimension(18) :: neighborSubCellVolume   
      doubleprecision, dimension(6)    :: centerSubCellFaceFlows
      ! 8 corners with information of nc1, nc2, nc3, fN
      integer, dimension(8,6)  :: cornerPorosityIndexes ! 6 Initialized subcells + self + missing = 8
      logical :: twoDimensionsDomain
      real    :: flowContributionFactor
      integer :: cid, scid
      doubleprecision :: dx, dy, dz, dz0, dz1
      doubleprecision :: areaFlowX0, areaFlowX1
      doubleprecision :: areaFlowY0, areaFlowY1
      doubleprecision :: areaFlowZ
      doubleprecision :: volumeCenter
      doubleprecision :: volumeLower 
      doubleprecision :: volumeUpper
      !-----------------------------------------------------------------


      !print *, '** TrackSubCell:ComputeCornerPorosity: entered, constant value of', this%SubCellData%Porosity



      ! This is the fallback for spatially 
      ! constant porosity
      this%porosity000 = this%SubCellData%Porosity   
      this%porosity100 = this%SubCellData%Porosity
      this%porosity010 = this%SubCellData%Porosity       
      this%porosity110 = this%SubCellData%Porosity
      this%porosity001 = this%SubCellData%Porosity
      this%porosity101 = this%SubCellData%Porosity       
      this%porosity011 = this%SubCellData%Porosity
      this%porosity111 = this%SubCellData%Porosity


      return 




      !!-------------------! THIS CODE IS FOR VARIABLE POROSITY CASE ! --------------------

      !! Following indexes should be computed only once 
      !! subRow and subColumn changes. If kept at this point, 
      !! then computation is performed several times for the same 
      !! subCell

      !! subRow and subColumn determine indexes from NeighborCellData to be used 
      !! in interpolation process
      !if ( this%SubCellData%Row .eq. 1 ) then
      !    ! If subcell from first subRow 
      !    if ( this%SubCellData%Column .eq. 1 ) then
      !        ! If subcell from first subColumn
      !        ! Review
      !        neighborSubCellIndexes(1,:)  = [ 1,1,2]
      !        neighborSubCellIndexes(2,:)  = [ 1,2,2]
      !        neighborSubCellIndexes(3,:)  = [ 3,2,2]
      !        neighborSubCellIndexes(4,:)  = [ 0,1,2]
      !        neighborSubCellIndexes(5,:)  = [ 0,2,2]
      !        neighborSubCellIndexes(6,:)  = [10,2,2]
      !        neighborSubCellIndexes(7,:)  = [ 0,2,1]
      !        neighborSubCellIndexes(8,:)  = [13,2,1]
      !        neighborSubCellIndexes(9,:)  = [16,2,1]
      !        neighborSubCellIndexes(10,:) = [10,2,1]
      !        neighborSubCellIndexes(11,:) = [11,2,1]
      !        neighborSubCellIndexes(12,:) = [12,2,1]
      !        neighborSubCellIndexes(13,:) = [13,1,1]
      !        neighborSubCellIndexes(14,:) = [14,1,2]
      !        neighborSubCellIndexes(15,:) = [13,1,2]
      !        neighborSubCellIndexes(16,:) = [16,1,1]
      !        neighborSubCellIndexes(17,:) = [17,1,2]
      !        neighborSubCellIndexes(18,:) = [16,1,2]
      !    else
      !        ! If subcell from second subColumn
      !        ! Review
      !        neighborSubCellIndexes(1,:)  = [ 0,1,1]
      !        neighborSubCellIndexes(2,:)  = [ 0,2,1]
      !        neighborSubCellIndexes(3,:)  = [10,2,1]
      !        neighborSubCellIndexes(4,:)  = [ 4,1,1]
      !        neighborSubCellIndexes(5,:)  = [ 4,2,1]
      !        neighborSubCellIndexes(6,:)  = [ 6,2,1]
      !        neighborSubCellIndexes(7,:)  = [ 0,2,2]
      !        neighborSubCellIndexes(8,:)  = [13,2,2]
      !        neighborSubCellIndexes(9,:)  = [16,2,2]
      !        neighborSubCellIndexes(10,:) = [10,2,2]
      !        neighborSubCellIndexes(11,:) = [11,2,2]
      !        neighborSubCellIndexes(12,:) = [12,2,2]
      !        neighborSubCellIndexes(13,:) = [13,1,2]
      !        neighborSubCellIndexes(14,:) = [13,1,1]
      !        neighborSubCellIndexes(15,:) = [15,1,1]
      !        neighborSubCellIndexes(16,:) = [16,1,2]
      !        neighborSubCellIndexes(17,:) = [16,1,1]
      !        neighborSubCellIndexes(18,:) = [18,1,1]
      !    end if        
      !else
      !    ! If subcell from second subRow 
      !    if ( this%SubCellData%Column .eq. 1 ) then 
      !        ! If subcell from first subColumn
      !        ! Review
      !        neighborSubCellIndexes(1,:)  = [ 1,2,2]
      !        neighborSubCellIndexes(2,:)  = [ 2,1,2]
      !        neighborSubCellIndexes(3,:)  = [ 1,1,2]
      !        neighborSubCellIndexes(4,:)  = [ 0,2,2]
      !        neighborSubCellIndexes(5,:)  = [ 7,1,2]
      !        neighborSubCellIndexes(6,:)  = [ 0,1,2]
      !        neighborSubCellIndexes(7,:)  = [ 7,1,1]
      !        neighborSubCellIndexes(8,:)  = [ 8,1,1]
      !        neighborSubCellIndexes(9,:)  = [ 9,1,1]
      !        neighborSubCellIndexes(10,:) = [ 0,1,1]
      !        neighborSubCellIndexes(11,:) = [13,1,1]
      !        neighborSubCellIndexes(12,:) = [16,1,1]
      !        neighborSubCellIndexes(13,:) = [13,2,1]
      !        neighborSubCellIndexes(14,:) = [14,2,2]
      !        neighborSubCellIndexes(15,:) = [13,2,2]
      !        neighborSubCellIndexes(16,:) = [16,2,1]
      !        neighborSubCellIndexes(17,:) = [17,2,2]
      !        neighborSubCellIndexes(18,:) = [16,2,2]
      !    else
      !        ! If subcell from second subColumn
      !        ! Review
      !        neighborSubCellIndexes(1,:)  = [ 0,2,1]
      !        neighborSubCellIndexes(2,:)  = [ 7,1,1]
      !        neighborSubCellIndexes(3,:)  = [ 0,1,1]
      !        neighborSubCellIndexes(4,:)  = [ 4,2,1]
      !        neighborSubCellIndexes(5,:)  = [ 5,1,1]
      !        neighborSubCellIndexes(6,:)  = [ 4,1,1]
      !        neighborSubCellIndexes(7,:)  = [ 7,1,2]
      !        neighborSubCellIndexes(8,:)  = [ 8,1,2]
      !        neighborSubCellIndexes(9,:)  = [ 9,1,2]
      !        neighborSubCellIndexes(10,:) = [ 0,1,2]
      !        neighborSubCellIndexes(11,:) = [13,1,2]
      !        neighborSubCellIndexes(12,:) = [16,1,2]
      !        neighborSubCellIndexes(13,:) = [13,2,2]
      !        neighborSubCellIndexes(14,:) = [13,2,1]
      !        neighborSubCellIndexes(15,:) = [15,2,1]
      !        neighborSubCellIndexes(16,:) = [16,2,2]
      !        neighborSubCellIndexes(17,:) = [16,2,1]
      !        neighborSubCellIndexes(18,:) = [18,2,1]
      !    end if        
      !end if 



      !! Compute cell sizes for areas computation
      !dx  = this%SubCellData%DX
      !dy  = this%SubCellData%DY
      !dz  = this%SubCellData%DZ
      !dz0 = neighborCellData(13)%GetDZ() ! query a lower layer cell and getdz 
      !dz1 = neighborCellData(16)%GetDZ() ! query an upper layer cell and getdz 


      !volumeCenter = dx*dy*dz
      !volumeLower  = dx*dy*dz0
      !volumeUpper  = dx*dy*dz1


      !do cid = 1,18

      !    ! Fill porosities 
      !    if ( neighborSubCellIndexes( cid, 1 ) .gt. 0 ) then
      !        ! A no center cell

      !        if ( neighborCellData( neighborSubCellIndexes( cid, 1 ) )%CellNumber .ne. -999 ) then
      !            ! the buffer is a normal cell
      !            neighborSubCellPorosity( cid ) = neighborCellData( neighborSubCellIndexes( cid,  1 ) )%Porosity
      !        else
      !            ! if the buffer is a custom cell
      !            ! Compute subcell index
      !            ! HERE WE HAD A COLUMNCOUNT
      ! 
      !     scid = ( neighborSubCellIndexes( cid, 2 ) - 1 )*2  + neighborSubCellIndexes( cid, 3 )

      !            ! Get porosity from subcell index. This
      !            ! should work fine in cases that 
      !            ! interpolation grid has nested subcells
      !            ! as required subcells for porosity are direct connections
      !            ! and for a smoothed usg, direct connections
      !            ! will be regular/real cells with defined porosity and not
      !            ! custom/linked cells. 
      !            neighborSubCellPorosity( cid ) = &
      !                neighborCellData( neighborSubCellIndexes( cid,  1 ) )%SubCellDataBuffer( scid )%Porosity
      !        end if
      !    else
      !        ! The center cell
      !        neighborSubCellPorosity( cid ) = currentCellData%Porosity
      !    end if

      !    
      !    ! Fill cell volumes 
      !    if ( neighborCellData( neighborSubCellIndexes( cid, 1 ))%Layer .eq. currentCellData%Layer ) then
      !        neighborSubCellVolume( cid ) = volumeCenter
      !    else if ( neighborCellData( neighborSubCellIndexes( cid, 1 ))%Layer .gt. currentCellData%Layer ) then
      !        ! Remember that deeper layers have higher index
      !        neighborSubCellVolume( cid ) = volumeLower
      !    else
      !        ! Is an upper layer
      !        neighborSubCellVolume( cid ) = volumeUpper
      !    end if

      !end do


      !! Convention for cell indexes of 
      !! coordinate components
      !! 1: 000
      !! 2: 100
      !! 3: 010
      !! 4: 110
      !! 5: 001
      !! 6: 101
      !! 7: 011
      !! 8: 111
      !cornerPorosityIndexes(1,:) = [ 1, 2,  7,  8, 13, 14 ] 
      !cornerPorosityIndexes(2,:) = [ 4, 5,  7,  8, 13, 15 ] 
      !cornerPorosityIndexes(3,:) = [ 1, 3, 10, 11, 13, 14 ] 
      !cornerPorosityIndexes(4,:) = [ 4, 6, 10, 12, 16, 18 ] 
      !cornerPorosityIndexes(5,:) = [ 1, 2,  7,  9, 16, 17 ] 
      !cornerPorosityIndexes(6,:) = [ 4, 5,  7,  9, 16, 18 ] 
      !cornerPorosityIndexes(7,:) = [ 1, 3, 10, 12, 16, 17 ] 
      !cornerPorosityIndexes(8,:) = [ 4, 6, 10, 11, 13, 15 ] 


      !! Now access required indexes for each corner
      !! In this step it is required initialization
      !! of porosity and volume of corner cells 
      !! not initialized in buffer because are 
      !! not required for flow purposes
      !! (related to volumeUpper or volumeLower according to corner).
      !! Remember that it is also required to consider 
      !! contribution of the center cell


      !! So the missing parameter it is only the porosity
      !! of the surrounding cube corners


      !! Functions for building cells 
      !! are defined at ParticleTrackingEngine.f90
      !! Thus building corner cells makes sense when 
      !! working at that file. 
      !! Rather than define the missing things at this file.







      !! Missing areas, should divide sum of face flows
      !!! Corner porosities follow similar indexation 
      !!! as corner discharge
      !!! Compute equivalent corner volume
      !!cornerVolumeZ0 = 1/8*( 4*dx*dy*dz + 4*dx*dy*dz0 )
      !!cornerVolumeZ1 = 1/8*( 4*dx*dy*dz + 4*dx*dy*dz1 )

      !!! Corner porosities follow similar indexation 
      !!! as corner discharge
      !!! Compute equivalent corner volume
      !!cornerVolumeZ0 = 1/8*( 4*dx*dy*dz + 4*dx*dy*dz0 )
      !!cornerVolumeZ1 = 1/8*( 4*dx*dy*dz + 4*dx*dy*dz1 )


      !!! Porosity is used as prototype of new variable structure
      !!this%cornerPorosity(0,0,0) = 

  end subroutine pr_ComputeCornerPorosity


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
      
      ! Initialize
      output = 0d0

      v00    = ( 1.0d0 - x )*v000 + x*v100
      v01    = ( 1.0d0 - x )*v001 + x*v101
      v10    = ( 1.0d0 - x )*v010 + x*v110
      v11    = ( 1.0d0 - x )*v011 + x*v111
      v0     = ( 1.0d0 - y )*v00  + y*v10
      v1     = ( 1.0d0 - y )*v01  + y*v11
      output = ( 1.0d0 - z )*v0   + z*v1

  end subroutine pr_Trilinear



  subroutine pr_FillNeighborSubCellFaceFlowsAreas( this, centerCellData, neighborCellData, &
                neighborSubCellIndexes, neighborSubCellFaceFlows, neighborSubCellFaceAreas )
      !-----------------------------------------------------------
      ! WHERE ?
      !-----------------------------------------------------------
      ! Specifications
      !-----------------------------------------------------------
      implicit none
      class (TrackSubCellType) :: this
      class (ModpathCellDataType), intent(in) :: centerCellData 
      class (ModpathCellDataType), dimension(2,18), intent(in) :: neighborCellData 
      integer, dimension(18,3), intent(in) ::  neighborSubCellIndexes ! nCellBuffer, subRow, subCol
      doubleprecision, dimension(18,6), intent(inout) :: neighborSubCellFaceFlows ! faceFlows
      doubleprecision, dimension(18,3), intent(inout) :: neighborSubCellFaceAreas ! faceAreas
      logical :: skipSubCells = .false.
      integer :: subConnectionIndex, neighborSubRow, neighborSubColumn
      integer :: n, m
      !-----------------------------------------------------------


      ! Reset arrays
      neighborSubCellFaceFlows = 0d0
      neighborSubCellFaceAreas = 0d0


      ! If current cell is not refined
      if ( centerCellData%GetSubCellCount() .eq. 1 ) then 
        ! DEBUG/DEV
        print *, 'CURRENT CELL IS NOT REFINED'
        ! Cell is connected to non refined cells, 
        ! it does not uses sub cells from itself.
        ! Neighbors without subcells
        ! Actually it does not need neighborSubCellIndexes
        do n = 1,18
            call neighborCellData( &
                1, neighborSubCellIndexes( n, 1 ) )%FillSubCellFaceFlowsBuffer( 1, 1, &
                                                     neighborSubCellFaceFlows( n, : ) )
            call neighborCellData( & 
                1, neighborSubCellIndexes( n, 1 ) )%FillSubCellFaceAreas( neighborSubCellFaceAreas( n, : ), skipSubCells ) 

            if ( neighborCellData( 1, neighborSubCellIndexes( n, 1 ) )%GetSubCellCount() .gt. 1 ) then 
                print *, n , '!!!!!!!!!! NEIGHBOR is REFINED !!!!!!!!!!!'
            else
                if ( neighborCellData( 1, neighborSubCellIndexes( n, 1 ) )%fromSubSubCell ) then 
                    print *, n, 'NEIGHBOR is GRANDCHILDREN OF  ', &
                        neighborCellData( 1, neighborSubCellIndexes( n, 1 ) )%parentCellNumber
                else if ( neighborCellData( 1, neighborSubCellIndexes( n, 1 ) )%fromSubCell ) then 
                    print *, n, 'NEIGHBOR is CHILDREN OF  ', &
                        neighborCellData( 1, neighborSubCellIndexes( n, 1 ) )%parentCellNumber
                else 
                    print *, n, 'NEIGHBOR is CELL NUMBER  ', &
                        neighborCellData( 1, neighborSubCellIndexes( n, 1 ) )%CellNumber

                end if


            end if
        end do

        ! Done 
        return

      end if  


      ! DEBUG/DEV
      print *, 'CURRENT CELL IS REFINED'

      ! DEBUG/DEV
      ! Print a report
      print *, 'TrackSubCell: will compute discharges'
      print *, '********************************************************************************************'
      print *, '** ComputeCornerDischarge:  neighbors information for TrackCell', centerCellData%CellNumber
      print *, '** ComputeCornerDischarge: subcell indexes',  this%SubCellData%Row, this%SubCellData%Column
      print *, '********************************************************************************************'

      do n = 1, 18
          print *, '-------------------------------------------------------------------'
          print *, '   -- Should fill flows with indexes: ', neighborSubCellIndexes( n, : )
          print *, '   -- NeighborCellBuffer ', neighborSubCellIndexes( n, 1 )
          skipSubCells = .false.

          ! Fill face flows with one of its own subcells
          if ( neighborSubCellIndexes( n, 1 ) .eq. 0 ) then
              print *, '   -- Will fill with subcell data from itself '
              call centerCellData%FillSubCellFaceFlowsBuffer( &
                                neighborSubCellIndexes( n, 2 ), neighborSubCellIndexes( n, 3 ), & 
                                                               neighborSubCellFaceFlows( n, : ) )

              call centerCellData%FillSubCellFaceAreas( neighborSubCellFaceAreas( n, : ), skipSubCells )

              cycle
          end if

          ! Fill face flows with data obtained from neighbors
          print *, '   -- Will fill with data from neighbor '
          if ( neighborCellData( 2, neighborSubCellIndexes( n, 1 ) )%CellNumber .gt. 0 ) then 
              print *, '   -- This is a DOUBLE BUFFER'
              print *, '   -- Based on requested face direction should determine INDEX'

        
              if ( & 
                  ( neighborCellData(1, neighborSubCellIndexes( n, 1 ) )%GetSubCellCount() .gt. 1 ) .or. &
                  ( neighborCellData(2, neighborSubCellIndexes( n, 1 ) )%GetSubCellCount() .gt. 1 ) ) then 
                  print *, 'IS REFINED HOW WHAT'
                  skipSubCells = .true.
              end if 

              ! By design, when the buffer is double, cells are not refined.
              ! Detect from which direction where requested.
              if ( neighborCellData( 2, neighborSubCellIndexes( n, 1 ) )%requestedFromDirection .eq. 1 ) then 
                  ! If x-direction, location in buffer is given by subcell row index 
                  subConnectionIndex = neighborSubCellIndexes( n, 2 )
                  neighborSubRow     = 1 
                  neighborSubColumn  = 1

                  print *, '   -- X DIRECTION THEN LOCATION IN BUFFER IS GIVEN BY SUBCELL ROW',&
                  neighborCellData( 2, neighborSubCellIndexes( n, 1 ) )%requestedFromDirection
                  !call neighborCellData(&
                  !    neighborSubCellIndexes( n, 2 ),&
                  !    neighborSubCellIndexes( n, 1 ) )%FillSubCellFaceFlowsBuffer( &
                  !                          1, 1, neighborSubCellFaceFlows( n, : ) )

                  print *, '   -- FILLING FROM CELLNUMBER',&
                  neighborCellData(&
                      neighborSubCellIndexes( n, 2 ), &
                      neighborSubCellIndexes( n, 1 ) )%CellNumber

              else if ( neighborCellData( 2, neighborSubCellIndexes( n, 1 ) )%requestedFromDirection .eq. 2 ) then
                  ! If y-direction, location in buffer is given by subcell column index 
                  subConnectionIndex = neighborSubCellIndexes( n, 3 )
                  neighborSubRow     = 1 
                  neighborSubColumn  = 1

                  print *, '   -- Y DIRECTION THEN LOCATION IN BUFFER IS GIVEN BY SUBCELL COLUMN',&
                  neighborCellData( 2, neighborSubCellIndexes( n, 1 ) )%requestedFromDirection
                  !call neighborCellData(&
                  !    neighborSubCellIndexes( n, 3 ), &
                  !    neighborSubCellIndexes( n, 1 ) )%FillSubCellFaceFlowsBuffer( &
                  !                          1, 1, neighborSubCellFaceFlows( n, : ) )
                  print *, '   -- FILLING FROM CELLNUMBER',&
                  neighborCellData(&
                      neighborSubCellIndexes( n, 3 ), &
                      neighborSubCellIndexes( n, 1 ) )%CellNumber

              else 
                  ! If z-direction, location in buffer is also given 
                  ! by subcell column index. This is true because double 
                  ! buffers requested from a vertical face are possible 
                  ! only if horizontal buffer is double too. This means 
                  ! that this case is related to indirect buffers, which 
                  ! in the current protocol of neighbor cells initialization 
                  ! are only possible after initialization of a direct connection 
                  ! through the y-direction. 
                  subConnectionIndex = neighborSubCellIndexes( n, 3 )
                  neighborSubRow     = 1 
                  neighborSubColumn  = 1

                  print *, '   -- Z DIRECTION THEN LOCATION IN BUFFER IS GIVEN BY SUBCELL COLUMN TOO'

                  !call neighborCellData(&
                  !    neighborSubCellIndexes( n, 3 ), &
                  !    neighborSubCellIndexes( n, 1 ) )%FillSubCellFaceFlowsBuffer( &
                  !                          1, 1, neighborSubCellFaceFlows( n, : ) )

                  print *, '   -- FILLING FROM CELLNUMBER',&
                  neighborCellData(&
                      neighborSubCellIndexes( n, 3 ), &
                      neighborSubCellIndexes( n, 1 ) )%CellNumber


              end if 

          else if ( &
              ( neighborCellData( 1, neighborSubCellIndexes( n, 1 ) )%CellNumber .gt. 0 ) .or. & 
              ( neighborCellData( 1, neighborSubCellIndexes( n, 1 ) )%fromSubCell ) ) then 
              print *, '   -- This is a SINGLE BUFFER'

              if ( neighborCellData(1, neighborSubCellIndexes( n, 1 ) )%GetSubCellCount() .gt. 1 ) then 
                  ! If its refined request the sub cell indexes
                  subConnectionIndex = 1
                  neighborSubRow     = neighborSubCellIndexes( n, 2 )
                  neighborSubColumn  = neighborSubCellIndexes( n, 3 ) 

                  print *, '         -- IS REFINED  ', ' GET THE REQUESTED SUB CELL INDEXES'

                  !call neighborCellData( 1, neighborSubCellIndexes( n, 1 ) )%FillSubCellFaceFlowsBuffer( &
                  !                       neighborSubCellIndexes( n, 2 ), neighborSubCellIndexes( n, 3 ), & 
                  !                                                      neighborSubCellFaceFlows( n, : ) )

              else
                  ! If its not refined, all ones
                  subConnectionIndex = 1
                  neighborSubRow     = 1
                  neighborSubColumn  = 1

                  print *, '         -- NOT REFINED  '
                  if ( neighborCellData( 1, neighborSubCellIndexes( n, 1 ) )%fromSubCell )  then 
                    print *, '              -- WAS FILLED FROM SUBCELL  '

                    !call neighborCellData( 1, neighborSubCellIndexes( n, 1 ) )%FillSubCellFaceFlowsBuffer( &
                    !                                                1, 1, neighborSubCellFaceFlows( n, : ) )

                    !subConnectionIndex = 1
                    !neighborSubRow     = 1
                    !neighborSubColumn  = 1

                  else
                    print *, '              -- WAS FILLED FROM NORMAL CELL  '

                    !call neighborCellData( 1, neighborSubCellIndexes( n, 1 ) )%FillSubCellFaceFlowsBuffer( &
                    !                                                1, 1, neighborSubCellFaceFlows( n, : ) )
                  end if 
              end if  

          else
              ! Buffer is empty, continue
              print *, '   -- This is a EMPTY BUFFER'
              cycle
          end if  

          ! Fill them 
          !call neighborCellData( subConnectionIndex, neighborSubCellIndexes( n, 1 ) )%FillSubCellFaceFlowsBuffer( &
          !                                 neighborSubRow, neighborSubColumn, neighborSubCellFaceFlows( n, : ) )
          call neighborCellData( subConnectionIndex, neighborSubCellIndexes( n, 1 ) )%FillSubCellFaceFlows( &
                          neighborSubRow, neighborSubColumn, neighborSubCellFaceFlows( n, : ), skipSubCells )

          call neighborCellData( subConnectionIndex, neighborSubCellIndexes( n, 1 ) )%FillSubCellFaceAreas( &
                                                             neighborSubCellFaceAreas( n, : ), skipSubCells )  


      end do


      ! Done
      return


  end subroutine pr_FillNeighborSubCellFaceFlowsAreas



end module TrackSubCellModule



!!! THRASH 

!      call exit(0)
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!      ! Following indexes are fixed,
!      ! do not depend on subRow, subColumn
!
!      ! Convention for cell indexes of 
!      ! coordinate components
!      ! 1: 000
!      ! 2: 100
!      ! 3: 010
!      ! 4: 110
!      ! 5: 001
!      ! 6: 101
!      ! 7: 011
!      ! 8: 111
!
!      ! Cell indexes for x-component
!      ! Complete connection id
!      ! 000: 7-8-13   : fn1 : 
!      ! 100: 7-8-13   : fn2 : 
!      ! 010: 10-11-13 : fn1 : 
!      ! 110: 10-11-13 : fn2 :  
!      ! 001: 7-9-16   : fn1 :
!      ! 101: 7-9-16   : fn2 :
!      ! 011: 10-12-16 : fn1 :
!      ! 111: 10-12-16 : fn2 :
!      cornerXComponentIndexes(1,:) = [ 7, 8,13,1] 
!      cornerXComponentIndexes(2,:) = [ 7, 8,13,2] 
!      cornerXComponentIndexes(3,:) = [10,11,13,1] 
!      cornerXComponentIndexes(4,:) = [10,11,13,2] 
!      cornerXComponentIndexes(5,:) = [ 7, 9,16,1] 
!      cornerXComponentIndexes(6,:) = [ 7, 9,16,2] 
!      cornerXComponentIndexes(7,:) = [10,12,16,1] 
!      cornerXComponentIndexes(8,:) = [10,12,16,2] 
!
!      ! Cell indexes for y-component
!      ! Complete connection id
!      !000: 1-13-14 : fn3 :
!      !100: 4-13-15 : fn3 :
!      !010: 1-13-14 : fn4 :
!      !110: 4-13-15 : fn4 :
!      !001: 1-16-17 : fn3 :
!      !101: 4-16-18 : fn3 :
!      !011: 1-16-17 : fn4 :
!      !111: 4-16-18 : fn4 :
!      cornerYComponentIndexes(1,:) = [ 1,13,14,3] 
!      cornerYComponentIndexes(2,:) = [ 4,13,15,3] 
!      cornerYComponentIndexes(3,:) = [ 1,13,14,4] 
!      cornerYComponentIndexes(4,:) = [ 4,13,15,4] 
!      cornerYComponentIndexes(5,:) = [ 1,16,17,3] 
!      cornerYComponentIndexes(6,:) = [ 4,16,18,3] 
!      cornerYComponentIndexes(7,:) = [ 1,16,17,4] 
!      cornerYComponentIndexes(8,:) = [ 4,16,18,4]
!
!      ! Cell indexes for z-component
!      ! Complete connection id
!      !000: 1-2-7   : fn5 :
!      !100: 4-5-7   : fn5 :
!      !010: 1-3-10  : fn5 :
!      !110: 4-6-10  : fn5 :
!      !001: 1-2-7   : fn6 :
!      !101: 4-5-7   : fn6 :
!      !011: 1-3-10  : fn6 :
!      !111: 4-6-10  : fn6 :
!      cornerZComponentIndexes(1,:) = [ 1, 2, 7,5] 
!      cornerZComponentIndexes(2,:) = [ 4, 5, 7,5] 
!      cornerZComponentIndexes(3,:) = [ 1, 3,10,5] 
!      cornerZComponentIndexes(4,:) = [ 4, 6,10,5] 
!      cornerZComponentIndexes(5,:) = [ 1, 2, 7,6] 
!      cornerZComponentIndexes(6,:) = [ 4, 5, 7,6] 
!      cornerZComponentIndexes(7,:) = [ 1, 3,10,6] 
!      cornerZComponentIndexes(8,:) = [ 4, 6,10,6]



    
      !!! subcell faceflows buffer is filled, now use it for interpolation
      !!! to corner values

      !!! Compute cell sizes for areas computation
      !!dx  = this%SubCellData%DX
      !!dy  = this%SubCellData%DY
      !!dz  = this%SubCellData%DZ
      !!dz0 = neighborCellData(13)%GetDZ() ! query a lower layer cell and getdz 
      !!dz1 = neighborCellData(16)%GetDZ() ! query an upper layer cell and getdz 


      !!! Missing areas, should divide sum of face flows
      !!flowContributionFactor = 0.25 ! fraction of total face flow contributed to corner values
      !!areaFlowX0             = flowContributionFactor*( dy*dz + dy*dz + dy*dz0 + dy*dz0 ) 
      !!areaFlowX1             = flowContributionFactor*( dy*dz + dy*dz + dy*dz1 + dy*dz1 ) 
      !!areaFlowY0             = flowContributionFactor*( dx*dz + dx*dz + dx*dz0 + dx*dz0 ) 
      !!areaFlowY1             = flowContributionFactor*( dx*dz + dx*dz + dx*dz1 + dx*dz1 ) 
      !!areaFlowZ              = dy*dx


      ! NO
      !!! function pr_ComputeCornerFaceArea( this, neighborSubCellFaceAreas, cornerIndex, faceDirection ) result( faceArea )
      !!! select case( faceDirection )
      !!!     case (1)
      !!!         faceArea = this%SubCellData%DY*this%SubCellData%DZ
      !!!         do m = 1, 3
      !!!             faceArea = faceArea + &
      !!!                 neighborSubCellFaceAreas( this%cornerXComponentIndexes( cornerIndex, m ), 1 )
      !!!         end do
      !!!     case (2)
      !!!         faceArea = this%SubCellData%DX*this%SubCellData%DZ
      !!!         do m = 1, 3
      !!!             faceArea = faceArea + &
      !!!                 neighborSubCellFaceAreas( this%cornerYComponentIndexes( cornerIndex, m ), 2 )
      !!!         end do
      !!!     case (3)
      !!!         faceArea = this%SubCellData%DX*this%SubCellData%DY
      !!!         do m = 1, 3
      !!!             faceArea = faceArea + &
      !!!                 neighborSubCellFaceAreas( this%cornerZComponentIndexes( cornerIndex, m ), 3 )
      !!!         end do
      !!! end select
      !!! end function pr_ComputeCornerFaceArea






      !!! Compute x-components
      !!this%qCorner000(1) = flowContributionFactor*(&
      !!        centerSubCellFaceFlows(   cornerXComponentIndexes(1,4) ) + &
      !!        neighborSubCellFaceFlows( cornerXComponentIndexes(1,1), cornerXComponentIndexes(1,4) ) + & 
      !!        neighborSubCellFaceFlows( cornerXComponentIndexes(1,2), cornerXComponentIndexes(1,4) ) + &
      !!        neighborSubCellFaceFlows( cornerXComponentIndexes(1,3), cornerXComponentIndexes(1,4) )   &
      !!    )/areaFlowX0
      !!this%qCorner100(1) = flowContributionFactor*(&
      !!        centerSubCellFaceFlows(   cornerXComponentIndexes(2,4) ) + &
      !!        neighborSubCellFaceFlows( cornerXComponentIndexes(2,1), cornerXComponentIndexes(2,4) ) + & 
      !!        neighborSubCellFaceFlows( cornerXComponentIndexes(2,2), cornerXComponentIndexes(2,4) ) + &
      !!        neighborSubCellFaceFlows( cornerXComponentIndexes(2,3), cornerXComponentIndexes(2,4) )   &
      !!    )/areaFlowX0
      !!this%qCorner010(1) = flowContributionFactor*(&
      !!        centerSubCellFaceFlows(   cornerXComponentIndexes(3,4) ) + &
      !!        neighborSubCellFaceFlows( cornerXComponentIndexes(3,1), cornerXComponentIndexes(3,4) ) + & 
      !!        neighborSubCellFaceFlows( cornerXComponentIndexes(3,2), cornerXComponentIndexes(3,4) ) + &
      !!        neighborSubCellFaceFlows( cornerXComponentIndexes(3,3), cornerXComponentIndexes(3,4) )   &
      !!    )/areaFlowX0
      !!this%qCorner110(1) = flowContributionFactor*(&
      !!        centerSubCellFaceFlows(   cornerXComponentIndexes(4,4) ) + &
      !!        neighborSubCellFaceFlows( cornerXComponentIndexes(4,1), cornerXComponentIndexes(4,4) ) + & 
      !!        neighborSubCellFaceFlows( cornerXComponentIndexes(4,2), cornerXComponentIndexes(4,4) ) + &
      !!        neighborSubCellFaceFlows( cornerXComponentIndexes(4,3), cornerXComponentIndexes(4,4) )   &
      !!    )/areaFlowX0
      !!this%qCorner001(1) = flowContributionFactor*(&
      !!        centerSubCellFaceFlows(   cornerXComponentIndexes(5,4) ) + &
      !!        neighborSubCellFaceFlows( cornerXComponentIndexes(5,1), cornerXComponentIndexes(5,4) ) + & 
      !!        neighborSubCellFaceFlows( cornerXComponentIndexes(5,2), cornerXComponentIndexes(5,4) ) + &
      !!        neighborSubCellFaceFlows( cornerXComponentIndexes(5,3), cornerXComponentIndexes(5,4) )   &
      !!    )/areaFlowX1
      !!this%qCorner101(1) = flowContributionFactor*(&
      !!        centerSubCellFaceFlows(   cornerXComponentIndexes(6,4) ) + &
      !!        neighborSubCellFaceFlows( cornerXComponentIndexes(6,1), cornerXComponentIndexes(6,4) ) + & 
      !!        neighborSubCellFaceFlows( cornerXComponentIndexes(6,2), cornerXComponentIndexes(6,4) ) + &
      !!        neighborSubCellFaceFlows( cornerXComponentIndexes(6,3), cornerXComponentIndexes(6,4) )   &
      !!    )/areaFlowX1
      !!this%qCorner011(1) = flowContributionFactor*(&
      !!        centerSubCellFaceFlows(   cornerXComponentIndexes(7,4) ) + &
      !!        neighborSubCellFaceFlows( cornerXComponentIndexes(7,1), cornerXComponentIndexes(7,4) ) + & 
      !!        neighborSubCellFaceFlows( cornerXComponentIndexes(7,2), cornerXComponentIndexes(7,4) ) + &
      !!        neighborSubCellFaceFlows( cornerXComponentIndexes(7,3), cornerXComponentIndexes(7,4) )   &
      !!    )/areaFlowX1
      !!this%qCorner111(1) = flowContributionFactor*(&
      !!        centerSubCellFaceFlows(   cornerXComponentIndexes(8,4) ) + &
      !!        neighborSubCellFaceFlows( cornerXComponentIndexes(8,1), cornerXComponentIndexes(8,4) ) + & 
      !!        neighborSubCellFaceFlows( cornerXComponentIndexes(8,2), cornerXComponentIndexes(8,4) ) + &
      !!        neighborSubCellFaceFlows( cornerXComponentIndexes(8,3), cornerXComponentIndexes(8,4) )   & 
      !!    )/areaFlowX1


      !!! Compute y-components
      !!this%qCorner000(2) = flowContributionFactor*(&
      !!        centerSubCellFaceFlows(   cornerYComponentIndexes(1,4) ) + &
      !!        neighborSubCellFaceFlows( cornerYComponentIndexes(1,1), cornerYComponentIndexes(1,4) ) + & 
      !!        neighborSubCellFaceFlows( cornerYComponentIndexes(1,2), cornerYComponentIndexes(1,4) ) + &
      !!        neighborSubCellFaceFlows( cornerYComponentIndexes(1,3), cornerYComponentIndexes(1,4) )   &
      !!    )/areaFlowY0
      !!this%qCorner100(2) = flowContributionFactor*(&
      !!        centerSubCellFaceFlows(   cornerYComponentIndexes(2,4) ) + &
      !!        neighborSubCellFaceFlows( cornerYComponentIndexes(2,1), cornerYComponentIndexes(2,4) ) + & 
      !!        neighborSubCellFaceFlows( cornerYComponentIndexes(2,2), cornerYComponentIndexes(2,4) ) + &
      !!        neighborSubCellFaceFlows( cornerYComponentIndexes(2,3), cornerYComponentIndexes(2,4) )   &
      !!    )/areaFlowY0
      !!this%qCorner010(2) = flowContributionFactor*(&
      !!        centerSubCellFaceFlows(   cornerYComponentIndexes(3,4) ) + &
      !!        neighborSubCellFaceFlows( cornerYComponentIndexes(3,1), cornerYComponentIndexes(3,4) ) + & 
      !!        neighborSubCellFaceFlows( cornerYComponentIndexes(3,2), cornerYComponentIndexes(3,4) ) + &
      !!        neighborSubCellFaceFlows( cornerYComponentIndexes(3,3), cornerYComponentIndexes(3,4) )   &
      !!    )/areaFlowY0
      !!this%qCorner110(2) = flowContributionFactor*(&
      !!        centerSubCellFaceFlows(   cornerYComponentIndexes(4,4) ) + &
      !!        neighborSubCellFaceFlows( cornerYComponentIndexes(4,1), cornerYComponentIndexes(4,4) ) + & 
      !!        neighborSubCellFaceFlows( cornerYComponentIndexes(4,2), cornerYComponentIndexes(4,4) ) + &
      !!        neighborSubCellFaceFlows( cornerYComponentIndexes(4,3), cornerYComponentIndexes(4,4) )   &
      !!    )/areaFlowY0
      !!this%qCorner001(2) = flowContributionFactor*(&
      !!        centerSubCellFaceFlows(   cornerYComponentIndexes(5,4) ) + &
      !!        neighborSubCellFaceFlows( cornerYComponentIndexes(5,1), cornerYComponentIndexes(5,4) ) + & 
      !!        neighborSubCellFaceFlows( cornerYComponentIndexes(5,2), cornerYComponentIndexes(5,4) ) + &
      !!        neighborSubCellFaceFlows( cornerYComponentIndexes(5,3), cornerYComponentIndexes(5,4) )   &
      !!    )/areaFlowY1
      !!this%qCorner101(2) = flowContributionFactor*(&
      !!        centerSubCellFaceFlows(   cornerYComponentIndexes(6,4) ) + &
      !!        neighborSubCellFaceFlows( cornerYComponentIndexes(6,1), cornerYComponentIndexes(6,4) ) + & 
      !!        neighborSubCellFaceFlows( cornerYComponentIndexes(6,2), cornerYComponentIndexes(6,4) ) + &
      !!        neighborSubCellFaceFlows( cornerYComponentIndexes(6,3), cornerYComponentIndexes(6,4) )   &
      !!    )/areaFlowY1
      !!this%qCorner011(2) = flowContributionFactor*(&
      !!        centerSubCellFaceFlows(   cornerYComponentIndexes(7,4) ) + &
      !!        neighborSubCellFaceFlows( cornerYComponentIndexes(7,1), cornerYComponentIndexes(7,4) ) + & 
      !!        neighborSubCellFaceFlows( cornerYComponentIndexes(7,2), cornerYComponentIndexes(7,4) ) + &
      !!        neighborSubCellFaceFlows( cornerYComponentIndexes(7,3), cornerYComponentIndexes(7,4) )   &
      !!    )/areaFlowY1
      !!this%qCorner111(2) = flowContributionFactor*(&
      !!        centerSubCellFaceFlows(   cornerYComponentIndexes(8,4) ) + &
      !!        neighborSubCellFaceFlows( cornerYComponentIndexes(8,1), cornerYComponentIndexes(8,4) ) + & 
      !!        neighborSubCellFaceFlows( cornerYComponentIndexes(8,2), cornerYComponentIndexes(8,4) ) + &
      !!        neighborSubCellFaceFlows( cornerYComponentIndexes(8,3), cornerYComponentIndexes(8,4) )   &
      !!    )/areaFlowY1


      !!! Compute z-components
      !!this%qCorner000(3) = flowContributionFactor*(&
      !!        centerSubCellFaceFlows(   cornerZComponentIndexes(1,4) ) + &
      !!        neighborSubCellFaceFlows( cornerZComponentIndexes(1,1), cornerZComponentIndexes(1,4) ) + & 
      !!        neighborSubCellFaceFlows( cornerZComponentIndexes(1,2), cornerZComponentIndexes(1,4) ) + &
      !!        neighborSubCellFaceFlows( cornerZComponentIndexes(1,3), cornerZComponentIndexes(1,4) )   &
      !!    )/areaFlowZ
      !!this%qCorner100(3) = flowContributionFactor*(&
      !!        centerSubCellFaceFlows(   cornerZComponentIndexes(2,4) ) + &
      !!        neighborSubCellFaceFlows( cornerZComponentIndexes(2,1), cornerZComponentIndexes(2,4) ) + & 
      !!        neighborSubCellFaceFlows( cornerZComponentIndexes(2,2), cornerZComponentIndexes(2,4) ) + &
      !!        neighborSubCellFaceFlows( cornerZComponentIndexes(2,3), cornerZComponentIndexes(2,4) )   &
      !!    )/areaFlowZ
      !!this%qCorner010(3) = flowContributionFactor*(&
      !!        centerSubCellFaceFlows(   cornerZComponentIndexes(3,4) ) + &
      !!        neighborSubCellFaceFlows( cornerZComponentIndexes(3,1), cornerZComponentIndexes(3,4) ) + & 
      !!        neighborSubCellFaceFlows( cornerZComponentIndexes(3,2), cornerZComponentIndexes(3,4) ) + &
      !!        neighborSubCellFaceFlows( cornerZComponentIndexes(3,3), cornerZComponentIndexes(3,4) )   &
      !!    )/areaFlowZ
      !!this%qCorner110(3) = flowContributionFactor*(&
      !!        centerSubCellFaceFlows(   cornerZComponentIndexes(4,4) ) + &
      !!        neighborSubCellFaceFlows( cornerZComponentIndexes(4,1), cornerZComponentIndexes(4,4) ) + & 
      !!        neighborSubCellFaceFlows( cornerZComponentIndexes(4,2), cornerZComponentIndexes(4,4) ) + &
      !!        neighborSubCellFaceFlows( cornerZComponentIndexes(4,3), cornerZComponentIndexes(4,4) )   &
      !!    )/areaFlowZ
      !!this%qCorner001(3) = flowContributionFactor*(&
      !!        centerSubCellFaceFlows(   cornerZComponentIndexes(5,4) ) + &
      !!        neighborSubCellFaceFlows( cornerZComponentIndexes(5,1), cornerZComponentIndexes(5,4) ) + & 
      !!        neighborSubCellFaceFlows( cornerZComponentIndexes(5,2), cornerZComponentIndexes(5,4) ) + &
      !!        neighborSubCellFaceFlows( cornerZComponentIndexes(5,3), cornerZComponentIndexes(5,4) )   &
      !!    )/areaFlowZ
      !!this%qCorner101(3) = flowContributionFactor*(&
      !!        centerSubCellFaceFlows(   cornerZComponentIndexes(6,4) ) + &
      !!        neighborSubCellFaceFlows( cornerZComponentIndexes(6,1), cornerZComponentIndexes(6,4) ) + & 
      !!        neighborSubCellFaceFlows( cornerZComponentIndexes(6,2), cornerZComponentIndexes(6,4) ) + &
      !!        neighborSubCellFaceFlows( cornerZComponentIndexes(6,3), cornerZComponentIndexes(6,4) )   &
      !!    )/areaFlowZ
      !!this%qCorner011(3) = flowContributionFactor*(&
      !!        centerSubCellFaceFlows(   cornerZComponentIndexes(7,4) ) + &
      !!        neighborSubCellFaceFlows( cornerZComponentIndexes(7,1), cornerZComponentIndexes(7,4) ) + & 
      !!        neighborSubCellFaceFlows( cornerZComponentIndexes(7,2), cornerZComponentIndexes(7,4) ) + &
      !!        neighborSubCellFaceFlows( cornerZComponentIndexes(7,3), cornerZComponentIndexes(7,4) )   &
      !!    )/areaFlowZ
      !!this%qCorner111(3) = flowContributionFactor*(&
      !!        centerSubCellFaceFlows(   cornerZComponentIndexes(8,4) ) + &
      !!        neighborSubCellFaceFlows( cornerZComponentIndexes(8,1), cornerZComponentIndexes(8,4) ) + & 
      !!        neighborSubCellFaceFlows( cornerZComponentIndexes(8,2), cornerZComponentIndexes(8,4) ) + &
      !!        neighborSubCellFaceFlows( cornerZComponentIndexes(8,3), cornerZComponentIndexes(8,4) )   &
      !!    )/areaFlowZ

