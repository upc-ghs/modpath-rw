module TrackSubCellModule
  use ParticleLocationModule,only : ParticleLocationType
  use TrackSubCellResultModule,only : TrackSubCellResultType
  use ModpathSubCellDataModule,only : ModpathSubCellDataType
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

  contains
    procedure,private :: CalculateDT=>pr_CalculateDT
    procedure,private :: NewXYZ=>pr_NewXYZ
    procedure :: ExecuteTracking=>pr_ExecuteTracking

    ! RWPT
    procedure :: ComputeCornerVelocities=>pr_ComputeCornerVelocities
    procedure :: Trilinear=>pr_Trilinear
    procedure :: TrilinearDerivative=>pr_TrilinearDerivative
    procedure :: DispersionDivergence=>pr_DispersionDivergence
    procedure :: DisplacementRandom=>pr_DisplacementRandom
    procedure :: DetectExitFaceAndUpdateTimeStep=>pr_DetectExitFaceAndUpdateTimeStep

  end type
  
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

  ! RWPT
  doubleprecision :: alphaT, alphaL, Dmol
  doubleprecision :: divDx, divDy, divDz
  doubleprecision :: dBx, dBy, dBz
  doubleprecision :: dx, dy, dz
  doubleprecision :: nx, ny, nz
  doubleprecision :: dxrw, dyrw, dzrw
  logical         :: continueTimeLoop
  logical         :: reachedMaximumTime
  doubleprecision, dimension(3) :: dts
  !------------------------------------------------------------

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

  ! Local copies of cell size
  dx = this%SubCellData%DX
  dy = this%SubCellData%DY
  dz = this%SubCellData%DZ

  ! Initialize positions
  x = initialLocation%LocalX
  y = initialLocation%LocalY
  z = initialLocation%LocalZ


  ! RWPT
  ! COMPUTE A DELTA T FOR THE TRACK CELL
  ! - NOTES:
  !     - CalculateDT performs the linear interpolation of velocities
  !     - Before computing dispersion, requires velocities x,y,z
  !     - CalculateDT defines three values of dt, one for each axis
  !       which are then compared and used in computing the final 
  !       position of particles 
  ! RWPT
  ! In the case with hydrodynamic dispersion 
  ! particles should be moved in several time 
  ! steps, in each tracked cell.
  ! This is essentially due to the fact that 
  ! the exit face of the particle is not know beforehand
  ! as there is random motion of particles and dispersion.
  ! Alternatives for computing delta t
  !     - Constant
  !     - Constant Courant 
  !     - Constant Peclet
  !     - Several alternatives 
  ! The thing is, that given a deltat, 
  ! there should a loop cycle that should run 
  ! until the particle leaves the cell,
  ! which requires a detection mechanism
  ! At each time iteration, quantities related to 
  ! RWPT should be recomputed 
  statusVX = this%CalculateDT(vx1, vx2, this%SubCellData%DX, initialX, vx, dvxdx, dtx)
  statusVY = this%CalculateDT(vy1, vy2, this%SubCellData%DY, initialY, vy, dvydy, dty)
  statusVZ = this%CalculateDT(vz1, vz2, this%SubCellData%DZ, initialZ, vz, dvzdz, dtz)

  dts(1) = dtx
  dts(2) = dty
  dts(3) = dtz
  dt = 0.001*minval( dts, dts > 0 )
  

  ! At this point we have defined a reasonable dt
  t = initialTime + dt

  alphaL = 0
  alphaT = 0
  Dmol   = 0

  !------- 
  ! RWPT
  !-------
  exitFace = 0
  continueTimeLoop = .true.
  reachedMaximumTime = .false.
  !LOCAL TIME LOOP
  do while( continueTimeLoop )

      ! Recompute dt for maximumTime 
      if (maximumTime .lt. t) then
          dt = t - maximumTime
          reachedMaximumTime = .true.
      end if 


      ! Move particle

      ! Interpolate velocities and
      ! compute RWPT movement
      vx = ( 1.0d0 - x )*vx1 + x*vx2
      vy = ( 1.0d0 - y )*vy1 + y*vy2
      vz = ( 1.0d0 - z )*vz1 + z*vz2
      call this%DispersionDivergence( x, y, z, alphaL, alphaT, Dmol, divDx, divDy, divDz )
      call this%DisplacementRandom( x, y, z, alphaL, alphaT, Dmol, dBx, dBy, dBz )
      dxrw = ( vx + divDx )*dt + dBx*sqrt( dt )
      dyrw = ( vy + divDy )*dt + dBy*sqrt( dt )
      dzrw = ( vz + divDz )*dt + dBz*sqrt( dt )
      nx   = x + dxrw/dx
      ny   = y + dyrw/dy
      nz   = z + dzrw/dz

      !! Very unlikely for particles to fall exactly
      !! on the interface

      ! Detect if particle leaving the cell
      ! and force the particle into exactly one
      ! interface by computing the required dt
      if (                                             &
          ( nx .gt. 1.0d0 ) .or. ( nx .lt. 0d0 )  .or. &
          ( ny .gt. 1.0d0 ) .or. ( ny .lt. 0d0 )  .or. &
          ( nz .gt. 1.0d0 ) .or. ( nz .lt. 0d0 )       & 
      ) then                                           
      
          call this%DetectExitFaceAndUpdateTimeStep(              &
                  x, y, z, nx, ny, nz,                            &
                  ( vx + divDx ), ( vy + divDy ), ( vz + divDz ), &
                  dBx, dBy, dBz, t, dt, exitFace                  &
              )

          ! Find new displacements using the recomputed dt
          ! and update coordinates
          if ( ( exitFace .eq. 1 ) .or. ( exitFace .eq. 2 ) ) then 
              dyrw = ( vy + divDy )*dt + dBy*sqrt( dt )
              dzrw = ( vz + divDz )*dt + dBz*sqrt( dt )
              ny = y + dyrw/dy
              nz = z + dzrw/dz
          else if ( ( exitFace .eq. 3 ) .or. ( exitFace .eq. 4 ) ) then 
              dxrw = ( vx + divDx )*dt + dBx*sqrt( dt )
              dzrw = ( vz + divDz )*dt + dBz*sqrt( dt )
              nx = x + dxrw/dx
              nz = z + dzrw/dz
          else if ( ( exitFace .eq. 5 ) .or. ( exitFace .eq. 6 ) ) then
              dxrw = ( vx + divDx )*dt + dBx*sqrt( dt )
              dyrw = ( vy + divDy )*dt + dBy*sqrt( dt )
              nx = x + dxrw/dx
              ny = y + dyrw/dy
          else
              ! Something wrong
              trackingResult%ExitFace = exitFace
              trackingResult%Status = trackingResult%Status_Undefined()
              return
          end if
     
      end if

      ! Update particle positions and time
      x = nx
      y = ny
      z = nz
      t = t + dt

      ! Report and leave
      if ( (reachedMaximumTime) .or. ( exitFace .ne. 0 ) ) then
          if (reachedMaximumTime) then
              ! In this case it is allowed for a particle
              ! to reach the maximum time and eventually have an exit face
              trackingResult%Status = trackingResult%Status_ReachedMaximumTime()
          else
              trackingResult%ExitFaceConnection = this%SubCellData%Connection(exitFace)
              if(this%SubCellData%Connection(exitFace) .lt. 0) then
                  ! This is the case for unstructured grid, NOT IMPLEMENTED YET
                  trackingResult%Status = trackingResult%Status_ExitAtInternalFace()
              else
                  trackingResult%Status = trackingResult%Status_ExitAtCellFace()
              end if
          end if
          trackingResult%ExitFace = exitFace
          trackingResult%FinalLocation%CellNumber = cellNumber
          trackingResult%FinalLocation%LocalX = x
          trackingResult%FinalLocation%LocalY = y
          trackingResult%FinalLocation%LocalZ = z
          trackingResult%FinalLocation%TrackingTime = t
          continueTimeLoop = .false.         
      end if
      
  end do 
  ! END LOCAL TIME LOOP 


  return






  !!!----------------------
  !!! ORIGINAL MODPATH 
  !!!---------------------- 
  !!call trackingResult%Reset()
  !!
  !!cellNumber = initialLocation%CellNumber
  !!initialX = initialLocation%LocalX
  !!initialY = initialLocation%LocalY
  !!initialZ = initialLocation%LocalZ
  !!initialTime = initialLocation%TrackingTime
  !!
  !!trackingResult%CellNumber = cellNumber
  !!trackingResult%Row = this%SubCellData%Row
  !!trackingResult%Column = this%SubCellData%Column
  !!trackingResult%InitialLocation%CellNumber = cellNumber
  !!trackingResult%InitialLocation%LocalX = initialX
  !!trackingResult%InitialLocation%LocalY = initialY
  !!trackingResult%InitialLocation%LocalZ = initialZ
  !!trackingResult%InitialLocation%TrackingTime = initialTime
  !!trackingResult%FinalLocation%LocalX = initialX
  !!trackingResult%FinalLocation%LocalY = initialY
  !!trackingResult%FinalLocation%LocalZ = initialZ
  !!trackingResult%FinalLocation%TrackingTime = initialTime
  !!trackingResult%MaximumTime = maximumTime
  !!trackingResult%Status = trackingResult%Status_Undefined()
  !!
  !!if(stopIfNoExit) then
  !!  if(.not. this%SubCellData%HasExitFace()) then
  !!    trackingResult%Status = trackingResult%Status_NoExitPossible()
  !!    return
  !!  end if
  !!end if
  !!
  !!! Make local copies of face velocities for convenience
  !!vx1 = this%SubCellData%VX1
  !!vx2 = this%SubCellData%VX2
  !!vy1 = this%SubCellData%VY1
  !!vy2 = this%SubCellData%VY2
  !!vz1 = this%SubCellData%VZ1
  !!vz2 = this%SubCellData%VZ2

  !!statusVX = this%CalculateDT(vx1, vx2, this%SubCellData%DX, initialX, vx, dvxdx, dtx)
  !!statusVY = this%CalculateDT(vy1, vy2, this%SubCellData%DY, initialY, vy, dvydy, dty)
  !!statusVZ = this%CalculateDT(vz1, vz2, this%SubCellData%DZ, initialZ, vz, dvzdz, dtz)

  !!! This snippet selects dt and determines the exitFace
  !!exitFace = 0
  !!dt = 1.0d+30
  !!if((statusVX .lt. 2) .or. (statusVY .lt. 2) .or. (statusVZ .lt. 2)) then
  !!  dt = dtx
  !!  if(vx .lt. 0d0) then
  !!    exitFace = 1
  !!  else if(vx .gt. 0) then
  !!    exitFace = 2
  !!  end if
  !!  
  !!  if(dty .lt. dt) then
  !!    dt = dty
  !!    if(vy .lt. 0d0) then
  !!      exitFace = 3
  !!    else if(vy .gt. 0d0) then
  !!      exitFace = 4
  !!    end if
  !!  end if
  !!  
  !!  if(dtz .lt. dt) then
  !!    dt = dtz
  !!    if(vz .lt. 0d0) then
  !!      exitFace = 5
  !!    else if(vz .gt. 0d0) then
  !!      exitFace = 6
  !!    end if
  !!  end if
  !!else
  !!end if

  !!! Increment t
  !!t = initialTime + dt
  !!
  !!! If the maximum time is less than the computed exit time, then calculate the
  !!! particle location at the maximum time and set the final time equal to
  !!! the maximum time. Set status to 2 (MaximumTimeReached) and return.
  !!if(maximumTime .lt. t) then
  !!  dt = maximumTime - initialTime
  !!  t = maximumTime
  !!  x = this%NewXYZ(vx, dvxdx, vx1, vx2, dt, initialX, this%SubCellData%DX, statusVX)
  !!  y = this%NewXYZ(vy, dvydy, vy1, vy2, dt, initialY, this%SubCellData%DY, statusVY)
  !!  z = this%NewXYZ(vz, dvzdz, vz1, vz2, dt, initialZ, this%SubCellData%DZ, statusVZ)

  !!  trackingResult%ExitFace = 0
  !!  trackingResult%FinalLocation%CellNumber = cellNumber
  !!  trackingResult%FinalLocation%LocalX = x
  !!  trackingResult%FinalLocation%LocalY = y
  !!  trackingResult%FinalLocation%LocalZ = z
  !!  trackingResult%FinalLocation%TrackingTime = t
  !!  trackingResult%Status = trackingResult%Status_ReachedMaximumTime()
  !!else
  !!    ! Otherwise, if the computed exit time is less than or equal to the maximum time,
  !!    ! then calculate the exit location and set the final time equal to the computed
  !!    ! exit time.
  !!    if((exitFace .eq. 1) .or. (exitFace .eq.2)) then
  !!        x = 0d0
  !!        y = this%NewXYZ(vy, dvydy, vy1, vy2, dt, initialY, this%SubCellData%DY, statusVY)
  !!        z = this%NewXYZ(vz, dvzdz, vz1, vz2, dt, initialZ, this%SubCellData%DZ, statusVZ)
  !!        if(exitFace .eq. 2) x = 1.0d0
  !!    else if((exitFace .eq. 3) .or. (exitFace .eq.4)) then
  !!        x = this%NewXYZ(vx, dvxdx, vx1, vx2, dt, initialX, this%SubCellData%DX, statusVX)
  !!        y = 0d0
  !!        z = this%NewXYZ(vz, dvzdz, vz1, vz2, dt, initialZ, this%SubCellData%DZ, statusVZ)
  !!        if(exitFace .eq. 4) y = 1.0d0
  !!    else if((exitFace .eq. 5) .or. (exitFace .eq.6)) then
  !!        x = this%NewXYZ(vx, dvxdx, vx1, vx2, dt, initialX, this%SubCellData%DX, statusVX)
  !!        y = this%NewXYZ(vy, dvydy, vy1, vy2, dt, initialY, this%SubCellData%DY, statusVY)
  !!        z = 0d0
  !!        if(exitFace .eq. 6) z = 1.0d0
  !!    else
  !!        ! If it gets this far, something went wrong. Signal an error condition by
  !!        ! setting Status = 0 (undefined) and then return.
  !!        trackingResult%ExitFace = exitFace
  !!        trackingResult%Status = trackingResult%Status_Undefined()
  !!        return
  !!    end if
  !!    
  !!    ! Assign tracking result data
  !!   trackingResult%ExitFaceConnection = this%SubCellData%Connection(exitFace)
  !!    if(this%SubCellData%Connection(exitFace) .lt. 0) then
  !!        trackingResult%Status = trackingResult%Status_ExitAtInternalFace()
  !!    else
  !!        trackingResult%Status = trackingResult%Status_ExitAtCellFace()
  !!    end if
  !!    trackingResult%ExitFace = exitFace
  !!    trackingResult%FinalLocation%CellNumber = cellNumber
  !!    trackingResult%FinalLocation%LocalX = x
  !!    trackingResult%FinalLocation%LocalY = y
  !!    trackingResult%FinalLocation%LocalZ = z
  !!    trackingResult%FinalLocation%TrackingTime = t
  !!end if


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
  subroutine pr_DetectExitFaceAndUpdateTimeStep( this, x, y, z, nx, ny, nz, & 
                                    Ax, Ay, Az, Bx, By, Bz, t, dt, exitFace )
      !----------------------------------------------------------------
      ! Given the initial and updated coordinates, and rwpt terms, 
      ! detects the exit face and recomputes time step to move the 
      ! particle exactly into the corresponding interface 
      ! 
      ! Params:
      !     - x, y, z    | doubleprecision |: initial local cell coordinates
      !     - nx, ny, nz | doubleprecision |: coordinates after RWPT
      !     - Ai, Bi     | doubleprecision |: RWPT terms 
      !     - t, dt      | doubleprecision |: current time and time step
      !     - exitFace   |     integer     |: exit face index
      !----------------------------------------------------------------
      ! Specifications
      !----------------------------------------------------------------
      implicit none
      class (TrackSubCellType) :: this
      ! input
      doubleprecision :: x, y, z
      doubleprecision :: nx, ny, nz
      doubleprecision :: Ax, Ay, Az, Bx, By, Bz
      doubleprecision, intent(inout) :: t, dt
      integer, intent(inout)         :: exitFace
      ! local
      doubleprecision :: AFace, BFace, z1, z2, zsqrt, dInterface
      doubleprecision :: dx, dy, dz
      !----------------------------------------------------------------

      ! Initialize
      AFace = 0
      BFace = 0
      z1    = 0
      z2    = 0
      zsqrt = 0 
      dInterface = 0

      ! Local copies of cell size
      dx = this%SubCellData%DX
      dy = this%SubCellData%DY
      dz = this%SubCellData%DZ

      ! There should be something that analyzes precedence
      ! what happened first, maybe something involving dt
      if ( ( nx .gt. 1.0d0 ) .or. ( nx .lt. 0d0 )  ) then ! leaving through x face
          AFace = Ax
          BFace = Bx
          if ( nx .gt. 1.0d0 ) then 
              dInterface = dx*( 1.0d0 - x )
              nx         = 1.0d0
              exitFace   = 2
          else 
              dInterface = -dx*x
              nx         = 0d0
              exitFace   = 1
          end if 
      else if ( ( ny .gt. 1.0d0 ) .or. ( ny .lt. 0d0 )  ) then ! leaving through y face
          AFace = Ay
          BFace = By
          if ( ny .gt. 1.0d0 ) then 
              dInterface = dy*( 1.0d0 - y )
              ny         = 1.0d0
              exitFace   = 4
          else 
              dInterface = -dy*y
              ny         = 0d0
              exitFace   = 3
          end if
      else if ( ( nz .gt. 1.0d0 ) .or. ( nz .lt. 0d0 )  ) then ! leaving through z face
          AFace = Az
          BFace = Bz
          if ( nz .gt. 1.0d0 ) then
              dInterface = dz*( 1.0d0 - z )
              nz         = 1.0d0
              exitFace   = 6
          else
              dInterface = -dz*z
              nz         = 0d0
              exitFace   = 5
          end if 
      end if

      ! Reset t, dt will be replaced
      t = t - dt

      ! Given dInterface, compute new dt
      zsqrt = sqrt( BFace**2 + 4*dInterface*AFace )
      z1    = (-BFace + zsqrt )/( 2*AFace )
      z2    = (-BFace - zsqrt )/( 2*AFace )
      if ( z1 .gt. 0d0 ) then
          dt = z1**2
      else if ( z2 .gt. 0d0 ) then
          dt = z2**2
      end if

  end subroutine pr_DetectExitFaceAndUpdateTimeStep


  ! RWPT
  subroutine pr_DispersionDivergence( this, x, y, z, alphaL, alphaT, Dmol, divDx, divDy, divDz )
      !----------------------------------------------------------------
      ! Compute dispersion divergence terms 
      ! 
      ! Params:
      !     - x, y, z, doubleprecision: local cell coordinates
      !     - alphaL, doubleprecision: longidutinal dispersivity
      !     - alphaT, doubleprecision: transverse dispersivity
      !     - Dmol, doubleprecision: molecular diffusion
      !     - divDx, divDy, divDz, doubleprecision: output terms
      !----------------------------------------------------------------
      ! Specifications
      !----------------------------------------------------------------
      implicit none
      class (TrackSubCellType)    :: this
      ! input
      doubleprecision, intent(in) :: x, y, z
      doubleprecision, intent(in) :: alphaL, alphaT, Dmol
      ! output
      doubleprecision             :: divDx, divDy, divDz
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

      !! RWPT
      !alphaL  = 1
      !alphaT  = 0.1
      !Dmol    = 1e-8
    
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
  subroutine pr_DisplacementRandom( this, x, y, z, alphaL, alphaT, Dmol, dBx, dBy, dBz ) 
      !----------------------------------------------------------------
      ! Computes the product between displacement matrix and random 
      ! vector
      !
      ! Params:
      !     - x, y, z, doubleprecision: local cell coordinates
      !     - alphaL, doubleprecision: longidutinal dispersivity
      !     - alphaT, doubleprecision: transverse dispersivity
      !     - Dmol, doubleprecision: molecular diffusion
      !     - dBx, dBy, dBz, doubleprecision: output terms
      !----------------------------------------------------------------
      ! Specifications
      !----------------------------------------------------------------
      implicit none
      class (TrackSubCellType) :: this
      ! input
      doubleprecision, intent(in) :: x, y, z
      doubleprecision, intent(in) :: alphaL, alphaT, Dmol
      ! output
      doubleprecision             :: dBx, dBy, dBz
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
    
      ! Displacement terms (Salamon et al. 2006)
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
      call random_number( rdmx )
      call random_number( rdmy )
      call random_number( rdmz )

      ! Compute displacement times random
      dBX = B11*rdmx + B12*rdmy + B13*rdmz 
      dBY = B21*rdmx + B22*rdmy + B23*rdmz 
      dBZ = B31*rdmx + B32*rdmy 

  end subroutine pr_DisplacementRandom






  ! RWPT
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
      ! 
      this%vCorner000(1) = 0.25*( this%SubCellData%VX1 + neighborSubCellData(7)%VX1  + &
          neighborSubCellData(13)%VX1 + neighborSubCellData(8)%VX1 )
      this%vCorner100(1) = 0.25*( this%SubCellData%VX2 + neighborSubCellData(7)%VX2  + &
          neighborSubCellData(13)%VX2 + neighborSubCellData(8)%VX2 )
      this%vCorner010(1) = 0.25*( this%SubCellData%VX1 + neighborSubCellData(10)%VX1 + &
          neighborSubCellData(13)%VX1 + neighborSubCellData(11)%VX1 )
      this%vCorner110(1) = 0.25*( this%SubCellData%VX2 + neighborSubCellData(10)%VX2 + &
          neighborSubCellData(13)%VX2 + neighborSubCellData(11)%VX2 )
      this%vCorner001(1) = 0.25*( this%SubCellData%VX1 + neighborSubCellData(7)%VX1  + &
          neighborSubCellData(16)%VX1 + neighborSubCellData(9)%VX1 )
      this%vCorner101(1) = 0.25*( this%SubCellData%VX2 + neighborSubCellData(7)%VX2  + &
          neighborSubCellData(16)%VX2 + neighborSubCellData(9)%VX2 )
      this%vCorner011(1) = 0.25*( this%SubCellData%VX1 + neighborSubCellData(10)%VX1 + &
          neighborSubCellData(16)%VX1 + neighborSubCellData(12)%VX1 )
      this%vCorner111(1) = 0.25*( this%SubCellData%VX2 + neighborSubCellData(10)%VX2 + &
          neighborSubCellData(16)%VX2 + neighborSubCellData(12)%VX2 )


      this%vCorner000(2) = 0.25*( this%SubCellData%VY1 + neighborSubCellData(1)%VY1  + &
          neighborSubCellData(13)%VY1 + neighborSubCellData(14)%VY1 )
      this%vCorner010(2) = 0.25*( this%SubCellData%VY2 + neighborSubCellData(1)%VY2  + &
          neighborSubCellData(13)%VY2 + neighborSubCellData(14)%VY2 )
      this%vCorner100(2) = 0.25*( this%SubCellData%VY1 + neighborSubCellData(4)%VY1  + &
          neighborSubCellData(13)%VY1 + neighborSubCellData(15)%VY1 )
      this%vCorner110(2) = 0.25*( this%SubCellData%VY2 + neighborSubCellData(4)%VY2  + &
          neighborSubCellData(13)%VY2 + neighborSubCellData(15)%VY2 )
      this%vCorner001(2) = 0.25*( this%SubCellData%VY1 + neighborSubCellData(1)%VY1  + &
          neighborSubCellData(16)%VY1 + neighborSubCellData(17)%VY1 )
      this%vCorner011(2) = 0.25*( this%SubCellData%VY2 + neighborSubCellData(1)%VY2  + &
          neighborSubCellData(16)%VY2 + neighborSubCellData(17)%VY2 )
      this%vCorner101(2) = 0.25*( this%SubCellData%VY1 + neighborSubCellData(4)%VY1  + &
          neighborSubCellData(16)%VY1 + neighborSubCellData(18)%VY1 )
      this%vCorner111(2) = 0.25*( this%SubCellData%VY2 + neighborSubCellData(4)%VY2  + &
          neighborSubCellData(16)%VY2 + neighborSubCellData(18)%VY2 )


      this%vCorner000(3) = 0.25*( this%SubCellData%VZ1 + neighborSubCellData(1)%VZ1  + &
          neighborSubCellData(7)%VZ1 + neighborSubCellData(2)%VZ1 )
      this%vCorner001(3) = 0.25*( this%SubCellData%VZ2 + neighborSubCellData(1)%VZ2  + &
          neighborSubCellData(7)%VZ2 + neighborSubCellData(2)%VZ2 )
      this%vCorner100(3) = 0.25*( this%SubCellData%VZ1 + neighborSubCellData(4)%VZ1  + &
          neighborSubCellData(7)%VZ1 + neighborSubCellData(5)%VZ1 )
      this%vCorner101(3) = 0.25*( this%SubCellData%VZ2 + neighborSubCellData(4)%VZ2  + &
          neighborSubCellData(7)%VZ2 + neighborSubCellData(5)%VZ2 )
      this%vCorner010(3) = 0.25*( this%SubCellData%VZ1 + neighborSubCellData(1)%VZ1  + &
          neighborSubCellData(10)%VZ1 + neighborSubCellData(3)%VZ1 )
      this%vCorner011(3) = 0.25*( this%SubCellData%VZ2 + neighborSubCellData(1)%VZ2  + &
          neighborSubCellData(10)%VZ2 + neighborSubCellData(3)%VZ2 )
      this%vCorner110(3) = 0.25*( this%SubCellData%VZ1 + neighborSubCellData(4)%VZ1  + &
          neighborSubCellData(10)%VZ1 + neighborSubCellData(6)%VZ1 )
      this%vCorner111(3) = 0.25*( this%SubCellData%VZ2 + neighborSubCellData(4)%VZ2  + &
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


  subroutine pr_TrilinearDerivative( this, direction, x, y, z, v000, v100, v010, v110, v001, v101, v011, v111, output )
      !-----------------------------------------------------------
      ! Compute derivative of a trilinear interpolation in a given
      ! direction
      ! 
      ! Params
      !     - x, y, z, doubleprecision: local cell coordinates 
      !     - direction, integer: specifies x=1, y=2 or z=3 direction
      !     - vijk , doubleprecision: values at corresponding corners
      !     - output, doubleprecision: the output variable
      !-----------------------------------------------------------
      ! Specifications
      !-----------------------------------------------------------
      implicit none
      class (TrackSubCellType) :: this
      integer :: direction
      doubleprecision :: x, y, z
      doubleprecision :: v000, v100, v010, v110, v001, v101, v011, v111
      doubleprecision :: output
      doubleprecision :: v0, v1, v00, v10, v01, v11
      !-----------------------------------------------------------

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
      ! Compute trilinear interpolation of corner velocities
      ! for the given coordinates 
      !  
      ! Params
      !     - x, y, z, doubleprecision: coordinates to interpolate
      !     - vijk , doubleprecision: values at corresponding corners
      !     - output, doubleprecision: the output variable
      !-----------------------------------------------------------
      ! Specifications
      !-----------------------------------------------------------
      implicit none
      class (TrackSubCellType) :: this
      doubleprecision :: x, y, z
      doubleprecision :: v000, v100, v010, v110, v001, v101, v011, v111
      doubleprecision :: output
      doubleprecision :: v0, v1, v00, v10, v01, v11
      !-----------------------------------------------------------

      v00    = ( 1.0d0 - x )*v000 + x*v100
      v01    = ( 1.0d0 - x )*v001 + x*v101
      v10    = ( 1.0d0 - x )*v010 + x*v110
      v11    = ( 1.0d0 - x )*v011 + x*v111
      v0     = ( 1.0d0 - y )*v00  + y*v10
      v1     = ( 1.0d0 - y )*v01  + y*v11
      output = ( 1.0d0 - z )*v0   + z*v1

  end subroutine pr_Trilinear



end module TrackSubCellModule
