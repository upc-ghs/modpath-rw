module ParticleManagerModule
  use ParticleGroupModule,only : ParticleGroupType
  use ParticleModule,only : ParticleType
  use ParticleCoordinateModule,only : ParticleCoordinateType
  use ModpathSimulationDataModule,only : ModpathSimulationDataType
  use GeoReferenceModule,only : GeoReferenceType
  use ModflowRectangularGridModule,only : ModflowRectangularGridType
#ifdef _OPENMP
  use omp_lib
#endif
  implicit none
 
  ! For parallel output in timeseries
  abstract interface
    subroutine TimeseriesWriter( sequenceNumber, particleID, groupIndex, timeStep, & 
             timePointIndex, pCoord, geoRef, outUnit, tsRecordCounts, tsTempUnits)
      !-------------------------------------------------------------
      import ParticleCoordinateType
      import GeoReferenceType
      !-------------------------------------------------------------
      type(ParticleCoordinateType),intent(in) :: pCoord
      type(GeoReferenceType),intent(in) :: geoRef
      integer,intent(in) :: outUnit, particleID, timePointIndex, & 
                            groupIndex, timeStep, sequenceNumber
      integer, dimension(:), intent(inout) :: tsRecordCounts
      integer, dimension(:) :: tsTempUnits
    end subroutine TimeseriesWriter
  end interface


  ! Write resident observation records 
  abstract interface
    subroutine ResidentObsWriter( timeStep, timePointIndex, particle, &
                      pCoord, soluteID, rFactor, waterVolume, outUnit )
      !-------------------------------------------------------------
      import ParticleType
      import ParticleCoordinateType
      !-------------------------------------------------------------
      integer,intent(in)                      :: timeStep, timePointIndex
      type(ParticleType),intent(in)           :: particle
      type(ParticleCoordinateType),intent(in) :: pCoord
      integer, intent(in)                     :: soluteID
      doubleprecision, intent(in)             :: rFactor, waterVolume 
      integer, intent(in)                     :: outUnit
    end subroutine ResidentObsWriter
  end interface


  ! Write sink observation records 
  abstract interface
    subroutine SinkObsWriter( timeStep, timePointIndex, particle, & 
                                      soluteID, flowRate, outUnit )
      !-------------------------------------------------------------
      import ParticleType
      !-------------------------------------------------------------
      ! input
      integer,intent(in)              :: timeStep, timePointIndex
      type(ParticleType),intent(in)   :: particle
      integer, intent(in)             :: soluteID
      doubleprecision, intent(in)     :: flowRate
      integer, intent(in)             :: outUnit
      !----------------------------------------------
    end subroutine SinkObsWriter
  end interface


  ! Write endpoint file
  abstract interface
    subroutine EndpointWriter(simData, grid, geoRef, outUnit)
      !-------------------------------------------------------------
      import ParticleType
      import ModflowRectangularGridType
      import ModpathSimulationDataType
      import GeoReferenceType
      !-------------------------------------------------------------
      ! input
      type(ModpathSimulationDataType), intent(in), target :: simData
      class(ModflowRectangularGridType), intent(in)       :: grid
      type(GeoReferenceType), intent(in)                  :: geoRef
      integer, intent(in)                                 :: outUnit
      !-------------------------------------------------------------
    end subroutine EndpointWriter
  end interface

contains

  
  subroutine WriteEndpoints(simData, grid, geoRef, outUnit)
  !--------------------------------------------------------------------------------
  ! Specifications
  !--------------------------------------------------------------------------------
  implicit none
  integer,intent(in) :: outUnit
  integer :: version, subversion, pgIndex, pIndex, n, face, zone, initialZone,  &
    totalCount, releaseCount, maximumID
  doubleprecision :: initialGlobalZ, globalZ, initialModelX, initialModelY, modelX, modelY
  integer,dimension(0:9) :: statusSums
  type(ParticleType),pointer :: p
  class(ModflowRectangularGridType),intent(in) :: grid
  type(ModpathSimulationDataType),intent(in),target :: simData
  type(GeoReferenceType),intent(in) :: geoRef
  !--------------------------------------------------------------------------------
  
    version = 7
    subversion = 2
    
    ! Compute status sums for all particles
    totalCount = 0
    releaseCount = 0
    maximumID = 0
    do n = 0, 9
        statusSums(n) = 0
    end do
    
    do pgIndex = 1, simData%ParticleGroupCount
      do pIndex = 1, simData%ParticleGroups(pgIndex)%TotalParticleCount
        totalCount = totalCount + 1
        p => simData%ParticleGroups(pgIndex)%Particles(pIndex)
        if((p%Status .ge. 0) .and. (p%Status .le. 9)) then
          statusSums(p%Status) = statusSums(p%Status) + 1
          if((p%Status .ne. 0) .and. (p%Status .ne. 8) .and. (p%ID .gt. maximumID)) maximumID = p%ID
        end if
      end do
    end do
    releaseCount = totalCount - statusSums(0) - statusSums(8)
    
    write(outUnit, '(a,2i10)') 'MODPATH_ENDPOINT_FILE', version, subversion
    write(outUnit, '(i5,3i10,1x,4e18.10)') simData%TrackingDirection, totalCount, releaseCount, maximumID, &
        simData%ReferenceTime, grid%OriginX, grid%OriginY, grid%RotationAngle
    write(outUnit, '(10i8)') (statusSums(n), n = 0, 9)
    write(outUnit, '(i10)') simData%ParticleGroupCount
    do n = 1, simData%ParticleGroupCount
        write(outUnit, '(a)') simData%ParticleGroups(n)%Name
    end do
    write(outUnit, '(a)') 'END HEADER'
    
    do pgIndex = 1, simData%ParticleGroupCount
      do pIndex = 1, simData%ParticleGroups(pgIndex)%TotalParticleCount
        p => simData%ParticleGroups(pgIndex)%Particles(pIndex)
        if((p%Status .ne. 0) .and. (p%Status .ne. 8) ) then
          call grid%ConvertToModelXYZ(p%InitialCellNumber,                 &
            p%InitialLocalX, p%InitialLocalY, p%InitialLocalZ,             &
            initialModelX, initialModelY, initialGlobalZ)
          call grid%ConvertToModelXYZ(p%CellNumber, p%LocalX, p%LocalY,    &
            p%LocalZ, modelX, modelY, globalZ)
          face = FindFace(p%LocalX, p%LocalY, p%LocalZ)
          initialZone = simData%Zones(p%InitialCellNumber)
          zone = simData%Zones(p%CellNumber)
          write(outUnit, '(3i10,i5,2es18.9e3,i10,i5,6es18.9e3,2i5,i10,i5,6es18.9e3,2i5)') &
            p%SequenceNumber, p%Group, p%ID, p%Status,                    &
            p%InitialTrackingtime, p%TrackingTime, p%InitialCellNumber,   &
            p%InitialLayer, p%InitialLocalX, p%InitialLocalY,             &
            p%InitialLocalZ, initialModelX, initialModelY,                &
            initialGlobalZ, initialZone, p%InitialFace, p%CellNumber,     &
            p%Layer, p%LocalX, p%LocalY, p%LocalZ, modelX, modelY,        &
            globalZ, zone, face
        end if
      end do
    end do
  
  
  end subroutine WriteEndpoints


  subroutine WriteEndpointsMassParticles(simData, grid, geoRef, outUnit)
  !--------------------------------------------------------------------------------
  ! Is a modified form of WriteEndpoints which includes a solute id 
  ! and the particle mass. It removes the local cell coordinates and the face
  ! locations
  !
  !--------------------------------------------------------------------------------
  ! Specifications
  !--------------------------------------------------------------------------------
  implicit none
  integer,intent(in) :: outUnit
  integer :: version, subversion, pgIndex, pIndex, n, zone, initialZone,  &
    totalCount, releaseCount, maximumID
  !integer :: face
  doubleprecision :: initialGlobalZ, globalZ, initialModelX, initialModelY, modelX, modelY
  integer,dimension(0:9) :: statusSums
  type(ParticleType),pointer :: p
  class(ModflowRectangularGridType),intent(in) :: grid
  type(ModpathSimulationDataType),intent(in),target :: simData
  type(GeoReferenceType),intent(in) :: geoRef
  !--------------------------------------------------------------------------------
  
    version = 1
    subversion = 0
    
    ! Compute status sums for all particles
    totalCount = 0
    releaseCount = 0
    maximumID = 0
    do n = 0, 9
        statusSums(n) = 0
    end do
    
    do pgIndex = 1, simData%ParticleGroupCount
      do pIndex = 1, simData%ParticleGroups(pgIndex)%TotalParticleCount
        totalCount = totalCount + 1
        p => simData%ParticleGroups(pgIndex)%Particles(pIndex)
        if((p%Status .ge. 0) .and. (p%Status .le. 9)) then
          statusSums(p%Status) = statusSums(p%Status) + 1
          if((p%Status .ne. 0) .and. (p%Status .ne. 8) .and. (p%ID .gt. maximumID)) maximumID = p%ID
        end if
      end do
    end do
    releaseCount = totalCount - statusSums(0) - statusSums(8)
    
    write(outUnit, '(a,2i10)') 'MODPATHRW_ENDPOINT_FILE', version, subversion
    write(outUnit, '(i5,3i10,1x,4e18.10)') simData%TrackingDirection, totalCount, releaseCount, maximumID, &
        simData%ReferenceTime, grid%OriginX, grid%OriginY, grid%RotationAngle
    write(outUnit, '(10i8)') (statusSums(n), n = 0, 9)
    write(outUnit, '(i10)') simData%ParticleGroupCount
    do n = 1, simData%ParticleGroupCount
        write(outUnit, '(a)') simData%ParticleGroups(n)%Name
    end do
    write(outUnit, '(a)') 'END HEADER'
    
    do pgIndex = 1, simData%ParticleGroupCount
      do pIndex = 1, simData%ParticleGroups(pgIndex)%TotalParticleCount
        p => simData%ParticleGroups(pgIndex)%Particles(pIndex)
        if((p%Status .ne. 0) .and. (p%Status .ne. 8) ) then
          call grid%ConvertToModelXYZ(p%InitialCellNumber,                 &
            p%InitialLocalX, p%InitialLocalY, p%InitialLocalZ,             &
            initialModelX, initialModelY, initialGlobalZ)
          call grid%ConvertToModelXYZ(p%CellNumber, p%LocalX, p%LocalY,    &
            p%LocalZ, modelX, modelY, globalZ)
          !face = FindFace(p%LocalX, p%LocalY, p%LocalZ)
          initialZone = simData%Zones(p%InitialCellNumber)
          zone = simData%Zones(p%CellNumber)
          write(outUnit, '(4i10,i5,3es18.9e3,i10,i5,3es18.9e3,1i5,i10,i5,3es18.9e3,1i5)')            &
            p%SequenceNumber, p%Group, simData%ParticleGroups(pgIndex)%Solute, p%ID, p%Status,       &
            p%Mass, p%InitialTrackingtime, p%TrackingTime, p%InitialCellNumber,                      &
            p%InitialLayer, initialModelX, initialModelY, initialGlobalZ, initialZone, p%CellNumber, &
            p%Layer, modelX, modelY, globalZ, zone
        end if
      end do
    end do
  
  
  end subroutine WriteEndpointsMassParticles


  subroutine WriteTimeseriesHeader(outUnit, trackingDirection, referenceTime,   &
      originX, originY, rotationAngle)
  !--------------------------------------------------------------------------------
  ! Specifications
  !--------------------------------------------------------------------------------
  implicit none
  integer,intent(in) :: outUnit, trackingDirection
  doubleprecision,intent(in) :: referenceTime, originX, originY, rotationAngle
  integer :: version, subversion
  !--------------------------------------------------------------------------------
  
    version = 7
    subversion = 2
    write(outUnit, '(a,2i10)') 'MODPATH_TIMESERIES_FILE', version, subversion
    write(outUnit, '(i5,1x,4e18.10)') trackingDirection, referenceTime, originX, originY, & 
        rotationAngle
    write(outUnit, '(a)') 'END HEADER'
  
  end subroutine WriteTimeseriesHeader
 

  subroutine WriteTimeseriesRecord(sequenceNumber, particleID, groupIndex,      &
    timeStep, timePointIndex, pCoord, geoRef, outUnit)
  !--------------------------------------------------------------------------------
  ! Specifications
  !--------------------------------------------------------------------------------
  implicit none
  type(ParticleCoordinateType),intent(in) :: pCoord
  type(GeoReferenceType),intent(in) :: geoRef
  integer,intent(in) :: outUnit, particleID, timePointIndex, groupIndex,        &
    timeStep, sequenceNumber
  doubleprecision :: modelX, modelY
  !--------------------------------------------------------------------------------
  
    modelX = pCoord%GlobalX
    modelY = pCoord%GlobalY
    
    write(outUnit, '(2I8,es18.9e3,i10,i5,2i10,6es18.9e3,i10)')                    &
      timePointIndex, timeStep, pCoord%TrackingTime, sequenceNumber, groupIndex,  &
      particleID, pCoord%CellNumber, pCoord%LocalX, pCoord%LocalY, pCoord%LocalZ, &
      modelX, modelY, pCoord%GlobalZ, pCoord%Layer

  end subroutine WriteTimeseriesRecord


  subroutine WriteTimeseriesRecordStatus(sequenceNumber, particleID, groupIndex,      &
    timeStep, timePointIndex, pCoord, geoRef, particleStatus, outUnit)
  !--------------------------------------------------------------------------------
  ! Specifications
  !--------------------------------------------------------------------------------
  implicit none
  type(ParticleCoordinateType),intent(in) :: pCoord
  type(GeoReferenceType),intent(in) :: geoRef
  integer,intent(in) :: outUnit, particleID, timePointIndex, groupIndex,        &
    timeStep, sequenceNumber, particleStatus
  doubleprecision :: modelX, modelY
  !--------------------------------------------------------------------------------
  
    modelX = pCoord%GlobalX
    modelY = pCoord%GlobalY
    
    write(outUnit, '(2I8,es18.9e3,i10,i5,2i10,6es18.9e3,2i10)')                   &
      timePointIndex, timeStep, pCoord%TrackingTime, sequenceNumber, groupIndex,  &
      particleID, pCoord%CellNumber, pCoord%LocalX, pCoord%LocalY, pCoord%LocalZ, &
      modelX, modelY, pCoord%GlobalZ, pCoord%Layer, particleStatus
  
  end subroutine WriteTimeseriesRecordStatus


  subroutine WriteBinaryTimeseriesRecordId(sequenceNumber, particleID, groupIndex, &
    timeStep, timePointIndex, pCoord, geoRef, recordID, outUnit)
  !--------------------------------------------------------------------------------
  ! Write timeseries record plus recordID for consolidation 
  !--------------------------------------------------------------------------------
  implicit none
  type(ParticleCoordinateType),intent(in) :: pCoord
  type(GeoReferenceType),intent(in) :: geoRef
  integer,intent(in) :: outUnit, particleID, timePointIndex, groupIndex,&
    timeStep, sequenceNumber
  integer, intent(in) :: recordID
  doubleprecision :: modelX, modelY
  !--------------------------------------------------------------------------------
    
    modelX = pCoord%GlobalX
    modelY = pCoord%GlobalY
    
    !'(2I8,es18.9e3,i10,i5,2i10,6es18.9e3,i10)+recordID'
    write(outUnit) &
      timePointIndex, timeStep, pCoord%TrackingTime, sequenceNumber, groupIndex,  &
      particleID, pCoord%CellNumber, pCoord%LocalX, pCoord%LocalY, pCoord%LocalZ, &
      modelX, modelY, pCoord%GlobalZ, pCoord%Layer, recordID
  
  end subroutine WriteBinaryTimeseriesRecordId


  subroutine WriteTimeseriesRecordSerial(sequenceNumber, particleID, groupIndex, timeStep, &
                      timePointIndex, pCoord, geoRef, outUnit, tsRecordCounts, tsTempUnits )
  !--------------------------------------------------------------------------------
  ! WriteTimeseriesRecord classical, this method is only for interface consistency
  !--------------------------------------------------------------------------------
  implicit none
  type(ParticleCoordinateType),intent(in) :: pCoord
  type(GeoReferenceType),intent(in) :: geoRef
  integer,intent(in) :: outUnit, particleID, timePointIndex, groupIndex,&
    timeStep, sequenceNumber
  integer, dimension(:), intent(inout) :: tsRecordCounts
  integer, dimension(:) :: tsTempUnits
  !--------------------------------------------------------------------------------

    call WriteTimeseriesRecord(sequenceNumber, particleID, groupIndex, &
                      timeStep, timePointIndex, pCoord, geoRef, outUnit)


  end subroutine WriteTimeseriesRecordSerial


  subroutine WriteTimeseriesRecordCritical(sequenceNumber, particleID, groupIndex, timeStep, &
                        timePointIndex, pCoord, geoRef, outUnit, tsRecordCounts, tsTempUnits )
  !--------------------------------------------------------------------------------
  ! WriteTimeseriesRecord with OpenMP critical directive
  !--------------------------------------------------------------------------------
  implicit none
  type(ParticleCoordinateType),intent(in) :: pCoord
  type(GeoReferenceType),intent(in) :: geoRef
  integer,intent(in) :: outUnit, particleID, timePointIndex, groupIndex, &
    timeStep, sequenceNumber
  integer, dimension(:), intent(inout) :: tsRecordCounts
  integer, dimension(:) :: tsTempUnits
  !--------------------------------------------------------------------------------


    !$omp critical( timeseries )
    call WriteTimeseriesRecord(sequenceNumber, particleID, groupIndex, &
                      timeStep, timePointIndex, pCoord, geoRef, outUnit)
    !$omp end critical( timeseries )


  end subroutine WriteTimeseriesRecordCritical


  subroutine WriteTimeseriesRecordConsolidate(sequenceNumber, particleID, groupIndex, timeStep, & 
                           timePointIndex, pCoord, geoRef, outUnit, tsRecordCounts, tsTempUnits )
  !--------------------------------------------------------------------------------
  ! WriteBinaryTimeseriesRecordId for each parallel unit 
  !--------------------------------------------------------------------------------
  implicit none
  type(ParticleCoordinateType),intent(in) :: pCoord
  type(GeoReferenceType),intent(in) :: geoRef
  integer,intent(in) :: outUnit, particleID, timePointIndex, groupIndex, &
    timeStep, sequenceNumber
  integer, dimension(:), intent(inout) :: tsRecordCounts
  integer, dimension(:) :: tsTempUnits
  integer :: ompThreadId
  !--------------------------------------------------------------------------------

#ifdef _OPENMP
    ompThreadId = omp_get_thread_num() + 1 ! Starts at zero
#else
    ompThreadId = 1
#endif

    tsRecordCounts( ompThreadId ) = tsRecordCounts( ompThreadId ) + 1
    call WriteBinaryTimeseriesRecordId(sequenceNumber, particleID, groupIndex, & 
      timeStep, timePointIndex, pCoord, geoRef, tsRecordCounts( ompThreadId ), &
      tsTempUnits( ompThreadId ))

  end subroutine WriteTimeseriesRecordConsolidate


  subroutine WriteTimeseriesRecordThread(sequenceNumber, particleID, groupIndex, timeStep, & 
                      timePointIndex, pCoord, geoRef, outUnit, tsRecordCounts, tsTempUnits )
  !--------------------------------------------------------------------------------
  ! WriteBinaryTimeseriesRecord for each parallel unit 
  !--------------------------------------------------------------------------------
  implicit none
  type(ParticleCoordinateType),intent(in) :: pCoord
  type(GeoReferenceType),intent(in) :: geoRef
  integer,intent(in) :: outUnit, particleID, timePointIndex, groupIndex, &
    timeStep, sequenceNumber
  integer, dimension(:), intent(inout) :: tsRecordCounts
  integer, dimension(:) :: tsTempUnits
  integer :: ompThreadId
  !--------------------------------------------------------------------------------

#ifdef _OPENMP
    ompThreadId = omp_get_thread_num() + 1 ! Starts at zero
#else
    ompThreadId = 1
#endif

    call WriteTimeseriesRecord(sequenceNumber, particleID, groupIndex, timeStep, &
                      timePointIndex, pCoord, geoRef, tsTempUnits( ompThreadId ) )


  end subroutine WriteTimeseriesRecordThread


  subroutine ConsolidateParallelTimeseriesRecords(inUnits, outUnit, recordCounts, lastRecord)
  !--------------------------------------------------------------------------------------
  ! Read data from binary temporal inUnits and consolidate timeseries into outUnit
  !--------------------------------------------------------------------------------------
  implicit none
  integer, dimension(:), intent(in) :: inUnits, recordCounts
  integer, intent(in)    :: outUnit
  integer, intent(inout) :: lastRecord
  type(ParticleCoordinateType) :: pCoord
  integer  :: timePointIndex, timeStep, sequenceNumber,  groupIndex, particleID
  doubleprecision :: modelX, modelY
  integer  :: n, i
  integer  :: nThreads, recordID
  integer  :: startFromRecord(size(recordCounts)+1)
  !--------------------------------------------------------------------------------------

    nThreads = size( inUnits ) 

    !! ALL THREADS PROCESS THE SAME UNIT 
    !! Read from temporal units and dump
    !do n = 1, nThreads
    !  if( recordCounts(n) .ge. 1 ) then
    !      rewind( inUnits(n) )
    !      !$omp parallel do schedule(dynamic,1)                         &
    !      !$omp default( none )                                         &
    !      !$omp shared( n, inUnits, recordCounts, outUnit, lastRecord ) &
    !      !$omp private( timePointIndex, timeStep, pCoord )             &
    !      !$omp private( sequenceNumber, groupIndex, particleID )       &
    !      !$omp private( modelX, modelY, recordID, reclen, reclencumm )
    !      do i = 1, recordCounts(n)
    !          ! Read from temporal binary
    !          read( inUnits(n) ) &
    !            timePointIndex, timeStep, pCoord%TrackingTime, sequenceNumber, groupIndex,   &
    !            particleID, pCoord%CellNumber, pCoord%LocalX, pCoord%LocalY, pCoord%LocalZ,  &
    !            modelX, modelY, pCoord%GlobalZ, pCoord%Layer, recordID
    !          ! Write to consolidated direct access unit
    !          write(outUnit, rec=lastRecord+recordID, fmt='(2I8,es18.9e3,i10,i5,2i10,6es18.9e3,i10,a1)') & 
    !            timePointIndex, timeStep, pCoord%TrackingTime, sequenceNumber, groupIndex,               &
    !            particleID, pCoord%CellNumber, pCoord%LocalX, pCoord%LocalY, pCoord%LocalZ,              &
    !            modelX, modelY, pCoord%GlobalZ, pCoord%Layer, char(10)
    !      end do
    !      !$omp end parallel do
    !      lastRecord = lastRecord + recordCounts(n)
    !  end if
    !end do
 
    ! EACH THREAD PROCESS ITS OWN UNIT
    startFromRecord     = 0
    startFromRecord(2:) = [ (sum(recordCounts(1:i)), i=1, size(recordCounts)) ]
    startFromRecord     = startFromRecord + lastRecord 
    ! Read from temporal units and dump
    !$omp parallel do                                             &
    !$omp default( none )                                         &
    !$omp shared( inUnits, recordCounts, outUnit, lastRecord )    &
    !$omp shared( startFromRecord, nThreads )                     &
    !$omp private( timePointIndex, timeStep, pCoord )             &
    !$omp private( sequenceNumber, groupIndex, particleID )       &
    !$omp private( modelX, modelY, recordID )
    do n = 1, nThreads
      if( recordCounts(n) .ge. 1 ) then
          rewind( inUnits(n) )
          do i = 1, recordCounts(n)
              ! Read from temporal binary
              read( inUnits(n) ) &
                timePointIndex, timeStep, pCoord%TrackingTime, sequenceNumber, groupIndex,   &
                particleID, pCoord%CellNumber, pCoord%LocalX, pCoord%LocalY, pCoord%LocalZ,  &
                modelX, modelY, pCoord%GlobalZ, pCoord%Layer, recordID
              ! Write to consolidated direct access unit
              write(outUnit, rec=startFromRecord(n)+recordID, fmt='(2I8,es18.9e3,i10,i5,2i10,6es18.9e3,i10,a1)') & 
                timePointIndex, timeStep, pCoord%TrackingTime, sequenceNumber, groupIndex,               &
                particleID, pCoord%CellNumber, pCoord%LocalX, pCoord%LocalY, pCoord%LocalZ,              &
                modelX, modelY, pCoord%GlobalZ, pCoord%Layer, char(10)
          end do
      end if
    end do
    !$omp end parallel do

    lastRecord = startFromRecord(size(recordCounts)+1)

  end subroutine ConsolidateParallelTimeseriesRecords


  subroutine WritePathlineHeader(outUnit, trackingDirection, referenceTime,     &
      originX, originY, rotationAngle)
  !---------------------------------------------------------------------------------
  ! Specifications
  !---------------------------------------------------------------------------------
  implicit none
  integer,intent(in)  :: outUnit, trackingDirection
  doubleprecision,intent(in) :: referenceTime, originX, originY, rotationAngle
  integer :: version, subversion
  !---------------------------------------------------------------------------------

    version = 7
    subversion = 2
    write(outUnit, '(a,2i10)') 'MODPATH_PATHLINE_FILE', version, subversion
    write(outUnit, '(i5,1x,4e18.10)') trackingDirection, referenceTime, originX, originY, &
        rotationAngle
    write(outUnit, '(a)') 'END HEADER'
  
  end subroutine WritePathlineHeader
 

  subroutine WritePathlineRecord(tpResult, outUnit, stressPeriod, timeStep, geoRef)
  use TrackPathResultModule,only : TrackPathResultType
  use ParticleCoordinateModule,only : ParticleCoordinateType
  !---------------------------------------------------------------------------------
  ! Specifications
  !---------------------------------------------------------------------------------
  implicit none
  integer,intent(in) :: outUnit, stressPeriod, timeStep
  type(TrackPathResultType),intent(in),target :: tpResult
  type(ParticleCoordinateType),pointer :: c
  type(GeoReferenceType) :: geoRef
  integer  :: n, count
  doubleprecision :: modelX, modelY
  !---------------------------------------------------------------------------------
  
    count = tpResult%ParticlePath%Pathline%GetItemCount()
    if(count .lt. 2) return
    
    write(outUnit, '(4i10)') tpResult%SequenceNumber, tpResult%Group, tpResult%ParticleID, count
    do n = 1, count
      c => tpResult%ParticlePath%Pathline%Items(n)
      modelX = c%GlobalX
      modelY = c%GlobalY
      write(outUnit, "(i10,7es16.7e3,3I10)")                            &
        c%CellNumber, modelX, modelY, c%GlobalZ, c%TrackingTime,        &
        c%LocalX, c%LocalY, c%LocalZ, c%Layer, stressPeriod, timeStep
    end do
 

  end subroutine WritePathlineRecord
 

  subroutine WriteBinaryPathlineRecord(tpResult, outUnit, stressPeriod, timeStep, geoRef)
  use TrackPathResultModule,only : TrackPathResultType
  use ParticleCoordinateModule,only : ParticleCoordinateType
  !---------------------------------------------------------------------------------
  ! Specifications
  !---------------------------------------------------------------------------------
  implicit none
  integer,intent(in) :: outUnit, stressPeriod, timeStep
  type(TrackPathResultType),intent(in),target :: tpResult
  type(ParticleCoordinateType),pointer :: c
  type(GeoReferenceType) :: geoRef
  integer  :: n, count, currentPosition
  doubleprecision :: modelX, modelY
  !---------------------------------------------------------------------------------
  
    count = tpResult%ParticlePath%Pathline%GetItemCount()
    if(count .lt. 2) return

    inquire(unit=outUnit, pos=currentPosition)
    write(outUnit) tpResult%SequenceNumber, count, tpResult%Group, tpResult%ParticleID
    do n = 1, count
      c => tpResult%ParticlePath%Pathline%Items(n)
      modelX = c%GlobalX
      modelY = c%GlobalY
      write(outUnit) c%CellNumber, modelX, modelY, c%GlobalZ,           &
        c%TrackingTime, c%LocalX, c%LocalY, c%LocalZ, c%Layer, stressPeriod,  &
        timeStep
    end do
  

  end subroutine WriteBinaryPathlineRecord
  

  subroutine ConsolidatePathlines(inUnit, outUnit, recordCount, particleCount)
  !---------------------------------------------------------------------------------
  ! Specifications
  !---------------------------------------------------------------------------------
  implicit none
  integer, intent(in) :: inUnit, outUnit, recordCount, particleCount
  integer, dimension(:), allocatable :: sequenceNumbers, recordPointCounts,     &
    particlePointCounts, particleRecordCounts
  integer (kind=8), dimension(:), allocatable :: recordPointers
  integer :: dataOffset, n, m, i, count, group, id,               &
    stressPeriod, timeStep
  integer (kind=8) :: pos, ptr
  integer  :: cellNumber, layer, pointCount
  doubleprecision :: modelX, modelY, globalZ, trackingTime, localX, localY, localZ
  !---------------------------------------------------------------------------------
  
    allocate(sequenceNumbers(recordCount))
    allocate(recordPointCounts(recordCount))
    allocate(recordPointers(recordCount))
    allocate(particlePointCounts(particleCount))
    allocate(particleRecordCounts(particleCount))
    
    ! Initialize all particle point count and record count array elements to 0
    do n = 1, particleCount
      particlePointCounts(n) = 0
      particleRecordCounts(n) = 0
    end do
    
    do n = 1, recordCount
      sequenceNumbers(n) = 0
      recordPointCounts(n) = 0
      recordPointers(n) = 0
    end do 
   
    ! Read and save pathline record header index information
    pos = 1
    do n = 1, recordCount
      recordPointers(n) = pos
      read(inUnit, pos=pos) sequenceNumbers(n), recordPointCounts(n)
      m = sequenceNumbers(n)
      if(m .gt. 0) then
        particlePointCounts(m) = particlePointCounts(m) + recordPointCounts(n)
        particleRecordCounts(m) = particleRecordCounts(m) + 1
      end if
      dataOffset = 16 + (recordPointCounts(n)*72)
      pos = recordPointers(n) + dataOffset
    end do
    
    ! Loop through particles. For each particle, loop through the pathline record headers and read and 
    ! consolidate pathline points for all of the segments belonging to that particle. Then write the 
    ! consolidated pathline for the particle. Skip particles for which the particlePointCount equals 0.
    do n = 1, particleCount
      if(particlePointCounts(n) .gt. 0) then
        count = 0
        do m = 1, recordCount
          if(sequenceNumbers(m) .eq. n) then
            ! Increment the record count for this particle
            count = count + 1
            
            ! If this is the first pathline record for this particle, write the record header for the
            ! consolidated pathline.
            if(count .eq. 1) then
              ! Read the rest of the record header
              ptr = recordPointers(m) + 8
              read(inUnit, pos=ptr) group, id
              pointCount = particlePointCounts(n) - ParticleRecordCounts(n) + 1
              write(outUnit, '(4i10)') n, group, id, pointCount
            else
              ptr = recordPointers(m) + 16
              read(inUnit, pos=ptr)
            end if
            
            do i = 1, recordPointCounts(m)
              read(inUnit) cellNumber, modelX, modelY, globalZ, &
                trackingTime, localX, localY, localZ, layer,    &
                stressPeriod, timeStep
              if((count .gt. 1) .and. (i .eq. 1)) cycle
              write(outUnit, "(i10,7es16.7e3,3I10)") cellNumber, modelX,&
                modelY, globalZ, trackingTime, localX, localY, localZ,  &
                layer, stressPeriod, timeStep
            end do
          end if
          ! Check to see if this was the last pathline record for this particle. If so, exit the record loop because there is no need 
          ! to read through the rest of the records for this particle.
          if(count .eq. particleRecordCounts(n)) exit
        end do
      end if
    end do
  
  end subroutine ConsolidatePathlines


  function FindFace(x, y, z) result(face)
  !---------------------------------------------------------------------------------
  ! Specifications
  !---------------------------------------------------------------------------------
  implicit none
  doubleprecision, intent(in) :: x, y, z
  integer :: face
  !---------------------------------------------------------------------------------
  
    face = 0
    if(x .eq. 0d0) then
      face = 1
    else if(x .eq. 1d0) then
      face = 2
    else if(y .eq. 0d0) then
      face = 3
    else if(y .eq. 1d0) then
      face = 4
    else if(z .eq. 0d0) then
      face = 5
    else if(z .eq. 1d0) then
      face = 6
    end if
 

  end function FindFace


  subroutine WriteResidentObsRecord(timeStep, timePointIndex, particle, pCoord, &
                                        soluteID, rFactor, waterVolume, outUnit )
  !--------------------------------------------------------------------------------
  ! Specifications
  !---------------------------------------------------------------------------------
  implicit none
  integer,intent(in)                      :: timeStep, timePointIndex
  type(ParticleType),intent(in)           :: particle
  type(ParticleCoordinateType),intent(in) :: pCoord
  integer, intent(in)                     :: soluteID
  doubleprecision, intent(in)             :: rFactor, waterVolume 
  integer, intent(in)                     :: outUnit
  doubleprecision :: modelX, modelY
  !---------------------------------------------------------------------------------
  
    modelX = pCoord%GlobalX
    modelY = pCoord%GlobalY
    
    write(outUnit, '(2I8,es18.9e3,i10,es18.9e3,2i5,2i10,5es18.9e3)')             &
      timePointIndex, timeStep, pCoord%TrackingTime, particle%ID, particle%Mass, & 
                      particle%Group, soluteID, pCoord%CellNumber, pCoord%Layer, &
                           rFactor, waterVolume, modelX, modelY, pCoord%GlobalZ


  end subroutine WriteResidentObsRecord


  subroutine WriteResidentObsRecordBinary(timeStep, timePointIndex, particle, pCoord, &
                                              soluteID, rFactor, waterVolume, outUnit )
  !--------------------------------------------------------------------------------
  ! Specifications
  !---------------------------------------------------------------------------------
  implicit none
  integer,intent(in)                      :: timeStep, timePointIndex
  type(ParticleType),intent(in)           :: particle
  type(ParticleCoordinateType),intent(in) :: pCoord
  integer, intent(in)                     :: soluteID
  doubleprecision, intent(in)             :: rFactor, waterVolume 
  integer, intent(in)                     :: outUnit
  doubleprecision :: modelX, modelY
  !---------------------------------------------------------------------------------
  
    modelX = pCoord%GlobalX
    modelY = pCoord%GlobalY
    
    write(outUnit) &
      timePointIndex, timeStep, pCoord%TrackingTime, particle%ID, particle%Mass, & 
                      particle%Group, soluteID, pCoord%CellNumber, pCoord%Layer, &
                           rFactor, waterVolume, modelX, modelY, pCoord%GlobalZ


  end subroutine WriteResidentObsRecordBinary


  subroutine WriteSinkObsRecord( timeStep, timePointIndex, particle, &
                                         soluteID, flowRate, outUnit )
  !---------------------------------------------------------------------------------
  ! Write observation sink cell record
  !---------------------------------------------------------------------------------
  ! Specifications
  !---------------------------------------------------------------------------------
  implicit none
  ! input
  integer,intent(in)              :: timeStep, timePointIndex
  type(ParticleType),intent(in)   :: particle
  integer, intent(in)             :: soluteID
  doubleprecision, intent(in)     :: flowRate
  integer, intent(in)             :: outUnit
  !---------------------------------------------------------------------------------


    write(outUnit, '(2I8,es18.9e3,i10,es18.9e3,2i5,2i10,es18.9e3)')                &
      timePointIndex, timeStep, particle%TrackingTime, particle%ID, particle%Mass, & 
        particle%Group, soluteID, particle%CellNumber, particle%Layer, flowRate


  end subroutine WriteSinkObsRecord


  subroutine WriteSinkObsRecordBinary( timeStep, timePointIndex, particle, & 
                                               soluteID, flowRate, outUnit )
  !---------------------------------------------------------------------------------
  ! Write observation sink cell record
  !---------------------------------------------------------------------------------
  ! Specifications
  !---------------------------------------------------------------------------------
  implicit none
  ! input
  integer,intent(in)              :: timeStep, timePointIndex
  type(ParticleType),intent(in)   :: particle
  integer, intent(in)             :: soluteID
  doubleprecision, intent(in)     :: flowRate
  integer, intent(in)             :: outUnit
  !---------------------------------------------------------------------------------


    write(outUnit) &
      timePointIndex, timeStep, particle%TrackingTime, particle%ID, particle%Mass, & 
        particle%Group, soluteID, particle%CellNumber, particle%Layer, flowRate


  end subroutine WriteSinkObsRecordBinary


end module ParticleManagerModule
