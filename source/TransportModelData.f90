module TransportModelDataModule
  use ModflowRectangularGridModule,only : ModflowRectangularGridType
  use ParticleTrackingOptionsModule,only : ParticleTrackingOptionsType
  use ParticleGroupModule,only : ParticleGroupType
  use StartingLocationReaderModule,only : ReadAndPrepareLocations,    &
                                   pr_CreateParticlesAsInternalArray, &
                                   CreateMassParticlesAsInternalArray
  use FlowModelDataModule,only : FlowModelDataType
  use ModpathSimulationDataModule, only: ModpathSimulationDataType
  implicit none
  !---------------------------------------------------------------------------------------------------------------

  ! Set default access status to private
  private

    type,public :: TransportModelDataType

      logical :: Initialized = .false.
      doubleprecision,dimension(:),allocatable :: AlphaLong
      doubleprecision,dimension(:),allocatable :: AlphaTrans
      integer,allocatable,dimension(:)         :: ICBound
      integer,allocatable,dimension(:)         :: ICBoundTS
      doubleprecision :: DMol

      ! grid
      class(ModflowRectangularGridType),pointer :: Grid => null()

      ! Private variables
      integer,private :: CurrentStressPeriod = 0
      integer,private :: CurrentTimeStep = 0

    contains

      procedure :: Initialize=>pr_Initialize
      procedure :: Reset=>pr_Reset
      procedure :: ReadData=>pr_ReadData
      procedure :: LoadTimeStep=>pr_LoadTimeStep

    end type


contains


    subroutine pr_Initialize( this, grid )
    !***************************************************************************************************************
    !
    !***************************************************************************************************************
    ! Specifications
    !---------------------------------------------------------------------------------------------------------------
    implicit none
    class(TransportModelDataType) :: this
    class(ModflowRectangularGridType),intent(inout),pointer :: grid
    integer :: cellCount, gridType
    !---------------------------------------------------------------------------------------------------------------

        this%Initialized = .false.
        
        ! Call Reset to make sure that all arrays are initially unallocated
        call this%Reset()
        
        ! Return if the grid cell count equals 0
        cellCount = grid%CellCount
        if(cellCount .le. 0) return
        
        ! Check budget reader and grid data for compatibility and allocate appropriate cell-by-cell flow arrays
        gridType = grid%GridType
        this%Grid => grid

        if(allocated(this%AlphaLong)) deallocate(this%AlphaLong)
        if(allocated(this%AlphaTrans)) deallocate(this%AlphaTrans)
        if(allocated(this%ICBoundTS)) deallocate(this%ICBoundTS)
        if(allocated(this%ICBound)) deallocate(this%ICBound)
        allocate(this%AlphaTrans(grid%CellCount))
        allocate(this%AlphaLong(grid%CellCount))
        allocate(this%ICBoundTS(grid%CellCount))
        allocate(this%ICBound(grid%CellCount))


        this%Initialized = .true.


    end subroutine pr_Initialize


    subroutine pr_ReadData(this, inUnit, inFile, outUnit, simulationData, flowModelData, ibound, grid, trackingOptions )
    !***************************************************************************************************************
    !
    !***************************************************************************************************************
    ! Specifications
    !---------------------------------------------------------------------------------------------------------------
    use utl7module,only : urdcom, upcase
    use UTL8MODULE,only : urword, ustop, u1dint, u1drel, u1ddbl, u8rdcom, &
      u3dintmp, u3dintmpusg, u3ddblmp, u3ddblmpusg
    implicit none
    class(TransportModelDataType) :: this
    class(FlowModelDataType), intent(in) :: flowModelData
    class(ModpathSimulationDataType), intent(inout) :: simulationData
    integer,intent(in) :: inUnit, outUnit
    character(len=200), intent(in)    :: inFile 
    type(ParticleTrackingOptionsType),intent(inout) :: trackingOptions
    integer,dimension(:),allocatable :: cellsPerLayer
    integer,dimension(:),allocatable :: layerTypes
    !class(ModflowRectangularGridType),intent(inout) :: grid
    class(ModflowRectangularGridType),intent(in) :: grid
    integer,dimension(grid%CellCount),intent(in) :: ibound
    character(len=200) :: line
    character(len=16) :: txt
    integer :: n, m, nn, o, p, length, iface, errorCode, layer
    integer :: icol, istart, istop
    integer :: iodispersion = 0
    doubleprecision :: r
    character(len=24),dimension(4) :: aname
    data aname(1) /'          BOUNDARY ARRAY'/
    data aname(2) /'            DISPERSIVITY'/
    data aname(3) /'            ICBOUND'/
    data aname(4) /'            IC'/
    integer :: tempAlphaUnit = 666
    character(len=200) :: tempAlphaFile

    ! Initial conditions
    integer :: nInitialConditions, particleCount, seqNumber, slocUnit
    integer :: releaseOption, releaseTimeCount
    doubleprecision :: initialReleaseTime, releaseInterval
    doubleprecision,dimension(:),allocatable :: releaseTimes
    doubleprecision :: frac, tinc
    type(ParticleGroupType),dimension(:),allocatable :: particleGroups
    type(ParticleGroupType),dimension(:),allocatable :: newParticleGroups
    integer :: initialConditionFormat, massProcessingFormat
    doubleprecision, dimension(:), allocatable :: densityDistribution
    doubleprecision, dimension(:), allocatable :: rawNParticles
    doubleprecision, dimension(:), allocatable :: nParticles
    doubleprecision, dimension(:), allocatable :: cellVolumes
    doubleprecision, dimension(:), allocatable :: delZ
    ! FOR DETERMINATION OF SUBDIVISIONS
    doubleprecision, dimension(:), allocatable :: shapeFactorX
    doubleprecision, dimension(:), allocatable :: shapeFactorY
    doubleprecision, dimension(:), allocatable :: shapeFactorZ
    doubleprecision, dimension(:), allocatable :: nParticlesX
    doubleprecision, dimension(:), allocatable :: nParticlesY
    doubleprecision, dimension(:), allocatable :: nParticlesZ

    integer :: nic
    integer :: newParticleGroupCount
    doubleprecision :: particleMass
    doubleprecision :: minOneParticleDensity

    ! FROM READLOCATIONS3
    type(ParticleGroupType) :: pGroup
    integer :: totalParticleCount,templateCount,templateCellCount,nc,nr,nl,row,column
    integer :: count,np,face,i,j,k,layerCount,rowCount,columnCount,       &    
               subCellCount,cell,offset,npcell
    doubleprecision :: dx,dy,dz,x,y,z,faceCoord,rowCoord,columnCoord,dr,dc 
    integer,dimension(:),allocatable :: templateSubDivisionTypes,         &
      templateCellNumbers,templateCellCounts, drape
    integer,dimension(:,:),allocatable :: subDiv
    integer,dimension(12) :: sdiv

    ! DIMENSIONALITY
    integer :: nDim
    integer, dimension(3) :: dimensionMask

    ! INJECTION BC's
    integer :: nInjectionConditions, nInjectionTimes
    integer :: injectionFormat
    integer :: injectionCellNumber
    doubleprecision, dimension(:), allocatable :: injectionTimes, injectionDensity
    character(len=200) :: injectionSeriesFile
    integer :: tempInjectionUnit = 667
    integer :: it, res
    doubleprecision :: totalMass, deltaTRelease

    !---------------------------------------------------------------------------------------------------------------

        
        ! Open dispersion unit 
        open( inUnit, file=inFile, status='old', access='sequential')


        ! Required for u3d
        allocate(cellsPerLayer(grid%LayerCount))
        do n = 1, grid%LayerCount
            cellsPerLayer(n) = grid%GetLayerCellCount(n)
        end do
    

        ! Write header to the listing file
        write(outUnit, *)
        write(outUnit, '(1x,a)') 'MODPATH RWPT dispersion data file'
        write(outUnit, '(1x,a)') '----------------------------'


        ! Read dispersivity variable, eventually file
        ! These methods follor OPEN/CLOSE, CONSTANT input format
        ! and variables are expected to be defined for each layer


        ! ALPHALONG
        if((grid%GridType .eq. 1) .or. (grid%GridType .eq. 3)) then
            call u3ddblmp(inUnit, outUnit, grid%LayerCount, grid%RowCount,      &
              grid%ColumnCount, grid%CellCount, this%AlphaLong, ANAME(2))                      
        else if((grid%GridType .eq. 2) .or. (grid%GridType .eq. 4)) then
            call u3ddblmpusg(inUnit, outUnit, grid%CellCount, grid%LayerCount,  &
              this%AlphaLong, aname(2), cellsPerLayer)
        else
            write(outUnit,*) 'Invalid grid type specified when reading ALPHALONG array data.'
            write(outUnit,*) 'Stopping.'
            call ustop(' ')          
        end if
        

        ! ALPHATRANS
        if((grid%GridType .eq. 1) .or. (grid%GridType .eq. 3)) then
            call u3ddblmp(inUnit, outUnit, grid%LayerCount, grid%RowCount,      &
              grid%ColumnCount, grid%CellCount, this%AlphaTrans, ANAME(2))                      
        else if((grid%GridType .eq. 2) .or. (grid%GridType .eq. 4)) then
            call u3ddblmpusg(inUnit, outUnit, grid%CellCount, grid%LayerCount,  &
              this%AlphaTrans, aname(2), cellsPerLayer)
        else
              write(outUnit,*) 'Invalid grid type specified when reading ALPHATRANS array data.'
              write(outUnit,*) 'Stopping.'
              call ustop(' ')          
        end if


        ! Dmol
        read( inUnit, * ) line
        icol = 1
        call urword(line,icol,istart,istop,3,n,r,0,0)
        this%DMol = r
        trackingOptions%Dmol = r

        ! Time Step kind 
        read(inUnit, '(a)') line
        icol = 1
        call urword(line,icol,istart,istop,1,n,r,0,0)
        line = line(istart:istop)
        ! Give a Courant
        if ( line .eq. 'CONSTANT_CU' ) then
            trackingOptions%timeStepKind = 1
            read( inUnit, * ) line
            icol = 1
            call urword(line,icol,istart,istop,3,n,r,0,0)
            trackingOptions%timeStepParameters(1) = r
        ! Give a Peclet
        else if ( line .eq. 'CONSTANT_PE' ) then
            trackingOptions%timeStepKind = 2
            read( inUnit, * ) line
            icol = 1
            call urword(line,icol,istart,istop,3,n,r,0,0)
            trackingOptions%timeStepParameters(2) = r
        ! Minimum between Courant and Peclet 
        else if ( line .eq. 'MIN_CU_PE' ) then
            trackingOptions%timeStepKind = 3
            read( inUnit, * ) line
            icol = 1
            call urword(line,icol,istart,istop,3,n,r,0,0)
            trackingOptions%timeStepParameters(1) = r
            read( inUnit, * ) line
            icol = 1
            call urword(line,icol,istart,istop,3,n,r,0,0)
            trackingOptions%timeStepParameters(2) = r
        else
            call ustop('RWPT: Invalid options for time step selection. Stop.')
        end if


        ! Advection Integration Kind
        read(inUnit, '(a)', iostat=iodispersion) line
        if ( iodispersion .lt. 0 ) then 
            ! end of file
            trackingOptions%advectionKind = 1
            write(outUnit,'(A)') 'RWPT: Advection integration not specified. Default to EXPONENTIAL.'
        else
            icol = 1
            call urword(line,icol,istart,istop,1,n,r,0,0)
            line = line(istart:istop)
            select case(line)
                case('EXPONENTIAL')
                    trackingOptions%advectionKind = 1
                    write(outUnit,'(A)') 'RWPT: Advection integration is EXPONENTIAL.'
                case('EULERIAN')
                    trackingOptions%advectionKind = 2
                    write(outUnit,'(A)') 'RWPT: Advection integration is EULERIAN.'
                case default
                    trackingOptions%advectionKind = 1
                    write(outUnit,'(A)') 'RWPT: Advection integration not specified. Default to EXPONENTIAL.'
            end select
        end if


        ! Domain dimensions option
        read(inUnit, '(a)', iostat=iodispersion) line
        if ( iodispersion .lt. 0 ) then 
            ! end of file
            trackingOptions%twoDimensions = .false.
            write(outUnit,'(A)') 'RWPT: Number of dimensions not specified. Defaults to 3D.'
        else
            icol = 1
            call urword(line,icol,istart,istop,1,n,r,0,0)
            line = line(istart:istop)
            select case(line)
                case('2D')
                    trackingOptions%twoDimensions = .true.
                    write(outUnit,'(A)') 'RWPT: Selected 2D domain solver.'
                case('3D')
                    trackingOptions%twoDimensions = .false.
                    write(outUnit,'(A)') 'RWPT: Selected 3D domain solver.'
                case default
                    trackingOptions%twoDimensions = .false.
                    write(outUnit,'(A)') 'RWPT: Invalid option for domain solver, defaults to 3D.'
            end select
        end if 


        ! ICBOUND
        if((grid%GridType .eq. 1) .or. (grid%GridType .eq. 3)) then
            call u3dintmp(inUnit, outUnit, grid%LayerCount, grid%RowCount,      &
              grid%ColumnCount, grid%CellCount, this%ICBound, ANAME(3))
        else if((grid%GridType .eq. 2) .or. (grid%GridType .eq. 4)) then
            call u3dintmpusg(inUnit, outUnit, grid%CellCount, grid%LayerCount,  &
              this%ICBound, aname(3), cellsPerLayer)
        else
              write(outUnit,*) 'Invalid grid type specified when reading ICBOUND array data.'
              write(outUnit,*) 'Stopping.'
              call ustop(' ')          
        end if
        ! END ICBOUND


        ! THIS COULD BE REUSED LATER FOR DIMENSIONALITY ADAPTIVITY
        ! Determine model dimensions
        dimensionMask = 0
        if ( grid%ColumnCount .gt. 1 ) dimensionMask(1) = 1 ! this dimension is active
        if ( grid%RowCount .gt. 1 )    dimensionMask(2) = 1 ! this dimension is active
        if ( grid%LayerCount .gt. 1 )  dimensionMask(3) = 1 ! this dimension is active
        nDim = sum(dimensionMask)


        ! INITIAL CONDITIONS
        ! MORE LIKE SPECIES AND EACH WITH AN INITIAL DISTRIBUTION
        ! THESE ARE TRANSFORMED INTO MASS PARTICLES 
        read(inUnit, *) nInitialConditions
        write(outUnit,'(/A,I5)') 'Number of initial conditions = ', nInitialConditions
        particleCount = 0
        slocUnit = 0
        seqNumber = 0


        ! ADD THEM TO EXISTING PARTICLE GROUPS
        ! LOOP IN MPATH IS ALL OVER PARTICLEGROUPS STORED IN SIMULATIONDATA
        if(nInitialConditions .gt. 0) then

            ! Would be like a realloc to extend 
            allocate(particleGroups(nInitialConditions))

            ! Loop over initial conditions
            do nic = 1, nInitialConditions

                particleGroups(nic)%Group = simulationData%ParticleGroupCount + nic

                ! Set release time for initial condition.
                ! In the meantime force 0d0 although it could 
                ! be extracted from a previuosly defined stage
                initialReleaseTime = 0d0
                call particleGroups(nic)%SetReleaseOption1(initialReleaseTime)

                ! Initial condition format 
                read(inUnit, *) initialConditionFormat ! 1: concentration, 2: particles (classic)

                read(inUnit, '(a)') particleGroups(nic)%Name

                select case ( initialConditionFormat )
                  case (1) 
                           
                      read(inUnit, *) massProcessingFormat ! 1:

                      select case ( massProcessingFormat )
                          case (1)
                              ! Read as mass concentration (ML^-3)

                              ! Given a value for the mass of particles, 
                              ! use flowModelData to compute cell volume 
                              ! and estimated number of particles per cell

                              ! Density/concentration arrays are expected 
                              ! to be consistent with flow model grid.
                              if(allocated(densityDistribution)) deallocate(densityDistribution)
                              allocate(densityDistribution(grid%CellCount))
                              if(allocated(rawNParticles)) deallocate(rawNParticles)
                              allocate(rawNParticles(grid%CellCount))
                              if(allocated(nParticles)) deallocate(nParticles)
                              allocate(nParticles(grid%CellCount))
                              if(allocated(cellVolumes)) deallocate(cellVolumes)
                              allocate(cellVolumes(grid%CellCount))
                              if(allocated(delZ)) deallocate(delZ)
                              allocate(delZ(grid%CellCount))
                              if(allocated(shapeFactorX)) deallocate(shapeFactorX)
                              allocate(shapeFactorX(grid%CellCount))
                              if(allocated(shapeFactorY)) deallocate(shapeFactorY)
                              allocate(shapeFactorY(grid%CellCount))
                              if(allocated(shapeFactorZ)) deallocate(shapeFactorZ)
                              allocate(shapeFactorZ(grid%CellCount))
                              if(allocated(nParticlesX)) deallocate(nParticlesX)
                              allocate(nParticlesX(grid%CellCount))
                              if(allocated(nParticlesY)) deallocate(nParticlesY)
                              allocate(nParticlesY(grid%CellCount))
                              if(allocated(nParticlesZ)) deallocate(nParticlesZ)
                              allocate(nParticlesZ(grid%CellCount))

                              ! Now convert to particles using mass
                              read(inUnit, *) particleMass

                              ! READ AS DENSITY/CONCENTRATION
                              if((grid%GridType .eq. 1) .or. (grid%GridType .eq. 3)) then
                                  call u3ddblmp(inUnit, outUnit, grid%LayerCount, grid%RowCount,      &
                                    grid%ColumnCount, grid%CellCount, densityDistribution, ANAME(4))
                              else if((grid%GridType .eq. 2) .or. (grid%GridType .eq. 4)) then
                                  call u3ddblmpusg(inUnit, outUnit, grid%CellCount, grid%LayerCount,  &
                                    densityDistribution, aname(4), cellsPerLayer)
                              else
                                    write(outUnit,*) 'Invalid grid type specified when reading IC array ', & 
                                        particleGroups(nic)%Name, ' name.'
                                    write(outUnit,*) 'Stopping.'
                                    call ustop(' ')          
                              end if
                              
   
                              ! Initialize
                              nParticles   = 0
                              cellVolumes  = 0d0
                              shapeFactorX = 0
                              shapeFactorY = 0
                              shapeFactorZ = 0
                              nParticlesX  = 0
                              nParticlesY  = 0
                              nParticlesZ  = 0


                              ! Do it simple (there's a volume involved)
                              ! Following is the particles density 
                              ! absolute value is required for the case that 
                              ! densityDsitribution contains negative values
                              rawNParticles = abs(densityDistribution)/particleMass


                              ! Simple delZ 
                              delZ      = grid%Top-grid%Bottom
                              ! LayerType if .eq. 1 convertible
                              where( this%Grid%CellType .eq. 1 )
                                  ! ONLY IF HEADS .lt. TOP 
                                  delZ = flowModelData%Heads-grid%Bottom
                              end where
                              where ( delZ .le. 0d0 )
                                  delZ = 0d0 
                              end where

                              ! Compute cell volumes
                              cellVolumes = 1
                              do n =1, 3
                                  if ( dimensionMask(n) .eq. 1 ) then 
                                      select case( n ) 
                                          case(1)
                                              cellVolumes = cellVolumes*grid%DelX
                                          case(2)
                                              cellVolumes = cellVolumes*grid%DelY
                                          case(3)
                                              cellVolumes = cellVolumes*delZ
                                      end select
                                  end if 
                              end do 


                              ! nParticles
                              ! notice rawNParticles = mass/volume/mass = 1/volume: particle density
                              nParticles  = rawNParticles*flowModelData%Porosity*cellVolumes

                              ! The minimum measurable absolute concentration/density 
                              minOneParticleDensity = minval(particleMass/flowModelData%Porosity/cellVolumes)


                          case default
                              continue

                      end select


                      ! Read the number of cells
                      templateCount = 1 
                      templateCellCount = grid%CellCount 

                      ! Allocate temporary arrays
                      allocate(subDiv(templateCount,12))
                      allocate(templateSubDivisionTypes(templateCount))
                      allocate(templateCellCounts(templateCount))
                      allocate(drape(templateCount))
                      allocate(templateCellNumbers(templateCellCount))

                      ! TEMP
                      drape = 0 ! SHOULD COME FROM SOMEWHERE
                      templateSubDivisionTypes(1) = 1
                      templateCellCounts(1) = grid%CellCount

                      ! Read data into temporary arrays and count the number of particles
                      rowCount    = grid%RowCount
                      columnCount = grid%ColumnCount
                      layerCount  = grid%LayerCount
                      do n = 1, templateCount
                          do m = 1, 12
                              subDiv(n,m) = 0
                          end do
                      end do
                     
                      np = 0
                      npcell = 0
                      offset = 0
                      do n = 1, templateCount
                          ! SHOULD DETERMINE SUBDIVISIONS GIVEN 
                          ! THE NUMBER OF PARTICLES. 
                          ! AS A DESIGN CONSIDERATION; THIS 
                          ! MASS PROCESSING METHOD CONSIDERS 
                          ! UNIFORM DISTRIBUTION OF PARTICLES
                          ! SOMETHING WITH A CELL SHAPE FACTOR
                          ! cellVolumes is already known
                          ! RAW ESTIMATE
                          
                          if ( dimensionMask(1) .eq. 1 ) then 
                              ! Shape factors
                              do m =1, grid%CellCount
                                  if ( cellVolumes(m) .le. 0d0 ) cycle
                                  shapeFactorX(m) = grid%DelX(m)/(cellVolumes(m)**(1d0/nDim))
                              end do 
                          end if 

                          if ( dimensionMask(2) .eq. 1 ) then 
                              ! Shape factors
                              do m =1, grid%CellCount
                                  if ( cellVolumes(m) .le. 0d0 ) cycle
                                  shapeFactorY(m) = grid%DelY(m)/(cellVolumes(m)**(1d0/nDim))
                              end do 
                          end if 

                          if ( dimensionMask(3) .eq. 1 ) then 
                              ! Shape factors
                              do m =1, grid%CellCount
                                  if ( cellVolumes(m) .le. 0d0 ) cycle
                                  shapeFactorZ(m) = delZ(m)/(cellVolumes(m)**(1d0/nDim))
                              end do 
                          end if 

                          where ( nParticles .lt. 1d0 )
                             nParticles = 0d0 
                          end where 
                          nParticlesX = shapeFactorX*(nParticles)**(1d0/nDim) 
                          nParticlesY = shapeFactorY*(nParticles)**(1d0/nDim)
                          nParticlesZ = shapeFactorZ*(nParticles)**(1d0/nDim)


                          ! THESE LINE READS THE CELL IDS THAT NEED TO BE PROCESSED
                          ! read(inUnit, *) (templateCellNumbers(offset + m), m = 1, templateCellCounts(n))
                          ! ALL ARE PROCESSED HERE
                          ! Loop through cells and count the number of particles
                          do cell = 1, templateCellCounts(n)
                              if(ibound(cell) .ne. 0) then
                                  ! Could be a user defined threshold
                                  if ( nParticles(cell) .lt. 1d0 ) cycle
                                  ! Cell subdivisions
                                  subDiv(n,1) = int( nParticlesX(cell) ) + 1
                                  subDiv(n,2) = int( nParticlesY(cell) ) + 1
                                  subDiv(n,3) = int( nParticlesZ(cell) ) + 1
                                  npcell      = subDiv(n,1)*subDiv(n,2)*subDiv(n,3)
                                  np = np + npcell
                              end if
                          end do
                          
                          ! Increment the offset
                          offset = offset + templateCellCounts(n)

                      end do


                      ! Calculate the total number of particles for all release time points.
                      totalParticleCount = 0
                      totalParticleCount = np*particleGroups(nic)%GetReleaseTimeCount()
                      if(allocated(particleGroups(nic)%Particles)) deallocate(particleGroups(nic)%Particles)
                      allocate(particleGroups(nic)%Particles(totalParticleCount))
                      particleGroups(nic)%TotalParticleCount = totalParticleCount
                     

                      ! UPDATE PARTICLES MASS
                      ! If for a given cell the value of density or concentration 
                      ! is negative, then assign the number of particles
                      ! and modify the sign of particles mass for that cell 
                      ! with a negative sign 
                      particleMass = sum( abs(densityDistribution)*cellVolumes*flowModelData%Porosity )/totalParticleCount 

                      ! Set the data for particles at the first release time point
                      m = 0
                      offset = 0
                      if(templateCount .gt. 0) then
                          do n = 1, templateCount
                              do cell = 1, templateCellCounts(n)
                                  if(ibound(cell) .ne. 0) then
                                      if ( nParticles(cell) .lt. 1d0 ) cycle
                                      sdiv(1) = int( nParticlesX(cell) ) + 1
                                      sdiv(2) = int( nParticlesY(cell) ) + 1
                                      sdiv(3) = int( nParticlesZ(cell) ) + 1
                                      if ( densityDistribution(cell) .gt. 0 ) then 
                                          particleMass = abs( particleMass ) 
                                      else 
                                          particleMass = -1*abs( particleMass )
                                      end if 
                                      call CreateMassParticlesAsInternalArray(& 
                                          particleGroups(nic), cell, m, sdiv(1), sdiv(2), sdiv(3), & 
                                          drape(n), particleMass, particleGroups(nic)%GetReleaseTime(1) )
                                  end if
                              end do
                              
                              ! Increment the offset
                              offset = offset + templateCellCounts(n)

                          end do
                      end if
                      

                      ! Deallocate temporary arrays
                      deallocate(subDiv)
                      deallocate(templateSubDivisionTypes)
                      deallocate(drape)
                      deallocate(templateCellCounts)
                      deallocate(templateCellNumbers)
                      !!  END  SIMPLIFIED FROM READLOATIONS3 !!!!!!!!!!!!!!!!!!!!!!!!!!


                  case (2) ! Read particles, the classical way + mass parameters
                   
                      ! releaseOption
                      read(inUnit, *) releaseOption
                      select case (releaseOption)
                          case (1)
                              read(inUnit, *) initialReleaseTime
                              call particleGroups(nic)%SetReleaseOption1(initialReleaseTime)
                          case (2)
                              read(inUnit, *) releaseTimeCount, initialReleaseTime, releaseInterval
                              call particleGroups(nic)%SetReleaseOption2(initialReleaseTime, &
                                releaseTimeCount, releaseInterval)
                          case (3)
                              read(inUnit, *) releaseTimeCount
                              if(allocated(releaseTimes)) deallocate(releaseTimes)
                              allocate(releaseTimes(releaseTimeCount))
                              read(inUnit, *) (releaseTimes(nn), nn = 1, releaseTimeCount)
                              call particleGroups(nic)%SetReleaseOption3(releaseTimeCount,   &
                                releaseTimes)
                          case default
                          ! write error message and stop
                      end select
                      
                      read(inUnit, '(a)') line
                      icol = 1
                      call urword(line,icol,istart,istop,1,n,r,0,0)
                      if(line(istart:istop) .eq. 'EXTERNAL') then
                          call urword(line,icol,istart,istop,0,n,r,0,0)
                          particleGroups(nic)%LocationFile = line(istart:istop)
                          slocUnit = 0
                      else if(line(istart:istop) .eq. 'INTERNAL') then
                          particleGroups(nic)%LocationFile = ''
                          slocUnit = inUnit
                      else
                          call ustop('Invalid starting locations file name. stop.')
                      end if
                      
                      ! LOAD LOCATION FILES 
                      call ReadAndPrepareLocations(slocUnit, outUnit, particleGroups(nic),   &
                        ibound, grid%CellCount, grid, seqNumber)

                      ! Read particles mass
                      read(inUnit, *) particleMass


                end select

                ! REPORT NUMBER OF PARTICLES
                write(outUnit, '(a,i4,a,i10,a)') 'Initial condition ', n, ' contains ',   &
                  particleGroups(nic)%TotalParticleCount, ' particles.'
                particleCount = particleCount + particleGroups(nic)%TotalParticleCount


            end do ! loop over particle groups

            write(outUnit, '(a,i10)') 'Total number of particles on initial conditions = ', particleCount
            write(outUnit, *)


            ! EXTEND SIMULATIONDATA TO INCLUDE THESE PARTICLE GROUPS
            ! OR RETURN THESE PARTICLEGROUPS AND REALLOCATE OUTSIDE
            newParticleGroupCount = simulationData%ParticleGroupCount + nInitialConditions
            allocate(newParticleGroups(newParticleGroupCount))
            do n = 1, simulationData%ParticleGroupCount
                newParticleGroups(n) = simulationData%ParticleGroups(n)
            end do 
            do n = 1, nInitialConditions
                newParticleGroups(n+simulationData%ParticleGroupCount) = particleGroups(n)
            end do 

            call move_alloc( newParticleGroups, simulationData%ParticleGroups )
            simulationData%ParticleGroupCount = newParticleGroupCount

            deallocate( particleGroups )
        

        end if ! if nInitialConditions .gt. 0



        ! INJECTION BOUNDARY CONDITIONS
        ! SHOULD BE RELATED TO ICBOUND ?
        read(inUnit, *) nInjectionConditions
        write(outUnit,'(/A,I5)') 'Number of mass injection cells = ', nInjectionConditions

        ! ADD THEM TO EXISTING PARTICLE GROUPS
        ! LOOP IN MPATH IS ALL OVER PARTICLEGROUPS STORED IN SIMULATIONDATA
        if(nInjectionConditions .gt. 0) then

            ! Allocate particlegroups  
            allocate(particleGroups(nInjectionConditions))

            ! Allocate these quantities to single size
            if(allocated(cellVolumes)) deallocate(cellVolumes)
            allocate(cellVolumes(1))
            if(allocated(delZ)) deallocate(delZ)
            allocate(delZ(1))
            if(allocated(shapeFactorX)) deallocate(shapeFactorX)
            allocate(shapeFactorX(1))
            if(allocated(shapeFactorY)) deallocate(shapeFactorY)
            allocate(shapeFactorY(1))
            if(allocated(shapeFactorZ)) deallocate(shapeFactorZ)
            allocate(shapeFactorZ(1))

            ! Loop over initial conditions
            do nic = 1, nInjectionConditions

                particleGroups(nic)%Group = simulationData%ParticleGroupCount + nic

                ! Injection condition format 
                ! 1: cell, constant concentration, start-end times
                ! 2: cell, injection time series
                read(inUnit, *) injectionFormat

                read(inUnit, '(a)') particleGroups(nic)%Name

                read(inUnit, *) injectionCellNumber ! SOMETHING TO VERIFY THAT IS VALID (IBOUND!=0)

                read(inUnit, *) particleMass

                read(inUnit, *) nInjectionTimes

                select case ( injectionFormat ) 
                    case (1)
                        continue 
                    case (2)
                        ! Density/concentration arrays are expected 
                        ! to be consistent with flow model grid.
                        if(allocated(injectionTimes)) deallocate(injectionTimes)
                        allocate(injectionTimes(nInjectionTimes))
                        if(allocated(injectionDensity)) deallocate(injectionDensity)
                        allocate(injectionDensity(nInjectionTimes))


                        ! Allocate arrays in time for nparticles
                        if(allocated(rawNParticles)) deallocate(rawNParticles)
                        allocate(rawNParticles(nInjectionTimes))
                        if(allocated(nParticles)) deallocate(nParticles)
                        allocate(nParticles(nInjectionTimes))
                        if(allocated(nParticlesX)) deallocate(nParticlesX)
                        allocate(nParticlesX(nInjectionTimes))
                        if(allocated(nParticlesY)) deallocate(nParticlesY)
                        allocate(nParticlesY(nInjectionTimes))
                        if(allocated(nParticlesZ)) deallocate(nParticlesZ)
                        allocate(nParticlesZ(nInjectionTimes))


                        ! Read file and load 
                        read(inUnit, '(a)') injectionSeriesFile 
                        open(tempInjectionUnit, file=injectionSeriesFile, &
                            access='sequential', form="formatted", iostat=res) 
                        do it = 1, nInjectionTimes
                            read(tempInjectionUnit,*) injectionTimes( it ), injectionDensity( it )
                        end do


                        ! Set release option
                        call particleGroups(nic)%SetReleaseOption3(nInjectionTimes, injectionTimes)


                        ! Initialize
                        nParticles   = 0
                        cellVolumes  = 0d0
                        shapeFactorX = 0
                        shapeFactorY = 0
                        shapeFactorZ = 0
                        nParticlesX  = 0
                        nParticlesY  = 0
                        nParticlesZ  = 0

                        ! Simple delZ 
                        delZ = grid%Top( injectionCellNumber ) - grid%Bottom(injectionCellNumber )
                        ! LayerType if .eq. 1 convertible
                        if( this%Grid%CellType( injectionCellNumber ) .eq. 1 ) then
                            if( flowModelData%Heads( injectionCellNumber ) .lt. & 
                                                grid%Top( injectionCellNumber ) ) then 
                                delZ = flowModelData%Heads( injectionCellNumber ) - grid%Bottom( injectionCellNumber )
                            end if
                        end if
                        if ( delZ(1) .le. 0d0 ) then
                            delZ(1) = 0d0 
                        end if

                        ! Compute cell volume
                        cellVolumes = 1 
                        do n =1,3
                            if ( dimensionMask(n) .eq. 1 ) then 
                                select case( n ) 
                                    case(1)
                                        cellVolumes = cellVolumes*grid%DelX( injectionCellNumber )
                                    case(2)
                                        cellVolumes = cellVolumes*grid%DelY( injectionCellNumber )
                                    case(3)
                                        cellVolumes = cellVolumes*delZ
                                end select
                            end if 
                        end do 


                        ! Allocate temporary arrays
                        ! TEMP
                        templateCount = 1 
                        templateCellCount = 1
                        allocate(subDiv(templateCount,12))
                        allocate(templateSubDivisionTypes(templateCount))
                        allocate(templateCellCounts(templateCount))
                        allocate(drape(templateCount))
                        allocate(templateCellNumbers(templateCellCount))
                        drape = 0 ! SHOULD COME FROM SOMEWHERE

                        np = 0
                        npcell = 0
                        offset = 0

                        ! THINGS THAT REMAIN EQUAL IN TIME
                        do n = 1, templateCount
                            do m = 1, 12
                                subDiv(n,m) = 0
                            end do
                        end do
                     
                        if ( dimensionMask(1) .eq. 1 ) then 
                            ! Shape factors
                            if ( cellVolumes(1) .le. 0d0 ) cycle
                            shapeFactorX(1) = grid%DelX(injectionCellNumber)/(cellVolumes(1)**(1d0/nDim))
                        end if 

                        if ( dimensionMask(2) .eq. 1 ) then 
                            ! Shape factors
                            if ( cellVolumes(1) .le. 0d0 ) cycle
                            shapeFactorY(1) = grid%DelY(injectionCellNumber)/(cellVolumes(1)**(1d0/nDim))
                        end if 

                        if ( dimensionMask(3) .eq. 1 ) then 
                            ! Shape factors
                            if ( cellVolumes(1) .le. 0d0 ) cycle
                            shapeFactorZ(1) = delZ(1)/(cellVolumes(1)**(1d0/nDim))
                        end if 

                        do it =1, nInjectionTimes

                            ! For each injectionDensity 
                            rawNParticles(it) = abs(injectionDensity(it))/particleMass

                            ! nParticles
                            ! notice rawNParticles = mass/volume/mass = 1/volume: particle density
                            nParticles(it)  = rawNParticles(it)*flowModelData%Porosity( injectionCellNumber )*cellVolumes(1)

                            if ( nParticles(it) .lt. 1d0 ) then
                               nParticles(it) = 0d0
                               cycle
                            end if 

                            nParticlesX(it) = shapeFactorX(1)*(nParticles(it))**(1d0/nDim) 
                            nParticlesY(it) = shapeFactorY(1)*(nParticles(it))**(1d0/nDim)
                            nParticlesZ(it) = shapeFactorZ(1)*(nParticles(it))**(1d0/nDim)

                            ! Cell subdivisions
                            subDiv(1,1) = int( nParticlesX(it) ) + 1
                            subDiv(1,2) = int( nParticlesY(it) ) + 1
                            subDiv(1,3) = int( nParticlesZ(it) ) + 1
                            npcell      = subDiv(1,1)*subDiv(1,2)*subDiv(1,3)
                            np          = np + npcell
                            
                        end do

                        ! AFTER LAST LOOP THE TOTAL NUMBER OF PARTICLES 
                        ! TO BE USED IN INJECTION TIMESERIES
                        ! Calculate the total number of particles for all release time points.
                        totalParticleCount = 0
                        totalParticleCount = np
                        if(allocated(particleGroups(nic)%Particles)) deallocate(particleGroups(nic)%Particles)
                        allocate(particleGroups(nic)%Particles(totalParticleCount))
                        particleGroups(nic)%TotalParticleCount = totalParticleCount

                        !print *, 'TOTAL PARTICLE COUNT INJECTION ', totalParticleCount

                        ! INTEGRATE MASS IN TIME
                        !totalMass = 0d0
                        !deltaTRelease = 0d0
                        !do n =2,nInjectionTimes
                        !    deltaTRelease = injectionTimes( n ) - injectionTimes( n-1 )
                        !end do

                        ! UPDATE PARTICLES MASS
                        ! THIS SHOULD BE DONE IN TIME
                        particleMass = sum( &
                            abs(injectionDensity)*cellVolumes(1)*flowModelData%Porosity( injectionCellNumber ) &
                        )/totalParticleCount 

                        m = 0
                        do it =1, nInjectionTimes
                            ! Set the data for particles at the first release time point
                            if ( nParticles(it) .lt. 1d0 ) cycle
                            sdiv(1) = int( nParticlesX(it) ) + 1
                            sdiv(2) = int( nParticlesY(it) ) + 1
                            sdiv(3) = int( nParticlesZ(it) ) + 1
                            if ( injectionDensity(it) .gt. 0 ) then 
                                particleMass = abs( particleMass ) 
                            else 
                                particleMass = -1*abs( particleMass )
                            end if 
                            call CreateMassParticlesAsInternalArray(     & 
                                particleGroups(nic), injectionCellNumber,&
                                m, sdiv(1), sdiv(2), sdiv(3),            &
                                drape(1), particleMass,                  &
                                injectionTimes( it )                     & 
                            )
                        end do 

                        ! Deallocate temporary arrays
                        deallocate(subDiv)
                        deallocate(templateSubDivisionTypes)
                        deallocate(drape)
                        deallocate(templateCellCounts)
                        deallocate(templateCellNumbers)

                end select 
    

            end do ! Loop over nInjectionConditions


            ! EXTEND SIMULATIONDATA TO INCLUDE THESE PARTICLE GROUPS
            ! OR RETURN THESE PARTICLEGROUPS AND REALLOCATE OUTSIDE
            newParticleGroupCount = simulationData%ParticleGroupCount + nInjectionConditions
            allocate(newParticleGroups(newParticleGroupCount))
            do n = 1, simulationData%ParticleGroupCount
                newParticleGroups(n) = simulationData%ParticleGroups(n)
            end do 
            do n = 1, nInjectionConditions
                newParticleGroups(n+simulationData%ParticleGroupCount) = particleGroups(n)
            end do 

            call move_alloc( newParticleGroups, simulationData%ParticleGroups )
            simulationData%ParticleGroupCount = newParticleGroupCount

            deallocate( particleGroups )


        end if ! if nInjectionConditions .gt. 0


    
        ! NONLINEARDISPERSIONDEV
        trackingOptions%dispersionModel = 2
        trackingOptions%mediumDistance  = 1
        trackingOptions%mediumDelta     = 5.5
        trackingOptions%betaLong        = 1
        trackingOptions%betaTrans       = 0.5






        ! Close dispersion data file
        close( inUnit )


    end subroutine pr_ReadData


    subroutine pr_Reset(this)
    !***************************************************************************************************************
    !
    !***************************************************************************************************************
    ! Specifications
    !---------------------------------------------------------------------------------------------------------------
    implicit none
    class(TransportModelDataType) :: this
    !---------------------------------------------------------------------------------------------------------------
       
        this%Grid => null()
        this%DMol = 0
        if(allocated(this%AlphaLong)) deallocate( this%AlphaLong )
        if(allocated(this%AlphaTrans)) deallocate( this%AlphaTrans )
        if(allocated(this%ICBoundTS)) deallocate( this%ICBoundTS )
        if(allocated(this%ICBound)) deallocate( this%ICBound )

    end subroutine pr_Reset



    subroutine pr_LoadTimeStep(this, stressPeriod, timeStep)
    !***************************************************************************************************************
    !
    !***************************************************************************************************************
    ! Specifications
    !---------------------------------------------------------------------------------------------------------------
    implicit none
    class(TransportModelDataType) :: this
    integer,intent(in) :: stressPeriod, timeStep
    integer :: firstRecord, lastRecord, n, m, firstNonBlank, lastNonBlank, &
      trimmedLength
    integer :: spaceAssigned, status,cellCount, iface, index,              &
      boundaryFlowsOffset, listItemBufferSize, cellNumber, layer
    !type(BudgetRecordHeaderType) :: header
    character(len=16) :: textLabel
    doubleprecision :: top 
    real :: HDryTol, HDryDiff
    !---------------------------------------------------------------------------------------------------------------
     

        !! Readers of time variable data 
        !call this%ClearTimeStepBudgetData()
        !call this%BudgetReader%GetRecordHeaderRange(stressPeriod, timeStep, firstRecord, lastRecord)
        !if(firstRecord .eq. 0) return
    
        cellCount = this%Grid%CellCount
        !listItemBufferSize = size(this%ListItemBuffer)
        
        ! Set steady state = true, then change it if the budget file contains storage
        !this%SteadyState = .true.
        
        ! MASS TRANSPORT; RWPT: could be recycled for changes in boundary conditions
        ! Load heads for this time step
        !call this%HeadReader%FillTimeStepHeadBuffer(stressPeriod, timeStep, &
        !  this%Heads, cellCount, spaceAssigned)
        
        ! Fill ICBoundTS array and set the SaturatedTop array for the Grid.
        ! The saturated top is set equal to the top for confined cells and water table cells 
        ! where the head is above the top or below the bottom.
        !HDryTol = abs(epsilon(HDryTol)*sngl(this%HDry))
        if(this%Grid%GridType .gt. 2) then
            do n = 1, cellCount
                this%ICBoundTS(n) = this%ICBound(n)
                ! DO SOMETHING WITH IBOUNDTS 
                ! DEPENDNING ON HEADS OR OTHERS
                !layer = this%Grid%GetLayer(n)
                !if(this%Grid%CellType(n) .eq. 1) then
                !    HDryDiff = sngl(this%Heads(n)) - sngl(this%HDry)
                !    if(abs(HDryDiff) .lt. HDryTol) then
                !        this%IBoundTS(n) = 0
                !        if(this%Heads(n) .lt. this%Grid%Bottom(n)) then
                !            this%IBoundTS(n) = 0
                !            this%Grid%SaturatedTop(n) = this%Grid%Bottom(n)
                !        end if
                !    end if
                !    if(this%IBoundTS(n) .ne. 0) then
                !        if((this%Heads(n) .le. this%Grid%Top(n)) .and. &
                !          (this%Heads(n) .ge. this%Grid%Bottom(n))) then
                !            this%Grid%SaturatedTop(n) = this%Heads(n)
                !        end if
                !    end if
                !end if
            end do
            
        else
            do n = 1, cellCount
                this%ICBoundTS(n) = this%ICBound(n)
                !this%Grid%SaturatedTop(n) = this%Grid%Top(n)
                !this%StorageFlows(n) = 0.0
                !this%IBoundTS(n) = this%IBound(n)
                !layer = this%Grid%GetLayer(n)
                !if(this%Grid%CellType(n) .eq. 1) then
                !    HDryDiff = sngl(this%Heads(n)) - sngl(this%HDry)
                !    if((abs(HDryDiff) .lt. HDryTol) .or. (this%Heads(n) .gt. 1.0d+6)) then
                !        this%IBoundTS(n) = 0
                !    end if
                !    if(this%IBoundTS(n) .ne. 0) then
                !        if((this%Heads(n) .le. this%Grid%Top(n)) .and. &
                !          (this%Heads(n) .ge. this%Grid%Bottom(n))) then
                !            this%Grid%SaturatedTop(n) = this%Heads(n)
                !        end if
                !    end if
                !end if
            end do
        end if
       

        !! MASS TRANSPORT, RWPT: COULD BE USEFUL FOR READING SOMETHING FROM EXTERNAL FILE 
        !! Loop through record headers
        !do n = firstRecord, lastRecord
        !     header = this%BudgetReader%GetRecordHeader(n)
        !     textLabel = header%TextLabel
        !     call TrimAll(textLabel, firstNonBlank, lastNonBlank, trimmedLength)
        !     
        !     select case(textLabel(firstNonBlank:lastNonBlank))
        !     case('CONSTANT HEAD', 'CHD')
        !          ! Read constant head flows into the sinkFlows and sourceFlows arrays.
        !          ! For a standard budget file, Method = 0. For a compact budget file,
        !          ! Method = 2.
        !          if(header%Method .eq. 0) then
        !              call this%BudgetReader%FillRecordDataBuffer(header,             &
        !                this%ArrayBufferDbl, cellCount, spaceAssigned, status)
        !              if(cellCount .eq. spaceAssigned) then
        !                  do m = 1, spaceAssigned
        !                      if(this%ArrayBufferDbl(m) .gt. 0.0d0) then
        !                          this%SourceFlows(m) = this%SourceFlows(m) +         &
        !                            this%ArrayBufferDbl(m)
        !                      else if(this%ArrayBufferDbl(m) .lt. 0.0d0) then
        !                          this%SinkFlows(m) = this%SinkFlows(m) +             &
        !                            this%ArrayBufferDbl(m)
        !                      end if
        !                  end do
        !              end if
        !          else if(header%Method .eq. 2) then
        !              call this%BudgetReader%FillRecordDataBuffer(header,             &
        !                this%ListItemBuffer, listItemBufferSize, spaceAssigned, status)
        !              if(spaceAssigned .gt. 0) then
        !                  do m = 1, spaceAssigned
        !                      cellNumber = this%ListItemBuffer(m)%CellNumber
        !                      if(this%ListItemBuffer(m)%BudgetValue .gt. 0.0d0) then
        !                          this%SourceFlows(cellNumber) =                      &
        !                            this%SourceFlows(cellNumber) + this%ListItemBuffer(m)%BudgetValue
        !                      else if(this%ListItemBuffer(m)%BudgetValue .lt. 0.0d0) then
        !                          this%SinkFlows(cellNumber) =                        &
        !                            this%SinkFlows(cellNumber) + this%ListItemBuffer(m)%BudgetValue
        !                      end if
        !                  end do
        !              end if
        !          else if((header%Method .eq. 5) .or. (header%Method .eq. 6)) then
        !              call this%BudgetReader%FillRecordDataBuffer(header,             &
        !                this%ListItemBuffer, listItemBufferSize, spaceAssigned,       &
        !                status)
        !              if(spaceAssigned .gt. 0) then
        !                  do m = 1, spaceAssigned
        !                      call this%CheckForDefaultIface(header%TextLabel, iface)
        !                      index = header%FindAuxiliaryNameIndex('IFACE')
        !                      if(index .gt. 0) then
        !                          iface = int(this%ListItemBuffer(m)%AuxiliaryValues(index))
        !                      end if
        !                      
        !                      cellNumber = this%ListItemBuffer(m)%CellNumber
        !                      if(iface .gt. 0) then
        !                          boundaryFlowsOffset = 6 * (cellNumber - 1)
        !                          this%BoundaryFlows(boundaryFlowsOffset + iface) =   &
        !                            this%BoundaryFlows(boundaryFlowsOffset + iface) + &
        !                            this%ListItemBuffer(m)%BudgetValue
        !                      else
        !                          if(this%ListItemBuffer(m)%BudgetValue .gt. 0.0d0) then
        !                              this%SourceFlows(cellNumber) =                  &
        !                                this%SourceFlows(cellNumber) +                &
        !                                this%ListItemBuffer(m)%BudgetValue
        !                          else if(this%ListItemBuffer(m)%BudgetValue .lt. 0.0d0) then
        !                              this%SinkFlows(cellNumber) =                    &
        !                                this%SinkFlows(cellNumber) +                  &
        !                                this%ListItemBuffer(m)%BudgetValue
        !                          end if
        !                      end if
        !                  end do
        !              end if
        !          end if
        !         
        !     case('STORAGE', 'STO-SS', 'STO-SY')
        !          ! Read storage for all cells into the StorageFlows array.
        !          ! Method should always be 0 or 1, but check anyway to be sure.
        !          if((header%Method .eq. 0) .or. (header%Method .eq. 1)) then
        !              if(header%ArrayItemCount .eq. cellCount) then
        !                  call this%BudgetReader%FillRecordDataBuffer(header,         &
        !                    this%ArrayBufferDbl, cellCount, spaceAssigned, status)
        !                  if(cellCount .eq. spaceAssigned) then
        !                      do m = 1, spaceAssigned
        !                          this%StorageFlows(m) = this%StorageFlows(m) + this%ArrayBufferDbl(m)
        !                          if(this%StorageFlows(m) .ne. 0.0) this%SteadyState = .false.
        !                      end do
        !                  end if
        !              end if
        !          end if
        !         
        !     case('FLOW JA FACE', 'FLOW-JA-FACE')
        !          ! Read connected face flows into the FlowsJA array for unstructured grids.
        !          if((header%Method .eq. 0) .or. (header%Method .eq. 1)) then
        !              ! Method should always be 0 or 1 for flow between grid cells. 
        !              if(header%ArrayItemCount .eq. this%BudgetReader%GetFlowArraySize()) then
        !                  call this%BudgetReader%FillRecordDataBuffer(header,         &
        !                    this%FlowsJA, header%ArrayItemCount, spaceAssigned,       &
        !                    status)
        !              end if
        !          else if(header%Method .eq. 6) then
        !              ! Method code 6 indicates flow to or from cells in the current model grid
        !              ! and another connected model grid in a multi-model MODFLOW-6 simulation. 
        !              ! Treat flows to or from connected model grids as distributed source/sink flows 
        !              ! for the current grid.
        !              call this%BudgetReader%FillRecordDataBuffer(header,             &
        !                this%ListItemBuffer, listItemBufferSize, spaceAssigned,       &
        !                status)
        !              if(spaceAssigned .gt. 0) then
        !                  do m = 1, spaceAssigned
        !                      cellNumber = this%ListItemBuffer(m)%CellNumber
        !                      if(this%ListItemBuffer(m)%BudgetValue .gt. 0.0d0) then
        !                          this%SourceFlows(cellNumber) =                  &
        !                              this%SourceFlows(cellNumber) +                &
        !                              this%ListItemBuffer(m)%BudgetValue
        !                      else if(this%ListItemBuffer(m)%BudgetValue .lt. 0.0d0) then
        !                          this%SinkFlows(cellNumber) =                    &
        !                              this%SinkFlows(cellNumber) +                  &
        !                              this%ListItemBuffer(m)%BudgetValue
        !                      end if
        !                  end do
        !              end if
        !          end if
        !         
        !     case('FLOW RIGHT FACE')
        !          ! Read flows across the right face for structured grids.
        !          ! Method should always be 0 or 1, but check anyway to be sure.
        !          if((header%Method .eq. 0) .or. (header%Method .eq. 1)) then
        !              if(header%ArrayItemCount .eq. this%BudgetReader%GetFlowArraySize()) then
        !                  call this%BudgetReader%FillRecordDataBuffer(header,         &
        !                    this%FlowsRightFace, header%ArrayItemCount, spaceAssigned,&
        !                    status)
        !              end if
        !          end if
        !         
        !     case('FLOW FRONT FACE')
        !          ! Read flows across the front face for structured grids.
        !          ! Method should always be 0 or 1, but check anyway to be sure.
        !          if((header%Method .eq. 0) .or. (header%Method .eq. 1)) then
        !              if(header%ArrayItemCount .eq. this%BudgetReader%GetFlowArraySize()) then
        !                  call this%BudgetReader%FillRecordDataBuffer(header,         &
        !                    this%FlowsFrontFace, header%ArrayItemCount, spaceAssigned,&
        !                    status)
        !              end if
        !          end if
        !         
        !     case('FLOW LOWER FACE')
        !          ! Read flows across the lower face for structured grids.
        !          ! Method should always be 0 or 1, but check anyway to be sure.
        !          if((header%Method .eq. 0) .or. (header%Method .eq. 1)) then
        !              if(header%ArrayItemCount .eq. this%BudgetReader%GetFlowArraySize()) then
        !                  call this%BudgetReader%FillRecordDataBuffer(header,         &
        !                    this%FlowsLowerFace, header%ArrayItemCount, spaceAssigned,&
        !                    status)
        !              end if
        !          end if
        !     
        !      case default
        !          ! Now handle any other records in the budget file.
        !           if((header%Method .eq. 0) .or. (header%Method .eq. 1)) then
        !              if(header%ArrayItemCount .eq. cellCount) then
        !                  call this%BudgetReader%FillRecordDataBuffer(header,         &
        !                    this%ArrayBufferDbl, cellCount, spaceAssigned, status)
        !                  if(cellCount .eq. spaceAssigned) then
        !                      call this%CheckForDefaultIface(header%TextLabel, iface)
        !                      if(iface .gt. 0) then
        !                          do m = 1, spaceAssigned
        !                              boundaryFlowsOffset = 6 * (m - 1)
        !                              this%BoundaryFlows(boundaryFlowsOffset + iface) =   &
        !                                this%BoundaryFlows(boundaryFlowsOffset + iface) + &
        !                                this%ArrayBufferDbl(m)
        !                          end do
        !                      else
        !                          do m = 1, spaceAssigned
        !                              if(this%ArrayBufferDbl(m) .gt. 0.0d0) then
        !                                  this%SourceFlows(m) = this%SourceFlows(m) +     &
        !                                    this%ArrayBufferDbl(m)
        !                              else if(this%ArrayBufferDbl(m) .lt. 0.0d0) then
        !                                  this%SinkFlows(m) = this%SinkFlows(m) +         &
        !                                    this%ArrayBufferDbl(m)
        !                              end if
        !                          end do
        !                      end if
        !                  end if
        !              end if
        !           else if(header%Method .eq. 3) then
        !              call this%BudgetReader%FillRecordDataBuffer(header,             &
        !                this%ArrayBufferDbl, this%ArrayBufferInt,                     &
        !                header%ArrayItemCount, spaceAssigned, status)
        !              if(header%ArrayItemCount .eq. spaceAssigned) then
        !                  call this%CheckForDefaultIface(header%TextLabel, iface)
        !                  if(iface .gt. 0) then
        !                      do m = 1, spaceAssigned
        !                          cellNumber = this%ArrayBufferInt(m)
        !                          boundaryFlowsOffset = 6 * (cellNumber - 1)
        !                          this%BoundaryFlows(boundaryFlowsOffset + iface) =   &
        !                            this%BoundaryFlows(boundaryFlowsOffset + iface) + &
        !                            this%ArrayBufferDbl(m)
        !                      end do
        !                  else            
        !                      do m = 1, spaceAssigned
        !                          cellNumber = this%ArrayBufferInt(m)
        !                          if(this%ArrayBufferDbl(m) .gt. 0.0d0) then
        !                              this%SourceFlows(cellNumber) =                  &
        !                                this%SourceFlows(cellNumber) +                &
        !                                this%ArrayBufferDbl(m)
        !                          else if(this%ArrayBufferDbl(m) .lt. 0.0d0) then
        !                              this%SinkFlows(cellNumber) =                    &
        !                                this%SinkFlows(cellNumber) +                  &
        !                                this%ArrayBufferDbl(m)
        !                          end if
        !                      end do
        !                  end if
        !              end if
        !           else if(header%Method .eq. 4) then
        !              call this%BudgetReader%FillRecordDataBuffer(header,             &
        !                this%ArrayBufferDbl, header%ArrayItemCount, spaceAssigned,    &
        !                status)
        !              if(header%ArrayItemCount .eq. spaceAssigned) then
        !                  call this%CheckForDefaultIface(header%TextLabel, iface)
        !                  if(iface .gt. 0) then
        !                      do m = 1, spaceAssigned
        !                          boundaryFlowsOffset = 6 * (m - 1)
        !                          this%BoundaryFlows(boundaryFlowsOffset + iface) =   &
        !                            this%BoundaryFlows(boundaryFlowsOffset + iface) + &
        !                            this%ArrayBufferDbl(m)
        !                      end do
        !                  else            
        !                      do m = 1, spaceAssigned
        !                          if(this%ArrayBufferDbl(m) .gt. 0.0d0) then
        !                              this%SourceFlows(m) = this%SourceFlows(m) +     &
        !                                this%ArrayBufferDbl(m)
        !                          else if(this%ArrayBufferDbl(m) .lt. 0.0d0) then
        !                              this%SinkFlows(m) = this%SinkFlows(m) +         &
        !                                this%ArrayBufferDbl(m)
        !                          end if
        !                      end do
        !                  end if
        !              end if
        !          else if(header%Method .eq. 2) then
        !              call this%BudgetReader%FillRecordDataBuffer(header,             &
        !                this%ListItemBuffer, listItemBufferSize, spaceAssigned,       &
        !                status)
        !              if(spaceAssigned .gt. 0) then
        !                  call this%CheckForDefaultIface(header%TextLabel, iface)
        !                  if(iface .gt. 0) then
        !                      do m = 1, spaceAssigned
        !                          cellNumber = this%ListItemBuffer(m)%CellNumber
        !                          boundaryFlowsOffset = 6 * (cellNumber - 1)
        !                          this%BoundaryFlows(boundaryFlowsOffset + iface) =   &
        !                            this%BoundaryFlows(boundaryFlowsOffset + iface) + &
        !                            this%ListItemBuffer(m)%BudgetValue
        !                      end do
        !                  else            
        !                      do m = 1, spaceAssigned
        !                          cellNumber = this%ListItemBuffer(m)%CellNumber
        !                          if(this%ListItemBuffer(m)%BudgetValue .gt. 0.0d0) then
        !                              this%SourceFlows(cellNumber) =                  &
        !                                this%SourceFlows(cellNumber) +                &
        !                                this%ListItemBuffer(m)%BudgetValue
        !                          else if(this%ListItemBuffer(m)%BudgetValue .lt. 0.0d0) then
        !                              this%SinkFlows(cellNumber) =                    &
        !                                this%SinkFlows(cellNumber) +                  &
        !                                this%ListItemBuffer(m)%BudgetValue
        !                          end if
        !                      end do
        !                  end if
        !              end if
        !          else if((header%Method .eq. 5) .or. (header%Method .eq. 6)) then
        !              call this%BudgetReader%FillRecordDataBuffer(header,             &
        !                this%ListItemBuffer, listItemBufferSize, spaceAssigned,       &
        !                status)
        !              if(spaceAssigned .gt. 0) then
        !                  do m = 1, spaceAssigned
        !                      call this%CheckForDefaultIface(header%TextLabel, iface)
        !                      index = header%FindAuxiliaryNameIndex('IFACE')
        !                      if(index .gt. 0) then
        !                          iface = int(this%ListItemBuffer(m)%AuxiliaryValues(index))
        !                      end if
        !                      
        !                      cellNumber = this%ListItemBuffer(m)%CellNumber
        !                      if(iface .gt. 0) then
        !                          boundaryFlowsOffset = 6 * (cellNumber - 1)
        !                          this%BoundaryFlows(boundaryFlowsOffset + iface) =   &
        !                            this%BoundaryFlows(boundaryFlowsOffset + iface) + &
        !                            this%ListItemBuffer(m)%BudgetValue
        !                      else
        !                          if(this%ListItemBuffer(m)%BudgetValue .gt. 0.0d0) then
        !                              this%SourceFlows(cellNumber) =                  &
        !                                this%SourceFlows(cellNumber) +                &
        !                                this%ListItemBuffer(m)%BudgetValue
        !                          else if(this%ListItemBuffer(m)%BudgetValue .lt. 0.0d0) then
        !                              this%SinkFlows(cellNumber) =                    &
        !                                this%SinkFlows(cellNumber) +                  &
        !                                this%ListItemBuffer(m)%BudgetValue
        !                          end if
        !                      end if
        !                  end do
        !              end if
        !          end if
        !     
        !     end select
        !     
        !end do
    
        this%CurrentStressPeriod = stressPeriod
        this%CurrentTimeStep = timeStep
    
    end subroutine pr_LoadTimeStep


end module TransportModelDataModule
