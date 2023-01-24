module TransportModelDataModule
  use ModflowRectangularGridModule,only : ModflowRectangularGridType
  use FlowModelDataModule,only : FlowModelDataType
  use SoluteModule,only : SoluteType
  use DispersionDataModule,only : DispersionDataType
  use ModpathSimulationDataModule, only: ModpathSimulationDataType
  !---------------------------------------------------------------------------------
  use ParticleTrackingOptionsModule,only : ParticleTrackingOptionsType  ! DEPRECATED
  use ParticleGroupModule,only : ParticleGroupType                      ! DEPRECATED
  use StartingLocationReaderModule,only : ReadAndPrepareLocations,    & ! DEPRECATED
                                          ReadAndPrepareLocationsMass,& ! DEPRECATED
                                   pr_CreateParticlesAsInternalArray, & ! DEPRECATED
                                   CreateMassParticlesAsInternalArray   ! DEPRECATED
  implicit none
  !---------------------------------------------------------------------------------------------------------------

  ! Set default access status to private
  private

    type,public :: TransportModelDataType

      logical :: Initialized = .false.

      ! Local dispersion parameters
      doubleprecision,dimension(:),pointer :: AlphaLong => null()
      doubleprecision,dimension(:),pointer :: AlphaTran => null()
      doubleprecision,dimension(:),pointer :: AlphaL => null()
      doubleprecision,dimension(:),pointer :: AlphaTH => null()
      doubleprecision,dimension(:),pointer :: AlphaTV => null()

      doubleprecision,dimension(:),allocatable :: MediumDistance ! Nonlinear model

      ! Local ICBound
      integer,allocatable,dimension(:)         :: ICBound
      integer,allocatable,dimension(:)         :: ICBoundTS


      doubleprecision :: DMol

      ! Simulation data
      class( ModpathSimulationDataType ), pointer :: simulationData


      ! Solutes
      type( SoluteType ), allocatable, dimension(:) :: Solutes
      type( SoluteType ) :: BaseSolute
      integer            :: nSolutes

      ! Dispersion models 
      type( DispersionDataType ), allocatable, dimension(:) :: DispersionData 
      type( DispersionDataType ) :: BaseDispersionData
      integer                    :: nDispersion

      ! grid
      class(ModflowRectangularGridType),pointer :: Grid => null()

      ! Private variables
      integer,private :: CurrentStressPeriod = 0
      integer,private :: CurrentTimeStep = 0

    contains

      procedure :: Initialize=>pr_Initialize
      procedure :: Reset=>pr_Reset
      procedure :: ReadData=>pr_ReadData ! DEPRECATION WARNING
      procedure :: ReadSPCData=>pr_ReadSPCData
      procedure :: ReadDSPData=>pr_ReadDSPData
      procedure :: LoadTimeStep=>pr_LoadTimeStep
      procedure :: LoadSoluteDispersion=>pr_LoadSoluteDispersion
      procedure :: SetSoluteDispersion=>pr_SetSoluteDispersion

    end type


contains


  subroutine pr_Initialize( this, grid, simulationData )
  !---------------------------------------------------------------------------------------------------------------
  ! Specifications
  !---------------------------------------------------------------------------------------------------------------
  implicit none
  class(TransportModelDataType) :: this
  class(ModflowRectangularGridType),intent(inout),pointer :: grid
  class(ModpathSimulationDataType),intent(in),pointer :: simulationData

  integer :: cellCount, gridType
  !---------------------------------------------------------------------------------------------------------------

      ! Call Reset to make sure that all arrays are initially unallocated
      this%Initialized = .false.
      call this%Reset()
      
      ! Return if the grid cell count equals 0
      cellCount = grid%CellCount
      if(cellCount .le. 0) return
      
      ! Check budget reader and grid data for compatibility and allocate appropriate cell-by-cell flow arrays
      gridType = grid%GridType
      this%Grid => grid

      ! simulationData
      this%simulationData => simulationData

      ! IBounds
      if(allocated(this%ICBoundTS)) deallocate(this%ICBoundTS)
      if(allocated(this%ICBound)) deallocate(this%ICBound)
      allocate(this%ICBoundTS(grid%CellCount))
      allocate(this%ICBound(grid%CellCount))


      this%Initialized = .true.


  end subroutine pr_Initialize


  ! REQUIRES REVIEW
  subroutine pr_Reset(this)
  !---------------------------------------------------------------------------------------------------------------
  ! Specifications
  !---------------------------------------------------------------------------------------------------------------
  implicit none
  class(TransportModelDataType) :: this
  !---------------------------------------------------------------------------------------------------------------
     
      this%Grid => null()
      this%DMol = 0
      !if(allocated(this%AlphaLong)) deallocate( this%AlphaLong )
      !if(allocated(this%AlphaTran)) deallocate( this%AlphaTran )
      this%AlphaLong  => null()
      this%AlphaTran => null()
      if(allocated(this%ICBoundTS)) deallocate( this%ICBoundTS )
      if(allocated(this%ICBound)) deallocate( this%ICBound )

  end subroutine pr_Reset


  ! Read specific SPC data
  subroutine pr_ReadSPCData( this, spcFile, spcUnit, outUnit )
    use UTL8MODULE,only : urword,ustop
    !--------------------------------------------------------------
    ! Specifications
    !--------------------------------------------------------------
    implicit none
    class(TransportModelDataType), target :: this
    character(len=200), intent(in)           :: spcFile
    integer, intent(in)                      :: spcUnit
    integer, intent(in)                      :: outUnit
    ! local
    class(ModpathSimulationDataType), pointer :: simulationData
    integer :: isThisFileOpen = -1
    integer :: icol,istart,istop,n
    doubleprecision    :: r
    character(len=200) :: line
    integer :: nSolutes, ns, ncount, npg
    !--------------------------------------------------------------

    write(outUnit, *)
    write(outUnit, '(1x,a)') 'MODPATH-RW SPC file data'
    write(outUnit, '(1x,a)') '------------------------'

    ! Local pointer to class pointer 
    simulationData => this%simulationData

    ! Verify if unit is open 
    inquire( file=spcFile, number=isThisFileOpen )
    if ( isThisFileOpen .lt. 0 ) then 
      ! No spc file
      write(outUnit,'(A)') 'SPC package was not specified in name file.'

      ! If no SPC package was specified, then verify 
      ! if ParticlesMassOption .eq. 2. In said case, 
      ! soluteIds were specified in particle groups.
      ! Follow this specification, all with the same dispersion

      if (simulationData%ParticlesMassOption .eq. 2) then
        ! Create different solutes for id purposes
        ! However, they are all displaced with the 
        ! same transport properties
        ! Extracted from particles groups defined until 
        ! this very moment. 
        continue
      else
        ! If not solute id, then assumes single solute and
        ! all pgroups are of the same kind 
        this%nSolutes = 1
        if ( allocated(this%Solutes) ) deallocate( this%Solutes )
        allocate( this%Solutes(this%nSolutes) )
        !this%Solutes = this%BaseSolute
      end if

      return

    end if

    ! If the SPC package was specified, 
    ! then read it. 
    ! Interpret solutesoption to determined 
    ! whether the simulation is multidispersion or not

    ! Define dispersion according to solutes option 
    select case( simulationData%SolutesOption )
      ! 0: Single solute dispersion 
      case (0)
        write(outUnit,'(A)') 'Simulation is single dispersion for all solutes, will not read SPC specs.'
        ! A single virtual solute storing the same 
        ! transport parameters for all pgroups
        this%BaseSolute%id = 0
        this%BaseSolute%stringid = 'SPC0' 

        ! All pgroups to base solute 
        if ( allocated(this%BaseSolute%pGroups) ) deallocate(this%BaseSolute%pGroups)
        allocate(this%BaseSolute%pGroups(simulationData%ParticleGroupCount))
        do n=1,simulationData%ParticleGroupCount
          this%BaseSolute%pGroups(n) = n 
        end do 
        this%BaseSolute%nParticleGroups = simulationData%ParticleGroupCount

      ! 1: Multiple solute, multidispersion
      case (1)
        write(outUnit,'(A)') 'Simulation is multidispersion, will read species/solutes specs.'

        ! How many ?
        read(spcUnit, '(a)') line
        icol = 1
        call urword(line, icol, istart, istop, 2, n, r, 0, 0)
        this%nSolutes = n 
       
        ! Process SPC's
        write(outUnit,'(A,I5)') 'Given number of species/solutes = ', this%nSolutes
        if(this%nSolutes .le. 0) then
          ! No spc's
          write(outUnit,'(A)') 'Number of given species/solutes is .le. 0.'
          ! Shall initialize BaseSolute?
          return
        end if

        ! Allocate solutes
        if ( allocated(this%Solutes) ) deallocate( this%Solutes )
        allocate( this%Solutes(this%nSolutes) )

        ! Initialize these guys
        do ns =1,this%nSolutes
          write(outUnit,'(A,I4)') 'Processing species/solute ', ns
          ! Of course needs some attributes 
          ! of the solute

          ! Read the integer id 
          read(spcUnit, '(a)') line
          icol = 1
          call urword(line, icol, istart, istop, 2, n, r, 0, 0)
          this%Solutes(ns)%id = n 
          
          ! Read the string id
          read(spcUnit, '(a)') line
          icol = 1
          call urword(line, icol, istart, istop, 0, n, r, 0, 0)
          this%Solutes(ns)%stringid = line(istart:istop)

          ! Assign pgroups related to the solute
          if ( simulationData%ParticlesMassOption .ne. 2 ) then
            ! Read pgroups related to the solute from the pkg
            write(outUnit,'(A)') 'Read pgroups related to this species/solutes from the pkg.'

            ! How many pgroups
            read(spcUnit, '(a)') line
            icol = 1
            call urword(line, icol, istart, istop, 2, n, r, 0, 0)
            this%Solutes(ns)%nParticleGroups = n 
           
            ! Some health check
            if ( allocated( this%Solutes(ns)%pGroups ) ) deallocate( this%Solutes(ns)%pGroups )
            allocate(this%Solutes(ns)%pGroups(this%Solutes(ns)%nParticleGroups))
            
            ! Read related particle groups
            do npg =1,this%Solutes(ns)%nParticleGroups
              ! Read the pgroups
              read(spcUnit, '(a)') line
              icol = 1
              call urword(line, icol, istart, istop, 2, n, r, 0, 0)
              this%Solutes(ns)%pGroups(npg) = n

              ! It needs to assign the soluteId back to the 
              ! corresponding pgroup for simulations 
              ! where the solute is not specified in the pgroup
              simulationData%ParticleGroups(&
                  this%Solutes(ns)%pGroups(npg) )%Solute = ns
            end do

          else if ( simulationData%ParticlesMassOption .eq. 2 ) then
            ! Read pgroups related to the solute ids at pgroups
            write(outUnit,'(A)') 'Interpret pgroups related to this species/solutes from pgroups list.'

            ! Read the soluteId's from simulationData%ParticleGroup list
            ! and count how many for this solute
            ncount = 0
            do npg =1,simulationData%ParticleGroupCount
              if ( simulationData%ParticleGroups( npg )%Solute .eq. ns ) then
                ncount = ncount + 1
              end if 
            end do

            ! If no pgroups, let it continue ?
            if ( ncount .eq. 0 ) then 
              write(outUnit,*) 'Warning: no particle groups associated to solute id ', ns
            else
              ! Initialize pgroup ids for the solute
              this%Solutes(ns)%nParticleGroups = ncount
          
              if ( allocated( this%Solutes(ns)%pGroups ) ) deallocate( this%Solutes(ns)%pGroups )
              allocate( this%Solutes(ns)%pGroups( &
                 this%Solutes(ns)%nParticleGroups ) )

              ncount = 0 
              do npg =1,simulationData%ParticleGroupCount
                if ( simulationData%ParticleGroups( npg )%Solute .eq. ns ) then
                  ncount = ncount + 1
                  this%Solutes(ns)%pGroups(ncount) = npg
                end if 
              end do
            end if 

          end if 
         
          ! Read the dispersion id
          read(spcUnit, '(a)') line
          icol = 1
          call urword(line, icol, istart, istop, 2, n, r, 0, 0)
          this%Solutes(ns)%dispersionId = n

          ! Read the dispersion string id
          read(spcUnit, '(a)') line
          icol = 1
          call urword(line, icol, istart, istop, 0, n, r, 0, 0)
          this%Solutes(ns)%dispersionStringid = line(istart:istop)

          write(outUnit,'(A,I4,A)')&
          'Specie/solute ',ns, ' is related to dispersion:',trim(adjustl(this%Solutes(ns)%dispersionStringid))
        end do

      case default
        ! Some error handling
        call ustop('Invalid specified solutes option. Stop. ') 
    end select
    


    ! Close spc data file
    close( spcUnit )


  end subroutine pr_ReadSPCData


  ! Read specific DSP data
  subroutine pr_ReadDSPData( this, dspFile, dspUnit, outUnit )
    use UTL8MODULE,only : urword,ustop,u3ddblmpusg, u3ddblmp
    !--------------------------------------------------------------
    ! Specifications
    !--------------------------------------------------------------
    implicit none
    ! input 
    class(TransportModelDataType), target :: this
    character(len=200), intent(in)           :: dspFile
    integer, intent(in)                      :: dspUnit
    integer, intent(in)                      :: outUnit
    ! local
    class(ModpathSimulationDataType), pointer :: simulationData
    class(ModflowRectangularGridType),pointer :: grid
    type(DispersionDataType),pointer          :: disp
    integer :: isThisFileOpen = -1
    integer :: icol,istart,istop,n
    doubleprecision    :: r
    character(len=200) :: line
    integer :: nDispersion, ndis
    integer, dimension(:), allocatable :: cellsPerLayer
    character(len=24),dimension(4) :: anamelin
    data anamelin(1) /'          DMEFF'/
    data anamelin(2) /'         ALPHAL'/
    data anamelin(3) /'        ALPHATH'/
    data anamelin(4) /'        ALPHATV'/
    character(len=24),dimension(6) :: anamenlin
    data anamenlin(1) /'          DMEFF'/
    data anamenlin(2) /'          BETAL'/
    data anamenlin(3) /'         BETATH'/
    data anamenlin(4) /'         BETATV'/
    data anamenlin(5) /'          DELTA'/
    data anamenlin(6) /'         DGRAIN'/
    !--------------------------------------------------------------

    write(outUnit, *)
    write(outUnit, '(1x,a)') 'MODPATH-RW DSP file data'
    write(outUnit, '(1x,a)') '------------------------'

    ! Verify if unit is open 
    inquire( file=dspFile, number=isThisFileOpen )
    if ( isThisFileOpen .lt. 0 ) then 
      ! No spc file
      write(outUnit,'(A)') 'DSP package was not specified in name file and is required for RW.'
      call ustop('DSP package was not specified in name file and is required for RW. Stop.')
    end if

    ! Local pointers to class pointers
    simulationData => this%simulationData
    grid => this%Grid

    !! Maybe it helps in reducing the amount of parameters 
    !! to be read
    !! RW dimensionality vars
    !dimensionMask => this%TrackingOptions%dimensionMask
    !nDim => this%TrackingOptions%nDim

    ! Process DSP's
    read(dspUnit, *) nDispersion
    write(outUnit,'(A,I5)') 'Given number of dispersion parameters/models = ', nDispersion

    if(nDispersion .le. 0) then
      ! No dsp
      write(outUnit,'(A)') 'Number of given dispersion parameters/models is .le. 0.'
      call ustop('Number of given dispersion parameters/models is .le. 0. Stop.')
    end if

    ! In case the simulation is not multidispersion, 
    ! read only the first dispersion model at DSP pkg
    if ( simulationData%SolutesOption .eq. 0 ) then
      ! 0: Single solute dispersion 
      write(outUnit,'(A)') 'Simulation is single dispersion, will read only first dispersion data.'
      nDispersion = 1
    end if   
    this%nDispersion = nDispersion


    ! Allocate dispersion data array
    if ( allocated(this%DispersionData) ) deallocate(this%DispersionData)
    allocate( this%DispersionData(this%nDispersion) )

    ! CellsPerLayer, required for u3d reader
    allocate(cellsPerLayer(grid%LayerCount))
    do n = 1, grid%LayerCount
      cellsPerLayer(n) = grid%GetLayerCellCount(n)
    end do

    ! Loop over dispersion data
    do ndis=1,this%nDispersion
      ! Report which DSP will be processed
      write(outUnit,'(A,I5)') 'Processing dispersion data: ', ndis

      ! Assing local pointer 
      disp => this%DispersionData(ndis)

      ! Read the integer id 
      read(dspUnit, '(a)') line
      icol = 1
      call urword(line, icol, istart, istop, 2, n, r, 0, 0)
      disp%id = n 
      
      ! Read the string id
      read(dspUnit, '(a)') line
      icol = 1
      call urword(line, icol, istart, istop, 0, n, r, 0, 0)
      disp%stringid = line(istart:istop)

      ! Read the dispersion modelKind
      read(dspUnit, '(a)') line
      icol = 1
      call urword(line, icol, istart, istop, 2, n, r, 0, 0)
      disp%modelKind = n 

      ! Allocate
      call disp%InitializeByModelKind( grid%CellCount )

      ! Read dispersion data
      select case ( disp%modelKind )
        ! Linear
        case(1)
          write(outUnit,'(A)') 'Dispersion model is linear, will read dmeff and dispersivities'

          ! Read DMEFF
          if((grid%GridType .eq. 1) .or. (grid%GridType .eq. 3)) then
            call u3ddblmp(dspUnit, outUnit, grid%LayerCount, grid%RowCount,      &
              grid%ColumnCount, grid%CellCount, disp%DMEff, anamelin(1))                      
          else if((grid%GridType .eq. 2) .or. (grid%GridType .eq. 4)) then
            call u3ddblmpusg(dspUnit, outUnit, grid%CellCount, grid%LayerCount,  &
              disp%DMEff, anamelin(1), cellsPerLayer)
          else
            write(outUnit,*) 'Invalid grid type specified when reading DMEFF array data.'
            write(outUnit,*) 'Stopping.'
            call ustop(' ')          
          end if

          ! Read ALPHAL
          if((grid%GridType .eq. 1) .or. (grid%GridType .eq. 3)) then
            call u3ddblmp(dspUnit, outUnit, grid%LayerCount, grid%RowCount,      &
              grid%ColumnCount, grid%CellCount, disp%AlphaL, anamelin(2))                      
          else if((grid%GridType .eq. 2) .or. (grid%GridType .eq. 4)) then
            call u3ddblmpusg(dspUnit, outUnit, grid%CellCount, grid%LayerCount,  &
              disp%AlphaL, anamelin(2), cellsPerLayer)
          else
            write(outUnit,*) 'Invalid grid type specified when reading ALPHAL array data.'
            write(outUnit,*) 'Stopping.'
            call ustop(' ')          
          end if

          ! Read ALPHATH
          if((grid%GridType .eq. 1) .or. (grid%GridType .eq. 3)) then
            call u3ddblmp(dspUnit, outUnit, grid%LayerCount, grid%RowCount,      &
              grid%ColumnCount, grid%CellCount, disp%AlphaTH, anamelin(3))                      
          else if((grid%GridType .eq. 2) .or. (grid%GridType .eq. 4)) then
            call u3ddblmpusg(dspUnit, outUnit, grid%CellCount, grid%LayerCount,  &
              disp%AlphaTH, anamelin(3), cellsPerLayer)
          else
            write(outUnit,*) 'Invalid grid type specified when reading ALPHATH array data.'
            write(outUnit,*) 'Stopping.'
            call ustop(' ')          
          end if

          ! Read ALPHATV
          if((grid%GridType .eq. 1) .or. (grid%GridType .eq. 3)) then
            call u3ddblmp(dspUnit, outUnit, grid%LayerCount, grid%RowCount,      &
              grid%ColumnCount, grid%CellCount, disp%AlphaTV, anamelin(4))                      
          else if((grid%GridType .eq. 2) .or. (grid%GridType .eq. 4)) then
            call u3ddblmpusg(dspUnit, outUnit, grid%CellCount, grid%LayerCount,  &
              disp%AlphaTV, anamelin(4), cellsPerLayer)
          else
            write(outUnit,*) 'Invalid grid type specified when reading ALPHATV array data.'
            write(outUnit,*) 'Stopping.'
            call ustop(' ')          
          end if

        ! Nonlinear
        case(2)
          write(outUnit,'(A)') 'Dispersion model is nonlinear'

          ! Aqueous molecular diffusion 
          read( dspUnit, * ) line
          icol = 1
          call urword(line,icol,istart,istop,3,n,r,0,0)
          disp%dmaqueous = r

          ! Read DMEFF
          if((grid%GridType .eq. 1) .or. (grid%GridType .eq. 3)) then
            call u3ddblmp(dspUnit, outUnit, grid%LayerCount, grid%RowCount,      &
              grid%ColumnCount, grid%CellCount, disp%DMEff, anamenlin(1))                      
          else if((grid%GridType .eq. 2) .or. (grid%GridType .eq. 4)) then
            call u3ddblmpusg(dspUnit, outUnit, grid%CellCount, grid%LayerCount,  &
              disp%DMEff, anamenlin(1), cellsPerLayer)
          else
            write(outUnit,*) 'Invalid grid type specified when reading DMEFF array data.'
            write(outUnit,*) 'Stopping.'
            call ustop(' ')          
          end if

          ! Read BETAL
          if((grid%GridType .eq. 1) .or. (grid%GridType .eq. 3)) then
            call u3ddblmp(dspUnit, outUnit, grid%LayerCount, grid%RowCount,      &
              grid%ColumnCount, grid%CellCount, disp%BetaL, anamenlin(2))
          else if((grid%GridType .eq. 2) .or. (grid%GridType .eq. 4)) then
            call u3ddblmpusg(dspUnit, outUnit, grid%CellCount, grid%LayerCount,  &
              disp%BetaL, anamenlin(2), cellsPerLayer)
          else
            write(outUnit,*) 'Invalid grid type specified when reading BETAL array data.'
            write(outUnit,*) 'Stopping.'
            call ustop(' ')          
          end if

          ! Read BETATH
          if((grid%GridType .eq. 1) .or. (grid%GridType .eq. 3)) then
            call u3ddblmp(dspUnit, outUnit, grid%LayerCount, grid%RowCount,      &
              grid%ColumnCount, grid%CellCount, disp%BetaTH, anamenlin(3))
          else if((grid%GridType .eq. 2) .or. (grid%GridType .eq. 4)) then
            call u3ddblmpusg(dspUnit, outUnit, grid%CellCount, grid%LayerCount,  &
              disp%BetaTH, anamenlin(3), cellsPerLayer)
          else
            write(outUnit,*) 'Invalid grid type specified when reading BETATH array data.'
            write(outUnit,*) 'Stopping.'
            call ustop(' ')          
          end if

          ! Read BETATV
          if((grid%GridType .eq. 1) .or. (grid%GridType .eq. 3)) then
            call u3ddblmp(dspUnit, outUnit, grid%LayerCount, grid%RowCount,      &
              grid%ColumnCount, grid%CellCount, disp%BetaTV, anamenlin(4))
          else if((grid%GridType .eq. 2) .or. (grid%GridType .eq. 4)) then
            call u3ddblmpusg(dspUnit, outUnit, grid%CellCount, grid%LayerCount,  &
              disp%BetaTV, anamenlin(4), cellsPerLayer)
          else
            write(outUnit,*) 'Invalid grid type specified when reading BETATV array data.'
            write(outUnit,*) 'Stopping.'
            call ustop(' ')          
          end if

          ! Read DELTA
          if((grid%GridType .eq. 1) .or. (grid%GridType .eq. 3)) then
            call u3ddblmp(dspUnit, outUnit, grid%LayerCount, grid%RowCount,      &
              grid%ColumnCount, grid%CellCount, disp%Delta, anamenlin(5)) 
          else if((grid%GridType .eq. 2) .or. (grid%GridType .eq. 4)) then
            call u3ddblmpusg(dspUnit, outUnit, grid%CellCount, grid%LayerCount,  &
              disp%Delta, anamenlin(5), cellsPerLayer)
          else
            write(outUnit,*) 'Invalid grid type specified when reading DELTA array data.'
            write(outUnit,*) 'Stopping.'
            call ustop(' ')          
          end if

          ! Read DGRAIN
          if((grid%GridType .eq. 1) .or. (grid%GridType .eq. 3)) then
            call u3ddblmp(dspUnit, outUnit, grid%LayerCount, grid%RowCount,      &
              grid%ColumnCount, grid%CellCount, disp%DGrain, anamenlin(6)) 
          else if((grid%GridType .eq. 2) .or. (grid%GridType .eq. 4)) then
            call u3ddblmpusg(dspUnit, outUnit, grid%CellCount, grid%LayerCount,  &
              disp%DGrain, anamenlin(6), cellsPerLayer)
          else
            write(outUnit,*) 'Invalid grid type specified when reading DGRAIN array data.'
            write(outUnit,*) 'Stopping.'
            call ustop(' ')          
          end if
        case default
          write(outUnit,*) 'Invalid dispersion model. Accepts 1 or 2, given ', disp%modelKind
          write(outUnit,*) 'Stopping.'
          call ustop('Invalid dispersion model. Stop')          
      end select

    end do ! end do this%nDispersion


    ! Close dsp data file
    close( dspUnit )


  end subroutine pr_ReadDSPData




    subroutine pr_LoadTimeStep(this, stressPeriod, timeStep)
    !---------------------------------------------------------------------------------------------------------------
    ! Specifications:
    !---------------------------------------------------------------------------------------------------------------
    !   - Huge simplification from flowModelData%LoadTimeStep
    !   - It could be employed for loading time variable data for dispersion
    !   - In the meantime, only update time references and ICBOUND 
    !---------------------------------------------------------------------------------------------------------------
    implicit none
    class(TransportModelDataType) :: this
    integer,intent(in) :: stressPeriod, timeStep
    integer :: firstRecord, lastRecord, n, m, firstNonBlank, lastNonBlank, &
      trimmedLength
    integer :: spaceAssigned, status,cellCount, iface, index,              &
      boundaryFlowsOffset, listItemBufferSize, cellNumber, layer
    character(len=16) :: textLabel
    doubleprecision :: top 
    real :: HDryTol, HDryDiff
    !---------------------------------------------------------------------------------------------------------------

      cellCount = this%Grid%CellCount
      if(this%Grid%GridType .gt. 2) then
          do n = 1, cellCount
              this%ICBoundTS(n) = this%ICBound(n)
          end do
      else
          do n = 1, cellCount
              this%ICBoundTS(n) = this%ICBound(n)
          end do
      end if

      this%CurrentStressPeriod = stressPeriod
      this%CurrentTimeStep = timeStep
   

    end subroutine pr_LoadTimeStep

    
    subroutine pr_LoadSoluteDispersion( this, inUnit, outUnit, solute, &
                                   grid, cellsPerLayer, trackingOptions) 
    !-----------------------------------------------------------------------------------
    ! Specifications
    !-----------------------------------------------------------------------------------
    use UTL8MODULE,only : urword, ustop, u1dint, u1drel, u1ddbl, u8rdcom, &
                          u3dintmp, u3dintmpusg, u3ddblmp, u3ddblmpusg
    implicit none
    class(TransportModelDataType) :: this
    class(SoluteType), intent(inout) :: solute
    class(ModflowRectangularGridType),intent(in) :: grid
    type(ParticleTrackingOptionsType),intent(inout) :: trackingOptions
    integer,intent(in) :: inUnit, outUnit 
    integer, dimension(:), intent(in) :: cellsPerLayer
    character(len=300) :: line
    character(len=24),dimension(3) :: aname
    data aname(1) /'          LONG DISP'/
    data aname(2) /'          TRAN DISP'/
    data aname(3) /'          TRAN DISP'/
    doubleprecision :: r
    integer :: icol, istart, istop, n , ns
    !-----------------------------------------------------------------------------------
    
      ! Initialize dispersion properties 
      call solute%Initialize( grid%CellCount ) 

      ! Aqueous diffusion 
      read( inUnit, * ) line
      icol = 1
      call urword(line,icol,istart,istop,3,n,r,0,0)
      solute%dAqueous = r

      ! Dispersivities
      select case( solute%dispersionModel )

        case( 1 ) ! Linear
          ! Read dispersivities
          ! These methods follow OPEN/CLOSE, CONSTANT input format
          ! and variables are expected to be defined for each layer
          ! ALPHALONG
          if((grid%GridType .eq. 1) .or. (grid%GridType .eq. 3)) then
            call u3ddblmp(inUnit, outUnit, grid%LayerCount, grid%RowCount,      &
              grid%ColumnCount, grid%CellCount, solute%AlphaLong, aname(1))                      
          else if((grid%GridType .eq. 2) .or. (grid%GridType .eq. 4)) then
            call u3ddblmpusg(inUnit, outUnit, grid%CellCount, grid%LayerCount,  &
              solute%AlphaLong, aname(1), cellsPerLayer)
          else
            write(outUnit,*) 'Invalid grid type specified when reading ALPHALONG array data.'
            write(outUnit,*) 'Stopping.'
            call ustop(' ')          
          end if
          
          ! ALPHATRANS
          if((grid%GridType .eq. 1) .or. (grid%GridType .eq. 3)) then
            call u3ddblmp(inUnit, outUnit, grid%LayerCount, grid%RowCount,      &
              grid%ColumnCount, grid%CellCount, solute%AlphaTran, aname(2))                      
          else if((grid%GridType .eq. 2) .or. (grid%GridType .eq. 4)) then
            call u3ddblmpusg(inUnit, outUnit, grid%CellCount, grid%LayerCount,  &
              solute%AlphaTran, aname(2), cellsPerLayer)
          else
            write(outUnit,*) 'Invalid grid type specified when reading ALPHATRANS array data.'
            write(outUnit,*) 'Stopping.'
            call ustop(' ')          
          end if


        case( 2 ) ! Nonlinear
          ! NONLINEAR

          ! Nonlinear model
          if(allocated(this%MediumDistance)) deallocate(this%MediumDistance)
          allocate(this%MediumDistance(grid%CellCount))

          ! betaLong
          read( inUnit, * ) line
          icol = 1
          call urword(line,icol,istart,istop,3,n,r,0,0)
          trackingOptions%betaLong  = r
          solute%betaLong  = r

          ! betaTrans
          read( inUnit, * ) line
          icol = 1
          call urword(line,icol,istart,istop,3,n,r,0,0)
          trackingOptions%betaTrans = r
          solute%betaTrans  = r

          ! mediumDelta
          read( inUnit, * ) line
          icol = 1
          call urword(line,icol,istart,istop,3,n,r,0,0)
          trackingOptions%mediumDelta = r

          ! MEDIUMDISTANCE
          if((grid%GridType .eq. 1) .or. (grid%GridType .eq. 3)) then
              call u3ddblmp(inUnit, outUnit, grid%LayerCount, grid%RowCount,      &
                grid%ColumnCount, grid%CellCount, this%MediumDistance, aname(5)) 
          else if((grid%GridType .eq. 2) .or. (grid%GridType .eq. 4)) then
              call u3ddblmpusg(inUnit, outUnit, grid%CellCount, grid%LayerCount,  &
                this%MediumDistance, aname(5), cellsPerLayer)
          else
              write(outUnit,*) 'Invalid grid type specified when reading MEDIUMDISTANCE array data.'
              write(outUnit,*) 'Stopping.'
              call ustop(' ')          
          end if

          ! TEMPORAL 
          trackingOptions%mediumDistance = this%MediumDistance(1)

          ! END NONLINEAR

        case default
          write(outUnit,*) 'Invalid dispersion model. Accepts 1 or 2, given ', solute%dispersionModel
          write(outUnit,*) 'Stopping.'
          call ustop('Invalid dispersion model. Stop')          
      end select


    end subroutine pr_LoadSoluteDispersion 


    subroutine pr_SetSoluteDispersion( this, soluteId )
    !-----------------------------------------------------------------------------------
    ! Specifications:
    !   - Assigns transport model data dispersivities
    !     with values from the corresponding soluteId
    !   - Runs only if SolutesOption .eq. 1, Multidispersion
    !-----------------------------------------------------------------------------------
    implicit none
    class(TransportModelDataType), target :: this
    integer, intent(in) :: soluteId
    !-----------------------------------------------------------------------------------

      ! Some sanity check

      this%AlphaLong => this%Solutes(soluteId)%AlphaLong
      this%AlphaTran => this%Solutes(soluteId)%AlphaTran

      ! Note this, needs clarification of the actual values 
      ! of diffusion to be used in displacements
      this%DMol = this%Solutes(soluteId)%dAqueous 


    end subroutine pr_SetSoluteDispersion


    ! DEPRECATION WARNING !
    subroutine pr_ReadData(this, inUnit, inFile, outUnit, simulationData, flowModelData, ibound, grid, trackingOptions )
    !***************************************************************************************************************
    !
    !***************************************************************************************************************
    ! Specifications
    !   - Reading of dispersion related input parameters
    !   - This function reads what comes immediately after particle groups (pgroups) in classical modpath
    !   - Logic:
    !     - Initial conditions (increase pgroups) 
    !     - Prescribed flux boundaries (increase pgroups)
    !     - Solutes and dispersion parameters, dispersivities
    !     - Time step configuration
    !     - Advection
    !     - ICBound
    !     - Dimension mask 
    !---------------------------------------------------------------------------------------------------------------
    use utl7module,only : urdcom, upcase
    use UTL8MODULE,only : urword, ustop, u1dint, u1drel, u1ddbl, u8rdcom, &
      u3dintmp, u3dintmpusg, u3ddblmp, u3ddblmpusg, ugetnode
    implicit none
    class(TransportModelDataType), target:: this
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
    character(len=24),dimension(5) :: aname
    data aname(1) /'          BOUNDARY ARRAY'/
    data aname(2) /'            DISPERSIVITY'/
    data aname(3) /'            ICBOUND'/
    data aname(4) /'            IC'/
    data aname(5) /'          MEDIUMDISTANCE'/
    integer :: tempAlphaUnit = 666
    character(len=200) :: tempAlphaFile

    ! Initial conditions
    integer :: nInitialConditions, particleCount, seqNumber, slocUnit
    integer :: nValidInitialConditions
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
    integer :: totalParticleCount,templateCount,templateCellCount,nc,nr,nl,row,column,idmax
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
    integer :: nInjectionConditions, nInjectionTimes, nValidInjectionConditions
    integer :: injectionFormat
    integer :: injectionCellNumber
    doubleprecision, dimension(:), allocatable :: injectionTimes, injectionDensity
    character(len=200) :: injectionSeriesFile
    integer :: tempInjectionUnit = 667
    integer :: it, res, npg, ns, ncount, readStyle
    doubleprecision :: totalMass, deltaTRelease
    integer :: cellNumber, ioInUnit

    ! Solutes 
    integer :: soluteId
    !----------------------------------------------------------------------------------------------------

        ! Open dispersion unit 
        open( inUnit, file=inFile, status='old', access='sequential')

        ! Required for u3d
        allocate(cellsPerLayer(grid%LayerCount))
        do n = 1, grid%LayerCount
            cellsPerLayer(n) = grid%GetLayerCellCount(n)
        end do

        ! Cell size
        rowCount    = grid%RowCount
        columnCount = grid%ColumnCount
        layerCount  = grid%LayerCount


        ! Determine model dimensions (READIT!)
        dimensionMask = 0
        if ( grid%ColumnCount .gt. 1 ) dimensionMask(1) = 1 ! this dimension is active
        if ( grid%RowCount .gt. 1 )    dimensionMask(2) = 1 ! this dimension is active
        if ( grid%LayerCount .gt. 1 )  dimensionMask(3) = 1 ! this dimension is active
        nDim = sum(dimensionMask)


        ! Write header to the listing file
        write(outUnit, *)
        write(outUnit, '(1x,a)') '------------------------------------'
        write(outUnit, '(1x,a)') ' MODPATH-RW configuration file      '
        write(outUnit, '(1x,a)') '------------------------------------'

        ! Process initial conditions 
        read(inUnit, *) nInitialConditions
        write(outUnit,'(/A,I5)') ' Number of initial conditions = ', nInitialConditions
        nValidInitialConditions = 0

        ! What is this ?
        particleCount = 0
        slocUnit      = 0
        seqNumber     = 0

        ! Extend particle groups
        if(nInitialConditions .gt. 0) then

          ! Carrier of particle groups
          allocate(particleGroups(nInitialConditions))

          ! Loop over initial conditions
          do nic = 1, nInitialConditions

            ! Increase pgroup counter
            particleGroups(nic)%Group = simulationData%ParticleGroupCount + nic

            ! Set release time for initial condition.
            ! It is an initial condition, then
            ! assumes release at referencetime
            initialReleaseTime = simulationData%ReferenceTime
            call particleGroups(nic)%SetReleaseOption1(initialReleaseTime)

            ! Read id 
            read(inUnit, '(a)') particleGroups(nic)%Name

            ! Initial condition format
            ! 1: concentration
            read(inUnit, *) initialConditionFormat

            select case ( initialConditionFormat )

            ! Read initial condition as concentration times 
            ! porosity  (ML^-3)
            case (1) 

              ! Given a value for the mass of particles, 
              ! use flowModelData to compute cellvolume
              ! and a shape factor from which the number 
              ! of particles per cell is estimated

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

              ! Read particles mass
              read(inUnit, *) particleMass
        
              if ( ( simulationData%ParticlesMassOption .eq. 2 ) .or. & 
                   ( simulationData%SolutesOption .eq. 1 ) ) then 
                ! Read solute id
                ! Some validation 
                read(inUnit, *) soluteId
              end if 

              ! Read as density/concentration
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
           
              ! Validity of initial density distribution 
              if( all(abs(densityDistribution).eq.0d0) ) then
                write(outUnit,*) 'Warning: initial condition ',&
                    particleGroups(nic)%Name,' has a densityDistribution only with zeros'
                write(outUnit,*) 'It will not create a particle group. Continue to the next.'

                ! Process the next one
                cycle

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

              ! Simple delZ 
              delZ = grid%Top-grid%Bottom
              ! LayerType if .eq. 1 convertible
              where( this%Grid%CellType .eq. 1 )
                  ! ONLY IF HEADS .lt. TOP 
                  delZ = flowModelData%Heads-grid%Bottom
              end where
              where ( delZ .le. 0d0 )
                  delZ = 0d0 
              end where

              ! Compute cell volumes
              cellVolumes = 1d0
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

              ! Particles density 
              ! absolute value is required for the case that 
              ! densityDsitribution contains negative values
              rawNParticles = abs(densityDistribution)/particleMass

              ! nParticles
              ! notice rawNParticles = mass/volume/mass = 1/volume: particle density
              nParticles  = rawNParticles*cellVolumes

              ! The minimum measurable absolute concentration/density 
              minOneParticleDensity = minval(particleMass/cellVolumes)

              ! Allocate subdivision arrays
              templateCount     = 1 
              templateCellCount = grid%CellCount 
              allocate(subDiv(templateCount,12))
              allocate(templateSubDivisionTypes(templateCount))
              allocate(templateCellCounts(templateCount))
              allocate(drape(templateCount))
              allocate(templateCellNumbers(templateCellCount))

              subDiv(:,:) = 0
              drape       = 0 ! Should come from somewhere
              templateSubDivisionTypes(1) = 1
              templateCellCounts(1)       = grid%CellCount

              np     = 0
              npcell = 0
              offset = 0
              ! There is only one templateCount
              do n = 1, templateCount
                ! Determine subdivisions
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
                nParticlesX = shapeFactorX*( (nParticles)**(1d0/nDim) ) 
                nParticlesY = shapeFactorY*( (nParticles)**(1d0/nDim) )
                nParticlesZ = shapeFactorZ*( (nParticles)**(1d0/nDim) )

                ! Loop through cells and count the number of particles
                do cell = 1, templateCellCounts(n)
                  if(ibound(cell) .ne. 0) then
                    ! Could be a user defined threshold ?
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
              particleGroups(nic)%TotalParticleCount = totalParticleCount
              if ( totalParticleCount .eq. 0 ) then 
                write(outUnit,*) ' Warning: initial condition ',&
                    particleGroups(nic)%Name,' has zero particles, it will skip this group.'

                ! Process the next one
                cycle

              end if 

              if(allocated(particleGroups(nic)%Particles)) deallocate(particleGroups(nic)%Particles)
              allocate(particleGroups(nic)%Particles(totalParticleCount))

              ! Update particles mass
              ! If for a given cell the value of density or concentration 
              ! is negative, then assign the number of particles
              ! and modify the sign of particles mass for that cell 
              ! with a negative sign 
              particleMass = sum( abs(densityDistribution)*cellVolumes )/totalParticleCount 
              write(outUnit,'(/A,es18.9e3)') ' Effective particleMass for initial condition = ', particleMass

              ! Assign to the particle group 
              particleGroups(nic)%Mass = particleMass

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
                      ! For the weird requirement where density
                      ! might be negative...
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
            
              ! Assign layer value to each particle
              idmax = 0
              seqNumber = 0
              do m = 1, totalParticleCount
                  seqNumber = seqNumber + 1
                  if(particleGroups(nic)%Particles(m)%Id .gt. idmax) idmax = particleGroups(nic)%Particles(m)%Id
                  particleGroups(nic)%Particles(m)%Group = particleGroups(nic)%Group
                  particleGroups(nic)%Particles(m)%SequenceNumber = seqNumber
                  particleGroups(nic)%Particles(m)%InitialLayer =                                   &
                    grid%GetLayer(particleGroups(nic)%Particles(m)%InitialCellNumber)
                  particleGroups(nic)%Particles(m)%Layer =                                          &
                    grid%GetLayer(particleGroups(nic)%Particles(m)%CellNumber)
              end do
  

              ! Deallocate temporary arrays
              deallocate(subDiv)
              deallocate(templateSubDivisionTypes)
              deallocate(drape)
              deallocate(templateCellCounts)
              deallocate(templateCellNumbers)

            case default
              write(outUnit,*) ' Invalid initial condition kind ', initialConditionFormat 
              write(outUnit,*) ' Stop.'
              call ustop(' Error: Invalid initial condition kind ')
            end select

            ! Increment valid couter
            nValidInitialConditions = nValidInitialConditions + 1 

            if ( simulationData%ParticlesMassOption .eq. 2 ) then 
              ! Assign the solute id 
              particleGroups(nic)%Solute = soluteId
            end if 

            ! Report number of particles
            write(outUnit, '(a,i4,a,i10,a)') ' Initial condition ', nic, ' contains ',   &
              particleGroups(nic)%TotalParticleCount, ' particles.'
            particleCount = particleCount + particleGroups(nic)%TotalParticleCount

          end do ! loop over initial conditions 
          write(outUnit, '(a,i10)') ' Total number of particles on initial conditions = ', particleCount
          write(outUnit, *)


          ! Extend simulationdata to include these particle groups
          if ( nValidInitialConditions .gt. 0 ) then 
            newParticleGroupCount = simulationData%ParticleGroupCount + nValidInitialConditions
            allocate(newParticleGroups(newParticleGroupCount))
            ! If some particle groups existed previously
            if( simulationData%ParticleGroupCount .gt. 0 ) then 
              do n = 1, simulationData%ParticleGroupCount
                newParticleGroups(n) = simulationData%ParticleGroups(n)
              end do
            end if 
            ncount = 0
            do n = 1, nInitialConditions
              if ( particleGroups(n)%TotalParticleCount .eq. 0 ) cycle
              ncount = ncount + 1 
              newParticleGroups(ncount+simulationData%ParticleGroupCount) = particleGroups(n)
            end do 
            if( simulationData%ParticleGroupCount .gt. 0 ) then 
              call move_alloc( newParticleGroups, simulationData%ParticleGroups )
              simulationData%ParticleGroupCount = newParticleGroupCount
              simulationData%TotalParticleCount = simulationData%TotalParticleCount + particleCount
            else
              simulationData%ParticleGroupCount = newParticleGroupCount
              simulationData%TotalParticleCount = particleCount
              allocate(simulationData%ParticleGroups(simulationData%ParticleGroupCount))
              call move_alloc( newParticleGroups, simulationData%ParticleGroups )
            end if
          end if

          ! Will process the next at some point
          deallocate( particleGroups )
       
        end if ! if nInitialConditions .gt. 0


        !-- FLUX BOUNDARIES NEED VERIFICATION REVIEW, REMOVE ? --!
            
        ! Prescribed flux boundary conditions
        read(inUnit, *) nInjectionConditions
        write(outUnit,'(/A,I5)') ' Number of mass injection (flux) conditions = ', nInjectionConditions
        nValidInjectionConditions = 0


        ! Extends particle groups 
        if(nInjectionConditions .gt. 0) then

          ! Allocate particlegroups  
          allocate(particleGroups(nInjectionConditions))

          ! Allocate these quantities to single size
          if(allocated(cellVolumes)) deallocate(cellVolumes)
          allocate(cellVolumes(1))
          if(allocated(delZ)) deallocate(delZ)
          allocate(delZ(1))
          if(allocated(shapeFactorX)) deallocate(shapeFactorX)
      
          if(allocated(shapeFactorY)) deallocate(shapeFactorY)
          allocate(shapeFactorY(1))
          if(allocated(shapeFactorZ)) deallocate(shapeFactorZ)
          allocate(shapeFactorZ(1))

          ! Loop over injection conditions
          do nic = 1, nInjectionConditions

            ! Increase pgroup counter
            particleGroups(nic)%Group = simulationData%ParticleGroupCount + nic

            ! Injection name/id
            read(inUnit, '(a)') particleGroups(nic)%Name

            ! Injection condition format 
            ! 1: cell, constant concentration, release times
            ! 2: cell, injection time series
            read(inUnit, *) injectionFormat

            select case ( injectionFormat ) 
              ! Cell id and some information 
              ! related to particles release
              case (1)

                ! Are these ids as (lay,row,col) or (cellid) ?
                read(inUnit, '(a)', iostat=ioInUnit) line
                icol = 1
                call urword(line, icol, istart, istop, 2, n, r, 0, 0)
                readStyle = n

                ! Read the injection cell
                if( readStyle .eq. 1) then
                  ! Read as layer, row, column
                  read(inUnit, *) layer, row, column
                  call ugetnode(layerCount, rowCount, columnCount, layer, row, column,cellNumber)
                  injectionCellNumber = cellNumber
                else if ( readStyle .eq. 2 ) then 
                  read(inUnit,*) injectionCellNumber
                else
                    call ustop(' Invalid read style for injection cell. Stop.')
                end if

                ! Read mass and solute id
                read(inUnit, *) particleMass
                if ( ( simulationData%ParticlesMassOption .eq. 2 ) .or. & 
                     ( simulationData%SolutesOption .eq. 1 ) ) then 
                  ! Read solute id 
                  read(inUnit, *) soluteId
                end if 

                ! Read release times data
                read(inUnit, *) releaseTimeCount, initialReleaseTime, releaseInterval
                call particleGroups(nic)%SetReleaseOption2(initialReleaseTime, &
                  releaseTimeCount, releaseInterval)

                ! It needs something to allocate particles, 
                ! the read and prepare locations or something
                allocate(subDiv(1,12))
                subDiv(:,:) = 0
                read(inUnit, *) (subDiv(n,i), i = 1, 3)
                npcell = subDiv(1,1)*subDiv(1,2)*subDiv(1,3)

                ! Calculate the total number of particles for all release time points.
                totalParticleCount = npcell*particleGroups(nic)%GetReleaseTimeCount()
                if(allocated(particleGroups(nic)%Particles)) deallocate(particleGroups(nic)%Particles)
                allocate(particleGroups(nic)%Particles(totalParticleCount))
                particleGroups(nic)%TotalParticleCount = totalParticleCount

                call pr_CreateParticlesAsInternalArray(    & 
                  particleGroups(nic), injectionCellNumber,&
                  m, subDiv(1,1), subDiv(1,2), subDiv(1,3),&
                  drape(1) )
                ! Assign mass
                particleGroups(nic)%Mass = particleMass
                particleGroups(nic)%Particles(:)%Mass = particleMass
                if ( simulationData%ParticlesMassOption .eq. 2 ) then 
                  ! Assign the solute id 
                  particleGroups(nic)%Solute = soluteId
                end if 

                ! Next

              ! Cell id and a timeseries for concentrations
              ! Note: Times should be consistent with referenceTime
              case (2)

                ! Are these ids as (lay,row,col) or (cellid) ?
                read(inUnit, '(a)', iostat=ioInUnit) line
                icol = 1
                call urword(line, icol, istart, istop, 2, n, r, 0, 0)
                readStyle = n

                ! Read the injection cell
                if( readStyle .eq. 1) then
                  ! Read as layer, row, column
                  read(inUnit, *) layer, row, column
                  call ugetnode(layerCount, rowCount, columnCount, layer, row, column,cellNumber)
                  injectionCellNumber = cellNumber
                else if ( readStyle .eq. 2 ) then 
                  read(inUnit,*) injectionCellNumber
                else
                    call ustop(' Invalid read style for injection cell. Stop.')
                end if

                ! Read mass and solute id
                read(inUnit, *) particleMass
                if ( ( simulationData%ParticlesMassOption .eq. 2 ) .or. & 
                     ( simulationData%SolutesOption .eq. 1 ) ) then 
                  ! Read solute id 
                  read(inUnit, *) soluteId
                end if 

                ! N Timeseries 
                read(inUnit, *) nInjectionTimes

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

                ! Compute delZ
                ! Simple and correct if convertible 
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

                ! Compute cell volume ( is this a vector ?)
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
                templateCount = 1 
                templateCellCount = 1
                allocate(subDiv(templateCount,12))
                allocate(templateSubDivisionTypes(templateCount))
                allocate(templateCellCounts(templateCount))
                allocate(drape(templateCount))
                allocate(templateCellNumbers(templateCellCount))
                drape  = 0 ! SHOULD COME FROM SOMEWHERE
                np     = 0
                npcell = 0
                offset = 0
                subDiv = 0
              
                ! Cell shape factors 
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


                ! Compute release data for each time
                ! Relate injected concentration to 
                ! a number of particles released. Is there 
                ! a frequency ?
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


                ! Calculate the total number of particles for all release time points.
                totalParticleCount = np
                if(allocated(particleGroups(nic)%Particles)) deallocate(particleGroups(nic)%Particles)
                allocate(particleGroups(nic)%Particles(totalParticleCount))
                particleGroups(nic)%TotalParticleCount = totalParticleCount

                ! Update particles mass using the actual injected number
                particleMass = sum( &
                    abs(injectionDensity)*cellVolumes(1)*flowModelData%Porosity( injectionCellNumber ) &
                )/totalParticleCount 
                write(outUnit,'(/A,es18.9e3)') ' Effective particleMass for injection condition = ', particleMass

                m = 0
                do it =1, nInjectionTimes
                  ! Set the data for particles at the first release time point
                  if ( nParticles(it) .lt. 1d0 ) cycle
                  sdiv(1) = int( nParticlesX(it) ) + 1
                  sdiv(2) = int( nParticlesY(it) ) + 1
                  sdiv(3) = int( nParticlesZ(it) ) + 1
                  ! For said rare cases where the injected quantity is negative
                  if ( injectionDensity(it) .gt. 0 ) then 
                      particleMass = abs( particleMass ) 
                  else 
                      particleMass = -1*abs( particleMass )
                  end if 
                  call CreateMassParticlesAsInternalArray(   & 
                    particleGroups(nic), injectionCellNumber,&
                    m, sdiv(1), sdiv(2), sdiv(3),            &
                    drape(1), particleMass,                  &
                    injectionTimes( it )                     & 
                  )
                end do 
                if ( simulationData%ParticlesMassOption .eq. 2 ) then 
                  ! Assign the solute id 
                  particleGroups(nic)%Solute = soluteId
                end if 

                ! Deallocate temporary arrays
                deallocate(subDiv)
                deallocate(templateSubDivisionTypes)
                deallocate(drape)
                deallocate(templateCellCounts)
                deallocate(templateCellNumbers)

            case default
              write(outUnit,*) ' Invalid flux condition kind ', injectionFormat 
              write(outUnit,*) 'Stop.'
              call ustop('Error: Invalid flux condition kind ')          
            end select 

            if ( simulationData%ParticlesMassOption .eq. 2 ) then 
              ! Assign the solute id 
              particleGroups(nic)%Solute = soluteId
            end if 

            ! Increment valid injection conditions
            nValidInjectionConditions = nValidInjectionConditions

            ! Report number of particles
            write(outUnit, '(a,i4,a,i10,a)') ' Flux/injection condition ', n, ' contains ',   &
              particleGroups(nic)%TotalParticleCount, ' particles.'
            particleCount = particleCount + particleGroups(nic)%TotalParticleCount

          end do ! Loop over nInjectionConditions
          write(outUnit, '(a,i10)') ' Total number of particles on injection conditions = ', particleCount
          write(outUnit, *)

          if ( nValidInjectionConditions .gt. 0 ) then 
            ! Extend simulationdata to include new particle groups
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
          end if

          deallocate( particleGroups )

        end if ! if nInjectionConditions .gt. 0


        ! ICBound
        if((grid%GridType .eq. 1) .or. (grid%GridType .eq. 3)) then
          call u3dintmp(inUnit, outUnit, grid%LayerCount, grid%RowCount,      &
            grid%ColumnCount, grid%CellCount, this%ICBound, ANAME(3))
        else if((grid%GridType .eq. 2) .or. (grid%GridType .eq. 4)) then
          call u3dintmpusg(inUnit, outUnit, grid%CellCount, grid%LayerCount,  &
            this%ICBound, aname(3), cellsPerLayer)
        else
          write(outUnit,*) ' Invalid grid type specified when reading ICBOUND array data.'
          write(outUnit,*) ' Stop.'
          call ustop('Error: Invalid grid type when reading ICBOUND data.')          
        end if


        ! Define dispersion according to solutes option 
        select case( simulationData%SolutesOption )
          ! 0: Single solute dispersion 
          case (0)
            ! A single virtual solute storing the same 
            ! transport parameters for all pgroups
            this%BaseSolute%id = 0
            this%BaseSolute%stringid = 'S0' 

            ! Trust at least one particle group
            if ( allocated(this%BaseSolute%pGroups) ) deallocate(this%BaseSolute%pGroups)
            allocate(this%BaseSolute%pGroups(simulationData%ParticleGroupCount))
            do n=1,simulationData%ParticleGroupCount
              this%BaseSolute%pGroups(n) = n 
            end do 
            this%BaseSolute%nParticleGroups = simulationData%ParticleGroupCount

            ! Read dispersion model kind
            read(inUnit, '(a)') line
            icol = 1
            call urword(line, icol, istart, istop, 2, n, r, 0, 0)
            this%BaseSolute%dispersionModel = n 
            ! In the meantime this is neede because 
            ! this variable determines the dispersion 
            ! function 
            trackingOptions%dispersionModel = n 

            ! Read dispersion data
            call this%LoadSoluteDispersion(&
                inUnit, outUnit, this%BaseSolute, grid, &
                          cellsPerLayer, trackingOptions )

            ! Assign dispersion pointers
            ! Something different if nonlinear dispersion ?
            this%AlphaLong => this%BaseSolute%AlphaLong
            this%AlphaTran => this%BaseSolute%AlphaTran

            ! Needs clarification
            this%DMol = this%BaseSolute%dAqueous

            if (simulationData%ParticlesMassOption .eq. 2) then
              ! Create different solutes for id purposes
              ! However, they are all displaced with the 
              ! same transport properties
              ! Extracted from particles groups defined until 
              ! this very moment. 
              continue
            else
              ! If not solute id, then assumes
              ! all pgroups are the same solute
              this%nSolutes = 1
              if ( allocated(this%Solutes) ) deallocate( this%Solutes )
              allocate( this%Solutes(this%nSolutes) )
              this%Solutes = this%BaseSolute
            end if

          ! 1: Multiple solute, multidispersion
          case (1)
            ! How many ?
            read(inUnit, '(a)') line
            icol = 1
            call urword(line, icol, istart, istop, 2, n, r, 0, 0)
            this%nSolutes = n 
            
            ! Lets trust the user, but some checking should be done
            ! on the value of n 
            if ( allocated(this%Solutes) ) deallocate( this%Solutes )
            allocate( this%Solutes(this%nSolutes) )

            ! Initialize these guys
            do ns =1,this%nSolutes
              ! Of course needs some attributes 
              ! of the solute

              ! Read the integer id 
              read(inUnit, '(a)') line
              icol = 1
              call urword(line, icol, istart, istop, 2, n, r, 0, 0)
              this%Solutes(ns)%id = n 
              
              ! Read the string id
              read(inUnit, '(a)') line
              icol = 1
              call urword(line, icol, istart, istop, 0, n, r, 0, 0)
              this%Solutes(ns)%stringid = line(istart:istop)

              ! If the soluteId was assigned to each 
              ! particle group due to particlessmassoption 2 
              ! then respect those assignments
              if ( simulationData%ParticlesMassOption .ne. 2 ) then

                ! Read how many pgroups
                read(inUnit, '(a)') line
                icol = 1
                call urword(line, icol, istart, istop, 2, n, r, 0, 0)
                this%Solutes(ns)%nParticleGroups = n 
                
                if ( allocated( this%Solutes(ns)%pGroups ) ) deallocate( this%Solutes(ns)%pGroups )
                allocate( this%Solutes(ns)%pGroups( &
                   this%Solutes(ns)%nParticleGroups ) )
                
                ! Read related particle groups
                do npg =1,this%Solutes(ns)%nParticleGroups
                  ! Read the pgroups
                  read(inUnit, '(a)') line
                  icol = 1
                  call urword(line, icol, istart, istop, 2, n, r, 0, 0)
                  this%Solutes(ns)%pGroups(npg) = n

                  ! It needs to assign the soluteId back to the 
                  ! corresponding pgroup for simulations 
                  ! where the solute is not specified in the pgroup
                  simulationData%ParticleGroups(&
                      this%Solutes(ns)%pGroups(npg) )%Solute = ns
                end do

              else if ( simulationData%ParticlesMassOption .eq. 2 ) then

                ! Read the INDEXES in simulationData%ParticleGroup list
                ! and count how many for this solute
                ncount = 0
                do npg =1,simulationData%ParticleGroupCount
                  if ( simulationData%ParticleGroups( npg )%Solute .eq. ns ) then
                    ncount = ncount + 1
                  end if 
                end do

                ! If no pgroups, let it continue
                if ( ncount .eq. 0 ) then 
                  write(outUnit,*) 'Warning: no particle groups associated to solute id ', ns
                else
                  ! Initialize pgroups for the solute
                  this%Solutes(ns)%nParticleGroups = ncount
              
                  if ( allocated( this%Solutes(ns)%pGroups ) ) deallocate( this%Solutes(ns)%pGroups )
                  allocate( this%Solutes(ns)%pGroups( &
                     this%Solutes(ns)%nParticleGroups ) )

                  ncount = 0 
                  do npg =1,simulationData%ParticleGroupCount
                    if ( simulationData%ParticleGroups( npg )%Solute .eq. ns ) then
                      ncount = ncount + 1
                      this%Solutes(ns)%pGroups(ncount) = npg
                    end if 
                  end do
                end if 

              end if 
             
              ! Read dispersion model kind
              read(inUnit, '(a)') line
              icol = 1
              call urword(line, icol, istart, istop, 2, n, r, 0, 0)
              this%Solutes(ns)%dispersionModel = n

              ! Read dispersion data
              call this%LoadSoluteDispersion(&
                  inUnit, outUnit, this%Solutes(ns), grid, &
                            cellsPerLayer, trackingOptions )

            end do  

            !! At some point something need to done 
            !! with particlemassoption ? 
            !else if (simulationData%ParticlesMassOption .eq. 2) then
            !   ! Link particle groups to recently created solutes
            !   ! Already done by requesting pgroups 
            !   continue
            !end if

            ! In this case, dispersivities and other 
            ! parameters are asigned during the loop 
            ! over particle groups
          case default
            ! Some error handling
            print *, 'NOT IMPLEMENTED'
        end select


        !! Time Step kind 
        !read(inUnit, '(a)') line
        !icol = 1
        !call urword(line,icol,istart,istop,1,n,r,0,0)
        !line = line(istart:istop)
        !! Give a Courant
        !if ( line .eq. 'CONSTANT_CU' ) then
        !    trackingOptions%timeStepKind = 1
        !    read( inUnit, * ) line
        !    icol = 1
        !    call urword(line,icol,istart,istop,3,n,r,0,0)
        !    trackingOptions%timeStepParameters(1) = r
        !! Give a Peclet
        !else if ( line .eq. 'CONSTANT_PE' ) then
        !    trackingOptions%timeStepKind = 2
        !    read( inUnit, * ) line
        !    icol = 1
        !    call urword(line,icol,istart,istop,3,n,r,0,0)
        !    trackingOptions%timeStepParameters(2) = r
        !! Minimum between Courant and Peclet 
        !else if ( line .eq. 'MIN_ADV_DISP' ) then
        !    trackingOptions%timeStepKind = 3
        !    read( inUnit, * ) line
        !    icol = 1
        !    call urword(line,icol,istart,istop,3,n,r,0,0)
        !    trackingOptions%timeStepParameters(1) = r
        !    read( inUnit, * ) line
        !    icol = 1
        !    call urword(line,icol,istart,istop,3,n,r,0,0)
        !    trackingOptions%timeStepParameters(2) = r
        !else
        !    call ustop('RWPT: Invalid options for time step selection. Stop.')
        !end if


        !! Advection Integration Kind
        !read(inUnit, '(a)', iostat=iodispersion) line
        !if ( iodispersion .lt. 0 ) then 
        !    ! end of file
        !    trackingOptions%advectionKind = 1
        !    write(outUnit,'(A)') 'RWPT: Advection integration not specified. Default to EXPONENTIAL.'
        !else
        !    icol = 1
        !    call urword(line,icol,istart,istop,1,n,r,0,0)
        !    line = line(istart:istop)
        !    select case(line)
        !        case('EXPONENTIAL')
        !            trackingOptions%advectionKind = 1
        !            write(outUnit,'(A)') 'RWPT: Advection integration is EXPONENTIAL.'
        !        case('EULERIAN')
        !            trackingOptions%advectionKind = 2
        !            write(outUnit,'(A)') 'RWPT: Advection integration is EULERIAN.'
        !        case default
        !            trackingOptions%advectionKind = 2
        !            write(outUnit,'(A)') 'RWPT: Advection integration not specified. Default to EULERIAN.'
        !    end select
        !end if


        !! Domain dimensions option
        !read(inUnit, '(a)', iostat=iodispersion) line
        !if ( iodispersion .lt. 0 ) then 
        !    ! end of file
        !    trackingOptions%twoDimensions = .false.
        !    write(outUnit,'(A)') 'RWPT: Number of dimensions not specified. Defaults to 3D.'
        !else
        !    icol = 1
        !    call urword(line,icol,istart,istop,1,n,r,0,0)
        !    line = line(istart:istop)
        !    select case(line)
        !        case('2D')
        !            trackingOptions%twoDimensions = .true.
        !            write(outUnit,'(A)') 'RWPT: Selected 2D domain solver.'
        !        case('3D')
        !            trackingOptions%twoDimensions = .false.
        !            write(outUnit,'(A)') 'RWPT: Selected 3D domain solver.'
        !        case default
        !            trackingOptions%twoDimensions = .false.
        !            write(outUnit,'(A)') 'RWPT: Invalid option for domain solver, defaults to 3D.'
        !    end select
        !end if 
    

        ! Close dispersion data file
        close( inUnit )


    end subroutine pr_ReadData



    
end module TransportModelDataModule
