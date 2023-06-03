module TransportModelDataModule
  use ModflowRectangularGridModule,only : ModflowRectangularGridType
  use FlowModelDataModule,only : FlowModelDataType
  use SoluteModule,only : SoluteType
  use DispersionDataModule,only : DispersionDataType
  use ModpathSimulationDataModule, only: ModpathSimulationDataType
  implicit none
  !---------------------------------------------------------------------------------------------------------------

  ! Set default access status to private
  private

  type,public :: TransportModelDataType

    logical :: Initialized = .false.

    ! Local dispersion parameters
    doubleprecision,dimension(:),pointer :: AlphaL    => null()
    doubleprecision,dimension(:),pointer :: AlphaT    => null()
    doubleprecision,dimension(:),pointer :: AlphaLH   => null()
    doubleprecision,dimension(:),pointer :: AlphaLV   => null()
    doubleprecision,dimension(:),pointer :: AlphaTH   => null()
    doubleprecision,dimension(:),pointer :: AlphaTV   => null()
    doubleprecision,dimension(:),pointer :: DMEff     => null()
    logical :: uniformDispersionParameters

    ! ICBound
    integer,allocatable, dimension(:) :: ICBound
    integer, pointer, dimension(:)    :: ICBoundTS => null()
    integer                           :: defaultICBound ! applied to model boundaries
    logical                           :: followIBoundTS = .false.
    logical                           :: followIBound   = .false.

    ! Simulation data
    class( ModpathSimulationDataType ), pointer :: simulationData

    ! Solutes
    type( SoluteType ), allocatable, dimension(:) :: Solutes
    type( SoluteType ) :: BaseSolute
    integer            :: nSolutes
    logical            :: readSPCPkg

    ! Dispersion models 
    type( DispersionDataType ), allocatable, dimension(:) :: DispersionData 
    type( DispersionDataType ) :: BaseDispersionData
    integer                    :: nDispersion
    integer                    :: currentDispersionModelKind

    ! grid
    class(ModflowRectangularGridType),pointer :: Grid => null()

    ! Private variables
    integer,private :: CurrentStressPeriod = 0
    integer,private :: CurrentTimeStep = 0

    ! Flow model data
    type( FlowModelDataType ), pointer :: flowModelData => null()

    ! Budget reader
    !type(BudgetReaderType), pointer :: budgetReader => null()
    !this%BudgetReader => budgetReader

    ! TO BE DEPRECATED
    doubleprecision :: DMol

  contains

    procedure :: Initialize=>pr_Initialize
    procedure :: Reset=>pr_Reset
    procedure :: ReadIMPData=>pr_ReadIMPData
    procedure :: ReadSPCData=>pr_ReadSPCData
    procedure :: ReadDSPData=>pr_ReadDSPData
    procedure :: ValidateDataRelations => pr_ValidateDataRelations
    procedure :: LoadTimeStep=>pr_LoadTimeStep
    procedure :: SetSoluteDispersion=>pr_SetSoluteDispersion

  end type


contains


  subroutine pr_Initialize( this, grid, simulationData, flowModelData )
  !---------------------------------------------------------------------------------------------------------------
  ! Specifications
  !---------------------------------------------------------------------------------------------------------------
  implicit none
  class(TransportModelDataType) :: this
  class(ModflowRectangularGridType), intent(in), target :: grid
  type(ModpathSimulationDataType), intent(in), target :: simulationData
  type(FlowModelDataType), intent(in), target :: flowModelData

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

      ! flowModelData
      this%flowModelData => flowModelData

      ! IBounds
      !if(allocated(this%ICBoundTS)) deallocate(this%ICBoundTS)
      !allocate(this%ICBoundTS(grid%CellCount))
      !if(allocated(this%ICBound)) deallocate(this%ICBound)
      !allocate(this%ICBound(grid%CellCount))
      !this%ICBound => simulationData%ICBound

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
     
    this%Grid   => null()
    this%AlphaL => null()
    this%AlphaT => null()
    this%DMEff  => null()
    this%DMol   = 0

    if(allocated(this%ICBound)) deallocate( this%ICBound )
    this%ICBoundTS => null()

  end subroutine pr_Reset


  ! Read specific IMP data
  subroutine pr_ReadIMPData( this, impFile, impUnit, outUnit, grid )
    use UTL8MODULE,only : urword,ustop,u3dintmpusg,u3dintmp,ugetnode
    !--------------------------------------------------------------
    ! Specifications
    !--------------------------------------------------------------
    implicit none
    ! input 
    class(TransportModelDataType), target        :: this
    character(len=200), intent(in)               :: impFile
    integer, intent(in)                          :: impUnit
    integer, intent(in)                          :: outUnit
    class(ModflowRectangularGridType),intent(in) :: grid
    ! local
    integer            :: isThisFileOpen
    integer            :: icol,istart,istop,n
    integer            :: impFormat
    doubleprecision    :: r
    character(len=200) :: line
    integer, dimension(:), allocatable :: cellsPerLayer
    ! icbound
    character(len=24),dimension(1) :: aname
    data aname(1)   /'      IMPCELLS'/
    !--------------------------------------------------------------

    write(outUnit, *)
    write(outUnit, '(1x,a)') 'MODPATH-RW IMP file data'
    write(outUnit, '(1x,a)') '------------------------'

    ! CellsPerLayer, required for u3d reader
    allocate(cellsPerLayer(grid%LayerCount))
    do n = 1, grid%LayerCount
      cellsPerLayer(n) = grid%GetLayerCellCount(n)
    end do

    ! Initialize 
    this%followIBound = .false.
    this%followIBoundTS = .false.

    ! Allocate ICBound array, is needed downstream
    if(allocated(this%ICBound)) deallocate(this%ICBound)
    allocate(this%ICBound(grid%CellCount))

    ! Verify if unit is open 
    isThisFileOpen = -1
    inquire( file=impFile, number=isThisFileOpen )
    if ( isThisFileOpen .lt. 0 ) then 
      ! No file
      write(outUnit,'(A)') 'IMP package was not specified in name file.'
      ! Initialize ICBound with only zeroes
      this%ICBound(:) = 0
      ! Link pointer
      this%ICBoundTS => this%ICBound
      ! Assign defaultICBound as rebound
      this%defaultICBound = 1
      ! And leave
      return
    end if

    write(outUnit,'(A)') 'Will interpret specification of impermeable cells.'

    ! Read defaultICBound
    read(impUnit, '(a)') line
    icol = 1
    call urword(line, icol, istart, istop, 2, n, r, 0, 0)
    if ( n.ge.1) then 
     this%defaultICBound = 1
     write(outUnit,'(A,I1)') 'Default ICBOUND set to rebound when no cell connection: ', 1
    else
     this%defaultICBound = 0
     write(outUnit,'(A,I1)') 'Default ICBOUND set to open when no cell connection: ', 0
    end if
    
    ! Read impFormat
    read(impUnit, '(a)') line
    icol = 1
    call urword(line, icol, istart, istop, 2, n, r, 0, 0)
    ! Validate
    select case(n)
    case(0,1,2)
      continue
    case default
      write(outUnit,'(A)') 'Invalid format for IMP package. Stop.'
      call ustop('Invalid format for IMP package. Stop.')          
    end select
    impFormat = n

    select case(impFormat)
    case(0)
     ! From IBound assign into ICBound and assign pointer, does not 
     ! change in time
     this%followIBound = .true.
     ! Initialize with zeroes
     this%ICBound(:) = 0
     ! If flow-model inactive, transport rebound
     do n=1,grid%CellCount
       if ( this%flowModelData%IBound(n) .lt. 1 ) this%ICBound(n) = 1
     end do 
     this%ICBoundTS => this%ICBound
     write(outUnit,'(A)') 'Impermeable cells follow inactive cells at flowModelData%IBound.'
    case(1)
     ! Activate flag, ICBoundTS is updated every time step
     this%followIBoundTS = .true.
     ! Initialize with zeroes
     this%ICBound(:) = 0
     write(outUnit,'(A)') 'Impermeable cells follow inactive cells at flowModelData%IBoundTS, includes dry cells.'
    case(2)
     write(outUnit,'(A)') 'Impermeable cells read from IMPCELLS.'
     ! Read IMPCELLS
     if((grid%GridType .eq. 1) .or. (grid%GridType .eq. 3)) then
       call u3dintmp(impUnit, outUnit, grid%LayerCount, grid%RowCount,      &
         grid%ColumnCount, grid%CellCount, this%ICBound, aname(1))
     else if((grid%GridType .eq. 2) .or. (grid%GridType .eq. 4)) then
       call u3dintmpusg(impUnit, outUnit, grid%CellCount, grid%LayerCount,  &
         this%ICBound, aname(1), cellsPerLayer)
     else
       write(outUnit,*) 'Invalid grid type specified when reading IMPCELLS array data.'
       write(outUnit,*) 'Stopping.'
       call ustop(' ')          
     end if

     ! Check that not all cells are impermeable
     if ( all(this%ICBound(:).gt.0) ) then 
       write(outUnit,'(A)')'Invalid IMP specification. All cells were marked as impermable. Stop.'
       call ustop('Invalid IMP specification. All cells were marked as impermable. Stop.')          
     end if 

     ! Assign pointer, does not change in time.
     this%ICBoundTS => this%ICBound

    end select


    ! Close file
    close( impUnit )


  end subroutine pr_ReadIMPData


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
    integer :: isThisFileOpen 
    integer :: icol,istart,istop,n
    doubleprecision    :: r
    character(len=200) :: line
    integer :: ns, ncount, npg, pgcount, pgindex
    integer, allocatable, dimension(:) :: pGroupsIds
    integer, allocatable, dimension(:) :: soluteIds
    integer :: minSolId, maxSolId, countSolId
    character(len=20) :: tempChar1
    !--------------------------------------------------------------

    write(outUnit, *)
    write(outUnit, '(1x,a)') 'MODPATH-RW SPC file data'
    write(outUnit, '(1x,a)') '------------------------'

    ! Restart flag
    this%readSPCPkg = .false.

    ! Local pointer to class pointer 
    simulationData => this%simulationData

    ! Verify if unit is open 
    isThisFileOpen = -1 
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
        write(outUnit,'(A)') 'ParticlesMassOption.eq.2. Will read soluteIds from particle groups and leave.'

        if ( allocated( soluteIds ) ) deallocate( soluteIds )
        allocate( soluteIds(simulationData%ParticleGroupCount) )
        do npg=1,simulationData%ParticleGroupCount
          soluteIds(npg) = simulationData%ParticleGroups(npg)%Solute
        end do
        minSolId = minval(soluteIds)-1
        maxSolId = maxval(soluteIds)
        ns = 0
        ! Count how many different ids
        do while (minSolId<maxSolId)
          ns = ns+1
          minSolId = minval(soluteIds, mask=soluteIds>minSolId)
        end do
        this%nSolutes = ns
        if ( allocated(this%Solutes) ) deallocate( this%Solutes )
        allocate( this%Solutes(this%nSolutes) )

        ! Now assign in increasing id order 
        minSolId = minval(soluteIds)-1
        ns = 0
        do while (minSolId<maxSolId)
          ns = ns+1
          minSolId = minval(soluteIds, mask=soluteIds>minSolId)
          countSolId = count(soluteIds.eq.minSolId)
          if ( allocated(this%Solutes(ns)%pGroups) ) deallocate(this%Solutes(ns)%pGroups)
          allocate(this%Solutes(ns)%pGroups(countSolId))
          pgcount = 0
          do npg=1,simulationData%ParticleGroupCount
            if ( simulationData%ParticleGroups(npg)%Solute.eq.minSolId ) then 
              pgcount = pgcount + 1 
              this%Solutes(ns)%pGroups(pgcount) = simulationData%ParticleGroups(npg)%Group
            end if 
          end do 
          this%Solutes(ns)%nParticleGroups = pgcount
          this%Solutes(ns)%id = minSolId
          ! Assigned increasingly ?
          !this%Solutes(ns)%id = ns
          this%Solutes(ns)%userid = minSolId
          tempChar1 = ''
          write( unit=tempChar1, fmt=* )ns 
          this%Solutes(ns)%stringid = 'SPC'//trim(adjustl(tempChar1))
        end do
      else
        write(outUnit,'(A)') 'ParticlesMassOption.ne.2. Will create a single solute and leave.'

        ! If not solute id, then assumes single solute and
        ! all pgroups are of the same kind 
        this%nSolutes = 1
        if ( allocated(this%Solutes) ) deallocate( this%Solutes )
        allocate( this%Solutes(this%nSolutes) )

        ! A single virtual solute storing the same 
        ! transport parameters for all pgroups
        this%Solutes(1)%id = 1
        this%Solutes(1)%userid = 1
        this%Solutes(1)%stringid = 'SPC1'

        ! All pgroups to base solute
        ! The way to identify the pGroup is by their position in the pgroups array 
        if ( allocated(this%Solutes(1)%pGroups) ) deallocate(this%Solutes(1)%pGroups)
        allocate(this%Solutes(1)%pGroups(simulationData%ParticleGroupCount))
        do n=1,simulationData%ParticleGroupCount
          this%Solutes(1)%pGroups(n) = n 
        end do 
        this%Solutes(1)%nParticleGroups = simulationData%ParticleGroupCount

      end if

      ! Leave
      return

    end if

    ! If the SPC package was specified, 
    ! then read it. 
    this%readSPCPkg = .true.
    if( allocated(pGroupsIds) ) deallocate( pGroupsIds )
    ! Interpret solutesoption to determined 
    ! whether the simulation is multidispersion or not

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
      call ustop('Number of given species/solutes is .le. 0.')
    end if

    ! Report
    if ( simulationData%ParticlesMassOption .ne. 2 ) then
     write(outUnit,'(A)') 'ParticleMassOption.ne.2: Read pgroups related to each species from the pkg.'
    else if ( simulationData%ParticlesMassOption .eq. 2 ) then
     write(outUnit,'(A)') 'ParticleMassOption.eq.2: Interpret pgroups related to each species from particle groups.'
    end if 

    ! Allocate solutes
    if ( allocated(this%Solutes) ) deallocate( this%Solutes )
    allocate( this%Solutes(this%nSolutes) )

    ! Initialize these guys
    do ns =1,this%nSolutes
      write(outUnit,'(A,I4)') 'Processing species specification ', ns

      ! Internal id assigned increasingly
      this%Solutes(ns)%id = ns 
      this%Solutes(ns)%userid = ns
      
      ! Read the string id
      read(spcUnit, '(a)') line
      icol = 1
      call urword(line, icol, istart, istop, 0, n, r, 0, 0)
      this%Solutes(ns)%stringid = line(istart:istop)

      ! Assign pgroups related to the solute
      if ( simulationData%ParticlesMassOption .ne. 2 ) then
        ! Read pgroups related to the solute from the pkg

        if ( .not. allocated( pGroupsIds ) ) then 
          allocate( pGroupsIds(simulationData%ParticleGroupCount) ) 
          do npg=1,simulationData%ParticleGroupCount
            pGroupsIds(npg) = simulationData%ParticleGroups(npg)%Group
          end do
        end if

        ! How many pgroups
        read(spcUnit, '(a)') line
        icol = 1
        call urword(line, icol, istart, istop, 2, n, r, 0, 0)
        if ( n.lt.1 ) then 
          write(outUnit,'(A)') 'Given number of particle groups for species is less than 1. Stop.'
          call ustop('Given number of particle groups for specie is less than 1. Stop.')
        end if 
        this%Solutes(ns)%nParticleGroups = n 
        
        ! Allocate species particle groups
        if ( allocated( this%Solutes(ns)%pGroups ) ) deallocate( this%Solutes(ns)%pGroups )
        allocate(this%Solutes(ns)%pGroups(this%Solutes(ns)%nParticleGroups))
        
        ! Read particle groups 
        read(spcUnit, *) (this%Solutes(ns)%pGroups(npg), npg = 1,this%Solutes(ns)%nParticleGroups)

        ! It needs to assign the soluteId back to the 
        ! corresponding pgroup for simulations 
        ! where the solute is not specified in the pgroup
        do npg = 1,this%Solutes(ns)%nParticleGroups
          pgindex = findloc( pGroupsIds, this%Solutes(ns)%pGroups(npg), 1 )
          if ( pgindex .eq. 0 ) then 
            write(outUnit,'(a)')'Given particle group id was not found in the list of avaiable particle groups. Stop.'
            call ustop('Given particle group id was not found in the list of avaiable particle groups. Stop.')
          end if 
          simulationData%ParticleGroups(&
              this%Solutes(ns)%pGroups(npg) )%Solute = ns
        end do


      else if ( simulationData%ParticlesMassOption .eq. 2 ) then
        ! Read pgroups related to the solute ids at pgroups

        ! Read the soluteId's from simulationData%ParticleGroup list
        ! and count how many for this solute
        ncount = 0
        do npg =1,simulationData%ParticleGroupCount
          if ( (simulationData%ParticleGroups( npg )%Solute .lt. 1) .or. &
               (simulationData%ParticleGroups( npg )%Solute .gt. this%nSolutes ) ) then 
            write(outUnit,'(a)')'Given SpeciesID for particle group is outside the valid range. Stop.'
            call ustop('Given SpeciesID for particle group is outside the valid range. Stop.')
          end if 
          if ( simulationData%ParticleGroups( npg )%Solute .eq. ns ) then
            ncount = ncount + 1
          end if 
        end do

        ! If no pgroups, let it continue ?
        if ( ncount .gt. 0 ) then 
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
        else
          ! Will it happen ?
          write(outUnit,*) 'Warning: no particle groups associated to SpeciesID ', ns
          this%Solutes(ns)%nParticleGroups = 0
        end if
        
      end if

      ! If simulation is multispecies dispersiom
      if ( simulationData%SolutesOption.eq.1 ) then
       ! Read the dispersion string id
       read(spcUnit, '(a)') line
       icol = 1
       call urword(line, icol, istart, istop, 0, n, r, 0, 0)
       this%Solutes(ns)%dispersionStringid = line(istart:istop)

       write(outUnit,'(A,A,A,A)')&
       'Species ',trim(adjustl(this%Solutes(ns)%stringid)),&
       ' will be related to dispersion ',trim(adjustl(this%Solutes(ns)%dispersionStringid))
      end if 

    end do ! ns =1,this%nSolutes


    ! Some verification that not all species 
    ! remained without particle groups
    ncount = 0
    do ns=1,this%nSolutes
      if( this%Solutes(ns)%nParticleGroups .eq. 0 ) ncount = ncount + 1
    end do
    if ( ncount .eq. this%nSolutes ) then 
      write(outUnit,'(a)') 'Error: all the species remained without particle groups. Stop.'
      call ustop('Error: all the species remained without particle groups. Stop.')
    end if  


    ! Or some report on relations
        

    ! Close spc data file
    close( spcUnit )


  end subroutine pr_ReadSPCData


  ! Read specific DSP data
  subroutine pr_ReadDSPData( this, dspFile, dspUnit, outUnit )
    use UTL8MODULE,only : urword,ustop,u3ddblmpusg, u3ddblmp
    use UtilMiscModule,only : TrimAll
    !--------------------------------------------------------------
    ! Specifications
    !--------------------------------------------------------------
    implicit none
    ! input 
    class(TransportModelDataType), target    :: this
    character(len=200), intent(in)           :: dspFile
    integer, intent(in)                      :: dspUnit
    integer, intent(in)                      :: outUnit
    ! local
    type(ModpathSimulationDataType), pointer  :: simulationData
    class(ModflowRectangularGridType),pointer :: grid
    type(DispersionDataType),pointer          :: disp
    integer :: isThisFileOpen 
    integer :: icol,istart,istop,n
    doubleprecision    :: r
    character(len=200) :: line
    integer :: nDispersion, ndis
    integer, dimension(:), allocatable :: cellsPerLayer
    character(len=24),dimension(5) :: anamelin
    data anamelin(1) /'          DMEFF'/ 
    data anamelin(2) /'         ALPHAL'/
    data anamelin(3) /'        ALPHALV'/
    data anamelin(4) /'        ALPHATH'/
    data anamelin(5) /'        ALPHATV'/
    character(len=24),dimension(6) :: anamenlin
    data anamenlin(1) /'          DMEFF'/
    data anamenlin(2) /'          BETAL'/
    data anamenlin(3) /'         BETATH'/
    data anamenlin(4) /'         BETATV'/
    data anamenlin(5) /'          DELTA'/
    data anamenlin(6) /'         DGRAIN'/
    character(len=20), allocatable, dimension(:) :: dispnames
    integer :: dnindex 
    !--------------------------------------------------------------

    write(outUnit, *)
    write(outUnit, '(1x,a)') 'MODPATH-RW DSP file data'
    write(outUnit, '(1x,a)') '------------------------'

    ! Verify if unit is open
    isThisFileOpen = -1 
    inquire( file=dspFile, number=isThisFileOpen )
    if ( isThisFileOpen .lt. 0 ) then 
      ! No spc file
      write(outUnit,'(A)') 'DSP package was not specified in name file and is required for RW.'
      call ustop('DSP package was not specified in name file and is required for RW. Stop.')
    end if

    ! Local pointers to class pointers
    simulationData => this%simulationData
    grid => this%Grid


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

    if ( allocated( dispnames ) ) deallocate( dispnames )
    allocate( dispnames( this%nDispersion ) )
    dispnames(:) = ''

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

      ! Dispersion id assigned increasingly
      disp%id = ndis 
      
      ! Read the string id
      read(dspUnit, '(a)') line
      icol = 1
      call urword(line, icol, istart, istop, 0, n, r, 0, 0)
      ! Validate that it was not already given
      dnindex = 0
      dnindex = findloc( dispnames, line(istart:istop), 1 )
      if ( dnindex.gt.0 ) then 
       write(outUnit,'(a,a,a)') 'Dispersion name ', trim(adjustl(line(istart:istop))),' was already given  Stop.'
       call ustop('Dispersion name was already given. Stop.')
      end if
      dispnames(ndis) = line(istart:istop)
      disp%stringid   = line(istart:istop)

      ! Read the dispersion modelKind
      read(dspUnit, '(a)') line
      icol = 1
      call urword(line, icol, istart, istop, 2, n, r, 0, 0)
      disp%modelKind = n 
      ! Validate dispersion model 
      select case ( disp%modelKind )
        ! Linear, isotropic and axisymmetric
        case(0,1)
          continue
        case default
          write(outUnit,*) 'Invalid dispersion model, not implemented. Given ', disp%modelKind
          call ustop('Invalid dispersion model, not implemented. Stop')          
      end select

      ! Read the parameters input format 
      ! 0: uniform 
      ! 1: distributed, u3d reader
      call urword(line, icol, istart, istop, 2, n, r, 0, 0)
      select case(n) 
      case(0,1)
        continue
      case default
        write(outUnit,*) 'Invalid input format for dispersion parameters. Given ', n
        call ustop('Invalid input format for dispersion parameters. Stop')
      end select
      disp%parametersFormat = n 

      ! Allocate
      call disp%InitializeParameters( grid%CellCount )

      ! Read dispersion data
      select case ( disp%modelKind )
       ! Linear, isotropic
       case(0)
        write(outUnit,'(A)') 'Dispersion model is linear, isotropic: will read dmeff and dispersivities.'

        ! Read parameters, according to the given input format
        select case(disp%parametersFormat)
        case (0)
         ! The simplest uniform constant
         write(outUnit,'(A)') 'Dispersion parameters are spatially uniform.'

         ! Read DMEFF
         read(dspUnit, '(a)') line
         icol = 1
         call urword(line, icol, istart, istop, 3, n, r, 0, 0)
         if ( r.lt.0d0 ) then 
          call ustop('Invalid value for effective molecular diffusion, is .lt. 0. Stop.')
         end if 
         disp%DMEff(1) = r

         ! Read ALPHAL
         read(dspUnit, '(a)') line
         icol = 1
         call urword(line, icol, istart, istop, 3, n, r, 0, 0)
         if ( r.lt.0d0 ) then 
          call ustop('Invalid value for alphaL, is .lt. 0. Stop.')
         end if 
         disp%AlphaL(1) = r

         ! Read ALPHATH
         read(dspUnit, '(a)') line
         icol = 1
         call urword(line, icol, istart, istop, 3, n, r, 0, 0)
         if ( r.lt.0d0 ) then 
          call ustop('Invalid value for alphaT, is .lt. 0. Stop.')
         end if 
         disp%AlphaTH(1) = r

        case (1)
         ! Distributed, u3d reader
         write(outUnit,'(A)') 'Dispersion parameters are spatially distributed.'

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
             grid%ColumnCount, grid%CellCount, disp%AlphaTH, anamelin(4))
         else if((grid%GridType .eq. 2) .or. (grid%GridType .eq. 4)) then
           call u3ddblmpusg(dspUnit, outUnit, grid%CellCount, grid%LayerCount,  &
             disp%AlphaTH, anamelin(4), cellsPerLayer)
         else
           write(outUnit,*) 'Invalid grid type specified when reading ALPHATH array data.'
           write(outUnit,*) 'Stopping.'
           call ustop(' ')          
         end if

        end select

       ! Linear, axisymmetric
       case(1)
        write(outUnit,'(A)') 'Dispersion model is linear, axisymmetric: will read dmeff and oriented dispersivities.'

        ! Read parameters, according to the given input format
        select case(disp%parametersFormat)
        case (0)
         ! The simplest uniform constant
         write(outUnit,'(A)') 'Dispersion parameters are spatially uniform.'

         ! Read DMEFF
         read(dspUnit, '(a)') line
         icol = 1
         call urword(line, icol, istart, istop, 3, n, r, 0, 0)
         if ( r.lt.0d0 ) then 
          call ustop('Invalid value for effective molecular diffusion, is .lt. 0. Stop.')
         end if 
         disp%DMEff(1) = r

         ! Read ALPHALH
         read(dspUnit, '(a)') line
         icol = 1
         call urword(line, icol, istart, istop, 3, n, r, 0, 0)
         if ( r.lt.0d0 ) then 
          call ustop('Invalid value for alphaLH, is .lt. 0. Stop.')
         end if 
         disp%AlphaL(1) = r

         ! Read ALPHALV
         read(dspUnit, '(a)') line
         icol = 1
         call urword(line, icol, istart, istop, 3, n, r, 0, 0)
         if ( r.lt.0d0 ) then 
          call ustop('Invalid value for alphaLV, is .lt. 0. Stop.')
         end if 
         disp%AlphaLV(1) = r

         ! Read ALPHATH
         read(dspUnit, '(a)') line
         icol = 1
         call urword(line, icol, istart, istop, 3, n, r, 0, 0)
         if ( r.lt.0d0 ) then 
          call ustop('Invalid value for alphaTH, is .lt. 0. Stop.')
         end if 
         disp%AlphaTH(1) = r

         ! Read ALPHATV
         read(dspUnit, '(a)') line
         icol = 1
         call urword(line, icol, istart, istop, 3, n, r, 0, 0)
         if ( r.lt.0d0 ) then 
          call ustop('Invalid value for alphaTV, is .lt. 0. Stop.')
         end if 
         disp%AlphaTV(1) = r

        case (1)
         ! Distributed, u3d reader
         write(outUnit,'(A)') 'Dispersion parameters are spatially distributed.'

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

         ! Read ALPHALV
         if((grid%GridType .eq. 1) .or. (grid%GridType .eq. 3)) then
           call u3ddblmp(dspUnit, outUnit, grid%LayerCount, grid%RowCount,      &
             grid%ColumnCount, grid%CellCount, disp%AlphaLV, anamelin(3))
         else if((grid%GridType .eq. 2) .or. (grid%GridType .eq. 4)) then
           call u3ddblmpusg(dspUnit, outUnit, grid%CellCount, grid%LayerCount,  &
             disp%AlphaLV, anamelin(3), cellsPerLayer)
         else
           write(outUnit,*) 'Invalid grid type specified when reading ALPHALV array data.'
           write(outUnit,*) 'Stopping.'
           call ustop(' ')          
         end if

         ! Read ALPHATH
         if((grid%GridType .eq. 1) .or. (grid%GridType .eq. 3)) then
           call u3ddblmp(dspUnit, outUnit, grid%LayerCount, grid%RowCount,      &
             grid%ColumnCount, grid%CellCount, disp%AlphaTH, anamelin(4))
         else if((grid%GridType .eq. 2) .or. (grid%GridType .eq. 4)) then
           call u3ddblmpusg(dspUnit, outUnit, grid%CellCount, grid%LayerCount,  &
             disp%AlphaTH, anamelin(4), cellsPerLayer)
         else
           write(outUnit,*) 'Invalid grid type specified when reading ALPHATH array data.'
           write(outUnit,*) 'Stopping.'
           call ustop(' ')          
         end if

         ! Read ALPHATV
         if((grid%GridType .eq. 1) .or. (grid%GridType .eq. 3)) then
           call u3ddblmp(dspUnit, outUnit, grid%LayerCount, grid%RowCount,      &
             grid%ColumnCount, grid%CellCount, disp%AlphaTV, anamelin(5))
         else if((grid%GridType .eq. 2) .or. (grid%GridType .eq. 4)) then
           call u3ddblmpusg(dspUnit, outUnit, grid%CellCount, grid%LayerCount,  &
             disp%AlphaTV, anamelin(5), cellsPerLayer)
         else
           write(outUnit,*) 'Invalid grid type specified when reading ALPHATV array data.'
           write(outUnit,*) 'Stopping.'
           call ustop(' ')          
         end if

        end select

       !! Nonlinear
       !case(2)
       !  write(outUnit,'(A)') 'Dispersion model is nonlinear.'

       !  ! Aqueous molecular diffusion 
       !  read( dspUnit, * ) line
       !  icol = 1
       !  call urword(line,icol,istart,istop,3,n,r,0,0)
       !  disp%dmaqueous = r

       !  ! Read DMEFF
       !  if((grid%GridType .eq. 1) .or. (grid%GridType .eq. 3)) then
       !    call u3ddblmp(dspUnit, outUnit, grid%LayerCount, grid%RowCount,      &
       !      grid%ColumnCount, grid%CellCount, disp%DMEff, anamenlin(1))                      
       !  else if((grid%GridType .eq. 2) .or. (grid%GridType .eq. 4)) then
       !    call u3ddblmpusg(dspUnit, outUnit, grid%CellCount, grid%LayerCount,  &
       !      disp%DMEff, anamenlin(1), cellsPerLayer)
       !  else
       !    write(outUnit,*) 'Invalid grid type specified when reading DMEFF array data.'
       !    write(outUnit,*) 'Stopping.'
       !    call ustop(' ')          
       !  end if

       !  ! Read BETAL
       !  if((grid%GridType .eq. 1) .or. (grid%GridType .eq. 3)) then
       !    call u3ddblmp(dspUnit, outUnit, grid%LayerCount, grid%RowCount,      &
       !      grid%ColumnCount, grid%CellCount, disp%BetaL, anamenlin(2))
       !  else if((grid%GridType .eq. 2) .or. (grid%GridType .eq. 4)) then
       !    call u3ddblmpusg(dspUnit, outUnit, grid%CellCount, grid%LayerCount,  &
       !      disp%BetaL, anamenlin(2), cellsPerLayer)
       !  else
       !    write(outUnit,*) 'Invalid grid type specified when reading BETAL array data.'
       !    write(outUnit,*) 'Stopping.'
       !    call ustop(' ')          
       !  end if

       !  ! Read BETATH
       !  if((grid%GridType .eq. 1) .or. (grid%GridType .eq. 3)) then
       !    call u3ddblmp(dspUnit, outUnit, grid%LayerCount, grid%RowCount,      &
       !      grid%ColumnCount, grid%CellCount, disp%BetaTH, anamenlin(3))
       !  else if((grid%GridType .eq. 2) .or. (grid%GridType .eq. 4)) then
       !    call u3ddblmpusg(dspUnit, outUnit, grid%CellCount, grid%LayerCount,  &
       !      disp%BetaTH, anamenlin(3), cellsPerLayer)
       !  else
       !    write(outUnit,*) 'Invalid grid type specified when reading BETATH array data.'
       !    write(outUnit,*) 'Stopping.'
       !    call ustop(' ')          
       !  end if

       !  ! Read BETATV
       !  if((grid%GridType .eq. 1) .or. (grid%GridType .eq. 3)) then
       !    call u3ddblmp(dspUnit, outUnit, grid%LayerCount, grid%RowCount,      &
       !      grid%ColumnCount, grid%CellCount, disp%BetaTV, anamenlin(4))
       !  else if((grid%GridType .eq. 2) .or. (grid%GridType .eq. 4)) then
       !    call u3ddblmpusg(dspUnit, outUnit, grid%CellCount, grid%LayerCount,  &
       !      disp%BetaTV, anamenlin(4), cellsPerLayer)
       !  else
       !    write(outUnit,*) 'Invalid grid type specified when reading BETATV array data.'
       !    write(outUnit,*) 'Stopping.'
       !    call ustop(' ')          
       !  end if

       !  ! Read DELTA
       !  if((grid%GridType .eq. 1) .or. (grid%GridType .eq. 3)) then
       !    call u3ddblmp(dspUnit, outUnit, grid%LayerCount, grid%RowCount,      &
       !      grid%ColumnCount, grid%CellCount, disp%Delta, anamenlin(5)) 
       !  else if((grid%GridType .eq. 2) .or. (grid%GridType .eq. 4)) then
       !    call u3ddblmpusg(dspUnit, outUnit, grid%CellCount, grid%LayerCount,  &
       !      disp%Delta, anamenlin(5), cellsPerLayer)
       !  else
       !    write(outUnit,*) 'Invalid grid type specified when reading DELTA array data.'
       !    write(outUnit,*) 'Stopping.'
       !    call ustop(' ')          
       !  end if

       !  ! Read DGRAIN
       !  if((grid%GridType .eq. 1) .or. (grid%GridType .eq. 3)) then
       !    call u3ddblmp(dspUnit, outUnit, grid%LayerCount, grid%RowCount,      &
       !      grid%ColumnCount, grid%CellCount, disp%DGrain, anamenlin(6)) 
       !  else if((grid%GridType .eq. 2) .or. (grid%GridType .eq. 4)) then
       !    call u3ddblmpusg(dspUnit, outUnit, grid%CellCount, grid%LayerCount,  &
       !      disp%DGrain, anamenlin(6), cellsPerLayer)
       !  else
       !    write(outUnit,*) 'Invalid grid type specified when reading DGRAIN array data.'
       !    write(outUnit,*) 'Stopping.'
       !    call ustop(' ')          
       !  end if

       case default
         write(outUnit,*) 'Invalid dispersion model. Given ', disp%modelKind
         call ustop('Invalid dispersion model. Stop')          
      end select

    end do ! end do this%nDispersion


    ! Close dsp data file
    close( dspUnit )


  end subroutine pr_ReadDSPData


  subroutine pr_ValidateDataRelations(this, outUnit)
    use UTL8MODULE,only : ustop
    use UtilMiscModule,only : TrimAll
    !--------------------------------------------------------------
    ! Specifications
    !--------------------------------------------------------------
    ! - Verify whether SPC and DSP data relations are consistent,
    !   meaning, that given dispersionIds while defining solutes 
    !   exist in this%DispersionData
    ! - Assign solute%dispersion pointer to corresponding 
    !   dispersiondata model
    !--------------------------------------------------------------
    implicit none
    ! input
    class(TransportModelDataType), target    :: this
    integer, intent(in)                      :: outUnit
    ! local 
    class(ModpathSimulationDataType), pointer :: simulationData
    integer :: ns, ndis
    integer :: fnb1, fnb2, tl1, tl2, lnb1, lnb2
    character(len=300) :: dsid
    logical :: validRelation 
    !--------------------------------------------------------------

    write(outUnit, *)
    write(outUnit, '(1x,a)') 'TransportModelData: verify SPC-DSP relations '
    write(outUnit, '(1x,a)') '---------------------------------------------'

    simulationData => this%simulationData

    if ( simulationData%SolutesOption .eq. 0 ) then
     ! If simulation is single dispersion 
     write(outUnit,'(A)') 'Simulation is single dispersion. Is consistent.'

     ! Pointer names may change
     this%Solutes(1)%Dispersion => this%DispersionData(1)
     this%currentDispersionModelKind = this%Solutes(1)%Dispersion%modelKind

     select case(this%currentDispersionModelKind)
     case(0)
       this%AlphaL => this%Solutes(1)%Dispersion%AlphaL
       this%AlphaT => this%Solutes(1)%Dispersion%AlphaTH
       this%DMEff  => this%Solutes(1)%Dispersion%DMEff
     case(1)
       this%AlphaL  => this%Solutes(1)%Dispersion%AlphaL
       this%AlphaLV => this%Solutes(1)%Dispersion%AlphaLV
       this%AlphaTH => this%Solutes(1)%Dispersion%AlphaTH
       this%AlphaTV => this%Solutes(1)%Dispersion%AlphaTV
       this%DMEff   => this%Solutes(1)%Dispersion%DMEff
     end select

     ! Inform if parameters are spatially uniform
     this%uniformDispersionParameters = .false.
     if ( this%Solutes(1)%Dispersion%uniformParameters ) then
       this%uniformDispersionParameters = .true.
     end if 
     this%DMol = this%Solutes(1)%Dispersion%dmaqueous ! TO BE DEPRECATED

     ! Is valid
     return

    else if ( simulationData%SolutesOption .eq. 1 ) then

     ! If it was specified to be multi species dispersion
     ! but no SPC package was read, then defaults to 
     ! single dispersion and go on, or kill ?
     if ( .not. this%readSPCPkg ) then
      write(outUnit,'(A)')'Simulation is marked to be multispecies dispersion but SPC pkg was not given. Stop.'
      call ustop('Simulation is marked to be multispecies dispersion but SPC pkg was not given. Stop.')
     end if

     ! If simulation is multidispersion
     write(outUnit,'(A)') 'Simulation is multispecies dispersion, will validate specifications.'

     ! Should confirm that dispersion models 
     ! specified for the different solutes
     ! exist in dispersiondata array

     do ns = 1, this%nSolutes

      dsid = this%Solutes(ns)%dispersionStringId
      call TrimAll(dsid, fnb1, lnb1, tl1)

      validRelation = .false.
      ! Loop over dispersion models
      do ndis = 1, this%nDispersion
       call TrimAll( this%DispersionData(ndis)%stringid, fnb2, lnb2, tl2 )
       if (&
       dsid(fnb1:lnb1) .eq. &
       this%DispersionData(ndis)%stringid(fnb2:lnb2) ) then
        write(outUnit,'(A,A,A,A,A)')&
        'Species ', trim(adjustl(this%Solutes(ns)%stringid)),&
        ' relation with dispersion ', trim(adjustl(dsid)), ' is consistent.'
        validRelation = .true.
        ! If relation is valid then assign dispersion pointer
        this%Solutes(ns)%dispersion => this%DispersionData(ndis)
        exit      
       end if 
      end do 

      if ( .not. validRelation ) then 
       write(outUnit,'(A,A,A,A)') 'Invalid dispersion ',&
               trim(adjustl(dsid)),' given to species ', trim(adjustl(this%Solutes(ns)%stringid))
       call ustop('Invalid relation between dispersion model and solute. Stop.')
      end if 

     end do
     
     ! Consistent relations, report 
     write(outUnit,'(A)') 'Relations between dispersion data and species/solutes are consistent.'

    end if
   
  end subroutine pr_ValidateDataRelations


  subroutine pr_SetSoluteDispersion( this, soluteId )
  !-----------------------------------------------------------------------------------
  !   - Assigns transportmodeldata dispersion pointers
  !     to dispersion parameters specified for the solute
  !   - Runs only if simulation is multidispersion and RWPT
  !-----------------------------------------------------------------------------------
  ! Specifications:
  !-----------------------------------------------------------------------------------
  implicit none
  class(TransportModelDataType), target :: this
  integer, intent(in) :: soluteId
  !-----------------------------------------------------------------------------------

    ! Temp: only isotropic model 
    select case(this%Solutes(soluteId)%Dispersion%modelKind)
    case(0)
      this%AlphaL => this%Solutes(soluteId)%Dispersion%AlphaL
      this%AlphaT => this%Solutes(soluteId)%Dispersion%AlphaTH
      this%DMEff  => this%Solutes(soluteId)%Dispersion%DMEff
    case(1)
      this%AlphaL  => this%Solutes(soluteId)%Dispersion%AlphaL
      this%AlphaLV => this%Solutes(soluteId)%Dispersion%AlphaLV
      this%AlphaTH => this%Solutes(soluteId)%Dispersion%AlphaTH
      this%AlphaTV => this%Solutes(soluteId)%Dispersion%AlphaTV
      this%DMEff   => this%Solutes(soluteId)%Dispersion%DMEff
    end select
    this%uniformDispersionParameters = .false.
    if ( this%Solutes(soluteId)%Dispersion%uniformParameters ) then  
      this%uniformDispersionParameters = .true.
    end if 

    ! TO BE DEPRECATED
    this%DMol = this%Solutes(soluteId)%Dispersion%dmaqueous 

  end subroutine pr_SetSoluteDispersion


  subroutine pr_LoadTimeStep(this, stressPeriod, timeStep)
  !-----------------------------------------------------------------------------------
  !   - Huge simplification from flowModelData%LoadTimeStep
  !   - It could be employed for loading time variable data for dispersion
  !   - In the meantime, only update time references and ICBOUND
  !   - It considers execution after flowModelData%LoadTimeStep  
  !-----------------------------------------------------------------------------------
  ! Specifications:
  !-----------------------------------------------------------------------------------
  implicit none
  ! input
  class(TransportModelDataType), target :: this
  integer,intent(in) :: stressPeriod, timeStep
  integer :: n 
  !-----------------------------------------------------------------------------------

    ! If following iboundts, update icbound
    ! for inactive flow-model cells 
    if ( this%followIBoundTS ) then
      ! Restart values to open
      this%ICBound(:) = 0
      ! Assign to rebound if flow-model inactive
      do n=1,this%Grid%CellCount
        if ( this%flowModelData%IBoundTS(n) .lt. 1 ) this%ICBound(n) = 1
      end do 
      this%ICBoundTS => null()
      this%ICBoundTS => this%ICBound
    end if

    this%CurrentStressPeriod = stressPeriod
    this%CurrentTimeStep = timeStep


  end subroutine pr_LoadTimeStep

    
end module TransportModelDataModule
