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
      doubleprecision,dimension(:),pointer :: AlphaLong => null()
      doubleprecision,dimension(:),pointer :: AlphaTran => null()
      doubleprecision,dimension(:),pointer :: AlphaL    => null()
      doubleprecision,dimension(:),pointer :: AlphaTH   => null()
      doubleprecision,dimension(:),pointer :: AlphaTV   => null()
      doubleprecision,dimension(:),pointer :: DMEff     => null()

      ! ICBound
      integer,allocatable, dimension(:) :: ICBound
      integer, pointer, dimension(:)    :: ICBoundTS
      !integer,allocatable, dimension(:) :: ICBoundTS
      integer                           :: defaultICBound ! applied to model bounadaries
      logical                           :: followIBoundTS
      logical                           :: followIBound

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
     
      this%Grid => null()
      this%DMol = 0
      !if(allocated(this%AlphaLong)) deallocate( this%AlphaLong )
      !if(allocated(this%AlphaTran)) deallocate( this%AlphaTran )

      this%AlphaLong  => null()
      this%AlphaTran => null()
      this%DMEff => null()

      !if(allocated(this%ICBoundTS)) deallocate( this%ICBoundTS )
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
    integer            :: icol,istart,istop,n,nd,currentDim
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

    this%followIBound = .false.
    this%followIBoundTS = .false.
    select case(impFormat)
    case(0)
     ! Copy IBound into ICBound and assign pointer, does not 
     ! change in time
     this%followIBound = .true.
     this%ICBound(:) = this%flowModelData%IBound(:)
     this%ICBoundTS => this%ICBound
     write(outUnit,'(A)') 'Impermeable cells follow flowModelData%IBound.'
    case(1)
     ! Activate flag, ICBoundTS is updated every time step
     this%followIBoundTS = .true.
     write(outUnit,'(A)') 'Impermeable cells follow flowModelData%IBoundTS, includes dry cells.'
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

    ! Allocate solutes
    if ( allocated(this%Solutes) ) deallocate( this%Solutes )
    allocate( this%Solutes(this%nSolutes) )

    ! Initialize these guys
    do ns =1,this%nSolutes
      write(outUnit,'(A,I4)') 'Processing species specification ', ns

      ! Internal id assigned increasingly
      this%Solutes(ns)%id = n 
      this%Solutes(ns)%userid = n 
      
      ! Read the string id
      read(spcUnit, '(a)') line
      icol = 1
      call urword(line, icol, istart, istop, 0, n, r, 0, 0)
      this%Solutes(ns)%stringid = line(istart:istop)

      ! Assign pgroups related to the solute
      if ( simulationData%ParticlesMassOption .ne. 2 ) then
        ! Read pgroups related to the solute from the pkg
        write(outUnit,'(A)') 'ParticleMassOption.ne.2: '
        write(outUnit,'(A)') 'Read pgroups related to this species/solutes from the pkg.'

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
        write(outUnit,'(A)') 'ParticleMassOption.eq.2: '
        write(outUnit,'(A)') 'Interpret pgroups related to this species/solutes from pgroups list.'

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

       write(outUnit,'(A,A,A)')&
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
    type(ModpathSimulationDataType), pointer :: simulationData
    class(ModflowRectangularGridType),pointer :: grid
    type(DispersionDataType),pointer          :: disp
    integer :: isThisFileOpen 
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

      ! Dispersion id assigned increasingly
      disp%id = ndis 
      
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

      ! Validate dispersion model 
      select case ( disp%modelKind )
        ! Linear, isotropic
        case(1)
          continue
        case default
          write(outUnit,*) 'Invalid dispersion model. Given ', disp%modelKind
          call ustop('Invalid dispersion model. Stop')          
      end select

      ! Allocate
      call disp%InitializeByModelKind( grid%CellCount )

      ! Read dispersion data
      select case ( disp%modelKind )
        ! Linear, isotropic
        case(1)
          write(outUnit,'(A)') 'Dispersion model is linear, isotropic: will read dmeff and dispersivities'

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

          !! Read ALPHATV
          !if((grid%GridType .eq. 1) .or. (grid%GridType .eq. 3)) then
          !  call u3ddblmp(dspUnit, outUnit, grid%LayerCount, grid%RowCount,      &
          !    grid%ColumnCount, grid%CellCount, disp%AlphaTV, anamelin(4))                      
          !else if((grid%GridType .eq. 2) .or. (grid%GridType .eq. 4)) then
          !  call u3ddblmpusg(dspUnit, outUnit, grid%CellCount, grid%LayerCount,  &
          !    disp%AlphaTV, anamelin(4), cellsPerLayer)
          !else
          !  write(outUnit,*) 'Invalid grid type specified when reading ALPHATV array data.'
          !  write(outUnit,*) 'Stopping.'
          !  call ustop(' ')          
          !end if

        !! Nonlinear
        !case(2)
        !  write(outUnit,'(A)') 'Dispersion model is nonlinear'

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
    integer :: ns, ndis, did
    character(len=300) :: dsid
    logical :: validRelation 
    !--------------------------------------------------------------

    write(outUnit, *)
    write(outUnit, '(1x,a)') 'TransportModelData: verify SPC-DSP relations '
    write(outUnit, '(1x,a)') '---------------------------------------------'


    simulationData => this%simulationData

    if ( simulationData%SolutesOption .eq. 0 ) then
     ! If simulation is single dispersion 
     write(outUnit,'(A)') 'Simulation is single dispersion. Is consistent. '

     ! TEMPORARY
     this%Solutes(1)%dispersion => this%DispersionData(1)
     this%AlphaLong => this%Solutes(1)%Dispersion%AlphaL
     this%AlphaTran => this%Solutes(1)%Dispersion%AlphaTH
     this%DMEff     => this%Solutes(1)%Dispersion%DMEff


     this%DMol = this%Solutes(1)%Dispersion%dmaqueous ! TO BE DEPRECATED

     ! Is valid
     return

    else if ( simulationData%SolutesOption .eq. 1 ) then
     ! If simulation is multidispersion
     write(outUnit,'(A)') 'Simulation is multidispersion... '

     ! Should confirm that dispersion models 
     ! specified for the different solutes
     ! exist in dispersiondata array

     do ns = 1, this%nSolutes

      did  = this%Solutes(ns)%dispersionId
      dsid = this%Solutes(ns)%dispersionStringId

      validRelation = .false.
      ! Loop over dispersion models
      do ndis = 1, this%nDispersion
       if ( this%DispersionData(ndis)%id .eq. did ) then
        write(outUnit,'(A,I3,A,I3)') 'Solute ', ns, ' relation with dispersion ', did, ' is consistent.'
        validRelation = .true.
        ! If relation is valid then assign dispersion pointer
        this%Solutes(ns)%dispersion => this%DispersionData(ndis)
        exit      
       end if 
      end do 

      if ( .not. validRelation ) then 
       write(outUnit,'(A,I4,A,A)') 'Invalid dispersion id ',did,' for solute ', this%Solutes(ns)%stringid
       call ustop('Invalid relation between dispersion model and solute. Stop.')
      end if 

     end do
     
     ! Consistent relations, report 
     write(outUnit,'(A)') 'Relations between dispersion data and species/solutes are consistent.'

    end if
   
    
    ! Same thing for PGROUPS and SPC (?)


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

    ! TEMPORARY ( should understand al, ath, athv, and so on )
    this%AlphaLong => this%Solutes(soluteId)%Dispersion%AlphaL
    this%AlphaTran => this%Solutes(soluteId)%Dispersion%AlphaTH
    this%DMEff     => this%Solutes(soluteId)%Dispersion%DMEff

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
  class(TransportModelDataType) :: this
  integer,intent(in) :: stressPeriod, timeStep
  ! local
  integer :: n, cellCount
  !-----------------------------------------------------------------------------------


    if ( this%followIBoundTS ) then
      this%ICBoundTS => null()
      this%ICBoundTS => this%flowModelData%IBoundTS
    end if

    this%CurrentStressPeriod = stressPeriod
    this%CurrentTimeStep = timeStep


  end subroutine pr_LoadTimeStep

    
end module TransportModelDataModule
