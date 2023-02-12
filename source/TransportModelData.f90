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
      doubleprecision,dimension(:),pointer :: AlphaL  => null()
      doubleprecision,dimension(:),pointer :: AlphaTH => null()
      doubleprecision,dimension(:),pointer :: AlphaTV => null()
      doubleprecision,dimension(:),pointer :: DMEff   => null()


      ! Local ICBound
      integer,dimension(:), pointer            :: ICBound
      !integer,allocatable,dimension(:)         :: ICBound
      integer,allocatable,dimension(:)         :: ICBoundTS

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

      ! Flow model data
      !type( FlowModelDataType ), pointer :: flowModelData => null()

      ! Budget reader
      !type(BudgetReaderType), pointer :: budgetReader => null()
      !this%BudgetReader => budgetReader

      ! TO BE DEPRECATED
      doubleprecision :: DMol

    contains

      procedure :: Initialize=>pr_Initialize
      procedure :: Reset=>pr_Reset
      procedure :: ReadSPCData=>pr_ReadSPCData
      procedure :: ReadDSPData=>pr_ReadDSPData
      procedure :: ValidateDataRelations => pr_ValidateDataRelations
      procedure :: LoadTimeStep=>pr_LoadTimeStep
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
      allocate(this%ICBoundTS(grid%CellCount))

      ! OLD
      !if(allocated(this%ICBound)) deallocate(this%ICBound)
      !allocate(this%ICBound(grid%CellCount))
      this%ICBound => simulationData%ICBound

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

      if(allocated(this%ICBoundTS)) deallocate( this%ICBoundTS )
      !if(allocated(this%ICBound)) deallocate( this%ICBound )
      this%ICBound => null()


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

        ! A single virtual solute storing the same 
        ! transport parameters for all pgroups
        this%Solutes(1)%id = 1
        this%Solutes(1)%stringid = 'SPC1'

        ! All pgroups to base solute 
        if ( allocated(this%Solutes(1)%pGroups) ) deallocate(this%Solutes(1)%pGroups)
        allocate(this%Solutes(1)%pGroups(simulationData%ParticleGroupCount))
        do n=1,simulationData%ParticleGroupCount
          this%Solutes(1)%pGroups(n) = n 
        end do 
        this%Solutes(1)%nParticleGroups = simulationData%ParticleGroupCount

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
        this%nSolutes = 1
        if ( allocated(this%Solutes) ) deallocate( this%Solutes )
        allocate( this%Solutes(this%nSolutes) )

        ! A single virtual solute storing the same 
        ! transport parameters for all pgroups
        this%Solutes(1)%id = 1
        this%Solutes(1)%stringid = 'SPC1'

        ! All pgroups to base solute 
        if ( allocated(this%Solutes(1)%pGroups) ) deallocate(this%Solutes(1)%pGroups)
        allocate(this%Solutes(1)%pGroups(simulationData%ParticleGroupCount))
        do n=1,simulationData%ParticleGroupCount
          this%Solutes(1)%pGroups(n) = n 
        end do 
        this%Solutes(1)%nParticleGroups = simulationData%ParticleGroupCount
        

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
  ! Specifications:
  !   - Assigns transportmodeldata dispersion pointers
  !     to dispersion parameters specified for the solute
  !   - Runs only if simulation is multidispersion and RWPT
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

    
end module TransportModelDataModule
