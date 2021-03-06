#if defined( DEVEL )
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: tendencies_mod.F90
!
! !DESCRIPTION: Module tendencies\_mod.F90 is a module to define, archive and
! write tracer tendencies. This module is still under development and only
! active if the DEVEL compiler flag is enabled. Also, the tracers to be
! diagnosed as well as the processes for which tendencies shall be calculated
! are currently hardcoded below. If enabled, this module will calculate 
! concentration tendencies for all defined species and processes and write 
! these into netCDF diagnostics. All tendencies are given as v/v/s.
!\\
!\\
! Tracer tendencies can be archived for as many processes are desired. For 
! each tendency process, an own instance of a tendency class has to be 
! defined (subroutine Tend\_CreateClass) and all active tracers for this class
! need be defined via subroutine Tend\_Add. Once a tendencies class is defined,
! the entry and exit 'checkpoint' of this tendency need be manually set in the
! code (subroutines Tend\_Stage1 and Tend\_Stage2).
!\\
!\\
! For example, suppose there is routine PROCESS in module example_mod.F90 and 
! we are interested in the tracer tendencies of O3 and CO by this process. We
! can then define a new tendency class (named 'PROCESS') during initialization
! of the tendencies (i.e. in tend\_init):
!
!    ! Create new class
!    CALL Tend_CreateClass( am_I_Root, Input_Opt, 'PROCESS', RC )
!    IF ( RC /= GIGC_SUCCESS ) RETURN
!
! The second step is to assign the tracers of interest to this tendency class:
!
!    ! Add species to classes
!    CALL Tend_Add ( am_I_Root, Input_Opt, 'PROCESS', IDTO3, RC )
!    IF ( RC /= GIGC_SUCCESS ) RETURN
!    CALL Tend_Add ( am_I_Root, Input_Opt, 'PROCESS', IDTCO, RC )
!    IF ( RC /= GIGC_SUCCESS ) RETURN
!
! The last step then involves the definition of the entry and exit points of
! the tendencies, e.g. the interval in between the tendencies shall be 
! calculated. To do so, we can wrap the Tend\_Stage routines around the 
! process of interest, e.g. in module example_mod.F90:
!
! CALL Tend\_Stage1 ( ... TendName='PROCESS', ... )
! CALL PROCESS ( ... )
! CALL Tend\_Stage2 ( ..., TendName='PROCESS', ... ) 
!\\
!\\
! The following six tendency classes are implemented by default: ADV (transport),
! CONV (convection), CHEM (chemistry), WETD (wet deposition), PBLMIX (PBL mixing, 
! includes emissions and dry deposition below PBL if non-local PBL is enabled),
! FLUX (emissions and dry depositions not coverd in PBLMIX).
! Subroutine Tend\_Init contains some example tendencies that are calculated 
! if flag 'DoTend' (subroutine Tend\_Init) is enabled. 
!
! !INTERFACE:
!
MODULE Tendencies_Mod 
!
! !USES:
!
  USE Precision_Mod
  USE HCO_Error_Mod
  USE HCO_Diagn_Mod
  USE GIGC_ErrCode_Mod
  USE Error_Mod,          ONLY : Error_Stop
  USE GIGC_Input_Opt_Mod, ONLY : OptInput
  USE GIGC_State_Met_Mod, ONLY : MetState
  USE GIGC_State_Chm_Mod, ONLY : ChmState

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: Tend_Init
  PUBLIC :: Tend_CreateClass
  PUBLIC :: Tend_Add
  PUBLIC :: Tend_Stage1
  PUBLIC :: Tend_Stage2
  PUBLIC :: Tend_Get
  PUBLIC :: Tend_Cleanup
!
! !PRIVATE MEMBER FUNCTIONS:
!
  PRIVATE :: Tend_FindClass
!
! !MODULE VARIABLES
!
  ! maximum string length of tendency name
  INTEGER, PARAMETER   :: MAXSTR = 31 

  ! Number of GEOS-Chem tracers
  INTEGER              :: nSpc = 0

  ! vector of tendency arrays
  TYPE :: TendArr
     REAL(f4),       POINTER   :: Arr(:,:,:) => NULL()
  END TYPE TendArr

  ! type holding tendencies of one type (class). Will
  ! be linked together via linked list.
  TYPE :: TendClass
     CHARACTER(LEN=MAXSTR)     :: TendName
     INTEGER                   :: Stage
     INTEGER,         POINTER  :: SpcInUse(:) => NULL()
     TYPE(TendArr),   POINTER  :: Tendency(:) => NULL()
     TYPE(TendClass), POINTER  :: NextTend    => NULL()
  END TYPE TendClass

  ! Tendency class linked list
  TYPE(TendClass),   POINTER   :: TendList    => NULL()
!
! !REVISION HISTORY:
!  14 Jul 2015 - C. Keller   - Initial version. 
!  26 Oct 2015 - C. Keller   - Now organize in linked list for more flexibility.
!EOP
!------------------------------------------------------------------------------
!BOC
CONTAINS
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Tend_Init
!
! !DESCRIPTION: Subroutine Tend\_Init is the wrapper routine to initialize the
! tendencies. At the moment, all tendencies are hardcoded and tendencies will 
! only be written if the manual flag `DoTend` is enabled. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Tend_Init ( am_I_Root, Input_Opt, State_Met, State_Chm, RC ) 
!
! !USES:
!
    USE TRACERID_MOD,  ONLY : IDTO3, IDTCO
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   ) :: am_I_Root  ! Are we on the root CPU?
    TYPE(OptInput),   INTENT(IN   ) :: Input_Opt  ! Input opts
    TYPE(MetState),   INTENT(IN   ) :: State_Met  ! met. state 
    TYPE(ChmState),   INTENT(IN   ) :: State_Chm  ! chm. state 
!
! !OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(OUT)   :: RC         ! Failure or success
!
! !REVISION HISTORY: 
!  26 Oct 2015 - C. Keller   - Initial version 
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    CHARACTER(LEN=255)       :: MSG
    CHARACTER(LEN=255)       :: LOC = 'Tend_Init (tendencies_mod.F)' 
   
    ! Set this to .TRUE. to enable some test diagnostics
    LOGICAL, PARAMETER       :: DoTend = .TRUE.
 
    !=======================================================================
    ! Tend_Init begins here!
    !=======================================================================

    ! Assume successful return
    RC = GIGC_SUCCESS

    ! Execute only if DoTend is enabled
    IF ( DoTend ) THEN

    ! Define classes
    CALL Tend_CreateClass( am_I_Root, Input_Opt, 'CHEM', RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    CALL Tend_CreateClass( am_I_Root, Input_Opt, 'ADV' , RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    CALL Tend_CreateClass( am_I_Root, Input_Opt, 'CONV', RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    CALL Tend_CreateClass( am_I_Root, Input_Opt, 'WETD', RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    CALL Tend_CreateClass( am_I_Root, Input_Opt, 'FLUX', RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    CALL Tend_CreateClass( am_I_Root, Input_Opt, 'PBLMIX', RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN

    ! Add species to classes
    CALL Tend_Add ( am_I_Root, Input_Opt, 'CHEM', IDTO3, RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    CALL Tend_Add ( am_I_Root, Input_Opt, 'CHEM', IDTCO, RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    CALL Tend_Add ( am_I_Root, Input_Opt, 'CONV', IDTCO, RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    CALL Tend_Add ( am_I_Root, Input_Opt, 'DRYD', IDTO3, RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    CALL Tend_Add ( am_I_Root, Input_Opt, 'ADV', IDTO3, RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN

    ENDIF ! test toggle

  END SUBROUTINE Tend_Init
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Tend_CreateClass 
!
! !DESCRIPTION: Subroutine Tend\_CreateClass creates a new tendency class. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Tend_CreateClass ( am_I_Root, Input_Opt, TendName, RC ) 
!
! !USES:
!
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   ) :: am_I_Root  ! Are we on the root CPU?
    TYPE(OptInput),   INTENT(IN   ) :: Input_Opt  ! Input opts
    CHARACTER(LEN=*), INTENT(IN   ) :: TendName   ! tendency class name
!
! !OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(OUT)   :: RC         ! Failure or success
!
! !REVISION HISTORY: 
!  26 Oct 2015 - C. Keller   - Initial version 
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    TYPE(TendClass), POINTER :: NewTend => NULL()
    LOGICAL                  :: FOUND
    CHARACTER(LEN=255)       :: MSG
    CHARACTER(LEN=255)       :: LOC = 'Tend_CreateClass (tendencies_mod.F)' 
    
    !=======================================================================
    ! Tend_CreateClass begins here!
    !=======================================================================

    ! Assume successful return
    RC = GIGC_SUCCESS

    ! Check if class already exists
    CALL Tend_FindClass( am_I_Root, TendName, FOUND, RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN

    IF ( .NOT. FOUND ) THEN 

       ! Eventually set local # of tracers variable
       IF ( nSpc <= 0 ) THEN
          nSpc = Input_Opt%N_TRACERS
       ENDIF

       ! Initialize class
       ALLOCATE(NewTend)

       ! Set tendency class name
       NewTend%TendName = TRIM(TendName)

       ! Initialize stage
       newTend%Stage    = -1

       ! Initialize vector with species flags
       ALLOCATE(NewTend%SpcInUse(nSpc))
       NewTend%SpcInUse(:) = 0

       ! Initialize tendency arrays (only allocated when needed)
       ALLOCATE(NewTend%Tendency(nSpc))

       ! Add tendency class to linked list
       NewTend%NextTend => TendList
       TendList         => NewTend 
    ENDIF

  END SUBROUTINE Tend_CreateClass
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Tend_FindClass
!
! !DESCRIPTION: Subroutine Tend\_FindClass searches for a tendency class. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Tend_FindClass ( am_I_Root, TendName, FOUND, RC, ThisTend ) 
!
! !USES:
!
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   )          :: am_I_Root  ! Are we on the root CPU?
    CHARACTER(LEN=*), INTENT(IN   )          :: TendName   ! tendency class name
!
! !OUTPUT PARAMETERS:
!
    LOGICAL,          INTENT(  OUT)          :: FOUND      ! class found
    INTEGER,          INTENT(  OUT)          :: RC         ! Failure or success
    TYPE(TendClass),  POINTER,      OPTIONAL :: ThisTend   ! Pointer to this class
!
! !REVISION HISTORY: 
!  26 Oct 2015 - C. Keller   - Initial version 
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    TYPE(TendClass), POINTER :: TmpTend => NULL()
    
    !=======================================================================
    ! Tend_FindClass begins here!
    !=======================================================================

    ! Assume successful return
    RC = GIGC_SUCCESS

    ! Init
    FOUND = .FALSE.

    ! Loop through linked list and search for class with same name
    TmpTend => TendList
    DO WHILE ( ASSOCIATED(TmpTend) ) 
    
       ! Is this the tendency of interest?  
       IF ( TRIM(TmpTend%TendName) == TRIM(TendName) ) THEN
          FOUND = .TRUE.
          EXIT
       ENDIF
 
       ! Advance in list
       TmpTend => TmpTend%NextTend
    END DO

    ! Eventually 
    IF ( PRESENT(ThisTend) ) ThisTend => TmpTend

    ! Cleanup
    TmpTend => NULL()

  END SUBROUTINE Tend_FindClass
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Tend_Cleanup
!
! !DESCRIPTION: Subroutine Tend\_Cleanup cleans up the tendencies linked list. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Tend_Cleanup ( ) 
!
! !USES:
!
!
! !INPUT PARAMETERS:
!
!
! !REVISION HISTORY: 
!  26 Oct 2015 - C. Keller   - Initial version 
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER                  :: I
    TYPE(TendClass), POINTER :: ThisTend => NULL()
    TYPE(TendClass), POINTER :: NextTend => NULL()
    
    !=======================================================================
    ! Tend_Cleanup begins here!
    !=======================================================================

    ! Loop through linked list and search for class with same name
    ThisTend => TendList
    DO WHILE ( ASSOCIATED(ThisTend) ) 
   
       ! Get pointer to next tendency
       NextTend => ThisTend%NextTend

       ! Cleanup every array
       DO I = 1, nSpc
          IF ( ASSOCIATED(ThisTend%Tendency(I)%Arr) ) THEN
             DEALLOCATE( ThisTend%Tendency(I)%Arr)
          ENDIF
       ENDDO
       DEALLOCATE(ThisTend%Tendency)
       DEALLOCATE(ThisTend%SpcInUse)

       ! Cleanup
       ThisTend%NextTend => NULL()
       NULLIFY(ThisTend)
 
       ! Advance in list
       ThisTend => NextTend
    END DO

    ! Cleanup
    ThisTend => NULL()
    NextTend => NULL()
    nSpc = 0

  END SUBROUTINE Tend_Cleanup
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Tend_Add
!
! !DESCRIPTION: Subroutine Tend\_Add adds a species to a tendency class. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Tend_Add ( am_I_Root, Input_Opt, TendName, TrcID, RC, CreateClass )
!
! !USES:
!
    USE CMN_SIZE_MOD,       ONLY : IIPAR, JJPAR, LLPAR
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   )           :: am_I_Root   ! Are we on the root CPU?
    TYPE(OptInput),   INTENT(IN   )           :: Input_Opt   ! Input opts
    CHARACTER(LEN=*), INTENT(IN   )           :: TendName    ! Tendency class name 
    INTEGER,          INTENT(IN   )           :: TrcID       ! Tracer ID 
    LOGICAL,          INTENT(IN   ), OPTIONAL :: CreateClass ! Create class if missing?
!
! !OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(  OUT)           :: RC          ! Failure or success
!
! !REVISION HISTORY: 
!  14 Jul 2015 - C. Keller   - Initial version 
!  26 Oct 2015 - C. Keller   - Update for linked list 
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER                  :: AS
    INTEGER                  :: Collection
    LOGICAL                  :: FOUND 
    TYPE(TendClass), POINTER :: ThisTend => NULL()
    CHARACTER(LEN=63)        :: DiagnName
    CHARACTER(LEN=255)       :: MSG
    CHARACTER(LEN=255)       :: LOC = 'Tend_Add (tendencies_mod.F)' 
    
    !=======================================================================
    ! Tend_Add begins here!
    !=======================================================================

    ! Assume successful return
    RC = GIGC_SUCCESS

    ! Ignore invalid tracer IDs
    IF ( TrcID <= 0 ) RETURN

    ! Search for diagnostics class
    CALL Tend_FindClass( am_I_Root, TendName, FOUND, RC, ThisTend=ThisTend ) 

    ! Eventually create this class if it does not exist yet
    IF ( .NOT. FOUND .AND. PRESENT( CreateClass ) ) THEN
       IF ( CreateClass ) THEN
          ! Create class
          CALL Tend_CreateClass( am_I_Root, Input_Opt, TendName, RC )
          IF ( RC /= GIGC_SUCCESS ) RETURN
          ! Get pointer to class object 
          CALL Tend_FindClass( am_I_Root, TendName, FOUND, RC, ThisTend=ThisTend ) 
       ENDIF
    ENDIF

    ! Return here if class not found
    IF ( .NOT. FOUND ) RETURN

    ! Tracer ID must not exceed # of tendency tracers
    IF ( TrcID > nSpc ) THEN
       WRITE(MSG,*) 'Tracer ID exceeds number of tendency tracers: ', TrcID, ' > ', nSpc
       CALL ERROR_STOP( MSG, LOC ) 
       RETURN
    ENDIF 
 
    ! Name of diagnostics 
    DiagnName = 'TEND_' // TRIM(TendName) // '_' //   &
                TRIM( Input_Opt%TRACER_NAME( TrcID ) )

    ! Mark tracer as being used
    ThisTend%SpcInUse(TrcID) = 1

    ! Make sure array is allocated
    IF ( .NOT. ASSOCIATED(ThisTend%Tendency(TrcID)%Arr) ) THEN
       ALLOCATE(ThisTend%Tendency(TrcID)%Arr(IIPAR,JJPAR,LLPAR),STAT=AS)
       IF ( AS /= 0 ) THEN
          MSG = 'Tendency allocation error: ' // TRIM(DiagnName)
          CALL ERROR_STOP( MSG, LOC ) 
          RETURN
       ENDIF
    ENDIF

    ! Get diagnostic parameters from the Input_Opt object
    Collection = Input_Opt%DIAG_COLLECTION

    ! Create container for tendency
    CALL Diagn_Create( am_I_Root,                     &
                       Col       = Collection,        & 
!                       cID       = cID,               &
                       cName     = TRIM( DiagnName ), &
                       AutoFill  = 0,                 &
                       ExtNr     = -1,                &
                       Cat       = -1,                &
                       Hier      = -1,                &
                       HcoID     = -1,                &
                       SpaceDim  =  3,                &
                       OutUnit   = 'v/v/s',           &
                       OutOper   = 'Mean',            &
                       OkIfExist = .TRUE.,            &
                       RC        = RC )
   
    IF ( RC /= HCO_SUCCESS ) THEN
       MSG = 'Cannot create diagnostics: ' // TRIM(DiagnName)
       CALL ERROR_STOP( MSG, LOC ) 
    ENDIF

  END SUBROUTINE Tend_Add
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Tend_Stage1
!
! !DESCRIPTION: Subroutine Tend\_Stage1 archives the current tracer 
! concentrations into the local tendency arrays.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Tend_Stage1( am_I_Root, Input_Opt, State_Met, &
                          State_Chm, TendName,  IsInvv,    RC ) 
!
! !USES:
!
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   ) :: am_I_Root  ! Are we on the root CPU?
    TYPE(OptInput),   INTENT(IN   ) :: Input_Opt  ! Input opts
    TYPE(MetState),   INTENT(IN   ) :: State_Met  ! Met state
    TYPE(ChmState),   INTENT(IN   ) :: State_Chm  ! Chemistry state 
    CHARACTER(LEN=*), INTENT(IN   ) :: TendName   ! tendency name 
    LOGICAL,          INTENT(IN   ) :: IsInvv     ! Is STT in v/v? 
!
! !OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(OUT)   :: RC         ! Failure or success
!
! !REVISION HISTORY: 
!  14 Jul 2015 - C. Keller   - Initial version 
!  26 Oct 2015 - C. Keller   - Update to include tendency classes
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER                  :: I
    LOGICAL                  :: FOUND 
    REAL(f4), POINTER        :: Ptr3D(:,:,:) => NULL()
    TYPE(TendClass), POINTER :: ThisTend => NULL()
    CHARACTER(LEN=255)       :: MSG
    CHARACTER(LEN=255)       :: LOC = 'TEND_STAGE1 (tendencies_mod.F)' 
    
    !=======================================================================
    ! TEND_STAGE1 begins here!
    !=======================================================================

    ! Assume successful return
    RC = GIGC_SUCCESS

    ! Find tendency class
    CALL Tend_FindClass( am_I_Root, TendName, FOUND, RC, ThisTend=ThisTend )
    IF ( .NOT. FOUND .OR. .NOT. ASSOCIATED(ThisTend) ) RETURN

    ! Loop over # of tendencies species
    DO I = 1, nSpc 

       ! Skip if tracer is not in use    
       IF ( ThisTend%SpcInUse(I) <= 0 ) CYCLE

       ! Get pointer to 3D array to be filled 
       Ptr3D => ThisTend%Tendency(I)%Arr

       ! Fill 3D array with current values. Make sure it's in v/v
       IF ( IsInvv ) THEN
          Ptr3D = State_Chm%Tracers(:,:,:,I)
       ELSE
          Ptr3D = State_Chm%Tracers(:,:,:,I) &
                * Input_Opt%TCVV(I) / State_Met%AD(:,:,:)
       ENDIF

       ! Cleanup
       Ptr3D => NULL()
    ENDDO !I

    ! Update stage 
    ThisTend%Stage = 1

    ! Cleanup
    ThisTend => NULL()

  END SUBROUTINE Tend_Stage1
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Tend_Stage2
!
! !DESCRIPTION: Subroutine Tend_Stage2 calculates the tendencies and 
! writes them into the diagnostics arrays. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Tend_Stage2( am_I_Root, Input_Opt, State_Met, &
                          State_Chm, TendName,  IsInvv, DT, RC ) 
!
! !USES:
!
    USE CMN_SIZE_MOD,      ONLY : IIPAR, JJPAR, LLPAR
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   ) :: am_I_Root  ! Are we on the root CPU?
    TYPE(OptInput),   INTENT(IN   ) :: Input_Opt  ! Input opts
    TYPE(MetState),   INTENT(IN   ) :: State_Met  ! Met state
    TYPE(ChmState),   INTENT(IN   ) :: State_Chm  ! Chemistry state 
    CHARACTER(LEN=*), INTENT(IN   ) :: TendName   ! tendency name 
    LOGICAL,          INTENT(IN   ) :: IsInvv     ! Is tracer in v/v? 
    REAL(fp),         INTENT(IN   ) :: DT         ! delta time, in seconds 
!
! !OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(OUT)   :: RC         ! Failure or success
!
! !REVISION HISTORY: 
!  14 Jul 2015 - C. Keller   - Initial version 
!  26 Oct 2015 - C. Keller   - Update to include tendency classes
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    LOGICAL                  :: ZeroTend
    LOGICAL                  :: FOUND
    INTEGER                  :: cID, I
    REAL(f4), POINTER        :: Ptr3D(:,:,:) => NULL()
    REAL(f4)                 :: TEND(IIPAR,JJPAR,LLPAR)
    TYPE(TendClass), POINTER :: ThisTend => NULL()
    CHARACTER(LEN=63)        :: DiagnName
    CHARACTER(LEN=255)       :: MSG
    CHARACTER(LEN=255)       :: LOC = 'TEND_STAGE2 (tendencies_mod.F)' 
   
    !=======================================================================
    ! TEND_STAGE2 begins here!
    !=======================================================================

    ! Assume successful return
    RC = GIGC_SUCCESS

    ! Find tendency class
    CALL Tend_FindClass( am_I_Root, TendName, FOUND, RC, ThisTend=ThisTend )
    IF ( .NOT. FOUND .OR. .NOT. ASSOCIATED(ThisTend) ) RETURN

    ! Error check: stage 2 must be called after stage 1
    ZeroTend = .FALSE.
    IF ( ThisTend%Stage /= 1 ) THEN
       IF ( am_I_Root ) THEN
          WRITE(*,*) 'Warning: cannot do tendency stage 2 - stage 1 not yet called: ', TRIM(TendName) 
       ENDIF
       ZeroTend = .TRUE.
    ENDIF

    ! Loop over # of tendencies species
    DO I = 1, nSpc 

       ! Skip if not used    
       IF ( ThisTend%SpcInUse(I) <= 0 ) CYCLE

       ! Name of diagnostics 
       DiagnName = 'TEND_' // TRIM(ThisTend%TendName) // '_' //   &
                   TRIM( Input_Opt%TRACER_NAME(I) )

       ! Get pointer to 3D array, define time interval
       Ptr3D => ThisTend%Tendency(I)%Arr 

       ! Calculate tendency in v/v/s
       IF ( ZeroTend ) THEN
          Tend = 0.0_f4
       ELSE
          IF ( IsInvv ) THEN
             Tend = ( State_Chm%Tracers(:,:,:,I) - Ptr3D(:,:,:) ) / DT
          ELSE
             Tend = ( ( State_Chm%Tracers(:,:,:,I)                &
                      * Input_Opt%TCVV(I) / State_Met%AD(:,:,:) ) &
                      - Ptr3D(:,:,:) ) / DT
          ENDIF
       ENDIF

       ! Update diagnostics array
       CALL Diagn_Update( am_I_Root, cName=DiagnName, Array3D=Tend, &
                          COL=Input_Opt%DIAG_COLLECTION, RC=RC       )
       IF ( RC /= HCO_SUCCESS ) THEN 
          WRITE(MSG,*) 'Error in updating diagnostics with ID ', cID
          CALL ERROR_STOP ( MSG, LOC )
          RC = GIGC_FAILURE
          RETURN
       ENDIF

       ! Update values 
       Ptr3D = Tend 
       Ptr3D => NULL()

    ENDDO !I

    ! Update stage 
    ThisTend%Stage = 2

    ! Cleanup
    ThisTend => NULL()

  END SUBROUTINE Tend_Stage2
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Tend_Get
!
! !DESCRIPTION: Subroutine Tend_Get returns the current tendency for the 
! given tracer and tendency class. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Tend_Get( am_I_Root, Input_Opt, TendName, TrcID, Stage, Tend, RC ) 
!
! !USES:
!
    USE CMN_SIZE_MOD,      ONLY : IIPAR, JJPAR, LLPAR
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   ) :: am_I_Root   ! Are we on the root CPU?
    TYPE(OptInput),   INTENT(IN   ) :: Input_Opt   ! Input opts
    CHARACTER(LEN=*), INTENT(IN   ) :: TendName    ! tendency name 
    INTEGER,          INTENT(IN   ) :: TrcID       ! Tracer ID 
!
! !OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(  OUT) :: Stage       ! Stage of tendency:
                                                   ! 0=does not exist; 1=stage 1; 2=stage 2 
    REAL(f4),         POINTER       :: Tend(:,:,:) ! Tendency array
    INTEGER,          INTENT(OUT)   :: RC          ! Failure or success
!
! !REVISION HISTORY: 
!  14 Jul 2015 - C. Keller   - Initial version 
!  26 Oct 2015 - C. Keller   - Update to include tendency classes
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    LOGICAL                  :: FOUND
    INTEGER                  :: cID, I
    TYPE(TendClass), POINTER :: ThisTend => NULL()
    CHARACTER(LEN=63)        :: DiagnName
    CHARACTER(LEN=255)       :: MSG
    CHARACTER(LEN=255)       :: LOC = 'TEND_GET (tendencies_mod.F)' 
   
    !=======================================================================
    ! TEND_GET begins here!
    !=======================================================================

    ! Assume successful return
    RC = GIGC_SUCCESS

    ! Init
    Tend  => NULL()
    Stage =  0

    ! Find tendency class
    CALL Tend_FindClass( am_I_Root, TendName, FOUND, RC, ThisTend=ThisTend )
    IF ( .NOT. FOUND .OR. .NOT. ASSOCIATED(ThisTend) ) RETURN

    ! Skip if not used    
    IF ( ThisTend%SpcInUse(TrcID) <= 0 ) RETURN 

    ! Get pointer to tendency
    Tend  => ThisTend%Tendency(TrcID)%Arr
    Stage =  ThisTend%Stage

    ! Cleanup
    ThisTend => NULL()

  END SUBROUTINE Tend_Get
!EOC
END MODULE Tendencies_Mod
#endif
