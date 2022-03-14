! experiments with polymorphism and MPI in Fortran by Rand Huso
! Sending the object using OpenMPI requires referencing the reference component of the base class. See:
! http://www.mcs.anl.gov/research/projects/mpi/mpi-standard/mpi-report-2.0/node236.htm for more description.
! debug with: mpiexec -n 2 xterm -e gdb bin/PolymorphicTest
!
! All these funny lines are because my monitor is wider than 132 characters.. and it helps separate sections for finding the section
! I'm seeking:
!--------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|-|
! I just threw this module together to simplify constructing MPI Datatypes
MODULE MpiAssist
    USE mpi
    IMPLICIT NONE

    TYPE MpiAssistType
        INTEGER :: MA_elements
        INTEGER( kind=MPI_ADDRESS_KIND ), ALLOCATABLE :: MA_address(:)
        INTEGER, ALLOCATABLE :: MA_len(:)
        INTEGER, ALLOCATABLE :: MA_type(:)
        INTEGER( kind=MPI_ADDRESS_KIND ) :: MA_baseAddress
    CONTAINS
        PROCEDURE :: MA_loadComponent
        PROCEDURE :: MA_getType
        FINAL :: MA_destructor
    END TYPE

    INTERFACE MpiAssistType
        PROCEDURE MA_constructor
    END INTERFACE

CONTAINS

    FUNCTION MA_constructor( elements ) RESULT( self )
        TYPE( MpiAssistType ) :: self
        INTEGER, INTENT( in ) :: elements
        self%MA_elements = elements
        ALLOCATE( self%MA_address( self%MA_elements ))
        ALLOCATE( self%MA_len( self%MA_elements ))
        ALLOCATE( self%MA_type( self%MA_elements ))
    END FUNCTION

    SUBROUTINE MA_destructor( self )
        TYPE( MpiAssistType ) :: self
        IF( allocated( self%MA_address )) DEALLOCATE( self%MA_address )
        IF( allocated( self%MA_len )) DEALLOCATE( self%MA_len )
        IF( allocated( self%MA_type )) DEALLOCATE( self%MA_type )
    END SUBROUTINE

    SUBROUTINE MA_loadComponent( self, thisItem, itemNumber, itemCount, itemType )
        CLASS( MpiAssistType ), INTENT( inout ) :: self
        CLASS( * ), INTENT( in ) :: thisItem
        INTEGER, INTENT( in ) :: itemNumber
        INTEGER, INTENT( in ) :: itemCount
        INTEGER, INTENT( in ) :: itemType
        INTEGER :: iErr
        CALL MPI_GET_ADDRESS( thisItem, self%MA_address( itemNumber ), iErr )
        IF( 1 == itemNumber ) self%MA_baseAddress = self%MA_address( itemNumber )
        self%MA_address( itemNumber ) = self%MA_address( itemNumber ) -self%MA_baseAddress
        self%MA_len( itemNumber ) = itemCount
        self%MA_type( itemNumber ) = itemType
    END SUBROUTINE

    FUNCTION MA_getType( self )
        CLASS( MpiAssistType ), INTENT( inout ) :: self
        INTEGER :: MA_getType
        INTEGER :: iErr
        CALL MPI_TYPE_CREATE_STRUCT( self%MA_elements, self%MA_len, self%MA_address, self%MA_type, MA_getType, iErr )
        CALL MPI_TYPE_COMMIT( MA_getType, iErr )
    END FUNCTION
END MODULE

!--------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|-|
! Base class for this particular test.
! This lets me use polymorphism where needed when referencing objects, and provides the base reference item for sending with OpenMPI
MODULE MsgBase
    USE MpiAssist
    IMPLICIT NONE

    TYPE, ABSTRACT :: MsgBaseType
        INTEGER :: baseId = 59 ! referenced by OpenMPI (must be first item)
        INTEGER :: mpiDatatype = MPI_DATATYPE_NULL ! not intended to be sent over MPI
    CONTAINS
        PROCEDURE( getMpiDatatype ), DEFERRED :: getMpiDatatype
        PROCEDURE( getStr ), DEFERRED :: str
    END TYPE

    ABSTRACT INTERFACE
        FUNCTION getMpiDatatype( self ) RESULT( response )
            IMPORT :: MsgBaseType
            CLASS( MsgBaseType ) :: self
            INTEGER :: response
        END FUNCTION

        FUNCTION getStr( self ) RESULT( identityStr )
            IMPORT :: MsgBaseType
            CLASS( MsgBaseType ) :: self
            CHARACTER( len=60 ) :: identityStr
        END FUNCTION
    END INTERFACE
END MODULE

!--------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|-|
! The full implementation will also get the available memory, the recent figures for use, the number of cores, &c.
MODULE MsgCapabilities
    USE MsgBase
    USE mpi
    IMPLICIT NONE

    TYPE, EXTENDS( MsgBaseType ) :: MsgCapabilitiesType
        INTEGER :: mRank
        INTEGER :: mSize
        CHARACTER( len=MPI_MAX_PROCESSOR_NAME ) :: hostname
    CONTAINS
        PROCEDURE :: getMpiDatatype => MCT_getMpiDatatype
        PROCEDURE :: str => MCT_str
    END TYPE

    INTERFACE MsgCapabilitiesType
        PROCEDURE MCT_constructor
    END INTERFACE

CONTAINS

    !----|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|-|
    FUNCTION MCT_str( self ) RESULT( identityStr )
        CLASS( MsgCapabilitiesType ) :: self
        CHARACTER( len=60 ) :: identityStr
        WRITE( identityStr, '("[", I0, ":", I0, " id:", I0, " Datatype:", I0, " Name:", A, "]")') &
            self%mRank, self%mSize, self%baseId, self%mpiDatatype, trim(self%hostname)
    END FUNCTION

    !----|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|-|
    FUNCTION MCT_getMpiDatatype( self ) RESULT( response )
        CLASS( MsgCapabilitiesType ) :: self
        INTEGER :: response
        TYPE( MpiAssistType ) :: maType
        IF( MPI_DATATYPE_NULL == self%mpiDatatype ) THEN ! ensure that we only do this once on each node
           maType = MpiAssistType( 4 )
           CALL maType%MA_loadComponent( thisItem=self%baseId, itemNumber=1, itemCount=1, itemType=MPI_INTEGER )
           CALL maType%MA_loadComponent( thisItem=self%mRank, itemNumber=2, itemCount=1, itemType=MPI_INTEGER )
           CALL maType%MA_loadComponent( thisItem=self%mSize, itemNumber=3, itemCount=1, itemType=MPI_INTEGER )
           CALL maType%MA_loadComponent(thisItem=self%hostname,itemNumber=4,itemCount=MPI_MAX_PROCESSOR_NAME,itemType=MPI_CHARACTER)
           self%mpiDatatype = maType%MA_getType()
        END IF
        response = self%mpiDatatype
    END FUNCTION

    !----|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|-|
    FUNCTION MCT_constructor() RESULT( self )
        TYPE( MsgCapabilitiesType ) :: self
        INTEGER :: iStatus
        INTEGER :: iErr
        INTEGER :: hostnameLength
        CALL MPI_Get_processor_name( self%hostname, hostnameLength, iErr )
        CALL MPI_COMM_RANK( MPI_COMM_WORLD, self%mRank, iErr )
        CALL MPI_COMM_SIZE( MPI_COMM_WORLD, self%mSize, iErr )
    END FUNCTION
END MODULE

!--------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|-|
MODULE MsgOther
    USE MsgBase
    USE mpi
    IMPLICIT NONE

    INTEGER, PARAMETER :: MSG_LEN = 20

    TYPE, EXTENDS( MsgBaseType ) :: MsgOtherType
        INTEGER :: mRank
        INTEGER :: mSize
        CHARACTER( len=MPI_MAX_PROCESSOR_NAME ) :: hostname
        CHARACTER( len=MSG_LEN ) :: message
    CONTAINS
        PROCEDURE :: getMpiDatatype => MOT_getMpiDatatype
        PROCEDURE :: str => MOT_str
    END TYPE

    INTERFACE MsgOtherType
        PROCEDURE MOT_constructor
    END INTERFACE

CONTAINS

    !----|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|-|
    FUNCTION MOT_str( self ) RESULT( identityStr )
        CLASS( MsgOtherType ) :: self
        CHARACTER( len=60 ) :: identityStr
        WRITE( identityStr, '("[", "id:", I0, " Datatype:", I0, 1X, I0, ":", I0, " Name:", A, 1X, A, "]")') &
            self%baseId, self%mpiDatatype, self%mRank, self%mSize, trim(self%hostname), trim(self%message)
    END FUNCTION

    !----|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|-|
    FUNCTION MOT_getMpiDatatype( self ) RESULT( response )
        CLASS( MsgOtherType ) :: self
        INTEGER :: response
        TYPE( MpiAssistType ) :: maType
        IF( MPI_DATATYPE_NULL == self%mpiDatatype ) THEN ! ensures that we only do this once on each node
           maType = MpiAssistType( 5 )
           CALL maType%MA_loadComponent( thisItem=self%baseId, itemNumber=1, itemCount=1, itemType=MPI_INTEGER )
           CALL maType%MA_loadComponent( thisItem=self%mRank, itemNumber=2, itemCount=1, itemType=MPI_INTEGER )
           CALL maType%MA_loadComponent( thisItem=self%mSize, itemNumber=3, itemCount=1, itemType=MPI_INTEGER )
           CALL maType%MA_loadComponent(thisItem=self%hostname,itemNumber=4,itemCount=MPI_MAX_PROCESSOR_NAME,itemType=MPI_CHARACTER)
           CALL maType%MA_loadComponent( thisItem=self%message, itemNumber=5, itemCount=MSG_LEN, itemType=MPI_CHARACTER )
           self%mpiDatatype = maType%MA_getType()
        END IF
        response = self%mpiDatatype
    END FUNCTION

    !----|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|-|
    FUNCTION MOT_constructor( message ) RESULT( self )
        TYPE( MsgOtherType ) :: self
        INTEGER :: iStatus
        INTEGER :: iErr
        INTEGER :: hostnameLength
        CHARACTER( len=* ) :: message
        CALL MPI_Get_processor_name( self%hostname, hostnameLength, iErr )
        CALL MPI_COMM_RANK( MPI_COMM_WORLD, self%mRank, iErr )
        CALL MPI_COMM_SIZE( MPI_COMM_WORLD, self%mSize, iErr )
        self%message = message
    END FUNCTION
END MODULE

!--------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|-|
PROGRAM PolymorphicTest
    USE mpi
    USE MsgBase
    USE MsgCapabilities
    USE MsgOther
    IMPLICIT NONE

    INTEGER :: iErr
    INTEGER :: mRank
    INTEGER :: mSize

    !----|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|-|
    ! the BLOCK allows for object cleanup - mostly for when using valgrind
    BLOCK
        TYPE( MsgCapabilitiesType ), TARGET :: myObjC
        TYPE( MsgOtherType ), TARGET :: myObjO
        CLASS( MsgBaseType ), POINTER :: myObjPtr

        CALL MPI_INIT( iErr )
        IF( 0 == iErr ) CALL MPI_COMM_RANK( MPI_COMM_WORLD, mRank, iErr )
        IF( 0 == iErr ) CALL MPI_COMM_SIZE( MPI_COMM_WORLD, mSize, iErr )

        myObjC = MsgCapabilitiesType()
        myObjPtr => myObjC

        IF( 0 == mRank ) myObjC%baseId = 89 ! just to have something to track
        IF( 0 == mRank ) myObjC%hostname = "Fred"
        IF( 1 == mRank ) myObjC%hostname = "Barney"
        CALL sendReceive( myObjPtr, 0, 1 )

        myObjO = MsgOtherType( 'Pebbles' )
        myObjPtr => myObjO

        IF( 1 == mRank ) myObjO%baseId = 71
        IF( 1 == mRank ) myObjO%hostname = "Wilma"
        IF( 0 == mRank ) myObjO%hostname = "Bam-bam"
        IF( 0 == mRank ) myObjO%message = "Betty"
        CALL sendReceive( myObjPtr, 1, 0 )

        CALL MPI_FINALIZE( iErr )
    END BLOCK

CONTAINS

    SUBROUTINE sendReceive( objPtr, srcNode, destNode )
        CLASS( MsgBaseType ), POINTER :: objPtr
        INTEGER, PARAMETER :: tag = 23
        INTEGER, DIMENSION( MPI_STATUS_SIZE ) :: status
        INTEGER :: mpiDatatype, srcNode, destNode

        mpiDatatype = objPtr%getMpiDatatype()
        WRITE(*, '(I0, ":", I0, 1X, "     Has:", A)') mRank, mSize, objPtr%str()
        IF( srcNode == mRank ) THEN
            ! OpenMPI polymorphic datatypes reference the first item in the base class, mpich doesn't
            CALL MPI_SEND( objPtr%baseId, 1, mpiDatatype, destNode, tag, MPI_COMM_WORLD, iErr )
            WRITE(*, '(I0, ":", I0, 1X, "    Sent:", A)') mRank, mSize, objPtr%str()
        ELSE
            CALL MPI_RECV( objPtr%baseId, 1, mpiDatatype, MPI_ANY_SOURCE, tag, MPI_COMM_WORLD, status, iErr )
            WRITE(*, '(I0, ":", I0, 1X, "Received:", A)') mRank, mSize, objPtr%str()
        END IF
    END SUBROUTINE sendReceive
END PROGRAM PolymorphicTest

!--------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|-|
! The MPI for testing:
! $ mpiexec --help
!   mpiexec (OpenRTE) 1.10.2
!
! Execution on two nodes produces this print (your order of print may be slightly different)
! $ mpiexec -n 2 bin/PolymorphicFortran
! 0:2      Has:[0:2 id:89 Datatype:73 Name:Fred]
! 1:2      Has:[1:2 id:59 Datatype:73 Name:Barney]
! 0:2     Sent:[0:2 id:89 Datatype:73 Name:Fred]
! 0:2      Has:[id:59 Datatype:74 0:2 Name:Bam-bam Betty]
! 1:2 Received:[0:2 id:89 Datatype:73 Name:Fred]
! 1:2      Has:[id:71 Datatype:74 1:2 Name:Wilma Pebbles]
! 1:2     Sent:[id:71 Datatype:74 1:2 Name:Wilma Pebbles]
! 0:2 Received:[id:71 Datatype:74 1:2 Name:Wilma Pebbles]
!--------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|-|
