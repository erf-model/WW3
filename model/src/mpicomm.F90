!> @file
!> @brief Defines MPI communicator id as constants for global use.
!>
!> @author H. L. Tolman  @date 05-Jun-2018
!>
#include "w3macros.h"

!>
!> @brief Define some mpi constants for global use
!>
!> @author  H. L. Tolman  @date 05-Jun-2018
!>
!
#ifndef ENDIANNESS
#define ENDIANNESS "native"
#endif
!
!/ ------------------------------------------------------------------- /
MODULE MPICOMM
  !/
  !/                  +-----------------------------------+
  !/                  | WAVEWATCH III           NOAA/NCEP |
  !/                  |           H. L. Tolman            |
  !/                  |                        FORTRAN 90 |
  !/                  | Last update :         05-Jun-2018 |
  !/                  +-----------------------------------+
  !/
  !/    11-Nov-1999 : Fortran 90 version.                 ( version 2.00 )
  !/    29-May-2009 : Preparing distribution version.     ( version 3.14 )
  !/    25-Jun-2011 : Adding Kelvin functions.            ( version 4.05 )
  !/    03-Sep-2012 : Adding TSTOUT flag.                 ( version 4.10 )
  !/    28-Feb-2013 : Adding cap at 0.5 in FWTABLE        ( version 4.08 )
  !/    20-Jan-2017 : Add parameters for ESMF             ( version 6.02 )
  !/    01-Mar-2018 : Add UNDEF parameter                 ( version 6.02 )
  !/    05-Jun-2018 : Add PDLIB parameters                ( version 6.04 )
  !/
  !/    Copyright 2009-2012 National Weather Service (NWS),
  !/       National Oceanic and Atmospheric Administration.  All rights
  !/       reserved.  WAVEWATCH III is a trademark of the NWS.
  !/       No unauthorized use without permission.
  !/
  !  1. Purpose :
  !
  !     Define some mpi constants for global use
  !
  !  2. Variables and types :
  !
  !      Name      Type  Scope    Description
  !     ----------------------------------------------------------------
  !      UNDEF     Real  Global   Value for undefined variable in output
  !     ----------------------------------------------------------------
  !/ ------------------------------------------------------------------- /
  !/
  !
!#ifdef W3_MPI
!  INCLUDE "mpif.h"
!#endif

  INTEGER :: MPI_COMM_WW3=0 !< MPI_COMM_WW3
  !
  ! Parameters in support of running as ESMF component
  !
  ! --- Flag indicating whether or not the model has been invoked as an
  !     external Component.  This flag is set to true in the external
  !     module during initialization.
  LOGICAL :: IS_EXTERNAL_COMPONENT = .FALSE. !< IS_EXTERNAL_COMPONENT Flag for model invoked via external executable.
  !
CONTAINS

SUBROUTINE WW3_SEND_TO_ERF

    USE CONSTANTS
    USE W3GDATMD  ! HAS NX, NY
    USE W3ADATMD, ONLY: HS, WLM
    USE W3ODATMD, ONLY: NDST, UNDEF, IAPROC, NAPROC
    USE W3ADATMD, ONLY: NSEALM
    USE W3PARALL, ONLY : INIT_GET_ISEA, SYNCHRONIZE_GLOBAL_ARRAY
    ! USE W3IOGOMD,  ONLY: S2GRID 
#ifdef W3_PDLIB
    USE W3ODATMD, only : IAPROC, NAPROC, NTPROC
    USE W3ADATMD, ONLY: MPI_COMM_WCMP
    use yowDatapool, only: rtype, istatus
    USE yowNodepool, only: npa
    use yowNodepool, only: iplg
    use yowDatapool, only: rkind
#endif
    IMPLICIT NONE
#ifdef W3_MPI
  INCLUDE "mpif.h"
#endif


! PARAMETER LIST -------------------------------------------------- *
    INTEGER :: COUNTER
#define W3_MPMD
#ifdef W3_MPMD
  LOGICAL             :: FIRST_STEP = .TRUE., initialized, mpi_initialized_by_us
  integer             :: flag, myproc, nprocs, max_appnum, min_appnum, this_root, other_root, rank_offset, this_nboxes
  integer             :: p, appnum, all_appnum(10), napps, all_argc(10), IERR_MPI
  CHARACTER(LEN=80)   :: exename
  REAL, ALLOCATABLE       :: X1(:,:)

! MY EDITS
  INTEGER :: n_elements
   REAL(8), allocatable :: magnitude_values(:)
   REAL(8), allocatable :: theta_values(:)

#ifdef W3_PDLIB
  REAL(rkind)         :: XY_SEND(NX*NY)
  REAL(rkind)         :: XY_SYNCH_SEND(NSEA)
#else
  DOUBLE PRECISION    :: XY_SEND(NX*NY)
  DOUBLE PRECISION    :: XY_SYNCH_SEND(NSEA)
#endif
#endif

    INTEGER            :: JSEA, ISEA, IX, IY, I, J

! BEGIN SEND ---------------------------------------------------------- *

#ifdef W3_MPMD

#ifdef W3_MPI
  CALL MPI_COMM_SIZE ( MPI_COMM_WORLD, NPROCS, IERR_MPI )
#endif
#ifdef W3_MPI
  CALL MPI_COMM_RANK ( MPI_COMM_WORLD, MYPROC, IERR_MPI )
  MYPROC = MYPROC + 1
#endif

#ifdef W3_MPI
  print*, "My rank is ",MYPROC," out of ",NPROCS," total ranks in my part of MPI_COMM_WORLD communicator ",MPI_COMM_WORLD, "and my rank is ",IAPROC," out of ",NAPROC," total ranks in my part of the split communicator ", MPI_COMM_WW3

  rank_offset = MyProc - IAPROC;
  if (rank_offset .eq. 0) then ! First program
     this_root = 0
     other_root = NAPROC
  else
     this_root = rank_offset
     other_root = 0
  end if

  ALLOCATE(X1(NX+1,NY))
!  ALLOCATE(XY_SEND(NX*NY))
  if (MyProc-1 .eq. this_root) then
     if (rank_offset .eq. 0) then !  the first program
        CALL MPI_Send(NX, 1, MPI_INT, other_root, 0, MPI_COMM_WORLD, IERR_MPI)
        CALL MPI_Send(NY, 1, MPI_INT, other_root, 6, MPI_COMM_WORLD, IERR_MPI)
     else ! the second program
        CALL MPI_Send(NX, 1, MPI_INT, other_root, 1, MPI_COMM_WORLD, IERR_MPI)
        CALL MPI_Send(NY, 1, MPI_INT, other_root, 7, MPI_COMM_WORLD, IERR_MPI)
     end if
  end if

  if (MyProc-1 .eq. this_root) then
     if (rank_offset .eq. 0) then !  the first program
        X1     = UNDEF
        XY_SEND     = UNDEF
        ! CALL S2GRID(HS, X1)
        XY_SYNCH_SEND = HS
        CALL SYNCHRONIZE_GLOBAL_ARRAY(XY_SYNCH_SEND)
 
        DO JSEA=1, NSEA
           CALL INIT_GET_ISEA(ISEA, JSEA)
           IX     = MAPSF(ISEA,1)
           IY     = MAPSF(ISEA,2)
           XY_SEND((IX)+(IY-1)*NX)=XY_SYNCH_SEND(ISEA)
        END DO

        CALL MPI_Send(XY_SEND, NX*NY, MPI_DOUBLE, other_root, 2, MPI_COMM_WORLD, IERR_MPI)
        X1     = UNDEF
        XY_SYNCH_SEND = WLM
        CALL SYNCHRONIZE_GLOBAL_ARRAY(XY_SYNCH_SEND)
        DO JSEA=1, NSEA
           CALL INIT_GET_ISEA(ISEA, JSEA)
           IX     = MAPSF(ISEA,1)
           IY     = MAPSF(ISEA,2)
           XY_SEND((IX)+(IY-1)*NX)=XY_SYNCH_SEND(ISEA)
        END DO
        CALL MPI_Send(XY_SEND, NX*NY, MPI_DOUBLE, other_root, 4, MPI_COMM_WORLD, IERR_MPI)
     else ! the second program
        X1     = UNDEF
        XY_SEND     = UNDEF
        XY_SYNCH_SEND = HS
        CALL SYNCHRONIZE_GLOBAL_ARRAY(XY_SYNCH_SEND)
        DO JSEA=1, NSEA
           CALL INIT_GET_ISEA(ISEA, JSEA)
           IX     = MAPSF(ISEA,1)
           IY     = MAPSF(ISEA,2)
           XY_SEND((IX)+(IY-1)*NX)=XY_SYNCH_SEND(ISEA)
        END DO
        CALL MPI_Send(XY_SEND, NX*NY, MPI_DOUBLE, other_root, 3, MPI_COMM_WORLD, IERR_MPI)
        X1     = UNDEF
        XY_SYNCH_SEND = WLM
        CALL SYNCHRONIZE_GLOBAL_ARRAY(XY_SYNCH_SEND)
        DO JSEA=1, NSEA
           CALL INIT_GET_ISEA(ISEA, JSEA)
           IX     = MAPSF(ISEA,1)
           IY     = MAPSF(ISEA,2)
           XY_SEND((IX)+(IY-1)*NX)=XY_SYNCH_SEND(ISEA)
        END DO
        CALL MPI_Send(XY_SEND, NX*NY, MPI_DOUBLE, other_root, 5, MPI_COMM_WORLD, IERR_MPI)
     end if
  end if

! MY EDITS HERE
! CHECK XY_SYNCH_SEND, SYNCH_GLOBAL_ARRAY
    OPEN(5120, file='printmpi.txt', status='unknown', access='append', action="write")

    ! Write HS values to the new file
    DO JSEA=1, NSEAL
        CALL INIT_GET_ISEA(ISEA, JSEA)
        IX     = MAPSF(ISEA,1)
        IY     = MAPSF(ISEA,2)

        WRITE(5120, *) SIZE(XY_SEND), XY_SEND(ISEA), SIZE(XY_SYNCH_SEND), XY_SYNCH_SEND(ISEA)
    END DO
    CLOSE(5120)
  DEALLOCATE(X1)
#endif
#endif

! END SEND ----------------------------------------------------------------- *

END SUBROUTINE WW3_SEND_TO_ERF

SUBROUTINE WW3_RECEIVE_FROM_ERF

    USE CONSTANTS
    USE W3GDATMD
    USE W3ADATMD, ONLY: HS, WLM
    USE W3ODATMD, ONLY: NDST, UNDEF, IAPROC, NAPROC
    USE W3ADATMD, ONLY: NSEALM
    USE W3PARALL, ONLY : INIT_GET_ISEA, SYNCHRONIZE_GLOBAL_ARRAY
    ! USE W3IOGOMD,  ONLY: S2GRID
#ifdef W3_PDLIB
    USE W3ODATMD, only : IAPROC, NAPROC, NTPROC
    USE W3ADATMD, ONLY: MPI_COMM_WCMP
    use yowDatapool, only: rtype, istatus
    USE yowNodepool, only: npa
    use yowNodepool, only: iplg
    use yowDatapool, only: rkind
#endif
    IMPLICIT NONE
#ifdef W3_MPI
  INCLUDE "mpif.h"
#endif
! PARAMETER LIST -------------------------------------------------- *
    INTEGER :: COUNTER
#define W3_MPMD
#ifdef W3_MPMD
  LOGICAL             :: FIRST_STEP = .TRUE., initialized, mpi_initialized_by_us
  integer             :: flag, myproc, nprocs, max_appnum, min_appnum, this_root, other_root, rank_offset, this_nboxes
  integer             :: p, appnum, all_appnum(10), napps, all_argc(10), IERR_MPI
  CHARACTER(LEN=80)   :: exename
  REAL, ALLOCATABLE       :: X1(:,:)

! MY EDITS
  INTEGER :: n_elements
   REAL(8), allocatable :: magnitude_values(:)
   REAL(8), allocatable :: theta_values(:)

#ifdef W3_PDLIB
  REAL(rkind)         :: XY_SEND(NX*NY)
  REAL(rkind)         :: XY_SYNCH_SEND(NSEA)
#else
  DOUBLE PRECISION    :: XY_SEND(NX*NY)
  DOUBLE PRECISION    :: XY_SYNCH_SEND(NSEA)
#endif
#endif

    INTEGER            :: JSEA, ISEA, IX, IY, I, J

#ifdef W3_MPMD

#ifdef W3_MPI
  CALL MPI_COMM_SIZE ( MPI_COMM_WORLD, NPROCS, IERR_MPI )
#endif
#ifdef W3_MPI
  CALL MPI_COMM_RANK ( MPI_COMM_WORLD, MYPROC, IERR_MPI )
  MYPROC = MYPROC + 1
#endif

#ifdef W3_MPI
  print*, "My rank is ",MYPROC," out of ",NPROCS," total ranks in my part of MPI_COMM_WORLD communicator ",MPI_COMM_WORLD, "and my rank is ",IAPROC," out of ",NAPROC," total ranks in my part of the split communicator ", MPI_COMM_WW3

  rank_offset = MyProc - IAPROC;
  if (rank_offset .eq. 0) then ! First program
     this_root = 0
     other_root = NAPROC
  else
     this_root = rank_offset
     other_root = 0
  end if

! BEGIN RECEIVE FROM ERF ---------------------------------------------------- *

n_elements = NX * NY
PRINT *, "ABOUT TO RECEIVE FROM ERF"
  if (MyProc-1 .eq. this_root) then
     if (rank_offset .eq. 0) then !  the first program

        CALL MPI_RECV( n_elements, 1, MPI_INT, other_root, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE, IERR_MPI );


        allocate(magnitude_values(n_elements))
        allocate(theta_values(n_elements))

        CALL MPI_RECV(magnitude_values, n_elements, MPI_DOUBLE, other_root, 12, MPI_COMM_WORLD,MPI_STATUS_IGNORE, IERR_MPI)
        CALL MPI_RECV(theta_values, n_elements, MPI_DOUBLE, other_root, 14, MPI_COMM_WORLD, MPI_STATUS_IGNORE,IERR_MPI)
     else ! the second program

        CALL MPI_RECV( n_elements, 1, MPI_INT, other_root, 11, MPI_COMM_WORLD, MPI_STATUS_IGNORE,IERR_MPI );

        allocate(magnitude_values(n_elements))
        allocate(theta_values(n_elements))

        call MPI_RECV(magnitude_values, n_elements, MPI_DOUBLE, other_root, 13, MPI_COMM_WORLD,MPI_STATUS_IGNORE, IERR_MPI)
        call MPI_RECV(theta_values, n_elements, MPI_DOUBLE, other_root, 15, MPI_COMM_WORLD, MPI_STATUS_IGNORE, IERR_MPI)
     end if
  end if


    print*, "JUST RECEIVED AND ALLOCATED MAG_VALUES(n-elements)"! MPI RECEIVE TEST
#endif
#endif

! END RECEIVE ------------------------------------------------------------------ *
END SUBROUTINE WW3_RECEIVE_FROM_ERF

  !/
  !/ End of module MPICOMM ------------------------------------------- /
  !/
END MODULE MPICOMM
