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
!  CONTAINS

  !/
  !/ End of module MPICOMM ------------------------------------------- /
  !/
END MODULE MPICOMM
