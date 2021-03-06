
USE ATLAS_MODULE, ONLY : ATLAS_LIBRARY, ATLAS_LOG, ATLAS_LAMBERTREGIONALGRID, &
                       & ATLAS_GRID

USE FA_MOD, ONLY : FA_COM_DEFAULT, NEW_FA_DEFAULT, FACADR

USE PARKIND1, ONLY : JPRB, JPIM

IMPLICIT NONE

#include "abor1.intfb.h"

REAL (KIND=JPRB), PARAMETER :: RPI = 2.0_JPRB * ASIN (1.0_JPRB), &
                             & R2PI = 2._JPRB * RPI,             &
                             & RAD2DEG = 180._JPRB / RPI,        &
                             & DEG2RAD = RPI / 180._JPRB

CHARACTER (LEN=*), PARAMETER :: CLNOMC = 'c'
INTEGER (KIND=JPIM) :: ILUN, IREP, INBARP, INBARI

TYPE (ATLAS_LAMBERTREGIONALGRID) :: GRID

INTEGER, PARAMETER :: NX = 64, NY = 64;
INTEGER, PARAMETER :: NUX = 53, NCX = 8, NUY = 53, NCY = 8;

REAL (8), PARAMETER :: LADINDEGREES = 67.36, LOVINDEGREES = 26.64;
REAL (8), PARAMETER :: DXINMETRES = 2500, DYINMETRES = 2500;
REAL (8), PARAMETER :: LATIN1INDEGREES = 67.36, LATIN2INDEGREES = 67.36;

INTEGER :: IX, IY

CHARACTER (LEN=64) :: CLFILENAME

CALL GETARG (1, CLFILENAME)

CALL ATLAS_LIBRARY%INITIALISE()

ILUN = 77_JPIM
INBARP = 0
INBARI = 0
CALL FAITOU (IREP, ILUN, .TRUE., TRIM (CLFILENAME), 'OLD', &
           & .TRUE., .TRUE., 2_JPIM, INBARP, INBARI, CLNOMC)

BLOCK
  INTEGER (KIND=JPIM) :: IRANGC
  TYPE (FACADR), POINTER :: YLCADR
  REAL (KIND=JPRB) :: DXINMETRES, DYINMETRES
  REAL (KIND=JPRB) :: LOVINDEGREES, LADINDEGREES, LATIN1INDEGREES, LATIN2INDEGREES
  INTEGER (KIND=JPIM) :: NUX, NUY

  CALL FANUCA (CLNOMC, IRANGC, .FALSE.)

  YLCADR => FA_COM_DEFAULT%CADRE (IRANGC)

  DXINMETRES = YLCADR%SINLAT (7)
  DYINMETRES = YLCADR%SINLAT (8)
  NUX = YLCADR%NLOPAR (4)
  NUY = YLCADR%NLOPAR (6)

  LADINDEGREES    = YLCADR%SINLAT (4) * RAD2DEG
  LATIN1INDEGREES = YLCADR%SINLAT (4) * RAD2DEG
  LATIN2INDEGREES = YLCADR%SINLAT (4) * RAD2DEG
  LOVINDEGREES    = YLCADR%SINLAT (3) * RAD2DEG

  GRID = ATLAS_LAMBERTREGIONALGRID (YLCADR%NXLOPA, YLCADR%NLATIT, -NUX / 2 * DXINMETRES, -NUY / 2 * DYINMETRES, &
                                  & DXINMETRES, DYINMETRES, LOVINDEGREES, LADINDEGREES, LATIN1INDEGREES, LATIN2INDEGREES)

ENDBLOCK

CALL FAIRME (IREP, ILUN, 'KEEP')

PRINT *, " SIZE = ", GRID%SIZE ()

DO IY = 1, GRID%NY ()
  DO IX = 1, GRID%NX (IY)
    WRITE (*, '(2F20.10)') GRID%LONLAT (IX, IY)
  ENDDO
ENDDO


CALL GRID%FINAL ()

CALL ATLAS_LIBRARY%FINALISE()

END 

