#include "DOUBLEPLUME_OPTIONS.h"

CBOP
C     !ROUTINE: DOUBLEPLUME_SOLVE
C     !INTERFACE:
      SUBROUTINE DOUBLEPLUME_SOLVE( nPass, errOut, myThid )
C     !DESCRIPTION: \bv
C     *==========================================================*
C     | SUBROUTINE DOUBLEPLUME_SOLVE
C     | o Double plume solution for the current column (Python:
C     |   doublePlume.__init__, double_plume branch):
C     |   - if the ambient is not stably stratified anywhere, run
C     |     the inner plume only (no outer plume);
C     |   - otherwise iterate inner <-> outer until the inner
C     |     volume flux, the peel heights and the plume top change
C     |     by less than DP_convTol (at most DP_maxIter+1 passes).
C     | Requires the column set-up done in DOUBLEPLUME_CALC
C     | (taProf, saProf, dp_H, dp_nAmb, plumeLayers, Q_sg, r_sg).
C     *==========================================================*
C     nPass  :: number of inner/outer passes done
C     errOut :: final convergence error (Python err)
C     \ev
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "DOUBLEPLUME.h"
      INTEGER nPass, myThid
      DOUBLE PRECISION errOut

      INTEGER k, j, iter, nTop, nTopOld, nPk, nPkOld, kc
      LOGICAL stratified
      DOUBLE PRECISION QiOld(DP_maxFine)
      DOUBLE PRECISION hPk(DP_maxPeel), hPkOld(DP_maxPeel)
      DOUBLE PRECISION err, errQ, errZp, errTop, Hrise, sNum, sDen
      DOUBLE PRECISION wt, htmp
      CHARACTER*(MAX_LEN_MBUF) msgBuf
      DOUBLE PRECISION DP_RHO
      EXTERNAL DP_RHO
CEOP

      DO k = 1, plumeLayers
        QoProf_f(k) = 0.D0
        MoProf_f(k) = 0.D0
        HoProf_f(k) = 0.D0
        BoProf_f(k) = 0.D0
        CoProf_f(k) = 0.D0
      ENDDO
      dpo_nFine = plumeLayers
      nPass  = 0
      errOut = 0.D0

C--   Stratification (Python: nanmean(N) > 0, density at p = 0):
C     stratified if density increases with depth somewhere
      stratified = .FALSE.
      DO k = 1, dp_nAmb-1
        IF ( DP_RHO(taProf(k+1), saProf(k+1), 0.D0) .GT.
     &       DP_RHO(taProf(k),   saProf(k),   0.D0) )
     &       stratified = .TRUE.
      ENDDO

      IF ( .NOT. stratified ) THEN
        WRITE(msgBuf,'(A)') 'DOUBLEPLUME_SOLVE: no ambient '
     &   //'stratification, running the inner plume only'
        CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                      SQUEEZE_RIGHT, myThid )
        CALL DOUBLEPLUME_INNER_MODEL( myThid )
        nPass = 1
        RETURN
      ENDIF

      err     = 1.D3
      iter    = 0
      nTopOld = 0
      nPkOld  = 0
  100 CONTINUE
      IF ( err .GT. DP_convTol .AND. iter .LE. DP_maxIter ) THEN

        CALL DOUBLEPLUME_INNER_MODEL( myThid )
        CALL DOUBLEPLUME_OUTER_MODEL( myThid )
        nPass = nPass + 1

C--     Top of the inner plume and sorted peel heights
        nTop = 1
        DO k = 1, plumeLayers
          IF ( QiProf_f(k) .GT. 0.D0 ) nTop = k
        ENDDO
        nPk = j_peel
        DO j = 1, nPk
          hPk(j) = dp_H - z_peel(j)
        ENDDO
        DO j = 2, nPk
          DO k = j, 2, -1
            IF ( hPk(k) .LT. hPk(k-1) ) THEN
              htmp     = hPk(k)
              hPk(k)   = hPk(k-1)
              hPk(k-1) = htmp
            ENDIF
          ENDDO
        ENDDO

        IF ( iter .GT. 0 ) THEN
          Hrise = MAX( (nTop-1)*DP_dz, DP_dz )
C         integrated (L1) change of Qi on the common height range
          kc   = MIN( nTop, nTopOld )
          sNum = 0.D0
          sDen = 0.D0
          DO k = 1, kc
            wt = 1.D0
            IF ( k .EQ. 1 .OR. k .EQ. kc ) wt = 0.5D0
            sNum = sNum + wt*ABS( QiProf_f(k) - QiOld(k) )
            sDen = sDen + wt*ABS( QiOld(k) )
          ENDDO
          IF ( sDen .GT. 0.D0 ) THEN
            errQ = sNum/sDen
          ELSE
            errQ = 1.D0
          ENDIF
C         peel heights and plume top, relative to the rise height
          IF ( nPk .EQ. nPkOld ) THEN
            errZp = 0.D0
            DO j = 1, nPk
              errZp = MAX( errZp, ABS(hPk(j) - hPkOld(j))/Hrise )
            ENDDO
          ELSE
            errZp = 1.D0
          ENDIF
          errTop = ABS( nTop - nTopOld )*DP_dz/Hrise
          err = MAX( errQ, errZp, errTop )
          IF ( debugLevel .GE. debLevB ) THEN
            WRITE(msgBuf,'(A,I3,A,3ES10.2,A,I4)')
     &       'DOUBLEPLUME_SOLVE: iter', iter, ' dQ, dPeel, dTop =',
     &       errQ, errZp, errTop, '  peels =', nPk
            CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                          SQUEEZE_RIGHT, myThid )
          ENDIF
        ENDIF

C--     Keep this pass for the next comparison
        DO k = 1, plumeLayers
          QiOld(k) = QiProf_f(k)
        ENDDO
        DO j = 1, nPk
          hPkOld(j) = hPk(j)
        ENDDO
        nTopOld = nTop
        nPkOld  = nPk
        iter    = iter + 1
        GOTO 100
      ENDIF

      errOut = err
      IF ( iter .GT. DP_maxIter .AND. err .GT. DP_convTol ) THEN
        WRITE(msgBuf,'(A,ES10.2)') 'DOUBLEPLUME_SOLVE: max iterations'
     &   //' exceeded, keeping last pass; err =', err
        CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                      SQUEEZE_RIGHT, myThid )
      ENDIF

      RETURN
      END
