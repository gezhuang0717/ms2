c
c ----- LOADDALIPRM ------- read the parameters for DALI
c
c

      SUBROUTINE Loaddaliprm
      IMPLICIT NONE
      INCLUDE "daliprm.fh"
      REAL    tmp(2)
      INTEGER ier ! error flag
      CHARACTER*256 filedaliprm


c read parameters for dali2
      
      CALL GETENV('DALI_PRM',filedaliprm)
      OPEN(UNIT=80, FILE=filedaliprm, STATUS='old', ERR=110)  

c pedestal      
      CALL Read_flt_list(80,ped,nch_dali,ier)
      IF (ier.NE.0) GOTO 125

c gain
      CALL Read_flt_list(80,gain,nch_dali,ier)
      IF (ier.NE.0) GOTO 125

c Tcal
      CALL Read_flt_list(80,ch2ns,nch_dali,ier)
      IF (ier.NE.0) GOTO 125

c Toffset
      CALL Read_flt_list(80,toff,nch_dali,ier)
      IF (ier.NE.0) GOTO 125

c Doppler corrction
      CALL Read_flt_list(80,theta,nch_dali,ier)
      IF (ier.NE.0) GOTO 125
      CALL Read_flt(80,tmp,2,ier)
      IF (ier.NE.0) GOTO 125
      beta = tmp(1)
      gamma = tmp(2)
      CALL Read_flt(80,tmp,2,ier)
      IF(ier.NE.0) GOTO 125
      truemin = tmp(1)
      truemax = tmp(2)

c parameters for conversion from ID to Layer
      CALL Read_int_list(80,id2layer,nch_dali,ier)
      IF (ier.NE.0) GOTO 125

      WRITE(*,*) 'ENC_DALI-M: DALI_PRM read in.'
      WRITE(*,*) ' '
      CLOSE(80)
      GOTO 130

 110  WRITE(*,*) 'ENC_DALI-E: DALI_PRM is not found.'
      GOTO 130
      STOP

 125  WRITE(*,*)' '
      WRITE(*,*)' ENC_DALI-E: Error in DALI_PRM.'
      WRITE(*,*)' '
      CLOSE(80)
      STOP

 130  CONTINUE

      RETURN
      END
