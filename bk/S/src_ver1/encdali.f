      SUBROUTINE Encdali(val,nx,ny,naok,
     &     rawdata,nhitdata,hitdet,nhitdet,ndet,ndata)
c-----------------------------------------------------------------------
c ANALYZER 12 : Encdali
c-----------------------------------------------------------------------
c  ID=1,186 
c
c   W#: 1  2    3     4         5     6     7    8           9
c       ID Araw Traw  Traw-ref  Acal  Adop  Tcal Araw(Ttrue) Acal(Ttrue)
c       
c      11
c       LayerID
c
c  ID=200 Multiplicity (defined by 0 =< Araw =< 3840)
c   W#: 1  2
c       ID Multi
c  ID
c-----------------------------------------------------------------------

c
c
c
c
c     ihit_min                  : the number of required data (ihit_min0)
c     hitdet(1:nhitdet)         : hit detector id
c     nhitdata(1:ndet)          : the number of data for each id 
c     nhitdet                   : the number of hit detector
c     ndet                      : the number of total detectors
c     ndata                     : the length of rawdata array
      IMPLICIT NONE
      INCLUDE 'analyslogic.fh'
      INCLUDE 'daliprm.fh'

      INTEGER nx, ny, ndet, ndata, naok
      INTEGER rawdata(ndata,ndet)
      INTEGER nhitdata(ndet)
      INTEGER hitdet(ndet)
      INTEGER nhitdet
      REAL    val(nx,ny)

      CHARACTER ihitchara*4
      INTEGER ihit_min0

      INTEGER ihit, id
      INTEGER mult, tmult

      REAL    Cosd
      EXTERNAL Cosd

      SAVE ihit_min0

      IF (InitENCflag(12)) THEN
         CALL LOADDALIPRM
         CALL GETENV('IHIT_MIN0',ihitchara)
         READ(ihitchara,*) ihit_min0
         IF (ihit_min0.EQ.1) THEN
            WRITE(*,*) 'ENCDALI-M : IHIT_MIN0 = ', ihit_min0
         ELSE
            ihit_min0 = 3
            WRITE(*,*) 'ENCDALI-M : IHIT_MIN0 -> ', ihit_min0
         ENDIF
         InitENCFlag(12) = .FALSE.
c         DO i=1,186,1
c            write(*,*)i,ped(i)+gain(i)*3840.
c            write(*,*)i,(100.-ped(i))/gain(i),
c     &           int((100.-ped(i))/gain(i)/2.)
c            write(*,*)i, int((100.-ped(i))/gain(i)/2.)
c         ENDDO
      ENDIF

      mult = 0
      tmult = 0
      naok = 0

      DO ihit = 1,nhitdet
         id = hitdet(ihit)
         IF (id.GT.nch_dali) CYCLE
         IF (nhitdata(id).LT.ihit_min0) CYCLE ! ihit_min
         IF (rawdata(3,id).GT.1) CYCLE

         naok = naok + 1
         val(1,naok) = id
         val(2,naok) = rawdata(1,id) ! Araw
         val(3,naok) = rawdata(2,id)
         val(4,naok) = rawdata(2,id)
         val(5,naok) = ped(id) + rawdata(1,id) * gain(id) ! Acal
         val(6,naok) = val(5,naok) * gamma *(1.-beta*Cosd(theta(id)))! Adop

         val(7,naok) = val(4,naok) * ch2ns(id) + toff(id)! Tcal

         val(11,naok) = id2layer(id)

         IF (val(2,naok).GT.0. .AND. val(3,naok).GT.0.) THEN
            mult = mult + 1
         ENDIF

         val(15,naok) = rawdata(3,id) ! Multihit

         IF (val(7,naok).GT.truemin .AND. val(7,naok).LE.truemax .AND. 
     &       val(5,naok).GT.0.) THEN
            val(8,naok)=val(2,naok)
            val(9,naok)=val(5,naok)
            tmult = tmult + 1
         ELSE
            val(8,naok)=-1000.
            val(9,naok)=-1000.
         ENDIF

      ENDDO

ccc Multiplicity
      ID = 200

      naok = naok + 1
      val(1,naok) = id
      val(2,naok) = mult
      val(3,naok) = tmult

      RETURN
      END

c ======================================

      REAL FUNCTION Cosd(x)

      REAL x,pi

      pi = ACOS(-1.0)
      Cosd = COS(x/180.*pi)

      RETURN

      END
