      FUNCTION FUNCQ(X)
      COMMON/PAWPAR(4)

c     PAR(1)=alinment
c     PAR(2)=Q-moment
c     PAR(3)=offset

      real A2
      data A2/-0.3962/

      FUNCQ=(3.*(A2*PAR(1)*0.200+0.116*COS(PAR(2)*X)
     &+0.033*COS(2.*PAR(2)*X)+0.165*COS(3.*PAR(2)*X)
     &+0.097*COS(4.*PAR(2)*X)+0.192*COS(5.*PAR(2)*X)
     &+0.086*COS(6.*PAR(2)*X)+0.060*COS(7.*PAR(2)*X)
     &+0.036*COS(9.*PAR(2)*X)+0.014*COS(11.*PAR(2)*X)))/
     &(4+(A2*PAR(1)*0.200+0.116*COS(PAR(2)*X)
     &+0.033*COS(2.*PAR(2)*X)+0.165*COS(3.*PAR(2)*X)
     &+0.097*COS(4.*PAR(2)*X)+0.192*COS(5.*PAR(2)*X)
     &+0.086*COS(6.*PAR(2)*X)+0.060*COS(7.*PAR(2)*X)
     &+0.036*COS(9.*PAR(2)*X)+0.014*COS(11.*PAR(2)*X)))
     &+PAR(3)
      RETURN
      END
