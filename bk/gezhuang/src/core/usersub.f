      Subroutine USERSUB(GateFlag)
c
c                                             2000.Oct.20. S.Takeuchi
c                                             2003.Nov.09. S.Takeuchi modified
c
c Do not EDIT!! -->
      Integer*2      EvtData(2000)
      Integer        nx/500/,ny/500/,naok(50)
      Integer        USERFGATE(10)
      Integer        RNum
      Real           val(500,500,50)
      Logical        AnalyzerFlag(50),INITENCFLAG(50)
      Logical        USERFLAG(10),EVTERR
      Logical        GateFlag(0:2000)
      Common/ANALYSEVENT/ EvtData,RNum
      Common/ANALYSVALUE/ val,naok
      Common/ANALYSLOGIC/ AnalyzerFlag,INITENCFLAG,USERFLAG,EVTERR
      Common/ANALYSUSER/  USERFGATE
c Do not EDIT!! <--
c
c

c
c UFLAG : 1  --  Hogehoge  --
c

      If( (GateFlag(USERFGATE(1))) ) Then
c
c         Call USER_HogeHoge(val(1,1,3),nx,ny,naok(3))
c
      EndIf

      Return

      End

