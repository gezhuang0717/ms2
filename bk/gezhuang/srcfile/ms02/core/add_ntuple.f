c =========================================================================
      Subroutine NTBOOK

c -- from ANALYS Source !! Do Not Edit !! ------------------------------- c
      Logical  AnalyzerFlag(50),INITENCFLAG(50)
      Logical  USERFLAG(10),EVTERR
      Common/ANALYSLOGIC/  AnalyzerFlag,INITENCFLAG,USERFLAG,EVTERR
c ----------------------------------------------------------------------- c

c -- for ntuple defined by user ----------------------------------------- c

c -- Analyzer 1
      Integer   RunNum
      Integer   Coin
      Common/CoinReg/ RunNum,Coin
c -- Analyzer 2
      Integer   AnaRun1
      Real      F2Tl,F2Tr,F3Tl,F3Tr,RF1,RF2,
     &          F2Al,F2Ar,F3Al,F3Ar,SSDA
      Common/BeamLine/ AnaRun1,F2TL,F2TR,F3TL,F3TR,
     &                 RF1,RF2,F2AL,F2AR,F3AL,F3AR,SSDA
c -- Analyzer 3
      Integer   AnaRun2
      Real      P1X1,P1X2,P1Y1,P1Y2,P2X1,P2X2,P2Y1,P2Y2,
     &          PAT,PBT,PAA,PBA
      Common/PPAC/ AnaRun2,P1X1,P1X2,P1Y1,P1Y2,P2X1,P2X2,P2Y1,P2Y2,
     &             PAT,PBT,PAA,PBA


c == Booking Part ==

      If(AnalyzerFlag(1)) Then
         Call HBNAME(10, 'CoinReg', RunNum,'RunNum:I,Coin:I')
      EndIf
      If(AnalyzerFlag(2)) Then
         Call HBNAME(10, 'BeamLine', AnaRun1,
     &        'AnaRun1:I,F2TL, F2TR, F3TL, F3TR, '//
     &        'RF1, RF2 ,F2AL, F2AR, F3AL, F3AR,SSDA')
      EndIf
      If(AnalyzerFlag(3)) Then
         Call HBNAME(10, 'PPAC', AnaRun2,
     &        'AnaRun2:I,P1X1,P1X2,P1Y1,P1Y2,P2X1,P2X2,P2Y1,P2Y2, '//
     &        'PAT, PBT, PAA, PBA')
      EndIf

      Return

      End

c =========================================================================
      Subroutine Add_Ntuple(IOFLAG)

c -- for Local ---------------------------------------------------------- c
      Logical  IOFLAG   ! ( True : Next File )

c -- from ANALYS Source !! Do Not Edit !! ------------------------------- c
      Real     val(500,500,50)
      Integer  naok(50)
      Integer  RunNumber
      Logical  AnalyzerFlag(50),INITENCFLAG(50)
      Logical  USERFLAG(10),EVTERR
      Common/ANALYSVALUE/  val,naok
      Common/AnalysNtuple/ NTUPLE_FIRST,RunNumber
      Common/ANALYSLOGIC/  AnalyzerFlag,INITENCFLAG,USERFLAG,EVTERR
c ----------------------------------------------------------------------- c

c -- for ntuple defined by user ----------------------------------------- c

c -- Analyzer 1
      Integer   RunNum
      Integer   Coin
      Common/CoinReg/ RunNum,Coin

c -- Analyzer 2
      Integer   AnaRun1
      Real      F2Tl,F2Tr,F3Tl,F3Tr,RF1,RF2,
     &          F2Al,F2Ar,F3Al,F3Ar,SSDA
      Common/BeamLine/ AnaRun1,F2TL,F2TR,F3TL,F3TR,
     &                 RF1,RF2,F2AL,F2AR,F3AL,F3AR,SSDA

c -- Analyzer 3
      Integer   AnaRun2
      Real      P1X1,P1X2,P1Y1,P1Y2,P2X1,P2X2,P2Y1,P2Y2,
     &          PAT,PBT,PAA,PBA
      Common/PPAC/ AnaRun2,P1X1,P1X2,P1Y1,P1Y2,P2X1,P2X2,P2Y1,P2Y2,
     &             PAT,PBT,PAA,PBA

c ----------------------------------------------------------------------- c

      If(IOFLAG)Then                  !! Do Not Edit !!
         Call Ntuple_io(0)            !! Do Not Edit !!
         Call Ntuple_io(1)            !! Do Not Edit !!
      EndIf                           !! Do Not Edit !!

c ----------------------------------------------------------------------- c
c               Definition of val(Word,ID,Analyzer)                       c
c ----------------------------------------------------------------------- c

c -- Analyzer 1
      If(AnalyzerFlag(1)) Then

      EndIf
c -- Analyzer 2
      If(AnalyzerFlag(2)) Then
         AnaRun1 = RunNum
      EndIf

c -- Analyzer 3
      If(AnalyzerFlag(3)) Then
         AnaRun2 = RunNum
      EndIf

      Return
      
      End




