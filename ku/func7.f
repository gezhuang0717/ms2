      FUNCTION FUNC7(X)
      COMMON/PAWPAR/PAR(4)

c      PAR(1) = life1
c      PAR(2) = const1
c      PAR(3) = bg
c      PAR(4) = tau*Y Y=beam Yield

   real  LN 
   data  LN/0.693147/
   real  tau
   tau = PAR(1)/LN

    If (X.lt.500) Then
      FUNC7 = PAR(4) - PAR(2)*EXP(-X/tau) + PAR(3)
    Else
      FUNC7 = (PAR(4)-PAR(2)*EXP(-500/tau))*EXP(-(X-500)/tau) + PAR(3)
    EndIf
    RETURN
    END
