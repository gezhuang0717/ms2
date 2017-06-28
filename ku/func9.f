      FUNCTION FUNC9(X)
      COMMON/PAWPAR/PAR(3)

c      PAR(1) = life1
c      PAR(2) = bg
c      PAR(3) = tau*Y Y=beam Yield

   real  LN 
   data  LN/0.693147/
   real  tau
   real  A

   tau = PAR(1)/LN 
c   A= PAR(3)/( 1 + EXP(-500/tau) )
   A= PAR(3)/( 1 + EXP(-2900/tau) )

c    If (X.lt.500) Then
c      FUNC9 = PAR(3) - A*EXP(-X/tau) + PAR(2)
c    Else
c      FUNC9 = (PAR(3)-A*EXP(-500/tau))*EXP(-(X-500)/tau) + PAR(2)
c    EndIf

    If (X.lt.2900) Then
      FUNC9 = PAR(3) - A*EXP(-X/tau) + PAR(2)
    Else
      FUNC9 = (PAR(3)-A*EXP(-2900/tau))*EXP(-(X-2900)/tau) + PAR(2)
    EndIf
    RETURN
    END
