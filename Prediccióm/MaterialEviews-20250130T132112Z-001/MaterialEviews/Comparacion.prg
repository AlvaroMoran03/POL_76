' Comparacion.prg
' por Diego Winkelried (Universidad del Pacífico)
' Curso de Verano del BCRP - Economía Avanzada (Enero 2025)
' ----------------------------------------------------------------------------------------------------------------------
' Automatización del cómputo de pruebas de Diebold y Mariano para varias funciones de pérdida. Asimismo, se estiman la
' "regresión de combinación" y se realizan pruebas de hipótesis sobre "encompassing"
' ----------------------------------------------------------------------------------------------------------------------
' Este programa debe ser ejecutado después de abrir haber generado una secuencia de errores de predicción de nombres 
' e{%tag1}1, e{%tag1}2, ..., e{%tag1}H y e{%tag2}1, e{%tag2}2, ..., e{%tag2}H en el Workfile 
' ----------------------------------------------------------------------------------------------------------------------
' Nota al instructor: Previamente, discutir teoría sobre "encompassing" y Diebold y Mariano
' ----------------------------------------------------------------------------------------------------------------------

%tag1 = "a" ' La que se sospecha es la mejor
%tag2 = "c"

!maxH = Hmax

D Gy{%tag1}{%tag2}* dm*{%tag1}{%tag2}* comb*{%tag1}{%tag2}* tabla* Comp{%tag1}{%tag2}*

' Grupo 1a: Proyecciones para un horizonte seleccionado
GROUP Gy{%tag1}{%tag2}1a y{%tag1}1 y{%tag2}1  
Gy{%tag1}{%tag2}1a.SCAT LINEFIT

' Grupo 1b: Errores de predicción para un horizonte seleccionado
GROUP Gy{%tag1}{%tag2}1b e{%tag1}1 e{%tag2}1  
Gy{%tag1}{%tag2}1b.SCAT LINEFIT

' Grupo 2a: Proyecciones para otro horizonte seleccionado
GROUP Gy{%tag1}{%tag2}2a y{%tag1}{!maxH} y{%tag2}{!maxH}  
Gy{%tag1}{%tag2}2a.SCAT LINEFIT

' Grupo 2b: Errores de predicción para otro horizonte seleccionado
GROUP Gy{%tag1}{%tag2}2b e{%tag1}{!maxH} e{%tag2}{!maxH}  
Gy{%tag1}{%tag2}2b.SCAT LINEFIT

TABLE tabla

' Diebold y Mariano
' -----------------
SMPL @ALL
!h = 1
tabla(1,2) = %tag2
tabla(1,3) = %tag1
tabla(1,4) = "DM = (" + %tag2 + " - " + %tag1 + ")"
tabla(1,5) = "t-stat"
SETLINE(tabla, 2)
!q = 3
SERIES Ce1 = e{%tag1}{!h}^2
SERIES Ce2 = e{%tag2}{!h}^2
EQUATION dm.LS(n) Ce2-Ce1 C
tabla(!q,1) = "Cuadrática"
tabla(!q,2) = @MEAN(Ce2)
tabla(!q,3) = @MEAN(Ce1)
tabla(!q,4) = dm.C(1)
tabla(!q,5) = dm.@TSTATS(1)
RENAME dm dm1{%tag1}{%tag2}

!q = !q + 1
SERIES Ce1 = @ABS(e{%tag1}{!h})
SERIES Ce2 = @ABS(e{%tag2}{!h})
EQUATION dm.LS(n) Ce2-Ce1 C
tabla(!q,1) = "Absoluto"
tabla(!q,2) = @MEAN(Ce2)
tabla(!q,3) = @MEAN(Ce1)
tabla(!q,4) = dm.C(1)
tabla(!q,5) = dm.@TSTATS(1)
RENAME dm dm2{%tag1}{%tag2}

!q = !q + 1
SERIES Ce1 = @EXP(e{%tag1}{!h}/4) - e{%tag1}{!h}/4 - 1
SERIES Ce2 = @EXP(e{%tag2}{!h}/4) - e{%tag2}{!h}/4 - 1
EQUATION dm.LS(n) Ce2-Ce1 C
tabla(!q,1) = "LinEx"
tabla(!q,2) = @MEAN(Ce2)
tabla(!q,3) = @MEAN(Ce1)
tabla(!q,4) = dm.C(1)
tabla(!q,5) = dm.@TSTATS(1)
RENAME dm dm3{%tag1}{%tag2}
SETLINE(tabla, !q+1)

' Combination and Encompassing
' ----------------------------
!q = !q + 2

tabla(!q, 1) = "y = (1 - a)*y" + %tag1 + " + a*y" + %tag2 
!q = !q+1
EQUATION comb.LS e{%tag1}{!h} (e{%tag1}{!h}-e{%tag2}{!h})
tabla(!q,1) = "a estimado"
tabla(!q,2) = comb.C(1)
!q = !q+1
tabla(!q,1) = "H0: a = 0.0 [t-stat]" 
tabla(!q,3) = @ABS(comb.C(1)/comb.@STDERRS(1))
tabla(!q,4) = %tag1 + " domina"
!q = !q+1
tabla(!q,1) = "H0: a = 0.5 [t-stat]" 
tabla(!q,3) = @ABS((comb.C(1) - 0.5)/comb.@STDERRS(1))
!q = !q+1
tabla(!q,1) = "H0: a = 1.0 [t-stat]" 
tabla(!q,3) = @ABS((comb.C(1) - 1)/comb.@STDERRS(1))
tabla(!q,4) = %tag2 + " domina"
SETLINE(tabla, !q+1)

RENAME comb comb{%tag1}{%tag2}
RENAME tabla Comp{%tag1}{%tag2}


