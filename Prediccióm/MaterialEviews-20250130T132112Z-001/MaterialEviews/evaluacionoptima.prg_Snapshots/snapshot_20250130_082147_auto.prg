' EvaluacionOptima.prg
' por Diego Winkelried (Universidad del Pacífico)
' Curso de Verano del BCRP - Economía Avanzada (Enero 2025)
' ----------------------------------------------------------------------------------------------------------------------
' Automatización del cómputo de pruebas de hipótesis E[e*e(-i)] = 0 para determinar si una proyección es 
' "aproximadamente" ("estadísticamente") óptima
' ----------------------------------------------------------------------------------------------------------------------
' Este programa debe ser ejecutado después de abrir haber generado una secuencia de errores de predicción de nombre 
' e{%tag}1, e{%tag}2, ..., e{%tag}H en el Workfile 
' ----------------------------------------------------------------------------------------------------------------------
' Nota al instructor: Previamente, haga una demostración del uso de simples pruebas de hipótesis para la media/mediana
' ----------------------------------------------------------------------------------------------------------------------

SMPL @ALL
!maxH = Hmax

TABLE tabla
tabla(1,1) = "h"
tabla(1,2) = "j"
tabla(1,3) = "cov[e{h}, e{h}(-j)]"
tabla(1,4) = "t-stat"
tabla(1,5) = " "

!q = 1
FOR !h = 1 TO !maxH
	!q = !q + 1
	tabla(!q + 1, 1) = !H
	SETLINE(tabla, !q)
	FOR !j = 0 TO !maxH+1
		!q = !q + 1
		IF !j = !maxH+1 THEN
			EQUATION Eee.LS(n) e{%tag}{!h}*y(-!h) C 
			tabla(!q, 2) = " y"
		ELSE
			EQUATION Eee.LS(n) e{%tag}{!h}*e{%tag}{!h}(-!j) C 
			tabla(!q, 2) = !j
		ENDIF
		tabla(!q, 3) = Eee.C(1)
		!tstat  = Eee.@TSTATS(1)
		tabla(!q, 4) = @ABS(!tstat) 
		IF @ABS(!tstat) > 2 THEN
			tabla(!q, 5) = "*"
		ENDIF 
	NEXT !h
	
NEXT !j
tabla.SETFORMAT(A:B) f.0
tabla.SETFORMAT(C:D) f.4
SETCOLWIDTH(tabla, 1, 4)
SETCOLWIDTH(tabla, 2, 4)
SETCOLWIDTH(tabla, 3, 15)
SETCOLWIDTH(tabla, 5, 4)

D Eval{%tag}* Eee
RENAME tabla Eval{%tag}
'SHOW Eval{%tag}


