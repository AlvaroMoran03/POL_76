' Robusto.prg
' por Diego Winkelried (Universidad del Pacífico)
' Curso de Verano del BCRP - Economía Avanzada (Enero 2025)
' ----------------------------------------------------------------------------------------------------------------------
' Ilustración sobre estrategias para "robustificar" las proyecciones ante la posible presencia de cambios estrcuturales.
' En particular, la eliminación del "término de corrección de errores" y el uso de modelos en diferencias
' ----------------------------------------------------------------------------------------------------------------------

'CLOSE Ejemplo3
'WFOPEN Ejemplo3

' Los datos
y.LINE

' Parámetros
!horizonte = 20
!maxfin = 116
!minfin = 80

' Tabla con resultados sobre error cuadrático medio e importancia del sesgo
' Ojo: Estos cálculos, que son automáticos de EViews, son a lo largo de "h"
'      y no son comparables con los vistos en secciones anteriores, que eran a lo largo de "t"
TABLE RECM
RECM.SETFORMAT(A) f.0
RECM.SETFORMAT(B:E) f.4
RECM(1,1) = "Origen"

FOR !fin = !minfin TO !maxfin STEP 2
	GROUP G{!fin} y
NEXT

' Proyecciones del modelo AR(1) en niveles (reparametrizado tipo Dikcey-Fuller) 
' -----------------------------------------------------------------------------
%tag = "AR"
GROUP y_{%tag} y
GROUP e_{%tag} 0
!col = 2
!row = 2
RECM(1,!col) = %tag
FOR !fin = !minfin TO 150 STEP 2
	SMPL 1 !fin
	EQUATION eq.LS D(y) C y(-1)
	SMPL !fin+1 !fin+!horizonte
	IF !fin <= !maxfin THEN
		FREEZE(tab_{%tag}{!fin}) eq.FORECAST(F=NA,E) y_{%tag}{!fin}
		tab_{%tag}{!fin}.SETFORMAT(B) f.3
		!row = !row + 1
		RECM(!row, 1) = @STR(!fin)
		RECM(!row, !col) = tab_{%tag}{!fin}(6, 2)
		!row = !row + 1
		RECM(!row, !col) = "[" +  @STR(tab_{%tag}{!fin}(10, 2)) + "]"
		G{!fin}.ADD y_{%tag}{!fin}
	ELSE
		eq.FORECAST(F=NA) y_{%tag}{!fin}
	ENDIF
	y_{%tag}.ADD y_{%tag}{!fin}
	e_{%tag}.ADD y-y_{%tag}{!fin}
NEXT !fin
y_{%tag}.LINE
CLOSE y_{%tag}
e_{%tag}.LINE
CLOSE e_{%tag}

' Proyecciones tipo "random walk" o "no change"
' ---------------------------------------------
%tag = "RW"
GROUP y_{%tag} y
GROUP e_{%tag} 0
!col = 3
!row = 2
RECM(1,!col) = %tag
FOR !fin = !minfin TO 150 STEP 2
	SMPL 1 !fin
	SERIES intercepto = 1
	EQUATION eq.LS D(y) intercepto 
	SMPL !fin+1 !fin+!horizonte
	' Intercepto = 0 en el horizonte de proyección para asegurar que D(y) = 0
 	intercepto = 0 
	IF !fin <= !maxfin THEN
		FREEZE(tab_{%tag}{!fin}) eq.FORECAST(F=NA,E) y_{%tag}{!fin}
		tab_{%tag}{!fin}.SETFORMAT(B) f.3
		!row = !row + 1
		RECM(!row, !col) = tab_{%tag}{!fin}(6, 2)
		!row = !row + 1
		RECM(!row, !col) = "[" +  @STR(tab_{%tag}{!fin}(10, 2)) + "]"
		G{!fin}.ADD y_{%tag}{!fin}
	ELSE
		eq.FORECAST(F=NA) y_{%tag}{!fin}
	ENDIF
	y_{%tag}.ADD y_{%tag}{!fin}
	e_{%tag}.ADD y-y_{%tag}{!fin}
NEXT !fin
y_{%tag}.LINE
e_{%tag}.LINE

' Proyecciones del modelo AR(1) en diferencias
' --------------------------------------------
%tag = "DAR"
GROUP y_{%tag} y
GROUP e_{%tag} 0
!col = 4
!row = 2
RECM(1,!col) = %tag
FOR !fin = !minfin TO 150 STEP 2
	SMPL 1 !fin
	EQUATION eq.LS D(y) D(y(-1))
	SMPL !fin+1 !fin+!horizonte
	IF !fin <= !maxfin THEN
		FREEZE(tab_{%tag}{!fin}) eq.FORECAST(F=NA,E) y_{%tag}{!fin}
		tab_{%tag}{!fin}.SETFORMAT(B) f.3
		!row = !row + 1
		RECM(!row, !col) = tab_{%tag}{!fin}(6, 2)
		!row = !row + 1
		RECM(!row, !col) = "[" +  @STR(tab_{%tag}{!fin}(10, 2)) + "]"
		G{!fin}.ADD y_{%tag}{!fin}
	ELSE
		eq.FORECAST(F=NA) y_{%tag}{!fin}
	ENDIF
	y_{%tag}.ADD y_{%tag}{!fin}
	e_{%tag}.ADD y-y_{%tag}{!fin}
NEXT !fin
y_{%tag}.LINE
e_{%tag}.LINE

' Proyecciones del modelo AR(1) en diferencias
' --------------------------------------------
'%tag = "IC"
'GROUP y_{%tag} y
'GROUP e_{%tag} 0
'!col = 4
'!row = 2
'RECM(1,!col) = %tag
'FOR !fin = !minfin TO 150 STEP 2
'	SMPL @ALL 
'	SERIES dum = 0
'	SMPL !fin @LAST
'	dum = 1
'	SMPL 1 !fin
'	EQUATION eq.LS y C y(-1) dum
'	SMPL !fin+1 !fin+!horizonte
'	IF !fin <= !maxfin THEN
'		FREEZE(tab_{%tag}{!fin}) eq.FORECAST(F=NA,E) y_{%tag}{!fin}
'		tab_{%tag}{!fin}.SETFORMAT(B) f.3
'		!row = !row + 1
'		RECM(!row, !col) = tab_{%tag}{!fin}(6, 2)
'		!row = !row + 1
'		RECM(!row, !col) = "[" + tab_{%tag}{!fin}(10, 2) + "]"
'		G{!fin}.ADD y_{%tag}{!fin}
'	ELSE
'		eq.FORECAST(F=NA) y_{%tag}{!fin}
'	ENDIF
'	y_{%tag}.ADD y_{%tag}{!fin}
'	e_{%tag}.ADD y-y_{%tag}{!fin}
'NEXT !fin
'y_{%tag}.LINE
'e_{%tag}.LINE

' Gráficos comparativos por origen de la proyecciòn
%juntar = " "
FOR !fin = !minfin TO !maxfin STEP 2
	SMPL 90 !maxfin+!horizonte+10
	FREEZE(_g{!fin}) G{!fin}.LINE
	%juntar = %juntar + " _g" + @STR(!fin)
	D G{!fin}
NEXT
GRAPH __g.merge {%JUNTAR}
D _g* 

y.LINE
SMPL @ALL


